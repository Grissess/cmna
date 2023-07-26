#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmna/cmna.h>

#define GROUND ((size_t)-1)

struct component {
	struct component *next;
	char type;
	unsigned line;
	size_t node_a;
	size_t node_b;
	void *data;
};

#define AS(tp, ptr) ((tp *)ptr)

int node_to_idx(unsigned, char *s, size_t *res);

int main() {
	struct cmna_circuit circuit;
	struct component *components = NULL, *top;
	enum cmna_error error;
	char type, node_a[32], node_b[32];
	size_t max_node = 0;
	size_t vss = 0, vs = 0;
	double value, *potentials, *currents;
	size_t node;
	unsigned lineno = 0;

	while(scanf(" %c %31s%31s", &type, node_a, node_b) != EOF) {
		lineno++;
		assert(top = malloc(sizeof(*top)));
		top->type = type;
		top->line = lineno;
		if(
				node_to_idx(lineno, node_a, &top->node_a) || \
				node_to_idx(lineno, node_b, &top->node_b)
		) {
			free(top);
			continue;
		}

		switch(type) {
			case 'R': case 'V': case 'I':
				assert(scanf("%lf", &value) == 1);
				assert(top->data = malloc(sizeof(double)));
				*AS(double, top->data) = value;
				break;

			default:
				fprintf(stderr, "Line %u: Unknown component %c\n", type);
				free(top);
				continue;
		}

		if(type == 'V') vss++;
		if(top->node_a != GROUND && top->node_a > max_node) max_node = top->node_a;
		if(top->node_b != GROUND && top->node_b > max_node) max_node = top->node_b;
		top->next = components;
		components = top;
	}

#define pl(v) (v == 1 ? "" : "s")

	fprintf(stderr, "Read %u line%s, %zd node%s, %zd source%s\n", lineno, pl(lineno), max_node + 1, pl(max_node + 1), vss, pl(vss));

#define check(e) if((error = (e))) { \
	fprintf(stderr, #e ": Error %d: %s\n", (int)error, cmna_strerror(error)); \
	return 1; \
}

	check(cmna_circuit_init(&circuit, CMNA_PREC_DOUBLE, max_node + 1, vss, NULL));

	while(top) {
		switch(top->type) {
			case 'R':
				if(top->node_a == GROUND && top->node_b == GROUND) {
					fprintf(stderr, "Line %u: useless resistor removed\n", top->line);
				} else if(top->node_a == GROUND) {
					check(cmna_circuit_add_conductance_to_ground(&circuit, top->node_b, 1.0 / *AS(double, top->data)));
				} else if(top->node_b == GROUND) {
					check(cmna_circuit_add_conductance_to_ground(&circuit, top->node_a, 1.0 / *AS(double, top->data)));
				} else {
					check(cmna_circuit_add_conductance(&circuit, top->node_a, top->node_b, 1.0 / *AS(double, top->data)));
				}
				break;

			case 'V':
				if(top->node_a != GROUND) {
					check(cmna_circuit_add_source_pos(&circuit, vs, top->node_a));
				}
				if(top->node_b != GROUND) {
					check(cmna_circuit_add_source_neg(&circuit, vs, top->node_b));
				}
				check(cmna_circuit_add_source_potential(&circuit, vs, *AS(double, top->data)));
				vs++;
				break;

			case 'I':
				if(top->node_a != GROUND) {
					check(cmna_circuit_add_current(&circuit, top->node_a, *AS(double, top->data)));
				}
				if(top->node_b != GROUND) {
					check(cmna_circuit_add_current(&circuit, top->node_b, -*AS(double, top->data)));
				}
				break;

			default:
				assert(0);
				break;
		}

		top = top->next;
	}

	fprintf(stderr, "Simulating...\n");

	check(cmna_circuit_solve(&circuit));

	check(cmna_circuit_node_potentials(&circuit, (void **)&potentials));
	check(cmna_circuit_source_currents(&circuit, (void **)&currents));
	for(node = 0; node < circuit.nodes; node++) {
		printf("%lf ", potentials[node]);
	}
	putchar('\n');
	for(vs = 0; vs < circuit.sources; vs++) {
		printf("%lf ", -currents[vs]);
	}
	putchar('\n');

	cmna_circuit_cleanup(&circuit);

	return 0;
}

int node_to_idx(unsigned lineno, char *s, size_t *res) {
	char *ep;
	if(!strcmp(s, "G")) {
		*res = GROUND;
	} else {
		*res = strtoul(s, &ep, 10);
		if(*ep) {
			fprintf(stderr, "Line %u: Failed to parse: %s not a valid base 10 string or G\n", lineno, s);
			return 1;
		}
	}
	return 0;
}
		

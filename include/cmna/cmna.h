#ifndef CMNA_H
#define CMNA_H

#include <stdint.h>
#include <stddef.h>

enum cmna_prec {
	CMNA_PREC_SINGLE,
	CMNA_PREC_DOUBLE
};

struct cmna_alloc {
	void *(*calloc)(size_t, size_t);
	void (*free)(void *);
};
extern struct cmna_alloc CMNA_ALLOC_MALLOC;

enum cmna_solvestate {
	CMNA_SS_UNINIT,
	CMNA_SS_SOLVED,
	CMNA_SS_DIRTY_MATRIX,
	CMNA_SS_DIRTY_KNOWNS,
	CMNA_SS_ERROR
};

struct cmna_solve {
	enum cmna_solvestate state;
	void *factorized;
	int *pivots;
	void *unknowns;
};

struct cmna_circuit {
	enum cmna_prec precision;
	size_t nodes;
	size_t sources;
	struct cmna_alloc alloc;
	void *matrix;
	void *knowns;
	struct cmna_solve solve;
};
#define cmna_circuit_matrix_size(c) ((c)->nodes + (c)->sources)
#define cmna_circuit_dirty(c, d) ((c)->solve.state = (c)->solve.state == CMNA_SS_UNINIT ? CMNA_SS_UNINIT : (d))
#define cmna_circuit_matrix_add(c, m, i, j, v) ((m)[(i) * cmna_circuit_matrix_size(c) + (j)] += (v))

enum cmna_error {
	CMNA_E_SUCCESS,
	CMNA_E_NO_MEM,
	CMNA_E_INVALID,
	CMNA_E_OOB,
	CMNA_E_SINGULAR,
	CMNA_E_NOT_READY,
	_CMNA_E_MAX_ERROR
};

enum cmna_error cmna_circuit_init(
		struct cmna_circuit *circuit, 
		enum cmna_prec precision,
		size_t nodes, size_t sources,
		struct cmna_alloc *alloc
);

void cmna_circuit_cleanup(struct cmna_circuit *circuit);

enum cmna_error cmna_circuit_solve(struct cmna_circuit *circuit);

enum cmna_error cmna_circuit_add_conductance_to_ground(
		struct cmna_circuit *circuit,
		size_t node,
		double cond
);

enum cmna_error cmna_circuit_add_conductance(
		struct cmna_circuit *circuit,
		size_t node_a, size_t node_b,
		double cond
);

enum cmna_error cmna_circuit_add_source(
		struct cmna_circuit *circuit,
		size_t source, size_t node, double value
);
#define cmna_circuit_add_source_pos(c, s, n) cmna_circuit_add_source((c), (s), (n), 1.0)
#define cmna_circuit_add_source_neg(c, s, n) cmna_circuit_add_source((c), (s), (n), -1.0)

enum cmna_error cmna_circuit_add_source_potential(
		struct cmna_circuit *circuit,
		size_t source, double potential
);

enum cmna_error cmna_circuit_add_current(
		struct cmna_circuit *circuit,
		size_t node, double current
);

enum cmna_error cmna_circuit_node_potentials(
		struct cmna_circuit *circuit,
		void **potentials
);

enum cmna_error cmna_circuit_source_currents(
		struct cmna_circuit *circuit,
		void **currents
);

size_t cmna_element_size(enum cmna_prec precision);

const char *cmna_strerror(enum cmna_error err);

#endif

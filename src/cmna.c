#include <stdlib.h>
#include <string.h>

#include <lapack.h>

#include <cmna/cmna.h>

#define CHECK_NODE(c, n) if((n) >= (c)->nodes) return CMNA_E_OOB
#define CHECK_SOURCE(c, n) if((n) >= (c)->sources) return CMNA_E_OOB

struct cmna_alloc CMNA_ALLOC_MALLOC = {
	.calloc = calloc,
	.free = free,
};

enum cmna_error cmna_circuit_init(
		struct cmna_circuit *circuit,
		enum cmna_prec precision,
		size_t nodes, size_t sources,
		struct cmna_alloc *alloc
) {
	if(!alloc) alloc = &CMNA_ALLOC_MALLOC;

	circuit->precision = precision;
	circuit->nodes = nodes;
	circuit->sources = sources;
	circuit->solve.state = CMNA_SS_UNINIT;
	circuit->alloc = *alloc;

	size_t elem_size = cmna_element_size(precision);
	size_t mtx_size = cmna_circuit_matrix_size(circuit);

	circuit->matrix = circuit->alloc.calloc(
			mtx_size * mtx_size, elem_size
	);
	if(!circuit->matrix) return CMNA_E_NO_MEM;
	circuit->knowns = circuit->alloc.calloc(
			mtx_size, elem_size
	);
	if(!circuit->knowns) goto out_free_mtx;

	goto ok;

out_free_mtx:
	circuit->alloc.free(circuit->matrix);
	return CMNA_E_NO_MEM;
	
ok:
	return CMNA_E_SUCCESS;
}

void cmna_circuit_cleanup(struct cmna_circuit *circuit) {
	circuit->alloc.free(circuit->matrix);
	circuit->alloc.free(circuit->knowns);
	if(circuit->solve.state != CMNA_SS_UNINIT) {
		circuit->alloc.free(circuit->solve.factorized);
		circuit->alloc.free(circuit->solve.pivots);
		circuit->alloc.free(circuit->solve.unknowns);
	}
}

enum cmna_error cmna_circuit_solve(struct cmna_circuit *circuit) {
	static char trans = 'N';
	static int one = 1;

	if(circuit->solve.state == CMNA_SS_SOLVED) return CMNA_E_SUCCESS;

	size_t elem_size = cmna_element_size(circuit->precision);
	size_t mtx_size = cmna_circuit_matrix_size(circuit);

	if(circuit->solve.state == CMNA_SS_UNINIT) {
		circuit->solve.factorized = circuit->alloc.calloc(
				mtx_size * mtx_size, elem_size
		);
		if(!circuit->solve.factorized) return CMNA_E_NO_MEM;
		circuit->solve.pivots = circuit->alloc.calloc(
				mtx_size, sizeof(int)
		);
		if(!circuit->solve.pivots) goto out_free_fac;
		circuit->solve.unknowns = circuit->alloc.calloc(
				mtx_size, elem_size
		);
		if(!circuit->solve.unknowns) goto out_free_piv;

		goto ok;

out_free_piv:
		circuit->alloc.free(circuit->solve.pivots);
out_free_fac:
		circuit->alloc.free(circuit->solve.factorized);
		return CMNA_E_NO_MEM;

ok:
		circuit->solve.state = CMNA_SS_DIRTY_MATRIX;
	}

	int mtx_size_i = mtx_size;  /* XXX overflow */
	int info;

	if(circuit->solve.state == CMNA_SS_DIRTY_MATRIX) {
		size_t bytes = mtx_size * mtx_size * elem_size;  /* XXX idem */

		memcpy(circuit->solve.factorized, circuit->matrix, bytes);
		switch(circuit->precision) {
			case CMNA_PREC_SINGLE:
				LAPACK_sgetrf(
						&mtx_size_i, &mtx_size_i,  /* M, N */
						circuit->solve.factorized, &mtx_size_i,  /* A, LDA */
						circuit->solve.pivots,  /* IPIV */
						&info  /* INFO */
				);
				break;

			case CMNA_PREC_DOUBLE:
				LAPACK_dgetrf(
						&mtx_size_i, &mtx_size_i,  /* M, N */
						circuit->solve.factorized, &mtx_size_i,  /* A, LDA */
						circuit->solve.pivots,  /* IPIV */
						&info  /* INFO */
				);
				break;

			default:
				return CMNA_E_INVALID;
		}

		circuit->solve.state = CMNA_SS_ERROR;
		if(info < 0) return CMNA_E_INVALID;
		if(info > 0) return CMNA_E_SINGULAR;
		circuit->solve.state = CMNA_SS_DIRTY_KNOWNS;
	}

	if(circuit->solve.state == CMNA_SS_DIRTY_KNOWNS) {
		size_t bytes = elem_size * mtx_size;
		
		memcpy(circuit->solve.unknowns, circuit->knowns, bytes);
		switch(circuit->precision) {
			case CMNA_PREC_SINGLE:
				LAPACK_sgetrs(
						&trans,  /* TRANS */
						&mtx_size_i, &one,  /* N, NRHS */
						circuit->solve.factorized, &mtx_size_i,  /* A, LDA */
						circuit->solve.pivots,  /* IPIV */
						circuit->solve.unknowns, &mtx_size_i, /* B, LDB */
						&info  /* INFO */
				);
				break;

			case CMNA_PREC_DOUBLE:
				LAPACK_dgetrs(
						&trans,  /* TRANS */
						&mtx_size_i, &one,  /* N, NRHS */
						circuit->solve.factorized, &mtx_size_i,  /* A, LDA */
						circuit->solve.pivots,  /* IPIV */
						circuit->solve.unknowns, &mtx_size_i, /* B, LDB */
						&info  /* INFO */
				);
				break;

			default:
				return CMNA_E_INVALID;
		}

		circuit->solve.state = CMNA_SS_ERROR;
		if(info < 0) return CMNA_E_INVALID;
		circuit->solve.state = CMNA_SS_SOLVED;
	}
	return CMNA_E_SUCCESS;
}

enum cmna_error cmna_circuit_add_conductance_to_ground(
		struct cmna_circuit *circuit,
		size_t node,
		double cond
) {
	CHECK_NODE(circuit, node);

	switch(circuit->precision) {
		case CMNA_PREC_SINGLE:
			cmna_circuit_matrix_add(
					circuit,
					(float *) circuit->matrix,
					node, node,
					cond
			);
			break;

		case CMNA_PREC_DOUBLE:
			cmna_circuit_matrix_add(
					circuit,
					(double *) circuit->matrix,
					node, node,
					cond
			);
			break;

		default:
			return CMNA_E_INVALID;
	}
	cmna_circuit_dirty(circuit, CMNA_SS_DIRTY_MATRIX);
	return CMNA_E_SUCCESS;
}

enum cmna_error cmna_circuit_add_conductance(
		struct cmna_circuit *circuit,
		size_t node_a, size_t node_b,
		double cond
) {
	enum cmna_error error;
	CHECK_NODE(circuit, node_a);
	CHECK_NODE(circuit, node_b);

	if((error = cmna_circuit_add_conductance_to_ground(circuit, node_a, cond)))
		return error;
	if((error = cmna_circuit_add_conductance_to_ground(circuit, node_b, cond)))
		return error;

	switch(circuit->precision) {
		case CMNA_PREC_SINGLE:
			{
				float *mtx = (float *)circuit->matrix;
				cmna_circuit_matrix_add(
						circuit, mtx,
						node_a, node_b,
						-cond
				);
				cmna_circuit_matrix_add(
						circuit, mtx,
						node_b, node_a,
						-cond
				);
			}
			break;

		case CMNA_PREC_DOUBLE:
			{
				double *mtx = (double *)circuit->matrix;
				cmna_circuit_matrix_add(
						circuit, mtx,
						node_a, node_b,
						-cond
				);
				cmna_circuit_matrix_add(
						circuit, mtx,
						node_b, node_a,
						-cond
				);
			}
			break;

		default:
			return CMNA_E_INVALID;
	}
	cmna_circuit_dirty(circuit, CMNA_SS_DIRTY_MATRIX);
	return CMNA_E_SUCCESS;
}

enum cmna_error cmna_circuit_add_source(
		struct cmna_circuit *circuit,
		size_t source, size_t node, double value
) {
	size_t si = source + circuit->nodes;
	CHECK_NODE(circuit, node);
	CHECK_SOURCE(circuit, source);

	switch(circuit->precision) {
		case CMNA_PREC_SINGLE:
			{
				float *mtx = (float *)circuit->matrix;
				cmna_circuit_matrix_add(
						circuit, mtx,
						si, node,
						value
				);
				cmna_circuit_matrix_add(
						circuit, mtx,
						node, si,
						value
				);
			}
			break;

		case CMNA_PREC_DOUBLE:
			{
				double *mtx = (double *)circuit->matrix;
				cmna_circuit_matrix_add(
						circuit, mtx,
						si, node,
						value
				);
				cmna_circuit_matrix_add(
						circuit, mtx,
						node, si,
						value
				);
			}
			break;

		default:
			return CMNA_E_INVALID;
	}
	cmna_circuit_dirty(circuit, CMNA_SS_DIRTY_MATRIX);
	return CMNA_E_SUCCESS;
}

enum cmna_error cmna_circuit_add_source_potential(
		struct cmna_circuit *circuit,
		size_t source, double potential
) {
	size_t si = source + circuit->nodes;
	CHECK_SOURCE(circuit, source);

	switch(circuit->precision) {
		case CMNA_PREC_SINGLE:
			((float *)circuit->knowns)[si] += potential;
			break;

		case CMNA_PREC_DOUBLE:
			((double *)circuit->knowns)[si] += potential;
			break;

		default:
			return CMNA_E_INVALID;
	}
	cmna_circuit_dirty(circuit, CMNA_SS_DIRTY_KNOWNS);
	return CMNA_E_SUCCESS;
}

enum cmna_error cmna_circuit_add_current(
		struct cmna_circuit *circuit,
		size_t node, double current
) {
	CHECK_NODE(circuit, node);

	switch(circuit->precision) {
		case CMNA_PREC_SINGLE:
			((float *)circuit->knowns)[node] += current;
			break;

		case CMNA_PREC_DOUBLE:
			((double *)circuit->knowns)[node] += current;
			break;

		default:
			return CMNA_E_INVALID;
	}
	cmna_circuit_dirty(circuit, CMNA_SS_DIRTY_KNOWNS);
	return CMNA_E_SUCCESS;
}

enum cmna_error cmna_circuit_node_voltages(
		struct cmna_circuit *circuit,
		void **voltages
) {
	if(circuit->solve.state != CMNA_SS_SOLVED) return CMNA_E_NOT_READY;
	if(!voltages) return CMNA_E_INVALID;
	*voltages = circuit->solve.unknowns;
	return CMNA_E_SUCCESS;
}

enum cmna_error cmna_circuit_source_currents(
		struct cmna_circuit *circuit,
		void **currents
) {
	size_t elem_size = cmna_element_size(circuit->precision);
	if(circuit->solve.state != CMNA_SS_SOLVED) return CMNA_E_NOT_READY;
	if(!currents) return CMNA_E_INVALID;
	*currents = ((uint8_t*)circuit->solve.unknowns) + elem_size*circuit->nodes;
	return CMNA_E_SUCCESS;
}

size_t cmna_element_size(enum cmna_prec precision) {
	switch(precision) {
		case CMNA_PREC_SINGLE: return sizeof(float);
		case CMNA_PREC_DOUBLE: return sizeof(double);
	}
	/* This return value SHOULD cause a fault at invocation */
	return (size_t)-1;
}

static const char *_strerrors[] = {
	"Success",
	"Out of memory",
	"Invalid argument",
	"Out of bounds",
	"Singular matrix",
	"Not yet ready"
};

const char *cmna_strerror(enum cmna_error err) {
	if(err < 0 || err >= _CMNA_E_MAX_ERROR) return "(invalid error code!)";
	return _strerrors[err];
}

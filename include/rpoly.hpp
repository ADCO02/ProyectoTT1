/**
 * @file rpoly.h
 * 
 * @brief Provides functions to find roots of real-coefficient polynomials.
 * 
 * This module implements methods to find all roots (real and complex) of
 * polynomials with real coefficients using a stateful approach.
 * 
 * It allows allocation and release of internal state needed for root finding,
 * and provides functions to compute roots given polynomial coefficients.
 * The polynomial is represented as:
 *   p(x) = p[0]*x^degree + p[1]*x^(degree-1) + ... + p[degree].
 * 
 * Functions support polynomials up to a specified maximum degree.
 */

#ifndef BEGIN_C_DECLS
#define BEGIN_C_DECLS

/**
 * @brief Opaque state handle for polynomial root finding.
 */
// Opaque state handle
struct RPoly_State;

/**
 * @brief Allocates internal state for polynomial root finding.
 * 
 * @param max_degree Maximum degree of polynomials supported.
 * 
 * @return Pointer to allocated RPoly_State.
 */
// Allocate state for finding roots of a real polynomial.
// Degrees up to max_degree are supported by returned state.
struct RPoly_State *real_poly_alloc(int max_degree);

/**
 * @brief Releases allocated polynomial root finding state.
 * 
 * @param s Pointer to state to be released.
 */
// Release state
void real_poly_release(struct RPoly_State *s);


/**
 * @brief Computes roots of a polynomial with real coefficients.
 * 
 * @param p Array of polynomial coefficients (length degree+1), with
 *          p[0] as coefficient for x^degree.
 * @param degree Degree of the polynomial.
 * @param state Pointer to pre-allocated polynomial root finding state.
 * @param zeror Output array to hold real parts of roots.
 * @param zeroi Output array to hold imaginary parts of roots.
 * 
 * @return Number of roots found (degree).
 */
// Find roots of a polynomial with real coefficients.
// p[] holds coefficients of the polynomial:
//   p(x) = p[0]*x^degree + ... + p[degree]
//
// Caller must have already allocated 'state'.
// Returns number of roots stored in zeror[], zeroi[].
int real_poly_roots_compute(const double p[], int degree,
                            struct RPoly_State *state,
                            double zeror[], double zeroi[]);

                            
/**
 * @brief Convenience function to compute roots without managing state.
 * 
 * Allocates and releases state internally.
 * 
 * @param p Array of polynomial coefficients.
 * @param degree Degree of the polynomial.
 * @param zeror Output array for real parts of roots.
 * @param zeroi Output array for imaginary parts of roots.
 * 
 * @return Number of roots found.
 */
// Convenience function.
// Find roots of a polynomial with real coefficients.
// p[] holds coefficients of the polynomial:
//   p(x) = p[0]*x^degree + ... + p[degree]
//
// Returns number of roots stored in zeror[], zeroi[]
static inline int real_poly_roots(const double p[], int degree,
                                  double zeror[], double zeroi[])
{
    struct RPoly_State *state = real_poly_alloc(degree);
    int nr = real_poly_roots_compute(p, degree, state, zeror, zeroi);
    real_poly_release(state);

    return nr;
}

#endif

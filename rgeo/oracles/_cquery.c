/* Compiled query of CertifiedCutIndex: the same arithmetic as CertifiedCutIndex.bounds, in C.
 * Built on first use by rgeo/oracles/cquery.py (cc -O2 -ffp-contract=off, so the results match the
 * Python query to the last bits up to libm), loaded with ctypes. */
#include <math.h>
#include <stdint.h>

typedef struct {
    int64_t n, depth, exact, interp; /* interp: 0 off, 1 parabola, 2 model */
    double Rn2, lam, f32rel;
    const int32_t *anc, *row, *depth_of; /* n x depth, n x depth, n */
    const double *xy, *V3;               /* n x 2, n x 3 */
    const int64_t *near_ptr;             /* n + 1 */
    const int32_t *near_ids;
    const float *near_d;
    const int64_t *L_off, *P_off, *M_off; /* per cut: offset of its label matrix, portal arrays, member rows */
    const int64_t *K;                     /* per cut: portals */
    const float *L;                       /* label matrices, row-major |M| x K, concatenated */
    const double *gap, *arc;              /* per portal (gap[j] = arc gap to portal j+1, inf across runs and at the end) */
    const int64_t *run;
    const float *foot_r, *foot_arc;       /* per (cut, member row) */
    const int64_t *foot_run;
    const int64_t *E_off, *E_n, *end_idx; /* per cut: run-end stretches (portal index, 2 x stretch arc) */
    const double *end_gap;
    const uint8_t *lb_ok;                 /* per cut: interval lower bound valid */
} CIndex;

static inline double hyp(double u, double xc, double r) { double d = u - xc; return sqrt(r * r + d * d); }

static double model_min(const double *x, const double *y, double xs, double r1, double xt, double r2) {
    double g0 = y[0] - hyp(x[0], xs, r1) - hyp(x[0], xt, r2);
    double g1 = y[1] - hyp(x[1], xs, r1) - hyp(x[1], xt, r2);
    double g2 = y[2] - hyp(x[2], xs, r1) - hyp(x[2], xt, r2);
    double d1 = (g1 - g0) / (x[1] - x[0]), d2 = (g2 - g1) / (x[2] - x[1]);
    double a = (d2 - d1) / (x[2] - x[0]);
    double b = d1 - a * (x[0] + x[1]);
    double c = g0 - a * x[0] * x[0] - b * x[0];
    double lo = x[0], hi = x[2], u, D, DD, M1, M2;
#define DER(uu) do { M1 = hyp((uu), xs, r1); M2 = hyp((uu), xt, r2); \
        D = ((uu) - xs) / M1 + ((uu) - xt) / M2 + 2 * a * (uu) + b; \
        DD = r1 * r1 / (M1 * M1 * M1) + r2 * r2 / (M2 * M2 * M2) + 2 * a; } while (0)
    DER(lo);
    if (D >= 0) {
        u = lo;
    } else {
        DER(hi);
        if (D <= 0) {
            u = hi;
        } else {
            u = x[1];
            for (int it = 0; it < 40; it++) {
                DER(u);
                if (D < 0) lo = u; else hi = u;
                double un = DD > 0 ? u - D / DD : 0.5 * (lo + hi);
                if (!(lo < un && un < hi)) un = 0.5 * (lo + hi);
                if (fabs(un - u) < 1e-9 * fmax(1.0, fabs(u))) { u = un; break; }
                u = un;
            }
        }
    }
#undef DER
    double base = hyp(u, xs, r1) + hyp(u, xt, r2);
    double val = base + a * u * u + b * u + c;
    if (val < base) val = base;
    return y[1] < val ? y[1] : val;
}

static double interp(const CIndex *ix, int64_t C, const double *f, int64_t K, int64_t i, int64_t rs, int64_t rt) {
    const int64_t *run = ix->run + ix->P_off[C];
    if (!(0 < i && i < K - 1) || !(run[i - 1] == run[i] && run[i] == run[i + 1])) return f[i];
    const double *arc = ix->arc + ix->P_off[C];
    double x[3] = {arc[i - 1], arc[i], arc[i + 1]}, y[3] = {f[i - 1], f[i], f[i + 1]};
    if (ix->interp == 2) {
        int64_t mo = ix->M_off[C];
        if (ix->foot_run[mo + rs] == run[i] && ix->foot_run[mo + rt] == run[i]) {
            double xs = ix->foot_arc[mo + rs], r1 = ix->foot_r[mo + rs];
            double xt = ix->foot_arc[mo + rt], r2 = ix->foot_r[mo + rt];
            if (fmin(r1, r2) < 2.0 * (x[2] - x[0])) return model_min(x, y, xs, r1, xt, r2);
        }
    }
    double d1 = (y[1] - y[0]) / (x[1] - x[0]), d2 = (y[2] - y[1]) / (x[2] - x[1]);
    double a2 = (d2 - d1) / (x[2] - x[0]);
    if (a2 <= 0) return y[1];
    double xv = 0.5 * (x[0] + x[1]) - d1 / (2.0 * a2);
    if (!(x[0] <= xv && xv <= x[2])) return y[1];
    double fv = y[0] + d1 * (xv - x[0]) + a2 * (xv - x[0]) * (xv - x[1]);
    double cap = 0.25 * fmin(y[0] - y[1], y[2] - y[1]);
    return fmax(fv, y[1] - cap);
}

static double chord(const CIndex *ix, int64_t s, int64_t t) {
    const double *a = ix->V3 + 3 * s, *b = ix->V3 + 3 * t;
    double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return ix->lam * sqrt(dx * dx + dy * dy + dz * dz);
}

/* out = (LB, UB, estimate); scratch holds at least max K doubles */
void cq_bounds(const CIndex *ix, int64_t s, int64_t t, double *out, double *f) {
    if (s == t) { out[0] = out[1] = out[2] = 0.0; return; }
    double dx = ix->xy[2 * s] - ix->xy[2 * t], dy = ix->xy[2 * s + 1] - ix->xy[2 * t + 1];
    if (dx * dx + dy * dy <= ix->Rn2) {
        int64_t lo = ix->near_ptr[s], hi = ix->near_ptr[s + 1], end = hi;
        while (lo < hi) { int64_t m = (lo + hi) >> 1; if (ix->near_ids[m] < t) lo = m + 1; else hi = m; }
        if (lo < end && ix->near_ids[lo] == t) {
            double d = ix->near_d[lo];
            out[0] = ix->exact ? d * (1.0 - ix->f32rel) : chord(ix, s, t);
            out[1] = out[2] = d;
            return;
        }
    }
    const int32_t *as = ix->anc + s * ix->depth, *at = ix->anc + t * ix->depth;
    int64_t k = -1;
    for (int64_t j = 0; j < ix->depth; j++) if (as[j] != at[j]) { k = j; break; }
    if (k < 0) k = ix->depth_of[s] < ix->depth_of[t] ? ix->depth_of[s] : ix->depth_of[t];
    int64_t C = as[k - 1], rs = ix->row[s * ix->depth + k - 1], rt = ix->row[t * ix->depth + k - 1];
    int64_t K = ix->K[C];
    const float *Ls = ix->L + ix->L_off[C] + rs * K, *Lt = ix->L + ix->L_off[C] + rt * K;
    const double *gap = ix->gap + ix->P_off[C];
    double ub = INFINITY, lb = INFINITY;
    int64_t imin = 0;
    int use_lb = ix->exact && ix->lb_ok[C];
    for (int64_t j = 0; j < K; j++) {
        f[j] = (double)Ls[j] + (double)Lt[j];
        if (f[j] < ub) { ub = f[j]; imin = j; }
        if (use_lb && j > 0 && isfinite(gap[j - 1])) {
            double m = 0.5 * (f[j - 1] + f[j]) - gap[j - 1];
            if (m < lb) lb = m;
        }
    }
    if (use_lb) {
        for (int64_t e = ix->E_off[C]; e < ix->E_off[C] + ix->E_n[C]; e++) {
            double m = f[ix->end_idx[e]] - ix->end_gap[e];
            if (m < lb) lb = m;
        }
    }
    if (!use_lb || !isfinite(lb)) lb = 0.0;
    lb *= 1.0 - 2 * ix->f32rel;
    double ch = chord(ix, s, t);
    if (ch > lb) lb = ch;
    double est = ix->interp ? interp(ix, C, f, K, imin, rs, rt) : ub;
    if (est < lb) est = lb;
    if (est > ub) est = ub;
    out[0] = lb; out[1] = ub; out[2] = est;
}

void cq_bounds_batch(const CIndex *ix, int64_t m, const int64_t *S, const int64_t *T, double *out, double *f) {
    for (int64_t i = 0; i < m; i++) cq_bounds(ix, S[i], T[i], out + 3 * i, f);
}

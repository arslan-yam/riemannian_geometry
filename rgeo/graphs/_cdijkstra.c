/* Multi-seed Dijkstra with early termination at the targets (binary heap, lazy deletion).
 * Built on first use by rgeo/graphs/cdijkstra.py and loaded with ctypes.
 *
 * Graph: CSR (indptr[n+1], indices, data), nonnegative weights.
 * Source i is a set of seeds seed_node[seed_ptr[i] .. seed_ptr[i+1]) with initial distances
 * seed_dist[...] (a point inside a face is the set of that face's boundary nodes at their straight
 * distances). Targets: if tgt_ptr is NULL every source uses targets[0 .. nt), else source i uses
 * targets[tgt_ptr[i] .. tgt_ptr[i+1]). out receives one distance per (source, target) in that layout
 * (INFINITY when unreachable). A source's search stops as soon as all of its targets are settled or
 * the smallest open distance exceeds `limit`: the settled distances are final, so the answers are the
 * same as a full Dijkstra. Returns 0, or -1 if memory runs out. */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct { double d; int64_t v; } item;

static inline void push(item *h, int64_t *sz, double d, int64_t v) {
    int64_t i = (*sz)++;
    while (i > 0) {
        int64_t p = (i - 1) >> 1;
        if (h[p].d <= d) break;
        h[i] = h[p];
        i = p;
    }
    h[i].d = d; h[i].v = v;
}

static inline item pop(item *h, int64_t *sz) {
    item top = h[0], last = h[--(*sz)];
    int64_t i = 0, n = *sz;
    for (;;) {
        int64_t c = 2 * i + 1;
        if (c >= n) break;
        if (c + 1 < n && h[c + 1].d < h[c].d) c++;
        if (h[c].d >= last.d) break;
        h[i] = h[c];
        i = c;
    }
    if (n > 0) h[i] = last;
    return top;
}

int64_t cd_multi(int64_t n, const int64_t *indptr, const int32_t *indices, const double *data,
                 int64_t ns, const int64_t *seed_ptr, const int64_t *seed_node, const double *seed_dist,
                 int64_t nt, const int64_t *tgt_ptr, const int64_t *targets, double limit, double *out) {
    int64_t max_seeds = 0;
    for (int64_t i = 0; i < ns; i++)
        if (seed_ptr[i + 1] - seed_ptr[i] > max_seeds) max_seeds = seed_ptr[i + 1] - seed_ptr[i];
    double *dist = malloc((size_t)n * sizeof(double));
    uint8_t *done = malloc((size_t)n);
    int32_t *mark = calloc((size_t)n, sizeof(int32_t)); /* mark[v] = i+1 when v is a target of source i */
    int64_t cap = indptr[n] + max_seeds + 16;
    item *heap = malloc((size_t)cap * sizeof(item));
    int64_t *touched = malloc((size_t)n * sizeof(int64_t));
    if (!dist || !done || !mark || !heap || !touched) {
        free(dist); free(done); free(mark); free(heap); free(touched);
        return -1;
    }
    for (int64_t v = 0; v < n; v++) { dist[v] = INFINITY; done[v] = 0; }
    for (int64_t i = 0; i < ns; i++) {
        const int64_t *T = tgt_ptr ? targets + tgt_ptr[i] : targets;
        int64_t m = tgt_ptr ? tgt_ptr[i + 1] - tgt_ptr[i] : nt;
        double *o = tgt_ptr ? out + tgt_ptr[i] : out + i * nt;
        int64_t need = 0, ntouch = 0, sz = 0;
        for (int64_t j = 0; j < m; j++)
            if (mark[T[j]] != (int32_t)(i + 1)) { mark[T[j]] = (int32_t)(i + 1); need++; }
        for (int64_t k = seed_ptr[i]; k < seed_ptr[i + 1]; k++) {
            int64_t v = seed_node[k];
            double d = seed_dist[k];
            if (d < dist[v]) {
                if (dist[v] == INFINITY) touched[ntouch++] = v;
                dist[v] = d;
                push(heap, &sz, d, v);
            }
        }
        while (sz > 0 && need > 0) {
            item it = pop(heap, &sz);
            int64_t u = it.v;
            if (done[u] || it.d > dist[u]) continue;
            if (it.d > limit) break;
            done[u] = 1;
            if (mark[u] == (int32_t)(i + 1)) need--;
            for (int64_t e = indptr[u]; e < indptr[u + 1]; e++) {
                int64_t v = indices[e];
                double nd = it.d + data[e];
                if (nd < dist[v]) {
                    if (dist[v] == INFINITY) touched[ntouch++] = v;
                    dist[v] = nd;
                    if (sz >= cap) { /* cannot happen with lazy deletion bound, kept for safety */
                        item *h2 = realloc(heap, (size_t)(2 * cap) * sizeof(item));
                        if (!h2) { free(dist); free(done); free(mark); free(heap); free(touched); return -1; }
                        heap = h2; cap *= 2;
                    }
                    push(heap, &sz, nd, v);
                }
            }
        }
        for (int64_t j = 0; j < m; j++) o[j] = done[T[j]] ? dist[T[j]] : INFINITY;
        for (int64_t k = 0; k < ntouch; k++) { dist[touched[k]] = INFINITY; done[touched[k]] = 0; }
    }
    free(dist); free(done); free(mark); free(heap); free(touched);
    return 0;
}

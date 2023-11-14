#ifndef TRAVERSAL_HH
#define TRAVERSAL_HH

#include <assert.h>
#include <stdint.h>
#include <climits>
#include <vector>
#include <queue>          // std::priority_queue
#include <algorithm>      // std::min
#include <functional>

#include "mgraph.hh"

/**
 * Classical graph traversals: BFS, Dijkstra, ...
 * Example :
 *    mgraph<> g;
 *    traversal<graph> trav(n);
 *    int u=0; v=1;
 *    trav.dijkstra(g, u);
 *    int dist_u_v = trav.dist(v);
 *    trav.bfs(g, u);
 *    int hops_u_v = trav.dist(v);
 */

enum dfs_option { NONE = 0, SCC = 1, CHECK_DAG = 2 };

template<typename G,  // Graph type with G::vertex = int
         typename WL = int64_t, // long weight for summing weights, 
                     // with appropriate casting for the following constants:
         int64_t max_weight = INT64_MAX, int64_t zero_weight = 0>
class traversal {
public:
    typedef int V;
    typedef typename G::weight W;
    typedef WL long_weight;
private:
    typedef edge::dst_wgt<V,WL> wl_head;
    
    static inline bool wl_head_further (const wl_head &u, const wl_head &v) {
        return u.wgt > v.wgt;
    }

    std::vector<wl_head> q_vec_;
    std::vector<V> q_v_, scc_stack_;
    std::priority_queue<wl_head,
                        std::vector<wl_head>,
                        std::function<bool(const wl_head &, const wl_head &)>
                        > queue_;
    std::vector<V> visit_;
    std::vector<bool> visited_, in_scc_stack_;
    std::vector<int> lowlink_, visit_end_; // for DFS
    std::vector<V> parent_;
    std::vector<WL> dist_, tree_ecc_; // for Dijkstra and BFS
    
    std::vector<int> size_; // for tree centroid, and component sizes
    
    int n_, nvis_;
    
public:

    const V not_vertex = G::not_vertex;

    // Handle graphs with vertices in [0..n-1].
    traversal (int n)
        : queue_(wl_head_further, q_vec_), dist_(n, max_weight),
          visit_(n), visited_(n, false),
          size_(n, 0), parent_(n, n), n_(n), nvis_(0),
          in_scc_stack_(n, false), lowlink_(n, n), visit_end_(n, n),
          tree_ecc_(n, zero_weight)
    {
        q_vec_.reserve(n_); q_v_.reserve(n_);
    }

    int n() const { return n_; }
    int nvis() const { return nvis_; }
    V visit(int i) const { assert(i < nvis_); return visit_[i]; }
    bool visited(V u) { return visited_[u]; }
    V parent(V u) const { return parent_[u]; }
    WL dist(V u) const { return dist_[u]; }
    const std::vector<WL>& distances() const { return dist_; }
    int size(V u) const { return size_[u]; }
    WL tree_ecc(V u) const { return tree_ecc_[u]; }
    V first_visited (int i_start = 0) const { return visit_[i_start]; }
    V last_visited () const { return visit_[nvis_ - 1]; }
    
    void clear(WL dft_wgt = max_weight, int n = 0) {
        int up_to = std::max(n, nvis_);
        for (int i = 0; i < nvis_; ++i)  {
            V u = visit_[i];
            visit_[i] = not_vertex;
            visited_[u] = false;
            in_scc_stack_[u] = false;
            lowlink_[u] = n_;
            visit_end_[u] = n_;
            dist_[u] = dft_wgt;
            parent_[u] = not_vertex;
            size_[u] = 0;
            tree_ecc_[u] = zero_weight;
        }
        q_vec_.clear();
        scc_stack_.clear();
        nvis_ = 0;
    }

    static bool visit_all(V v, WL d, V p, WL dp) { return true; }
    
    // returns number of visited nodes
    int bfs(const G &g, V s,
            std::function<bool(V, WL, V, WL)> filtr
               = [](V v, WL d, V p, WL dp) { return true; },
            std::vector<V> more_sources = {}, WL step_weight = 1) {
        assert(g.n() <= n_);
        int nvis0 = nvis_;
        int head = 0, tail = 0;
        more_sources.push_back(s);
        for (int s : more_sources) {
            dist_[s] = zero_weight;
            parent_[s] = s;
            q_vec_[tail++] = wl_head(s, zero_weight);
        }

        while (head < tail) {
            const wl_head &u_du = q_vec_[head++];
            WL du = u_du.wgt;
            V u = u_du.dst;
            if ( ! visited_[u]) {
                visit_[nvis_++] = u;
                visited_[u] = true;
                WL dv = du + step_weight;
                for (V v : g[u]) {
                    if (dist_[v] == max_weight && filtr(v, dv, u, du)) {
                        dist_[v] = dv;
                        parent_[v] = u;
                        q_vec_[tail++] = wl_head(v, dv);
                    }
                }
            }
        }
        return nvis_ - nvis0;
    }

    // returns number of visited nodes
    int dijkstra (const G &g, V s,
                 std::function<bool(V, WL, V, WL)> filtr
                  = [](V v, WL d, V p, WL dp) { return true; },
                  std::vector<V> more_sources = {}) {
        assert(g.n() <= n_);
        int nvis0 = nvis_;
        more_sources.push_back(s);
        for (int s : more_sources) {
            dist_[s] = zero_weight;
            parent_[s] = s;
            queue_.push(wl_head(s, zero_weight));
        }
            
        while ( ! queue_.empty()) {
            const wl_head &u_du = queue_.top();
            WL du = u_du.wgt;
            V u = u_du.dst;
            queue_.pop();
            if (du == dist_[u] && ! visited_[u]) {
                visit_[nvis_++] = u;
                visited_[u] = true;
                for (auto e : g[u]) {
                    WL dv = du + e.wgt;
                    V v = e.dst;
                    if((! visited_[v]) && dv < dist_[v] && filtr(v, dv, u, du)){
                        dist_[v] = dv;
                        parent_[v] = u;
                        queue_.push(wl_head(v, dv));
                    }
                }
            }
        }
        return nvis_ - nvis0;
    }

    // returns number of visited nodes
    int dijkstra_i (const G &g, V s,
                    std::function<bool(V, WL, V, WL, int)> filtr
                    = [](V v, WL d, V p, WL dp, int iv) { return true; },
                  std::vector<V> more_sources = {}) {
        assert(g.n() <= n_);
        int nvis0 = nvis_;
        more_sources.push_back(s);
        for (int s : more_sources) {
            dist_[s] = zero_weight;
            parent_[s] = s;
            queue_.push(wl_head(s, zero_weight));
        }
            
        while ( ! queue_.empty()) {
            const wl_head &u_du = queue_.top();
            WL du = u_du.wgt;
            V u = u_du.dst;
            queue_.pop();
            if (du == dist_[u] && ! visited_[u]) {
                visit_[nvis_++] = u;
                visited_[u] = true;
                int iv = 0;
                for (auto e : g[u]) {
                    WL dv = du + e.wgt;
                    V v = e.dst;
                    if((! visited_[v]) && dv < dist_[v]
                       && filtr(v, dv, u, du, iv)){
                        dist_[v] = dv;
                        parent_[v] = u;
                        queue_.push(wl_head(v, dv));
                    }
                    ++iv;
                }
            }
        }
        return nvis_ - nvis0;
    }

    // returns number of visited nodes
    int a_star (const G &g, V s, V t,
                std::function<WL(V)> dist_to_t_lb,
                std::function<bool(V, WL, V, WL)> filtr
                = [](V v, WL d, V p, WL dp) { return true; }) {
        assert(g.n() <= n_);
        int nvis0 = nvis_;
        dist_[s] = zero_weight;
        parent_[s] = s;
        queue_.push(wl_head(s, dist_to_t_lb(s)));

        while ( ! queue_.empty()) {
            const wl_head &u_dut = queue_.top();
            WL d_ut = u_dut.wgt;
            V u = u_dut.dst;
            queue_.pop();
            if (d_ut == dist_[u] + dist_to_t_lb(u) && ! visited_[u]) {
                visit_[nvis_++] = u;
                visited_[u] = true;
                if (u == t) { break; }
                for (auto e : g[u]) {
                    WL dv = dist_[u] + e.wgt;
                    V v = e.dst;
                    if((! visited_[v]) && dv < dist_[v]
                       && filtr(v, dv, u, dist_[u])){
                        dist_[v] = dv;
                        parent_[v] = u;
                        queue_.push(wl_head(v, dv + dist_to_t_lb(v)));
                    }
                }
            }
        }
        return nvis_ - nvis0;
    }

    void clear_a_star(const G &g, WL dft_wgt = max_weight, int n = 0) {
        int up_to = std::max(n, nvis_);
        for (int i = 0; i < nvis_; ++i)  {
            V u = visit_[i];
            visit_[i] = not_vertex;
            visited_[u] = false;
            lowlink_[u] = n_;
            visit_end_[u] = n_;
            dist_[u] = dft_wgt;
            parent_[u] = not_vertex;
            size_[u] = 0;
            tree_ecc_[u] = zero_weight;
            for (V v : g[u]) {
                dist_[v] = dft_wgt;
                parent_[v] = not_vertex;
            }
        }
        q_vec_.clear();
        nvis_ = 0;
    }


    // returns number of visited nodes
    int dfs (const G &g, V s, dfs_option opt=NONE) {
        assert(g.n() <= n_);
        int nvis0 = nvis_, vis_end = nvis_;
        parent_[s] = s;
        int tail = 0;
        q_v_[tail++] = s;

        while (tail > 0) {
            V u = q_v_[--tail];
            /* lowlink_[u] (see Tarjan alg. for strongly connected components)
             * is the smallest visit number of a node accessible from u through
             * a sequence of forward edges plus eventually one backward edge. 
             * Here we consider more paths, possibly with several backward
             * edges; the important point is that all paths considered by
             * Tarjan algorithm are considered here, and the first node visited
             * in a strongly connected component (scc) will thus be the only 
             * node of the component that has its own visit  number as lowlink 
             * here also. (Not storing the visit number of nodes allows to
             * avoid allocating one more array).
             */
            if ( ! visited_[u]) { // begin visit of [u]
                lowlink_[u] = nvis_;
                visit_[nvis_++] = u;
                visited_[u] = true;
                if (opt == SCC) {
                    scc_stack_.push_back(u);
                    in_scc_stack_[u] = true;
                }
                
                // add u to stack for detecting end of visit when popped again
                if (tail >= q_v_.size()) q_v_.push_back(not_vertex);
                q_v_[tail++] = u;
                
                for (V v : g[u]) {
                    if ( ! visited_[v]) {
                        // possible forward edge
                        parent_[v] = u; // possibly overwritten
                        if (tail >= q_v_.size()) q_v_.push_back(not_vertex);
                        q_v_[tail++] = v;
                    } else if (visited_[v]) {
                        // backward edge
                        if (opt == SCC && in_scc_stack_[v]) {
                            //lowlink_[u] = std::min(lowlink_[u], visit_beg[v]);
                            lowlink_[u] = std::min(lowlink_[u], lowlink_[v]);
                        } else if (opt == CHECK_DAG && visit_end_[v] == n_) {
                            throw "cycle";
                        }
                    }
                }
            } else if (visit_end_[u] == n_) { // end visit of [u]
                visit_end_[u] = vis_end++;
                if (opt == SCC) {
                    V p = parent_[u];
                    int lowlink_u = lowlink_[u];
                    lowlink_[p] = std::min(lowlink_[p], lowlink_u);
                    assert(0 <= lowlink_u && lowlink_u < nvis_);
                    if (visit_[lowlink_u] == u) {
                        // u is the first visited node in its scc,
                        // pop scc and
                        // set all lowlinks in the scc to u's visit nb
                        while ( ! scc_stack_.empty() ) {
                            V v = scc_stack_.back();
                            scc_stack_.pop_back();
                            lowlink_[v] = lowlink_u;
                            in_scc_stack_[v] = false;
                            if (v == u) break; // all scc scanned
                        }
                    }
                }
            }
        }
        return nvis_ - nvis0;
    }

    // Compute strongly connected components of the input graph and returns
    // the number of components found.
    int strongly_connected_components(const G &g) {
        clear();
        for (V s = 0; s < g.n(); ++s) {
            if ( ! visited_[s]) {
                dfs(g, s, SCC);
            }
        }
        std::cerr <<"nvis="<< nvis_ <<"\n";
        int nb = 0;
        for (int i = 0; i < nvis_; ++i) {
            V u = visit_[i];
            int ll = lowlink_[u];
            assert(0 <= ll && ll < nvis_);
            ++(size_[ll]);
            if (ll == i) ++nb;
        }
        return nb;
    }

    // Returns the scc number of the strongly connected component of [v],
    // asserts that [strongly_connected_components(g)] has been called
    // and that [v <= g.n()].
    // The return number is in the range 0 .. g.n() - 1.
    int scc_number(V v) {
        assert(v < nvis_);
        return lowlink_[v];
    }
    int scc_largest() {
        assert(nvis_ > 0);
        int m = 0;
        for (int i = 1; i < nvis_; ++i) {
            if (size_[i] > size_[m]) m = i;
        }
        return m;
    }
    int scc_size(int i) { return size_[i]; }
    std::vector<int> scc_vector() {
        return lowlink_;
    }
    // Returns a node of scc number i.
    V scc_node(int i) {
        assert(i < nvis_);
        return visit_[i];
    }
    // Returns the jth pseudo-random node of scc number i.
    V scc_random_node(int i, int j) {
        assert(i < nvis_);
        int i_r = i + (127L + 1237L * j) % size_[i]; // pseudo random jth
        V v = visit_[i_r];
        assert(lowlink_[v] == i);
        return v;
    }
    
    // Returns a topological ordering [ord] of [g]
    // such that for all edge [uv] of [g], [u] appears before [v] in [ord].
    std::vector<int> topological_ordering (G &g) {
        clear();
        for (V u : g) {
            if ( ! visited_[u]) {
                dfs(g, u, true);
            }
        }
        int n = g.n();
        std::vector<int> ord;
        ord.reserve(n);
        for (V u : g) ord[n - 1 - visit_end_[u]] = u;
        return ord;
    }


    
    /** The following functions assume that WL is an int type. */
    
    // returns number of visited nodes, assumes 
    int max_card_search (const G &g, V s,
                 std::function<bool(V, V)> filtr
                 = [](V v, V p) { return true; }) {
        assert(g.n() <= n_);
        int nvis0 = nvis_;
        const WL nvis_max = n_;
        dist_[s] = nvis_max;
        // dist_[u] will contain nvis_max - (nb of visited neigh)
        parent_[s] = s;
        queue_.push(wl_head(s, dist_[s]));

        while ( ! queue_.empty()) {
            const wl_head &u_du = queue_.top();
            WL nb_neighb_vis = u_du.wgt;
            V u = u_du.dst;
            queue_.pop();
            if (nb_neighb_vis == dist_[u] && ! visited_[u]) {
                visit_[nvis_++] = u;
                visited_[u] = true;
                for (V v : g[u]) {
                    if ( (! visited_[v]) && filtr(v, u)) {
                        dist_[v] = std::min(dist_[v], nvis_max) - 1;
                        queue_.push(wl_head(v, dist_[v]));
                    }
                }
            }
        }
        return nvis_ - nvis0;
    }

    // Asserts that the graph is a tree (with symmetric links).
    int tree_size(const G &g, int root,
                  std::function<int(V)> node_size = [](V v) { return 1; }) {
        return tree_size(g, root, root);
    }
    int tree_size(const G &g, int u, int par,
                  std::function<int(V)> node_size = [](V v) { return 1; }) {
        int s = 1; // for u
        size_[u] = node_size(u);
        assert (size_[u] > 0);
        for (V v : g[u]) {
            if (size_[v] <= 0) {
                s += tree_size(g, v, u);
            } else assert(v == par || v == u); // Check that we are in a tree.
        }
        size_[u] = s;
        return s ;
    }

    void tree_size_clear(const G &g, int r) { tree_size_clear(g, r, r); }
    void tree_size_clear(const G &g, int u, int par) {
        if (size_[u] > 0) {
            size_[u] = 0;
            for (V v : g[u]) {
                if (size_[v] > 0) {
                    tree_size_clear(g, v, u);
                } else assert(v == par || v == u); // Check tree.
            }
        }
    }

    // Asserts that tree_size(g, root) has been called.
    V centroid(const G &g, int r, int n_threshold) {
        return centroid(g, r, r, n_threshold);
    }
    V centroid(const G &g, int r, int par, int n_threshold) {
        for (V v : g[r]) {
            if (v != par && v != r && size_[v] > n_threshold)
                return centroid (g, v, r,  n_threshold);
        }
        return r;
    }

    /** 
     * Eccentricities in the undirected tree given by parent_
     * and having root visit[i_r].
     *
     * We use the property that in a tree rooted at r,
     * each node u has eccentricity max(d(u,f), d(u,f2)) where f
     * is the fursthest node from r and f2 is the furthest node from f.
     */
    void tree_eccentricities(int i_r=0) {
        V r = visit_[i_r];
        V f = visit_[nvis_ - 1];
        
        // set d_f_f2[u] to max(d(u,f), d(u,f2)) and return node with max dist
        // assumes dist_ contains distances from r and that f is furthest from r
        auto set_dist = [this, i_r, r, f](std::vector<WL> &d_f_f2, V f2) {
            for (int i = i_r; i < nvis_; ++i) {
                d_f_f2[visit_[i]] = max_weight;
            }
            // dist from f in f--r path
            for (V u = f; true; u = parent_[u]) {
                d_f_f2[u] = dist_[f] - dist_[u];
                if (u == r) break;
            }
            if (f2 != f) {
                assert (dist_[f2] <= dist_[f]);
                // find the common ancestor a of f and f2
                V a = f2;
                while (d_f_f2[a] == max_weight) {
                    d_f_f2[a] = dist_[f2] - dist_[a]; // set dist from f2
                    a = parent_[a];
                }
                // update f--a path
                WL a_f2 = dist_[f2] - dist_[a];
                for (V u = f; u != a; u = parent_[u]) {
                    d_f_f2[u] = std::max(d_f_f2[u], a_f2 + dist_[u] - dist_[a]);
                }
                // update f2--a path
                WL  a_f = dist_[f] - dist_[a];
                for (V u = f2; u != a; u = parent_[u]) {
                    d_f_f2[u] = std::max(d_f_f2[u], a_f + dist_[u] - dist_[a]);
                }
            }
            // propagate to the tree
            V far = r;
            for (int i = i_r; i < nvis_; ++i) {
                V u = visit_[i];
                if (d_f_f2[u] == max_weight) {
                    V p = parent_[u];
                    d_f_f2[u] = d_f_f2[p] + dist_[u] - dist_[p];
                    if (d_f_f2[u] > d_f_f2[far]) far = u;
                }
            }
            return far;
        };

        V f2 = set_dist(tree_ecc_, f);
        set_dist(tree_ecc_, f2);
    }

    
    const std::vector<WL> & dist_from(V f, int i_r=0) {
        V r = visit_[i_r];
        std::vector<WL> &d_f = tree_ecc_;
        for (int i = i_r; i < nvis_; ++i) {
            d_f[visit_[i]] = max_weight;
        }
        // dist from f in f--r path
        for (V u = f; true; u = parent_[u]) {
            d_f[u] = dist_[f] - dist_[u];
            if (u == r) break;
        }
        // propagate to the tree
        for (int i = i_r; i < nvis_; ++i) {
            V u = visit_[i];
            if (d_f[u] == max_weight) {
                V p = parent_[u];
                d_f[u] = d_f[p] + dist_[u] - dist_[p];
            }
        }
        return tree_ecc_;
    }

    

    

}; // traversal

#endif // TRAVERSAL_HH

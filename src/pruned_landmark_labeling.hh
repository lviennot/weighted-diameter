#ifndef PRUNED_LANDMARK_LABELING_HH
#define PRUNED_LANDMARK_LABELING_HH

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <xmmintrin.h> // prefetch

#include "edge.hh"

template<typename G,
         typename WL = int64_t, // long weights
         int64_t max_weight = INT64_MAX, int64_t zero_weight = 0>
class pruned_landmark_labeling {

private:
    typedef typename G::vertex V;
    typedef typename G::weight W;
    typedef typename edge::src_dst_wgt<V, W> edge;
    
    struct label_t {
        std::vector<V> out_v; // out hubs
        std::vector<WL> out_d;   // dist to them
        std::vector<V> in_v;  // in hubs
        std::vector<WL> in_d;    // dist from them
    }; // __attribute__((aligned(64)));

    static inline void * memalign (size_t bytes) {
        void *ptr;
        assert(0 == posix_memalign(&ptr, 64, bytes));
        return ptr;
    }
    
    int n_; // number of nodes
    std::vector<label_t> index_;


public:
    pruned_landmark_labeling() : n_(0) {}

    inline WL distance(int u, int v) {
        const label_t &lab_u = index_[u];
        const label_t &lab_v = index_[v];

        WL d = max_weight;

        //__mm_prefetch(&(lab_u.out_v[0], _MM_HINT_T0));
        //__mm_prefetch(&(lab_v.in_v[0], _MM_HINT_T0));

        for (int iu = 0, iv = 0 ; ; ) {
            V xu = lab_u.out_v[iu], xv = lab_v.in_v[iv];
            if (xu == xv) { // common hub
                if (xu == n_) break; // sentinel
                WL d_uxv = lab_u.out_d[iu] + lab_v.in_d[iv];
                if (d_uxv < d) d = d_uxv;
                ++iu;
                ++iv;
            } else {
                if (xu < xv) ++iu;
                else ++iv;
            }
        }

        return d;
    }

    void print_stats(std::ostream &cout) {
        long int fwd = 0, bwd = 0, fmax = 0, bmax = 0;
        for (int u = 0; u < n_; ++u) {
            if (index_[u].out_v.size() > fmax) fmax = index_[u].out_v.size();
            if (index_[u].in_v.size() > bmax) bmax = index_[u].in_v.size();
            fwd += index_[u].out_v.size();
            bwd += index_[u].in_v.size();
        }
        cout << n_ << " labels\n";
        if (n_ > 0) cout << "  avg fwd = " << (fwd / n_)
                              << "  avg bwd = " << (bwd / n_)
                              << "  max fwd = " << fmax
                              << "  max bwd = " << bmax
                              << "\n";
    }

    pruned_landmark_labeling(const G &g,
                    std::vector<V> rank_order = std::vector<int>()) {
        construct_index(g, rank_order);
    }

    void construct_index(const G &g, std::vector<V> rank_order) {

        // Number of vertices
        n_ = g.n();

        // Allocate index
        index_ = std::vector<label_t>(n_);
        for (int u = 0; u < n_; ++u) {
            label_t &lab_u = index_[u];
            lab_u.in_v.push_back(n_);
            lab_u.in_d.push_back(max_weight);
            lab_u.out_v.push_back(n_);
            lab_u.out_d.push_back(max_weight);
        }

        // Rank of a vertex and graph with ranks as ID
        std::vector<int> rank_of(n_);
        std::vector<edge> edges(g.m());
        {
            if (rank_order.size() != n_) { // default order
                assert(rank_order.size() == 0);
                rank_order = std::vector<int>(n_);
                for (int i = 0; i < n_; ++i) rank_order[i] = i;
            }
            for (int i = 0; i < n_; ++i) {
                rank_of[rank_order[i]] = i;
            }
            for (int u = 0, i = 0; u < n_; ++u)
                for (auto e : g[u])
                    edges[i++] = edge(rank_of[u], rank_of[e.dst], e.wgt);
        }
        G succ(g.n(), edges);
        for (int e = 0; e < g.m(); ++e) std::swap(edges[e].src, edges[e].dst);
        G prec(g.n(), edges);

        // Pruned labeling : make two pruned Dijkstras (forward and backward)
        // from each root r (by rank order), and add r as hub to labels
        // of visited nodes.
        {
            traversal<G, WL, max_weight, zero_weight> trav(n_);
            
            int64_t sum_nvis = 0, last_nvis = 0, last_r = -1;
            
            // Visit by rank order:
            for (int r = 0; r < n_; ++r) {
                int ro = rank_order[r];

                // forward dijkstra:
                {
                    trav.dijkstra(succ, r,
                               [this,&rank_order,ro](V v, WL dv, V prt, WL dprt)
                               -> bool {
                                   return dv < distance(ro, rank_order[v]);
                               });

                    int nvis = trav.nvis();
                    for (int i = 0; i < nvis; ++i) {
                        int u = trav.visit(i);
                        int uo = rank_order[u];
                        label_t &lab_u = index_[uo];
                        lab_u.in_v.back() = r;
                        lab_u.in_d.back() = trav.dist(u);
                        lab_u.in_v.push_back(n_);
                        lab_u.in_d.push_back(max_weight);
                    }
                    sum_nvis += nvis; last_nvis += nvis;
                    trav.clear();
                }
                    
                // backward dijkstra:
                {
                    trav.dijkstra(prec, r,
                               [this,&rank_order,ro](V v, WL dv, V prt, WL dprt)
                               -> bool {
                                   return dv < distance(rank_order[v], ro);
                               });
                    
                    int nvis = trav.nvis();
                    for (int i = 0; i < nvis; ++i) {
                        int u = trav.visit(i);
                        int uo = rank_order[u];
                        label_t &lab_u = index_[uo];
                        lab_u.out_v.back() = r;
                        lab_u.out_d.back() = trav.dist(u);
                        lab_u.out_v.push_back(n_);
                        lab_u.out_d.push_back(max_weight);
                    }
                    sum_nvis += nvis; last_nvis += nvis;
                    trav.clear();
                }

                // progress:
                if (r <= 10 || (r <= 100 && r % 10 == 0) || r % 100 == 0) {
                    std::cerr << r  << " avg_nvis="<< (sum_nvis / (2*(r+1)))
                              <<" lst_nvis="<< (last_nvis / 2 / (r - last_r))
                              << " " << " ro=" << ro << " n=" << n_ <<" "
                              << "avg_hs=" << (sum_nvis / (2*n_))
                              <<"         \r";
                    std::cerr.flush();
                    last_r = r;
                    last_nvis = 0;
                }

            }

            std::cerr << "\n"; std::cerr.flush();
            
        } // Pruned labeling
    }
    
};

#endif //PRUNED_LANDMARK_LABELING_HH

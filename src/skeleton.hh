#ifndef SKELETON_HH
#define SKELETON_HH

#include <assert.h>
#include <vector>
#include <queue>          // std::priority_queue


template<typename T, // Traversal.
         typename WL = int64_t, // long weigth (not necessarily T::long_weigth)
         int64_t zero_weight = 0>
class skeleton {
private:

    typedef int V;
    typedef typename T::W W;

    std::vector<V> furthest_;
    std::vector<V> visit_;
    std::vector<WL> dist_;
    std::vector<int> nsons_;
    int n_; // max tree size
    int n_sk; // skeleton size

public:

    skeleton(int n)
        : furthest_(n), visit_(n), dist_(n), nsons_(n),
          n_(n), n_sk(0) { assert(n > 0); }

    void clear() {
        for (int i = 0; i < n_sk; ++i) nsons_[visit_[i]] = 0;
        n_sk = 0;
    }

    int size() { return n_sk; }

    static WL traversal_metric (V u, V v, W w) { return w; }
    static WL log_metric (V u, V v, W w) { return std::log(1 + 2*w); }
    static WL hop_count_metric (V u, V v, W w) { return 1; }
    
    void of_traversal (const T &trav,
                       WL alpha_p = 1, WL alpha_q = 2, // alpha = 1/2
                       std::function<WL(V, V, W)> edge_metr
                       = hop_count_metric) {
        // alpha = alpha_p / alpha_q
        int n_tr = trav.nvis(); // tree size
        assert(n_tr <= n_);

        // reach of u is distance of furhtest node from u
        for (int i = 0; i < n_tr; ++i) {
            V u = trav.visit(i);
            furthest_[u] = u;
            V p = trav.parent(u);
            W d_pu = trav.dist(u) - trav.dist(p);
            dist_[u] = i == 0 ? 0. : (dist_[p] + edge_metr(p, u, d_pu));
        }
        auto dist_par = [this,&trav](V u) -> WL {
            return dist_[u] - dist_[trav.parent(u)];
        };
        auto dist_furth = [this,&trav](V u) -> WL {
            return dist_[furthest_[u]] - dist_[u];
        };
        for (int i = n_tr-1; i >= 1; --i) {
            V u = trav.visit(i);
            V p = trav.parent(u);
            if (dist_par(u) + dist_furth(u) > dist_furth(p))
                furthest_[p] = furthest_[u];
            /* std::cout << u <<" "<< furthest_[u]
                      <<" "<< dist_[u] <<" "<< dist_furth(u)
                      << std::endl; */
        }
        //std::cout << "---\n";

        // skeleton is made of edges with sufficient reach
        auto long_reach = [this,&trav,alpha_p,alpha_q,
                           &dist_par,&dist_furth](V u) -> bool {
            // true if edge par(u)-->u has reach > alpha dist(par(u))
            /* Bad since we modify dist_:
               return (dist_par(u) + dist_furth(u)) * alpha_q
               > alpha_p * dist_[trav.parent(u)]; */
            return dist_[furthest_[u]] * alpha_q
            > (alpha_q + alpha_p) * dist_[trav.parent(u)]; 
        };
        // root is in the skeleton:
        V r = trav.visit(0);
        visit_[n_sk++] = r;
        // others :
        for (int i = 1; i < n_tr; ++i) {
            V u = trav.visit(i);
            if (long_reach(u)) {
                visit_[n_sk++] = u;
                WL d_prn = (dist_[u] + dist_furth(u)) // /(1 + alpha), i.e.:
                    * alpha_q / (alpha_p + alpha_q);
                dist_[u] = std::min(dist_[u], d_prn);
                nsons_[trav.parent(u)] += 1;
                /*std::cout << trav.parent(u) <<" "<< u <<" "<< dist_par(u)
                  << std::endl;*/
            }
        }
        // sort according to corrected distances :
        std::sort(visit_.begin() + 1, visit_.begin() + n_sk,
                  [this](int u, int v) {
                      // Important for width computation :
                      if (dist_[u] == dist_[v]) return nsons_[u] < nsons_[v];
                      return dist_[u] < dist_[v];
                  });
    }

    int width() {
        int wmax = 0;
        for (int i = 0, wdt = 1/* for root */; i < n_sk; ++i) {
            V u = visit_[i];
            wdt += nsons_[u] - 1;
            if (wdt > wmax) wmax = wdt;
        }
        return wmax;
    }

    double integrated_width(const T &trav) {
        double swdt = 0.;
        for (int i = 1; i < n_sk; ++i) {
            V u = visit_[i];
            V prt = trav.parent(u);
            if (dist_[u] == zero_weight) ; /* root */
            else if (dist_[prt] == zero_weight) swdt += std::log(dist_[u]);
            else swdt += std::log(dist_[u]) - std::log(dist_[prt]);
        }
        return swdt;
    }


};

#endif // SKELETON_HH

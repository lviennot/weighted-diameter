#ifndef TREEDEC_HH
#define TREEDEC_HH

#include <assert.h>
#include <stdint.h>
#include <climits>
#include <vector>

#include "mgraph.hh"
#include "dyn_graph.hh"
#include "traversal.hh"

template<typename G> // Dynamic symmetric graph class with G::vertex = int.

class treedec { // Heuristic for tree decompositon.

public:    
    typedef typename G::vertex V;

    const G &g;
    G chordal; // chordal completion
    std::vector<V> elim; // simplicial elimination ordering
    std::vector<int> elim_inv; // u -> pos of u in elim
    G clique_tree; // clique_tree as a graph
    std::vector<std::vector<V> > cliques;
    std::vector<int> parent; // parent of a clique in clique tree
    std::vector<std::vector<V> > separator_from_parent;
    std::vector<int> included_in;
    // u -> last clique containing it in elimination order.

    treedec(const G &g)
        : g(g), n(g.n()), cliques(g.n()), separator_from_parent(g.n()),
          included_in(g.n(),-1), parent(g.n()) {
        elim.reserve(n);
        completion_min_fill_in();
    }

    const std::vector<V> &elim_order() {
        return elim;
    }
    
private:

    const int n;

    struct q_node_ {
        V u;
        int level;
        int fill;
        int rnd;
        q_node_() {}
        q_node_(V u, int l, int f, int r) : u(u), level(l), fill(f), rnd(r) {}
    };

    static bool cmp_q_node(const q_node_ &qu, const q_node_ &qv) {
        if (qu.level != qv.level) return qu.level > qv.level;
        if (qu.fill != qv.fill) return qu.fill > qv.fill;
        return qu.u < qv.u;
    };

    typedef std::vector<q_node_> q_vec_;

    void completion_min_fill_in (int avg_deg_max=100 /* INT_MAX */,
                                 int max_bag_size=INT_MAX) {
        int n0 = g.n();
        elim.reserve(n0);
        if (max_bag_size == -1) max_bag_size = (int) sqrt(g.m());
        G h(g); 

        const int level_ceil = 8;
        // level of u is 1 + max lev of v for v visited neigh (initially 0),
        // smaller level nodes are visited first until ceil is reached
        
        std::vector<bool> vis(n0, false);
        std::vector<int> level(n0, 0);
        std::vector<int> neighb_vis(n0, 0);
        std::priority_queue<q_node_, std::vector<q_node_>,
                std::function<bool(const q_node_ &, const q_node_ &)> >
            queue(cmp_q_node);

        for (V u : g) {
            queue.push(q_node_(u, 0, h.fill_in(u) - h.degree(u), rand()));
        }
        std::cerr << "init fill in" << std::endl;
        std::cerr.flush();
        int64_t nq = 1, nvis = 0, madd = 0, degmax = 0, levelmoy = 0, m = h.m();
            
        while ( ! queue.empty() && (2*m / h.n() < avg_deg_max
                                    || h.n() > max_bag_size)) {
            q_node_ qu = queue.top();
            queue.pop();
            V u = qu.u;
            bool go_for_u = true;
            int l_now = std::min(level[u], level_ceil);
            int f_now = (qu.level == l_now)
                ? h.fill_in(u) - h.degree(u) : qu.fill;
            if (l_now > qu.level || (f_now >= 0 && f_now > qu.fill * 11 / 10)) {
                go_for_u = false;
                qu.level = l_now;
                qu.fill = f_now;
                queue.push(qu);
                ++nq;
            }
            if (go_for_u && ! vis[u]) {
                vis[u] = true;
                //rm std::cerr << nvis <<" "<< u <<" "<< l_now <<" "<< f_now <<"\n"; std::cerr.flush();
                ++nvis;
                elim.push_back(u);
                for (V v : h[u]) {
                    chordal.add_edge(u, v);
                    level[v] = std::max(level[v], level[u]+1);
                }

                int f = h.fill_in(u, true); // Neighborhood to clique:
                madd += f;
                int degu = h.degree(u);
                m += f - degu;
                h.del_vertex(u);
                //rm std::cout << nvis <<" "<< degu <<" "<< f << std::endl;

                if (degu > degmax) degmax = degu;
                levelmoy += level[u];
                if (nvis % 1000 == 0 || h.n() == 1) {
                    std::cerr << "td nvis=" << nvis
                              << " n=" << h.n() << " m=" << m
                              << " degmoy="<< (m / h.n()) 
                              << " maddd=" << madd 
                              << " lvl="<< level[u]
                              << " deg="<< degu
                              << " fill="<< f_now
                              << " nrq="<< nq
                              << "\r" ;
                    std::cerr.flush();
                }
            }
        }

        if ( ! queue.empty()) {
            std::cerr << "\n avg_deg limit, adding " << h.n() << std::endl;
            std::cerr.flush();
        }
        while ( ! queue.empty()) { queue.pop(); }
        for (V u : h) {
            queue.push(q_node_(u, 0, h.degree(u), rand()));
        }        
        while ( ! queue.empty()) {
            q_node_ qu = queue.top();
            queue.pop();
            V u = qu.u;
            vis[u] = true;
            ++nvis;
            elim.push_back(u);

            //for (V v : h) chordal.add_edge(u,v);
            /*
            int degu = 0;
            for (V v : h[u]) {
                chordal.add_edge(u, v);
                level[v] = std::max(level[v], level[u]+1);
                ++degu;
            }
            int f = h.fill_in(u, true);
            if (nvis % 1000 == 0 || h.n() == 1) {
                    std::cerr << "td nvis=" << nvis
                              << " n=" << h.n() << " m=" << m
                              << " degmoy="<< (m / h.n()) 
                              << " lvl="<< level[u]
                              << " deg="<< degu
                              << " fill="<< f
                              << "\r" ;
                    std::cerr.flush();
                }
            */
            h.del_vertex(u);
        }

        std::cerr << "\n maxdeg=" << degmax
                  << " levelmoy=" << (levelmoy / (nvis > 0 ? nvis : 1))
                  << " nvis=" << nvis << " ordsize=" << elim.size()
                  << " requeuings=" << nq << "\n";

    }

    
public:

    std::vector<V> mcs_clique_tree () {
 
        // take elimination ordering from MCS of chordal completion
        std::vector<V> ord(n), ord_inv(n);
        {
            traversal<G> trav(n);
            for (int r = n-1, ord_end = n-1, prev_vis = 0; r >= 0; --r) {
                if ( ! trav.visited(elim[r])) {
                    std::cerr << "mcs vis " << elim[r] <<" "<< r << std::endl; std::cerr.flush();
                    int nvis = trav.max_card_search(chordal, elim[r]);
                    std::cerr << "mcs vis 1" << std::endl; std::cerr.flush();
                    for (int i = 0, j = ord_end; i < nvis; ++i, --j) {
                        ord[j] = trav.visit(prev_vis + i);
                        ord_inv[ord[j]] = j;
                    }
                    prev_vis += nvis;
                    ord_end -= nvis;
                    std::cerr << "mcs vis 2" << std::endl; std::cerr.flush();
                }
            }
            for (int i = 0; i < n; ++i) {
                assert(ord[ord_inv[i]] == i);
                assert(ord_inv[ord[i]] == i);
                assert(0 <= ord[i] && ord[i] < n);
            }
        }

        std::cerr << "mcs 1" << std::endl; std::cerr.flush();
        
        // Since ord is a simplicial elimination ordering, [u U N(u)] is
        // a clique C_i where [u = ord[i]] and [N(u)] are neighbors of [u]
        // in [chordal] coming after [u] in ord. 
        {
            std::vector<bool> mark(n, false);
            std::vector<V> clq;
            //for (int i =0; i < n; ++i) parent[i] = i;
            for (int i = 0; i < n; ++i) {
                V u = ord[i];
                clq.push_back(u);
                
                int first = n, last = -1; // after in elim ordering
                for (V v : chordal[u]) {
                    int j = ord_inv[v];
                    if (i < j) {
                        mark[v] = true;
                        clq.push_back(v);
                        if (j < first) first = j;
                        if (j > last) last = j;
                    }
                }
                
                // If N(i) contains N(first), first is not a maximal clique
                if (first < n) {
                    bool included = true;
                    for (V x : chordal[ord[first]]) {
                        if (ord_inv[x] > first && ! mark[x]) included = false;
                    }
                    if (included) included_in[first] = i;
                }
                for (V v : chordal[u]) mark[v] = false;

                if (included_in[i] == -1) {
                    // [clq] is a maximal clique (because ord is a simplicial
                    // elimination ordering; checking first only is sufficient).
                    std::swap(cliques[i], clq);
                    included_in[i] = i;
                } else {  
                    clq.clear();
                    included_in[i] = included_in[included_in[i]]; // trans clos.
                }
                parent[included_in[i]] = first < n ? first : included_in[i];
            }
            for (int i =0; i < n; ++i) {
                assert(0 <= parent[i] && parent[i] < n);
                assert(0 <= included_in[i] && included_in[i] < n);
                if (included_in[i] == i) {
                    parent[i] = included_in[parent[i]];
                    clique_tree.add_edge(i, parent[i]);
                    // std::cerr << "clq " << i << " --> " << parent[i] << std::endl; std::cerr.flush();
                    mark[ord[i]] = true;
                    for (V v : chordal[ord[i]]) {
                        if (ord_inv[v] > i) mark[v] = true;
                    }
                    for (V v : chordal[ord[parent[i]]]) {
                        if (mark[v] && ord_inv[v] >= parent[i]) {
                            separator_from_parent[i].push_back(v);
                        }
                    }
                    V v = ord[parent[i]];
                    if (mark[v]) {
                        separator_from_parent[i].push_back(v);
                    }
                    for (V v : chordal[ord[i]]) mark[v] = false;
                    mark[ord[i]] = false;
                } else {
                    assert(cliques[i].size()==0);
                    clique_tree.add_edge(i, included_in[i]);
                    // std::cerr << "inc " << i << " --> " << included_in[i] << std::endl;
                    // Helps for finding a centroid bag.
                }
            }
        }

        return ord;

    } // mcs_clique_tree()
    
    std::vector<V> contraction_order (int nb_main_ldmks, int nb_per_sep=3) {
        std::vector<V> mcs = mcs_clique_tree(), ord;
        std::vector<V> elim_inv(n);
        for (int i = 0; i < n; ++i) elim_inv[elim[i]] = i;
        std::vector<bool> added(n, false);
        auto max_elim = [&elim_inv,&added] (std::vector<V> s) {
            assert(s.size() > 0);
            V m = s[0];
            for (V v : s) {
                if (( ! added[v]) && elim_inv[v] > elim_inv[m]) m = v;
            }
            return m;
        };
        auto add_to_ord = [&added,&ord] (V u) {
            if ( ! added[u]) {
                ord.push_back(u);
                added[u] = true;
            }
        };
        traversal<G> trav(n);
        std::priority_queue< std::pair<int, int> > queue;
        for (int i = n-1; i >= 0; --i) {
            if (parent[i] == i) {
                int sz = trav.tree_size(clique_tree, i);
                queue.push (std::make_pair(sz, i));

            }
        }
        for (int ldmks = 0; ! queue.empty() && ldmks < nb_main_ldmks; ++ldmks) {
            std::pair<int, int> sz_r = queue.top();
            queue.pop();
            int sz = sz_r.first, r = sz_r.second;
            int c = trav.centroid(clique_tree, r, sz/2);
            int cm = included_in[c];
            if (cm != r && r == included_in[r]) {
                //for (V v : separator_from_parent[cm]) {
                //    add_to_ord(v);
                //}
                if (separator_from_parent[cm].size() > 0) {
                    for (int i = 0; i < nb_per_sep; ++i)
                        add_to_ord(max_elim(separator_from_parent[cm]));
                }
                clique_tree.del_edge(cm, parent[cm]);
                queue.push(std::make_pair(trav.size(cm), cm));
                trav.tree_size_clear(clique_tree, r);
                sz = trav.tree_size(clique_tree, r);
                queue.push(std::make_pair(sz, r));
            } else {
                for (int chl : clique_tree[r]) {
                    if (chl != r) {
                        queue.push(std::make_pair(trav.size(chl), chl));
                    }
                }
                //for (V v : cliques[r]) {
                //    add_to_ord(v);
                //}
                for (int i = 0; i < nb_per_sep; ++i)
                    add_to_ord(max_elim(cliques[r]));
                clique_tree.del_vertex(r);
            }
        }

        // Add what's remaining while preserving elim order.
        for (int i = n-1; i >= 0; --i) {
            add_to_ord(elim[i]);
        }

        assert(ord.size() == n);
        
        return ord;
    }
    
}; // treedec

#endif // TREEDEC_HH

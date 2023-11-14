#ifndef DYN_GRAPH_HH
#define DYN_GRAPH_HH

#include <assert.h>
#include <vector>
#include <unordered_map>
#include <set>

#include "edge.hh"

/**
 * Dynamic graph implementation as hastable of maps.
 */

template<typename V, // vertex type (anything hashable)
         V not_vtx> // some value not used as vertex
         
class dyn_graph { // A (symmetric) dynamic unweighted graph.

private:
    std::unordered_map<V, std::set<V> > adj;
    
public:
    typedef V vertex;
    typedef edge::src_dst<V> edge;
    typedef int weight;
    static const V not_vertex = not_vtx;
    
    dyn_graph() {}

    dyn_graph(const dyn_graph &g) {
        for (V u : g)
            for (V v : g[u])
                add_edge(u, v);
    }

    bool add_edge(V u, V v, bool symmetric = true) {
        auto res = adj[u].insert(v);
        if (symmetric) adj[v].insert(u);
        return res.second; // true if edge was created (not present before)
    }

    void add_vertex (V u) { adj[u]; /* have u in table */ }

    size_t n() const { return adj.size(); }

    // Takes linear time :
    size_t m() const {
        size_t sum = 0;
        for (auto nu : adj)
            sum += nu.second.size();
        return sum;
    }

    size_t degree(V u) const {
        auto const nu = adj.find(u);
        return nu->second.size();
    }

    bool mem_edge (V u, V v) const {
        auto const nu = adj.find(u);
        return nu->second.find(v) != nu->second.end();
    }

    bool mem_vertex (V u) const { return adj.find(u) != adj.end(); }

    void del_edge (V u, V v, bool symmetric = true) {
        auto search = adj.find(u);
        if (search != adj.end()) search->second.erase(v);
        if (symmetric) {
            auto search = adj.find(v);
            if (search != adj.end()) search->second.erase(u);
        }
    }

    // Asserts that the graph is symmetric.
    void del_vertex (V u) {
        auto search = adj.find(u);
        if (search != adj.end()) {
            for (V v : search->second) {
                if (v != u) {
                    auto search = adj.find(v);
                    if (search != adj.end()) {
                        search->second.erase(u);
                    }
                }
            }
            adj.erase(u);
        }
    }

    std::vector<edge> all_edges() const {
        std::vector<edge> es;
        for (auto const &un : adj) {
            V u = un.first;
            for (auto v : un.second) {
                es.push_back (edge(u, v));
            }
        }
        return es;
    }

    int fill_in(V u, bool shortcut = false) {
        int fill = 0;
        for (auto e = adj[u].begin(), end = adj[u].end() ; e != end ; ) {
            V x = *e;
            ++e;
            if (x == u) continue;
            for (auto f = e ; f != end ; ++f) {
                V y = *f;
                if (y == u) continue;
                if ( ! mem_edge(x, y)) {
                    ++fill;
                    if (shortcut)
                        add_edge(x,y);
                }
            }
        }
        return fill;
    }

    
    //  --------------------- iterators : -----------------------
    //
    // for (V u : g)
    //     for (V v : g[u])
    //         f(u, v);
    //
    
    typedef
    typename std::unordered_map<V, std::set<V> >::const_iterator umap_it_;
    struct const_iterator : public umap_it_ {
        const_iterator(umap_it_ &&it) : umap_it_(it) {};
        V operator*() { return umap_it_::operator*().first; }
    };
    const_iterator begin() const { return const_iterator(adj.cbegin()); }
    const_iterator end() const { return const_iterator(adj.cend()); }

    const std::set<V> & operator[](V u) const {
        auto const nu = adj.find(u);
        return nu->second;
    }

    std::vector<V> vertices() {
        std::vector<V> vec;
        for (V u : *this) vec.push_back(u);
        std::sort(vec.begin(), vec.end());
        return vec;
    }

    
}; // dyn_graph

typedef dyn_graph<int, -1> i_dyn_graph;

#endif // DYN_GRAPH_HH

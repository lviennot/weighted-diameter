
// Compilation :
// without lex : g++ -std=c++11 -pthread -O3 src/graph_test.cc

#include <sys/time.h>
#include <stdint.h>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <thread>         // std::thread

#include "mgraph.hh"
#include "traversal.hh"
#include "dyn_graph.hh"
#include "treedec.hh"
#include "pruned_landmark_labeling.hh"
#include "skeleton.hh"

typedef mgraph<int> graph;

void bye (std::string msg) {
    std::cerr << msg << "\n";
    std::cerr.flush();
    exit(3);
}

double top (double t1, std::string msg) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double t2 = tv.tv_sec + tv.tv_usec * 1e-6;
    std::cerr << "-- time " << msg << " : " << (t2 - t1) << "s\n";
    std::cerr.flush();
    return t2;
}

/* g++ -std=c++11 -pthread -O3 src/lex.yy.c src/graph_test.cc
   extern FILE *yyin;
   extern int nlines;
   extern long long int yyint;
   int yylex();
*/

int main (int argc, char **argv) {
    double t = top (0., "start");

    // ------------------------- load ----------------------
    auto i_arg = [&argc,&argv](std::string a) {
        for (int i = 1; i < argc; ++i)
            if (a == argv[i])
                return i;
        return -1;
    };
    auto del_arg = [&argc,&argv,i_arg](std::string a) {
        int i = i_arg(a);
        if (i >= 0) {
            for (int j = i+1; j < argc; ++j)
                argv[j-1] = argv[j];
            --argc;
            return true;
        }
        return false;
    };
    auto get_arg = [&argc,&argv,i_arg](std::string a, std::string dft="") {
        int i = i_arg(a);
        if (i >= 0) {
            std::string b(argv[i+1]);
            for (int j = i+2; j < argc; ++j)
                argv[j-2] = argv[j];
            argc -= 2;
            return b;
        }
        return dft;
    };
    
    bool symmetrize = del_arg("-symmetrize");
    bool simplify = del_arg("-simple");
    bool output_order = del_arg("-output-order");
    std::string input_order = get_arg("-input-order");
    int src = 0, tgt = 100000;
    bool int_vertices = del_arg("-i");
    int nth = 1;
    if (argc > 2) nth = std::stoi (argv[2]) ;

    FILE *in;
    std::string dash = "-";
    if (argc == 1 || dash == argv[1]) { in = stdin; }
    else {
        in = fopen(argv[1], "r");
    }

    std::vector<graph::edge> edg;
    std::unordered_map<std::string,int> vi; // vertex index
    std::vector<std::string> lab;
    size_t n = 0;
    if (in && int_vertices) { // vertices are already ints from 0 to n-1
        long long int u, v, w;
        for ( ; fscanf(in, " %lld %lld %lld", &u, &v, &w) >= 3 ; ) {
            if (u >= n) n = u+1;
            if (v >= n) n = v+1;
            edg.push_back(graph::edge(u, v, w));
            if(symmetrize) edg.push_back(graph::edge(v, u, w));
        }
    } else if (in) { // vertices are any string
        char u[1024], v[1024];
        long long int w;
        for ( ; fscanf(in, " %s %s %lld", u, v, &w) >= 3 ; ) {
            if (vi[u] == 0) { lab.push_back(u); vi[u] = 1+n++; }
            if (vi[v] == 0) { lab.push_back(v); vi[v] = 1+n++; }
            edg.push_back(graph::edge(vi[u]-1, vi[v]-1, w));
            if(symmetrize) edg.push_back(graph::edge(vi[v]-1, vi[u]-1, w));
        }
    }
    size_t m = edg.size();
    std::cerr << "n=" << n << " m=" << m <<  std::endl;
    t = top (t, "load");

    // ------------------------- graph -----------------------
    graph g(edg);
    n = g.n(); m = g.m();
    std::cerr << "n=" << n << " m=" << m <<  std::endl;
    t = top (t, "graph");

    // ------------------------- simple -----------------------
    if (simplify) {
        g = g.simple();
        n = g.n(); m = g.m();
        std::cerr << "n=" << n << " m=" << m <<  std::endl;
        t = top (t, "simple");
    }

    // ------------------------ BFS ---------------------------
    {
        traversal<graph> trav(n);
        for (int i = 0; i < nth ; ++i) {
            int s = (src + i*12345) % n;
            trav.bfs(g, s);
            if (i <= 10) {
                int64_t dist = trav.dist(tgt % n);
                std::cerr << "bfs dist (" << i << "/" << nth << ") "
                          << s <<" -> "<< (tgt % n)
                          << " = " << dist <<  std::endl;
            }
            trav.clear();
        }
    }
    t = top (t, "bfs");
        
    /*
    // ----------------topological ordering ------------------------
    {
        traversal<graph> trav(n);
        try {
            std::vector<int> ord = trav.topological_ordering(g);
            for (int u : ord) std::cout << u << std::endl;
        } catch (std::string s) { std::cerr << "Cycle!\n"; }
    }
    std::cout.flush();
    t = top (t, "topological ordering");
        
    // ------------------------ dijkstra ---------------------------
    {
        traversal<graph> trav(n);
        trav.dijkstra(g, src % n);
        int64_t dist = trav.dist(tgt % n);
        std::cerr << "dist "<< (src % n) <<" -> "<< (tgt % n)
                  << " = " << dist <<  std::endl;
    }
    t = top (t, "dijkstra");
        
    // ------------------------ skels ---------------------------
    {
        int swdt = 0, mwdt = 0;
        double siw = 0., miw = 0.;
        const int batch = 1;

        auto dij_skel = [&g,n,tgt,&swdt,&mwdt,&siw,&miw](int src) {
            traversal<graph> trav(n);
            skeleton<traversal<graph> > sk(n);
            for (int i = 0; i < batch; ++i) {
                int s = src + i*i*7013;
                trav.dijkstra(g, s % n);
                int64_t dist = trav.dist(tgt % n);
                sk.of_traversal(trav, 1, 2, // 1/2
                                skeleton<traversal<graph> >::hop_count_metric);
                int width = sk.width();
                swdt += width ; mwdt = std::max(mwdt, width);
                double iwdt = sk.integrated_width(trav);
                siw += iwdt ; miw = std::max(miw, iwdt);            
                std::cerr << "dist "<< (s % n) <<" -> "<< (tgt % n)
                          << " = " << dist
                          << ", skel. width = " << width
                          << " " << iwdt
                          << ", skel. size = " << sk.size()
                          <<  std::endl;
                trav.clear();
                sk.clear();
            }
        };

        if (nth > 1) {
            std::cerr << "lauching " << nth << " threads " << std::endl;
            std::vector<std::thread> threads(nth);
            for (int i = 0; i < nth ; ++i)
                threads[i] = std::thread (dij_skel, src + i*10000);
            for (int i = 0 ; i < nth ; ++i) threads[i].join();
        } else {
            dij_skel(src);
        }

        std::cerr << "avg width=" << (swdt / (nth * batch))
                  << " max width=" << mwdt << std::endl; 
        std::cerr << "avg integr_width=" << (siw / (nth * batch))
                  << " max integer_width=" << miw << std::endl; 
    }
    t = top (t, "skels");
*/

    // ------------------------- dyn_graph -----------------------
    dyn_graph<int, -1> d;
    for (int u : g) {
        for (int v : g[u]) d.add_edge(u, v);
    }
    std::cerr << "n=" << d.n() << " m_sym=" << d.m() << std::endl;
    /*
      for (int u : d)
      for (int v : d[u])
      std::cerr << u <<" "<< v << std::endl;
      std::cerr << d.fill_in(3) << std::endl;
    */      
    t = top (t, "dyn_graph");

    // ------------------------- conn comps -----------------------
    {
        traversal<dyn_graph<int, -1> > trav(n);
        int s = 0;
        std::cerr << "cc sizes:" ;
        for (int u : d) {
            if ( ! trav.visited(u)) {
                int scc = trav.bfs(d, u);
                s += scc;
                std::cerr << " " << scc;
            }
        }
        std::cerr << ",   " << s << " in total" << std::endl;
    }
    t = top (t, "conn comps");

    // ------------------------- treedec -----------------------
    std::vector<int> ord; // = contraction_order(10, 4);
    if (input_order.compare("") != 0) {
        char u[1024];
        FILE *in = fopen(input_order.c_str(), "r");
        for ( ; fscanf(in, " %s", u) >= 1 ; ) {
            ord.push_back(vi[u]-1);
        }
        t = top (t, "read order"); 
        if (ord.size() != n)
            throw "number of elements in order do not correspond to graph size";
        std::vector<int> inv(n);
        for (int i=0; i<n; ++i) {
            assert(0 <= ord[i] && ord[i] < n);
            inv[ord[i]] = i;
        }
        bool perm = true;
        for (int j=0; j<n; ++j)
            if (ord[inv[j]] != j) perm = false;
        if ( ! perm) throw "not a permutation of vertices of the graph";
    } else {
        treedec<dyn_graph<int, -1> > td(d);
        t = top (t, "treedec"); 

        
        std::vector<int> mcs = td.mcs_clique_tree();
        t = top (t, "treedec mcs");
        for (int i =0; i < td.cliques.size() ; ++i) {
            if (td.cliques[i].size() > 0) {
                //std::cout << td.separator_from_parent[i].size() << std::endl;
                std::cout << td.cliques[i].size()
                          << " " << td.separator_from_parent[i].size()
                          << std::endl;
            }
        }
        //traversal<dyn_graph<int, -1> > trav(n);
        //std::cerr << trav.tree_size(td.clique_tree, 1000) << std::endl; 
        
        ord = td.elim_order();
        if (output_order) {
            for (int i = 0; i < n; ++i)
                std::cout << lab[ord[i]] << std::endl;
        }
    }
    t = top (t, "order");

    //exit(0);
    
    // ------------------------- pruned ll -----------------------
    std::vector<int> rev_ord(n);
    for (int i = 0, j = n-1; i < n; ++i, --j)
        rev_ord[i] = ord[j];
    pruned_landmark_labeling<mgraph<int> > pll(g, rev_ord);
    pll.print_stats(std::cerr);
    std::cerr << "dist "<< (src % n) <<" -> "<< (tgt % n)
              << " = " << pll.distance(src % n, tgt % n) <<  std::endl;
    t = top (t, "pruned ll");
    int64_t sd = 0;
    for (int i = 0; i < 1000 * 1000; ++i)
        sd = std::max(sd, pll.distance(rand() % n, rand() %n));
    std::cerr << sd << std::endl;
    t = top (t, "pruned ll 1M req dist");
    
    exit(0);
    
}

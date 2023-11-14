
/* Unweight a weighted graph by subdividing edges. */

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

    // ------------------------- arguments ----------------------
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
    bool subdivide = del_arg("-subdivide");
    float alpha = 0.4;
    
    FILE *in;
    std::string dash = "-";
    if (argc == 1 || dash == argv[1]) { in = stdin; }
    else {
        in = fopen(argv[1], "r");
    }
    
    // ------------------------- load ----------------------
    std::vector<graph::edge> edg;
    std::unordered_map<std::string,int> vi; // vertex index
    std::vector<std::string> lab;
    size_t n = 0;
    if (in) { // vertices are any string
        char u[1024], v[1024];
        long long int w;
        for ( ; fscanf(in, " %s %s %lld", u, v, &w) >= 3 ; ) {
            if (vi[u] == 0) { lab.push_back(u); vi[u] = 1+n++; }
            if (vi[v] == 0) { lab.push_back(v); vi[v] = 1+n++; }
            edg.push_back(graph::edge(vi[u]-1, vi[v]-1, w));
            if(symmetrize) edg.push_back(graph::edge(vi[v]-1, vi[u]-1, w));
        }
    } else { bye ("no graph ?"); }
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

    // ------------------------- conn comps -----------------------
    {
        traversal<graph> trav(n);
        int cc = 0, biggest=0, u_biggest=0;
        for (int u : g) {
            if ( ! trav.visited(u)) {
                int scc = trav.bfs(g,u);
                if (scc > biggest) {
                    biggest = scc;
                    u_biggest = u;
                }
                ++cc;
            }
        }
        std::cerr << cc << " connected components, biggest : "
                  << biggest << std::endl;
        t = top (t, "conn comps");

        // Restrict graph to biggest component :
        trav.clear();
        biggest = trav.bfs(g, u_biggest);
        n = 0;
        m = 0;
        edg.clear();
        vi.clear();
        std::vector<std::string> lab2;
        for (int i = 0; i < biggest; ++i) {
            int u = trav.visit(i);
            for (auto e : g[u]) {
                int v = e.dst;
                std::string su = lab[u], sv = lab[v];
                if (vi[su] == 0) { lab2.push_back(su); vi[su] = 1+n++; }
                if (vi[sv] == 0) { lab2.push_back(sv); vi[sv] = 1+n++; }
                edg.push_back(graph::edge(vi[su]-1, vi[sv]-1, e.wgt));
            }
        }
        graph g2(edg);
        g = g2;
        lab = lab2;
        n = g.n(); m = g.m();
        std::cerr << "Biggest comp : n=" << n << " m=" << m <<  std::endl;
    }
    t = top (t, "biggest comp");


    // ------------------- Weight stats -------------------------
    std::vector<int> weights;
    int weights_min = INT_MAX, weights_max = 0;
    long long int weights_sum = 0;
    int nw = 0;
    
    for (int u : g) {
        for (auto e : g[u]) {
            ++nw;
            weights.push_back(e.wgt);
            weights_sum += e.wgt;
            if (e.wgt < weights_min) weights_min = e.wgt;
            if (e.wgt > weights_max) weights_max = e.wgt;
        }
    }
    assert(nw == m);
    std::sort(weights.begin(), weights.end());
    float x = 0., y = 0.;
    std::cerr << "#_distr_weights ";
    for (int i = 0; i < m ; ++i) {
        float y2 = ((float)(i+1)) / ((float) m);
        float x2 = ((float)(weights[i] - weights_min))
            / ((float)(weights_max - weights_min));
        if (i == 0 || i == m-1 || i == m/2
            || x2 - x >= 0.1 || y2 - y >= 0.1) {
            std::cerr << weights[i] <<","<< (i+1) <<" ";
            x = x2; y = y2;
        }
    }
    std::cerr << std::endl;

    std::cerr << "nb weightsest (min/med/avg/max) : "
              << weights_min <<" / "<< weights[m/2]
              <<" / "<< (weights_sum / m) <<" / "<< weights_max
              << std::endl;    

    int one = (int)(weights_sum / m / 20);
    if (subdivide) one = weights_max;
    std::cerr << "unit = " << one << std::endl;
    t = top (t, "weight stats");


    // ----------------- Unweighted graph ------------------------
    for (int u : g) {
        for (auto e : g[u]) {
            int v = e.dst, w = e.wgt;
            if (u < v || ! symmetrize) {
                int len = w / one;
                ++m;
                if (len <= 1) std::cout << u <<" "<< v << std::endl;
                else {
                    int last = n++;
                    std::cout << u <<"\t"<< last << std::endl;
                    for (int i = 1; i < len; ++i) {
                        std::cout << last <<"\t"<< (last+1) << std::endl;
                        last = n++;
                        ++m;
                    }
                    std::cout << last <<"\t"<< v << std::endl;
                    ++m;
                }
            }
        }
    }
    std::cerr << "generated graph with n = "
              << n <<"  m = " << m << std::endl;
    t = top (t, "graph gen");
    
    exit(0);
    
}

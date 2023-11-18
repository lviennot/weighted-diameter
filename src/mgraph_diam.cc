
// Compilation :
// without lex : g++ -std=c++11 -pthread -O3 src/graph_test.cc

#include <sys/time.h>
#include <stdint.h>
#include <climits>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <algorithm>      // std::min
#include <thread>         // std::thread
#include <atomic>
#include <random>
#include <utility>        // std::pair

#include "verbose.hh"
#include "mgraph.hh"
#include "traversal.hh"
#include "eccentricity.hh"
#include "dyn_graph.hh"
#include "graphgen.hh"
//#include "treedec.hh"
//#include "pruned_landmark_labeling.hh"
//#include "skeleton.hh"

typedef mgraph<int> graph;


void usage_exit (char **argv) {
    auto paragraph = [](std::string s, int width=80) -> std::string {
        std::string acc;
        while (s.size() > 0) {
            int pos = s.size();
            if (pos > width) pos = s.rfind(' ', width);
            std::string line = s.substr(0, pos);
            acc += line + "\n";
            s = s.substr(pos);
        }
        return acc;
    };
    
    std::cerr <<"\nUsage: "
              << argv[0] <<" [options] [filename]\n\n"

              << paragraph ("Read or generate a graph and compute some of its "
                            "parameters such as radius, diameter, or all "
                            "eccentricities. The graph is considered directed "
                            "(use -symmetrize to obtain a symmetric graph "
                            "equivalent to an undirected graph) and can be"
                            " weighted (use option -weighted). "
                            "The graph is either "
                            "generated from options (see, e.g. -grid), or "
                            "read from a file (use filename - for reading from "
                            "the standard input). The input format consists of "
                            "a sequence of [u v] pairs (for arc u -> v) for "
                            "unweighted graphs or [u v wgt] triples for "
                            "weighted graphs (space and newline characters "
                            "are ignored). Nodes can be identified by any "
                            "string (excluding whitespaces).\n")

              << paragraph ("Distances: this program computes various "
                            "paramters related to distances in the input graph."
                            " The distance d(u,v) from u to v is defined as the"
                            " length of a shortest path from u to v "
                            "(non-negative weights are assumed). The length "
                            "of a path is its number of arcs if the graph is "
                            "unweighted, and the sum of weigts of its arcs "
                            "if it is weighted. It is "
                            "considered infinite if v is not reachable from u."
                            " Note that d(u,v) can be different from d(v,u) in "
                            "a directed graph.\n")

              << paragraph ("Efficiency: nodes are numbered as ints as they "
                            "are encountered, this limits to graphs of 2^31 "
                            "nodes at most (the code can easily be patched to "
                            "use larger integers). Some computations use "
                            "multithreading (see -n-thread).\n")

              << paragraph ("Possible options:\n")

              << paragraph ("  -eccentricity-all  Compute all eccentricities "
                            "using multithreading. "
                            "The eccentricity e(u) of a node u is max_v d(u,v) "
                            "where d(u,v).\n")
        
              << paragraph ("  -closeness-all  Compute the closeness "
                            "centrality of all nodes using multiple threads."
                            " Outputs for each node u"
                            " a quadruple [u r s h] where u is the label of u "
                            "(as given in input), r is the size |R(u)| where "
                            "R(u) denotes the set of nodes reachable from u,"
                            " s is the sum of distances s = sum_{v in R(u)} "
                            "d(u,v), and h is the harmonic sum h = "
                            "sum_{v in R(u)} 1 / d(u,v). Closeness centrality"
                            " is classically defined as 1/s, a normalized "
                            "value can be obtained with (r-1)^2 / (s(n-1)). "
                            "Use -columns-verb 0 to output only those lines.\n")
        
              << paragraph ("  -params-verb v   Computed parameters are "
                            "output to stdout on one line, separated by "
                            "spaces. The first four fields are [n m dir wgt]:"
                            " the number n of nodes, the number m of edges, "
                            "dir=1 for a directed graph (0 if -symmetrize is "
                            "used) and wgt=1 for a weighted graph (use "
                            "-weighted). Use -diameter, -radius for diameter, "
                            "and/or radius. Nothing is output with v=0, the "
                            "default is 1, while more information about the "
                            "computation is output with 2.\n")

              << paragraph ("  -diameter       Compute the diameter, that is "
                            "the maximum eccentricity of a node (see "
                            "-eccentricity-all). Consider -largest-scc for non "
                            "(strongly) connected graphs.\n")
        
              << paragraph ("  -grid l         Generate a l x l grid.\n")
        
              << paragraph ("  -n-thread       Specify the number of threads "
                            "to use. Options using multithreading are "
                            "-eccentricity-all, -closeness-all.\n")
        
              << paragraph ("  -power-law b    Generate a random graph "
                            "according to the configuration model such that "
                            "the degree sequence follows a power law with "
                            "parameter b.\n")
        
              << paragraph ("  -radius         Compute the radius, that is "
                            "the minimum eccentricity of a node (see "
                            "-eccentricity-all). Consider -largest-scc for non "
                            "(strongly) connected graphs.\n")
        
              << paragraph ("  -reverse        Reverse all arcs.\n")
        
              << paragraph ("  -largest-scc    Restrict the graph to its "
                            "largest strongly connected component.\n")

              << paragraph ("  -simple         Ensure that the graph is simple "
                            "by removing duplicate arcs: among "
                            
                            "all u -> v arcs, keep one with minimum weight.\n"
                            )
        
              << paragraph ("  -symmetrize     Add arc v -> u for each arc "
                            "u -> v in the original graph.\n")

              << paragraph ("  -verbosity v    Print on stderr various "
                            "information depending on v=0,1,2 (defaults to "
                            "0 for nothing).\n")

              << paragraph ("  -weighted       Read a weighted graph: the "
                            "sequence read is interpreted as a list of triples "
                            "[u v wgt] instead of pairs [u v].\n")

              <<"";

    exit(1);
}



void bye (std::string msg) {
    std::cerr << msg << "\n";
    std::cerr.flush();
    verb::end();
    exit(3);
}


void printerr_distribution(std::vector<int64_t> &v, std::string what="") {
    int n = v.size();
    std::sort(v.begin(), v.end());
    int64_t v_min = INT_MAX, v_max = 0, a_prev = -1;
    int nvals = 0;
    for (int64_t a : v) {
        if (a < v_min) v_min = a;
        if (a > v_max) v_max = a;
        if (nvals == 0 || a != a_prev) { ++nvals; a_prev = a; } 
    }

    double x = 0., y = 0.;
    verb::cerr() << "#_distr_" << what <<" ";
    for (int i = 0; i < n ; ++i) {
        double y2 = ((double)(i+1)) / ((double) n);
        double x2 = ((double)(v[i] - v_min))
            / ((double)(v_max - v_min));
        int i_back = i - (int)(0.03 * n);
        if (i < 3 || i >= n-3 || i == n/2
            || x2 - x >= 0.1 || y2 - y >= 0.1 
            || (v[i] != v[i-1] && (nvals <= 15
                                   || (i_back >= 0 && v[i_back] == v[i-1])) )) {
            verb::cerr() << v[i] <<","<< (i+1) <<" ";
            x = x2; y = y2;
        }
    }
    verb::cerr() << std::endl;
}

int64_t mean1000(std::vector<int64_t> &v) {
    int64_t sum = 0;
    for (int64_t a : v) sum += a;
    int n = v.size();
    if (n == 0) return -1;
    return (sum * 1000) / n;
}

class vector_int : public std::vector<int> {
 public:
    vector_int() : std::vector<int>(2) { }
};



int main (int argc, char **argv) {

    // ------------------------- arguments ----------------------
    auto i_arg = [&argc,&argv](std::string a) {
        //for (int i = 1; i < argc; ++i) std::cerr << argv[i] <<" "; std::cerr <<" for "<< a <<"\n";
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
            if (i+1 > argc) bye(a + " requires an argument, try -h");
            std::string b(argv[i+1]);
            for (int j = i+2; j < argc; ++j)
                argv[j-2] = argv[j];
            argc -= 2;
            //std::cerr << a <<" "<< b <<"\n";
            return b;
        }
        return dft;
    };
    auto get_iarg = [&get_arg](std::string a, int dft) {
        std::string s = get_arg(a, "");
        if (s != "") return std::stoi(s);
        else return dft;
    };
    auto get_darg = [&get_arg](std::string a, double dft) {
        std::string s = get_arg(a);
        if (s != "") return std::stod(s);
        else return dft;
    };
    auto no_more_options = [&argc,&argv]() {
        for (int i = 1; i < argc; ++i) {
            std::string b(argv[i]);
            if(b.size() > 1 && b[0] == '-')
                bye(b + ": unrecognize option ?????, try -h");
        }
    };
    
    bool do_vector = del_arg("-test-vector"); //notcom
    bool symmetrize = del_arg("-symmetrize");
    bool del_zero_edges = del_arg("-del-zero-edges"); //tocom
    bool do_reverse = del_arg("-reverse");
    bool do_scc = del_arg("-largest-scc");
    int cols_verb = get_iarg("-params-verb", 1);
    if (del_arg("-no-cols") || del_arg("-no-params")) cols_verb = 0;
    bool directed = ! symmetrize;
    bool simple = del_arg("-simple");
    bool weighted = del_arg("-weighted");
    //bool biggest_comp = del_arg("-bcc");
    bool do_all_ecc = del_arg("-eccentricity-all") || del_arg("-all-ecc");
    int many_bfs = get_iarg("-many-bfs", 0); //notcom
    int n_big_graph = get_iarg("-n-big-graph", 0); //notcom
    bool do_skel = del_arg("-skel"); //tocom
    bool do_rad = del_arg("-rad") || del_arg("-radius");
    int sample_size = get_iarg("-sample-size", 20); //notcom
    bool do_diam = del_arg("-diam") || del_arg("-diameter");
    bool do_diam_all = del_arg("-diam-all"); //notcom
    bool do_closeness = del_arg("-closeness"); //notcom
    int source_node = get_iarg("-source", 0);
    int loop_limit = get_iarg("-loop-limit", -1); //notcom
    bool do_closeness_all = del_arg("-closeness-all");
    bool all_bfs = del_arg("-all-bfs");
    bool very_low_cert = del_arg("-very-low-cert"); //very low mem usag for cert
    std::atomic<bool> low_cert(del_arg("-low-cert") || very_low_cert);
    std::string optim_certif = get_arg("-optim-certif");
    std::string algo = get_arg("-algo");
    int n_thread = get_iarg("-n-thread", std::thread::hardware_concurrency());
    int beta_hyp = get_iarg("-beta-hyp", INT_MAX);
    bool do_quad_antipode = del_arg("-quad-antipode");
    bool do_antipodes = del_arg("-antipodes");
    bool do_coballs = del_arg("-coballs");
    //double alpha = 0.4;
    int verbosity = get_iarg("-verbosity", 0);
    verbosity = get_iarg("-verb", verbosity);
    int mem_limit_mb = get_iarg("-mem-limit-mb", 48000);
    int rand_edge_del = get_iarg("-edge-del", 0); // percentage
    int rand_edge_weight = get_iarg("-rand-weight", 0); // percentage
    bool rand_orient = del_arg("-rand-orient");
    int losange = get_iarg("-losange", 0);
    int bow_tie = get_iarg("-bow-tie", 0);
    int cycle = get_iarg("-cycle", 0);
    int path = get_iarg("-path", 0);
    int grid = get_iarg("-grid", 0);
    int n_gen = get_iarg("-n-gen", 0);
    double udg_deg = get_darg("-udg-deg", 0);
    double power_law = get_darg("-power-law", 0);
    bool generate = losange > 0 || bow_tie > 0 || udg_deg > 0
        || cycle > 0 || grid > 0 || path > 0 || power_law > 0.;
    bool print_graph = del_arg("-print-graph");
    
    if (del_arg("-h") || del_arg("-help") || del_arg("--help")) {
        usage_exit(argv);
    }
    no_more_options();

    if (argc != generate ? 0 : 1) usage_exit(argv);

    
    verb::begin(verbosity);

    if (n_thread > 1 && n_thread != std::thread::hardware_concurrency()) {
        verb::cerr() <<"System recommendation for -n-thread: "
                     << std::thread::hardware_concurrency() <<"\n";
    }
    
    if (do_vector) {
        verb::lap("test empty");
        std::vector<vector_int> t(10*1000*1000);
        for (int i = 0; i < t.size(); ++i) {
            for (int j = 1; j <= 0; ++j) t[i].push_back(j);
        }
        verb::cerr() << "capcity: "<< t[0].capacity() <<"\n";
        verb::lap("test vect of vect");
    }
    
    // ------------- Printing results in space-separated columns ------------
    std::vector<std::string> col_name;
    std::unordered_map<std::string, int> col_i;
    std::vector<int64_t> col_val;
    auto set_col = [&col_name, &col_i, &col_val, cols_verb](std::string name,
                                                 int64_t val, int verb=2) {
        verb::cerr(2) << name <<" "<< val
                      <<" verb="<< verb <<" params-verb="<< cols_verb <<"\n";
        if (verb > cols_verb) return;
        int i = col_i[name];
        if (i == 0) {
            col_name.push_back(name);
            i = col_name.size();
            col_i[name] = i;
        }
        while (col_val.size() < i) col_val.push_back(0);
        col_val[i-1] = val;
    };
    auto set_distr_col = [&set_col, cols_verb](std::vector<int64_t> &v,
                                               std::string name,
                                    int verb=2) {
        printerr_distribution(v, name);
        if (verb > cols_verb) return;
        int n = v.size();
        set_col(name + "_min", n > 0 ? v[0] : -1, verb);
        set_col(name + "_n/4", n > 0 ? v[n/4] : -1, verb);
        set_col(name + "_n/2", n > 0 ? v[n/2] : -1, verb);
        set_col(name + "_3n/4", n > 0 ? v[3*n/4] : -1, verb);
        set_col(name + "_90%", n > 0 ? v[9*n/10] : -1, verb);
        set_col(name + "_99%", n > 0 ? v[99*n/100] : -1, verb);
        set_col(name + "_max", n > 0 ? v[n-1] : -1, verb);
        set_col(name + "_avg*1000", mean1000(v), verb);
    };
    auto print_cols = [&col_name, &col_val](std::string more="") {
        std::cout << "#" ;
        for (int i = 0; i < col_name.size() ; ++i) {
            std::cout <<","<< col_name[i] <<":"<< col_val[i];
        }
        std::cout << std::endl;
        
        for (int64_t v : col_val) std::cout << v <<" ";

        std::cout << more << std::endl << "#_";
        int i = 0;
        for (std::string c : col_name) {
            std::cout << ++i <<"_"<< c <<" ";
        }
        std::cout << std::endl;
    };

    
    /*
    std::cout << "#_n m rad rad_nbfs rad_cert diam diam_nbfs diam_cert" // 8
              <<" nallecc lb_cert ub_cert last_lb" // 4
              <<" periph center furthest tree_ctr" // 4
              <<" far_nodes tree_ctr_far" // 2
              <<" nb_furth_min nb_furth_med nb_furth_avg" // 3
              <<" nb_furth_3/4 nb_furth_9/10 nb_furth_99/100 nb_furth_max" // 4
              <<" pruned_min pruned_med pruned_avg" // 3
              <<" pruned_3/4 pruned_9/10 pruned_99/100 pruned_max" // 4
              <<" furth_ctr ctr2_size furth_ctr2"  // 3
              // 35
              <<" all_ecc_lb_cert rad_gdy_cert all_ecc_ub_cert diam_gdy_cert"
              // 39
        //<< " nb_corners alpha"
              << "\n";
    */
    
    
    // ------------------------- load ----------------------
    FILE *in = nullptr;
    std::string dash = "-";
    if ((argc == 1 || dash == argv[1]) && ! generate) {
        verb::cerr() << "open stdin" << std::endl;
        in = stdin;
    } else if (argc >= 1) {
        verb::cerr() << "open " << argv[1] << std::endl;
        in = fopen(argv[1], "r");
    }
    
    std::vector<graph::edge> edg;
    std::unordered_map<std::string,int> vi; // vertex index
    std::vector<std::string> lab;
    int n = 0;
    if (in) { // vertices are any string
        char u[1024], v[1024];
        if (weighted) {
            int64_t w;
            for ( ; fscanf(in, " %s %s %lld \n", u, v, &w) >= 3 ; ) {
                assert(w >= 0 && w <= INT_MAX);
                if (vi[u] == 0) { lab.push_back(u); vi[u] = ++n; }
                if (vi[v] == 0) { lab.push_back(v); vi[v] = ++n; }
                //w /= 1000;
                edg.push_back(graph::edge(vi[u]-1, vi[v]-1, w));
                if(symmetrize) edg.push_back(graph::edge(vi[v]-1, vi[u]-1, w));
            }
        } else {
            for ( ; fscanf(in, " %s %s \n", u, v) >= 2 ; ) {
                if (vi[u] == 0) { lab.push_back(u); vi[u] = ++n; }
                if (vi[v] == 0) { lab.push_back(v); vi[v] = ++n; }
                int w = 1; //1000000 + (rand() % 1000);
                edg.push_back(graph::edge(vi[u]-1, vi[v]-1, w));
                if(symmetrize) edg.push_back(graph::edge(vi[v]-1, vi[u]-1, w));
            }
        }
        fclose(in);
    } else if (generate) {
        verb::cerr() << "--- losange: " << losange <<" bow_tie: " << bow_tie
                  <<" cycle: "<< cycle
                  <<" grid: "<< grid <<" path: "<< path
                  <<" udg_deg: "<< udg_deg
                  <<" power_law: "<< power_law
                  <<"\n";
        if (udg_deg > 0) {
            edg = unit_disk_graph(n_gen, udg_deg);
            for (int u = 0; u < n_gen; ++u) {
                std::string su = std::to_string(u);
                assert(vi[su] == 0);
                lab.push_back(su);
                assert(u == n);
                vi[su] = ++n;
            }
        } else {
            dyn_graph<int, -1> g;
            if (losange > 0) g = losange_graph(losange);
            else if (bow_tie > 0) g = bow_tie_graph(bow_tie, bow_tie);
            else if (cycle > 0) g = cycle_graph(cycle);
            else if (grid > 0) g = grid_graph(grid, grid);
            else if (path > 0) g = path_graph(path);
            else if (power_law > 0.) g = power_law_random_graph(n_gen, power_law);
            else assert(false);
            for (int u = 0; true; ++u) {
                if ( ! g.mem_vertex(u)) break;
                std::string su = std::to_string(u);
                if (vi[su] == 0) { lab.push_back(su); vi[su] = ++n; }
            }
            for (int u : g.vertices()) {
                std::string su = std::to_string(u);
                if (vi[su] == 0) { lab.push_back(su); vi[su] = ++n; }
                //std::cerr << "g["<< u <<"] =";;
                //for (int v : g[u]) std::cerr <<" "<< v;
                //std::cerr << "\n";
            }
            for (int u : g) {
                std::string su = std::to_string(u);
                for (int v : g[u]) {
                    std::string sv = std::to_string(v);
                    assert (vi[sv] >  0);
                    edg.push_back(graph::edge(vi[su]-1, vi[sv]-1, 1));
                }
            }
        }
        directed = false;
    } else { bye ("no graph ?"); }
    size_t m = edg.size();
    verb::cerr() << "n=" << n << " m=" << m
              << " weighted=" << weighted << " symmetrize=" << symmetrize
              <<  std::endl;
    verb::lap("load");

    

    if (n <= 0) bye("Empty graph ???\n  (Check input format : two space-separated ints per line for an unweighted graph and 3 space-separated ints per line for a weighted graph.)");

    // ------------------------- graph -----------------------
    if (del_zero_edges) {
        // graph of zero weight edges :
        std::vector<graph::edge> edg_zero;
        for (graph::edge e : edg) {
            if (e.wgt == 0) {
                edg_zero.push_back(e);
                edg_zero.push_back(graph::edge(e.dst, e.src, e.wgt));
            }
        }
        graph gz(n, edg_zero);
        // connected components :
        traversal<graph> trav(n);
        trav.strongly_connected_components(gz);
        // Contract components :
        std::vector<graph::edge> edg_contr;
        for (graph::edge e : edg) {
            int c_src = trav.scc_number(e.src), c_dst = trav.scc_number(e.dst);
            int src = trav.scc_node(c_src), dst = trav.scc_node(c_dst);
            if (e.wgt != 0) edg_contr.push_back(graph::edge(src, dst, e.wgt));
        }
        edg = edg_contr;
    }
    
    if (rand_orient) {
        assert( ! directed);
        std::vector<graph::edge> edg_sel;
        for (graph::edge e : edg) {
            if (e.src < e.dst) {
                if (rand() >= RAND_MAX / 2) {
                    edg_sel.push_back(e);
                } else {
                    edg_sel.push_back(graph::edge(e.dst, e.src, e.wgt));
                }
            }
        }
        edg = edg_sel;
        directed = true;
    }
    
    if (rand_edge_del > 0 || rand_edge_weight > 0) {
        if (rand_edge_weight > 0) weighted = true;
        std::vector<graph::edge> edg_sel;
        for (graph::edge e : edg) {
            if ( (directed || e.src < e.dst) && rand() % 100 >= rand_edge_del) {
                if (rand_edge_weight > 0) {
                    e.wgt = 0 + (rand() % rand_edge_weight);
                }
                edg_sel.push_back(e);
                if ( ! directed)
                    edg_sel.push_back(graph::edge(e.dst, e.src, e.wgt));
            }
        }
        edg = edg_sel;
    }
    
    graph g(edg);
    n = g.n(); m = g.m();
    verb::cerr() << "n=" << n << " m=" << m <<  std::endl;
    verb::lap("graph");

    // ------------------------- simple -----------------------
    if (simple) {
        g = g.simple();
        n = g.n(); m = g.m();
        verb::cerr() << "n=" << n << " m=" << m <<  std::endl;
        verb::lap("simple");
    }

    // ------------------------- strong conn comps -----------------------
    int g_scc_nb = 0, u_biggest = 0;
    {
        traversal<graph> trav(n);
        g_scc_nb = trav.strongly_connected_components(g);
        int biggest = trav.scc_largest();
        u_biggest = trav.scc_node(biggest);
        verb::cerr() << g_scc_nb << " strongly connected component(s), biggest : "
                  << trav.scc_size(biggest) <<" (that of node"
                  << lab[u_biggest] <<"=lab["<< u_biggest <<"])\n";
        verb::lap("strong conn comps");

        if (do_scc) {
            // check scc comput is right (at least in largest comp.)
            graph h = g.reverse();
            traversal<graph> tg(n), th(n);
            tg.bfs(g, u_biggest); th.bfs(h, u_biggest);
            int n_biggest = 0;
            for (int u : g) if (tg.visited(u) && th.visited(u)) ++n_biggest;
            verb::cerr() << " check : " << n_biggest << " nodes in biggest\n";
            assert(trav.scc_size(biggest) == n_biggest);
            set_col("scc_of", u_biggest);
        
            // Restrict graph to biggest component :
            if (g_scc_nb > 1) {
                n = 0;
                m = 0;
                edg.clear();
                vi.clear();
                std::vector<std::string> lab2;
                for (int u : g) {
                    if (trav.scc_number(u) == biggest) {
                        const std::string &su = lab[u];
                        for (auto e : g[u]) {
                            int v = e.dst;
                            if (vi[su] == 0)
                                { lab2.push_back(su); vi[su] = ++n; }
                            if (trav.scc_number(v) == biggest) {
                                const std::string &sv = lab[v];
                                if (vi[sv] == 0)
                                    { lab2.push_back(sv); vi[sv] = ++n; }
                                edg.push_back(graph::edge(vi[su]-1,vi[sv]-1,e.wgt));
                            }
                        }
                    }
                }
                graph g2(n, edg);
                g = g2;
                lab = lab2;
                n = g.n(); m = g.m();
                trav.clear();
                g_scc_nb = trav.strongly_connected_components(g);
                verb::cerr() << "Graph scc : n=" << n << " m=" << m <<"\n";
                verb::lap("restrict to scc");
            }
            
        }
    }

    n_thread = std::min(n_thread, n);

    // ----------------------- reverse graph -----------------

    graph g_rev = g.reverse();
    {
        // check reverse
        for (int u : g) {
            for (auto e : g[u]) {
                assert(g_rev.has_edge(e.dst, u));
                assert(g_rev.has_edge(e.dst, u)
                       && e.wgt == g_rev.edge_weight(e.dst, u));
            }
        }
        if ( ! directed) { // check symmetric
            g = g_rev.reverse();
            for (int u : g_rev) {
                for (auto e : g_rev[u]) {
                    if ( ! (g.has_edge(e.dst, u) && g.has_edge(u, e.dst)
                            && e.wgt == g.edge_weight(e.dst, u)
                            && e.wgt == g.edge_weight(u, e.dst)) )
                        throw std::invalid_argument("graph is not symmetric");
                }
            }
        }
        verb::cerr() <<"reverse graph : n=" << g_rev.n()
                  <<" m="<< g_rev.m() <<"\n";
        verb::lap("reverse");
    }
    
    if(do_reverse) {
        std::swap(g, g_rev);
        verb::cerr() << "reverse graph!\n";
        //for (int u : g_rev) for (int v : g_rev[u]) std::cerr<< u<<" "<< v<<"\n";
    }


    // ------------------------- print graph -----------------------
    if (print_graph) {
        for (int u : g) {
            for (auto e : g[u]) {
                if ( ! weighted)
                    std::cout << lab[u] <<"\t"<< lab[e.dst] <<"\n";
                else
                    std::cout << lab[u] <<"\t"<< lab[e.dst]
                              <<"\t"<< e.wgt <<"\n";
            }
        }
        verb::lap("print-graph");
    }


    // ----------------------- Closeness -------------------------
    {
        traversal<graph> trav(n);

        auto print_sums = [&g,&trav,&lab,weighted](int source_node) {
            trav.clear();
            int nvis = 0;
            if (weighted) {
                nvis = trav.dijkstra(g, source_node);
            } else {
                nvis = trav.bfs(g, source_node);
            }
            int64_t sum = 0;
            double harm = 0.;
            for (int i = 0; i < nvis; ++i) {
                int v = trav.visit(i);
                int64_t d = trav.dist(v);
                sum += d;
                if (d > 0) harm += 1./double(d);
            }
            std::cout << lab[source_node] <<" "<< nvis
                     <<" "<< sum <<" "<< harm <<"\n";
        };

        if (do_closeness) {
            std::cout << "# v |Reach(v)| dist_sum harm_sum\n";
            print_sums(source_node);
            verb::lap("closeness");
        }

        if (do_closeness_all && n_thread <= 1) {
            
            std::cout << "# v |Reach(v)| dist_sum harm_sum\n";
            for (int s = 0; s < n && (loop_limit < 0 || s < loop_limit); ++s) {
                print_sums(s);
            }
            verb::lap("closeness-all");
            
        } else if (do_closeness_all) {
            
            std::cout << "# v |Reach(v)| dist_sum harm_sum\n";
            
            std::mutex mutex;

            auto go_clo = [n,loop_limit,n_thread,
                           &g,&lab,weighted,&mutex](int i_thd) {
                traversal<graph> trav(n);
                for (int s = 0;
                     s < n && (loop_limit < 0 || s < loop_limit); ++s) {
                    if ((1237L * s) % n_thread == i_thd) {
                        // print_sums(s); :
                        int source_node = s;
                        trav.clear();
                        int nvis = 0;
                        if (weighted) {
                            nvis = trav.dijkstra(g, source_node);
                        } else {
                            nvis = trav.bfs(g, source_node);
                        }
                        int64_t sum = 0;
                        double harm = 0.;
                        for (int i = 0; i < nvis; ++i) {
                            int v = trav.visit(i);
                            int64_t d = trav.dist(v);
                            sum += d;
                            if (d > 0) harm += 1./double(d);
                        }
                        mutex.lock();
                        std::cout << lab[source_node] <<" "<< nvis
                                  <<" "<< sum <<" "<< harm <<"\n";
                        mutex.unlock();
                        if (i_thd == 0 && verb::progress())
                            verb::cerr("closeness") << s <<" "
                                <<"todo: "<< (n - s) <<" / "<< n <<"\n";
                    }
                }
            };

            std::vector<std::thread> threads(n_thread);
            for (int i = 0; i < n_thread; ++i) {
                threads[i] = std::thread(go_clo, i);
            }
            verb::cerr() << n_thread <<" closeness threads\n";
            for (int i = 0; i < n_thread; ++i) threads[i].join();
            
            verb::lap("closeness-all");
            
        }
    }
    

    // ----------------------- Work on eccentricities -----------
    
    eccentricity<graph> ecc(g, g_rev, directed, weighted);
    verb::lap("setup ecc");

    // ----------------------- Pruned Dijkstra for skeleton -----------------

    if (do_skel) {

        traversal<graph> trav(n);
        
        /*
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> unif(0.0, 1.0);

        auto edge_rnd = [&gen, &unif](int len){ // min rand val of subdived edge
            int j = 0;
            float p = unif(gen);
            while (j <= len) {
                // dist to next rnd val < p follows geom distr:
                float dj = std::log(1. - unif(gen)) / std::log(1. - p);
                if (dj + j >= INT_MAX) break;
                j += (int) dj;
                p = p * unif(gen);
            }
            return std::pair<float, int>(p, j);
        };

        std::vector<float> edg_rnd(g.m());
        std::vector<int> edg_rnd_len(g.m());
        
        for (int e=0; e<g.m(); ++e) edg_rnd[e] = 1.;
        for (int u=0; u<n; ++u) {
            int iv = 0, degsum = g.degree_sum(u);
            for (auto e : g[u]) {
                auto r = edge_rnd(e.wgt);
                edg_rnd[degsum + iv] = r.first;
                edg_rnd_len[degsum + iv] = r.second;
                if(u % 300 == 0) {
                    std::cerr <<"  rnd "<< u <<" -> "<< e.dst
                              <<" : "<< edg_rnd[degsum + iv] <<"\n";
                }
                ++iv;
            }
        }
        */

        std::vector<int> edg_rnd(g.m());
        for (size_t e = g.m()-1; e != 0; --e) edg_rnd[e] = rand();

        std::vector<graph::edge> edg_perm;
        size_t e = 0;
        for (int u : g) {
            for (int v : g[u]) edg_perm.push_back(graph::edge(u,v,edg_rnd[e]));
            ++e;
        }
        std::sort(edg_perm.begin(), edg_perm.end(),
                  [](graph::edge e, graph::edge f) {
                      if (e.wgt != f.wgt) return e.wgt < f.wgt;
                      if (e.dst != f.dst) return e.dst < f.dst;
                      return e.src < f.src;
                  });

        verb::lap("rnd val for edges");
        
        std::vector<int64_t> rnd_stop_dist(n);

        int64_t nv = 0;
        for (size_t e = 0; e < g.m(); ++e) {

            int e_s = edg_perm[e].src, e_d = edg_perm[e].dst;
            int e_rnd = edg_perm[e].wgt;
            
            for (int v=0; v<n; ++v) {
                rnd_stop_dist[v] = INT64_MAX;
            }
            
            auto filter = [e_s, e_d, &trav, &g, &edg_rnd, //&edg_rnd_len,
                           &e_rnd, &rnd_stop_dist]
                    (int v, int64_t d, int par, int64_t dpar, int iv) {
                //d = dpar + edg_rnd_len[g.degree_sum(par) + iv];
                if (edg_rnd[g.degree_sum(par) + iv] < e_rnd) {
                    rnd_stop_dist[v] = std::min(rnd_stop_dist[par],
                                            d + /* alpha / (1 - alpha) */ d);
                }
                return d <= rnd_stop_dist[par] // prune
                       && (par != e_s || v == e_d); // dijkstra from edge e
            };

            trav.clear();
            trav.dijkstra_i(g, e_s, filter);
            
            if(e % 1000 == 0) verb::cerr() <<"---  nvis : "<< trav.nvis() <<"\n";
            nv += trav.nvis();
        }
        std::cout <<"avg nvis : "<< nv / n <<" / "<< n <<"\n";
        
        verb::lap("skel");
    }


    // ----------------------- Radius heuristic -----------------

    if (do_rad) {

        if(g_scc_nb != 1) bye("The graph is not strongly connected.");
        
        ecc.clear();
        verb::lap("rad_with_sumsw");
        int64_t r = ecc.radius(graph::not_vertex, true /*do not optimize cert*/,
                               false, -1, true);
        double t_r = verb::lap_time() * 1000;
        set_col("rad_smsw_t", (int64_t)t_r);
        
        verb::cerr() << "rad_smsw : R=" << r << " D>=" << ecc.diam_lb
                  <<" nbfs="<< ecc.rad_nsweep
                  <<" ncertif="<< ecc.rad_certif.size()
                  <<" time="<< t_r
                  << std::endl;
    
        set_col("rad_smsw", r);
        set_col("rad_smsw_diam_lb", ecc.diam_lb);
        set_col("rad_smsw_nbfs", ecc.rad_nsweep);
        set_col("rad_smsw_cert", ecc.P.size());
        set_col("rad_smsw_cert_opt", ecc.rad_certif.size());
        set_col("rad_smsw_cert_nbfs", ecc.nsweep - ecc.rad_nsweep);

        ecc.clear();
        verb::lap("rad");
        r = ecc.radius(graph::not_vertex, true);
        t_r = verb::lap_time() * 1000;
        set_col("rad_t", (int64_t)t_r);
        
        verb::cerr() << "rad : R=" << r << " D>=" << ecc.diam_lb
                  <<" nbfs="<< ecc.rad_nsweep
                  <<" ncertif="<< ecc.rad_certif.size()
                  <<" time="<< t_r
                  << std::endl;
    
        set_col("rad", r, 0);
        set_col("rad_diam_lb", ecc.diam_lb);
        set_col("rad_nbfs", ecc.rad_nsweep);
        set_col("rad_cert", ecc.P.size());
        set_col("rad_cert_opt", ecc.rad_certif.size());
        set_col("rad_cert_nbfs", ecc.nsweep - ecc.rad_nsweep);

        traversal<graph> trav(n);
        std::vector<std::vector<int>> coballs;
        for (int u : ecc.Rcoballs) {
            trav.clear();
            if (weighted) trav.dijkstra(g, u); else trav.bfs(g, u);
            std::vector<int> cob;
            for (int i = trav.nvis() - 1; i >= 0; --i) {
                int v = trav.visit(i);
                if (trav.dist(v) < r) break;
                // else in coball of u
                cob.push_back(v);
            }
            coballs.push_back(cob);
        }
        
        int disj = ecc.gdy_disjoint(coballs, n).size();
        disj = std::max(disj, 2);
        verb::cerr() << "rad_disj: "<< disj <<"\n";
        set_col("rad_disj", disj);

        ecc.clear();
        verb::lap("rad-cert-app");
        assert(r == ecc.radius_certif_approx(r, sample_size));
        t_r = verb::lap_time() * 1000;
        verb::cerr() << "rad_cert_app: P="<< ecc.P.size()
                  <<" nbfs="<< ecc.rad_nsweep
                  <<" ncertif="<< ecc.rad_certif.size()
                  <<" time="<< t_r <<"\n";
        set_col("rad_cert_app", ecc.P.size());
        set_col("rad_cert_app_opt", ecc.rad_certif.size());
        verb::cerr() << "spl_sz: "<< sample_size <<"\n";
        set_col("spl_sz", sample_size);

        verb::lap("rad");
    }
    
    // --------------------- Various diameter heuristic ---------------

    if (do_diam) {

        if(g_scc_nb != 1) bye("The graph is not strongly connected.");

        ecc.clear();
        int start = graph::not_vertex;
        start_diam do_rad = RADIUS;
        int64_t d = -1;
        verb::lap("diam");
        if (algo == "bare") {
            d = ecc.diameter_bare(start, do_rad, 0, false);
        } else if (algo == "bare-smsw") {
            d = ecc.diameter_bare(start, SUMSWEEP, 0, false);
        } else if (algo == "sample") {
            d = ecc.diameter_sample(start, do_rad);
        } else if (algo == "bi") {
            d = ecc.diameter_bi(start, do_rad);
        } else if (algo == "tk11") {
            d = ecc.diameter_tk11(start, do_rad);
        } else {
            d = ecc.diameter_bare1(start, do_rad);
        }
        double t_d = verb::lap_time() * 1000;
        set_col("diam_t", (int64_t)t_d);

        verb::cerr() << "D=" << d
                  <<" nbfs="<< ecc.diam_nsweep
                  <<" ncertif="<< ecc.diam_certif.size()
                  <<" time="<< t_d
                  << std::endl;

        set_col("diam", d, 0);
        set_col("diam_nbfs", ecc.diam_nsweep);
        set_col("diam_cert", ecc.C.size());
        set_col("diam_cert_opt", ecc.diam_certif.size());
        set_col("diam_cert_nbfs", ecc.nsweep - ecc.diam_nsweep);

        verb::lap("diam");
    }
    
    // --------------------- Various diameter heuristic ---------------

    if (do_diam_all) {

        if(g_scc_nb != 1) bye("The graph is not strongly connected.");
        
        int best_nbfs = INT_MAX, besti = -1;
        std::string best = "";
        std::mutex col_mutex;
        
        auto concl = [&set_col, &best_nbfs, &besti, &best, &col_mutex]
                     (std::string diam, int i, eccentricity<graph> &ecc,
                      start_diam rad, int loo, int64_t d) {
            col_mutex.lock();
            diam = diam + "_" + std::to_string(rad) + "_" + std::to_string(loo);
            verb::cerr() << diam <<" : D=" << d
               <<" nbfs="<< ecc.diam_nsweep
               <<" ncertif="<< ecc.diam_certif.size()
               << std::endl;

            if (ecc.diam_nsweep < best_nbfs) {
                best_nbfs = ecc.diam_nsweep;
                best = diam;
                besti = i;
            }
            
            set_col(diam, d);
            set_col(diam +"_nbfs", ecc.diam_nsweep);
            set_col(diam +"_rad_nbfs", ecc.rad_nsweep);
            set_col(diam +"_cert", ecc.C.size());
            set_col(diam +"_P", ecc.P.size());
            set_col(diam +"_cert_opt", ecc.diam_certif.size());
            set_col(diam +"_cert_nbfs", ecc.nsweep - ecc.diam_nsweep);
            
            verb::cerr() <<"last_lb_improve " << ecc.last_lb_improve <<std::endl;
            set_col("last_lb_improve", ecc.last_lb_improve);
            
            col_mutex.unlock();
        };

        auto go_D = [&concl, &g, &g_rev,
                     directed, weighted](int i_thd, int n_thd) {
            eccentricity<graph> ecc(g, g_rev, directed, weighted);
            int novtx = graph::not_vertex;

            bool ch_ub = true;
            int a = -1;

            for (start_diam rad : {STRAIGHT, SUMSWEEP, RADIUS}) {
            if (i_thd == ++a % n_thd ) {
                ecc.clear();
                concl("diam_sample", a, ecc, rad, 0,
                      ecc.diameter_sample(novtx, rad, false));
            }
            if (i_thd == ++a % n_thd ) {
                ecc.clear();
                concl("diam_sample", a, ecc, rad, 1,
                      ecc.diameter_sample(novtx, rad, true));
            }
            for (int loose : {0, 1, 2}) {
            if (i_thd == ++a % n_thd ) {
                ecc.clear();
                concl("diam_prune", a, ecc, rad, loose,
                      ecc.diameter_prune(novtx, rad, ch_ub, true, loose));
            }
            if (i_thd == ++a % n_thd ) {
                ecc.clear();
                concl("diam", a, ecc, rad, loose,
                      ecc.diameter(novtx, rad, ch_ub, loose));
            }
            if (i_thd == ++a % n_thd ) {
                ecc.clear();
                concl("diam_bare", a, ecc, rad, loose,
                      ecc.diameter_bare(novtx, rad, loose));
            }
            if (i_thd == ++a % n_thd ) {
                ecc.clear();
                concl("diam_bare_p", a, ecc, rad, loose,
                      ecc.diameter_bare1(novtx, rad, true, true, loose));
            }
            if (i_thd == ++a % n_thd ) {
                ecc.clear();
                concl("diam_path_lb", a, ecc, rad, loose,
                      ecc.diameter_path_to_lb(novtx, rad, ch_ub, false, loose));
            }
            }
            }
            
        };
        if (n_thread <= 1) go_D(0, 1);
        else {
            int n_thd = std::min(8, n_thread);
            std::vector<std::thread> threads(n_thd);
            for (int i = 0; i < n_thd; ++i) {
                threads[i] = std::thread(go_D, i, n_thread);
            }
            verb::cerr() << n_thd <<" D threads\n";
            for (int i = 0; i < n_thd; ++i) threads[i].join();
        }

        verb::cerr() << "Best: "<< best <<"("<< besti
                  <<") with "<< best_nbfs <<" sweeps\n";
        set_col("best_diam_nbfs", besti);
        
        verb::lap("diam-all");
    }

    // --------------------- All ecc heuristic ---------------

    if (many_bfs > 0) {
        // Do many_bfs Dij/BFS on pseudo-random nodes

        traversal<graph> trav(n);

        verb::lap("many bfs");
        for (int i = 1; i <= many_bfs; ++i) {
            int u = (int)((1237L * i + 12347L) % n);
            trav.clear();
            if (weighted) trav.dijkstra(g, u); else trav.bfs(g, u);
            if (verb::progress())
                verb::cerr("many")
                    << i <<" / "<< many_bfs
                    <<" t_avg = " << (verb::lap_time() * 1000 / i) <<"ms"
                    <<" edges/s: "<< (1.0 * i * m / verb::lap_time())
                    <<"\n";
        }
        double t_all = verb::lap_time() * 1000;
        set_col("t_all", (int64_t)t_all);
        verb::cerr() << "t_all: " << t_all << " t_avg: " << (t_all / many_bfs)
                  <<" edges/s: "<< (1.0 * many_bfs * m / t_all)
                  <<"\n";

        verb::lap("many-bfs");
    }
    
    // --------------------- All ecc heuristic ---------------

    if (do_all_ecc) {

        if(g_scc_nb != 1) bye("The graph is not strongly connected.");
        
        ecc.clear();
        if (n_thread <= 1) {
            ecc.all();
        } else {
            verb::cerr() << n_thread <<" eccentricity threads\n";
            ecc.all_threaded(n_thread);
        }
        verb::cerr() <<"last_lb_improve " << ecc.last_lb_improve <<std::endl;

        set_col("last_lb_improve", ecc.last_lb_improve);

        verb::cerr() << "all ecc: "
                  <<" nbfs="<< ecc.all_ecc_nsweep
                  <<" ncertif="<< ecc.P.size() <<","<< ecc.C.size()
                  << std::endl;
        
        set_col("nall_nbfs", ecc.all_ecc_nsweep);
        set_col("nall_lb_P", ecc.P.size());
        set_col("nall_ub_C", ecc.C.size());
        set_col("nall_cert_lb", ecc.all_lb_certif.size());
        set_col("nall_cert_ub", ecc.all_ub_certif.size());
        set_col("nall_cert_nbfs",
                ecc.nsweep + ecc.P.size() - ecc.all_ecc_nsweep);

        // quite expansive with false, false :
        int d_cert_final = ecc.optim_ub_certif_one_shot(false, true).size();
        verb::cerr() <<"diam_cert_final " << d_cert_final << std::endl;
        set_col("diam_cert_final", d_cert_final);

        // Centers :
        verb::cerr() <<"centers:";
        int nb_centers = 0;
        for (int u : g)
            if (ecc.ecc(u) == ecc.rad_ub)
                { verb::cerr() <<" "<< lab[u]; ++nb_centers; }
        verb::cerr() <<"\n";
        set_col("nb_centers", nb_centers);
        
        // R cert :
        verb::cerr() <<"R_cert:";
        for (int u : ecc.rad_certif) verb::cerr() <<" "<< lab[u];
        verb::cerr() <<"\n";
        // ub cert :
        verb::cerr() <<"D_cert:";
        for (int u : ecc.diam_certif) verb::cerr() <<" "<< lab[u];
        verb::cerr() <<"\n";
        // lb cert :
        verb::cerr() <<"lb_cert:";
        for (int u : ecc.all_lb_certif) verb::cerr() <<" "<< lab[u];
        verb::cerr() <<"\n";
        // ub cert :
        verb::cerr() <<"ub_cert:";
        for (int u : ecc.all_ub_certif) verb::cerr() <<" "<< lab[u];
        verb::cerr() <<"\n";

        // stats :
        std::vector<int64_t> e(n);
        double sum = 0.;
        for (int u = 0; u < n; ++u) {
            e[u] = ecc.ecc_lb(u);
            std::cout << lab[u] <<" "<< e[u] <<"\n";
            sum += (double)e[u];
        }
        set_distr_col(e, "ecc");
        
        verb::cerr() <<"avg ecc "<< (sum / (double)n) << std::endl;
        
        verb::lap("all_ecc");
    }

    // ------------------- All BFS -----------
    if (all_bfs) {

        /** How many nodes are furthest from some node u ?
         *         according to various tie break rules:
         * rnd_perm: last in a given random permutation
         * last_bfs: last in the BFS/Dijkstra from u
         * max_idx: last encountered node in the graph file
         * bfs_ctr: last in a BFS from a center
         * bfs_diam: last in a BFS from a diametral node
         * min_idx: first encountered node in the graph fie
         * rnd_bfs: last in a BFS from a random node
         * periph_bfs: last in a BFS from the last visited from the rnd node
         * max_lab: node with maximum label in the graph file
         * min_lab: node with minimum label in the graph file
         *
         * nb_furthest: number of nodes u at distance ecc(v) from some v
         */
        
        if(g_scc_nb != 1) bye("The graph is not strongly connected.");
        
        traversal<graph> trav(n), ctrav(n), ctrav_fwd(n);
        std::mutex mutex;

        ecc.clear();
        int pseudo_center = ecc.sum_sweep();

        ecc.clear();
        int64_t radius = ecc.radius();
        int center = ecc.rad_node;
        int64_t diameter = ecc.diameter_sample(graph::not_vertex, STRAIGHT);

        // Need all eccentricities for full stats
        ecc.clear();

        std::vector<int> certif; // lower certificate
        if (optim_certif != "") {
            if (optim_certif == "-") { in = stdin; }
            else {
                verb::cerr() << "open '" << optim_certif << "' to read certif\n";
                in = fopen(optim_certif.c_str(), "r");
            }
            if (in != nullptr) {
                char u[1024];
                for ( ; fscanf(in, " %s \n", u) >= 1 ; ) {
                    assert(vi[u] > 0);
                    certif.push_back(vi[u]-1);
                }
                fclose(in);
            }
            verb::cerr() << "certif size: " << certif.size() << std::endl;
            for (int u : certif) verb::cerr() << u <<" ";
            verb::cerr() << std::endl;

            verb::cerr() << "file t_all: " + (optim_certif + ".t_all") <<"\n";
            FILE *in = fopen((optim_certif + ".t_all").c_str(), "r");
            int64_t t_all = -1;
            if (in != nullptr) {
                fscanf(in, "%lld\n", &t_all);
                fclose(in);
            }
            set_col("t_all", t_all);
            verb::cerr() << "t_all: " << t_all <<"\n";
        }

        if (certif.size() > 0) {
            verb::cerr() <<"all ecc from lb certif -----------------------\n";
            ecc.all_from_lb_cert(certif, n_thread, false);
        } else {
            verb::lap("beg all ecc");
            ecc.all_threaded(n_thread, graph::not_vertex, false, false);
            double t_all = verb::lap_time() * 1000;
            set_col("t_all", (int64_t)t_all);
            verb::cerr() << "t_all: " << t_all <<"\n";
            
            if (optim_certif != "") { // save it
                FILE *out = fopen(optim_certif.c_str(), "w");
                for (int u : ecc.all_lb_certif) {
                    fprintf(out, "%s\n", lab[u].c_str());
                }
                fclose(out);
                out = fopen((optim_certif + ".t_all").c_str(), "w");
                fprintf(out, "%lld\n", (int64_t)t_all);
                fclose(out);
            }
        }
            
        verb::cerr() << "rad: "<< ecc.rad_ub <<"  diam: "<< ecc.diam_lb
                  << "  bfs: " << ecc.all_ecc_nsweep
                  << "  all cert: "
                  << ecc.P.size() <<", " << ecc.C.size()
                  << "  all cert opt: "
                  << ecc.all_lb_certif.size() <<", " << ecc.all_ub_certif.size()
                  <<" last_lb_improve: " << ecc.last_lb_improve
                  << std::endl;

        // Centers :
        std::cerr <<"centers:";
        int nb_centers = 0;
        for (int u : g) {
            if (ecc.ecc(u) == ecc.rad_ub) {
                ++nb_centers;
                if (nb_centers <= 100) std::cerr <<" "<< lab[u];
                if (nb_centers == 101) std::cerr <<",... ";
            }
        }
        std::cerr <<" ; nb_centers="<< nb_centers <<"\n";
        // R cert :
        std::cerr <<"R_cert:";
        //for (int u : ecc.rad_certif) std::cerr <<" "<< lab[u];
        std::cerr <<"\n";
        // D cert :
        std::cerr <<"D_cert:";
        //for (int u : ecc.diam_certif) std::cerr <<" "<< lab[u];
        std::cerr <<"\n";
        // lb cert :
        std::cerr <<"lb_cert:";
        for (int u : ecc.all_lb_certif) std::cerr <<" "<< lab[u];
        std::cerr <<"\n";

        set_col("last_lb_improve", ecc.last_lb_improve);
        set_col("all_ecc_nsweep", ecc.all_ecc_nsweep);
        set_col("all_ecc_cert_lb", ecc.P.size());
        set_col("all_ecc_cert_ub", ecc.C.size());
        set_col("all_ecc_cert_lb_opt", ecc.all_lb_certif.size());
        set_col("all_ecc_cert_ub_opt", ecc.all_ub_certif.size());
        set_col("all_ecc_cert_rad", ecc.rad_certif.size());
        set_col("all_ecc_cert_diam", ecc.diam_certif.size());
        int64_t rad = ecc.rad_ub;
        set_col("rad", rad);
        int64_t diam = ecc.diam_lb;
        set_col("diam", diam);

        std::vector<int> disj = ecc.gdy_disjoint(ecc.Pcoballs, n);
        std::cerr << "all_ecc_disj: "<< disj.size() <<"\n";
        set_col("all_ecc_disj", disj.size());

        verb::lap("all_bfs_eccs");

        int cert_center = center; // --------- use true center !!!!!!!!!!
        verb::cerr() << "cert_center lab: "<< lab[cert_center]
                      <<" ("<< cert_center <<")\n"; 
        ctrav.clear();
        if (weighted) ctrav.dijkstra(g_rev, cert_center);
        else ctrav.bfs(g_rev, cert_center);
        ctrav_fwd.clear();
        if (weighted) ctrav_fwd.dijkstra(g, cert_center);
        else ctrav_fwd.bfs(g, cert_center);
        auto in_c_scc = [&ctrav, &ctrav_fwd](int v) {
            return ctrav.visited(v) && ctrav_fwd.visited(v);
        };
        int64_t cert_rad = ctrav_fwd.dist(ctrav_fwd.last_visited());
        std::cerr <<"R, certR, D: "<< ecc.rad_ub <<", "<< cert_rad
                  <<", "<< ecc.diam_lb <<"\n";

        assert(rad == ecc.rad_ub);
        assert(diam == ecc.diam_lb);
        
        std::vector<int> ctr_uncov; // nodes outside B(c, D - R)
        for (int v = 0; v < n; ++v) {
            if (ctrav.dist(v) > diam - cert_rad || ! in_c_scc(v))
                ctr_uncov.push_back(v);
        }

        std::vector<int> C_opt;
        std::vector<bool> in_C(n, false);
        for (int c : ecc.C) in_C[c] = true;
        std::vector<bool> one_cert1(n);
        
        { // release some memory after
            
            std::vector<std::vector<int>> diam_cert_of(n);
            // Mesure packing if radius of balls of certificates is divd by 2/3
            std::vector<std::vector<int>> diam1_cert_of(n),
                diam2_cert_of(n), diam3_cert_of(n);
            double center1_fact = 1.0 / 0.8, center2_fact = 2.0,
                center3_fact = 3.0;
            int n_ctr = 0, n_ctr1 = 0, n_ctr2 = 0, n_ctr3 = 0;

            { // release diam_cert
            std::vector<std::vector<int>> diam_cert(n);
            // Mesure packing if radius of balls of certificates is divd by 2/3
            std::vector<std::vector<int>>  diam1_cert(n),
                diam2_cert(n), diam3_cert(n);

            int64_t nprune = 0, mprune = 0;
            // (use theorem here) :
            auto go_cert = [n, &one_cert1, &C_opt, &in_C, &ctrav, &ecc,
                            cert_rad, cert_center, &in_c_scc,
                            &nprune, &mprune, &mutex](int i_thd, int n_thd) {
                traversal<graph> travloc(n);
                int nsweep, mprune_loc, nprune_loc;
                for (int u = 0; u < n; ++u) {
                    if ((1237L * u) % n_thd == i_thd &&
                        (in_C[u] || ctrav.dist(u) > ecc.diam_lb - cert_rad
                         || ! in_c_scc(u))) {
                        ecc.pruned_sweep_to_certifiers_threaded
                            (u, travloc, nsweep, mprune_loc, nprune_loc);
                        one_cert1[u] = travloc.nvis() == 1;
                        if (one_cert1[u]) {
                            mutex.lock();
                            C_opt.push_back(u);
                            mutex.unlock();
                        }
                    }
                    if (i_thd == 0 && verb::progress())
                        verb::cerr("prune cert")
                            <<"todo: "<< (n - u) <<" / "<< n <<"\n";
                }
                mutex.lock();
                nprune += nprune_loc;
                mprune += mprune_loc;
                mutex.unlock();
            };
            auto go_diam = [n, &diam_cert, cert_rad, cert_center,
                            &diam1_cert, &diam2_cert, &diam3_cert,
                            center1_fact, center2_fact, center3_fact,
                            &in_C, &ctrav, &ecc, &low_cert, mem_limit_mb,
                            &in_c_scc, &very_low_cert,
                            &nprune, &mprune, &mutex](int i_thd, int n_thd) {
                traversal<graph> travloc(n);
                bool low_cert_loc = low_cert;
                int nsweep = 0, mprune_loc = 0, nprune_loc = 0,
                    nsweep_null = 0, mprune_null = 0, nprune_null = 0;
                for (int u = 0; u < n; ++u) {
                    if ((1237L * u) % n_thd == i_thd
                        && ((! in_c_scc(u))
                            || ctrav.dist(u) >= ecc.diam_lb - cert_rad
                            || (center3_fact * ctrav.dist(u)
                                   > ecc.diam_lb - cert_rad
                                && ! (low_cert_loc
                                      || verb::mem_now_kb()/1000 > mem_limit_mb)
                                ))) {
                        double dfact = center3_fact;
                        if (center2_fact * ctrav.dist(u)
                            > ecc.diam_lb - cert_rad)
                            dfact = center2_fact;
                        if (center1_fact * ctrav.dist(u)
                            > ecc.diam_lb - cert_rad)
                            dfact = center1_fact;

                        if (ctrav.dist(u) > ecc.diam_lb - cert_rad
                                || !  in_c_scc(u))
                            ecc.pruned_sweep_to_certifiers_threaded
                                (u, travloc, nsweep, mprune_loc, nprune_loc,
                                 ecc.diam_lb);
                        else
                            ecc.pruned_sweep_to_certifiers_threaded
                                (u, travloc, nsweep_null, mprune_null,
                                  nprune_null, ecc.diam_lb, dfact);

                        for (int i = travloc.nvis() - 1; i >= 0; --i) {
                            int c = travloc.visit(i);
                            if ((ctrav.dist(u) > ecc.diam_lb - cert_rad
                                 || ! in_c_scc(u))
                                && ! very_low_cert)
                                diam_cert[u].push_back(c);
                            if ((ctrav.dist(u) >= ecc.diam_lb - cert_rad
                                 || (! in_c_scc(u))
                                 || (center1_fact * ctrav.dist(u)
                                     > ecc.diam_lb - cert_rad))
                                &&
                                (travloc.dist(c) == 0
                                 || center1_fact * travloc.dist(c)
                                 < ecc.diam_lb - ecc.ecc(c))
                                )//&& ! very_low_cert)
                                diam1_cert[u].push_back(c);
                            if ((ctrav.dist(u) >= ecc.diam_lb - cert_rad
                                 || center2_fact * ctrav.dist(u)
                                     > ecc.diam_lb - cert_rad
                                 || ! in_c_scc(u))
                                &&
                                  (travloc.dist(c) == 0
                                   || center2_fact * travloc.dist(c)
                                   < ecc.diam_lb - ecc.ecc(c))
                                && ! very_low_cert)
                                diam2_cert[u].push_back(c);
                            if ((ctrav.dist(u) >= ecc.diam_lb - cert_rad
                                 || center3_fact * ctrav.dist(u)
                                     > ecc.diam_lb - cert_rad
                                 || ! in_c_scc(u))
                                &&
                                (travloc.dist(c) == 0
                                 || center3_fact * travloc.dist(c)
                                 < ecc.diam_lb - ecc.ecc(c)))
                                diam3_cert[u].push_back(c);
                        }
                        
                    }
                    if (( ! low_cert_loc) &&
                        (low_cert || verb::mem_now_kb()/1000 > mem_limit_mb)) {
                        if ( ! low_cert) {
                            verb::cerr()
                                << "Switching to -low-cert mode (low mem)!\n";
                            low_cert = true;
                        }
                        low_cert_loc = true;
                        // free memory
                        for (int v = 0; v <= u; ++v) {
                            if ((1237L * v) % n_thd == i_thd
                                && (in_c_scc(v)
                                    && ctrav.dist(v) < ecc.diam_lb - cert_rad)
                                ) {
                                diam1_cert[v] = {};
                                diam2_cert[v] = {};
                                diam3_cert[v] = {};
                            }
                        }
                    }
                    if (i_thd == 0 && verb::progress())
                        verb::cerr("prune diam")
                            <<"todo: "<< (n - u) <<" / "<< n <<"\n";
                }
                mutex.lock();
                nprune += nprune_loc;
                mprune += mprune_loc;
                mutex.unlock();
            };
            
            verb::cerr() << "prune cert\n";
            if (n_thread <= 1) go_cert(0, 1);
            else {
                std::vector<std::thread> threads(n_thread);
                for (int i = 0; i < n_thread; ++i) {
                    threads[i] = std::thread(go_cert, i, n_thread);
                }
                std::cerr << n_thread <<" cert threads\n";
                for (int i = 0; i < n_thread; ++i) threads[i].join();
            }
            
            while (mprune > g.m()) { mprune -= g.m()+1; ++nprune; }
            nprune += (mprune > 0 ? 1 : 0);
            std::cerr << "n_prune for cert :" << nprune << std::endl;
            set_col("n_prune_cert", nprune);
            
            verb::cerr() << "prune diam\n";
            nprune = 0; mprune = 0;
            if (n_thread <= 1) go_diam(0, 1);
            else {
                std::vector<std::thread> threads(n_thread);
                for (int i = 0; i < n_thread; ++i) {
                    threads[i] = std::thread(go_diam, i, n_thread);
                }
                std::cerr << n_thread <<" diam threads\n";
                for (int i = 0; i < n_thread; ++i) threads[i].join();
            }
            
            while (mprune > g.m()) { mprune -= g.m()+1; ++nprune; }
            nprune += (mprune > 0 ? 1 : 0);
            ecc.mprune = 0;
            std::cerr << "n_prune for diam :" << nprune << std::endl;
            set_col("n_prune_diam", nprune);
            
            verb::cerr() << "invert\n";
            assert(ctrav.first_visited() == cert_center);
            for (int u = 0; u < n; ++u) {
                for (int c : diam_cert[u]) diam_cert_of[c].push_back(u);
                diam_cert[u] = {}; // release mem
                for (int c : diam1_cert[u]) diam1_cert_of[c].push_back(u);
                for (int c : diam1_cert[u])
                diam1_cert[u] = {}; // release mem
                for (int c : diam2_cert[u]) diam2_cert_of[c].push_back(u);
                diam2_cert[u] = {}; // release mem
                for (int c : diam3_cert[u]) diam3_cert_of[c].push_back(u);
                diam3_cert[u] = {}; // release mem
                if (verb::progress())
                    verb::cerr("cert") <<"todo: "<< (n - u) <<" / "<< n <<"\n";
            }

            } // release diam_cert

            std::vector<int> tocov1, tocov2, tocov3;
            for (int u = 0; u < n; ++u) {
                if (ctrav.dist(u) <= ecc.diam_lb - cert_rad && in_c_scc(u)) {
                    //diam_cert_of[cert_center].push_back(u);
                    ++n_ctr;
                }
                if (center1_fact * ctrav.dist(u) <= ecc.diam_lb - cert_rad
                     && in_c_scc(u)) {
                    //diam1_cert_of[cert_center].push_back(u);
                    ++n_ctr1;
                } else tocov1.push_back(u);
                if (center2_fact * ctrav.dist(u) <= ecc.diam_lb - cert_rad
                     && in_c_scc(u)) {
                    //diam2_cert_of[cert_center].push_back(u);
                    ++n_ctr2;
                } else tocov2.push_back(u);
                if (center3_fact * ctrav.dist(u) <= ecc.diam_lb - cert_rad
                     && in_c_scc(u)) {
                    //diam3_cert_of[cert_center].push_back(u);
                    ++n_ctr3;
                } else tocov3.push_back(u);
            }

            verb::lap("all_bfs_pruning");

            std::cerr << "n_ctr1,2,3 : " << n_ctr
                      <<" "<< n_ctr2 <<" "<< n_ctr3 <<" / "<< n <<"\n";

            set_col("n_ctr", n_ctr);
            set_col("n_ctr2", n_ctr2);
            set_col("n_ctr3", n_ctr3);
        
            int one_ub_cert = 0;
            for (int u = 0; u < n; ++u) {
                if (one_cert1[u] == 1) ++one_ub_cert;
            }
            std::cerr << "one certif only: " << one_ub_cert << std::endl;
            set_col("one_ub_cert", one_ub_cert);


            int D1 = -1, D1_inf = -1, D1bis = -1;
            std::vector<int> D1cert = {};
            
            if ( ! very_low_cert) {
                D1cert = ecc.gdy_set_cov(diam_cert_of, n, ctr_uncov);
                D1 = 1 + D1cert.size();
                std::cerr <<"D1_inf, D1:   " << D1_inf <<", "<< D1 << "\n";
                D1_inf = ecc.gdy_packing(diam_cert_of, n, ctr_uncov).size();
                D1_inf = std::max(1, D1_inf);
                std::cerr <<"D1_inf, D1:   " << D1_inf <<", "<< D1 << "\n";
                std::cerr <<"low_cert=" << low_cert <<"\n";
                if ( ! low_cert) {
                    D1bis = 1 + ecc.gdy_set_cov(diam1_cert_of,
                                                n, tocov1).size();
                }

            }
            

            // D cert :
            std::cerr <<"D1_cert: ";
            std::cerr << lab[cert_center];
            for (int u : D1cert) std::cerr <<" "<< lab[u];
            std::cerr <<"\n";
            
            std::cerr <<"D1_inf, D1:   " << D1_inf <<", "<< D1 << "\n";

            set_col("D1_inf", D1_inf);
            set_col("D1", D1);

            std::cerr <<"D1.25:   " << D1bis <<"\n";
            
            int D2 = -1;
            if ( ! low_cert) {
                D2 = 1 + ecc.gdy_set_cov(diam2_cert_of, n, tocov2).size();
            }
            std::cerr <<"D2:   " << D2 <<"\n";
            int D3 = -1;
            if ( ! low_cert) {
                D3 = 1 + ecc.gdy_set_cov(diam3_cert_of, n, tocov3).size();
            }
        
            std::cerr <<"D1.25, D2, D3:   " << D1bis
                      <<", "<< D2 <<", "<< D3 << "\n";

            set_col("D1.25", D1bis);
            set_col("D2", D2);
            set_col("D3", D3);

            // Using center:
            std::vector<int> ctr_unstr; // nodes outside B(c,D-R-1)
            for (int v = 0; v < n; ++v) {
                if (ctrav.dist(v) >= diam - cert_rad || ! in_c_scc(v)) {
                    ctr_unstr.push_back(v);
                }
            }

            int D1ctr = -1, D1ctr_inf = -1;
            //if ( ! very_low_cert) {
                D1ctr = 1 + ecc.gdy_set_cov(diam1_cert_of, n, ctr_unstr).size();
                D1ctr_inf = ecc.gdy_packing(diam1_cert_of, n, ctr_unstr).size();
                D1ctr_inf = std::max(1, D1ctr_inf);
                //}
            std::cerr <<"D1.25ctr\n";
            
            int D2ctr = -1, D2ctr_inf = -1;
            if ( ! very_low_cert) {
                D2ctr = 1 + ecc.gdy_set_cov(diam2_cert_of, n, ctr_unstr).size();
                D2ctr_inf = ecc.gdy_packing(diam2_cert_of, n, ctr_unstr).size();
                D2ctr_inf = std::max(1, D2ctr_inf);
            }
            std::cerr <<"D2ctr\n";

            std::vector<int> D3ctr_cov =
                ecc.gdy_set_cov(diam3_cert_of, n, ctr_unstr);
            std::cerr <<"D3ctr_cert: "<< cert_center;
            for (int u : D3ctr_cov) std::cerr <<" "<< lab[u];
            std::cerr <<"\n";
            int D3ctr = 1 + ecc.gdy_set_cov(diam3_cert_of, n, ctr_unstr).size();
            int D3ctr_inf = ecc.gdy_packing(diam3_cert_of, n, ctr_unstr).size();
            D3ctr_inf = std::max(1, D3ctr_inf);
            std::cerr <<"D3ctr\n";
        
            std::cerr <<"D1.25ctr, D2ctr, D3ctr:   "
                      << D1ctr <<", "<< D2ctr <<", "<< D3ctr << "\n"
                      <<"D1.25ctr_inf, D2ctr_inf, D3ctr_inf:   "
                      << D1ctr_inf <<", "<< D2ctr_inf <<", "<< D3ctr_inf << "\n";

            set_col("D1.25ctr_inf", D1ctr_inf);
            set_col("D2ctr_inf", D2ctr_inf);
            set_col("D3ctr_inf", D3ctr_inf);

            set_col("D1.25ctr", D1ctr);
            set_col("D2ctr", D2ctr);
            set_col("D3ctr", D3ctr);

        
            std::cerr <<"upp. cert C_opt : "<< C_opt.size() << std::endl;
            set_col("Copt", C_opt.size());

       } // release big arrays

        
        // ------------ eccentricities stats --------------------
        std::vector<int64_t> e(n), e_one;
        double sum = 0., sum_one = 0.;
        int n_ecc_min = 0, n_ecc_min1 = 0, n_ecc_min2 = 0;
        for (int u = 0; u < n; ++u) {
            e[u] = ecc.ecc(u);
            if (e[u] <= rad) ++n_ecc_min;
            if (e[u] <= rad + 1) ++n_ecc_min1;
            if (e[u] <= rad + 2) ++n_ecc_min2;
            sum += (double)e[u];
            if (one_cert1[u]) {
                e_one.push_back(ecc.ecc(u));
                sum_one += (double)ecc.ecc(u);
            }
        }
        set_distr_col(e, "ecc");
        set_distr_col(e_one, "ecc_one");

        set_col("n_ecc_min", n_ecc_min);
        set_col("n_ecc_min1", n_ecc_min1);
        set_col("n_ecc_min2", n_ecc_min2);

        int m_e_one = e_one.size();

        std::cerr <<"avg ecc "<< (sum / (double)n)
                  <<"  avg ecc_one "<< (sum_one / (double)m_e_one)
                  << std::endl;

        // ------------- true center for ctrav ----------------
        
        ctr_uncov.clear(); // nodes outside B(c, D - R)
        for (int v = 0; v < n; ++v) {
            if (ctrav.dist(v) > diam - rad || ! in_c_scc(v))
                ctr_uncov.push_back(v);
        }

        // --- stats for pruning : nb of neighbors with eccentricity ecc - 1
        std::vector<int64_t> neighb_ecc_m1(n, 0), neighb_ecc_m1_far;
        // for nodes not certified by a center
        double sum_neighb = 0.;
        int neighb_n = 0;
        for (int u = 0; u < n; ++u) {
            int64_t e_u = ecc.ecc(u);
            bool far = ctrav.dist(u) > ecc.diam_lb - rad || ! in_c_scc(u);
            int nb = 0;
            for (auto e : g[u]) {
                if (ecc.ecc(e.dst) + e.wgt == e_u) ++nb;
            }
            neighb_ecc_m1[u] = nb;
            if (far) {
                neighb_n += 1;
                neighb_ecc_m1_far.push_back(nb);
            }
            sum_neighb += neighb_ecc_m1[u];
        }
        set_distr_col(neighb_ecc_m1, "neighb"); //"nb_neighb_ecc_m1");
        set_distr_col(neighb_ecc_m1_far, "neighb_far"); //nb_neighb_ecc_m1_far

        if (neighb_n == 0) {
            neighb_ecc_m1_far.push_back(-1);
            neighb_n = 1;
        }        
        
        // ---- orderings        

        std::unordered_map<std::string, std::vector<int>>order;

        // random permutation :
        std::vector<int> perm(n);
        for (int i = 0; i < n; ++i) perm[i] = i; 
        for (int i = n - 1; i > 0; --i) {
            int j = rand() % (i + 1);
            std::swap(perm[i], perm[j]);
        }
        order["rnd_perm"] = perm;

        std::vector<int> visit_inv(n,-1);
        auto trav_order = [&trav, n, &perm](std::vector<int> &ord) {
            int nv = trav.nvis();
            for (int i = 0; i < nv; ++i) ord[trav.visit(i)] = i;
            for (int u = 0; u < n; ++u) {
                int i = ord[u];
                if ( ! (i >= 0 && i < nv && trav.visit(i) == u)) {
                    ord[u] = nv + perm[u];
                }
            }
        };

        // Last BFS
        std::vector<int> last_bfs(n,0); // same pos, no one beats the last
        order["last_bfs"] = last_bfs;

        // Rand BFS
        std::vector<int> rnd_bfs(n);
        int u = rand() % n;
        trav.clear();
        if (weighted) trav.dijkstra(g, u); else trav.bfs(g, u);
        trav_order(rnd_bfs);
        order["rnd_bfs"] = rnd_bfs;

        // BFS from periph node:
        std::vector<int> periph_bfs(n);
        int p = trav.last_visited();
        trav.clear();
        if (weighted) trav.dijkstra(g, p); else trav.bfs(g, p);
        trav_order(periph_bfs);
        order["periph_bfs"] = periph_bfs;

        // BFS order from center :
        // ecc.radius();
        trav.clear();
        if (weighted) trav.dijkstra(g, center); else trav.bfs(g, center);
        std::vector<int> bfs_ctr(n);
        trav_order(bfs_ctr);
        order["bfs_ctr"] = bfs_ctr;

        // BFS from diametral node
        // ecc.diameter();
        int diam_node = ecc.diam_node;
        trav.clear();
        if (weighted) trav.dijkstra(g, diam_node); else trav.bfs(g, diam_node);
        std::vector<int> bfs_diam(n);
        trav_order(bfs_diam);
        order["bfs_diam"] = bfs_diam;

        // Min index:
        std::vector<int> min_idx(n);
        for (int u = 0; u < n; ++u) min_idx[u] = n - u;
        order["min_idx"] = min_idx;

        // Max index:
        std::vector<int> max_idx(n);
        for (int u = 0; u < n; ++u) max_idx[u] = u;
        order["max_idx"] = max_idx;

    std::cerr<<"6 end big\n"; std::cerr.flush();

        // Labels:
        std::vector<std::pair<int64_t, int>> lab_idx(n);
        try {
            for (int u = 0; u < n; ++u) lab_idx[u] = { std::stoll(lab[u]), u };
        } catch (...) {
            std::vector<std::pair<std::string,int>> li(n);
            for (int u = 0; u < n; ++u) li[u] = { lab[u], u };
            std::sort(li.begin(), li.end());
            for (int i = 0; i < n; ++i) {
                int u = li[i].second;
                lab_idx[u] = { i, u };
            }
        }
    std::cerr<<"7 end big\n"; std::cerr.flush();
        std::sort(lab_idx.begin(), lab_idx.end());
        std::vector<int> min_lab(n), max_lab(n);
        for (int i = 0; i < n; ++i) {
            int u = lab_idx[i].second;
            min_lab[u] = n - i;
            max_lab[u] = i;
        }
        order["min_lab"] = min_lab;
        order["max_lab"] = max_lab;

     std::cerr<<"8 end big\n"; std::cerr.flush();
        
        std::vector<std::atomic<int>> is_last_bfs(n);
        for (int u = 0; u < n; ++u) is_last_bfs[u] = false;
        
        std::unordered_map<std::string,
                           std::vector<std::atomic<int>>> is_last, is_cert;
        std::unordered_map<std::string, int> n_last, n_cert;
        for(const auto &entry : order) {
            is_last[entry.first] = std::vector<std::atomic<int>>(n);
            for (int u = 0; u < n; ++u) is_last[entry.first][u] = false;
            n_last[entry.first] = 0;
            is_cert[entry.first] = std::vector<std::atomic<int>>(n);
            for (int u = 0; u < n; ++u) is_cert[entry.first][u] = false;
            n_cert[entry.first] = 0;
        }

        std::cerr<<"9 end big\n"; std::cerr.flush();

            
        // ---- end orderings        
        
        if (n > n_big_graph) { // -------------- big graph

            // fill unset cols with -1
            for (const auto & entry : order) {
                set_col(entry.first + "_lst", -1);
            }
            for (const auto & entry : order) {
                set_col(entry.first + "_crt", -1);
            }
            set_col("nb_furthest", -1);
            set_col("nb_far_cert", -1);
            set_col("nb_one_cert", -1);
            set_col("one_ub_ext", -1);
            set_col("P", -1);
            set_col("Pext", -1);
            set_col("C", C_opt.size());
            set_col("Cext", -1);
            set_col("P_inf", -1);
            set_col("Pext_inf", -1);
            
            set_col("C_inf", C_opt.size());
            set_col("Cext_inf", -1);
            
        } else {          // -------------- small graph

            ecc.make_dist_lab();

            std::cerr<<"dist labels to P: "<< ecc.P.size() <<"\n";
            std::cerr.flush();
            
            // is at distance ecc(v) for some v
            std::vector<std::atomic<int>> is_furthest(n),
                low_ecc(n), is_last_low_ecc(n),
                is_far_cert(n), one_cert(n);
            int nb_furthest = 0, nb_far_cert = 0, nb_one_cert = 0;
            // furthest from nodes at distance > R from antipode :
            std::vector<int64_t> dist_far_antipode(n, 0), quad_antipode(n, 0),
                rad_triangle(n, 0), ecc_max_dipole(n, 0), n_inter_dipole(n, 0);
            std::vector<int64_t> antipode_hyperbol(n, 0), rad_coball_size(n, 0);
            std::vector<std::vector<int>> rad_coball(n);
            std::vector<std::vector<int64_t>> rad_coball_dist(n);

            std::cerr<<"9.5 vects\n"; std::cerr.flush();

            // --------- x % closest to centers :
            int64_t ecc5pct = e[5*n/100];
            std::vector<int> the_centers;
            for (int u = 0; u < n; ++u) {
                if (ecc.ecc(u) <= ecc5pct) low_ecc[u] = true; 
            }
            
            std::cerr<<"10 bef count\n"; std::cerr.flush();

            auto count = [n, &is_last, &n_last, &is_furthest, &nb_furthest,
                          &n_cert, &is_cert, &is_far_cert, &nb_far_cert,
                          &one_cert, &nb_one_cert, &is_last_bfs]() {
                for (auto &entry : is_last) {
                    int cnt = 0;
                    for (int u = 0; u < n; ++u) {
                        if (entry.second[u].load(std::memory_order_relaxed))
                            ++cnt;
                    }
                    n_last[entry.first] = cnt;
                    std::cerr << entry.first <<"_l:"<< cnt <<" ";
                }
                for (auto &entry : is_cert) {
                    int cnt = 0;
                    for (int u = 0; u < n; ++u) {
                        if (entry.second[u].load(std::memory_order_relaxed))
                            ++cnt;
                    }
                    n_cert[entry.first] = cnt;
                    std::cerr << entry.first <<"_c:"<< cnt <<" ";
                }
                nb_furthest = 0; nb_far_cert = 0; nb_one_cert = 0;
                int nb_last_bfs = 0;
                for (int u = 0; u < n; ++u){
                    if (is_furthest[u].load(std::memory_order_relaxed))
                        ++nb_furthest;
                    if (is_far_cert[u].load(std::memory_order_relaxed))
                        ++nb_far_cert;
                    if (one_cert[u].load(std::memory_order_relaxed))
                        ++nb_one_cert;
                    if (is_last_bfs[u].load(std::memory_order_relaxed))
                        ++nb_last_bfs;
                }
                std::cerr << "nb_furthest:" << nb_furthest
                <<" nb_far_cert:" << nb_far_cert
                <<" nb_one_cert:" << nb_one_cert
                <<" nb_last_bfs:" << nb_last_bfs
                << std::endl;
                std::cerr.flush();
            };

            // Certifiers
            std::vector<std::vector<int>> lb_cert_of(n), lb_ext_of(n),
                ub_cert_of(n), ub_ext_of(n), rad_cert_of(n), diam_cert_of(n);
            std::vector<int> n_ub_ext(n);
            std::vector<int64_t> n_furthest(n), n_visit(n); // to call prt_distr
            
        std::cerr<<"11 small graph\n"; std::cerr.flush();

            int r_lim = 80000; // store rad cert if not certified by diam pair

            traversal<graph> atrav(n), btrav(n);
            int diam_a = ecc.diam_node;
            atrav.clear();
            if (weighted) atrav.dijkstra(g, diam_a); else atrav.bfs(g, diam_a);
            int diam_b = atrav.last_visited();
            btrav.clear();
            if (weighted) btrav.dijkstra(g, diam_b); else btrav.bfs(g, diam_b);
            
            int d_lim = 400000; // compute diam cert for nodes outside B(c,R)
            int64_t m_dcert_lim = 500000000, m_dcert_max = 0;
            
            auto thread_for = [&](int i_thd, int n_thd) { 
                traversal<graph> trav(n), trav_reverse(directed ? n : 0);
                int64_t m_diam_cert = 0, m_lb_cert = 0,
                        m_lb_ext = 0, m_ub_ext = 0;
                std::unordered_map<std::string, int> f_max;
                for(const auto &entry : order) {
                    f_max[entry.first] = graph::not_vertex;
                }
                
                //std::cerr<<"12 go\n"; std::cerr.flush();

            for (int u = 0; u < n; ++u) {
              if (u % n_thd == i_thd) {
                rad_coball[u] = std::vector<int>();
                  
                trav.clear();
                if (weighted) trav.dijkstra(g, u); else trav.bfs(g, u);
                int nvis = trav.nvis();
                if (directed) {
                    trav_reverse.clear();
                    if (weighted) trav_reverse.dijkstra(g_rev, u);
                    else trav_reverse.bfs(g_rev, u);
                }
                traversal<graph> &trav_rev = directed ? trav_reverse : trav;
                n_visit[u] = nvis;
                //assert(trav.nvis() == n);

                int f = trav.last_visited();
                int64_t e_u = trav.dist(f);
                if(e_u != ecc.ecc(u)) {
                    mutex.lock();
                    std::cerr <<" u="<< u <<" e_rev_u="
                              << trav_rev.dist(trav_rev.last_visited())
                              <<" e_u="<< e_u <<" != "<< ecc.ecc(u) <<"\n";
                    mutex.unlock();
                }
                assert(e_u == ecc.ecc(u));
                assert(e_u <= ecc.diam_lb);
                assert(e_u <= diameter);
                assert(ecc.scc_nb[u] != ecc.scc_nb[center] || e_u >= radius);

                // Certifiers:
                trav.tree_eccentricities();
                for (int i = trav_rev.nvis()-1; i >= 0; --i) {
                    int v = trav_rev.visit(i);
                    int64_t e_v = ecc.ecc(v);
                    int64_t lb = trav_rev.dist(v);
                    int64_t ub = trav_rev.dist(v) + e_u;
                    bool in_scc_of_u = trav.visited(v) && trav_rev.visited(v);
                    if (lb == e_v) { lb_cert_of[u].push_back(v); }
                    if (ub == e_v && in_scc_of_u)
                        { ub_cert_of[u].push_back(v); }
                    if (n <= r_lim
                        || (atrav.dist(v) < rad && btrav.dist(v) < rad)) {
                        if (lb >= rad)
                            { rad_cert_of[u].push_back(v); ++m_lb_cert; }
                    }
                    if ((n <= d_lim && m_diam_cert <= m_dcert_lim && ! low_cert)
                        || ctrav.dist(v) > diam-rad || ! in_c_scc(u)) {
                       if (ub <= diam && in_scc_of_u)
                           { diam_cert_of[u].push_back(v); ++m_diam_cert; }
                    }
                    lb = in_scc_of_u ? trav_rev.dist(v)
                        : std::max(trav_rev.dist(v), e_u - trav.dist(v));
                    ub = trav.tree_ecc(v);
                    if (lb == e_v) { lb_ext_of[u].push_back(v); ++m_lb_ext; }
                    if (ub == e_v  && in_scc_of_u
                        && (! low_cert || n_ub_ext[v] <= 1)) {
                        ub_ext_of[u].push_back(v);
                        mutex.lock();
                        ++(n_ub_ext[v]); // mem race here
                        mutex.unlock();
                        ++m_ub_ext;
                    }
                }
            
                // trav_order(last_bfs); ** all zero will work too : f is last 

                // ------------ coball : complementary of B(u, rad)

                for (int i = nvis-1; i >= 0; --i) {
                    int v = trav.visit(i);
                    if (trav.dist(v) <= rad + beta_hyp) {
                        rad_coball_size[u] = n-1 - i;
                        break;
                    }
                    if (do_coballs) rad_coball[u].push_back(v);
                }

                std::sort(rad_coball[u].begin(), rad_coball[u].end());
                rad_coball_dist[u] = std::vector<int64_t>();
                for (int v : rad_coball[u]) {
                    rad_coball_dist[u].push_back(trav.dist(v));
                }
          
                // ------------ ecc. lower bound certificates:

                // Max for some order (init):
                for (auto & entry : f_max) entry.second = f;
                is_last_bfs[f].fetch_or(true, std::memory_order_relaxed);
            
                // Scan all of them:
                for (int i = nvis-1; i >= 0; --i) {
                    int v = trav.visit(i);
                    if (trav.dist(v) < e_u) { // not lb cert any more
                        n_furthest[u] = nvis-1 - i;
                        break;
                    }
                    //if (e_u <= rad + beta_hyp)
                    is_furthest[v].fetch_or(true, std::memory_order_relaxed);
                    for (auto & entry : f_max) {
                        const auto & ord = order[entry.first];
                        if (ord[v] > ord[entry.second]) entry.second = v;
                    }
                }

                // Max for some order:
                if (e_u <= rad + beta_hyp) {
                    for (const auto & entry : f_max) {
                        is_last[entry.first][entry.second]
                            .fetch_or(true, std::memory_order_relaxed);
                    }
                }

            
                // ------------ diam upper bound certificates (at max distance):

                int64_t diam = ecc.diam_lb;
                
                // last cert:
                int c = u, i_c = 0;
                for (int i = nvis-1; i >= 0; --i) {
                    int v = trav.visit(i);
                    if (trav.dist(v) + ecc.ecc(v) <= diam) {
                        c = v; i_c = i;
                        break;
                    }
                }

                if(low_ecc[u])
                    is_last_low_ecc[f].fetch_or(true,
                                                std::memory_order_relaxed);

                if (c == u) one_cert[u].fetch_or(true,
                                                 std::memory_order_relaxed);

                // Max for some order (init):
                for (auto & entry : f_max) entry.second = c;
            
                // Scan all of them:
                for (int i = i_c; i >= 0; --i) {
                    int v = trav.visit(i);
                    if (trav.dist(v) < trav.dist(c)) break; // not far any more
                    if (trav.dist(v) + ecc.ecc(v) <= diam) { // diam cert
                        is_far_cert[v].fetch_or(true,
                                                std::memory_order_relaxed);
                        for (auto & entry : f_max) {
                            const auto & ord = order[entry.first];
                            if (ord[v] > ord[entry.second]) entry.second = v;
                        }
                    }
                }

                // Max for some order:
                for (const auto & entry : f_max) {
                    is_cert[entry.first][entry.second]
                        .fetch_or(true, std::memory_order_relaxed);
                }

            
                // -------- log progress

                if (i_thd == 0 && verb::progress()) {
                    count();
                    verb::cerr() << "m_lb_ext:"<< m_lb_ext
                                 << " m_ub_ext:"<< m_ub_ext
                                 << " m_lb_cert:"<< m_lb_cert
                                 << " m_diam_cert:"<< m_diam_cert <<"\n";
                    verb::cerr("BFSs") << "todo: " << (n-u) <<" / "<< n <<"\n";
                }
              }
            }

            mutex.lock();
            m_dcert_max = std::max(m_dcert_max, m_diam_cert);
            mutex.unlock();
                
            };// thread_for

            
          auto thread_antipode = [&](int i_thd, int n_thd) { 
                traversal<graph> trav(n), trav_f(n);

            for (int u = 0; u < n; ++u) {
                if (u % n_thd == i_thd && (do_coballs || is_furthest[u])) {

                trav.clear();
                if (weighted) trav.dijkstra(g, u); else trav.bfs(g, u);

                int f = trav.last_visited();
                int64_t e_u = trav.dist(f);
                assert(e_u == ecc.ecc(u));

                  // antipode quadruples with a = f:
                int ecc_i_f = ecc.dist_lab_index(f); // dist to f in ecc.dlb ?
                if (ecc_i_f < 0) {
                    trav_f.clear();
                    if (weighted) trav_f.dijkstra(g, f); else trav_f.bfs(g, f);
                    assert(trav.nvis() == n);
                }
                // ------------ antipode-hyperbolicity:

                // Scan v,w,u,f antipode quadruples (vw, wu > R+beta, uf=e_u)
                //   check vu <= R || vf <= R

                // Easier : long triangles: max uv s.t. vf > R (+ beta)
                int64_t delta = e_u - rad, dmax = -rad;
                if (do_quad_antipode) {
                    for (int v : rad_coball[u]) {
                        int64_t d_vf = ecc_i_f < 0 ? trav_f.dist(v)
                                                 : ecc.dist_lab(v, ecc_i_f);
                        if (d_vf > std::max(rad + beta_hyp, dmax))
                            dmax = d_vf - rad;
                    }
                }
                dist_far_antipode[u] = dmax;

                // Midle : R+1 triangle: uvw s.t. uv, vw, wu > rad

                // Harder :
                // for all v in coball(f) cap coball(u),
                // calc max_{ w in coball(a) cap coball(v) } min(vw,uw)
                dmax = beta_hyp;
                if (do_quad_antipode && e_u > rad + beta_hyp) {
                    for (int iw = rad_coball[u].size() - 1; iw >= 0; --iw) {
                        int w = rad_coball[u][iw];
                        for (int iv = rad_coball[w].size() - 1; iv >= 0; --iv) {
                            int v = rad_coball[w][iv];
                            int64_t d_vu = trav.dist(v);
                            if (d_vu <= rad) {
                                int64_t d_vf = ecc_i_f < 0 ? trav_f.dist(v)
                                                   : ecc.dist_lab(v, ecc_i_f);
                                if (d_vf <= rad) { // violating quadruple
                                    int64_t dm =
                                        std::min(rad_coball_dist[u][iw],
                                                 rad_coball_dist[w][iv]) - rad;
                                    if (dm > dmax) dmax = dm;
                                }
                            } else { //R+1-triangle
                                ++(rad_triangle[u]);
                            }
                        }
                    }
                }
                quad_antipode[u] = dmax;

                
                // eccentricity of nodes at distance <= R from u,f
                int n_inter_balls = 0;
                if (is_furthest[u]) {
                    int64_t rad_dip = rad;
                    for (int i = 0; i < trav.nvis(); ++i) {
                        int v = trav.visit(i);
                        if (trav.dist(v) > rad) break; // else v in B(u,rad)
                        int64_t d_vf = ecc_i_f < 0 ? trav_f.dist(v)
                                                   : ecc.dist_lab(v, ecc_i_f);
                        int64_t r = std::max(trav.dist(v), d_vf);
                        if (r < rad_dip) rad_dip = r;
                    }
                    int dmax = 0;
                    for (int i = 0; i < trav.nvis(); ++i) {
                        int v = trav.visit(i);
                        if (trav.dist(v) > rad_dip) break; // else v in B(u,rad)
                        int64_t d_vf = ecc_i_f < 0 ? trav_f.dist(v)
                                                   : ecc.dist_lab(v, ecc_i_f);
                        if (d_vf <= rad_dip) {
                            ++n_inter_balls;
                            int64_t e_v = ecc.ecc(v);
                            if (e_v > dmax) dmax = e_v;
                        }
                    }
                    n_inter_dipole[u] = n_inter_balls;
                    ecc_max_dipole[u] = dmax;
                }
                

                if (i_thd == 0 && verb::progress()) {
                    verb::cerr("antipodes")
                        << "todo: " << (n-u) <<" / "<< n <<"\n";
                }

              }
            }

          };// thread_antipode
              
            std::vector<std::thread> threads(n_thread);

            if (n_thread <= 1) thread_for(0, 1);
            else {
                for (int i = 0; i < n_thread; ++i) {
                    threads[i] = std::thread(thread_for, i, n_thread);
                }
                std::cerr << n_thread <<" threads\n";
                for (int i = 0; i < n_thread; ++i) threads[i].join();
            }
            
            verb::lap("all_bfs_BFSs");

            if (do_antipodes) {
                if (n_thread <= 1) thread_antipode(0, 1);
                else {
                    for (int i = 0; i < n_thread; ++i) {
                        threads[i] = std::thread(thread_antipode, i, n_thread);
                    }
                    std::cerr << n_thread <<" threads\n";
                    for (int i = 0; i < n_thread; ++i) threads[i].join();
                }

                verb::lap("quad_antipode");
            }

            // ------------ nb lb cert per node ---------------

            set_distr_col(n_furthest, "n_furthest");
            set_distr_col(n_visit, "n_visit");
            
            // ------------ last ub/lb cert for various orders ---------------

            count();
            for (const auto & entry : n_last) {
                set_col(entry.first + "_lst", entry.second);
            }
            for (const auto & entry : n_cert) {
                set_col(entry.first + "_crt", entry.second);
            }
            set_col("nb_furthest", nb_furthest);
            set_col("nb_far_cert", nb_far_cert);
            set_col("nb_one_cert", nb_one_cert);

            int n_last_low_ecc = 0;
            for (int u = 0; u < n; ++u) {
                if(is_last_low_ecc[u]) ++n_last_low_ecc;
            }
            std::cerr << "n_last_low_ecc: " << n_last_low_ecc <<"\n";
            set_col("n_last_low_ecc", n_last_low_ecc);
            
            std::vector<int64_t> ecc_furthest;
            for (int u = 0; u < n; ++u) {
                if (is_furthest[u]) ecc_furthest.push_back(ecc.ecc(u));
            }

            set_distr_col(ecc_furthest, "ecc_furthest");
            
            std::vector<int64_t> ecc_idx_lst;
            for (int u = 0; u < n; ++u) {
                if (is_last["max_idx"][u]) ecc_idx_lst.push_back(ecc.ecc(u));
            }

            set_distr_col(ecc_idx_lst, "ecc_max_idx_lst");
            
            // ------------ max distance to rad-far from antipode ----

            set_distr_col(quad_antipode, "quad_antipode");

            std::cerr <<"quad-hyperbol: "<< quad_antipode[n-1]
                      <<" ,  rad:"<< rad <<"\n";
            
            set_distr_col(rad_triangle, "rad_triangle");

            std::cerr <<"rad-triangle: "<< rad_triangle[n-1]
                      <<" + rad:"<< rad <<" = "<< rad_triangle[n-1] + rad
                      <<"\n";
            
            set_distr_col(dist_far_antipode, "dist_far_antipode");
            
            std::cerr <<"bi-hyperbol: "<< dist_far_antipode[n-1]
                      <<" + rad:"<< rad <<" = "<< dist_far_antipode[n-1] + rad
                      <<"\n";

            set_distr_col(rad_coball_size, "rad_coball_size");

            std::cerr <<"rad_coball_size: "
                      << ((float) mean1000(rad_coball_size)) / 1000.0
                      <<"  max:"<< rad_coball_size[n-1]
                      <<"\n";

            std::vector<int64_t> ecc_dip, n_dip;
            for (int u = 0; u < n; ++u) {
                if (is_furthest[u]) {
                    ecc_dip.push_back(ecc_max_dipole[u]);
                    n_dip.push_back(n_inter_dipole[u]);
                }
            }
            set_distr_col(ecc_dip, "ecc_max_dipole");
            set_distr_col(n_dip, "n_inter_dipole");

            std::cerr <<"ecc_max_dipole: "<< ecc_dip[ecc_dip.size() - 1]
                      <<" n_inter_dipole: " << n_dip[n_dip.size() - 1]
                      <<" ,  +rad:"<< rad <<"\n";
            
            // ------------ certifiers --------------------

            int one_ub_ext = 0;
            for (int u = 0; u < n; ++u) {
                if (n_ub_ext[u] == 1) ++one_ub_ext;
            }
            std::cerr << "one certif_ext: " << one_ub_ext << std::endl;
            set_col("one_ub_ext", one_ub_ext);

            int P = ecc.gdy_set_cov(lb_cert_of, n).size();
            std::cerr <<"P done\n";
            int C = ecc.gdy_set_cov(ub_cert_of, n).size();
            std::cerr <<"C done\n";
            int Pext = directed ? -1 : ecc.gdy_set_cov(lb_ext_of, n).size();
            int Cext = -1;
            if ( ! low_cert) Cext = ecc.gdy_set_cov(ub_ext_of, n).size();

            std::cerr <<"P, Pext, C, Cext:  " << P <<", "<< Pext
                      <<", "<< C <<", "<< Cext << std::endl;

            set_col("P", P);
            set_col("Pext", Pext);
            set_col("C", C);
            set_col("Cext", Cext);
        
            int P_inf = ecc.gdy_packing(lb_cert_of, n).size();
            int C_inf = ecc.gdy_packing(ub_cert_of, n).size();
            int Pext_inf = ecc.gdy_packing(lb_ext_of, n).size();
            int Cext_inf = -1;
            if ( ! low_cert) Cext_inf = ecc.gdy_packing(ub_ext_of, n).size();

            std::cerr <<"P_inf, Pext_inf, C_inf, Cext_inf:  "
                      << P_inf <<", "<< Pext_inf
                      <<", "<< C_inf <<", "<< Cext_inf << std::endl;

            set_col("P_inf", P_inf);
            set_col("Pext_inf", Pext_inf);
            set_col("C_inf", C_inf);
            set_col("Cext_inf", Cext_inf);

            std::vector<int> tocov = {-1}; // default for gdycov/pack
            if (n > r_lim) {
                tocov.clear();
                for (int u = 0; u < n; ++u)
                    if (atrav.visited(u) && btrav.visited(u)
                        && atrav.dist(u) < rad && btrav.dist(u) < rad)
                        tocov.push_back(u);
            }
            
            int R_inf = ecc.gdy_packing(rad_cert_of, n, tocov).size();
            R_inf = std::max(2, R_inf);
            int R = (n <= r_lim ? 0 : 2)
                + ecc.gdy_set_cov(rad_cert_of, n, tocov).size();
            std::cerr <<"R_inf, R:   " << R_inf <<", "<< R << "\n";

            tocov.clear();
            for (int u = 0; u < n; ++u) {
                if (ecc.ecc(u) > rad + beta_hyp) tocov.push_back(u);
            }
            
            int Rbeta_inf = -1;
            if (do_coballs) ecc.gdy_packing(rad_coball, n, tocov).size();
            int Rbeta = -1;
            if (do_coballs ) ecc.gdy_set_cov(rad_coball, n, tocov).size();
            std::cerr <<"Rbeta_inf, Rbeta:   " << Rbeta_inf
                      <<", "<< Rbeta << "\n";

            tocov = {-1}; // default for gdycov/pack
            if (n > d_lim || m_dcert_max > m_dcert_lim || low_cert) {
                tocov = ctr_uncov;
            }
                
            int D_inf = ecc.gdy_packing(diam_cert_of, n, tocov).size();
            D_inf = std::max(1, D_inf);
            unsigned long D =
                ((n > d_lim || m_dcert_max > m_dcert_lim || low_cert) ? 1 : 0)
                + ecc.gdy_set_cov(diam_cert_of, n, tocov).size();
            verb::cerr() << "m_dcert max,lim : "
                         << m_dcert_max <<","<< m_dcert_lim << "\n";
            if (n <= d_lim && m_dcert_max <= m_dcert_lim && ! low_cert)
                D = std::min(D, 1 + ecc.gdy_set_cov(diam_cert_of, n,
                                                    ctr_uncov).size());
            std::cerr <<"D_inf, D:   "<< D_inf <<", "<< D << "\n";

            set_col("m_dcert_max", m_dcert_max);
            set_col("R_inf", R_inf);
            set_col("R", R);
            set_col("Rbeta_inf", Rbeta_inf);
            set_col("Rbeta", Rbeta);
            set_col("D_inf", D_inf);
            set_col("D", D);
            
       } // small graph

        verb::lap("all_bfs");
    }
    
    // ------------------- Optimize a given certificate -----------
    if (optim_certif != "" && ! all_bfs) {

        if(g_scc_nb != 1) bye("The graph is not strongly connected.");
        assert( ! directed); // not checked                
        
        std::vector<int> certif;

        if (optim_certif == "-") { in = stdin; }
        else {
            std::cerr << "open '" << optim_certif << "' to read certif\n";
            in = fopen(optim_certif.c_str(), "r");
        }
        char u[1024];
        for ( ; fscanf(in, " %s \n", u) >= 1 ; ) {
            certif.push_back(vi[u]-1);
        }
        std::cerr << "certif size: " << certif.size() << std::endl;
        for (int u : certif) std::cerr << u <<" ";
        std::cerr << std::endl;

        // We need radius and diameer.
        //eccentricity<graph> ecc(g, g_rev, directed, weighted);
        int rad = ecc.radius(), diam = ecc.diameter_sample();
        std::cerr << "R = " << rad <<" (cert";
        for (int p : ecc.P) std::cerr <<" "<< p; 
        std::cerr << ", opt: "
                  << ecc.rad_certif.size()
                  <<") center: "<< ecc.rad_node
                  <<"  D = " << diam <<" (cert";
        for (int c : ecc.C) std::cerr <<" "<< c;
        std::cerr << ", opt: "<< ecc.diam_certif.size()
                  <<") periph: "<< ecc.diam_node << "\n";

        set_col("ecc_rad_true", ecc.optim_lb_certif(ecc.P, rad, true).size());
        set_col("ecc_diam_true", ecc.optim_ub_certif(ecc.C, diam, true).size());
        set_col("ecc_rad_opt", ecc.rad_certif.size());
        set_col("ecc_diam_opt", ecc.diam_certif.size());
        
        
        ecc.clear();
        for (int u : certif) {
            ecc.sweepImproveBounds(u, true, true, rad, diam, true, true);
        }
        for (int u : g)
            if (ecc.ecc_lb(u) < rad || ecc.ecc_ub(u) > diam) {
                std::cerr << u << " not certified : "
                          << ecc.ecc_lb(u) <<" <= e_u <= "<< ecc.ecc_ub(u)
                          <<std::endl;
                bye("not a certificate!");
            }
        
        std::vector<int> rad_cert = ecc.optim_lb_certif(certif, rad, true);
        std::vector<int> diam_cert = ecc.optim_ub_certif(certif, diam, true);

        std::cerr << "rad_cert: " << rad_cert.size()
                  <<" diam_cert: "<< diam_cert.size() << std::endl;

        set_col("cert", certif.size());
        set_col("rad_cert_opt", rad_cert.size());
        set_col("diam_cert_opt", diam_cert.size());

        try {
            rad_cert = ecc.optim_lb_certif(certif, rad, false);
        } catch(const std::string& s) {
            std::cerr << s << std::endl;
            rad_cert = {};
        }
        try {
            diam_cert = ecc.optim_ub_certif(certif, diam, false);
        } catch(const std::string& s) {
            std::cerr << s << std::endl;
            diam_cert = {};
        }
        
        set_col("rad_cert__false", rad_cert.size());
        set_col("diam_cert_false", diam_cert.size());
        
        verb::lap("optim_certif");
    }

    
    // ------------------- Print columns ------------------------------
    
    set_col("n", n, 1);
    set_col("m", m, 1);
    set_col("dir", directed, 1);
    set_col("wgt", weighted, 1);
    set_col("nb_scc", g_scc_nb, 2);


    if (cols_verb > 0) { print_cols("end"); }
    verb::end();
    exit(0);
    
    /*** ------------------------- Smooth corners ----------------------
    {
        std::vector<int> corners(furthest_central2);
        std::vector<std::vector<int>> dist(n);
        traversal<graph> trav(n);
        for (int i = 0; i < corners.size() ; ++i) {
            int u = corners[i];
            if (weighted) trav.dijkstra(g, u);
            else trav.bfs(g, u);
            dist[u].reserve(n);
            for (int v : g) dist[u][v] = trav.dist(v);
            trav.clear();
        }

        auto dist_corners =
            [&dist,&corners,n](int v) -> int {
            int d = n;
            for (int f : corners)
                if (dist[f][v] < d) d = dist[f][v];
            return d;
        };

        for (int u : g) {
            int e = ecc[u];
            if (e > rad) {
                for (int v : furthest[u]) {
                    int d = dist_corners(v);
                    double a = ((double) d) / ((double) (e - rad));
                    if (a > 0.4) { // add v to corners
                        corners.push_back(v);
                        if (weighted) trav.dijkstra(g, v);
                        else trav.bfs(g, v);
                        dist[v].reserve(n);
                        for (int x : g) dist[v][x] = trav.dist(x);
                        trav.clear();
                    }
                }
            }
        }
        
        double alpha = 0.;
        
        for (int u : g) {
            int e = ecc[u], dmax = 0;
            for (int v : furthest[u]) {
                int d = dist_corners(v);
                if (d > dmax) dmax = d;
            }
            if (e > rad) {
                double a = ((double) dmax) / ((double) (e - rad));
                if (a > alpha) alpha = a;
            }
        }
        std::cerr << "alpha = " << alpha
                  << " with " << corners.size() << " corners\n";

        std::cout <<" "<< corners.size() <<" "<< alpha; 
    }
    verb::lap("smooth corners");
    ***/

}

#ifndef GRAPHGEN_HH
#define GRAPHGEN_HH

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <utility>
#include <random>

#include "mgraph.hh"
#include "dyn_graph.hh"

/**
 * Generate graphs
 */

typedef dyn_graph<int, -1> dgraph;

int add_path(int fresh_nb, dgraph &g, int from, int to, int len) {
    assert(len > 0);
    if (len == 1) {
        g.add_edge(from, to);
        return fresh_nb;
    }
    int u = fresh_nb++;
    g.add_edge(from, u);
    --len;
    while (len > 1) {
        int v = fresh_nb++;
        g.add_edge(u, v);
        --len;
        u = v;
    }
    g.add_edge(u, to);
    return fresh_nb;
}

dgraph path_graph(int len) {
    dgraph g;
    int from = 0, to = len, fresh = 1;
    add_path(fresh, g, from, to, len);
    return g;
}

dgraph cycle_graph(int len) {
    assert(len > 1);
    dgraph g = path_graph(len - 1);
    g.add_edge(len - 1, 0);
    return g;
}

dgraph product(const dgraph &gx, const dgraph &gy,
               std::function<bool(int, int)>
                 filtr = [](int x, int y) { return true; }
               ) {
    dgraph h;
    int nx = gx.n();
    auto num = [nx](int x, int y) { return nx * y + x; };
    for (int x : gx) {
        for (int y : gy) {
            for (int x2 : gx[x])
                if (filtr(x,y) && filtr(x2,y))
                    h.add_edge(num(x,y), num(x2,y));
            for (int y2 : gy[y])
                if (filtr(x,y) && filtr(x,y2))
                    h.add_edge(num(x,y), num(x,y2));
        }
    }
    return h;
}

dgraph grid_graph(int nx, int ny) {
    return product(path_graph(nx), path_graph(ny));
}

dgraph grid_ball_graph(int nx, int ny) {
    auto filtr = [nx, ny](int x, int y) {
        x -= nx / 2; y -= ny / 2;
        return std::max(x, -x) + std::max(y, -y) <= nx/2 + ny/2;
    };
    return product(path_graph(nx), path_graph(ny), filtr);
}

int skel_prefix_min(int len, int iter = 1000) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> unif(0.0, 1.0);
    int64_t s = 0;
    for (int i = 0; i < iter; ++i) {
        int j = 0;
        float p = unif(gen);
        ++s;
        while (j <= len) {
            ++s;
            float dj = std::log(1. - unif(gen)) / std::log(1. - p);
            if (dj + j >= INT_MAX) break;
            j += (int) dj;
            p = p * unif(gen);
        }
    }
    return (int) s / iter;
}

dgraph erdos_renyi_random_graph(int n, float p) {
    std::random_device rand_dev;
    std::mt19937 rand_gen(rand_dev());
    std::uniform_real_distribution<> unif(0.0, 1.0);
    // generate graph :
    dgraph g;
    for (int u=0; u<n; ++u) {
         for (int v=0; v<n; ++v) {
              if (unif(rand_gen) < p) {
                   g.add_edge(u, v);
              }
         }
    }
    return g;
}

dgraph power_law_random_graph(int n, float beta) {
    // power law with parameter beta :
    std::vector<double> p_deg(n);
    double nf = n;
    for (int d=1; d<n; ++d) p_deg[d] = nf / pow((float)d, beta);
    p_deg[0.] = 0.;
    std::discrete_distribution<> power_law(p_deg.begin(), p_deg.end());
    // generate degrees :
    std::random_device rand_dev;
    std::mt19937 rand_gen(rand_dev());
    std::vector<int> deg(n);
    int m = 0;
    for(int u=0; u<n; ++u) {
        deg[u] = power_law(rand_gen);
        m += deg[u];
    }
    // generate edges :
    std::vector<int> dst, src;
    dst.reserve(m); src.reserve(m);
    for (int u=0; u<n; ++u) {
        for (int e=0; e<deg[u]; ++e) {
            dst.push_back(u); src.push_back(u);
        }
    }
    std::uniform_int_distribution<> rnd_int(0, m-1);
    for (int e=m-1; e>0; --e) {
        int f = rnd_int(rand_gen) % (e+1);
        std::swap(dst[e], dst[f]);
    }
    // avoid self edge :
    for (int u=0, e=0; u<n; ++u) {
        for (int i=0; i<deg[u]; ++i, ++e) {
            while (u == dst[e]) {
                int f = rnd_int(rand_gen) % m;
                while (u == dst[f] || src[f] == dst[e]) {
                    f = rnd_int(rand_gen) % m;
                }
                std::swap(dst[e], dst[f]);
            }
        }
    }
    // generate graph :
    dgraph g;
    for (int u=0, e=0; u<n; ++u) {
        for (int i=0; i<deg[u]; ++i, ++e) g.add_edge(u, dst[e]);
    }
    return g;
}

typedef mgraph<int> wgraph;

std::vector<wgraph::edge> unit_disk_graph(int n, double deg,
                                          int max_weight = 1000) {
    // pi r^2 n = deg
    double r = sqrt(deg / (M_PI * n));
    std::vector<double> x(n), y(n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    auto dist2_center = [](double x, double y) {
        x -= 0.5; y -= 0.5;
        return (x * x + y * y);
    };
    for (int i = 0; i < n; ++i) {
        x[i] = dis(gen);
        y[i] = dis(gen);
        while (dist2_center(x[i], y[i]) >= 1) {
            x[i] = dis(gen);
            y[i] = dis(gen);
        }
    }
    // square rxr cells:
    int row_len = (int)(1/r) + 1;
    auto cell_index = [r, row_len](double x, double y) {
        int i = (int)(x / r);
        int j = (int)(y / r);
        return i + j * row_len;
    };
    std::vector<std::vector<int>> cell(row_len * row_len);
    for (int i = 0; i < n; ++i) {
        int c = cell_index(x[i], y[i]);
        cell[c].push_back(i);
    }
    // graph: r-disk is included in 9 neighboring cells
    std::vector<wgraph::edge> edges;
    auto dist = [&x, &y](int u, int v) {
        double dx = x[u] - x[v], dy = y[u] - y[v];
        return sqrt(dx * dx + dy * dy);
    };
    for (int u = 0; u < n; ++u) {
        int c = cell_index(x[u], y[u]);
        for (int i=-1; i<=1; ++i) {
            for (int j=-1; j<=1; ++j) {
                int d = c + i + j * row_len; // neighboring cell
                if (d >= 0 && d < row_len * row_len) {
                    for (int v : cell[d]) {
                        double d_uv = dist(u, v);
                        if (d_uv < r && u < v) {
                            int w_uv = (int)(max_weight * d_uv / r);
                            edges.push_back(wgraph::edge(u, v, w_uv));
                            edges.push_back(wgraph::edge(v, u, w_uv));
                        }
                    }
                }
            }
        }
    }
    return edges;
}


/** Bow tie :
 *             ...p
 *              2
 * .           q-1           .  
 * .  q 1  q+1  0  q+1  3 q  .
 * .                         .
 * p           q+2           p
 * 
 *              4
 */

dgraph bow_tie_graph(int p, int q) {
    dgraph g;
    int n = 5, ni = -1;
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 0 5
    n = add_path(n, g, 1, 0, q+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 1 55
    n = add_path(n, g, 0, 3, q+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 2 105
    n = add_path(n, g, 2, 0, q-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 3 153
    n = add_path(n, g, 2, 1, q-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 4 201
    n = add_path(n, g, 2, 3, q-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 5 249
    n = add_path(n, g, 4, 0, q+2);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 6 300
    n = add_path(n, g, 4, 1, q+2);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 7 351
    n = add_path(n, g, 4, 3, q+2);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 8 402

    // max degree node : 0
    for (int i = 1; i <= p; ++i) g.add_edge(0, n++);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 9 452

    // add p nodes attached to 2
    int u = n; n += p;
    std::cerr <<"top:";
    for (int v = u; v < n; ++v) {
        std::cerr <<" "<< v;
        g.add_edge(v, 2);
    }
    std::cerr <<"\n";
    // add_path attached to 1
    u = n++;
    int w = n++;
    n = add_path(n, g, u, w, p-1);
    std::cerr <<"left:";
    for (int v = u, m = n; v < m; ++v) {
        //g.add_edge(v, 1);
        n = add_path(n, g, v, 1, 2);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, q-2);
    }
    std::cerr <<"\n";
    // add_path attached to 3
    u = n++;
    w = n++;
    n = add_path(n, g, u, w, p-1);
    std::cerr <<"right:";
    for (int v = u, m = n; v < m; ++v) {
        //g.add_edge(v, 3);
        n = add_path(n, g, v, 3, 2);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, q-2);
    }
    std::cerr <<"\n";
    return g;
}

/** Losange1 :
 *       2
 *       6
 * 1     5      3
 *
 *       4
 */

dgraph losange_graph(int x) {
    dgraph g;
    int n = 10 * x * x, ni = -1;
    int y = 0;
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 0 26010
    n = add_path(n, g, 1, 2, x-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 1 26059
    n = add_path(n, g, 2, 3, x-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 2 26108
    n = add_path(n, g, 3, 5, x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 3 26159
    n = add_path(n, g, 5, 1, x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 4 26210
    n = add_path(n, g, 1, 4, x+2 + y);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 5 26262
    n = add_path(n, g, 4, 3, x+2 + y);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 6 26314
    n = add_path(n, g, 4, 5, x+2 + y);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 7 26366
    n = add_path(n, g, 5, 6, x-2);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 8 26414
    g.add_edge(6,2);

    for (int i = 1; i <= x+2; ++i) g.add_edge(5, n++);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 9 26467
    
    //g.add_edge(4,7);
    //g.add_edge(7,8);
    // add_path attached to 2
    int u = n++;
    int w = n++;
    //n = add_path(n, g, u, w, x-1);
    n += x-2;
    std::cerr <<"top:";
    for (int v = u; v < n; ++v) {
        std::cerr <<" "<< v;
        g.add_edge(v, 2);
    }
    std::cerr <<"\n";
    // add_path attached to 1
    u = n++;
    w = n++;
    n = add_path(n, g, u, w, x-1);
    int m = n;
    std::cerr <<"left:";
    for (int v = u; v < m; ++v) {
        //g.add_edge(v, 1);
        n = add_path(n, g, v, 1, 2);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, x-2 + y);
    }
    std::cerr <<"\n";
    // add_path attached to 3
    u = n++;
    w = n++;
    n = add_path(n, g, u, w, x-1);
    m = n;
    std::cerr <<"right:";
    for (int v = u; v < m; ++v) {
        //g.add_edge(v, 3);
        n = add_path(n, g, v, 3, 2);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, x-2 + y);
    }
    std::cerr <<"\n";
    return g;
}


/** Losange ok rad :
 *       2
 *       6
 * 1     5      3
 *
 *       4
 */

dgraph losange_graph_ok_rad(int x) {
    dgraph g;
    int n = 10 * x * x, ni = -1;
    int y = 0;
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 0 26010
    n = add_path(n, g, 1, 2, x-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 1 26059
    n = add_path(n, g, 2, 3, x-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 2 26108
    n = add_path(n, g, 3, 5, x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 3 26159
    n = add_path(n, g, 5, 1, x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 4 26210
    n = add_path(n, g, 1, 4, x+2 + y);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 5 26262
    n = add_path(n, g, 4, 3, x+2 + y);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 6 26314
    n = add_path(n, g, 4, 5, x+2 + y);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 7 26366
    n = add_path(n, g, 5, 6, x-2);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 8 26414
    g.add_edge(6,2);

    for (int i = 1; i <= x+2; ++i) g.add_edge(5, n++);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 9 26467
    
    //g.add_edge(4,7);
    //g.add_edge(7,8);
    // add_path attached to 2
    int u = n++;
    int w = n++;
    //n = add_path(n, g, u, w, x-1);
    n += x-2;
    std::cerr <<"top:";
    for (int v = u; v < n; ++v) {
        std::cerr <<" "<< v;
        g.add_edge(v, 2);
    }
    std::cerr <<"\n";
    // add_path attached to 1
    u = n++;
    w = n++;
    n = add_path(n, g, u, w, x-1);
    int m = n;
    std::cerr <<"left:";
    for (int v = u; v < m; ++v) {
        g.add_edge(v, 1);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, x-1 + y);
    }
    std::cerr <<"\n";
    // add_path attached to 3
    u = n++;
    w = n++;
    n = add_path(n, g, u, w, x-1);
    m = n;
    std::cerr <<"right:";
    for (int v = u; v < m; ++v) {
        g.add_edge(v, 3);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, x-1 + y);
    }
    std::cerr <<"\n";
    return g;
}

/** Losange hard diam :
 *    . . . . 
 *     5 4 6
 *       0
 *.             .
 *.  1         3.
 *.             .
 *       2
 */

dgraph losange_graph_diam(int x) {
    dgraph g;
    int n = 7, ni = -1;
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 0 7
    n = add_path(n, g, 1, 5, x-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 1
    n = add_path(n, g, 6, 3, x-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 2
    n = add_path(n, g, 3, 2, x-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 3
    n = add_path(n, g, 2, 1, x-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 4
    n = add_path(n, g, 2, 0, 2*x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 5
    g.add_edge(0, 4);

    // highest degree node : 0
    for (int i = 1; i <= x+5; ++i) g.add_edge(0, n++);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 9
    
    //g.add_edge(4,7);
    //g.add_edge(7,8);
    // add_path attached to 2
    int u = n++;
    int w = n++;
    //n = add_path(n, g, u, w, x-1);
    n += x-2;
    std::cerr <<"top:";
    for (int v = u; v < n; ++v) {
        std::cerr <<" "<< v;
        g.add_edge(v, 4); g.add_edge(v, 5); g.add_edge(v, 6);
    }
    std::cerr <<"\n";
    // add_path attached to 1
    u = n++;
    w = n++;
    n = add_path(n, g, u, w, x-1);
    //n += x-2;
    int m = n;
    std::cerr <<"left:";
    for (int v = u; v < m; ++v) {
        g.add_edge(v, 1);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, x-1);
    }
    std::cerr <<"\n";
    // add_path attached to 3
    u = n++;
    w = n++;
    n = add_path(n, g, u, w, x-1);
    //n += x-2;
    m = n;
    std::cerr <<"right:";
    for (int v = u; v < m; ++v) {
        g.add_edge(v, 3);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, x-1);
    }
    std::cerr <<"\n";
    return g;
}

/** Losange_hard_diam1 :
 *     8   9
 *      6 7
 *       2
 *       
 * 1     5      3
 *
 *       4
 */

dgraph losange_graph_diam1(int x) {
    dgraph g;
    int n = 10 * x * x, ni = -1;
    int y = 0;
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 0 26010
    n = add_path(n, g, 1, 8, x-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 1 26059
    n = add_path(n, g, 9, 3, x-1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 2 26108
    n = add_path(n, g, 3, 5, x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 3 26159
    n = add_path(n, g, 5, 1, x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 4 26210
    n = add_path(n, g, 1, 4, x-1 + y);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 5 26262
    n = add_path(n, g, 4, 3, x-1 + y);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 6 26314
    n = add_path(n, g, 4, 5, x+4 + y);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 7 26366
    n = add_path(n, g, 5, 2, x-3);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 8 26414
    g.add_edge(2, 6); g.add_edge(2, 7);
    g.add_edge(6, 8); g.add_edge(7, 9);

    // let 2 (not 5) be the node with highest degree :
    for (int i = 1; i <= x+2; ++i) g.add_edge(2, n++);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 9 26467

    //n = add_path(n, g, 1, 3, 2*x-1);
    
    //g.add_edge(4,7);
    //g.add_edge(7,8);
    // add_path attached to 2
    int u = n++;
    int w = n++;
    //n = add_path(n, g, u, w, x-1);
    n += x-2;
    std::cerr <<"top:";
    for (int v = u; v < n; ++v) {
        std::cerr <<" "<< v;
        g.add_edge(v, 8); g.add_edge(v, 9);
    }
    std::cerr <<"\n";
    // add_path attached to 1
    u = n++;
    w = n++;
    n = add_path(n, g, u, w, x-1);
    int m = n;
    std::cerr <<"left:";
    for (int v = u; v < m; ++v) {
        g.add_edge(v, 1);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, x-1 + y);
    }
    std::cerr <<"\n";
    // add_path attached to 3
    u = n++;
    w = n++;
    n = add_path(n, g, u, w, x-1);
    m = n;
    std::cerr <<"right:";
    for (int v = u; v < m; ++v) {
        g.add_edge(v, 3);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, x-1 + y);
    }
    std::cerr <<"\n";
    return g;
}

/** Losange2 :
 *       6
 *       7
 *       2
 *       
 * 1     5      3
 *
 *       4
 *       8
 */

dgraph losange_graph2(int x) {
    dgraph g;
    int n = 10 * x * x, ni = -1;
    int y = 0;
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 0 26010
    n = add_path(n, g, 1, 2, x-3);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 1 26059
    n = add_path(n, g, 2, 3, x-3);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 2 26108
    n = add_path(n, g, 3, 5, x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 3 26159
    n = add_path(n, g, 5, 1, x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 4 26210
    n = add_path(n, g, 1, 4, x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 5 26262
    n = add_path(n, g, 4, 3, x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 6 26314
    n = add_path(n, g, 4, 5, x+1);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 7 26366
    n = add_path(n, g, 5, 2, x-3);
    std::cerr << ++ni <<" n = "<< n <<"\n"; // 8 26414
    g.add_edge(6, 7); g.add_edge(7, 2);
    g.add_edge(4, 8);
    
    for (int i = 1; i <= x+2; ++i) g.add_edge(5, n++);
    
    //g.add_edge(4,7);
    //g.add_edge(7,8);
    // path attached to 6
    int u = n++;
    int w = n++;
    n = add_path(n, g, u, w, x-1);
    //n += x-2;
    std::cerr <<"top:";
    for (int v = u; v < n; ++v) {
        std::cerr <<" "<< v;
        g.add_edge(v, 6);
    }
    std::cerr <<"\n";
    // path attached to 1
    u = n++;
    w = n++;
    n = add_path(n, g, u, w, x-1);
    int m = n;
    std::cerr <<"left:";
    for (int v = u; v < m; ++v) {
        g.add_edge(v, 1);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, x-1 + y);
    }
    std::cerr <<"\n";
    // path attached to 3
    u = n++;
    w = n++;
    n = add_path(n, g, u, w, x-1);
    m = n;
    std::cerr <<"right:";
    for (int v = u; v < m; ++v) {
        g.add_edge(v, 3);
        int a = n++;
        std::cerr <<" "<< a;
        n = add_path(n, g, a, v, x-1 + y);
    }
    std::cerr <<"\n";
    return g;
}



#endif // GRAPHGEN_HH

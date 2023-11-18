#ifndef ECCENTRICITY_HH
#define ECCENTRICITY_HH

#include <assert.h>
#include <stdint.h>
#include <climits>
#include <vector>
#include <algorithm>      // std::min
#include <thread>
#include <mutex>
#include <atomic>

#include "mgraph.hh"
#include "traversal.hh"
#include "verbose.hh"

/**
 * Compute some or all eccentricities of a graph G=(V,E).
 * The excentricity of a node u is ecc(u) = max_{v in V} d(u,v).
 * The radius of G is R = min_{u in G} ecc(u).
 * The diameter of G is D = max_{u in G} ecc(u).
 * Example :
 *    mgraph<> g;
 *    eccentricity<graph> exc(g);
 *    int r = exc.radius();
 *    int d = exc.diameter();
 *    exc.all();
 *    int u = 1;
 *    int ecc_u = exc.excentricity(u);
 */

/* How to start diameter computation: find center first ? */
enum start_diam { STRAIGHT = 0, SUMSWEEP = 1, RADIUS = 2 };


template<typename G,  // Graph type with G::vertex = int
         typename WL = int64_t, // long weight for summing weights, 
                     // with appropriate casting for the following constants:
         int64_t max_weight = INT64_MAX, int64_t zero_weight = 0>
class eccentricity {
public:
    typedef int V;
    typedef typename G::weight W;
    typedef WL long_weight;

private:
    const int n;
    std::vector<WL> ecc_lb_, ecc_ub_, sum_, next_hop_dist, dist_C;
    std::vector<std::vector<WL>> dist_lab_lb, dist_lab_ub;
    std::vector<V> lb_node, // node p in P giving lower bound
        next_hop_to_lb, // first node on sh. path to p
        diam_c_to_lb; // furth. node c on path to p s.t. dist(c)+ecc(c)<=diam_lb
    traversal<G, WL, max_weight, zero_weight> trav;
    std::vector<WL> save_dist;
    std::vector<V> save_visit;
    std::vector<V> sample;
    std::vector<WL> sample_count;
 
public:
    const G &graph, &graph_rev;
    const bool directed, weighted;
    bool strongly_connected;
    WL rad_ub, diam_lb;
    V rad_node, diam_node, last_sweep_node;
    bool last_sweep_backward;
    int nsweep, mprune, nprune_sweep, rad_nsweep, diam_nsweep, all_ecc_nsweep;
    int rad_todo, diam_todo, all_ecc_todo;
    int last_lb_improve, last_ub_improve;
    std::vector<V> P /* pseudo-peripheral nodes giving ecc lower bounds */,
        C, /* pseudo-central nodes giving ecc upper bounds */
        Rcoballs;
    std::vector<std::vector<V>> Pcoballs;
    std::vector<V> rad_certif, diam_certif, all_lb_certif, all_ub_certif;
    std::vector<int> scc_nb; int scc_largest; V scc_node_largest;
    
    static const int not_nsweep = -1;
    
    eccentricity(const G &g, const G &g_rev,
                 bool directed = true, bool weighted = true) :
        n(g.n()),
        ecc_lb_(g.n(), zero_weight), ecc_ub_(g.n(), max_weight),
        dist_lab_lb(g.n()), dist_lab_ub(g.n()), dist_C(g.n(),max_weight),
        lb_node(n, -1), next_hop_to_lb(n, -1), diam_c_to_lb(n, -1),
        sum_(g.n(), zero_weight), next_hop_dist(g.n(), zero_weight),
        trav(g.n()), save_dist(g.n(),max_weight), save_visit(g.n()),
        graph(g), graph_rev(g_rev), directed(directed), weighted(weighted),
        rad_ub(max_weight), diam_lb(zero_weight),
        rad_node(G::not_vertex), diam_node(G::not_vertex),
        last_sweep_node(G::not_vertex), last_sweep_backward(false),
        nsweep(0), mprune(0), nprune_sweep(0),
        rad_nsweep(not_nsweep), diam_nsweep(not_nsweep),
        all_ecc_nsweep(not_nsweep),
        rad_todo(g.n()), diam_todo(g.n()), all_ecc_todo(g.n()),
        last_lb_improve(0), last_ub_improve(0),
        sample(), sample_count(g.n(), 0)
    {
        verb::cerr("eccentricity()", 1)
            << " ++++++ graph n="<< n <<" m="<< graph.m()
            << " directed: "<< directed <<" weighted: "<< weighted <<"\n";
        assert(graph.n() == graph_rev.n() && graph.m() == graph_rev.m() && n>0);
        int n_scc = trav.strongly_connected_components(graph);
        strongly_connected = n_scc <= 1;
        if (n_scc > 1) verb::cerr("eccentricity()", 1)
                           << n_scc <<" strongly connected components !\n";
        scc_nb = trav.scc_vector();
        scc_largest = trav.scc_largest();
        scc_node_largest = trav.scc_node(scc_largest);
        for (int u = 0; u < n; ++u) {
            lb_node[u] = u;
            next_hop_to_lb[u] = u;
            next_hop_dist[u] = zero_weight;
            diam_c_to_lb[u] = u;
        }
    }

    WL ecc_lb(V u) const { return ecc_lb_[u]; }
    WL ecc_lb_node(V u) const { return lb_node[u]; }
    WL ecc_ub(V u) const { return ecc_ub_[u]; }
    WL ecc(V u) const {
        if (ecc_lb_[u] != ecc_ub_[u]) {
            //throw std::invalid_argument ("this in excentricity.ecc():
            std::cerr << "bounds are not tight, run this.all()" //);
                      <<"\nu="<< u <<" lb="<< ecc_lb_[u]
                      <<" ub="<< ecc_ub_[u] <<"\n";
            assert(false);
        }
        return ecc_lb_[u];
    }
    
    V center() const {
        if (rad_todo > 0) radius();
        return rad_node;
    }
    
    int dist_lab_index(V p) {
        assert(P.size() == dist_lab_lb[p].size());
        for (int i = P.size()-1; i >= 0; --i) {
            if (P[i] == p) {
                assert(dist_lab_lb[p][i] == zero_weight);
                return i;
            }
        }
        return -1;
    }
    WL dist_lab(V u, int i) {
        return dist_lab_lb[u][i];
    }
    void make_dist_lab() {
        if (n > 0 && P.size() > dist_lab_lb[0].size()) { 
            clear_ecc_lb();
            for (int i = 0; i < P.size(); ++i) {
                V p = P[i];
                last_sweep_node = G::not_vertex;
                trav.clear();
                if (weighted) trav.dijkstra(graph_rev, p);
                else trav.bfs(graph_rev, p);
                improves_ecc_lb(max_weight, false, false, true);
            }
        }
    }

    // Smallest eccentricity of a node in largest strong. conn. comp.
    WL radius(V start = G::not_vertex,
              bool optim_diam = false, bool optim_all_ecc = false,
              int scc = -1, bool do_sum_sweep = false) {
        if (do_sum_sweep) {
            start = sum_sweep(start);
            update_todo();
        }
        if (rad_todo == 0) return rad_ub;
        if (start == G::not_vertex) start = max_degree_node();
        if (scc == -1) {
            scc = (start != G::not_vertex) ? scc_nb[start] : scc_largest;
        }
        
        V c = start, p = G::not_vertex;
        while (rad_todo > 0) {
            
            // Sweep from pseudo-center : node with minimal ecc_lb_
            for (V u = 0; u < n; ++u) {
                if (ecc_lb_[u] < ecc_lb_[c] && scc_nb[u] == scc)
                    c = u;
            }
            Rcoballs.push_back(c);
            if (c != p || directed) {// c == p at iter 2
                WL e_c = sweep(c);
                // not good for diam cert size :
                // if (optim_diam && c != start) improves_ecc_ub(diam_lb); 
                if (optim_all_ecc) {
                    if (directed) sweep(c, true);
                    improves_ecc_ub(e_c);
                }
            }

            // Sweep from pseudo-peripheral node : the furthest from c
            p = trav.last_visited();
            WL d_p = trav.dist(p);
            //Rcoballs.push_back(p);
            if ( ! is_in(p, P)) { // possibly false at iter 2 and last iter
                sweep(p, true);
                assert(d_p == trav.dist(c));
                improves_ecc_lb();
            }
            
            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.radius()", 1)
                    << nsweep <<" sweeps, todo: "<< rad_todo <<" / "<< n <<"\n";
        }

        verb::cerr("eccentricity.radius()", 1)
            << nsweep <<" sweeps,"
            << " R_ub=" << rad_ub <<" D_lb="<< diam_lb
            <<" P:"<< P.size() <<" C:" << C.size() << std::endl;


        rad_certif = (optim_diam || optim_all_ecc) ? P
                                                   : optim_lb_certif(P, rad_ub);
        // to check correct : rad_certif = optim_rad_certif(rad_certif);

        verb::cerr("eccentricity.radius()", 1)
            << "rad_cert: " << rad_certif.size() <<"\n";

        return rad_ub;
    }

    WL radius_certif_approx(WL rad_estim, int sample_size=1,
                            V start = G::not_vertex, int scc = -1) {
        if (scc == -1) {
            scc = (start != G::not_vertex) ? scc_nb[start] : scc_largest;
        }

        // recover scc and select random sample
        sample.clear();
        for (int u = 0; u < n; ++u) {
            if (scc_nb[u] == scc) sample.push_back(u);
        }
        for (int i = 0; i < sample.size()-1 && i+1 < n; ++i) {
            int r = i + (127L + 1237L * i) % (sample.size() - i);
            std::swap(sample[i], sample[r]);
        }
        
        rad_ub = rad_estim;
        int total_sampled = 0;

        while (rad_todo > 0) {
            count_far_from_sample(rad_ub, sample_size);
            V p = argmax(sample_count); // greedy set cover

            assert(sample_count[p] > 0); // otherwise we should resample
            assert( ! is_in(p, P));

            WL e_p = sweep(p, true);
            assert(e_p >= rad_estim); // rad_estim seems too low!
            improves_ecc_lb();
            
            update_todo();
            if (true || verb::progress())
                verb::cerr("eccentricity.radius()", 1)
                    << sample_count[p] <<" covered, "
                    << nsweep <<" sweeps, todo: "<< rad_todo <<" / "<< n <<"\n";
        }

        rad_certif = optim_lb_certif(P, rad_ub);

        verb::cerr("eccentricity.radius_certif_approx()", 1)
            << "total_sampled: " << total_sampled
            << " P: " << P.size()
            << " rad_cert: " << rad_certif.size() <<"\n";
        
        return rad_ub;
    }

    // sweeps from sample, returns the total number of nodes used for sweeps
    int count_far_from_sample(WL rad_estim, int sample_size) {
        std::fill(sample_count.begin(), sample_count.end(), 0);
        int i = 0;
        for (int nb = 0; i < sample.size() && nb < sample_size / 2; ++i) {
            V u = sample[i];
            if (ecc_lb_[u] < rad_ub) {
                ++nb;
                sweep(u);
                for (int v = 0; v < n; ++v) {
                    if (trav.dist(v) >= rad_ub) ++(sample_count[v]);
                }
            }
        }
        return i;
    }
    
    /*
    WL diameter_set_cov(V start = G::not_vertex) {
        if (rad_todo > 0) radius(start, true);
        int sweep_budget = rad_nsweep;
        while (diam_todo > 0) {
            // heur. improving all lower bounds, & upper bounds of central nodes
            all_basic(sweep_budget, sweep_budget, false);
            // if lower bounds are not good enough, avoid too much pruned sweeps
            optim_ub_when_lb_tight(false, sweep_budget);
            sweep_budget *= 2;
        }
        diam_certif = optim_ub_certif_one_shot(false, true);
        //diam_certif = optim_diam_certif(C);
        return diam_lb;
    }
    
    WL diameter_1(V start = G::not_vertex) {
        if (rad_todo > 0) radius(start, true);
        if (diam_todo > 0) {
            // heur. improving all lower bounds, & upper bounds of central nodes
            all_basic(rad_nsweep, 5*rad_nsweep, false); 
            optim_ub_when_lb_tight(false);
        }
        diam_certif = optim_ub_certif_one_shot(false, true);
        //diam_certif = optim_diam_certif(C);
        return diam_lb;
    }
    
    WL diameter_basic(V start = G::not_vertex, start_diam do_rad = SUMSWEEP,
            bool choose_ub = true, bool decr_ecc = true, int loose_ineq = 0) {

        if (do_rad == RADIUS || rad_todo == 0) {
            radius(start, true);
            sweep(rad_node);
        } else if (do_rad == SUMSWEEP) {
            sweep(sum_sweep(start));
        } else {
            sweep(start == G::not_vertex ? max_degree_node() : start);
        }
        improves_ecc_ub(diam_lb);
        update_todo();

        std::vector<int> trav_nsons(n);

        int iter = 0;
        while (diam_todo > 0) {

            // Sweep from pseudo-peripheral node : node with maximal ecc_ub_
            V p = G::not_vertex;
            for (V u = 0; u < n; ++u) {
                if (ecc_ub_[u] > diam_lb
                    && (p == G::not_vertex ||
                        (choose_ub ? ecc_ub_[u] > ecc_ub_[p]
                                   : ecc_lb_[u] > ecc_lb_[p])))
                    p = u;
            }
            assert(p != G::not_vertex);
            
            WL e_p = sweep(p);
            //faster but bigger C: improves_ecc_ub(diam_lb);
            
            // Find furthest certif for p
            V c = p; //furthest_certif_for_root(diam_lb);
                int f_e = 2 - loose_ineq, f_D = loose_ineq;
                WL e_D = (f_e * e_p + f_D * diam_lb) / 2;
            if ( ! decr_ecc) {
                for (int i = 0; i < trav.nvis() -1; ++i) {
                    V u = trav.visit(i);
                    if (trav.dist(u) + ecc_lb_[u] <= e_D
                        )// &&  ecc_lb_[u] < ecc_lb_[c])
                        c = u;
                }
            } else {
                for (int i = 0; i < trav.nvis() -1; ++i) {
                    V u = trav.visit(i);
                    if (trav.dist(u) + ecc_lb_[u] <= e_D
                        &&  ecc_lb_[u] < ecc_lb_[c])
                        c = u;
                }
            }

            WL e_c = sweep(c);
            improves_ecc_ub(diam_lb);
            
            // Check that lower bound was tight (used for choosing certif) :
            if (ecc_lb_[c] < e_c) {
                V f = trav.last_visited();
                sweep(f);
                improves_ecc_lb();
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter_basic()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }
    */

    // returns a pseudo-center 
    V sum_sweep(V start = G::not_vertex, int iter = 3) {
        if (directed) return sum_sweep_dir(start);

        for (V u = 0; u < n; ++u)
            assert(ecc_lb_[u] == zero_weight && sum_[u] == zero_weight);

        if (start == G::not_vertex) start = max_degree_node();

        verb::cerr(2) << "sumsweep 1 : " << start <<"\n";
        sweep(start);
        improves_ecc_lb();
        update_sums();

        while (--iter > 0) {
            V u = argmax_sum();
            verb::cerr(2) << "sumsweep   : " << u <<"\n";
            sweep(u);
            improves_ecc_lb();
            update_sums();
        }
        
        V c = argmin_lb();
        assert(c != G::not_vertex);
        return c;
    }

    V sum_sweep_dir(V start = G::not_vertex, int iter = 3) {
        assert(directed);
        for (V u = 0; u < n; ++u)
            assert(ecc_lb_[u] == zero_weight && sum_[u] == zero_weight);

        if (start == G::not_vertex) start = max_degree_node();
        
        sweep(start);
        save_trav(); // store rev distances

        while (--iter > 0) {
            V v = argmax_save(); // max rev sum
            sweep(v, true); // backward
            improves_ecc_lb();
            update_sums();

            if (iter > 0) { // used only for [v] selection
                V u = argmax_sum();
                sweep(u);
                save_trav_add();
            }
        }
        
        V c = argmin_lb();
        assert(c != G::not_vertex);
        return c;
    }
        

    WL diameter_sample(V start = G::not_vertex, start_diam do_rad = SUMSWEEP,
                       bool find_rad = false) {

        diameter_start(start, do_rad);

        std::vector<std::vector<V>> certif_for(n);
        std::vector<bool> in_sample(n, false);
        //for (V x = 0; x < n; ++x) certif_for[x].clear();
        
        int sum_certif = 0, n_sample = 0, wait_rounds = 3;
        int sample_round = 0, maxub_round = 0, clear_round = 0, total_sampled=0;
        int maxub_ncov = 1, x_ncov = 0, avg_deg = (graph.m() + n - 1) / n;
        int smooth_ncov = n;

        int iter = 0;
        while (diam_todo > 0) {
            ++iter;
            
            WL diam_prev_lb = diam_lb;

            auto clear_cert = [&certif_for, &in_sample,
                               this, iter, &clear_round, &n_sample,
                               &sum_certif, &diam_prev_lb](std::string msg=""){
                clear_round = iter;
                for (V y = 0; y < n; ++y) {
                    certif_for[y].clear();
                }
                sum_certif = 0;
                for (V v = 0; v < n; ++v) {
                    in_sample[v] = false;
                }
                n_sample = 0;
                diam_prev_lb = diam_lb;
                //verb::cerr() << iter << " ---- clear --- " << msg <<"\n";
                update_todo();
            };
            
            if (find_rad && rad_todo > 0 && iter % wait_rounds == 1) {
                V c = G::not_vertex;
                for (V v = 0; v < n; ++v) {
                    if (c == G::not_vertex || ecc_lb_[v] < ecc_lb_[c]) c = v;
                }
                assert(c != G::not_vertex);

                WL e_c = sweep(c);

                verb::cerr() << iter << " radius sweep : "<< c 
                             <<"  e_c: "<< ecc_lb_[c] <<" <= "<< e_c <<"\n";

                if(ecc_lb_[c] < e_c) {
                    V p = trav.last_visited();
                    sweep(p, true);
                    assert(improves_ecc_lb());
                    clear_cert();
                    if (rad_todo > 0) continue;
                    for (V v = 0; v < n; ++v) {
                        if (ecc_lb_[v] < ecc_lb_[c]) c = v;
                    }
                } else clear_cert("rad");
                // else c is a center
                e_c = sweep(c);
                if (directed) pruned_sweep_from_certifier(c, trav, diam_lb,e_c);
                improves_ecc_ub(e_c, diam_lb);
                update_todo();
                /*
                if (diam_lb > diam_prev_lb) clear_cert("diam");
                for (V v = 0; v < n; ++v) {
                    in_sample[v] = false;
                }
                */
                continue;
            }

            // sample random uncover nodes
            int in_this_round = 0, n_tot_round = 0, mprune_rnd = mprune;

            V vi = 0;
            while (iter >= wait_rounds
                   // limit mem usage compared to graph size:
                   && sum_certif < (32 + 3 * avg_deg) * n
                   // increase slowly:
                   && n_sample <= iter
                   // bound variance under smooth_ncov (Chernoff):
                   && n_sample < (diam_todo < total_sampled / 3
                         ? diam_todo
                         : std::max(16,
                           16 * diam_todo / (smooth_ncov<=0 ? 1 : smooth_ncov)))
                   && n_sample < diam_todo
                   // at most one sweep equiv
                   && (mprune - mprune_rnd) < n * avg_deg
                   ) {

                V r = G::not_vertex; // pseudo-random uncovered node
                while (true) {
                    V v = (127L * iter + 1237L * vi) % n; // pseudo rnd order
                    ++vi;
                    if (verb::progress(0.2))
                        verb::cerr() << iter
                                     <<" v, vi "<< v <<", "<< vi
                                     <<" todo:"<< diam_todo <<"\n";
                    if (ecc_ub_[v] > diam_lb && ! in_sample[v]) {
                        r = v;
                        break;
                    }
                    if (vi > n) {
                        int ntodo = 0, nsmpl = 0, ntodsmp = 0;
                        for (V w = 0; w < n; ++w) {
                            if (ecc_ub_[w] > diam_lb) ++ntodo;
                            if (in_sample[w]) ++nsmpl;
                            if (ecc_ub_[w] > diam_lb && in_sample[w]) ++ntodsmp;
                        }
                        verb::cerr() << "ntodo:" << ntodo
                                     << " nsmpl:" << nsmpl
                                     << " ntodsmp:" << ntodsmp
                                     << " todo:" << diam_todo
                                     << " sample:"<< n_sample
                                     <<"\n";
                    }
                    assert(vi <= n);
                }
                assert(r != G::not_vertex);

                pruned_sweep_to_certifiers(r, trav, diam_lb, true);

                ++n_sample;
                ++total_sampled;
                ++n_tot_round;
                in_sample[r] = true;
                for (int i = trav.nvis() - 1; i >= 0; --i) {
                    V x = trav.visit(i);
                    certif_for[x].push_back(r);
                    ++sum_certif;
                    ++in_this_round;
                }

                if (verb::progress(0.2))
                    verb::cerr() << iter << " sample   sum : " << sum_certif
                                 << " nsample: "<< n_sample
                                 << " total_sampled: " << total_sampled
                                 <<" r, vi "<< r <<", "<< vi
                                 <<" todo:"<< diam_todo <<"\n";
            }
            
            V x = G::not_vertex;
            for (V v = 0; v < n; ++v) {
                if (x == G::not_vertex
                    || certif_for[v].size() > certif_for[x].size()) {
                    x = v;
                }
            }
            assert(x != G::not_vertex);
            int cert_max = certif_for[x].size();
            for (V v = 0; v < n; ++v) {
                if (4 * certif_for[v].size() >= 3 * cert_max
                    && ecc_lb_[v] < ecc_lb_[x]) {
                    x = v;
                }
            }

            // try node covering max of uncovered according to sampling
            int n_tocov = diam_todo;
            int expected_cov = certif_for[x].size() * (n_tocov+1)/(n_sample+1);
            x_ncov = 0;
            
            if (expected_cov * 2 >= maxub_ncov
                && iter > wait_rounds) {

                WL e_x = sweep(x);
                V p = trav.last_visited();
                if (directed) pruned_sweep_from_certifier(x, trav, diam_lb,e_x);

                x_ncov = 0;
                for (int i = 0; i < trav.nvis(); ++i) {
                    V v = trav.visit(i);
                    if (ecc_ub_[v] > diam_lb && scc_nb[v] == scc_nb[x]
                        && trav.dist(v) + e_x <= diam_lb)
                        ++x_ncov;
                }

                /*verb::cerr() << iter <<" "<< nsweep
                             << " x : " << x
                             <<" "<< certif_for[x].size()
                             <<" n_sample:"<< n_sample<<" tot:"<< total_sampled
                             <<" sum_cert:"<< sum_certif <<" tocov "<< n_tocov
                             <<" expt:"<< expected_cov
                             <<" ncov:"<< x_ncov <<"\n";*/

                if (x_ncov * 4 >= expected_cov && x_ncov > 0) {
                    // try in certificate anyway
                    assert(improves_ecc_ub(e_x, diam_lb));
                    update_todo();
                } else {
                    /*verb::cerr() << " ------------ bad cov "
                      << "expected "<< expected_cov <<" "<< x_ncov<<"\n";*/
                }

                if (x_ncov > 2 * expected_cov) {
                    /* verb::cerr() << " ------ too much cov "
                       << "expected "<< expected_cov <<" "<< x_ncov<<"\n";*/
                    smooth_ncov = (3 * smooth_ncov + x_ncov + 1) / 4;
                }
                
                if (e_x > ecc_lb_[x]) {
                    WL e_p = sweep(p, true);
                    assert(improves_ecc_lb());
                    clear_cert();
                    continue;
                } // else :

                if (x_ncov * 2 < expected_cov) {
                    /* if (x_ncov * 4 >= expected_cov)
                         verb::cerr() << " ------ low cov "
                           << "expected "<< expected_cov <<" "<< x_ncov<<"\n";*/
                    smooth_ncov = (2 * smooth_ncov + x_ncov + 1) / 3;
                }

                // in all cases:
                sum_certif -= certif_for[x].size();
                for (V v : certif_for[x]) {
                    if (ecc_ub_[v] <= diam_lb && in_sample[v]) {
                        in_sample[v] = false;
                        --n_sample;
                    }
                }
                certif_for[x].clear();
                if (diam_lb > diam_prev_lb) {
                    clear_cert("diam");
                    continue;
                }

            } else { // max ub node (basic heuristic)

                maxub_round = iter;
                n_tocov = diam_todo;
                V w = G::not_vertex;
                for (V v = 0; v < n; ++v) {
                    if (w == G::not_vertex || ecc_ub_[v] > ecc_ub_[w])
                        w = v;
                }
                assert(w != G::not_vertex && ecc_ub_[w] > diam_lb);

                V x = pruned_sweep_to_certifiers(w, trav, ecc_lb_[w]);

                WL e_x = sweep(x);
                V p_x = trav.last_visited();
                if (directed)
                    pruned_sweep_from_certifier(x, trav, diam_lb, e_x);

                assert(improves_ecc_ub(e_x, diam_lb) || e_x > ecc_lb_[x]);
                update_todo();
                maxub_ncov = (2*(n_tocov - diam_todo) + maxub_ncov + 1) / 3;

                // in all cases:
                sum_certif -= certif_for[x].size();
                for (V v : certif_for[x]) {
                    if (ecc_ub_[v] <= diam_lb && in_sample[v]) {
                        in_sample[v] = false;
                        --n_sample;
                    }
                }
                certif_for[x].clear();
                if (diam_lb > diam_prev_lb) {
                    clear_cert("diam");
                    continue;
                }

                /* verb::cerr() << iter <<" "<< nsweep << " max ub : " << x
                             <<" e_x :"<< ecc_lb_[x] <<", "<< e_x
                             <<" tsmpled:"<< total_sampled
                             <<" todo:"<< diam_todo
                             <<" maxubcov:"<< maxub_ncov
                             <<" ncov:"<< n_tocov - diam_todo <<"\n"; */

                    
                if (e_x > ecc_lb_[x]) {
                    WL e_p = sweep(p_x, true);
                    assert(improves_ecc_lb());
                    clear_cert();
                    continue;
                } // else:
            }
                
            // remove covered nodes in certif_for[y] for all y
            if (diam_todo > 0) {
                for (V v = 0; v < n; ++v) {
                    if (ecc_ub_[v] <= diam_lb && in_sample[v]) {
                        in_sample[v] = false;
                        --n_sample;
                    }
                }
                for (V y = 0; y < n; ++y) {
                    int n_uncov = 0, size = certif_for[y].size();
                    for (int i = 0; i < size; ++i) {
                        V v = certif_for[y][i];
                        if (ecc_ub_[v] > diam_lb) {
                            certif_for[y][n_uncov] = v;
                            ++n_uncov;
                        }
                    }
                    sum_certif += n_uncov - size;
                    certif_for[y].resize(n_uncov);
                    // count this scan in terms of sweeps:
                    mprune += size;
                    if (mprune > graph.m()) {
                        mprune -= graph.m()+1; ++nsweep; ++nprune_sweep;
                    }
                }
            }
                
            if (verb::progress())
                verb::cerr("eccentricity.diameter_sample()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n
                    << " n_smpl=" << n_sample
                    <<" sum_cerrt: "<< sum_certif <<"< "<< (32 + 3 * avg_deg)*n
                    <<"\n";
        }

        if (mprune > 0) {
            mprune = 0; ++nsweep; ++nprune_sweep;
        }
        
        verb::cerr("eccentricity.diameter_sample()", 1)
                <<" total_sampled " << total_sampled
                <<",  "<< nsweep <<" sweeps\n";

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }

    /*
    WL diameter_diff(V start = G::not_vertex, start_diam do_rad = SUMSWEEP,
                     float f_ub = 1.1) {

        if (do_rad == RADIUS || rad_todo == 0) {
            radius(start, true);
            sweep(rad_node);
        } else if (do_rad == SUMSWEEP) {
            sweep(sum_sweep(start));
        } else {
            sweep(start == G::not_vertex ? max_degree_node() : start);
        }
        improves_ecc_ub(diam_lb);
        update_todo();

        int iter = 0;
        while (diam_todo > 0 && ++iter < 5000) {

            // Sweep from pseudo-peripheral node : node with maximal ecc_ub_
            V p = G::not_vertex;
            for (V u = 0; u < n; ++u) {
                if (ecc_ub_[u] > diam_lb
                    && (p == G::not_vertex ||
                        (f_ub * ecc_ub_[u] - ecc_lb_[u] >
                         f_ub * ecc_ub_[p] - ecc_lb_[p])))
                    p = u;
            }
            assert(p != G::not_vertex);
            
            WL e_p = sweep(p);
            improves_ecc_ub(diam_lb);
            
            // Check that lower bound was tight (used for choosing certif) :
            if (ecc_lb_[p] < e_p) {
                V f = trav.last_visited();
                sweep(f);
                improves_ecc_lb();
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter_diff()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }

        if (diam_todo > 0) {
            verb::cerr("eccentricity.diameter_diff()", 0) << "failed !!!!!!\n";
        } else {
            diam_certif = optim_ub_certif(C, diam_lb);
        }

        return diam_lb;
    }
    */

    void diameter_start(V start = G::not_vertex, start_diam do_rad = SUMSWEEP) {
        WL e_start = max_weight;
        if (do_rad == RADIUS || rad_todo == 0) {
            radius(start, true);
            e_start = rad_ub;
            sweep(rad_node, true);
        } else if (do_rad == SUMSWEEP) {
            start = sum_sweep(start);
            sweep(start); // asserts undirected
            e_start = trav.dist(trav.last_visited());
        } else {
            if (start == G::not_vertex) start = max_degree_node();
            sweep(start);
            e_start = trav.dist(trav.last_visited());
            if (directed) sweep(start, true);
        }
        assert(improves_ecc_ub(e_start, diam_lb));
        update_todo();
    }

    void diameter_update(std::string fname = "") {
        update_todo();
        if (verb::progress())
            verb::cerr("eccentricity." + fname + "()", 1)
                << nsweep <<" sweeps, |C|="<< C.size()
                <<" todo: "<< diam_todo <<" / "<< n <<"\n";
    }

    WL diameter_bare(V start = G::not_vertex, start_diam do_rad = SUMSWEEP,
                     int loose_ineq = 0, bool optim_certif = true) {
        diameter_start(start, do_rad);

        int iter = 0;
        while (diam_todo > 0 && ++iter < 11000) {
            V u = argmax_ub();
            assert(u != G::not_vertex && ecc_ub_[u] > diam_lb);
            
            WL e_u = sweep(u);
            save_trav();
            if (diam_todo == 0) break; // can happen if diam_lb improved
            
            V x = save_argmin_ecc([this, u, e_u, loose_ineq](V v, WL e_v) { 
                    int f_e = 2 - loose_ineq, f_D = loose_ineq;
                    WL e_u_D = (f_e * e_u + f_D * diam_lb) / 2;
                    if (scc_nb[v] == scc_nb[u] && save_dist[v] + e_v <= e_u_D)
                        return e_v;
                    else
                        return max_weight;
                });
            assert(last_sweep_node == x && ! last_sweep_backward);
            WL e_x = ecc_lb_[x]; // this is tight (= ecc(x))
            if (directed) sweep(x, true);
            assert(improves_ecc_ub(e_x, diam_lb) || e_u == diam_lb); 

            diameter_update("diameter_bare");
        }

        if (diam_todo > 0) {
            verb::cerr("eccentricity.diameter_bare()", 0) << "failed !!!!!!\n";
        } else {
            if (optim_certif) diam_certif = optim_ub_certif(C, diam_lb);
            else diam_certif = C;
        }
        
        return diam_lb;
    }

    WL diameter_bare1(V start = G::not_vertex, start_diam do_rad = SUMSWEEP,
                     bool choose_p = true, bool min_ecc = true,
                     int loose_ineq = 0) {
        diameter_start(start, do_rad);

        int iter = 0;
        while (diam_todo > 0 && ++iter < 11000) {

            // Sweep from pseudo-peripheral node : node with maximal ecc_ub_
            V p = G::not_vertex;
            for (V u = 0; u < n; ++u) {
                if (ecc_ub_[u] > diam_lb
                    && (p == G::not_vertex || ecc_ub_[u] > ecc_ub_[p]))
                    p = u;
            }
            assert(p != G::not_vertex);
            
            WL e_p = sweep(p);
            
            // Certif for p
            V c = p;
            if (choose_p) {

                if (directed) sweep(p, true);
                improves_ecc_ub(e_p, diam_lb);

            } else { // exact certif with smallest ecc
                int f_e = 2 - loose_ineq, f_D = loose_ineq;
                WL e_p_D = (f_e * e_p + f_D * diam_lb) / 2;

                WL e_c = e_p;
                for (int i = 0; i < trav.nvis(); ++i) {
                    V u = trav.visit(i);
                    if (scc_nb[u] == scc_nb[p]
                        && trav.dist(u) + ecc_lb_[u] <= e_p_D
                        &&  (min_ecc ? ecc_lb_[u] < ecc_lb_[c]
                             : ecc_lb_[u] <= ecc_lb_[c])) {
                        c = u;
                    }
                }

                e_c = sweep(c);

                // Check lower bound was tight (used for choosing certif) :
                if (ecc_lb_[c] < e_c) {
                    V f = trav.last_visited();
                    sweep(f, true);
                    improves_ecc_lb();
                    // try again
                } else {
                    if (directed) sweep(c, true);
                    improves_ecc_ub(e_c, diam_lb);
                }
            
            }

            diameter_update("diameter_bare1");
        }
        
        /*
        sweep(rad_node);
        int in_ball_ctr = 0;
        for (int c : C) {
            if (trav.dist(c) + rad_ub <= diam_lb) ++in_ball_ctr;
        }
        verb::cerr() << "in_ball_ctr: "<< in_ball_ctr
                     <<" choose_p: " << choose_p <<"\n";
        */

        if (diam_todo > 0) {
            verb::cerr("eccentricity.diameter_bare1()", 0) << "failed !!!!!!\n";
        } else {
            diam_certif = optim_ub_certif(C, diam_lb);
        }
        
        return diam_lb;
    }


    WL diameter_bi(V start = G::not_vertex, start_diam do_rad = SUMSWEEP,
                     int loose_ineq = 0) {
        assert(strongly_connected && ! directed);
        
        diameter_start(start, do_rad);

        int iter = 0;
        while (diam_todo > 0) {
            V u = argmax_ub();
            assert(u != G::not_vertex && ecc_ub_[u] > diam_lb);
            
            V x = pruned_sweep_to_certifiers(u, trav, ecc_lb_[u]);
            WL e_x = sweep(x);
            improves_ecc_ub(e_x, diam_lb);
            save_trav(); // dist to x
            
            V y = trav.last_visited();
            //if (ecc_lb_[x] < e_x) {            
            WL e_y = sweep(y);
            improves_ecc_lb();
            //    continue;
            //} // else ecc_lb_[x] == e_x :
            
            for (int i = trav.nvis() - 1; i >= 0; --i) {
                V v = trav.visit(i);
                if (trav.dist(v) <= diam_lb / 2) { y = v; break; }
            }
            
            improves_ecc_ub(e_y, diam_lb);

            // improve ub:
            for (V v = 0; v < n; ++v) {
                WL vx = save_dist[v], vy = trav.dist(v);
                WL ub_max = zero_weight;
                for (V w = 0; w < n; ++w) {
                    WL xw = save_dist[w], yw = trav.dist(w);
                    WL ub = std::min(vx + xw, vy + yw);
                    if (ub > ub_max) { ub_max = ub;  }
                }
                if (ub_max < ecc_ub_[v]) ecc_ub_[v] = ub_max;
            }
            
            diameter_update("diameter_bi");
        }

        return diam_lb;
    }



    WL diameter_tk11(V start = G::not_vertex, start_diam do_rad = SUMSWEEP) {
        diameter_start(start, do_rad);

        auto ecc_untight = [this](V w) { return ecc_lb_[w] < ecc_ub_[w]; };

        int iter = 0;
        while (diam_todo > 0 && ++iter < 11000) {
            
            V v = argmax_ub(ecc_untight);
            WL e_v = sweep(v);
            improves_ecc_lb(max_weight, true);
            improves_ecc_ub(e_v);
            
            update_todo();
            if (diam_todo == 0) break;

            V u = argmin_lb(ecc_untight);
            WL e_u = sweep(u);
            improves_ecc_lb(max_weight, true);
            improves_ecc_ub(e_v);
            
            diameter_update("diameter_tk11");
        }
        
        return diam_lb;
    }


    WL diameter_sumsweep(V start = G::not_vertex,
                         start_diam do_rad = SUMSWEEP) {
        diameter_start(start, do_rad);

        auto ecc_untight = [this](V w) { return ecc_lb_[w] < ecc_ub_[w]; };

        int iter = 0;
        while (diam_todo > 0 && ++iter < 11000) {
            
            V v = argmax_ub(ecc_untight);
            WL e_v = sweep(v);
            improves_ecc_lb(max_weight, true);
            improves_ecc_ub(e_v);
            
            update_todo();
            if (diam_todo == 0) break;

            V u = argmin_lb(ecc_untight);
            WL e_u = sweep(u);
            improves_ecc_lb(max_weight, true);
            improves_ecc_ub(e_v);
            
            diameter_update("diameter_tk11");
        }
        
        return diam_lb;
    }


    
    /*
    WL diameter_basic_lb(V start = G::not_vertex, start_diam do_rad = SUMSWEEP,
                         bool choose_ub = true) {

        if (do_rad == RADIUS || rad_todo == 0) {
            radius(start, true);
            sweep(rad_node);
        } else if (do_rad == SUMSWEEP) {
            sweep(sum_sweep(start));
        } else {
            sweep(start == G::not_vertex ? max_degree_node() : start);
        }
        improves_ecc_ub(diam_lb);
        update_todo();

        std::vector<int> trav_nsons(n);

        int iter = 0;
        while (diam_todo > 0) {

            // Sweep from pseudo-peripheral node : node with maximal ecc_ub_
            V p = G::not_vertex;
            for (V u = 0; u < n; ++u) {
                if (ecc_ub_[u] > diam_lb
                    && (p == G::not_vertex ||
                        (choose_ub ? ecc_ub_[u] > ecc_ub_[p]
                                   : ecc_lb_[u] > ecc_lb_[p])))
                    p = u;
            }
            assert(p != G::not_vertex);
            
            WL e_p = sweep(p);
            //faster but bigger C: improves_ecc_ub(diam_lb);
            if (ecc_lb_[p] < e_p) {
                V f = trav.last_visited();
                sweep(f);
                improves_ecc_lb();
            } else {
            
                // Find furthest certif for p
                V c = p; //furthest_certif_for_root(diam_lb);
                for (int i = 0; i < trav.nvis() -1; ++i) {
                    V u = trav.visit(i);
                    if (trav.dist(u) + ecc_lb_[u] <= diam_lb
                        )// &&  ecc_lb_[u] < ecc_lb_[c])
                        c = u;
                }
                
                WL e_c = sweep(c);
                improves_ecc_ub(diam_lb);
                
                // Check that lower bound was tight (used for choosing certif) :
                if (ecc_lb_[c] < e_c) {
                    V f = trav.last_visited();
                    sweep(f);
                    improves_ecc_lb();
                }

            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter_basic_lb()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }

    WL diameter_stopcov(V start = G::not_vertex, start_diam do_rad = SUMSWEEP,
                        bool choose_ub = true, int div = 2) {

        if (do_rad == RADIUS || rad_todo == 0) {
            radius(start, true);
            sweep(rad_node);
        } else if (do_rad == SUMSWEEP) {
            sweep(sum_sweep(start));
        } else {
            sweep(start == G::not_vertex ? max_degree_node() : start);
        }
        improves_ecc_ub(diam_lb);
        update_todo();

        std::vector<int> trav_nsons(n);

        int iter = 0;
        while (diam_todo > 0) {

            // Sweep from pseudo-peripheral node : node with maximal ecc_ub_
            V p = G::not_vertex;
            for (V u = 0; u < n; ++u) {
                if (ecc_ub_[u] > diam_lb
                    && (p == G::not_vertex ||
                        (choose_ub ? ecc_ub_[u] > ecc_ub_[p]
                                   : ecc_lb_[u] > ecc_lb_[p])))
                    p = u;
            }
            assert(p != G::not_vertex);
            
            WL e_p = sweep(p);
            //faster but bigger C: improves_ecc_ub(diam_lb);

            WL dist_cov = max_weight;
            for (int i = 0; i < trav.nvis() - 1; ++i) {
                V u = trav.visit(i);
                if (ecc_ub_[u] <= diam_lb) break;
                else dist_cov = trav.dist(u);
            }
            
            // Find furthest certif for p
            V c = p; //furthest_certif_for_root(diam_lb);
            for (int i = 0; i < trav.nvis() - 1; ++i) {
                V u = trav.visit(i);
                if (dist_cov / div + trav.dist(u) + ecc_lb_[u] <= diam_lb)
                    c = u;
            }

            WL e_c = sweep(c);
            improves_ecc_ub(diam_lb);
            
            // Check that lower bound was tight (used for choosing certif) :
            if (ecc_lb_[c] < e_c) {
                V f = trav.last_visited();
                sweep(f);
                improves_ecc_lb();
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter_stopcov()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }

    WL diameter_basic_minecc(V start = G::not_vertex,
                             start_diam do_rad = SUMSWEEP,
                             bool choose_ub = true) {

        if (do_rad == RADIUS || rad_todo == 0) {
            radius(start, true);
            sweep(rad_node);
        } else if (do_rad == SUMSWEEP) {
            sweep(sum_sweep(start));
        } else {
            sweep(start == G::not_vertex ? max_degree_node() : start);
        }
        improves_ecc_ub(diam_lb);
        update_todo();
        
        std::vector<int> trav_nsons(n);

        int iter = 0;
        while (diam_todo > 0) {

            // Sweep from pseudo-peripheral node : node with maximal ecc_ub_
            V p = G::not_vertex;
            for (V u = 0; u < n; ++u) {
                if (ecc_ub_[u] > diam_lb
                    && (p == G::not_vertex ||
                        (choose_ub ? ecc_ub_[u] > ecc_ub_[p]
                                   : ecc_lb_[u] > ecc_lb_[p])))
                    p = u;
            }
            assert(p != G::not_vertex);
            
            WL e_p = sweep(p);
            //faster but bigger C: improves_ecc_ub(diam_lb);
            
            // Find furthest certif for p
            V c = p; //furthest_certif_for_root(diam_lb);
            for (int i = 0; i < trav.nvis() -1; ++i) {
                V u = trav.visit(i);
                if (trav.dist(u) + ecc_lb_[u] <= diam_lb
                     &&  ecc_lb_[u] <= ecc_lb_[c])
                    c = u;
            }

            WL e_c = sweep(c);
            improves_ecc_ub(diam_lb);
            
            // Check that lower bound was tight (used for choosing certif) :
            if (ecc_lb_[c] < e_c) {
                V f = trav.last_visited();
                sweep(f);
                improves_ecc_lb();
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter_basic_minecc()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }

    WL diameter_pack_approx(V start = G::not_vertex, bool do_rad = true, bool choose_ub = true, double slack = 2.) {
        if (do_rad) {
            radius(start, true);
        } else {
            if (start == G::not_vertex) start = max_degree_node();
            sweep(start);
            improves_ecc_ub();
        }

        std::vector<int> trav_nsons(n);

        int iter = 0;
        while (diam_todo > 0) {

            // Sweep from pseudo-peripheral node : node with maximal ecc_ub_
            V p = G::not_vertex;
            for (V u = 0; u < n; ++u) {
                if (ecc_ub_[u] > diam_lb
                    && (p == G::not_vertex ||
                        (choose_ub ? ecc_ub_[u] > ecc_ub_[p]
                                   : ecc_lb_[u] > ecc_lb_[p])))
                    p = u;
            }
            assert(p != G::not_vertex);
            
            WL e_p = sweep(p);
            if (choose_ub && ecc_lb_[p] < e_p) {
                V f = trav.last_visited();
                sweep(f);
                improves_ecc_lb();
                update_todo();
                continue;
            }
            
            // Find furthest certif c for p that covers also nodes in B(u,uc)
            V c = p; //furthest_certif_for_root(diam_lb);
            for (V u = 0; u < n; ++u) trav_nsons[u] = 0;
            for (int i = 0; i < trav.nvis() -1; ++i) {
                V u = trav.visit(i);
                trav_nsons[trav.parent(u)] += 1;
            }
            V i_one = 0;
            while (i_one < n && trav_nsons[trav.visit(i_one)] == 1) ++i_one;
            WL dist_one = trav.dist(std::min(n-1, i_one));
            for (int i = i_one; i < trav.nvis() -1; ++i) {
                V u = trav.visit(i);
                if (trav.dist(u) + ecc_lb_[u] <= diam_lb
                    && dist_one + slack * (trav.dist(u) - dist_one)
                       + ecc_lb_[u] <= diam_lb
                    )// &&  ecc_lb_[u] < ecc_lb_[c])
                    c = u;
            }
            
            WL e_c = sweep(c);
            improves_ecc_ub(diam_lb);
            
            // Check that lower bound was tight (used for choosing certif) :
            if (ecc_lb_[c] < e_c) {
                V f = trav.last_visited();
                sweep(f);
                improves_ecc_lb();
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter_pack_approx()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }

    WL diameter_dist_C(V start = G::not_vertex, bool do_rad = true, bool choose_ub = true) {
        if (do_rad) {
            radius(start, true);
        } else {
            if (start == G::not_vertex) start = max_degree_node();
            sweep(start);
            improves_ecc_ub();
        }

        std::vector<int> trav_nsons(n);

        int iter = 0;
        while (diam_todo > 0) {

            // Sweep from pseudo-peripheral node : node with maximal ecc_ub_
            V p = G::not_vertex;
            for (V u = 0; u < n; ++u) {
                if (ecc_ub_[u] > diam_lb
                    && (p == G::not_vertex ||
                        (choose_ub ? ecc_ub_[u] > ecc_ub_[p]
                                   : ecc_lb_[u] > ecc_lb_[p])))
                    p = u;
            }
            assert(p != G::not_vertex);
            
            WL e_p = sweep(p);
            //faster but bigger C: improves_ecc_ub(diam_lb);
            
            // Find certif for p which is closest to C
            V c = p; 
            for (int i = 0; i < trav.nvis() -1; ++i) {
                V u = trav.visit(i);
                if (trav.dist(u) + ecc_lb_[u] <= diam_lb
                    && dist_C[u] <= dist_C[c]
                    )// &&  ecc_lb_[u] < ecc_lb_[c])
                    c = u;
            }
            
            WL e_c = sweep(c);
            improves_ecc_ub(diam_lb);
            
            // Check that lower bound was tight (used for choosing certif) :
            if (ecc_lb_[c] < e_c) {
                V f = trav.last_visited();
                sweep(f);
                improves_ecc_lb();
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter_dist_C()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }

        WL diameter_basic_f(V start = G::not_vertex, bool do_rad = true, bool choose_ub = true) {
        if (do_rad) {
            radius(start, true);
        } else {
            if (start == G::not_vertex) start = max_degree_node();
            sweep(start);
            improves_ecc_ub();
        }

        std::vector<int> trav_nsons(n);

        int iter = 0;
        while (diam_todo > 0) {

            // Sweep from pseudo-peripheral node : node with maximal ecc_ub_
            V p = G::not_vertex;
            for (V u = 0; u < n; ++u) {
                if (ecc_ub_[u] > diam_lb
                    && (p == G::not_vertex ||
                        (choose_ub ? ecc_ub_[u] > ecc_ub_[p]
                                   : ecc_lb_[u] > ecc_lb_[p])))
                    p = u;
            }
            assert(p != G::not_vertex);
            
            WL e_p = sweep(p);
            //faster but bigger C: improves_ecc_ub(diam_lb);

            // Find furthest certif for p on path to furthest node
            V c = trav.last_visited();
            while (c != p) {
                if (trav.dist(c) + ecc_lb_[c] <= diam_lb) break;
                c = trav.parent(c);
            }

            WL e_c = sweep(c);
            improves_ecc_ub(diam_lb);
            
            // Check that lower bound was tight (used for choosing certif) :
            if (ecc_lb_[c] < e_c) {
                V f = trav.last_visited();
                sweep(f);
                improves_ecc_lb();
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter_basic_f()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }

        WL diameter_pack_f(V start = G::not_vertex, bool do_rad = true, bool choose_ub = true, double slack = 1.2) {
        if (do_rad) {
            radius(start, true);
        } else {
            if (start == G::not_vertex) start = max_degree_node();
            sweep(start);
            improves_ecc_ub();
        }

        std::vector<int> deep_desc_1(n), deep_desc_2(n);

        int iter = 0;
        while (diam_todo > 0) {

            // Sweep from pseudo-peripheral node : node with maximal ecc_ub_
            V p = G::not_vertex;
            for (V u = 0; u < n; ++u) {
                if (ecc_ub_[u] > diam_lb
                    && (p == G::not_vertex ||
                        (choose_ub ? ecc_ub_[u] > ecc_ub_[p]
                                   : ecc_lb_[u] > ecc_lb_[p])))
                    p = u;
            }
            assert(p != G::not_vertex);
            
            WL e_p = sweep(p);
            if (choose_ub && ecc_lb_[p] < e_p) {
                V f = trav.last_visited();
                sweep(f);
                improves_ecc_lb();
                update_todo();
                continue;
            }

            for (int u = 0; u < n; ++u) {
                deep_desc_1[u] = u;
                deep_desc_2[u] = u;
            }
            for (int i = trav.nvis() - 1; i > 0; --i) {
                V u = trav.visit(i);
                V pu = trav.parent(u);
                if (trav.dist(deep_desc_1[u]) > trav.dist(deep_desc_1[pu]))
                    deep_desc_1[pu] = deep_desc_1[u];
            }
            for (int i = trav.nvis() - 1; i > 0; --i) {
                V u = trav.visit(i);
                V pu = trav.parent(u);
                if (trav.dist(deep_desc_1[u]) > trav.dist(deep_desc_2[pu])
                    && deep_desc_1[u] != deep_desc_1[pu])
                    deep_desc_2[pu] = deep_desc_1[u];
            }

            V c = deep_desc_1[p];
            while (c != p) {
                if (trav.dist(c) + ecc_lb_[c] <= diam_lb) { // c cert for p
                    // require also that B(c,D-ecc(c)) covers any u s.t.
                    // p and u belong to B(c',(D-ecc(c)) / 3) for any c'
                    // approx: cover B(p, slack * uc) (slack = 3 to be sure)
                    // but optimize test using tree structure
                    bool covers = true;
                    V v = trav.parent(c); // node on path to p
                    while (covers && v != p) {
                        WL dist_branch =
                            trav.dist(deep_desc_2[v]) - trav.dist(v);
                        WL c_rad = diam_lb - ecc_lb_[c];
                        if (std::min((double)dist_branch,
                                     slack * trav.dist(c) - trav.dist(v)) 
                            + (trav.dist(c) - trav.dist(v)) > c_rad) {
                            covers = false;
                        }
                        v = trav.parent(v);
                    }
                    if (covers) break;
                }
                c = trav.parent(c);
            }

            WL e_c = sweep(c);
            improves_ecc_ub(diam_lb);
            
            // Check that lower bound was tight (used for choosing certif) :
            if (ecc_lb_[c] < e_c) {
                V f = trav.last_visited();
                sweep(f);
                improves_ecc_lb();
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }
    */

    

    void all(V start = G::not_vertex, bool opt_certif=true, bool prune=true) {

        bool first = true;
        
        while (all_ecc_todo > 0) {

            V c = G::not_vertex;
            if (first) {
                first = false;
                if (start != G::not_vertex) c = start;
                else c = max_degree_node();
            } else {
                // Select a node with minimal untight lower-bound:
                for (V v = 0; v < n; ++v) {
                    if (ecc_lb_[v] < ecc_ub_[v] &&
                        (c == G::not_vertex || ecc_lb_[v] < ecc_lb_[c]))
                    c = v;
                }
            }
            assert(c != G::not_vertex);

            WL e_c = sweep(c);

            if (ecc_lb_[c] == e_c) {
                if (directed) {
                    if ( ! prune) sweep(c, true);
                    else pruned_sweep_from_certifier(c, trav, zero_weight, e_c);
                }
                assert(improves_ecc_ub(e_c));
            } else {
                V p = trav.last_visited();
                WL e_p = sweep(p, true);
                assert(improves_ecc_lb());
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.all()", 1)
                    << nsweep <<" sweeps,"
                    <<" P:"<< P.size()
                    <<" todo: "<< all_ecc_todo
                    <<" / "<< n << std::endl;
        }
        if (mprune > 0) { ++nsweep; mprune = 0; ++nprune_sweep; }
        
        verb::cerr("eccentricity.all()", 1)
            << nsweep <<" sweeps ("<< nprune_sweep <<" pruned),   "
            << "R_ub=" << rad_ub <<" D_lb="<< diam_lb
            <<" P:"<< P.size() <<" C:" << C.size() << std::endl;

        all_lb_certif = opt_certif ? optim_lb_certif(P, max_weight) : P;
        all_ub_certif = C; // already optimal !
        rad_certif = opt_certif ? optim_lb_certif(P, rad_ub) : P;
        diam_certif = opt_certif ? diam_certif_when_lb_tight() : C;

        verb::cerr("eccentricity.all()", 1)
            << nsweep <<" sweeps ("<< nprune_sweep <<" pruned),   "
            << "R_cert:" << rad_certif.size()
            <<" D_cert:"<< diam_certif.size()
            <<" P:"<< all_lb_certif.size() << std::endl;
    }

    void all_from_lb_cert(const std::vector<V> & certif, int n_thd = 0,
                          bool opt_certif = true) {
        std::vector<V> last;
        
        for (V p : certif) {
            sweep(p, true);
            improves_ecc_lb();
            V l = trav.last_visited();
            if ( ! is_in(l, certif)) last.push_back(l);
            if (verb::progress())
                verb::cerr("eccentricity.all_from_lb_cert()", 1)
                        << nsweep <<" sweeps\n";
        }

        // pseudo-check
        verb::cerr()<< "nb last: "<< last.size() <<"\n";
        for (V u : last) {
            WL e = sweep(u); assert(e == ecc_lb_[u]);
        }

        // center:
        V c = G::not_vertex;
        for (V v = 0; v < n; ++v) {
            if (ecc_lb_[v] < ecc_ub_[v] &&
                (c == G::not_vertex || ecc_lb_[v] < ecc_lb_[c]))
                c = v;
        }
        assert(c != G::not_vertex);
        WL e_c = sweep(c);
        if (directed) pruned_sweep_from_certifier(c, trav, zero_weight, e_c);
        improves_ecc_ub(e_c);

        // random sampling for upper bounds and more checking
        for (int n_smpl=0, i=0; i < n && n_smpl < 3*(20+certif.size()); ++i) {
            V u = (127L * i + 1237L) % n;
            if (ecc_lb_[u] < ecc_ub_[u]) {
                ++n_smpl;
                V x = pruned_sweep_to_certifiers(u, trav, ecc_lb_[u], true);
                WL e_x = sweep(x);
                V p_x = trav.last_visited();
                if (directed)
                    pruned_sweep_from_certifier(x, trav, zero_weight, e_x);
                improves_ecc_ub(e_x);
                // pseudo-check:
                assert(e_x == ecc_lb_[x]); 
                if (( ! is_in(p_x, certif)) && ! is_in(p_x, last)) {
                    WL e_p = sweep(p_x); assert(e_p == ecc_lb_[p_x]);
                }
            }
        }
        // include in ub cert all uncovered nodes
        int n_uncov = 0;
        for (V v = 0; v < n; ++v) {
            if (ecc_lb_[v] < ecc_ub_[v]) {
                ++n_uncov;
                C.push_back(v);
                ecc_ub_[v] = ecc_lb_[v]; // assume certif is a tight-lb cert
            }
        }
        verb::cerr() << "n_uncov: "<< n_uncov <<" (eccs unchecked)\n";
        
        update_todo();
        all_lb_certif = P; // already optimized a priori
        all_ub_certif = optim_ub_certif_one_shot(true, true);
        // checking that certif is a tight-lb cert requires |all_ub_cert| sweeps
        rad_certif = opt_certif ? optim_lb_certif(P, rad_ub) : P;
        diam_certif = opt_certif ? (n_thd == 0 ? diam_certif_when_lb_tight()
                                    : diam_certif_when_lb_tight_threaded(n_thd))
            : all_ub_certif;
            
        verb::cerr("eccentricity.all_from_lb_cert()", 1)
            << nsweep <<" sweeps ("<< nprune_sweep <<" pruned),   "
            << "R_cert:" << rad_certif.size()
            <<" D_cert:"<< diam_certif.size()
            <<" lb_cert:"<< all_lb_certif.size()
            <<" ub_cert:"<< all_ub_certif.size() << std::endl;
    }

    
    void all_threaded(int n_thd, V start = G::not_vertex,
                      bool do_rad = false, bool opt_certif = true) {
        
        if (n_thd < 1) return all(start, opt_certif);

        if (do_rad) {
            radius(start, false, true);
        } else {
            if (start == G::not_vertex) start = max_degree_node();
            WL e_start = sweep(start);
            if (directed)
                pruned_sweep_from_certifier(start, trav, zero_weight, e_start);
            improves_ecc_ub(e_start);
            update_todo();
        }

        verb::cerr() << "all_ecc_todo: " << all_ecc_todo <<"\n";
        
        std::atomic<int> could_not_p(0);
                
     {
        const std::memory_order relaxed = std::memory_order_relaxed,
                                acquire = std::memory_order_acquire,
                                release = std::memory_order_release,
                                acq_rel = std::memory_order_acq_rel;
        
        //std::vector<std::atomic_flag> sweeping(n);
        //for (V u = 0; u < n; ++u) sweeping[u].clear();
        std::vector<std::atomic_bool> sweeping(n);
        for (V u = 0; u < n; ++u) sweeping[u].store(false, release);
        std::vector<std::atomic_bool> sweeping_P(n);
        for (V u = 0; u < n; ++u) sweeping_P[u].store(false, release);
        std::mutex this_mutex;
        // copies of ecc_lb and ecc_ub, thread i writes only its own copy
        std::vector<std::vector<std::atomic<WL>>> lb(n_thd), ub(n_thd);
        for (int i = 0; i < n_thd; ++i) {
            lb[i] = std::vector<std::atomic<WL>>(n);
            for (V u = 0; u < n; ++u) lb[i][u] = ecc_lb_[u];
            ub[i] = std::vector<std::atomic<WL>>(n);
            for (V u = 0; u < n; ++u) ub[i][u] = ecc_ub_[u];
        }

        auto go = [this, start, &sweeping, &sweeping_P, n_thd,
                   &this_mutex, &lb, &ub, &could_not_p,
                   acquire, release](const int i_thd) {
            
            traversal<G, WL, max_weight, zero_weight> trav(n);
            int64_t v_seed = 12347L * i_thd;
            
            int iter = 0, local_todo = all_ecc_todo;
            int nsweep_loc = 0, mprune_loc = 0, nprune_loc = 0;


            while (all_ecc_todo > 0) {
                ++iter;
                
                // Select a node c with minimal ecc_lb
                V c = G::not_vertex, c_locked = G::not_vertex;
                WL e_c_lb = max_weight;
                int n_sweeping = 0;
              
                for (int i = 0; i < n; ++i) {
                    V v = (v_seed + i /* * 127L */) % n;
                    if (lb[i_thd][v] < ub[i_thd][v]
                        && (c == G::not_vertex || lb[i_thd][v] < e_c_lb)) {
                        bool not_sweeping = false;
                        if (sweeping[v].compare_exchange_strong
                                       (not_sweeping, true, acq_rel, acquire)) {
                            V prev_c = c;
                            c = v;
                            e_c_lb = lb[i_thd][v].load(acquire);
                            if (prev_c != G::not_vertex)
                                sweeping[prev_c].store(false, release);
                        } else {
                            ++n_sweeping;
                            if (c_locked == G::not_vertex || v < c_locked)
                                c_locked = v;
                        }
                    }
                }

                WL e_m = max_weight;
                for (V v : graph)
                    if (lb[i_thd][v] < e_m && lb[i_thd][v] < ub[i_thd][v])
                        e_m = lb[i_thd][v];
                V vc = c == G::not_vertex ? 0 : c;
                
                auto sweep = [this, &trav,
                              &this_mutex](V u, bool backward = false) {
                    trav.clear();
                    const G &g = backward ? graph_rev : graph;
                    int nvis =
                      weighted ? trav.dijkstra(g, u) : trav.bfs(g, u);
                    V last = trav.last_visited();
                    V ecc = trav.dist(last);
                    this_mutex.lock();
                    ++nsweep;
                    if (ecc < rad_ub && ! (backward && directed))
                        { rad_ub = ecc; rad_node = u; }
                    if (ecc > diam_lb)
                        { diam_lb = ecc; diam_node = backward ? last : u; }
                    this_mutex.unlock();
                    return ecc;
                };

                // Get other thread bounds
                auto update_bounds = [this, &lb, &ub, i_thd, n_thd,
                                     &local_todo, acquire, release]() {
                    local_todo = 0;                    
                    for (V v = 0; v < n; ++v) {
                        WL l0 = lb[i_thd][v].load(acquire);
                        WL l = l0;
                        for (int j = 0; j < n_thd; ++j) {
                            if (j != i_thd)
                                l = std::max(l, lb[j][v].load(acquire));
                        }
                        if (l > l0)
                            lb[i_thd][v].store(l, release);
                        
                        WL u0 = ub[i_thd][v].load(acquire);
                        WL u = u0;
                        for (int j = 0; j < n_thd; ++j) {
                            if (j != i_thd)
                                u = std::min(u, ub[j][v].load(acquire));
                        }
                        if (u < u0) ub[i_thd][v].store(u, release);
                        if (l < u) ++local_todo;
                    }
                    if (i_thd == 0) all_ecc_todo = local_todo;
                };

                bool force_update = false;

                if (c == G::not_vertex) {
                    if (iter % 500 == 50 * i_thd) {
                        verb::cerr() << i_thd <<": no c found, "
                                     << n_sweeping << " locked as "
                                     << c_locked <<", todo: "
                                     << local_todo <<"  ";
                        if (c_locked >= 0) 
                            verb::cerr() << lb[i_thd][c_locked]
                                         << " <= "
                                         << sweep(c_locked)
                                         << " <= "
                                         << ub[i_thd][c_locked]
                                         <<"\n";
                    }
                    update_bounds();
                    continue;
                }
                
                if (i_thd == 0 && verb::progress())
                    verb::cerr("eccentricity.all_threaded()", 1)
                        << (iter * n_thd) <<" sweeps,"
                        <<" P:"<< P.size()
                        <<" todo: "<< local_todo <<" / "<< n << "\n";

                
                WL e_c = sweep(c);
                
                auto improve_ecc = [this, &trav, &lb, &ub, i_thd,
                                    acquire, release, &this_mutex]
                                   (bool do_lb, bool do_ub, WL ecc, int scc) {
                    if (do_lb) {
                        bool improves = false;
                        for (int i = 0; i < trav.nvis(); ++i) {
                            V v = trav.visit(i);
                            WL l = trav.dist(v);
                            if (l > lb[i_thd][v].load(acquire)) {
                                lb[i_thd][v].store(l, release);
                                improves = true;
                            }
                        }
                        if (improves) {
                            this_mutex.lock();
                            last_lb_improve = nsweep; // race
                            this_mutex.unlock();
                        }
                    }
                    if (do_ub) {
                        for (int i = 0; i < trav.nvis(); ++i) {
                            V v = trav.visit(i);
                            WL u = trav.dist(v) + ecc;
                            if (scc_nb[v] == scc
                                && u < ub[i_thd][v].load(acquire))
                                ub[i_thd][v].store(u, release);
                        }
                    }
                };

                if (e_c == e_c_lb) {
                    if (directed)
                        pruned_sweep_from_certifier_threaded
                            (c, trav, lb[i_thd],
                             nsweep_loc, mprune_loc, nprune_loc,
                             zero_weight, e_c);
                    if (nprune_loc > 0) {
                        this_mutex.lock();
                        nprune_sweep += nprune_loc;
                        nprune_loc = 0;
                        this_mutex.unlock();
                    }
                    improve_ecc(false, true, e_c, scc_nb[c]); // impr upper bnds
                    
                    this_mutex.lock();
                    C.push_back(c);
                    if (nsweep_loc > 0) { nsweep += nsweep_loc; nsweep_loc = 0;}
                    if (nprune_loc > 0)
                        { nprune_sweep += nprune_loc; nprune_loc = 0;}
                    this_mutex.unlock();
                } else {                    
                    V p = trav.last_visited();
                    std::vector<V> cob;
                    for (int i = trav.nvis() - 1; i >= 0; --i) {
                        V v = trav.visit(i);
                        if (trav.dist(v) < e_c) break;
                        // else v in coball ovB(c,e_c)
                        cob.push_back(v);
                    }
                    this_mutex.lock();
                    Pcoballs.push_back(cob);
                    this_mutex.unlock();
                    bool not_sweeping = false;
                    if (sweeping_P[p].compare_exchange_strong
                                      (not_sweeping, true, acq_rel, acquire)) {
                        verb::cerr() << "thrd "<< i_thd <<": "
                                     << "  p: " << p << std::endl;
                        WL e_p = sweep(p, true);
                        improve_ecc(true, false, max_weight, scc_nb[p]); //l.bnd

                        this_mutex.lock();
                        P.push_back(p);
                        this_mutex.unlock();
                    } else {
                        could_not_p.fetch_add(1);
                        verb::cerr() <<i_thd <<": could not p on " << p
                                     <<" c=" << c <<" e_c="<< e_c <<"\n";
                        force_update = true;
                    }

                    // Here to prevent others from doing c before bound update
                    sweeping[c].store(false, release);
                }
                // following assert may fail if interrupted here
                assert((! sweeping[c].load(acquire))
                       || lb[i_thd][c].load(acquire) == e_c);

                if (force_update || iter % 4 == i_thd % 4) 
                    update_bounds();
                /*
                    else { // update from 0
                        for (V v = 0; v < n; ++v) {
                            WL l0 = lb[i_thd][v];
                            WL l = lb[0][v];
                            if (l > l0) lb[i_thd][v] = l; 
                            WL u0 = ub[i_thd][v];
                            WL u = ub[0][v];
                            if (u < u0) ub[i_thd][v] = u;
                        }                    
                    }
                */
            } // end while
            
            this_mutex.lock();
            if (nsweep_loc > 0) { nsweep += nsweep_loc; nsweep_loc = 0;}
            if (nprune_loc > 0)
                { nprune_sweep += nprune_loc; nprune_loc = 0;}
            if (mprune_loc > 0) { ++nsweep; ++nprune_sweep; }
            this_mutex.unlock();            
        };

        std::vector<std::thread> threads(n_thd);
        for (int i = 0; i < n_thd; ++i) {
            threads[i] = std::thread(go, i);
        }
        for (int i = 0; i < n_thd; ++i) threads[i].join();

        for (V u = 0; u < n; ++u) {
            ecc_lb_[u] = lb[0][u].load(acquire);
            ecc_ub_[u] = ub[0][u].load(acquire);
            if (ecc_lb_[u] != ecc_ub_[u])
                std::cerr <<"u="<< u <<" lb="<< ecc_lb_[u]
                             <<" ub="<< ecc_ub_[u] <<"\n";
            assert(ecc_lb_[u] == ecc_ub_[u]);
        }
        verb::cerr("all_threaded()", 1) <<"checked !";

    }

        update_todo(); assert(all_ecc_todo == 0);

        verb::cerr("eccentricity.all_threaded()", 1)
            << nsweep <<" sweeps ("<< could_not_p <<" useless parallel), "
            << "R_ub=" << rad_ub <<" D_lb="<< diam_lb
            <<" P:"<< P.size() <<" C:" << C.size() << std::endl;

        all_lb_certif = opt_certif ? optim_lb_certif(P, max_weight) : P;
        verb::cerr() << "all_lb_certif "<< all_lb_certif.size() <<"\n";
        all_ub_certif = C; // already optimal if 1 thread
        rad_certif = opt_certif ? optim_lb_certif(P, rad_ub) : P;
        verb::cerr() << "rad_certif "<< rad_certif.size() <<"\n";
        diam_certif = opt_certif ? diam_certif_when_lb_tight_threaded(n_thd)
                                 : C;

        verb::cerr("eccentricity.all_threaded()", 1)
            << nsweep <<" sweeps ("<< nprune_sweep <<" pruned), "
            << "R_cert:" << rad_certif.size()
            <<" D_cert:"<< diam_certif.size()
            <<" P:"<< all_lb_certif.size() << std::endl;
    }

    int center_on_path_to_lb(V u, WL ub_goal = zero_weight) {
        ub_goal = std::max(ecc_lb_[u], ub_goal);
        V c = u, p = next_hop_to_lb[u];
        WL d = zero_weight;
        while (p != c && scc_nb[p] == scc_nb[u] &&
               d + next_hop_dist[c] + ecc_lb_[p] <= ub_goal) {
            d += next_hop_dist[c];
            c = p;
            p = next_hop_to_lb[c];
        }
        return c;
    }

    int center_on_path_to_lb1(V u, WL ub_goal = zero_weight) {
        ub_goal = std::max(ecc_lb_[u], ub_goal);
        V c = u, p = next_hop_to_lb[u];
        WL d = zero_weight;
        std::vector<V> stack = {c};
        while (p != c && scc_nb[p] == scc_nb[u]
               && d + next_hop_dist[c] + ecc_lb_[p] <= ub_goal
               /*&& ecc_lb_[p] < ecc_lb_[c]*/) {
            d += next_hop_dist[c];
            c = p;
            stack.push_back(c);
            p = next_hop_to_lb[c];
        }
        c = u;
        for (int i = 0; i < stack.size(); ++i) {
            V v = stack[i];
            if (ecc_lb_[v] <= ecc_lb_[c]) c = v;
        }
        return c;
    }

    /*
    WL diameter_a_star(V start = G::not_vertex, start_diam do_rad = SUMSWEEP) {
        
        if (do_rad == RADIUS || rad_todo == 0) {
            radius(start, true);
            sweep(rad_node);
        } else if (do_rad == SUMSWEEP) {
            sweep(sum_sweep(start));
        } else {
            sweep(start == G::not_vertex ? max_degree_node() : start);
        }
        improves_ecc_ub(diam_lb);
        update_todo();

        int iter = 0;
        while (diam_todo > 0) {
            ++iter;
            
            V u = G::not_vertex;
            for (V v = 0; v < n; ++v) {
                if (ecc_ub_[v] > diam_lb
                    && (u == G::not_vertex || ecc_ub_[v] > ecc_ub_[u]))
                    u = v;
            }
            assert(u != G::not_vertex);
            V c = pruned_sweep_to_certifiers(u, trav, diam_lb);

            int n_far = 0;
            for (V v = 0; v < n; ++v) {
                if (//ecc_ub_[v] > ecc_lb_[c] &&
                    dist_ub(u, v) > diam_lb)
                    ++n_far;
            }
            bool c_certif_of_u = false;
            if (iter >= 10 && n_far < 100) {
                c_certif_of_u = true;
                int loc_prune = 0, n_loc = 0;
                V lst = u;
                for (V v = 0; v < n; ++v) {
                    if (//ecc_ub_[v] > ecc_lb_[c] &&
                        dist_ub(u, v) > diam_lb) {
                        lst = v;
                        ++n_loc;
                        auto d_c_lb = [this, c](V u) { return dist_lb(u, c); };
                        WL e_c = ecc_lb_[c];
                        auto filter = [this, c, e_c, &loc_prune]
                            (V p, WL dp, V u, WL du) {
                            ++mprune; ++loc_prune;
                            return du + dist_lb(u,c) <= e_c;
                        };
                        last_sweep_node = G::not_vertex;
                        trav.clear();
                        trav.a_star(graph, v, c, d_c_lb, filter);
                        WL u_c = trav.dist(c);
                        trav.clear_a_star(graph);
                        if (u_c > ecc_lb_[c]) { c_certif_of_u = false; break; }
                        if (mprune > graph.m()) {
                            mprune -= graph.m()+1; ++nsweep; ++nprune_sweep;
                        }
                    }
                }
                if (n_loc > 1)
                verb::cerr() <<"u="<< u <<"  c="<< c
                          << "  lpn: " << c_certif_of_u << "  nloc="<< n_loc
                          <<" lst=" << lst
                          <<" avg:"<< loc_prune / (n_loc == 0 ? 1 : n_loc)<<"   "
                          << loc_prune <<" "<< graph.m()
                          << " : " << (loc_prune * 100 / graph.m()) <<" %\n";;
            }
            

            if (c_certif_of_u) {
                ecc_ub_[u] = diam_lb;
                if (is_in(c, C)) verb::cerr() <<"."; else C.push_back(c);
            } else {

                WL e_c = sweep(c);

                // If the lower bound for u is tight, either c is a certif of u
                //   or the lower bound for c is untight.
                if (ecc_lb_[c] == e_c) {
                    assert(improves_ecc_ub());
                } else {
                    V p = trav.last_visited();
                    verb::cerr() <<"  p="<< p << std::endl;
                    WL e_p = sweep(p);
                    assert(improves_ecc_lb());
                }
            }

            update_todo();
        }
        if (mprune > 0) { ++nsweep; mprune = 0; ++nprune_sweep; }
        verb::cerr("eccentricity.diameter_a_star()", 1)
            << nsweep <<" sweeps ("<< nprune_sweep <<" pruned),"
            << " R_ub=" << rad_ub <<" D_lb="<< diam_lb
            <<" P:"<< P.size() <<" C:" << C.size() << std::endl;

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }
    */

    WL diameter(V start = G::not_vertex, start_diam do_rad = SUMSWEEP,
                bool choose_ub = true, int loose_ineq = 0) {
        diameter_start(start, do_rad);

        std::vector<std::vector<V>> certif_for(n);
        std::vector<V> todo;

        const int prunePeriod = 8, pruneMax = 3, pruneMaxFail = 3;
        int iter = 0, lastPrunes = 0;
        while (diam_todo > 0) {
            ++iter;

            // do random sampling and choose certifiers that cover the sample
            if (iter >= lastPrunes + prunePeriod) {
                int ns = nsweep;
                for (V c = 0; c < n; ++c) certif_for[c].clear();
                todo.clear();
                V vi = 0;
                for ( ; vi < n && nsweep < ns + pruneMax; ++vi) {
                    V v = (127L * iter + 1237L * vi) % n; // pseudo rnd order
                    if (ecc_ub_[v] > diam_lb) {
                        todo.push_back(v);
                int f_e = 2 - loose_ineq, f_D = loose_ineq;
                WL e_D = (f_e * ecc_lb_[v] + f_D * diam_lb) / 2;
                        pruned_sweep_to_certifiers(v, trav, e_D);
                        for (int i = trav.nvis() - 1; i >= 0; --i) {
                            V c = trav.visit(i);
                            certif_for[c].push_back(v);
                        }
                    }
                }
                std::vector<V> cert = gdy_set_cov(certif_for, n, todo);
                verb::cerr() << iter << "  cert: " << cert.size()
                          << "  todo: " << todo.size()
                          <<"  nsweeps: "<< (nsweep - ns) <<"\n";
                lastPrunes = iter;

                int c_done = 0, c_fail = 0;
                for(V c : cert) {
                    WL e_c = sweep(c);
                    ++c_done;
                    if (ecc_lb_[c] == e_c) {
                        if (directed)
                            pruned_sweep_from_certifier(c, trav, diam_lb, e_c);
                        improves_ecc_ub(e_c, diam_lb);
                    } else {
                        V p = trav.last_visited();
                        sweep(p, true);
                        assert(improves_ecc_lb());
                        ++c_fail;
                        if(c_fail >= pruneMaxFail) break;
                    }
                }
                verb::cerr() <<"  nsweep: " << nsweep
                          <<"  c_done: " << c_done <<"  c_fail: " << c_fail
                          << std::endl;
            }
            
            V u = G::not_vertex;
            for (V v = 0; v < n; ++v) {
                if (ecc_ub_[v] > diam_lb
                    && (u == G::not_vertex ||
                        (choose_ub ? ecc_ub_[v] > ecc_ub_[u]
                                   : ecc_lb_[v] > ecc_lb_[u])))
                    u = v;
            }
            if (u == G::not_vertex) { update_todo(); break; }
            V c = pruned_sweep_to_certifiers(u, trav, ecc_lb_[u], true); 
            
            WL e_c = sweep(c);

            // If the lower bound for u is tight, either c is a certifier for u
            //   or the lower bound for c is untight.
            if (ecc_lb_[c] == e_c) {
                if (directed)
                    pruned_sweep_from_certifier(c, trav, diam_lb, e_c);
                assert(improves_ecc_ub(e_c, diam_lb));
            } else {
                V p = trav.last_visited();
                WL e_p = sweep(p, true);
                assert(improves_ecc_lb());
                lastPrunes = iter; // postpone prun as long as low. bnds improve
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }
        diam_nsweep += nprune_sweep; // account set cover computations
        if (mprune > 0) { ++nsweep; mprune = 0; ++nprune_sweep; }
        nsweep += nprune_sweep; // account set cover computations
        verb::cerr("eccentricity.diameter()", 1)
            << nsweep <<" sweeps ("<< nprune_sweep <<" pruned),"
            << " R_ub=" << rad_ub <<" D_lb="<< diam_lb
            <<" P:"<< P.size() <<" C:" << C.size() << std::endl;

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }

    /*
    WL diameter_sweep(V start = G::not_vertex, bool do_rad = true, bool choose_ub = true) {
        
        if (do_rad) {
            radius(start, true);
        } else {
            if (start == G::not_vertex) start = max_degree_node();
            sweep(start);
            improves_ecc_ub();
        }

        std::vector<std::vector<V>> certif_for(n);
        std::vector<V> todo;

        const int prunePeriod = 8, pruneMax = 3, pruneMaxFail = 3;
        int iter = 0, lastPrunes = 0;
        while (diam_todo > 0) {
            ++iter;

            // do random sampling and choose certifiers that cover the sample
            if (iter >= lastPrunes + prunePeriod) {
                int ns = nsweep;
                for (V c = 0; c < n; ++c) certif_for[c].clear();
                todo.clear();
                V vi = 0;
                for ( ; vi < n && nsweep < ns + pruneMax; ++vi) {
                    V v = (127L * iter + 1237L * vi) % n; // pseudo rnd order
                    if (ecc_ub_[v] > diam_lb) {
                        todo.push_back(v);
                        pruned_sweep_to_certifiers(v, trav, diam_lb);
                        for (int i = trav.nvis() - 1; i >= 0; --i) {
                            V c = trav.visit(i);
                            certif_for[c].push_back(v);
                        }
                    }
                }
                std::vector<V> cert = gdy_set_cov(certif_for, n, todo);
                verb::cerr() << iter << "  cert: " << cert.size()
                          << "  todo: " << todo.size()
                          <<"  nsweeps: "<< (nsweep - ns) <<"\n";
                lastPrunes = iter;

                int c_done = 0, c_fail = 0;
                for(V c : cert) {
                    WL e_c = sweep(c);
                    ++c_done;
                    if (ecc_lb_[c] == e_c) improves_ecc_ub();
                    else {
                        V p = trav.last_visited();
                        sweep(p);
                        assert(improves_ecc_lb());
                        ++c_fail;
                        if(c_fail >= pruneMaxFail) break;
                    }
                }
                verb::cerr() <<"  nsweep: " << nsweep
                          <<"  c_done: " << c_done <<"  c_fail: " << c_fail
                          << std::endl;
            }
            
            V u = G::not_vertex;
            for (V v = 0; v < n; ++v) {
                if (ecc_ub_[v] > diam_lb
                    && (u == G::not_vertex ||
                        (choose_ub ? ecc_ub_[v] > ecc_ub_[u]
                                   : ecc_lb_[v] > ecc_lb_[u])))
                    u = v;
            }
            if (u == G::not_vertex) { update_todo(); break; }
            
            WL e_u = sweep(u);
            if (ecc_lb_[u] < e_u) {
                V p = trav.last_visited();
                WL e_p = sweep(p);
                assert(improves_ecc_lb());
                update_todo();
                continue;
            }
            
            // Find furthest certif for p
            V c = u; //furthest_certif_for_root(diam_lb);
            for (int i = 0; i < trav.nvis() -1; ++i) {
                V v = trav.visit(i);
                if (trav.dist(v) + ecc_lb_[v] <= diam_lb
                    )// &&  ecc_lb_[v] < ecc_lb_[c])
                    c = v;
            }
            
            WL e_c = sweep(c);

            // If the lower bound for u is tight, either c is a certifier for u
            //   or the lower bound for c is untight.
            if (ecc_lb_[c] == e_c) {
                assert(improves_ecc_ub(diam_lb));
            } else {
                V p = trav.last_visited();
                WL e_p = sweep(p);
                assert(improves_ecc_lb());
                lastPrunes = iter; // postpone prun as long as low. bnds improve
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter_sweep()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }
        diam_nsweep += nprune_sweep; // account set cover computations
        if (mprune > 0) { ++nsweep; mprune = 0; ++nprune_sweep; }
        nsweep += nprune_sweep; // account set cover computations
        verb::cerr("eccentricity.diameter_sweep()", 1)
            << nsweep <<" sweeps ("<< nprune_sweep <<" pruned),"
            << " R_ub=" << rad_ub <<" D_lb="<< diam_lb
            <<" P:"<< P.size() <<" C:" << C.size() << std::endl;

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }
    */

    WL diameter_path_to_lb(V start = G::not_vertex,
                           start_diam do_rad = SUMSWEEP,
                           bool choose_ub = true, bool version1 = false,
                           int loose_ineq = 0) {
        diameter_start(start, do_rad);

        while (diam_todo > 0) {

            V u = G::not_vertex;
            for (V v = 0; v < n; ++v) {
                if (ecc_ub_[v] > diam_lb
                    && (u == G::not_vertex ||
                        (choose_ub ? ecc_ub_[v] > ecc_ub_[u]
                                   : ecc_lb_[v] > ecc_lb_[u])))
                    u = v;
            }
            assert(u != G::not_vertex);
                int f_e = 2 - loose_ineq, f_D = loose_ineq;
                WL e_D = (f_e * ecc_lb_[u] + f_D * diam_lb) / 2;
            V c = version1 ? center_on_path_to_lb1(u, e_D)
                           : center_on_path_to_lb(u, e_D);
            
            WL e_c = sweep(c);

            // If the lower bound for u is tight, either c is a certifier for u
            //   or the lower bound for c is untight.
            if (ecc_lb_[c] == e_c) {
                if (directed)
                    pruned_sweep_from_certifier(c, trav, diam_lb, e_c);
                assert(improves_ecc_ub(e_c, diam_lb));
            } else {
                V p = trav.last_visited();
                WL e_p = sweep(p, true);
                assert(improves_ecc_lb());
            }

            diameter_update("diameter_path_to_lb");
        }

        verb::cerr("eccentricity.diameter_path_to_lb()", 1)
            << nsweep <<" sweeps,"
            << " R_ub=" << rad_ub <<" D_lb="<< diam_lb
            <<" P:"<< P.size() <<" C:" << C.size() << std::endl;

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }

    /*
    WL diameter_c_to_lb(V start = G::not_vertex, bool do_rad = true, bool choose_ub = true) {
        
        if (do_rad) {
            radius(start, true);
        } else {
            if (start == G::not_vertex) start = max_degree_node();
            sweep(start);
            improves_ecc_ub();
        }

        int resweep = 0;
        
        while (diam_todo > 0) {

            V u = G::not_vertex;
            for (V v = 0; v < n; ++v) {
                if (ecc_ub_[v] > diam_lb
                    && (u == G::not_vertex ||
                        (choose_ub ? ecc_ub_[v] > ecc_ub_[u]
                                   : ecc_lb_[v] > ecc_lb_[u])))
                    u = v;
            }
            assert(u != G::not_vertex);
            V c = diam_c_to_lb[u];
            
            WL e_c = sweep(c);

            if (verb::progress())
                verb::cerr("diam",1)
                    <<"e_c="<< e_c <<" todo: "<< diam_todo <<"\n";

            improves_ecc_ub();
            V p = trav.last_visited();
            if (ecc_ub_[u] > diam_lb) {
                ++resweep; //std::cerr << "resweep lb_node\n";
                sweep(lb_node[u]);
                improves_ecc_lb();
            }
            if (ecc_lb_[c] == e_c) {
                //assert(improves_ecc_ub());
            } else {
                WL e_p = sweep(p);
                assert(improves_ecc_lb());
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter_c_to_lb()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }
        verb::cerr("eccentricity.diameter_c_to_lb()", 1)
            << nsweep <<" sweeps ("<< resweep <<" resweeps),"
            << " R_ub=" << rad_ub <<" D_lb="<< diam_lb
            <<" P:"<< P.size() <<" C:" << C.size() << std::endl;

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }
    */
    

    WL diameter_prune(V start = G::not_vertex, start_diam do_rad = SUMSWEEP,
            bool choose_ub = true, bool decr_ecc = true, int loose_ineq = 0) {
        
        diameter_start(start, do_rad);

        while (diam_todo > 0) {

            V u = G::not_vertex;
            for (V v = 0; v < n; ++v) {
                if (ecc_ub_[v] > diam_lb
                    && (u == G::not_vertex ||
                        (choose_ub ? ecc_ub_[v] > ecc_ub_[u]
                                   : ecc_lb_[v] > ecc_lb_[u])))
                    u = v;
            }
            assert(u != G::not_vertex);
                int f_e = 2 - loose_ineq, f_D = loose_ineq;
                WL e_D = (f_e * ecc_lb_[u] + f_D * diam_lb) / 2;
            V c = pruned_sweep_to_certifiers(u, trav, e_D, decr_ecc);
            
            WL e_c = sweep(c);

            // If the lower bound for u is tight, either c is a certifier for u
            //   or the lower bound for c is untight.
            if (ecc_lb_[c] == e_c) {
                if (directed)
                    pruned_sweep_from_certifier(c, trav, diam_lb, e_c);
                assert(improves_ecc_ub(e_c, diam_lb));
            } else {
                V p = trav.last_visited();
                WL e_p = sweep(p, true);
                assert(improves_ecc_lb());
            }

            diameter_update("diameter_prune");
        }
        if (mprune > 0) { ++nsweep; mprune = 0; ++nprune_sweep; }
        verb::cerr("eccentricity.diameter_prune()", 1)
            << nsweep <<" sweeps ("<< nprune_sweep <<" pruned),"
            << " R_ub=" << rad_ub <<" D_lb="<< diam_lb
            <<" P:"<< P.size() <<" C:" << C.size() << std::endl;

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }

    /*
    WL diameter_prune_diff(V start = G::not_vertex, bool do_rad = true, bool choose_ub = true, bool decr_ecc = false) {
        
        if (do_rad) {
            radius(start, true);
        } else {
            if (start == G::not_vertex) start = max_degree_node();
            sweep(start);
            improves_ecc_ub();
        }

        while (diam_todo > 0) {

            V u = G::not_vertex;
            for (V v = 0; v < n; ++v) {
                if (ecc_ub_[v] > diam_lb
                    && (u == G::not_vertex ||
                        (choose_ub ? ecc_ub_[v] > ecc_ub_[u]
                                   : ecc_lb_[v] > ecc_lb_[u])))
                    u = v;
            }
            assert(u != G::not_vertex);
            V c = pruned_sweep_to_certifiers(u, trav, diam_lb, decr_ecc);
            for (int i = trav.nvis()-1; i >= 0; --i) {
                int v = trav.visit(i);
                if (ecc_ub_[v] - ecc_lb_[v] > ecc_ub_[v] - ecc_lb_[c])
                    c = v;
            }
            
            WL e_c = sweep(c);

            // If the lower bound for u is tight, either c is a certifier for u
            //   or the lower bound for c is untight.
            if (ecc_lb_[c] == e_c) {
                assert(improves_ecc_ub());
            } else {
                V p = trav.last_visited();
                WL e_p = sweep(p);
                assert(improves_ecc_lb());
            }

            update_todo();
            if (verb::progress())
                verb::cerr("eccentricity.diameter_prune_diff()", 1)
                    << nsweep <<" sweeps, todo: "<< diam_todo <<" / "<< n<<"\n";
        }
        if (mprune > 0) { ++nsweep; mprune = 0; ++nprune_sweep; }
        verb::cerr("eccentricity.diameter_prune()", 1)
            << nsweep <<" sweeps ("<< nprune_sweep <<" pruned),"
            << " R_ub=" << rad_ub <<" D_lb="<< diam_lb
            <<" P:"<< P.size() <<" C:" << C.size() << std::endl;

        diam_certif = optim_ub_certif(C, diam_lb);

        return diam_lb;
    }
    */


    // Find nodes that c certifies
    V pruned_sweep_from_certifier(V c,  traversal<G> &trav,
                         WL ub_goal = zero_weight, WL e_c = zero_weight) {
        e_c = std::max(e_c, ecc_lb_[c]);
        int scc_c = scc_nb[c];
        auto filter = [this, c, e_c, scc_c,
                       ub_goal](V v, WL d, V par, WL dpar) {
            ++mprune; // count number of visited edges
            return scc_nb[v] == scc_c
                     && d + e_c <= std::max(ecc_lb_[v], ub_goal); 
        };
        last_sweep_node = G::not_vertex;
        trav.clear();
        if (weighted) trav.dijkstra(graph_rev, c, filter);
        else trav.bfs(graph_rev, c, filter);
        if (mprune > graph.m()) {
            mprune -= graph.m()+1; ++nsweep; ++nprune_sweep;
        }
        return trav.last_visited();
    }

    // Find nodes that c certifies
    V pruned_sweep_from_certifier_threaded(V c,  traversal<G> &travloc,
                         std::vector<std::atomic<WL>> &lb,
                         int &nsweep, int &mprune, int &nprune_sweep,
                         WL ub_goal = zero_weight, WL e_c = zero_weight) {
        const std::memory_order acquire = std::memory_order_acquire;
        e_c = std::max(e_c, lb[c].load(acquire));
        int scc_c = scc_nb[c];
        auto filter = [this, c, e_c, scc_c, &nsweep, &mprune, &nprune_sweep,
                       &lb, ub_goal, acquire](V v, WL d, V par, WL dpar) {
            ++mprune; // count number of visited edges
            return scc_nb[v] == scc_c
                     && d + e_c <= std::max(lb[v].load(acquire), ub_goal); 
        };
        last_sweep_node = G::not_vertex;
        travloc.clear();
        if (weighted) travloc.dijkstra(graph_rev, c, filter);
        else travloc.bfs(graph_rev, c, filter);
        if (mprune > graph.m()) {
            mprune -= graph.m()+1; ++nsweep; ++nprune_sweep;
        }
        return travloc.last_visited();
    }

    // returns furthest certifier found
    V pruned_sweep_to_certifiers(V u, traversal<G> &trav,
                               WL ub_goal = zero_weight, bool decr_ecc = false,
                               WL e_u = zero_weight, WL e_c = zero_weight) {
        ub_goal = std::max(std::max(ecc_lb_[u], e_u), ub_goal);
        int scc_u = scc_nb[u];
        auto filter = [this, ub_goal, u, scc_u,
                       decr_ecc](V v, WL d, V par, WL dpar) {
            ++mprune; // count number of visited edges
            return (decr_ecc ? (ecc_lb_[v] <= ecc_lb_[par]) : true)
                   && scc_nb[v] == scc_u
                   && d + ecc_lb_[v] <= ub_goal;
        };
        last_sweep_node = G::not_vertex;
        trav.clear();
        if (weighted) trav.dijkstra(graph, u, filter);
        else trav.bfs(graph, u, filter);
        if (mprune > graph.m()) {
            mprune -= graph.m()+1; ++nsweep; ++nprune_sweep;
        }
        return trav.last_visited();
    }

    // returns furthest certifier found
    V pruned_sweep_to_certifiers_threaded(V u, traversal<G> &trav,
                                 int &nsweep, int &mprune, int &nprune_sweep,
                                 WL ub_goal = zero_weight, double dfact = 1.0) {
        ub_goal = std::max(ecc_lb_[u], ub_goal);
        int scc_u = scc_nb[u];
        auto filter = [this, ub_goal, u, scc_u,
                       &nsweep, &mprune, &nprune_sweep, dfact]
               (V v, WL d, V par, WL dpar) {
            ++mprune; // count number of visited edges
            return scc_nb[v] == scc_u && dfact * d + ecc_lb_[v] <= ub_goal;
        };
        last_sweep_node = G::not_vertex;
        trav.clear();
        if (weighted) trav.dijkstra(graph, u, filter);
        else trav.bfs(graph, u, filter);
        if (mprune > graph.m()) {
            mprune -= graph.m()+1; ++nsweep; ++nprune_sweep;
        }
        return trav.last_visited();
    }

    std::vector<V> certifiers(V u, bool check = true) {
        pruned_sweep_to_certifiers(u, trav);
        std::vector<V> cert;
        for (int i = 0; i < trav.nvis(); ++i) {
            V c = trav.visit(i);
            // to be sure c is an ub cert for u:
            assert( (! check) || ecc_lb_[c] == ecc_ub_[c]);
            cert.push_back(c);
        }
        return cert;
    }
    
    void all_set_cov(V start = G::not_vertex) {
        if (rad_todo > 0) radius(start, false, true);

        // heur. to improve all lower bounds, and upper bounds of central nodes
        all_basic(rad_nsweep, 5*rad_nsweep); 

        optim_ub_when_lb_tight(true);

        all_lb_certif = optim_lb_certif(P, max_weight);
        all_ub_certif = optim_ub_certif(C, zero_weight);
    }

    /** Basic heuristic for computing all eccentricities :
     * Take a node with distinct bounds having smallest lower bound
     * and add certificates for that node (both upper and lower).
     * The choice of nodes with small lower bound helps finding
     * good upper bound certificates and new lower bound certificates.
     * This appears to be efficient in graphs where very few vertices
     * happen to be last visited in a sweep.
    */     
    void all_basic(int break_no_lb_improved = INT_MAX,
                   int sweep_budget = INT_MAX,
                   bool all_ecc = true) {
        int iter = 0;
        while ((all_ecc ? all_ecc_todo > 0 : diam_todo > 0)
               && --sweep_budget >= 0) {
            ++iter;
            
            // Node with smallest untight lower bound
            V s = G::not_vertex;
            for (V u = 0; u < n; ++u) {
                if (ecc_lb_[u] < ecc_ub_[u]
                         && (s == G::not_vertex || ecc_lb_[u] < ecc_lb_[s]))
                    s = u;
            }
            WL e_s = sweep(s);
            V p = trav.last_visited();
            V c = furthest_certif_for_root();
            
            // add a certificate for lower bound
            if (ecc_lb_[s] < e_s) { 
                sweep(p, true);
                improves_ecc_lb();
            }

            // add a certificate for upper bound
            if (ecc_ub_[s] > e_s) { 
                WL e_c = sweep(c);
                sweep(c, true);
                improves_ecc_ub(e_c);
            }
            
            update_todo();
            if (nsweep - last_lb_improve > break_no_lb_improved) break;
        }        
    }

    void clear() {
        clear_ecc_lb();
        clear_ecc_ub();
        clear_sum();
        for (V u = 0; u < n; ++u) {
            dist_lab_lb[u].clear();
            dist_lab_ub[u].clear();
            dist_C[u] = max_weight;
            save_dist[u] = max_weight;
            sample_count[u] = 0;
        }
        rad_ub = max_weight; diam_lb = zero_weight;
        rad_node = G::not_vertex; diam_node = G::not_vertex;
        last_sweep_node = G::not_vertex;
        nsweep = 0; mprune = 0; nprune_sweep = 0;
        rad_nsweep = not_nsweep; diam_nsweep = not_nsweep;
        all_ecc_nsweep = not_nsweep;
        rad_todo = n; diam_todo = n; all_ecc_todo = n;
        last_lb_improve = 0; last_ub_improve = 0;
        P.clear(); C.clear(); Pcoballs.clear(); Rcoballs.clear();
        rad_certif.clear(); diam_certif.clear();
        all_lb_certif.clear(); all_ub_certif.clear();
        sample.clear();
    }
    
    std::vector<V> optim_lb_certif(const std::vector<V> & P,
                                   WL lb_required = max_weight,
                                   bool optim_lb = false) {
        if (lb_required == zero_weight
            || (diam_todo == 0 && diam_lb == zero_weight)) { // all eccs are 0
            return {};
        }
        std::vector<std::vector<V>> certif_for(n);
        last_sweep_node = G::not_vertex;
        for (V p : P) {
            trav.clear();
            if (weighted) trav.dijkstra(graph_rev, p);
            else trav.bfs(graph_rev, p);
            ++nsweep;
            WL e_p = trav.dist(trav.last_visited());
            for (V u = 0; u < n; ++u) {
                WL lb = optim_lb ? std::max(trav.dist(u), e_p - trav.dist(u))
                                 : trav.dist(u);
                if (lb >= std::min(lb_required, ecc_ub_[u]))
                    certif_for[p].push_back(u);
            }
        }
        if (  ! strongly_connected) {
            std::vector<V> tocov;
            for (V u : graph) if (scc_nb[u] == scc_largest) tocov.push_back(u);
            return gdy_set_cov(certif_for, n, tocov);
        } else {
            return gdy_set_cov(certif_for, n);
        }
    }
    
    std::vector<V> optim_ub_certif(const std::vector<V> & C,
                                   WL ub_required = zero_weight,
                                   bool optim_ub = false) {
        std::vector<std::vector<V>> certif_for(n);
        last_sweep_node = G::not_vertex;
        for (V c : C) {
            WL e_c = ecc_ub_[c]; // since c is in C
            if (optim_ub) {
                assert( ! directed);
                trav.clear();
                if (weighted) trav.dijkstra(graph, c); else trav.bfs(graph, c);
                ++nsweep;
                trav.tree_eccentricities();
            } else {
                pruned_sweep_from_certifier(c, trav, diam_lb, ecc_ub_[c]);
                if (mprune > graph.m()) {
                    mprune -= graph.m()+1; ++nsweep; ++nprune_sweep;
                }
            }
            for (int i = 0; i < trav.nvis(); ++i) {
                V u = trav.visit(i);            
                WL ub = optim_ub ? trav.tree_ecc(u) : trav.dist(u) + e_c;
                if (ub <= std::max(ub_required, ecc_lb_[u]))
                    certif_for[c].push_back(u);
            }
        }
        if (mprune > 0) { mprune = 0; ++nsweep; ++nprune_sweep; }
        return gdy_set_cov(certif_for, n);
    }
    
    WL sweepImproveBounds(V u,
                          bool improve_lb = true, bool improve_ub = true,
                          WL desir_lb = max_weight, WL desir_ub = zero_weight,
                          bool optim_lb = false, bool optim_ub = false) {
        WL e_u = sweep(u);
        if (directed) sweep(u, true);
        if(improve_lb) improves_ecc_lb(desir_lb, optim_lb);
        if(improve_ub) improves_ecc_ub(e_u, desir_ub, optim_ub);
        return e_u;
    }

    WL dist_lb(V u, V v) { // works with directed graphs also
        WL d = zero_weight;
        for (int i = dist_lab_lb[u].size() - 1; i >= 0; --i) {
            WL u_pi = dist_lab_lb[u][i], v_pi = dist_lab_lb[v][i];
            if (u_pi < max_weight && v_pi < max_weight) {
                WL diff = u_pi < v_pi ? v_pi - u_pi : u_pi - v_pi;
                if (diff > d) d = diff;
            }
        }
        return d;
    }

    WL dist_ub(V u, V v) {
        WL d = max_weight;
        for (int i = dist_lab_ub[u].size() - 1; i >= 0; --i) {
            WL u_ci = dist_lab_ub[u][i], v_ci = dist_lab_ub[v][i];
            if (u_ci + v_ci < d) d = u_ci + v_ci;
        }
        return d;
    }

private:

    
    /* Sweep from u, i.e. perform a BFS from u if g is unweighted, and
     * a Dijkstra traversal otherwise. 
     */
    WL sweep(V u, bool backward = false) {
        if (u == last_sweep_node
            && (backward == last_sweep_backward || ! directed )) {
            assert(u == trav.first_visited());
            return trav.dist(trav.last_visited());
        }
        last_sweep_node = u;
        last_sweep_backward = backward;
        ++nsweep;
        trav.clear();
        const G &g = backward ? graph_rev : graph;
        int nvis = weighted ? trav.dijkstra(g, u) : trav.bfs(g, u);
        V last = trav.last_visited();
        V ecc = trav.dist(last);
        /*  verb::cerr() << "sweep from "<< u <<" to "<< last
                  <<" at dist "<< ecc
                  <<" in ["<< ecc_lb_[u] <<", "<< ecc_ub_[u] <<"]"
                  <<" diam_todo "<< diam_todo
                  <<std::endl;
        */
        if (ecc < rad_ub && ! (backward && directed))
            { rad_ub = ecc; rad_node = u; update_todo(); }
        if (ecc > diam_lb || diam_node == G::not_vertex) {
            diam_lb = ecc;
            diam_node = backward ? last : u;
            update_todo();
        }
        return ecc;
    }

    void save_trav() {
        int nv = trav.nvis();
        save_visit.resize(nv);
        for (int i = 0; i < nv; ++i) {
            V v = trav.visit(i);
            save_visit[i] = v;
            save_dist[v] = trav.dist(v);
        }
    }

    void save_trav_add() {
        int nv = trav.nvis();
        save_visit.resize(nv);
        for (int i = 0; i < nv; ++i) {
            V v = trav.visit(i);
            save_visit[i] = v;
            if (save_dist[v] < max_weight) save_dist[v] += trav.dist(v);
        }
    }

    static WL f_eccentricity(V v, WL e) { return e; }

    //  asserts f(v,.) non-decr for all v (f(v, e) <= f(v, e') when e <= e')
    V save_argmin_ecc(std::function<WL(V, WL)> f_ecc = f_eccentricity) {
        while (true) {
            V x = G::not_vertex;
            WL f_x = max_weight;
            for (V v : save_visit) {
                if (x == G::not_vertex || f_ecc(v, ecc_lb_[v]) < f_x) {
                    x = v;
                    f_x = f_ecc(v, ecc_lb_[v]);
                }
            }
            assert (x != G::not_vertex);

            WL e_x = sweep(x);

            if (e_x == ecc_lb_[x]) { // ok
                return x;
            } else { // try again
                sweep(trav.last_visited(), true);
                improves_ecc_lb();
            }
        }
        return G::not_vertex;
    }

    V argmin_f(std::function<bool(V)> filter, std::function<WL(V)> f) {
        V x = G::not_vertex;
        WL f_x = max_weight;
        for (V v = 0; v < n; ++v) {
            if (filter(v) && (x == G::not_vertex || f(v) < f_x)) {
                x = v;
                f_x = f(v);
            }
        }
        return x;
    }

    static bool filter_all(V v) { return true; }
    
    V argmin_lb(std::function<bool(V)> filter = filter_all) {
        return argmin_f(filter, [this](V v){ return ecc_lb_[v]; });
    }

    V argmax_ub(std::function<bool(V)> filter = filter_all) {
        return argmin_f(filter, [this](V v){ return max_weight - ecc_ub_[v]; });
    }

    V argmax_sum(std::function<bool(V)> filter = filter_all) {
        return argmin_f(filter, [this](V v){
                assert(sum_[v] >= 0);
                return max_weight - sum_[v];
            });
    }

    V argmax_save(std::function<bool(V)> filter = filter_all) {
        return argmin_f(filter, [this](V v){
                assert(save_dist[v] >= 0);
                return max_weight - save_dist[v];
            });
    }

    V max_degree_node(bool in_largest_scc = true) {
        V md = G::not_vertex;
        for (V u = 0; u < n; ++u) {
            if (scc_nb[u] == scc_largest
                && (md == G::not_vertex || graph.degree(u) > graph.degree(md))){
                md = u;
            }
        }
        return (md == G::not_vertex) ? 0 : md;
    }

    static bool is_in(V u, const std::vector<V> &set) {
        for (V v : set) if (u == v) return true;
        return false;
    }
    
    V furthest_certif_for_root_opt_cov(WL desired_ub) {
        V r = trav.visit(0);
        V f = trav.last_visited();
        if (desired_ub == max_weight) desired_ub = trav.dist(f);
        assert(desired_ub >= trav.dist(f));

        auto ecc_estim = [this](V u) { return ecc_lb_[u]; };

        auto n_cov = [this, desired_ub, ecc_estim](V c) {
            int count = 0;
            const std::vector<WL> & dist_c = trav.dist_from(c);
            for (int i = 0; i < trav.nvis(); ++i) {
                V u = trav.visit(i);
                if (ecc_ub_[u] > desired_ub
                    && dist_c[u] + ecc_estim(c) <= desired_ub) ++count;
            }
            return count;
        };
        
        V c = f;
        while (true) {
            if (trav.dist(c) + ecc_estim(c) <= desired_ub) break;
            c = trav.parent(c);
            if (c == r) break;
        }
        
        int h = 0; for (V u = c; u != r; u = trav.parent(u)) h++;
        int c_cov = n_cov(c);
        for (V u = c; true; u = trav.parent(u)) {
            if (trav.dist(c) + ecc_estim(c) > desired_ub) break;
            int u_cov = n_cov(u);
            if (u_cov > c_cov) { c = u; c_cov = u_cov; }
            if (u == r) break;
        }
        
        return c;
    }

    V furthest_certif_for_root(WL desired_ub = max_weight) {
        V r = trav.visit(0);
        V f = trav.last_visited();
        if (desired_ub == max_weight) desired_ub = trav.dist(f);
        assert(desired_ub >= trav.dist(f));

        auto ecc_estim = [this](V u) { return ecc_lb_[u]; };
        
        V c = f;
        while(true) {
            if (trav.dist(c) + ecc_estim(c) <= desired_ub) break;
            c = trav.parent(c);
            if (c == r) break;
        }
        /* Helps for diameter :
        WL d_c_max = trav.dist(c);
        while (c != r && trav.dist(trav.parent(c)) + ecc_estim(c) <= desired_ub
               && trav.dist(trav.parent(c)) > d_c_max > 2)
            c = trav.parent(c);
        */     
        return c;
    }


    static std::vector<V> empty_vect() { static std::vector<V> v; return v; }

    void clear_ecc_lb(const std::vector<V> &P = empty_vect()) {
        for (V u = 0; u < n; ++u) {
            ecc_lb_[u] = zero_weight;
            lb_node[u] = u;
            next_hop_to_lb[u] = u;
            next_hop_dist[u] = zero_weight;
            diam_c_to_lb[u] = u;
        }
        for (V p : P) {
            last_sweep_node = G::not_vertex;
            trav.clear();
            if (weighted) trav.dijkstra(graph, p); else trav.bfs(graph, p);
            improves_ecc_lb(max_weight, false, false);
        }
    }
    void clear_ecc_ub(const std::vector<V> &C = empty_vect()) {
        for (V u = 0; u < n; ++u) {
            ecc_ub_[u] = max_weight;
        }
        for (V c : C) {
            last_sweep_node = G::not_vertex;
            trav.clear();
            if (weighted) trav.dijkstra(graph, c); else trav.bfs(graph, c);
            improves_ecc_ub(zero_weight, false, false);
        }
    }
    void clear_sum() {
        for (V u = 0; u < n; ++u) sum_[u] = zero_weight;
    }

    bool dist_at_most(V s, V t, WL ub) {
        return false;
    }
    
    bool improves_ecc_lb(WL lb_desired=max_weight,
                         bool opt_lb=false, bool add_to_P = true,
                         bool force_lab=false) {
        bool improves = false;
        V src = trav.first_visited();
        WL ecc_src = trav.dist(trav.last_visited());
        for (int i = 0; i < trav.nvis(); ++i) {
            V u = trav.visit(i);
            if (ecc_lb_[u] < lb_desired) {
                WL lb = opt_lb ? std::max(trav.dist(u), ecc_src - trav.dist(u))
                               : trav.dist(u);
                if (lb > ecc_lb_[u]
                    || lb_node[u] == src) { // diam_c_to_lb can resweep
                    improves = true;
                    break;
                }
            }
        }
        if (improves || force_lab) {
            if ( ! force_lab) last_lb_improve = nsweep;
            for (int i = 0; i < trav.nvis(); ++i) {
                V u = trav.visit(i);
                if (add_to_P || force_lab) {
                    dist_lab_lb[u].push_back(trav.dist(u));
                }
                WL lb = opt_lb ? std::max(trav.dist(u), ecc_src - trav.dist(u))
                               : trav.dist(u);
                if (lb > ecc_lb_[u]
                    || lb_node[u] == src) { // diam_c_to_lb can resweep
                    ecc_lb_[u] = lb;
                    lb_node[u] = src;
                    //next_hop_to_lb[u] = trav.parent(u);
                    V p = trav.parent(u);
                    next_hop_to_lb[u] = p;
                    next_hop_dist[u] = trav.dist(u) - trav.dist(p);
                    while (trav.dist(u) - trav.dist(p)
                           + std::max(ecc_lb_[p], trav.dist(p)) <= diam_lb){
                        // dist(u,p)+ecc(p) <= D if bounds are tight
                        diam_c_to_lb[u] = p;
                        if (trav.parent(p) == p) break;
                        p = trav.parent(p);
                    }
                }
            }
            if (add_to_P) P.push_back(src);
        }
        return improves || ecc_src == 0;
    }
    
    bool improves_ecc_ub(WL ecc_src, WL ub_desired=zero_weight,
                         bool use_tree_ecc=false, bool add_to_C = true,
                         bool use_save_trav=false) {
        V scc_src = scc_nb[trav.first_visited()];
        bool improves = false;
        if (use_tree_ecc) { assert( ! directed); trav.tree_eccentricities(); }
        for (int i = 0; i < trav.nvis(); ++i) {
            V u = trav.visit(i);
            if (ecc_ub_[u] > ub_desired && scc_nb[u] == scc_src) {
                WL d_u = use_save_trav ? save_dist[u] : trav.dist(u);
                WL ub = use_tree_ecc
                    ? trav.tree_ecc(u) : d_u + ecc_src;
                if (ub < ecc_ub_[u]) {
                    improves = true;
                    break;
                }
            }
        }
        if (improves) {
            last_ub_improve = nsweep;
            for (int i = 0; i < trav.nvis(); ++i) {
                V u = trav.visit(i);
                if (scc_nb[u] == scc_src) {
                    WL d_u = use_save_trav ? save_dist[u] : trav.dist(u);
                    if (add_to_C && d_u < dist_C[u])
                        dist_C[u] = d_u; 
                    WL ub = use_tree_ecc ? trav.tree_ecc(u) : d_u + ecc_src;
                    //if (add_to_C) dist_lab_ub[u].push_back(d_u);
                    if (ub < ecc_ub_[u]) ecc_ub_[u] = ub;
                }
            }
            V src = trav.visit(0);
            if (add_to_C) C.push_back(src);
        }
        return improves;
    }

    void update_sums() {
        for (V u = 0; u < n; ++u)  {
            if (trav.visited(u)) sum_[u] += trav.dist(u);
        }
    }

    void update_todo() {
        rad_todo = 0; diam_todo = 0; all_ecc_todo = 0;
        for (V u = 0; u < n; ++u)  {
            if (ecc_lb_[u] < rad_ub && scc_nb[u] == scc_largest)
                ++rad_todo; // restricted to largest scc
            if (ecc_ub_[u] > diam_lb) ++diam_todo;
            if (ecc_lb_[u] < ecc_ub_[u]) ++all_ecc_todo;
        }
        int nswp = nsweep + (mprune > 0 ? 1 : 0);
        if (rad_todo == 0 && rad_nsweep == not_nsweep) {
            while (mprune > graph.m())
                { mprune -= graph.m()+1; ++nsweep; ++nprune_sweep; }
            int incr = (mprune > 0) ? 1 : 0;
            rad_nsweep = nswp + incr;
            verb::cerr() << "radius: " << rad_ub <<" in "<< rad_nsweep
                         <<" sweeps ("<< (nprune_sweep+incr) <<" pruned)\n";
        }
        if (diam_todo == 0 && diam_nsweep == not_nsweep) {
            while (mprune > graph.m())
                { mprune -= graph.m()+1; ++nsweep; ++nprune_sweep; }
            int incr = (mprune > 0) ? 1 : 0;
            diam_nsweep = nswp + incr;
            verb::cerr() << "diameter: " << diam_lb <<" in "<< diam_nsweep
                         <<" sweeps ("<< (nprune_sweep + incr) <<" pruned)\n";
        }
        if (all_ecc_todo == 0 && all_ecc_nsweep == not_nsweep)
                all_ecc_nsweep = nswp;
    }

    std::vector<V> optim_rad_certif(const std::vector<V> & P,
                                    bool all_ecc = false) {
        clear_ecc_lb();
        update_todo();
        std::vector<V> Popt;
        last_sweep_node = G::not_vertex;
        for (int i = P.size() - 1; i >= 0; --i) {
            V p = P[i];
            trav.clear();
            if (weighted) trav.dijkstra(graph, p); else trav.bfs(graph, p);
            if (improves_ecc_lb(all_ecc ? max_weight : rad_ub, false, false)) {
                Popt.push_back(p);
                update_todo();
                if (all_ecc ? all_ecc_todo == 0 : rad_todo == 0) break;
            }
        }
        assert(all_ecc ? all_ecc_todo == 0 : rad_todo == 0);
        return Popt;
    }
    
    std::vector<V> optim_diam_certif(const std::vector<V> &C,
                                     bool all_ecc = false) {
        clear_ecc_ub();
        update_todo();
        std::vector<V> Copt;
        last_sweep_node = G::not_vertex;
        for (int i = C.size() - 1; i >= 0; --i) {
            V c = C[i];
            trav.clear();
            if (weighted) trav.dijkstra(graph, c); else trav.bfs(graph, c);
            if (improves_ecc_ub(all_ecc ? zero_weight : diam_lb, false, false)){
                Copt.push_back(c);
                update_todo();
                if (all_ecc ? all_ecc_todo == 0 : diam_todo == 0) break;
            }
        }
        assert(all_ecc ? all_ecc_todo == 0 : diam_todo == 0);
        return Copt;
    }

    /* Compute an even smaller certificate for upper_bounds.
     * Assumes lower bounds are tight to compute efficiently
     * what nodes a given c can certify (in graphs with small all ecc. certif)
     * and then apply a greedy-set-cover like heuristic.
     * @param all_ecc specifies if certificates for all eccentricities are
     *   to be found or only necessary certificates for showing diam <= diam_lb.
     * @param prune_budget allows to avoid too much pruned sweeps.
     *   Too much pruned sweeps occur when lower bounds are not good enough
     *   (especially when certifying diameter only), it is then better 
     *    to call all_basic() to improve them (see diameter()).
     */
    void optim_ub_when_lb_tight(bool all_ecc,
                                int prune_budget = INT_MAX) {
        std::vector<V> todo;
        std::vector<int> certif_deg(n,0); // how much a node can certify
        std::vector<std::vector<V>> certif(n); // certifiers of a node

        last_sweep_node = G::not_vertex;
        bool do_prune = true;
        int iter = 0, m_vis = 0;
        while ((all_ecc ? all_ecc_todo : diam_todo) > 0
               && prune_budget >= 0) {
            ++iter;
            
            todo.clear();
            for (V u = 0; u < n; ++u)  {
                if (ecc_ub_[u] > (all_ecc ? ecc_lb_[u] : diam_lb)) {
                    todo.push_back(u);
                }
            }

            // Pruned sweeps
            int pruned_sweeps = 0;
            if (do_prune) {
                do_prune = false;
                for (V u = 0; u < n; ++u) {
                    certif[u].clear();
                    certif_deg[u] = 0;
                }
                for (V p : todo) {
                    auto filter = [this, p, &m_vis, &pruned_sweeps,
                                   &prune_budget, all_ecc]
                        (V u, WL d, V par, WL dpar) {
                        ++m_vis;
                        if (m_vis >= graph.m()) {
                            m_vis = 0;
                            ++pruned_sweeps;
                            --prune_budget;
                        }
                        return prune_budget >= 0 &&
                        d + ecc_lb_[u] <= (all_ecc ? ecc_lb_[p] : diam_lb);
                    };
                    trav.clear();
                    int nv = weighted ? trav.dijkstra(graph, p, filter)
                        : trav.bfs(graph, p, filter);
                    // visited nodes are certificates for p:
                    for (int i = 0; i < nv; ++i) {
                        V c = trav.visit(i);
                        certif[p].push_back(c);
                        certif_deg[c] += 1;
                    }
                }
            }
            nsweep += pruned_sweeps;

            // Select c that is a certificate for a maximum number of todo nodes
            int c = G::not_vertex;
            for (V u = 0; u < n; ++u)
                if (c == G::not_vertex || certif_deg[u] > certif_deg[c])
                    c = u;
            assert(c != G::not_vertex);
            
            WL e_c = sweep(c);

            int n_certified = 0;
            V p_false_lb = G::not_vertex;
            for (V u = 0; u < n; ++u)
                if (ecc_ub_[u] > (all_ecc ? ecc_lb_[u] : diam_lb)
                    && trav.dist(u) + e_c <= (all_ecc ? ecc_lb_[u] : diam_lb)) {
                    ++n_certified;
                    for (V c : certif[u]) --(certif_deg[c]);
                    certif[u].clear();
                }
            improves_ecc_ub(all_ecc ? zero_weight : diam_lb);
            if (ecc_lb_[c] < e_c) { // possible failure of certification
                V p = trav.last_visited();
                sweep(p);
                improves_ecc_lb();
                do_prune = true;
                certif_deg[c] = 0;
            } else {
                assert(certif_deg[c] <= 0);
            }

            update_todo();

            if (verb::progress() || iter <= 3)
                verb::cerr("optim_lb", 1) <<"nsweep: " << nsweep
                      <<" pruned: " << pruned_sweeps
                      <<"   c=" << c
                      <<" todo: "<< (all_ecc ? all_ecc_todo : diam_todo)
                      <<" certified: "<< n_certified
                      <<" C:" << C.size()
                      << "\n";
        }
        if (m_vis > 0) ++nsweep;
    }


public:
    
    /* Compute a small certificate for upper_bounds similarly to
     * optim_ub_when_lb_tight().
     * Assumes that lower bounds are tight enough to proceed in one shot.
     */
    std::vector<V> optim_ub_certif_one_shot(bool all_ecc, bool use_center) {

        std::vector<V> todo; // nodes to certify
        std::vector<std::vector<V>> certif_for(n); // nodes certified by a node
        
        // Center is usually a certificate for many nodes:
        if (use_center) {
            assert(rad_todo == 0);
            sweep(rad_node, true);
            for (V u = 0; u < n; ++u)
                if (trav.dist(u) + rad_ub > (all_ecc ? ecc_lb_[u] : diam_lb))
                    todo.push_back(u);
                else
                    certif_for[rad_node].push_back(u);
        } else {
            for (V u = 0; u < n; ++u) todo.push_back(u);
        }
        
        // Pruned sweeps
        int m_vis = 0, pruned_sweeps = 0;
        last_sweep_node = G::not_vertex;
        for (V p : todo) {
            auto filter = [this, p, &m_vis, &pruned_sweeps, all_ecc]
                (V u, WL d, V par, WL dpar) {
                ++m_vis;
                if (m_vis >= graph.m()) { m_vis = 0; ++pruned_sweeps; }
                return d + ecc_lb_[u] <= (all_ecc ? ecc_lb_[p] : diam_lb);
            };
            trav.clear();
            int nv = weighted ? trav.dijkstra(graph, p, filter)
                : trav.bfs(graph, p, filter);
            // visited nodes are certificates for p:
            for (int i = 0; i < nv; ++i) {
                V c = trav.visit(i);
                certif_for[c].push_back(p);
            }
        }
        if (m_vis > 0) ++pruned_sweeps; 
        nsweep += pruned_sweeps;
        
        std::vector<V> cov = gdy_set_cov(certif_for, n);
        return cov;
    }

    std::vector<V> diam_certif_when_lb_tight(bool use_center=true) {
        std::vector<V> todo; // nodes to certify
        std::vector<std::vector<V>> certif_for(n); // nodes certified by a node
        
        // Center is often a certificate for many nodes:
        if (use_center) {
            sweep(rad_node, true);
            for (V u = 0; u < n; ++u)
                if (trav.dist(u) + rad_ub > diam_lb) {
                    todo.push_back(u);
                } else {
                    //certif_for[rad_node].push_back(u);
                }
        } else {
            for (V u = 0; u < n; ++u) todo.push_back(u);
        }

        verb::cerr() << "diam_certif_when_lb_tight todo: "
                     << todo.size() << "\n";

        for (V u : todo) {
            pruned_sweep_to_certifiers(u, trav, diam_lb);
            for (int i = trav.nvis() - 1; i >= 0; --i) {
                V c = trav.visit(i);
                certif_for[c].push_back(u);
            }
        }
        if (mprune > 0) { ++nsweep; mprune = 0; ++nprune_sweep; }

        //int lb = gdy_packing(certif_for, n).size();
        //verb::cerr() << "lower bound : " << lb << std::endl;
        
        std::vector<V> C = gdy_set_cov(certif_for, n, todo);
        if (use_center) C.push_back(rad_node);
        return C;
        //return gdy_set_cov(certif_for, n);
    }
    
    std::vector<V> diam_certif_when_lb_tight_threaded(int n_thd,
                                                      bool use_center=true) {
        std::vector<V> todo = { -1 }; // nodes to certify
        std::vector<std::vector<V>> certif_for(n); // nodes certified by a node
        std::vector<std::vector<V>> certif(n); // certifiers of a node
        std::mutex this_mutex;
        
        // Center is often a certificate for many nodes:
        if (use_center) {
            todo.clear();
            sweep(rad_node, true);
            for (V u = 0; u < n; ++u)
                if (trav.dist(u) + rad_ub > diam_lb) {
                    todo.push_back(u);
                }
        }

        verb::cerr() << "diam_certif_when_lb_tight_threaded(): "
                     << todo.size() << "\n"; std::cerr.flush();

        auto go = [&todo, &certif, this, &this_mutex, n_thd](int i_thd) {
            int _nsweep = 0, _mprune = 0, _nprune_sweep = 0;
            traversal<G> trav(n);
            
            for (int i = todo.size() - 1; i >= 0; i--) {
                if (i % n_thd == i_thd) {
                    V u = todo[i];
                    pruned_sweep_to_certifiers_threaded
                        (u, trav, _nsweep, _mprune, _nprune_sweep, diam_lb);
                    for (int i = trav.nvis() - 1; i >= 0; --i) {
                        V c = trav.visit(i);
                        certif[u].push_back(c);
                    }
                }
                if (i_thd == 0 && verb::progress())
                    verb::cerr("diam_certif_when_lb_tight_threaded()", 1)
                        <<"todo: " << i <<" / "
                        << todo.size() <<"\n";
            }

            this_mutex.lock();
            nsweep += _nsweep;
            mprune += _mprune;
            nprune_sweep += _nprune_sweep;
            this_mutex.unlock();
        };

        std::vector<std::thread> threads(n_thd);
        for (int i = 0; i < n_thd; ++i) {
            threads[i] = std::thread(go, i);
        }
        for (int i = 0; i < n_thd; ++i) threads[i].join();

        while (mprune > graph.m()) {
            ++nsweep; mprune -= graph.m()+1; ++nprune_sweep;
        }
        if (mprune > 0) { ++nsweep; mprune = 0; ++nprune_sweep; }
        verb::cerr("diam_certif_when_lb_tight_threaded()", 1)
            << "nprune_sweep: " << nprune_sweep <<"\n"; std::cerr.flush();
        
        for (V u = 0; u < n; ++u) {
            for (V c : certif[u])
                certif_for[c].push_back(u);
            certif[u] = {}; // release mem
        }
        
        verb::cerr("diam_certif_when_lb_tight_threaded()", 1)
            << "gdy cover...\n"; std::cerr.flush();

        //int lb = gdy_packing(certif_for, n).size();
        //verb::cerr() << "lower bound : " << lb << std::endl;

        std::vector<V> C = gdy_set_cov(certif_for, n, todo);
        if (use_center) C.push_back(rad_node);
        return C;
    }
    
public:

    static std::vector<int>
    gdy_set_cov(const std::vector<std::vector<int>> &set,
                const int n, // sets are subsets of 0..n-1
                std::vector<int> to_be_covered = {-1}) {
        int m = set.size(); // number of sets
        std::vector<std::vector<int>> elt_of(n);

        // what can be covered: the union of all sets in set.
        std::vector<bool> can_be_covered(n);
        for (int s = 0; s < m; ++s)
            for (int u : set[s])
                can_be_covered[u] = true;
        // what nodes should we cover: 
        std::vector<bool> covered(n, true);
        int ncov = 0, n_to_cov = 0;
        if (to_be_covered.size() == 1 && to_be_covered[0] == -1) {
            // default : cover all
            to_be_covered.clear();
            to_be_covered.reserve(n);
            for (int u = 0; u < n; ++u) to_be_covered.push_back(u);
        }
        for (int u : to_be_covered) {
            if ( ! can_be_covered[u]) {
                std::cerr << u <<" cannot be covered!\n";
                throw (std::string("a node cannot be covered: ")
                       + std::to_string(u));
            }
            if (covered[u]) {
                covered[u] = false;
                ++n_to_cov;
            }
        }
            
        std::vector<int> degree(n, 0);
        std::vector<int> cov, // solution (sets that cover nodes)
            onecov; // nodes covered by only one set
        std::vector<bool> in_cov(n, false); // selected in solution
        for (int s = 0; s < m; ++s) {
            int deg = 0;
            for (int u : set[s]) {
                if ( ! covered[u]) { // consider only elts we have to cover
                    elt_of[u].push_back(s);
                    ++deg;
                }
            }
            degree[s] = deg;
        }
        for (int u = 0; u < n; ++u) {
            if (elt_of[u].size() == 1) { // only one set covers it, take it
                int s = elt_of[u][0];
                if ( ! in_cov[s]) {
                    onecov.push_back(s); // add it to cov later on
                    in_cov[s] = true; 
                }
            }
        }
        verb::cerr("gdy_cov") << m <<" sets "<< onecov.size() << " mandatory\n";

        auto cover = [&covered, &ncov, &elt_of, &degree](int u, int s) {
            if ( ! covered[u]) { // cover u with s
                covered[u] = true;
                ++ncov;
                for (int s : elt_of[u]) degree[s] -= 1;
            }
        };
        // mandatory ones : 
        for (int s : onecov) {
            for (int u : set[s]) cover(u, s);
        }
        std::vector<int> onecov_inter;
        while (ncov < n_to_cov) {
            int md = 0;
            // max degree set (covering most uncovered nodes)
            for (int s = 0; s < m; ++s) if (degree[s] > degree[md]) md = s;
            // add mandatory sets with degree at least that of md :
            onecov_inter.clear();
            for (int t : onecov) {
                if (degree[t] + 1 >= degree[md]) cov.push_back(t);
                else onecov_inter.push_back(t);
            }
            std::swap(onecov, onecov_inter);

            if (verb::progress())
                verb::cerr("gdy_cov", 1) << "deg " << degree[md] <<"\n";
            // add max degree set :
            cov.push_back(md);
            for (int u : set[md]) cover(u, md);
        }
        // add remaining mandatory ones :
        for (int s : onecov) {
            cov.push_back(s);
        }
        return cov;
    }

    /** Returns a maximal packing giving a lower bound on the size of
     * any set covering (a packing is a set of elements such that no set
     * contains two of them, any covering must thus use at least a distinct
     * set for covering each element of the packing).
     */
    static std::vector<int>
    gdy_packing_disj(const std::vector<std::vector<int>> &orig_set,
                const int n /* sets are subsets of 0..n-1 */,
                const std::vector<int> to_be_covered = {-1}) {
        
        // to_be_covered as a bool array:
        std::vector<bool> elt_remains(n, true);
        if (to_be_covered.size() != 1 || to_be_covered[0] != -1) {// not default
            for (int i = 0; i < n; ++i) elt_remains[i] = false;
            for (int i : to_be_covered) elt_remains[i] = true;
        }

        const int m = orig_set.size();
        std::vector<std::vector<int>> elt(n);
        for (int s = 0; s < m; ++s) {
            for (int i : orig_set[s])
                if (elt_remains[i]) elt[i].push_back(s);
        }

        return gdy_disjoint(elt, m);
    }
    
    /** Returns a maximal packing giving a lower bound on the size of
     * any set covering (a packing is a set of elements such that no set
     * contains two of them, any covering must thus use at least a distinct
     * set for covering each element of the packing).
     */
    static std::vector<int>
    gdy_packing(const std::vector<std::vector<int>> &orig_set,
                const int n /* sets are subsets of 0..n-1 */,
                const std::vector<int> to_be_covered = {-1},
                int n_thd = std::thread::hardware_concurrency()) {
        // Technically, we pack the sets in elt (see constr. bellow),
        // which gives a lower bound on a set cover for sets in set.

        std::vector<int> pck; // elts belonging to disjoint sets

        // to_be_covered as a bool array:
        std::vector<bool> elt_remains(n, true);
        if (to_be_covered.size() != 1 || to_be_covered[0] != -1) {// not default
            for (int i = 0; i < n; ++i) elt_remains[i] = false;
            for (int i : to_be_covered) elt_remains[i] = true;
        }
        
        // sets that contain an element:
        const int m = orig_set.size();
        std::vector<std::vector<int>> elt(n);
        for (int s = 0; s < m; ++s) {
            for (int i : orig_set[s])
                if (elt_remains[i]) elt[i].push_back(s);
        }

        // check that elements to be covered by sets are indeed covered
        for (int i = 0; i < n; ++i)
            assert(( ! elt_remains[i]) || elt[i].size() > 0); // i is covered
        
        std::vector<int> pck_disj = gdy_disjoint(elt, m);

        // Some sets must be in any covering, treat them specifically :
        std::vector<bool> set_remains(m, true);
        for (int i = 0; i < n; ++i) {
            if (elt_remains[i] && elt[i].size() == 1) {
                /* elt[i] must be in any covering (to cover i),
                 * Any maximal packing must hit this set (if the set is not
                 * hit, then we can safely add i in the packing).  */ 
                pck.push_back(i);
                int s = elt[i][0];
                for (int j : orig_set[s]) elt_remains[j] = false;
                set_remains[s] = false;
            }
        }

        // Clean up sets
        std::vector<std::vector<int>> set(m);
        for (int i = 0; i < n; ++i) {
            std::vector<int> e_inter;
            std::swap(elt[i], e_inter);
            if (elt_remains[i]) {
                for (int s : e_inter)
                    if (set_remains[s]) {
                        elt[i].push_back(s);
                        set[s].push_back(i);
                    }
            }
        }

        // neighbors are elements belonging to a common set
        // we avoid to construct this graph because of memory explosion
        // (code optimized for large sets)
        std::vector<std::vector<int>> neighb(n);
        std::vector<int> deg(n, 0);
        std::vector<bool> in_union(m, false);
        int64_t sum_size = 0, orig_sum_size = 0, work_estim = 0, mil = 1000000;
        for (int s = 0; s < m; ++s) {
            work_estim += set[s].size() * set[s].size();
            sum_size += set[s].size();
            orig_sum_size += orig_set[s].size();
        }
        int64_t affordable_degree = 3 * orig_sum_size / n;
        verb::cerr() << " work_estim: " << work_estim
                  << " set sum: " << sum_size << " orig sum: " << orig_sum_size
                  <<"  afford_deg: " << affordable_degree << std::endl;
        if (work_estim > std::max(1000 * orig_sum_size, 10000*mil)) {
            verb::cerr() << "use faster heuristic!\n";
            std::vector<int> to_cov;
            for (int i = 0; i < n; ++i)
                if (elt_remains[i]) to_cov.push_back(i);
            std::vector<int> pck2 = gdy_packing_fast(set, elt, n, to_cov);
            pck.insert(pck.end(), pck2.begin(), pck2.end());

            if (pck_disj.size() >= pck.size()) {
                verb::cerr("packing", 2) << "use pck_disj\n";
                std::swap(pck, pck_disj);
            }
            
            // elts of orig_sets
            for (int i = 0; i < n; ++i) elt[i].clear();
            for (int s = 0; s < m; ++s) {
                for (int i : orig_set[s]) elt[i].push_back(s);
            }
            check_packing(orig_set, elt, n, pck, to_be_covered);
            return pck;
        }

        // Compute degrees of neighb graph (i and j are neighbors is some set
        //   contains them both).
        auto go = [n, m, &elt_remains, &elt, &set_remains, &set, &deg, &neighb,
                   affordable_degree] (int i_thd, int n_thd) {
            std::vector<bool> in_union(m, false);
            for (int i = 0; i < n; ++i) {
                if (i_thd == 0 && verb::progress())
                    verb::cerr("packing neighb", 1)
                        <<"todo: "<< (n - i) <<" / "<< n <<"\n";
                if (i % n_thd == i_thd && elt_remains[i]) {
                    int count_dupl = 0;
                    for (int s : elt[i]) {
                        if (set_remains[s])
                            for (int j : set[s])
                                if (elt_remains[j]) {
                                    ++count_dupl;
                                    if ( ! in_union[j]) {
                                        ++(deg[i]);
                                        if (deg[i] <= affordable_degree)
                                            neighb[i].push_back(j);
                                        in_union[j] = true;
                                    }
                                }
                    }
                    if (deg[i] > affordable_degree) {
                        neighb[i].clear(); // recompute it if needed
                        neighb[i].shrink_to_fit();
                    }
                    // reset union
                    if (count_dupl >= n) {
                        for (int j = 0; j < n; ++j) {
                            if (in_union[j]) {
                                in_union[j] = false;
                            }
                        }
                    } else {
                        for (int s : elt[i]) {
                            if (set_remains[s])
                                for (int j : set[s])
                                    if (elt_remains[j] && in_union[j]) {
                                        in_union[j] = false;
                                    }
                        }
                    }
                }
            }
        };
        if (work_estim < 1000 * mil) {
            go(0, 1);
        } else {
            std::vector<std::thread> threads(n_thd);
            for (int i = 0; i < n_thd; ++i) {
                threads[i] = std::thread(go, i, n_thd);
            }
            for (int i = 0; i < n_thd; ++i) threads[i].join();
        }
        

        /* go sequential : 
        for (int i = 0; i < n; ++i) {
            if (elt_remains[i]) {
                int count_dupl = 0;
                for (int s : elt[i]) {
                    if (set_remains[s])
                        for (int j : set[s])
                            if (elt_remains[j]) {
                                ++count_dupl;
                                if ( ! in_union[j]) {
                                    ++(deg[i]);
                                    if (deg[i] <= affordable_degree)
                                        neighb[i].push_back(j);
                                    in_union[j] = true;
                                }
                    }
                }
                if (deg[i] > affordable_degree) {
                    neighb[i].clear(); // recompute it if needed
                    neighb[i].shrink_to_fit();
                }
                // reset in_union
                if (count_dupl >= n) {
                    for (int j = 0; j < n; ++j) {
                        if (in_union[j]) {
                            in_union[j] = false;
                        }
                    }
                } else {
                    for (int s : elt[i]) {
                        if (set_remains[s])
                            for (int j : set[s])
                                if (elt_remains[j] && in_union[j]) {
                                    in_union[j] = false;
                                }
                    }
                }
            }
        }
        */

        std::vector<int> sets_to_remove;
        
        // Greedily add nodes minimal degree to the packing
        while (true) {
            // take min deg elt
            int mde = -1, n_remain = 0;
            for (int i = 0; i < n; ++i) {
                if (elt_remains[i]) {
                    ++n_remain;
                    if(mde == -1 || deg[i] < deg[mde])
                        mde = i;
                }
            }
            if (mde == -1) break;
            assert(deg[mde] > 0);
            
            pck.push_back(mde);

            if (verb::progress())
                verb::cerr("packing", 1)
                    << "mde: " << mde <<" "<< deg[mde]
                    <<" todo: "<< n_remain <<" / "<< n <<std::endl;

            // remove neighbors (that includes mde) :
            if (deg[mde] >= n_remain) break; // will remove all remaining
            sets_to_remove.clear();
            for (int s : elt[mde])
                if (set_remains[s]) {
                    set_remains[s] = false;
                    sets_to_remove.push_back(s);
                }
            for (int s : sets_to_remove)
                for (int j : set[s])
                    if (elt_remains[j]) {
                        elt_remains[j] = false;
                        if (neighb[j].size() > 0) { //deg_ok, we have neighb
                            for (int k : neighb[j])
                                if (elt_remains[k])
                                    --(deg[k]);
                        } else { // recompute neighb[j]
                            int count_dupl = 0;
                            for (int t : elt[j])
                                if (set_remains[t])
                                    for (int k : set[t]) {
                                        ++count_dupl;
                                        if (elt_remains[k]
                                            && ! in_union[k]) {
                                            in_union[k] = true;
                                            --(deg[k]);
                                        }
                                    }
                            // reset in_union
                            if (count_dupl >= n) {
                                for (int k = 0; k < n; ++k)
                                    in_union[k] = false;
                            } else {
                                for (int t : elt[j])
                                    if (set_remains[t])
                                        for (int k : set[t])
                                            if (elt_remains[k]) {
                                                in_union[k] = false;
                                            }
                            }
                        }
                    }
            assert( ! elt_remains[mde]); 
        }

        verb::cerr("packing", 2) << "full heur. vs pck_disj: "
                                 << pck.size() <<", "<< pck_disj.size() <<"\n";
        if (pck_disj.size() >= pck.size()) {
            verb::cerr("packing", 2) << "use pck_disj\n";
            std::swap(pck, pck_disj);
        }

        // elts of orig_sets
        for (int i = 0; i < n; ++i) elt[i].clear();
        for (int s = 0; s < m; ++s) {
            for (int i : orig_set[s]) elt[i].push_back(s);
        }
        check_packing(orig_set, elt, n, pck, to_be_covered);
        return pck;
    }

    // Check that each set contains at most one element in pck
    static
    void check_packing(const std::vector<std::vector<int>> &set,
                       const std::vector<std::vector<int>> &elt,
                       const int n /* sets are subsets of 0..n-1 */,
                       const std::vector<int> pck,
                       const std::vector<int> to_be_covered = {-1}) {
        const int m = set.size();
        assert(n == elt.size());

        // check that a set contains at most one element of the packing pck
        std::vector<bool> in_pck(n, false);
        for (int i : pck) in_pck[i] = true;
        std::vector<bool> is_hit(n, false), duplicate(n, false);
        for (int s = 0; s < m; ++s) {
            int hits = 0;
            for (int i : set[s])
                if (in_pck[i] && ! duplicate[i]) {
                    duplicate[i] = true;
                    ++hits; assert(hits <= 1);
                }
            for (int i : set[s]) duplicate[i] = false;
            if (hits > 0) is_hit[s] = true;
        }
        
        // Check no element can be added to the packing (maximality)
        auto check = [&elt, &is_hit](int i) {
            bool cannot_add = false;
            for (int s : elt[i]) if (is_hit[s]) cannot_add = true;
            assert(cannot_add);
        };
        if (to_be_covered.size() != 1 || to_be_covered[0] != -1) {
            for (int i : to_be_covered) check(i);
        } else { // default : all elements
            for (int i = 0; i < n; ++i) check(i);
        }
    }

    // return a maximal set of elts in to_be_covered
    // s.t. no set contains two or more elts
    static std::vector<int>
    gdy_packing_fast(const std::vector<std::vector<int>> &set,
                     const std::vector<std::vector<int>> &elt,
                     const int n /* sets are subsets of 0..n-1 */,
                     std::vector<int> to_be_covered = {-1}) {
        // Technically, we pack the sets in elt (see constr. bellow),
        // which gives a lower bound on a set cover for sets in set.

        std::vector<int> pck; // elts belonging to disjoint sets

        // sets that contain an element:
        const int m = set.size();
        assert(n == elt.size());

        // check that elements to be covered by sets are indeed covered
        std::vector<bool> elt_remains(n, true);
        if (to_be_covered.size() != 1 || to_be_covered[0] != -1) {// not default
            for (int i = 0; i < n; ++i) elt_remains[i] = false;
            for (int i : to_be_covered) elt_remains[i] = true;
        }
        for (int i = 0; i < n; ++i)
            assert(( ! elt_remains[i]) || elt[i].size() > 0);

        // Greedily add nodes minimal degree to the packing
        while (true) {
            // take min deg elt
            int mde = -1;
            for (int i = 0; i < n; ++i) {
                if (elt_remains[i] && (mde == -1
                                       || elt[i].size() < elt[mde].size()))
                    mde = i;
            }
            if (mde == -1) break;
            
            pck.push_back(mde);

            if (verb::progress())
                verb::cerr("fst_packing", 1) << "f mde: " << mde
                                             <<" "<< elt[mde].size() <<"\n";

            for (int s : elt[mde])
                for (int j : set[s])
                    elt_remains[j] = false;
        }

        verb::cerr() << "packing of "<< pck.size() <<"\n";
        
        check_packing(set, elt, n, pck, to_be_covered);
        return pck;
    }

    static
    std::vector<int> disjoint_basic(const std::vector<std::vector<int>> &set,
                              const int n /* sets are subsets of 0..n-1 */) {
        std::vector<int> disj;

        int m = set.size();
        std::vector<bool> seen(n, false), taken(m, false);

        int nset = m;
        while (nset > 0) {
            int sm = -1;
            for (int s = 0; s < m; ++s) {
                if ( ! taken[s]) {
                    if (sm == -1 || set[s].size() < set[sm].size()) sm = s;
                }
            }
            assert(sm != -1);
            taken[sm] = true;
            --nset;
            bool is_disj = true;
            for (int u : set[sm]) {
                if (seen[u]) { is_disj = false; break; }
            }
            if (is_disj) {
                disj.push_back(sm);
                for (int u : set[sm]) {
                    seen[u] = true;
                }
            }
        }
        
        check_disjoint(set, n, disj);
        return disj;
    }

    static
    std::vector<int> gdy_disjoint(const std::vector<std::vector<int>> &set,
                                  const int n /* sets are subsets of 0..n-1 */){
        std::vector<int> disj;

        int m = set.size(), nset = 0;
        std::vector<bool> seen(n, false), taken(m, false);

        std::vector<int> deg(n,0), set_size(m,0);
        for (int s = 0; s < m; ++s) {
            for (int u : set[s])
                ++(deg[u]);
        }
        for (int s = 0; s < m; ++s) {
            int sz = 0;
            for (int u : set[s]) {
                if (deg[u] > 1) ++sz;
            }
            set_size[s] = sz;
            if (set[s].size() > 0) ++nset;
        }
        
        while (nset > 0) {
            int sm = -1;
            for (int s = 0; s < m; ++s) {
                if ( set[s].size() > 0 && ! taken[s]) {
                    if (sm == -1 || set_size[s] < set_size[sm]) sm = s;
                }
            }
            assert(sm != -1);
            taken[sm] = true;
            --nset;
            bool is_disj = true;
            for (int u : set[sm]) {
                if (seen[u]) { is_disj = false; break; }
            }
            if (is_disj) {
                disj.push_back(sm);
                for (int u : set[sm]) {
                    seen[u] = true;
                }
            }
        }
        
        check_disjoint(set, n, disj);
        return disj;
    }

    static
    void check_disjoint(const std::vector<std::vector<int>> &set,
                        const int n /* sets are subsets of 0..n-1 */,
                        const std::vector<int> &disj) {
        std::vector<bool> seen(n, false);
        for (int s : disj) {
            for (int u : set[s]) {
                assert( ! seen[u]);
                seen[u] = false;
            }
        }
    }

    static size_t argmin(const std::vector<WL> &v) {
        size_t n = v.size();
        size_t a = n;
        for (size_t i = 0; i < n; ++i) {
            if (a == n || v[i] < v[a])
                a = i;
        }
        assert(a != n);
        return a;
    }

    static size_t argmax(const std::vector<WL> &v) {
        size_t n = v.size();
        size_t a = n;
        for (size_t i = 0; i < n; ++i) {
            if (a == n || v[i] > v[a])
                a = i;
        }
        assert(a != n);
        return a;
    }

}; // eccentricity



#endif // ECCENTRICITY_HH

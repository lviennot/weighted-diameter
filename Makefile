
# --------- cpp tool ------

main: mdiam

MHDRS:=$(addprefix src/, edge.hh mgraph.hh dyn_graph.hh traversal.hh\
  eccentricity.hh graphgen.hh treedec.hh pruned_landmark_labeling.hh skeleton.hh verbose.hh)

mgraph: src/mgraph_test.cc $(MHDRS)
	@grep -e 'BUG\|TODO\|RM' $^ || true
	g++ -std=c++11 -O3 -g -rdynamic -pthread -o $@ $<

mdiam: src/mgraph_diam.cc $(MHDRS)
	@grep -e 'BUG\|TODO\|RM' $^ || true
	g++ -std=c++11 -O3 -g -rdynamic -pthread -o $@ $<

munweight: src/mgraph_unweight.cc $(MHDRS)
	@grep -e 'BUG\|TODO\|RM' $^ || true
	g++ -std=c++11 -O3 -g -rdynamic -pthread -o $@ $<





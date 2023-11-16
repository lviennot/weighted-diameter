
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




# ------ graph drawing -----

GRAPHVIZ:=neato -Ksfdp 
GRAPHVIZ2:=neato -Ksfdp -Goverlap=scale -Gsplines=curved -Nlabel="" -Earrowhead=none -Nshape=circle -Nstyle=filled -Nwidth=.1 -Nfixed-size=true -Nfontsize=15 -Ncolor="\#00000060" -Ecolor="\#00000020"

%.dot: %.txt
	(echo "graph g {"; \
	cat $< | awk '{if($$1 < $$2) print("  ", $$1, " -- ", $$2, " ;");}'; \
	echo "}") > $@

%.pdf: %.dot
	$(GRAPHVIZ) -o $@ -Tpdf $<

%.force:
	rm -f $*; make $*



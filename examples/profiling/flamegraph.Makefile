# This Makefile profiles stf using perf to produce flame graphs
.PHONEY: info record_prereq fold flamegraph run display

#script that produces lots of qsub scripts to run velociraptor on simulation output
info: 
	echo "This script profiles a VR run."

# PERF
_PERF_OPT_RECORD=--call-graph lbr -g -s -o
_PERF_PATH_ROOT=`pwd`/fg
_PERF_FILE_RAW=$(_PERF_PATH_ROOT)/perf.raw.`date +%Y%m%d`.`hostname`.data

# FLAMEGRAPH
_FLAMEGRAPH_GIT_URL=https://github.com/brendangregg/FlameGraph.git
_FLAMEGRAPH_PATH=$(_PERF_PATH_ROOT)/FlameGraph
_FLAMEGRAPH_FILE_FOLDED=$(_PERF_PATH_FOLDED)/perf.processed.`date +%Y%m%d`.`hostname`.folded
_FLAMEGRAPH_SCRIPT_STACKCOLLAPSE=$(_FLAMEGRAPH_PATH)/stackcollapse-perf.pl
_FLAMEGRAPH_SCRIPT_MAIN=$(_FLAMEGRAPH_PATH)/flamegraph.pl
_FLAMEGRAPH_FILE_SVG=$(_PERF_PATH_ROOT)/perf.html.`date +%Y%m%d`.`hostname`.svg

# STF
# Dynamic stf input parameters through environment variables to be set in an input payload "./payloadname.env "
_STF_PARAMETERS="DEFAULT"
# include ./payloadname.env 

record_prereq:
	@mkdir -p $(_PERF_PATH_ROOT)
	@cd $(_PERF_PATH_ROOT) && git clone $(_FLAMEGRAPH_GIT_URL) && cd..

record: record_prereq
	perf record $(_PERF_OPT_RECORD) -o $(_PERF_FILE_RAW) ./stf $(_STF_PARAMETERS)
	@echo "To manually test recorded data, execute: "
	@echo "perf report -i $(_PERF_FILE_RAW)"

fold:
	perf script -i $(_PERF_FILE_RAW) | $(_FLAMEGRAPH_SCRIPT_STACKCOLLAPSE) > $(_FLAMEGRAPH_FILE_FOLDED)

flamegraph:
	cat $(_FLAMEGRAPH_FILE_FOLDED) | $(_FLAMEGRAPH_SCRIPT_MAIN) > $(_FLAMEGRAPH_FILE_SVG)

run: record fold flamegraph

# Display output for proofreading
# It can be appended to job scripts or test run in CLI
display:
	@echo "========================================================================"
	@echo "Record execution with perf"
	@echo "perf record $(_PERF_OPT_RECORD) $(_PERF_FILE_RAW) ./stf $(_STF_PARAMETERS)"
	@echo ""
	@echo "========================================================================"
	@echo "Convert perf raw data to flamegraph folded data"
	@echo "perf script -i $(_PERF_FILE_RAW) | $(_FLAMEGRAPH_SCRIPT_STACKCOLLAPSE) > $(_FLAMEGRAPH_FILE_FOLDED)"
	@echo ""
	@echo "========================================================================"
	@echo "Convert flamegraph folded data to html format"
	@echo "cat $(_FLAMEGRAPH_FILE_FOLDED) | $(_FLAMEGRAPH_SCRIPT_MAIN) > $(_FLAMEGRAPH_FILE_SVG)"
	@echo ""
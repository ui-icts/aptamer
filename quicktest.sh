set -e

readonly graph_file=data/Trunc_Test/Trunc_test.fa.struct.fa.xgmml

echo "Running predict structures"
./predict-structures samples/Trunc_Test.fa --calculate_stats > /dev/null

echo "Comparing create graph output (no seed)"
rm -f $graph_file
./create-graph data/Trunc_Test/Trunc_Test.fa.struct.fa > /dev/null
icdiff test/create_graph_golden_master_no_seed.xgmml $graph_file

echo "Comparing create graph output (seed)"
rm -f $graph_file
./create-graph data/Trunc_Test/Trunc_Test.fa.struct.fa --seed > /dev/null
icdiff test/create_graph_golden_master_seed.xgmml $graph_file


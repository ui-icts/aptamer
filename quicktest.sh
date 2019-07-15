./predict-structures samples/Trunc_Test.fa --calculate_stats
./create-graph data/Trunc_Test/Trunc_Test.fa.struct.fa
icdiff data/Trunc_Test/golden_master.xgmml data/Trunc_Test/Trunc_Test.fa.struct.fa.xgmml


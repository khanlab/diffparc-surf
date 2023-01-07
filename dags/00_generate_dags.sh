#!/bin/bash

test_data="../diffparc/test_data/bids_AP_PA"
snakedwi_dir="../diffparc/test_data/snakedwi_dir"
tmp_dir=`mktemp -d`


for dagtype in dag rulegraph filegraph
do

mkdir -p ${dagtype}
poetry run diffparc $test_data $tmp_dir/out participant --${dagtype} | dot -Tpdf > ${dagtype}/${dagtype}_default.pdf

poetry run diffparc $test_data $tmp_dir/out participant --${dagtype} --in-snakedwi-dir $snakedwi_dir --config methods=['mrtrix','fsl'] | dot -Tpdf > ${dagtype}/${dagtype}_fsl_mrtrix.pdf

poetry run diffparc $test_data $tmp_dir/out participant --${dagtype} --in-snakedwi-dir $snakedwi_dir --config methods=['mrtrix'] | dot -Tpdf > ${dagtype}/${dagtype}_mrtrix.pdf

poetry run diffparc $test_data $tmp_dir/out participant --${dagtype} --in-snakedwi-dir $snakedwi_dir --config methods=['fsl'] | dot -Tpdf > ${dagtype}/${dagtype}_fsl.pdf

poetry run diffparc $test_data $tmp_dir/out participant --${dagtype} --use-eddy | dot -Tpdf > ${dagtype}/${dagtype}_use_eddy.pdf

poetry run diffparc $test_data $tmp_dir/out participant --${dagtype} --use-topup --use-eddy | dot -Tpdf > ${dagtype}/${dagtype}_use_topup_eddy.pdf

poetry run diffparc $test_data $tmp_dir/out participant --${dagtype} --in-snakedwi-dir $snakedwi_dir | dot -Tpdf > ${dagtype}/${dagtype}_skip_dwi_preproc.pdf

for metric in surfarea inout indepconn FA 
do
    poetry run diffparc $test_data $tmp_dir/out participant --${dagtype} --in-snakedwi-dir $snakedwi_dir --config surface_metrics=[\'${metric}\'] | dot -Tpdf > ${dagtype}/${dagtype}_${metric}_skip_dwi_preproc.pdf
done


done




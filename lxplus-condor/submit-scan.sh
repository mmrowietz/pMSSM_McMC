#!/bin/bash

npoints=100

subfile=logs/${npoints}.sub
rm $subfile

echo "universe                = docker" >> $subfile
echo "docker_image            = /cvmfs/unpacked.cern.ch/registry.hub.docker.com/jdickins/pmssm-env:latest" >> subfile
echo "executable              = run-scan.sh" >> $subfile
echo "arguments               = ${npoints}" >> $subfile
echo "should_transfer_files   = YES" >> $subfile
echo "when_to_transfer_output = ON_EXIT" >> $subfile
echo "output                  = logs/run_${npoints}.out" >> $subfile
echo "error                   = logs/run_${npoints}.err" >> $subfile
echo "log                     = logs/run_${npoints}.log" >> $subfile
echo "queue" >> $subfile

condor_submit $subfile

done


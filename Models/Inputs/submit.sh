#!/bin/tcsh

#BSUB -J landis_run		# Job name
#BSUB -o output.%J		# Standard output file
#BSUB -e error.%J		# Standard error file
#BSUB -W 48:00			# Maximum time (hh:mm)
#BSUB -n 5			# Number of MPI processes
#BSUB -R "span[hosts=1]"	# Use n cores on 1 node
#BSUB -R "rusage[mem=80GB]"	# Memory requirement (per core)
#BSUB -q sif			# The queue that I want it to join

module load apptainer

apptainer exec --no-mount --containall --cleanenv --bind /share/tcsi/lagoodal/LANDIS_Data/Scenario_Runs/RESISTANCE/RCP45/LSLV/landis:/landis,/share/tcsi/lagoodal/LANDIS_Data/Scenario_Runs/RESISTANCE/RCP45/LSLV/landis/run1:${PWD} /usr/local/usrapps/tcsi/LouisLANDIS/sif8.sif dotnet /Core-Model-v7-LINUX-7/build/Release/Landis.Console.dll /landis/clippedTest_scenario.txt &

apptainer exec --no-mount --containall --cleanenv --bind /share/tcsi/lagoodal/LANDIS_Data/Scenario_Runs/RESISTANCE/RCP45/LSLV/landis:/landis,/share/tcsi/lagoodal/LANDIS_Data/Scenario_Runs/RESISTANCE/RCP45/LSLV/landis/run2:${PWD} /usr/local/usrapps/tcsi/LouisLANDIS/sif8.sif dotnet /Core-Model-v7-LINUX-7/build/Release/Landis.Console.dll /landis/clippedTest_scenario.txt &

apptainer exec --no-mount --containall --cleanenv --bind /share/tcsi/lagoodal/LANDIS_Data/Scenario_Runs/RESISTANCE/RCP45/LSLV/landis:/landis,/share/tcsi/lagoodal/LANDIS_Data/Scenario_Runs/RESISTANCE/RCP45/LSLV/landis/run3:${PWD} /usr/local/usrapps/tcsi/LouisLANDIS/sif8.sif dotnet /Core-Model-v7-LINUX-7/build/Release/Landis.Console.dll /landis/clippedTest_scenario.txt &

apptainer exec --no-mount --containall --cleanenv --bind /share/tcsi/lagoodal/LANDIS_Data/Scenario_Runs/RESISTANCE/RCP45/LSLV/landis:/landis,/share/tcsi/lagoodal/LANDIS_Data/Scenario_Runs/RESISTANCE/RCP45/LSLV/landis/run4:${PWD} /usr/local/usrapps/tcsi/LouisLANDIS/sif8.sif dotnet /Core-Model-v7-LINUX-7/build/Release/Landis.Console.dll /landis/clippedTest_scenario.txt &

apptainer exec --no-mount --containall --cleanenv --bind /share/tcsi/lagoodal/LANDIS_Data/Scenario_Runs/RESISTANCE/RCP45/LSLV/landis:/landis,/share/tcsi/lagoodal/LANDIS_Data/Scenario_Runs/RESISTANCE/RCP45/LSLV/landis/run5:${PWD} /usr/local/usrapps/tcsi/LouisLANDIS/sif8.sif dotnet /Core-Model-v7-LINUX-7/build/Release/Landis.Console.dll /landis/clippedTest_scenario.txt &

wait

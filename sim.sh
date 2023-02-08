#!/usr/bin/env bash

python sim_poiss_poiss.py --repeats 100 --n-events 1000000 --min-freq 1e-4 --max-freq 1e4 --seed 7267

python sim_poiss_pareto.py --repeats 100 --n-events 10000000 --mean-pulse 1e3 --max-gap 1e4 --power-gap 1.0 --min-freq 1e-5 --max-freq 1e1 --seed 1186
python sim_poiss_pareto.py --repeats 100 --n-events 1000000 --pulse-magnitude 1e3 --mean-pulse 1e-3 --max-gap 1e4 --power-gap 1.0 --min-freq 1e-5 --max-freq 1e1 --seed 3981

python sim_uniform_pareto.py --repeats 100 --n-events 1000000 --pulse-magnitude 3 --max-pulse 1e3 --max-gap 1e4 --power-gap 1.0 --min-freq 1e-5 --max-freq 1e1 --seed 15769
python sim_const_pareto.py --repeats 100 --n-events 1000000 --pulse-magnitude 1e-1 --fixed-pulse 1e2 --max-gap 1e4 --power-gap 1.0 --min-freq 1e-5 --max-freq 1e1 --seed 27082
python sim_pareto_pareto.py --repeats 100 --n-events 10000000 --pulse-magnitude 3 --max-pulse 1e4 --power-pulse 1.0 --max-gap 1e4 --power-gap 1.0 --min-freq 1e-5 --max-freq 1e1 --seed 20358

python sim_poiss_pareto_dur.py --repeats 1000 --duration 1e4 --mean-pulse 1e2 --max-gap -1 --power-gap 1.0 --min-freq 1e-4 --max-freq 1e1 --seed 4859
python sim_poiss_pareto_dur.py --repeats 1000 --duration 1e6 --mean-pulse 1e2 --max-gap -1 --power-gap 1.0 --min-freq 1e-6 --max-freq 1e1 --seed 1081
python sim_poiss_pareto_dur.py --repeats 300 --duration 1e8 --mean-pulse 1e2 --max-gap -1 --power-gap 1.0 --min-freq 1e-8 --max-freq 1e1 --seed 22675
python sim_poiss_pareto_dur.py --repeats 10 --duration 1e10 --mean-pulse 1e2 --max-gap -1 --power-gap 1.0 --min-freq 1e-10 --max-freq 1e1 --seed 13007

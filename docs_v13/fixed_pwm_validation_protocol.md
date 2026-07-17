# V13 fixed-PWM validation protocol

V13 resolves the input-identifiability failure observed in V12 by using an
independent experiment with a known fan command.  Guo et al. (2024) report that
the ramp-loading experiments in their Section 3.1 use a fixed 60% PWM.  Their
Fig. 5 gives separate cathode inlet, intermediate, and outlet temperatures for
two input sequences of equal duration.

The 2 A/50 s ramp is calibration data.  Only the three dynamic heat-storage
parameters may change.  The 4 A/100 s ramp is a final confirmation sequence
and cannot enter the objective.  The V12 airflow, heat-transfer, conduction,
heat-source, and fan-footprint parameters remain byte-for-byte traceable to the
V12 frozen JSON.

The delivered plant remains the 4800-state mesh and the fine check remains
19200 states.  Passing V13 establishes a high-fidelity thermal reference only
for 23 °C ambient, 2--40 A, and fixed 60% fan PWM.  It does not validate
variable-duty MPC.  The 55 °C value remains a reference temperature and is not
called a safety limit.

# Cardiac Dynamics
Used a two-variable system to study the voltage flux across cardiac cell membranes in case of Cardiac Arrhythmia and for a regularly functioning heart. The implementation was performed using MATLAB.

## Motivation
Cardiac Arrhythmia, or irregular heart beat, is a name for a large family of cardiac behaviors that show abnormalities in the electrical behavior of the hard. A heartbeat that is too fast (“tachycardia”) or too slow (“bradycardia”) can be fatally dangerous. Other examples of arrhythmias include heart palpitations, stroke, and embolism. The natural pacemaker of the heart is called the sinoatrial node. Pacemaker cells are polarized. When the cells generate an electrical impulse (a cardiac action potential), voltage gated channels open to allow charge to move through the cell, creating heartbeats. These cardiac dynamics can be modeled with the help of a system of nonlinear ordinary differential equations. With the help of MATLAB, required graphs are generated to visualize the problem.

## Initial Equations
We use a two-variable system to study the voltage flux across cardiac cell membranes, which can give insight to how and why arrhythmias form:

<img src="/eq1.png" alt="Eq1" height="200"/>

If h = 0, voltage can pass freely. <br>
If h > 0, then the gate reduces the voltage passing into the cell. <br>
If h is an extremely large value, the gate is essentially closed. <br>
S(t) represents the electrical impulse generated by the pacemaker. <br>
a represents the threshold excitation in the system. <br>
k controls the magnitude of the electric current across the cell membrane. <br>
vh describes the repolarization current in the recovery process. <br>

## The Basic Model
First, we assume that S(t) = 0, now the equation becomes: <br>
<img src="/eq2.png" alt="Eq2" height="130"/>

<img src="/basic_model_v_nullclines.png" alt="v nullclines" width="450"/>
<img src="/basic_model_h_nullclines.png" alt="h nullclines" width="450"/>

There is only one non-negative equilibrium solution (v0,h0). The only non-negative equilibrium point occurs when v0= 0 and h0=0, at the point (0, 0). We will now define the following equations: <br>
v′ = f(v, h) <br>
h′ = g(v, h) <br>

We then calculate the eigenvalues and see that both of them are negative. This implies that the point (0,0) is asymptotically stable. If the solution starts close enough to an asymptotically stable equilibrium point, then the solution will converge to that equilibrium point as t→∞. 

Biologically, this means that the voltage across the cell membrane has a tendency to go towards zero, and the gating variable also approaches zero, in which voltage will
able to pass freely through the gate. <br>

Upon putting in the values of the various constants, the final equations we get are: <br>
<img src="/basic_model_final_eq.png" alt="final eq" width="450"/>

The graph of the nullclines for these equations is as follows: <br>
<img src="/figure1.png" alt="Figure1" height="350"/>

A sample solution curve in the vector field is plotted with the starting point (v0,h0)=(0.5,0.2): <br>
<img src="/sample_curve_basic_model.png" alt="Sample Curve" height="360"/>

A sample solution curve in the vector field, with the starting point (v0,h0)=(0.1,0.2): <br>
<img src="/sample_curve2_basic_model.png" alt="Sample Curve2" height="360"/>

1. The solution curve with starting position (0.5, 0.2) initially moves in the positive v direction, before rotating counter-clockwise and eventually reaching the equilibrium
point (0,0).
2. On the other hand, the solution curve with starting position (0.1, 0.2) initially moves in the negative v direction, and reaches the equilibrium point (0,0) much more quickly. Given the nullclines and the vector field, these solutions makes sense.
3. The trajectories of the solution curves follow the directions of the vectors on the vector field.
4. When the initial value v is to the right of a, the trajectory initially moves in the positive v direction, before rotating counter-clockwise to approach the equilibrium
point.
5. On the other hand, when the initial value v is to the left of a, the trajectory moves in the negative v direction and directly approaches the equilibrium point.
6. As the initial value v gets closer to zero, the trajectory approaches the equilibrium point quicker.

## Model Improvement: Periodic Stimulation
1.  For maintaining a steady heartbeat, periodic stimulation of cardiac cells is required.
2.  The period of stimulation is denoted by the parameter T, which defines the number of time units that pass between each stimulation.
3.  The cell will not be stimulated for any time between stimulation times.
4.  The system is initially assumed to be at (v0,h0)=(0,0) and a large positive stimulus is added to the voltage.

Sample solution curve with initial point (0,0): <br>
<img src="/improved_model_sample_curve.png" alt="Sample Curve" height="360"/>
1.  The solution starts at (v0, h0)=(0,0) and a positive stimulus S=0.25 of voltage is added.
2.   The flow is counter-clockwise. If no more stimuli are given after the initial push, then as t → ∞ the system will approach the equilibrium solution (0,0).

Sample solution curve with initial point (0.25,0): <br>
<img src="/improved_model_sample_curve2.png" alt="Sample Curve2" height="360"/>
1. It basically starts at (0.25,0) and goes to max v of 1, then it rotates counterclockwise approaching the equilibrium point (0,0).
2. The solutions to dv/dt and dh/dt are simulated starting from the initial condition (v0,h0)=(β,0) in which β = 0.25 over the time interval t ∈ [0, 500] with a time-step of Δt=0.2.

Sample solution curve with initial point (0,0) with a stimulation period T=100 and a stimulus of S(t)=0: <br>
<img src="/improved_model_sample_curve3.png" alt="Sample Curve3" height="360"/>
1. The solutions to dv/dt and dh/dt are simulated starting from the initial condition (v0,h0) = (0,0) over the time interval t ∈ [0, 500] with a time-step of Δt = 0.2 with a stimulation period T = 100 and a stimulus of S(t) = 0.25.
2. In the system without constant stimulus, v and h both approach zero. In the system with constant stimulus, v and h approach the equilibrium point zero, but then proceed to increase in magnitude periodically.

## Cardiac Action Potential
The APD (Action Potential Duration) is the duration from the time a cell is stimulated to the time it repolarizes. The APD is calculated as follows:
<img src="/apd_formula.png" alt="APD formula" height="30"/>
1. In which tup is the time at which the voltage v passes a constant critical voltage vc on the way up, and tdown is the time at which the voltage v passes that same constant critical voltage vc on the way down.
2. We let the critical voltage be vc = 0.1.
3. The APD for the last full beat of the v(t) solution curve corresponds to the steady state APD, and is denoted APD0.
4. We proceed to use the initial condition (v0,h0) = (0,0) over the time interval ∈ [0,1000] with a time-step of Δt = 0.2 and stimulus S(t) = 0.25 to stimulate solutions to
dv/dt and dh/dt to find the APD0 for T1= 100, T2=90, T3=80, T4=70, T5=60 and T6=50.

The following is a table of the data points that we have gathered: <br>
<img src="/t_apd_table.png" alt="APD formula" height="260"/>
<img src="/t_apd_graph.png" alt="APD formula" height="260"/> <br>
It is seen that APD0 increases as T increases. From a biological viewpoint, as the heart cell is 0 stimulated more frequently, the time period between consecutive beats decreases.

T vs h: <br>
<img src="/t_h_table.png" alt="APD formula" height="260"/>
<img src="/t_h_graph.png" alt="APD formula" height="260"/> <br>
1. An important feature of cardiac tissue is that theAPD needs to be long enough, especially for large animals.
2. The heartbeats of large animals are longer, because they have larger hearts, and so the voltage gated channels in their hearts need more time to allow charge to move through the cell to create heartbeats.
3. It is seen that the steady-state h decreases as T increases. From a biological point of view, if the heart cell is stimulated less frequently, it will return closer to a resting state, thus will be closer to zero.

APD vs h: <br>
<img src="/apd_h_graph.png" alt="APD formula" height="260"/> <br>
As can be seen, the steady-state h decreases as APD increases. This is because there will be a longer resting period between consecutive heartbeats if a heartbeat takes a longer amount of time, which means that the steady-state h will decrease (since it has more time to reach a resting state).

## Conclusion
1. If a heartbeat takes longer, the heart cell also needs to relax more. This helps explain why a heartbeat that is too fast (“tachycardia”) or too slow (“bradycardia”) can be fatally dangerous.
2. If a heartbeat is too fast, the heart cell may not have enough time to relax. On the other hand, if a heartbeat is too slow, the heart cell may relax too much, since it will tend towards the equilibrium resting state, so it may not be able to recover quickly enough.

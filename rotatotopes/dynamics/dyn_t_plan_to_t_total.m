function dx = dyn_t_plan_to_t_total(tdummy,in2,udummy)
%DYN_T_PLAN_TO_T_TOTAL
%    DX = DYN_T_PLAN_TO_T_TOTAL(TDUMMY,IN2,UDUMMY)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    05-Oct-2020 21:32:46

cqi = in2(1,:);
dcqi = in2(3,:);
dsqi = in2(4,:);
kai = in2(5,:);
kvi = in2(6,:);
sqi = in2(2,:);
t = in2(7,:);
t2 = kai./2.0;
t3 = kvi.*2.0;
t4 = kai+t3;
t5 = t-1.0./2.0;
t7 = t4.*t5;
t6 = kvi+t2-t7;
dx = [-sqi.*t6;cqi.*t6;-dsqi.*t6+sqi.*t4;-cqi.*t4+dcqi.*t6;0.0;0.0;1.0];

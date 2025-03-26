
%{ 
This is the 1st of two MATLAB scripts for particle-based microkinetic 
model (PB-MKM) within a continous-stirred-tank reactor (CSTR) for the case 
of ethane dehydrogenation over Pt. This scripts couples the MKMs for 
Pt(100), Pt(111), and Pt(211) and models cross-facet communication.

This script; 
- provides input parameters for the 2nd script 
- solves the ODE expressions in the 2nd script to compute values of adsorbate species
coverage on the 3 facets and gas mole fraction in the CSTR (y[1-57])
- calculates all kinetic parameters.

Please refer to the particle-based microkinetic modeling paper and SI for 
more details on the derivation of the mathematical expressions, reaction 
network and input parameters including DFT-based energy, lateral interaction 
parameters, and site occupancy for each species on the 3 facets.
%}

format long

%% Reaction conditions
PP_gas = zeros(6,1);                % Partial pressure of 1.CH3CH3(g), 2.CH2CH2(g), 3.CHCH(g), 4.CH4(g), 5.H2(g), 6.Inert(g)
PP_gas(1) = 0.03 ;                  % H2 partial pressure (bar)
PP_gas(5) = 0.2 ;                   % CH3CH3 partial pressure (bar)
PP_gas(6) = 1 - sum(PP_gas(1:5));   % Inert partial pressure (bar)
P_tot = sum(PP_gas);                % Total pressure = 1 (bar)
T = 630.0;                          % Temp (K) 
tauS = 1;                           % Residence time (s)
kB = 8.6173324E-05;                 % Boltzmann's constant (eV/K)
hh = 4.135667516E-15;               % Planck's constant (eV.s)

%% Variable definition
%{
   y[1-17] = Coverage of free site and adsorbate species on (100) facet
   y[18-34] = Coverage of free site and adsorbate species on (111) facet
   y[35-51] = Coverage of free site and adsorbate species on (111) facet
   y[52-57] = Mole fraction of gas species in the CSTR

   Arrangement of free site coverages on all 3 facets.
   y(1)= Free site; y(2)= H; y(3)= CH3CH3; y(4)= CH3CH2; y(5)= CH3CH; y(6)=
   CH3C; y(7)= CH2CH2; y(8)= CH2CH; y(9)= CH2C; y(10)= CHCH; y(11)= CHC; 
   y(12)= CC; y(13)= CH4; y(14)= CH3; y(15)= CH2; y(16)= CH; y(17)= C

   Arrangement of gas mole fractions.
   y(52)= CH3CH3(g); y(53)= CH2CH2(g); y(54)= CHCH(g); y(55)= CH4(g);
   y(56)= H2(g); y(57)= Inert(g)
%}

%% Initial conditions (y0) for coverages and mole fractions (y[1-57])
y0 = zeros(57,1);
y0(1)  = 1;                         % Initial free site coverage on Pt(100)
y0(18) = 1;                         % Initial free site coverage on Pt(111)
y0(35) = 1;                         % Initial free site coverage on Pt(211)
y0(52) = PP_gas(1)/P_tot;           % Initial mole fraction of CH3CH3(g)
y0(56) = PP_gas(5)/P_tot;           % Initial mole fraction of H2(g)
y0(57) = PP_gas(6)/P_tot ;          % Initial mole fraction of Inert(g)

%% Facet fraction at different particle sizes
% facet_frac = [0.001 0.25 0.75];     % 1.0 nm particle size [Pt(100) Pt(111) Pt(211)]
% facet_frac = [0.05 0.46 0.49];      % 1.7 nm particle size
facet_frac = [0.09 0.56 0.35];        % 2.5 nm particle size
% facet_frac = [0.13 0.65 0.22];      % 4.1 nm particle size
% facet_frac = [0.15 0.70 0.14];      % 6.4 nm particle size
% facet_frac = [0.17 0.73 0.10];      % 9.5 nm particle size
% facet_frac = [0.18 0.76 0.06];      % 15.7 nm particle size

XPt100 = facet_frac(1); XPt111 = facet_frac(2); XPt211 = facet_frac(3);
x111 = XPt211/XPt111;  x100 = XPt211/XPt100; 

%% Number of active site per total number of gas molecules in the CSTR (numS)
% numS = 7.6 ;                      % numS for 1.0 nm;  Conversion = ~0.100                                     
% numS = 8.0 ;                      % numS for 1.7 nm;  Conversion = ~0.100                                         
numS = 8.3 ;                        % numS for 2.5 nm;  Conversion = ~0.100
% numS = 8.6 ;                      % numS for 4.1 nm;  Conversion = ~0.100
% numS = 8.8 ;                      % numS for 6.4 nm;  Conversion = ~0.100
% numS = 8.9 ;                      % numS for 9.5 nm;  Conversion = ~0.100
% numS = 9.0 ;                      % numS for 15.7 nm; Conversion = ~0.100

%% Adsorbate, transition state, and gas phase referenced free energies at 873K using FT entropy 
% Energies are as shown and ordered in the Supporting Information
ads = [ 0.7856	0.9769	1.0685	1.0811	1.0090	1.3366	1.7213	0.9339	1.6955	2.7740	-0.0887	0.2908	0.5050	0.5589	0.9893	-0.5731 ...
        0.7352	1.0710	1.4636	0.9705	1.2168	1.6258	1.7330	1.9441	2.8686	4.0231	-0.1447	0.3818	0.9114	0.6480	1.5085	-0.4541 ...
		0.7391	0.9582	0.9868	1.0722	0.9586	1.4835	1.7973	1.6947	1.7895	2.8369	-0.0851	0.2300	0.4664	0.7287	1.5055	-0.5563];

TS  = [2.7565	1.1557	2.1567	1.2513	1.3151	2.5543	1.3899	1.6656	2.5268	1.9401	2.6547	1.4990	2.3848	1.9961	1.2534	3.1481	1.9721	2.4189	1.9185	3.2126	2.8842	4.1465	0.4049	0.5738	0.9717	1.4467 ...
       3.1823	1.4264	2.7438	1.7790	1.6961	2.6367	1.6305	2.0242	2.8577	2.1128	3.1640	1.9118	3.2535	2.0546	2.3404	3.7666	3.0192	2.8949	3.2538	3.8802	4.2873	5.2847	0.6788	1.1243	1.0384	1.8333 ...
       2.4130	1.1357	2.0784	1.1673	1.2589	2.3897	1.5343	1.7836	2.7924	1.9250	2.4754	1.4801	2.7853	1.9913	2.1028	3.4614	2.1213	2.7762	2.7801	3.3274	3.0226	4.3000	0.3918	0.5316	1.0017	1.7722]	;

gas = [0.5734	1.7744	3.5165	-0.3475	-0.5286];
gas_corr = [0.035	0.044	-0.068	-0.023	-0.08];
Gas = gas + gas_corr;

%% Forward rate constant for adsorbate diffusion
kf_SD = zeros(16,2);
% diffusion from (211) to (111)
kf_SD(:,1) = [1.3127E+13	1.2726E+11	8.4577E+07	7.1206E+11	1.0989E+11	7.5051E+08	1.0614E+12	2.2962E+08	2.3880E+03	8.7711E+00	1.3127E+13	1.9031E+11	6.0357E+08	1.2115E+12	5.4932E+11	2.2419E+11];
% diffusion from (211) to (100)
kf_SD(:,2) = [1.3127E+13	8.64581E+11	1.8764E+11	1.31271E+13	1.31271E+13	1.31271E+13	3.64075E+12	1.31271E+13	1.70687E+12	1.94001E+11	1.31271E+13	9.70266E+11	2.40061E+12	1.31271E+13	1.31271E+13	3.22376E+12];

%% Solving ODE equations
% Setting ODE solver tolerance
optode = odeset('NonNegative',1:57,'Abstol',1E-10,'RelTol',1E-10);                      

% Calling ODE solver function (ode15s) and passing arguments (ODE function in script2, initial to final time, initial conditions, and tolerance)
% Solver output is y's at different t's.
[t,y]=ode15s(@(t,y)EH_HO_PBED3_630_PBMKM_func(t,y,Gas,ads,TS,kf_SD,PP_gas,T,facet_frac,tauS,numS),[0,1e+8],y0,optode);

% Plotting the computed coverages and mole fractions as a function of time
loglog(t, y,'b')
xlabel('time t'); ylabel('coverage y');

%% Extracting adsorbate coverages and gas mole fractions (y)
y = y(end,:);

% Free site and adsorbate coverage on the (100) facet
% Computing pseudo coverages
y_100 = [y(1:17) y(34+10)*x100 y(34+11)*x100 y(34+12)*x100] ;
% Computing real coverages (multiplication of pseudo coverages by site occupancies)
y_100_real = (y_100.*[1	 1	1	1	2	4	2	3	3	4	4	4	1	1	2	4	4  2  2  2]) ;
% Computing sum of real coverages
y_100_total = sum(y_100_real)

% Free site and adsorbate coverage on the (111) facet
y_111 = [y(18:34) y(34+6)*x111 y(34+8)*x111 y(34+9)*x111 y(34+16)*x111 y(34+17)*x111] ;
y_111_real = (y_111.*[1	 1	1	1	2	3	2	3	3	3	3	4	1	1	2	3	3  1  1  2  1  1]);
y_111_total = sum(y_111_real)

% Free site and adsorbate coverage on the (211) facet
y_211 = y(35:51) ;
y_211_real = (y_211.*[1	 1	1	1	2	2	2	2	1	2	2	2	1	1	2	2	2]);
y_211_total = sum(y_211_real)

y_site_balance = y_100_total + y_111_total + y_211_total ;
sum(y_site_balance)

% Gas mole fraction in CSTR
MolFraction_total = sum(y(52:57))

%% Substituting y values in script2 to extract rate
[F, r, kf, Kq] = EH_HO_PBED3_630_PBMKM_func(t,y,Gas,ads,TS,kf_SD,PP_gas,T,facet_frac,tauS,numS);	

%% Computing CSTR ethane conversion
R_gas_Pt100 = -(r(1) + r(2) + r(3) + r(4) + r(5));
R_gas_Pt111 = -(r(31+1) + r(31+2) + r(31+3) + r(31+4) + r(31+5));
R_gas_Pt211 = -(r(62+1) + r(62+2) + r(62+3) + r(62+4) + r(62+5));

X_x = 1-((y(52)/(PP_gas(1)/P_tot))*(1 + tauS*numS*(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211)))

%% Computing ethane TOF
TOF_ethane = XPt100*r(1) + XPt111*r(31+1) + XPt211*r(62+1) ;

%% Computing CSTR mole fraction for each gas species
YsS = zeros(5,1);
YsS(1) =  (PP_gas(1)/P_tot - tauS*numS*(XPt100*r(1) + XPt111*r(31+1) + XPt211*r(62+1)))/(1 + tauS*numS*(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211));
YsS(2) =  (PP_gas(2)/P_tot - tauS*numS*(XPt100*r(2) + XPt111*r(31+2) + XPt211*r(62+2)))/(1 + tauS*numS*(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211));
YsS(3) =  (PP_gas(3)/P_tot - tauS*numS*(XPt100*r(3) + XPt111*r(31+3) + XPt211*r(62+3)))/(1 + tauS*numS*(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211));
YsS(4) =  (PP_gas(4)/P_tot - tauS*numS*(XPt100*r(4) + XPt111*r(31+4) + XPt211*r(62+4)))/(1 + tauS*numS*(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211));
YsS(5) =  (PP_gas(5)/P_tot - tauS*numS*(XPt100*r(5) + XPt111*r(31+5) + XPt211*r(62+5)))/(1 + tauS*numS*(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211));

%% Computing yield of gas products
Yld = zeros(3,1);
Yld(1) = (y(53)/(PP_gas(1)/P_tot))*(1 + tauS*numS*(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211));  % CH2CH2(g)
Yld(2) = (y(54)/(PP_gas(1)/P_tot))*(1 + tauS*numS*(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211));  % CHCH(g)
Yld(3) = (y(55)/(PP_gas(1)/P_tot))*(1 + tauS*numS*(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211));  % CH4(g)

%% Computing selectivity to gas products
Sel = zeros(1,3);
Sel(1) = Yld(1)/X_x;        % CH2CH2(g)   
Sel(2) = Yld(2)/X_x;        % CHCH(g) 
Sel(3) = 0.5*Yld(3)/X_x;    % CH4(g) 


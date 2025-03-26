
%{ 
This is the 2nd of two MATLAB scripts for particle-based microkinetic model
(PB-MKM) within a continous-stirred-tank reactor (CSTR) for the case of 
ethane dehydrogenation over Pt. This scripts couples the MKMs for Pt(100), 
Pt(111), and Pt(211) and models cross-facet communication.

This script is a function that; 
- takes input arguments from the 1st script
- outputs ODE expressions (F[1:57]) for the rate of change in
adsorbate species coverage on the 3 facets and gas mole fractions in the 
CSTR (dy/dt).

Please refer to the particle-based microkinetic modeling paper and SI for 
more details on the derivation of the mathematical expressions, reaction 
network and input parameters including DFT-based energy, lateral interaction 
parameters, and site occupancy for each species on the 3 facets.
%}

function [F, r, kf, Kq] = EH_HO_PBED3_630_PBMKM_func(t,y,Gas,ads,TS,kf_SD,PP_gas,T,facet_frac,tauS,numS)
    
    %% Constants
    P_tot = sum(PP_gas);            % Total pressure (bar)
    kB = 8.6173324E-05;             % Boltzmann's constant (eV/K)
    hh = 4.135667516E-15;           % Planck's constant (eV.s)
    XPt100 = facet_frac(1);         % Pt(100) facet fraction
    XPt111 = facet_frac(2);         % Pt(111) facet fraction
    XPt211 = facet_frac(3);         % Pt(211) facet fraction
    x111 = XPt211/XPt111;  
    x100 = XPt211/XPt100;  

	%% Lateral interaction parameter
    LI_100_CHCH = [-0.0785	0.1851	0.2610	0.1616	0.2408	0.3111	0.2657	0.3677	0.3959	0.4738	0.0322	0.1555	0.2289	0.1474	0.0933	0.1061];
    LI_100_CHC = [0.0095	0.1262	0.1438	0.0245	0.1909	0.2331	0.1532	0.3083	0.3901	0.5879	-0.0130	0.1158	0.1500	0.0357	0.1083	0.0645];
    LI_100_H = [-0.026	0.135	0.166	-0.114	-0.002	0.002	0.039	0.049	0.162	0.278	-0.002	0.082	0.052	-0.076	-0.032	0.047] ;
    LI_111_CH3C = [0.0304	0.3055	0.4277	0.5996	0.3407	0.7770	0.3871	0.8205	0.3895	0.7610	-0.0037	0.211	0.3412	0.5586	0.5608	0.187];
    %LI_111_H = [0.479	0.532	0.757	0.757	0.614	0.768	0.857	0.894	0.938	1.032	0.288	0.435	0.590	0.783	0.999	0.470] ;
    LI_211_CHC = [0.0003	0.0790	0.1108	0.1313	0.0784	0.2540	0.3076	-0.0714	0.0238	-0.0266	-0.0079	0.0673	0.0736	0.0606	0.1194	0.0436];
    LI_211_CH3C = [0.1690	0.0189	0.0485	0.0787	0.0496	0.0598	0.0671	0.2075	0.1772	0.0706	0.0795	0.0314	0.0176	0.0196	-0.0456	0.0041];
    LI_211_CH3C_other = [0.3379	0.0378	0.0969	0.1574	0.0993	0.1197	0.1343	0.4150	0.3544	0.1413	0.1589	0.0629	0.0351	0.0393	-0.0913	0.0082];
    LI_211_CHC_other = [0.0007	0.1580	0.2215	0.2627	0.1568	0.5080	0.6152	-0.1428	0.0475	-0.0532	-0.0158	0.1347	0.1472	0.1212	0.2388	0.0872];
    LI_211_H = [0.00974	0.0282	0.06178	0.06632	0.19162	0.0343	0.00692	0.096	0.06778	0.07976	0.03592	0.03726	0.0928	0.09878	0.14682	-0.01788] ;

    %% Forward rate constant computed from collision theory (1/s)
	f_100 = [146987950.4 152188562.7 157953028.9 201254808.8 569945659.3];
    f_111 = [127309091.3 131813441.6 136806156.7 174310661.2 493640898.9];
    f_211 = [159077847.7 164706215.2 170944814.2 217808206.3 616824226.2];
	
	%% Computing adsorbate free energy with lateral interaction
	sp = zeros(48,1);
    for i = 1:16
	    sp(i) = ads(i) + LI_100_CHCH(i)*y(10)*4  + LI_100_CHC(i)*y(11)*4 + LI_100_H(i)*y(2) ;
	    sp(16+i) = ads(16+i) + LI_111_CH3C(i)*y(17+6)*3 ;%+ LI_111_H(i)*y(17+2) ;
 	    sp(32+i) = ads(32+i) + LI_211_H(i)*y(34+2) + LI_211_CHC(i)*y(34+11)*2 + LI_211_CH3C(i)*y(34+6)*2 ...
            + LI_211_CHC_other(i)*y(34+11)*x100 + LI_211_CH3C_other(i)*y(34+6)*x111 ;
    end

    %% Computing reaction free energy 
    % Elementary reactions are as shown and ordered in the Supporting Information 
	 function rxn = reactionEquations(chemSpe)
	    rexn = zeros(26,1);
	    rxt = chemSpe;
	    rexn(1) = rxt(12) + rxt(12) - rxt(1) ;
	    rexn(2) = rxt(2) + rxt(16) - rxt(1) ;
	    rexn(3) = rxt(12) + rxt(13) - rxt(2) ;
	    rexn(4) = rxt(3) + rxt(16) - rxt(2) ;
		rexn(5) = rxt(5) + rxt(16) - rxt(2) ;
		rexn(6) = rxt(12) + rxt(14) - rxt(3) ;
		rexn(7) = rxt(4) + rxt(16) - rxt(3) ;
		rexn(8) = rxt(6) + rxt(16) - rxt(3) ;
		rexn(9) = rxt(12) + rxt(15) - rxt(4) ;
		rexn(10) = rxt(7) + rxt(16) - rxt(4) ;
		rexn(11) = rxt(13) + rxt(13) - rxt(5) ;
		rexn(12) = rxt(6) + rxt(16) - rxt(5) ;
		rexn(13) = rxt(13) + rxt(14) - rxt(6) ;
		rexn(14) = rxt(7) + rxt(16) - rxt(6) ;
		rexn(15) = rxt(8) + rxt(16) - rxt(6) ;
		rexn(16) = rxt(13) + rxt(15) - rxt(7) ;
		rexn(17) = rxt(9) + rxt(16) - rxt(7) ;
		rexn(18) = rxt(14) + rxt(14) - rxt(8) ;
		rexn(19) = rxt(9) + rxt(16) - rxt(8) ;
		rexn(20) = rxt(14) + rxt(15) - rxt(9) ;
		rexn(21) = rxt(10) + rxt(16) - rxt(9) ;
		rexn(22) = rxt(15) + rxt(15) - rxt(10) ;
		rexn(23) = rxt(12) + rxt(16) - rxt(11) ;
		rexn(24) = rxt(13) + rxt(16) - rxt(12) ;
		rexn(25) = rxt(14) + rxt(16) - rxt(13) ;
		rexn(26) = rxt(15) + rxt(16) - rxt(14) ;
		
		rxn = rexn ;
	end

	rxn_SCF_100 = reactionEquations(ads(1:16));
	rxn_SCF_111 = reactionEquations(ads(17:32));
	rxn_SCF_211 = reactionEquations(ads(33:48));
	
	rxn_egy_100 = reactionEquations(sp(1:16));
	rxn_egy_111 = reactionEquations(sp(17:32));
	rxn_egy_211 = reactionEquations(sp(33:48));	
	
	rxn_SCF = [rxn_SCF_100; rxn_SCF_111; rxn_SCF_211];
	rxn_egy = [rxn_egy_100; rxn_egy_111; rxn_egy_211];
	
	%% Computing reaction free energy barriers with lateral interaction
    % Elementary reactions are as shown and ordered Supporting Information
	function rxnBarr = reactionBarriers(TS, rxn_egy, rxn_SCF, ads)
		rxn_barr = zeros(26,1);
		rxn_barr(1)  = TS(1) + 0.5*(rxn_egy(1) - rxn_SCF(1)) - ads(1);
		rxn_barr(2)  = TS(2) + 0.5*(rxn_egy(2) - rxn_SCF(2)) - ads(1);
		rxn_barr(3)  = TS(3) + 0.5*(rxn_egy(3) - rxn_SCF(3)) - ads(2);
		rxn_barr(4)  = TS(4) + 0.5*(rxn_egy(4) - rxn_SCF(4)) - ads(2);
		rxn_barr(5)  = TS(5) + 0.5*(rxn_egy(5) - rxn_SCF(5)) - ads(2);
		rxn_barr(6)  = TS(6) + 0.5*(rxn_egy(6) - rxn_SCF(6)) - ads(3);
		rxn_barr(7)  = TS(7) + 0.5*(rxn_egy(7) - rxn_SCF(7)) - ads(3);
		rxn_barr(8)  = TS(8) + 0.5*(rxn_egy(8) - rxn_SCF(8)) - ads(3);
		rxn_barr(9)  = TS(9) + 0.5*(rxn_egy(9) - rxn_SCF(9)) - ads(4);
		rxn_barr(10) = TS(10) + 0.5*(rxn_egy(10) - rxn_SCF(10)) - ads(4);
		rxn_barr(11) = TS(11) + 0.5*(rxn_egy(11) - rxn_SCF(11)) - ads(5);
		rxn_barr(12) = TS(12) + 0.5*(rxn_egy(12) - rxn_SCF(12)) - ads(5);
		rxn_barr(13) = TS(13) + 0.5*(rxn_egy(13) - rxn_SCF(13)) - ads(6);
		rxn_barr(14) = TS(14) + 0.5*(rxn_egy(14) - rxn_SCF(14)) - ads(6);
		rxn_barr(15) = TS(15) + 0.5*(rxn_egy(15) - rxn_SCF(15)) - ads(6);
		rxn_barr(16) = TS(16) + 0.5*(rxn_egy(16) - rxn_SCF(16)) - ads(7);
		rxn_barr(17) = TS(17) + 0.5*(rxn_egy(17) - rxn_SCF(17)) - ads(7);
		rxn_barr(18) = TS(18) + 0.5*(rxn_egy(18) - rxn_SCF(18)) - ads(8);
		rxn_barr(19) = TS(19) + 0.5*(rxn_egy(19) - rxn_SCF(19)) - ads(8);
		rxn_barr(20) = TS(20) + 0.5*(rxn_egy(20) - rxn_SCF(20)) - ads(9);
		rxn_barr(21) = TS(21) + 0.5*(rxn_egy(21) - rxn_SCF(21)) - ads(9);
		rxn_barr(22) = TS(22) + 0.5*(rxn_egy(22) - rxn_SCF(22)) - ads(10);
		rxn_barr(23) = TS(23) + 0.5*(rxn_egy(23) - rxn_SCF(23)) - ads(11);
		rxn_barr(24) = TS(24) + 0.5*(rxn_egy(24) - rxn_SCF(24)) - ads(12);
		rxn_barr(25) = TS(25) + 0.5*(rxn_egy(25) - rxn_SCF(25)) - ads(13);
		rxn_barr(26) = TS(26) + 0.5*(rxn_egy(26) - rxn_SCF(26)) - ads(14);
		
        for j = 1:26
			if (rxn_barr(j) < rxn_egy(j))
				rxn_barr(j) = rxn_egy(j) ; 
			end
			
			if (rxn_barr(j) < 0)
				rxn_barr(j) = 0.001 ;
			end
        
        end

		rxnBarr = rxn_barr ;
		
	end
	
	rxnBarrier100 = reactionBarriers(TS(1:26), rxn_egy(1:26), rxn_SCF(1:26), ads(1:16));
	rxnBarrier111 = reactionBarriers(TS(27:52), rxn_egy(27:52), rxn_SCF(27:52), ads(17:32));
	rxnBarrier211 = reactionBarriers(TS(53:78), rxn_egy(53:78), rxn_SCF(53:78), ads(33:48));
	
	rxnBar = [rxnBarrier100; rxnBarrier111; rxnBarrier211];
	
	%% Computing Equilibrium constant
	function K_surf = EquilibriumConstant(chemSpe, rxnEnergy)
	
		K_surf = zeros(31,1);
		sps = chemSpe ;
		rxnEgy = rxnEnergy;

	% Adsorption/desorption reaction	
		K_surf(1) = exp(-(sps(1) - Gas(1))/(T*kB)) ;
		K_surf(2) = exp(-(sps(5) - Gas(2))/(T*kB)) ;
		K_surf(3) = exp(-(sps(8) - Gas(3))/(T*kB)) ;
		K_surf(4) = exp(-(sps(11) - Gas(4))/(T*kB)) ;
		K_surf(5) = exp(-(2*sps(16) - Gas(5))/(T*kB)) ;

    % Surface reaction
		for k = 1:26
			K_surf(5+k) = exp(-rxnEgy(k)/(T*kB)) ;
		end
	end
	
	Eqconst_100 = EquilibriumConstant(sp(1:16), rxn_egy(1:26));
	Eqconst_111 = EquilibriumConstant(sp(17:32), rxn_egy(27:52));
	Eqconst_211 = EquilibriumConstant(sp(33:48), rxn_egy(53:78));
	
	Kq = [Eqconst_100; Eqconst_111; Eqconst_211];
	
    %% Computing forward rate constant for surface reactions
	function f_surf = Kforward(fCollision, rxnBarrier)
		f_surf = zeros(31,1);
		
		for m=1:5
			f_surf(m) = fCollision(m) ;
		end
		
		for n=1:26
			f_surf(5+n) = (kB*T/hh)*exp(-(rxnBarrier(n)/(kB*T))) ;
		end		
	
	end
	
	kf_100 = Kforward(f_100, rxnBar(1:26)) ;
	kf_111 = Kforward(f_111, rxnBar(27:52)) ;	
	kf_211 = Kforward(f_211, rxnBar(53:78)) ;	
	
	kf = [kf_100; kf_111; kf_211] ;

    %% Computing backward rate constant
	kb = kf./Kq ;
	
    %% Rate equations on Pt(100)
	r = zeros(93,1);
	
	r(1)  = kf(1)*P_tot*y(52)*y(1) 					- kb(1)*y(3); 						
    r(2)  = kf(2)*P_tot*y(53)*y(1)*y(1) 			- kb(2)*y(7); 					
	r(3)  = kf(3)*P_tot*y(54)*y(1)*y(1)*y(1)*y(1) 	- kb(3)*y(10); 		
	r(4)  = kf(4)*P_tot*y(55)*y(1) 					- kb(4)*y(13); 						
	r(5)  = kf(5)*P_tot*y(56)*y(1)*y(1) 			- kb(5)*y(2)*y(2); 					
	r(6)  = kf(6)*y(3)*y(1) 						- kb(6)*y(14)*y(14); 					
	r(7)  = kf(7)*y(3)*y(1) 						- kb(7)*y(4)*y(2);  						
	r(8)  = kf(8)*y(4)*y(1)*y(1) 					- kb(8)*y(14)*y(15); 				
    r(9)  = kf(9)*y(4)*y(1)*y(1) 					- kb(9)*y(5)*y(2);  				
    r(10) = kf(10)*y(4)*y(1)*y(1) 					- kb(10)*y(7)*y(2); 				
    r(11) = kf(11)*y(5)*y(1)*y(1)*y(1) 				- kb(11)*y(14)*y(16); 		
    r(12) = kf(12)*y(5)*y(1)*y(1)*y(1) 				- kb(12)*y(6)*y(2); 			
    r(13) = kf(13)*y(5)*y(1)*y(1) 					- kb(13)*y(8)*y(2); 				
    r(14) = kf(14)*y(6)*y(1) 						- kb(14)*y(14)*y(17); 					
    r(15) = kf(15)*y(6) 							- kb(15)*y(9)*y(2); 						
    r(16) = kf(16)*y(7)*y(1)*y(1) 					- kb(16)*y(15)*y(15); 			
    r(17) = kf(17)*y(7)*y(1)*y(1) 					- kb(17)*y(8)*y(2); 				
    r(18) = kf(18)*y(8)*y(1)*y(1)*y(1) 				- kb(18)*y(15)*y(16);  	    
    r(19) = kf(19)*y(8)*y(1) 						- kb(19)*y(9)*y(2); 				    
    r(20) = kf(20)*y(8)*y(1)*y(1) 					- kb(20)*y(10)*y(2); 		        
    r(21) = kf(21)*y(9)*y(1)*y(1)*y(1) 				- kb(21)*y(15)*y(17); 		
    r(22) = kf(22)*y(9)*y(1)*y(1) 					- kb(22)*y(11)*y(2); 				
    r(23) = kf(23)*y(10)*y(1)*y(1)*y(1)*y(1) 		- kb(23)*y(16)*y(16); 	
    r(24) = kf(24)*y(10)*y(1) 						- kb(24)*y(11)*y(2); 					
    r(25) = kf(25)*y(11)*y(1)*y(1)*y(1)*y(1) 		- kb(25)*y(16)*y(17); 	
    r(26) = kf(26)*y(11)*y(1) 						- kb(26)*y(12)*y(2); 					
    r(27) = kf(27)*y(12)*y(1)*y(1)*y(1)*y(1) 		- kb(27)*y(17)*y(17); 	
    r(28) = kf(28)*y(13)*y(1) 						- kb(28)*y(14)*y(2);    				
    r(29) = kf(29)*y(14)*y(1)*y(1) 					- kb(29)*y(15)*y(2); 			
    r(30) = kf(30)*y(15)*y(1)*y(1)*y(1) 			- kb(30)*y(16)*y(2);	    
    r(31) = kf(31)*y(16)*y(1) 						- kb(31)*y(17)*y(2); 						

	%% Rate equations on Pt(111)	
	r(31+1)  = kf(31+1)*P_tot*y(52)*y(17+1) 					- kb(31+1)*y(17+3); 					
	r(31+2)  = kf(31+2)*P_tot*y(53)*y(17+1)*y(17+1) 			- kb(31+2)*y(17+7); 				
	r(31+3)  = kf(31+3)*P_tot*y(54)*y(17+1)*y(17+1)*y(17+1)		- kb(31+3)*y(17+10); 		    
	r(31+4)  = kf(31+4)*P_tot*y(55)*y(17+1) 					- kb(31+4)*y(17+13); 					
	r(31+5)  = kf(31+5)*P_tot*y(56)*y(17+1)*y(17+1) 			- kb(31+5)*y(17+2)*y(17+2); 				
	r(31+6)  = kf(31+6)*y(17+3)*y(17+1) 						- kb(31+6)*y(17+14)*y(17+14); 				
	r(31+7)  = kf(31+7)*y(17+3)*y(17+1) 						- kb(31+7)*y(17+4)*y(17+2);  					
	r(31+8)  = kf(31+8)*y(17+4)*y(17+1)*y(17+1) 				- kb(31+8)*y(17+14)*y(17+15); 			
	r(31+9)  = kf(31+9)*y(17+4)*y(17+1)*y(17+1) 				- kb(31+9)*y(17+5)*y(17+2);  			
	r(31+10) = kf(31+10)*y(17+4)*y(17+1)*y(17+1) 				- kb(31+10)*y(17+7)*y(17+2); 			
	r(31+11) = kf(31+11)*y(17+5)*y(17+1)*y(17+1) 				- kb(31+11)*y(17+14)*y(17+16); 		
	r(31+12) = kf(31+12)*y(17+5)*y(17+1)*y(17+1) 				- kb(31+12)*y(17+6)*y(17+2); 			
	r(31+13) = kf(31+13)*y(17+5)*y(17+1)*y(17+1) 				- kb(31+13)*y(17+8)*y(17+2); 			
	r(31+14) = kf(31+14)*y(17+6)*y(17+1) 						- kb(31+14)*y(17+14)*y(17+17); 				
	r(31+15) = kf(31+15)*y(17+6)*y(17+1) 						- kb(31+15)*y(17+9)*y(17+2); 				
	r(31+16) = kf(31+16)*y(17+7)*y(17+1)*y(17+1) 				- kb(31+16)*y(17+15)*y(17+15); 		
	r(31+17) = kf(31+17)*y(17+7)*y(17+1)*y(17+1) 				- kb(31+17)*y(17+8)*y(17+2); 			
	r(31+18) = kf(31+18)*y(17+8)*y(17+1)*y(17+1) 				- kb(31+18)*y(17+15)*y(17+16);  	    
	r(31+19) = kf(31+19)*y(17+8)*y(17+1) 						- kb(31+19)*y(17+9)*y(17+2); 				
	r(31+20) = kf(31+20)*y(17+8)*y(17+1) 						- kb(31+20)*y(17+10)*y(17+2); 		        
	r(31+21) = kf(31+21)*y(17+9)*y(17+1)*y(17+1) 				- kb(31+21)*y(17+15)*y(17+17); 		
	r(31+22) = kf(31+22)*y(17+9)*y(17+1) 						- kb(31+22)*y(17+11)*y(17+2); 				
	r(31+23) = kf(31+23)*y(17+10)*y(17+1)*y(17+1)*y(17+1) 		- kb(31+23)*y(17+16)*y(17+16); 	
	r(31+24) = kf(31+24)*y(17+10)*y(17+1) 						- kb(31+24)*y(17+11)*y(17+2); 			    
	r(31+25) = kf(31+25)*y(17+11)*y(17+1)*y(17+1)*y(17+1) 		- kb(31+25)*y(17+16)*y(17+17); 	
	r(31+26) = kf(31+26)*y(17+11)*y(17+1)*y(17+1) 				- kb(31+26)*y(17+12)*y(17+2); 		
	r(31+27) = kf(31+27)*y(17+12)*y(17+1)*y(17+1) 				- kb(31+27)*y(17+17)*y(17+17); 	    
	r(31+28) = kf(31+28)*y(17+13)*y(17+1) 						- kb(31+28)*y(17+14)*y(17+2);    			
	r(31+29) = kf(31+29)*y(17+14)*y(17+1)*y(17+1) 				- kb(31+29)*y(17+15)*y(17+2); 		
	r(31+30) = kf(31+30)*y(17+15)*y(17+1)*y(17+1) 				- kb(31+30)*y(17+16)*y(17+2);	        
	r(31+31) = kf(31+31)*y(17+16)*y(17+1) 						- kb(31+31)*y(17+17)*y(17+2); 	

	%% Rate equations on Pt(211)	
	r(62+1)  = kf(62+1)*P_tot*y(52)*y(34+1) 					- kb(62+1)*y(34+3); 								
	r(62+2)  = kf(62+2)*P_tot*y(53)*y(34+1)*y(34+1) 			- kb(62+2)*y(34+7); 							
	r(62+3)  = kf(62+3)*P_tot*y(54)*y(34+1)*y(34+1)*y(1)*y(1) 	- kb(62+3)*y(34+10); 		    
	r(62+4)  = kf(62+4)*P_tot*y(55)*y(34+1) 					- kb(62+4)*y(34+13); 								
	r(62+5)  = kf(62+5)*P_tot*y(56)*y(34+1)*y(34+1) 			- kb(62+5)*y(34+2)*y(34+2); 					
	r(62+6)  = kf(62+6)*y(34+3)*y(34+1) 						- kb(62+6)*y(34+14)*y(34+14); 								
	r(62+7)  = kf(62+7)*y(34+3)*y(34+1) 						- kb(62+7)*y(34+4)*y(34+2);  								
	r(62+8)  = kf(62+8)*y(34+4)*y(34+1)*y(34+1) 				- kb(62+8)*y(34+14)*y(34+15); 							
	r(62+9)  = kf(62+9)*y(34+4)*y(34+1)*y(34+1) 				- kb(62+9)*y(34+5)*y(34+2);  							
	r(62+10) = kf(62+10)*y(34+4)*y(34+1)*y(34+1) 				- kb(62+10)*y(34+7)*y(34+2); 							
	r(62+11) = kf(62+11)*y(34+5)*y(34+1)*y(18) 					- kb(62+11)*y(34+14)*y(34+16); 		    			
	r(62+12) = kf(62+12)*y(34+5)*y(34+1)*y(18) 					- kb(62+12)*y(34+6)*y(34+2); 			    		
	r(62+13) = kf(62+13)*y(34+5)*y(34+1)*y(18) 					- kb(62+13)*y(34+8)*y(34+2); 						
	r(62+14) = kf(62+14)*y(34+6)*y(34+1) 						- kb(62+14)*y(34+14)*y(34+17); 							
	r(62+15) = kf(62+15)*y(34+6)*y(18) 						    - kb(62+15)*y(34+9)*y(34+2); 								
	r(62+16) = kf(62+16)*y(34+7)*y(34+1)*y(34+1) 				- kb(62+16)*y(34+15)*y(34+15); 						
	r(62+17) = kf(62+17)*y(34+7)*y(34+1)*y(18) 					- kb(62+17)*y(34+8)*y(34+2); 						
	r(62+18) = kf(62+18)*y(34+8)*y(34+1)*y(34+1) 				- kb(62+18)*y(34+15)*y(34+16);  	        			
	r(62+19) = kf(62+19)*y(34+8)*y(18) 							- kb(62+19)*y(34+9)*y(34+2); 				    			
	r(62+20) = kf(62+20)*y(34+8)*y(34+1)*y(1)*y(1) 				- kb(62+20)*y(34+10)*y(34+2)*y(18); 		    
	r(62+21) = kf(62+21)*y(34+9)*y(34+1)*y(34+1)*y(34+1) 		- kb(62+21)*y(34+15)*y(34+17)*y(18); 		    
	r(62+22) = kf(62+22)*y(34+9)*y(34+1)*y(34+1)*y(1)*y(1) 		- kb(62+22)*y(34+11)*y(34+2)*y(18)*y(18); 					
	r(62+23) = kf(62+23)*y(34+10)*y(34+1)*y(34+1)*y(18)*y(18) 	- kb(62+23)*y(34+16)*y(34+16)*y(1)*y(1);
	r(62+24) = kf(62+24)*y(34+10)*y(34+1) 			    		- kb(62+24)*y(34+11)*y(34+2); 			
	r(62+25) = kf(62+25)*y(34+11)*y(34+1)*y(34+1)*y(18)*y(18)	- kb(62+25)*y(34+16)*y(34+17)*y(1)*y(1); 	        	
	r(62+26) = kf(62+26)*y(34+11)*y(34+1) 						- kb(62+26)*y(34+12)*y(34+2); 							
	r(62+27) = kf(62+27)*y(34+12)*y(34+1)*y(34+1)*y(18)*y(18) 	- kb(62+27)*y(34+17)*y(34+17)*y(1)*y(1); 	        	
	r(62+28) = kf(62+28)*y(34+13)*y(34+1) 						- kb(62+28)*y(34+14)*y(34+2);    							
	r(62+29) = kf(62+29)*y(34+14)*y(34+1)*y(34+1) 				- kb(62+29)*y(34+15)*y(34+2); 						
	r(62+30) = kf(62+30)*y(34+15)*y(34+1)*y(18) 				- kb(62+30)*y(34+16)*y(34+2);	            		
	r(62+31) = kf(62+31)*y(34+16)*y(34+1) 						- kb(62+31)*y(34+17)*y(34+2); 			

 	%% Rate equations for surface diffusion from (211) to Pt(111)  
    r(94)  = kf_SD(1,1)*y(34+3)*XPt211*y(17+1)                               - kf_SD(1,1)/exp(-(sp(16+1)-sp(32+1))/(kB*T))*y(17+3)*XPt211*y(34+1);                               % CH3CH3
    r(95)  = kf_SD(2,1)*y(34+4)*XPt211*y(17+1)                               - kf_SD(2,1)/exp(-(sp(16+2)-sp(32+2))/(kB*T))*y(17+4)*XPt211*y(34+1);                               % CH3CH2
    r(96)  = kf_SD(3,1)*y(34+5)*XPt211*y(17+1)*y(17+1)                       - kf_SD(3,1)/exp(-(sp(16+3)-sp(32+3))/(kB*T))*y(17+5)*XPt211*y(34+1)*y(34+1);                    	% CH3CH 
    r(97)  = kf_SD(4,1)*y(34+6)*XPt211*y(17+1)*y(17+1)                       - kf_SD(4,1)/exp(-(sp(16+4)-sp(32+4))/(kB*T))*y(17+6)*XPt211*y(34+1)*y(34+1);                    	% CH3C  
    r(98)  = kf_SD(5,1)*y(34+7)*XPt211*y(17+1)*y(17+1)                       - kf_SD(5,1)/exp(-(sp(16+5)-sp(32+5))/(kB*T))*y(17+7)*XPt211*y(34+1)*y(34+1);                    	% CH2CH2
    r(99)  = kf_SD(6,1)*y(34+8)*XPt211*y(17+1)*y(17+1)                       - kf_SD(6,1)/exp(-(sp(16+6)-sp(32+6))/(kB*T))*y(17+8)*XPt211*y(34+1)*y(34+1);                    	% CH2CH 
    r(100) = kf_SD(7,1)*y(34+9)*XPt211*y(17+1)                               - kf_SD(7,1)/exp(-(sp(16+7)-sp(32+7))/(kB*T))*y(17+9)*XPt211*y(34+1)  ;             				% CH2C  
    r(101) = kf_SD(8,1)*y(34+10)*XPt211*y(17+1)*y(17+1)*y(17+1)          	 - kf_SD(8,1)/exp(-(sp(16+8)-sp(32+8))/(kB*T))*y(17+10)*XPt211*y(34+1)*y(34+1)*y(1)*y(1) ;     		% CHCH  
    r(102) = kf_SD(9,1)*y(34+11)*XPt211*y(17+1)*y(17+1)*y(17+1)          	 - kf_SD(9,1)/exp(-(sp(16+9)-sp(32+9))/(kB*T))*y(17+11)*XPt211*y(34+1)*y(34+1)*y(1)*y(1) ;     		% CHC   
    r(103) = kf_SD(10,1)*y(34+12)*XPt211*y(17+1)*y(17+1)*y(17+1)*y(17+1) 	 - kf_SD(10,1)/exp(-(sp(16+10)-sp(32+10))/(kB*T))*y(17+12)*XPt211*y(34+1)*y(34+1)*y(1)*y(1);      	% CC    
    r(104) = kf_SD(11,1)*y(34+13)*XPt211*y(17+1)                             - kf_SD(11,1)/exp(-(sp(16+11)-sp(32+11))/(kB*T))*y(17+13)*XPt211*y(34+1)  ;                       	% CH4   
    r(105) = kf_SD(12,1)*y(34+14)*XPt211*y(17+1)                             - kf_SD(12,1)/exp(-(sp(16+12)-sp(32+12))/(kB*T))*y(17+14)*XPt211*y(34+1)  ;                        	% CH3   
    r(106) = kf_SD(13,1)*y(34+15)*XPt211*y(17+1)*y(17+1)                     - kf_SD(13,1)/exp(-(sp(16+13)-sp(32+13))/(kB*T))*y(17+15)*XPt211*y(34+1)*y(34+1);                    % CH2   
    r(107) = kf_SD(14,1)*y(34+16)*XPt211*y(17+1)*y(17+1)                     - kf_SD(14,1)/exp(-(sp(16+14)-sp(32+14))/(kB*T))*y(17+16)*XPt211*y(34+1)*y(34+1);                    % CH    
    r(108) = kf_SD(15,1)*y(34+17)*XPt211*y(17+1)*y(17+1)                     - kf_SD(15,1)/exp(-(sp(16+15)-sp(32+15))/(kB*T))*y(17+17)*XPt211*y(34+1)*y(34+1);                    % C	   
    r(109) = kf_SD(16,1)*y(34+2)*XPt211 *y(17+1)                             - kf_SD(16,1)/exp(-(sp(16+16)-sp(32+16))/(kB*T))*y(17+2)*XPt211*y(34+1);                             % H     
    
    %% Rate equations for surface diffusion from (211) to Pt(100) 
    r(110) =  kf_SD(1,2)*y(34+3)*XPt211*y(1)                         -  kf_SD(1,2)/exp(-(sp(1)-sp(32+1))/(kB*T))*y(3)*XPt211*y(34+1);                        % CH3CH3
    r(111) =  kf_SD(2,2)*y(34+4)*XPt211*y(1)                         -  kf_SD(2,2)/exp(-(sp(2)-sp(32+2))/(kB*T))*y(4)*XPt211*y(34+1);                        % CH3CH2
    r(112) =  kf_SD(3,2)*y(34+5)*XPt211*y(1)*y(1)              	     -  kf_SD(3,2)/exp(-(sp(3)-sp(32+3))/(kB*T))*y(5)*XPt211*y(34+1)*y(34+1);          	    % CH3CH 
    r(113) =  kf_SD(4,2)*y(34+6)*XPt211*y(1)*y(1)*y(1)*y(1)          -  kf_SD(4,2)/exp(-(sp(4)-sp(32+4))/(kB*T))*y(6)*XPt211*y(34+1)*y(34+1)*y(18);    	    % CH3C  
    r(114) =  kf_SD(5,2)*y(34+7)*XPt211*y(1)*y(1)              	     -  kf_SD(5,2)/exp(-(sp(5)-sp(32+5))/(kB*T))*y(7)*XPt211*y(34+1)*y(34+1);           	    % CH2CH2
    r(115) =  kf_SD(6,2)*y(34+8)*XPt211*y(1)*y(1)*y(1)        	     -  kf_SD(6,2)/exp(-(sp(6)-sp(32+6))/(kB*T))*y(8)*XPt211*y(34+1)*y(34+1)*y(18);    	    % CH2CH 
    r(116) =  kf_SD(7,2)*y(34+9)*XPt211*y(1)*y(1)*y(1)        	     -  kf_SD(7,2)/exp(-(sp(7)-sp(32+7))/(kB*T))*y(9)*XPt211*y(34+1)*y(18)*y(18);    	    % CH2C  
    r(117) =  kf_SD(8,2)*y(34+10)*XPt211*y(1)*y(1)                   -  kf_SD(8,2)/exp(-(sp(8)-sp(32+8))/(kB*T))*y(10)*XPt211*y(34+1)*y(34+1);          	    % CHCH  
    r(118) =  kf_SD(9,2)*y(34+11)*XPt211*y(1)*y(1)                   -  kf_SD(9,2)/exp(-(sp(9)-sp(32+9))/(kB*T))*y(11)*XPt211*y(34+1)*y(34+1);          	    % CHC   
    r(119) =  kf_SD(10,2)*y(34+12)*XPt211*y(1)*y(1)                  -  kf_SD(10,2)/exp(-(sp(10)-sp(32+10))/(kB*T))*y(12)*XPt211*y(34+1)*y(34+1);          	% CC    
    r(120) =  kf_SD(11,2)*y(34+13)*XPt211*y(1)                       -  kf_SD(11,2)/exp(-(sp(11)-sp(32+11))/(kB*T))*y(13)*XPt211*y(34+1);                     % CH4   
    r(121) =  kf_SD(12,2)*y(34+14)*XPt211*y(1)                       -  kf_SD(12,2)/exp(-(sp(12)-sp(32+12))/(kB*T))*y(14)*XPt211*y(34+1);                     % CH3   
    r(122) =  kf_SD(13,2)*y(34+15)*XPt211*y(1)*y(1)                  -  kf_SD(13,2)/exp(-(sp(13)-sp(32+13))/(kB*T))*y(15)*XPt211*y(34+1)*y(34+1);          	% CH2   
    r(123) =  kf_SD(14,2)*y(34+16)*XPt211*y(1)*y(1)*y(1)*y(1)  	     -  kf_SD(14,2)/exp(-(sp(14)-sp(32+14))/(kB*T))*y(16)*XPt211*y(34+1)*y(34+1)*y(18);       % CH    
    r(124) =  kf_SD(15,2)*y(34+17)*XPt211*y(1)*y(1)*y(1)*y(1)  	     -  kf_SD(15,2)/exp(-(sp(15)-sp(32+15))/(kB*T))*y(17)*XPt211*y(34+1)*y(34+1) *y(18);   	% C	  
    r(125) =  kf_SD(16,2)*y(34+2)*XPt211*y(1)                     	 -  kf_SD(16,2)/exp(-(sp(16)-sp(32+16))/(kB*T))*y(2)*XPt211*y(34+1);                      % H  

	%% Rate of change in adsorbate coverages on Pt(100) 
	F = zeros(57,1);
	
	F(2) = 2*r(5) +r(7) +r(9) + r(10) + r(13) +r(12) + r(15) +r(17) + r(19) +r(20) +r(22) ...
        + r(24) +r(26) + r(28) +r(29) +r(30) +r(31) + r(125)/XPt100;  		            % H
	F(3) = r(1) - r(6) -r(7) + r(110)/XPt100; 								            % CH3CH3
	F(4) = r(7) -r(8) -r(9) -r(10) + r(111)/XPt100;						                % CH3CH2
	F(5) = r(9) -r(11) -r(12) -r(13) + r(112)/XPt100; 						            % CH3CH
	F(6) = r(12) - r(14)- r(15)  + r(113)/XPt100; 							            % CH3C
	F(7) = r(2) + r(10) - r(16) -r(17) + r(114)/XPt100; 					            % CH2CH2
	F(8) = r(13) +r(17) - r(18) - r(19) -r(20) + r(115)/XPt100; 			            % CH2CH
	F(9) = r(15) + r(19) - r(21) -r(22) + r(116)/XPt100; 					            % CH2C
	F(10) = r(3) + r(20) - r(23) - r(24) + r(117)/XPt100; 					            % CHCH
	F(11) = r(22) + r(24) -r(25) -r(26)  + r(118)/XPt100; 					            % CHC
	F(12) = r(26) - r(27) + r(119)/XPt100; 									            % CC
	F(13) = r(4) - r(28) + r(120)/XPt100; 									            % CH4
	F(14) = 2*r(6) + r(8) + r(11) +r(14) +r(28) -r(29) + r(121)/XPt100; 	            % CH3
	F(15) = r(8) +2*r(16) + r(18) + r(21) +r(29) - r(30) + r(122)/XPt100; 	            % CH2
	F(16) = r(11) + r(18) + 2*r(23) + r(25) + r(30) - r(31) + r(123)/XPt100;            % CH
	F(17) = r(14) + r(21) + r(25) + 2*r(27) + r(31) + r(124)/XPt100; 		            % C	
	
	%% Rate of change in adsorbate coverages on Pt(111) 
	F(17+2) = 2*r(31+5) +r(31+7) +r(31+9) + r(31+10) + r(31+13) +r(31+12) + r(31+15) +r(31+17) ...
        + r(31+19) +r(31+20) +r(31+22) + r(31+24) +r(31+26) + r(31+28) +r(31+29) +r(31+30) +r(31+31) + r(109)/XPt111; 	   % H
	F(17+3) = r(31+1) - r(31+6) -r(31+7) + r(94)/XPt111; 								                                   % CH3CH3
	F(17+4) = r(31+7) -r(31+8) -r(31+9) -r(31+10)  + r(95)/XPt111;						          	                       % CH3CH2
	F(17+5) = r(31+9) -r(31+11) -r(31+12) -r(31+13) + r(96)/XPt111; 						                               % CH3CH
	F(17+6) = r(31+12) - r(31+14)- r(31+15) + r(97)/XPt111; 							                                   % CH3C
	F(17+7) = r(31+2) + r(31+10) - r(31+16) -r(31+17) + r(98)/XPt111; 					                                   % CH2CH2
	F(17+8) = r(31+13) +r(31+17) - r(31+18) - r(31+19) -r(31+20) + r(99)/XPt111; 			                               % CH2CH
	F(17+9) = r(31+15) + r(31+19) - r(31+21) -r(31+22) + r(100)/XPt111; 					                               % CH2C
	F(17+10) = r(31+3) + r(31+20) - r(31+23) - r(31+24) + r(101)/XPt111; 					                               % CHCH
	F(17+11) = r(31+22) + r(31+24) -r(31+25) -r(31+26) + r(102)/XPt111; 					                               % CHC
	F(17+12) = r(31+26) - r(31+27) + r(103)/XPt111; 									                                   % CC
	F(17+13) = r(31+4) - r(31+28) + r(104)/XPt111; 									                                       % CH4
	F(17+14) = 2*r(31+6) + r(31+8) + r(31+11) +r(31+14) +r(31+28) -r(31+29)+ r(105)/XPt111; 	                           % CH3
	F(17+15) = r(31+8) +2*r(31+16) + r(31+18) + r(31+21) +r(31+29) - r(31+30) + r(106)/XPt111; 	                           % CH2
	F(17+16) = r(31+11) + r(31+18) + 2*r(31+23) + r(31+25) + r(31+30) - r(31+31) + r(107)/XPt111;                          % CH
	F(17+17) = r(31+14) + r(31+21) + r(31+25) + 2*r(31+27) + r(31+31) + r(108)/XPt111; 		                               % C
	
	%% Rate of change in adsorbate coverages on Pt(211) 
	F(34+2) = 2*r(62+5) +r(62+7) +r(62+9) + r(62+10) + r(62+13) +r(62+12) + r(62+15) +r(62+17) ...
        + r(62+19) +r(62+20) +r(62+22) + r(62+24) +r(62+26) + r(62+28) +r(62+29) +r(62+30) +r(62+31) -r(109)/XPt211 -r(125)/XPt211; 	% H
	F(34+3) = r(62+1) - r(62+6) -r(62+7) -r(94)/XPt211 -r(110)/XPt211; 								                                    % CH3CH3
	F(34+4) = r(62+7) -r(62+8) -r(62+9) -r(62+10) -r(95)/XPt211 -r(111)/XPt211;						                                    % CH3CH2
	F(34+5) = r(62+9) -r(62+11) -r(62+12) -r(62+13) -r(96)/XPt211 -r(112)/XPt211; 						                                % CH3CH
	F(34+6) = r(62+12) - r(62+14)- r(62+15) -r(97)/XPt211 -r(113)/XPt211; 							                                    % CH3C
	F(34+7) = r(62+2) + r(62+10) - r(62+16) -r(62+17) -r(98)/XPt211 -r(114)/XPt211; 					                                % CH2CH2
	F(34+8) = r(62+13) +r(62+17) - r(62+18) - r(62+19) -r(62+20) -r(99)/XPt211 -r(115)/XPt211; 			                                % CH2CH
	F(34+9) = r(62+15) + r(62+19) - r(62+21) -r(62+22) -r(100)/XPt211 -r(116)/XPt211; 					                                % CH2C
	F(34+10) = r(62+3) + r(62+20) - r(62+23) - r(62+24) -r(101)/XPt211 -r(117)/XPt211; 					                                % CHCH
	F(34+11) = r(62+22) + r(62+24) -r(62+25) -r(62+26) -r(102)/XPt211 -r(118)/XPt211; 					                                % CHC
	F(34+12) = r(62+26) - r(62+27) -r(103)/XPt211 -r(119)/XPt211; 									                                    % CC
	F(34+13) = r(62+4) - r(62+28) -r(104)/XPt211 -r(120)/XPt211; 									                                    % CH4
	F(34+14) = 2*r(62+6) + r(62+8) + r(62+11) +r(62+14) +r(62+28) -r(62+29) -r(105)/XPt211 -r(121)/XPt211; 	                            % CH3
	F(34+15) = r(62+8) +2*r(62+16) + r(62+18) + r(62+21) +r(62+29) - r(62+30) -r(106)/XPt211 -r(122)/XPt211; 	                        % CH2
	F(34+16) = r(62+11) + r(62+18) + 2*r(62+23) + r(62+25) + r(62+30) - r(62+31) -r(107)/XPt211 -r(123)/XPt211;                         % CH
	F(34+17) = r(62+14) + r(62+21) + r(62+25) + 2*r(62+27) + r(62+31) -r(108)/XPt211 -r(124)/XPt211; 		                            % C	
	
	%% Rate of change in Free Site Coverage on Pt(100)	
  	F(1) = -r(1) - 2*r(2) - 4*r(3) - r(4) - 2*r(5) - r(6) - r(7) - 2*r(8) - 2*r(9) - 2*r(10) - 3*r(11) - 3*r(12) - 2*r(13) - r(14) - 2*r(16) - 2*r(17) - 3*r(18) - r(19) ...
 		- 2*r(20) - 3*r(21) - 2*r(22) - 4*r(23) - r(24) - 4*r(25) - r(26) - 4*r(27) - r(28) - 2*r(29) - 3*r(30) - r(31) - 2*x100*r(62+3) - 2*x100*r(62+20) - 2*x100*r(62+22) ...
		+ 2*x100*r(62+23) + 2*x100*r(62+25) + 2*x100*r(62+27) + (-r(110) -r(111)  -2*r(112) -4*r(113) - 2*r(114) -3*r(115) -3*r(116) -2*r(117) -2*r(118) -2*r(119) -r(120) ...
		-r(121) -2*r(122) -4*r(123) -4*r(124) -r(125) + 2*r(101) + 2*r(102) + 2*r(103))/XPt100;
 	
    %% Rate of change in Free Site Coverage on Pt(111)
  	F(17+1) = -r(31+1) - 2*r(31+2) - 3*r(31+3) - r(31+4) - 2*r(31+5) - r(31+6) - r(31+7) - 2*r(31+8) - 2*r(31+9) - 2*r(31+10) - 2*r(31+11) - 2*r(31+12) - 2*r(31+13) - r(31+14) ...
 		- r(31+15) - 2*r(31+16) - 2*r(31+17) - 2*r(31+18) - r(31+19) - r(31+20) - 2*r(31+21) - r(31+22) - 3*r(31+23) - r(31+24) - 3*r(31+25) - 2*r(31+26) - 2*r(31+27) - r(31+28) ...
 		- 2*r(31+29) - 2*r(31+30) - r(31+31) - x111*r(62+11) - x111*r(62+12) - x111*r(62+13) - x111*r(62+15) - x111*r(62+17) - x111*r(62+19) - 2*x111*r(62+23) - 2*x111*r(62+25) ...
		- 2*x111*r(62+27) - x111*r(62+30) + x111*r(62+20) + x111*r(62+21) + 2*x111*r(62+22) + (-r(94) -r(95) -2*r(96) -2*r(97)  -2*r(98) -2*r(99) -r(100) -3*r(101) -3*r(102) -4*r(103) ...
		-r(104) -r(105) -2*r(106) - 2*r(107) -2*r(108) -r(109)  +r(113) +r(115) + 2*r(116) +r(123) +r(124))/XPt111;
 
    %% Rate of change in Free Site Coverage on Pt(211)
  	F(34+1) = -r(62+1) - 2*r(62+2) - 2*r(62+3) - r(62+4) - 2*r(62+5) - r(62+6) - r(62+7) - 2*r(62+8) - 2*r(62+9) - 2*r(62+10) - r(62+11) - r(62+12) - r(62+13) - r(62+14) - 2*r(62+16) ...
 		- r(62+17) - 2*r(62+18) - 0*r(62+19) - r(62+20) - 3*r(62+21) - 2*r(62+22) - 2*r(62+23) - r(62+24) - 2*r(62+25) - r(62+26) - 2*r(62+27) - r(62+28) - 2*r(62+29) - r(62+30) - r(62+31) ...
 		 + (r(94)  +r(95)  +2*r(96)  +2*r(97)  +2*r(98)  +2*r(99)  +r(100)  +2*r(101)  +2*r(102)  +2*r(103)  + r(104)  + r(105) +2*r(106)  +2*r(107)  +2*r(108)  +r(109) ...
 		 +r(110) +r(111) +2*r(112) +2*r(113) +2*r(114) +2*r(115) +r(116)  +2*r(117)  +2*r(118)  +2*r(119)  +r(120)  +r(121)  +2*r(122)  +2*r(123)  +2*r(124) + r(125))/XPt211;
	
	%% Rate of change in gas mole fraction in CSTR
	R_gas_Pt100 = -(r(1) + r(2) + r(3) + r(4) + r(5));
	R_gas_Pt111 = -(r(31+1) + r(31+2) + r(31+3) + r(31+4) + r(31+5));
	R_gas_Pt211 = -(r(62+1) + r(62+2) + r(62+3) + r(62+4) + r(62+5));

	F(52) = PP_gas(1)/(tauS*P_tot) -y(52)*(1/tauS +(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211)*numS) -numS*(XPt100*r(1) + XPt111*r(31+1) + XPt211*r(62+1));
	F(53) = PP_gas(2)/(tauS*P_tot) -y(53)*(1/tauS +(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211)*numS) -numS*(XPt100*r(2) + XPt111*r(31+2) + XPt211*r(62+2));
	F(54) = PP_gas(3)/(tauS*P_tot) -y(54)*(1/tauS +(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211)*numS) -numS*(XPt100*r(3) + XPt111*r(31+3) + XPt211*r(62+3));
	F(55) = PP_gas(4)/(tauS*P_tot) -y(55)*(1/tauS +(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211)*numS) -numS*(XPt100*r(4) + XPt111*r(31+4) + XPt211*r(62+4));	
    F(56) = PP_gas(5)/(tauS*P_tot) -y(56)*(1/tauS +(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211)*numS) -numS*(XPt100*r(5) + XPt111*r(31+5) + XPt211*r(62+5));
    F(57) = PP_gas(6)/(tauS*P_tot) -y(57)*(1/tauS +(R_gas_Pt100*XPt100 + R_gas_Pt111*XPt111 + R_gas_Pt211*XPt211)*numS) ;

end	
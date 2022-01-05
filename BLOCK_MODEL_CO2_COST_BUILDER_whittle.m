%BLOCK MODEL CO2 COST BUILDER
%this matlab script reads the block model and computes the MCAF, PCAF and
%rehabilitation cost needed to include the carbon tax/price into the
%whittle NPV optimisation


clc; clear all; close all

T = readtable("marvin_reblocked.csv"); %import dataset (BM)

%% setup parameters

%parameters to compute energy and GWP
A = 116.9; %area of drilling section,cm2
L = 35; %drilling length,cm 
N = 10; %number of drillings per block
nD = 80; %drillers efficiency,%
LF = 6; %load factor(kg ANFO per ton of block)
Rs = 2; %rolling resistance of the truck,%
Ri = 1; %internal resistance of the truck,%
MT = 166; %mass of loaded truck
mT = 92.2; %filling capacity of the truck
rs = 10; %ramp slope,%
Nl = 70; %loader efficiency,%
t = 45; %loader cycle time,s
PI = 0.18; %loader power,mW
E_ANFO = 3.81581; %specific energy of the used explosive
Ev = 148*10^-6; %drilling specific energy of R1,MJ/cm3     TO BE DEFINED BASED ON THE ROCK TYPE
d1 = 0.09159134;%background energy value of diesel  kg co2 eq/ MJ
d2 = 0.51997136;%background energy value of diesel kg co2 eq/ MJ
b1 = 3.7;%background energy value of ANFO kg co2 eq/ MJ
b2 = 8.4;%background energy value of ANFO kg co2 eq/ MJ
wi = 15; %bond work index, kwh                             TO BE DEFINED BASED ON THE ROCK TYPE
crush_size_in = 75000; %paricle size, um
crush_size_out = 15000; %particle size, um
mill_size_in = 15000; %particle size, um
mill_size_out = 100; %particle size, um
p1 =0.56458524; %electricity impact, kg co2 eq/ KwH

%plant position                                            TO BE DEFINED 
X_p=74.452;%x-coord, m
Y_p=1101.133;%y-coord, m
Z_p=510;%z-coord, m

%dump position                                             TO BE DEFINED
X_d=35.665;%x-coord, m
Y_d=590.351;%y-coord, m
Z_d=510;%z-coord, m

%economic parameters
carbon_price = 24;%carbon price, $/tCO2-e
reference_mining_cost = 2; %reference mining cost, $/ton
reference_processing_cost = 15; %reference processing cost, $/ton
cost_increase_dept = 0.05; % increase per depth, $/ton/bench
start_cost_increase_Z = 442.5; %elevation where to start the cost increase with depth, m

%% computation of the specific energy [MJ/ton]

energy_drill = ((A*Ev*L*N)./((nD/100).*T.block_tonnage));
no_val = find(energy_drill == Inf);
energy_drill(no_val) = 0;
energy_blast = LF*E_ANFO;
energy_load = ((PI*t)/((Nl/100)*mT));

S_bp = sqrt((X_p-T.X).^2+(Y_p-T.Y).^2+(Z_p-T.Z).^2).*0.001;%distance block-plant [km]
S_bd = sqrt((X_d-T.X).^2+(Y_d-T.Y).^2+(Z_d-T.Z).^2).*0.001;%distance block-dump [km]

energy_haul_to_plant = ((9.81*S_bp*(mT*(rs/100)+((Rs+Ri)/100)*(2*MT-mT)))/mT);
energy_haul_to_dump = ((9.81*S_bd*(mT*(rs/100)+((Rs+Ri)/100)*(2*MT-mT)))/mT);

energy_crush = 3.6*wi*10*((1/sqrt(crush_size_out)-(1/sqrt(crush_size_in))));
energy_mill = 3.6*wi*10*((1/sqrt(mill_size_out)-(1/sqrt(mill_size_in))));

%% computations of GWP  [t CO2-e/ton]

GWP_drill = (d1+d2)*energy_drill/1000;
GWP_blast = (energy_blast/b1)*b2/1000; 
GWP_load = (d1+d2)*energy_load/1000;
GWP_haul_to_plant = (d1+d2)*energy_haul_to_plant/1000;
GWP_haul_to_dump = (d1+d2)*energy_haul_to_dump/1000;
GWP_proc = ((energy_crush+energy_mill)/3.6)*p1/1000;

%% computation GWP cost [$/ton]

cost_drill = carbon_price*GWP_drill;
cost_blast = carbon_price*GWP_blast; 
cost_load = carbon_price*GWP_load;
cost_haul_to_plant = carbon_price*GWP_haul_to_plant;
cost_haul_to_dump = carbon_price*GWP_haul_to_dump;
cost_proc = carbon_price*GWP_proc;

%% mining cost varying with dept

bottom_BM = min(T.Z); %find the bottom z coord of the bm
mining_cost_bottom = reference_mining_cost + ((start_cost_increase_Z - bottom_BM)/unique(T.size_Z_))*cost_increase_dept; %compute the mining cost at the bottom

%linear fit from top to bottom (the mining cost increases linearly)
[xData, yData] = prepareCurveData( [start_cost_increase_Z bottom_BM], [reference_mining_cost mining_cost_bottom] );
ft = fittype( 'poly1' );
[fitresult, gof] = fit( xData, yData, ft );
coeffvals= coeffvalues(fitresult);

%find mining cost for all coord Z
cost_increase_dept = coeffvals(1)*T.Z + coeffvals(2);

%where the Z coord >= start_cost_increase_Z set the cost to the reference
%mining cost
index_coord_costant_cost = find(T.Z >= start_cost_increase_Z);
cost_increase_dept(index_coord_costant_cost) = reference_mining_cost;

%% computation of MCAF, PCAF and rehab cost

CO2_cost_mining = cost_drill + cost_blast + cost_load;
total_mining_cost = cost_increase_dept + CO2_cost_mining;

CO2_cost_processing = cost_haul_to_plant + cost_proc;
total_processing_cost = reference_processing_cost + CO2_cost_processing;

MCAF = total_mining_cost/reference_mining_cost;
PCAF = total_processing_cost/reference_processing_cost;

max_MCAF = max(MCAF)
min_MCAF = min(MCAF)

max_PCAF = max(PCAF)
min_PCAF = min(PCAF)

rehab = cost_haul_to_dump; 

%IMPORTANT NOTE
%hauling to dump:
%this cost can be included within Whittle in two ways: 
%a) use the expression builder to define the entire rehabilitation
%expression
%b) import the variable 'rehab' as grade and point to this grade variable
%with the expression builder. If this option is chosen, the 'rehab' will
%appear in Whittle as a new ore, therefore the price of the 'rehab' must be
%set to 0 $

T = addvars(T,MCAF,PCAF,rehab);%add variable created to block model
writetable(T,"marvin_reblocked_CO2_COSTS.csv")




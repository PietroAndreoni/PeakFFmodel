$onText
This is a simple model to analize the role of oil, gas and coal as substitues
in various population and mitigation scenarios.
Mitigation (expressed in Carbon intensity trajectories from IPCC scenarios, and is
considered a proxy of both efficiency improvement and develompment and spread of
fossil fuels substitutes, both as primary energy and as final energy) is assumed
to be exogenous.
All emissions are assumed to be C02 and come from fossil fuels.
Fossil fuels are assumed to be perfect substitutes between themselves.
Ocean and land are assumed to sink half of C02 emissions at a constant rate.
No CCS is taken into consideration.
Unconventional natural resources such as shale are taken into consideration,
but not syntetic fuels or gasification.
The goal is to see if and what mix of fossil fuels can meet the residual fossil
demand (expressed in energy units) minimizing total emissions through different
emissions and socio-economic scenarios.
In particular, we want to verify how much fossil fuel (and of which type) must
be kept in the ground to meet various target and what influence the socio-economic
projections have in meeting certain goals.
Results in cumulative emissions are trasformed in atmospheric concentrations
in order to confront them with the policy targets.
Starting time is 2010.
RCP2.6 scenario in SSP3 is infeasible, therefore it is rapresented as another baseline
for coding reason
$offText

$onText
PROBLEMS E TODO
- set initial values for the system
- set lower and upper boundaries to stability
- create tables and data for multidimensional parameters
- increase robustness: what happens if FF cannot meet the demand?
$offText


set        t     Time periods (10 years per period)                    /1*10/;
set        SSP   Population scenario                              /SSP1*SSP5/;
set        RCP   Radiative forcing scenario                       /RCP26,BAU/;
set        i     Fossil fuels                                  /coal,oil,gas/;


scalar

         tstep   Years per period                                        /10/;

*multidimensional parameters
$include Datatables

parameters
         K(i)               Maximum recoverable resources by fossil source
                                                                         / coal
                                                                            oil
                                                                            gas  /
         EMFAC(i)           Emission factor by fossil source in
                                                                          / coal 18.1082399
                                                                            oil  35.7285436
                                                                            gas  25.79920907 /
         r(i)               Maximum growth rate of production by fossil source
                                                                          / coal 0.03
                                                                            oil  0.03
                                                                            gas  0.03 /
         tar(RCP)           Cumulative emission target by policy scenario /BAU
                                                                           RCP26/;

*initial values of variables
parameters

variables

         D(t)       Demand for fossil fuels energy
         P(t,i)     Production of fossil fuel i
         PTOT(t)    Total production of all fossil fuels (MWh)
         Np(t,i)    Cumulative production of fossil fuel i
         CIFF(t)    Carbon intensity of fossil fuel energy
         FFI(t)     Fossil fuel intensity of GDP in TWh over GDP
         EM(t)      Carbon emissions
         CUMEM      Cumulative emissions in 2100. Objective variable ;

equations

         DEQ             Demand equation
         PEQ             Production equation
         SPEQEQ          Supply and demand equilibrium equation
         CUMPEQ          Cumulative production equation
         PCOSTR          Geometrical progression constraint on avaliable resources
         EMEQ            Emission equations
         CIFFEQ          Carbon intensity of fossil fuels
         FFIEQ           Fossil fuel intensity equation
         CUMEMEQ         Cumulative emission equations ;


deq(t)..             D(t)       =E=   FFI(t) * GDP(t,SSP);
peq(t)..             PTOT(t)    =E=   sum(i,P(t,i));
speqeq(t)..          D(t)       =LE=  PTOT(t);
cumpeq(t,i)..        Np(t+1,i)  =E=   Np(t,i) + P(t,i);
pcostr(t,i)..        P(t,i)     =LE=  r(i) * Np(t,i) * ( 1 - Np(t,i)/K(i) );
emeq(t)              EM(t)      =E=   sum(i,P(t,i) * EMFAC(i));
ciffeq(t)..          CIFF(t)    =E=   EM(t) / PTOT(t);
ffieq(t)..           FFI(t)     =E=   CIFF(t) * CIGDP(t,SSP,RCP);
CUMEMEQ..            CUMEM      =E=   CUMEM(t0) + sum(t,EM(t));

*set initial values

*set boundaries

model peakoil /all/;

solve peakoil minimizing CUMEM using nlp;

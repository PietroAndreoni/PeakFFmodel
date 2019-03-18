$onText
This is a simple model to analize the role of oil, gas and coal as substitues
in various population and mitigation scenarios.
Mitigation (expressed in Carbon intensity trajectories from IPCC scenarios, and
considered a proxy of both efficiency improvement and develompment and spread of
fossil fuels substitutes) is assumed to be exogenous.
All emissions are assumed to be C02 and come from fossil fuels.
Fossil fuels are assumed to be perfect substitutes between themselves.
A certain maximum fuel switching rate is assumed. This is crearly highly arbitrary
and sector dependent, however as this model does not take into account costs, it
can be assumed that with greater economic effort substitution (between FF)
in all sector can be achieved.
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
- set initial values for the system  DONE
- set lower and upper boundaries to stability
- create tables and data for multidimensional parameters  DONE
- decide and add some parameters (K, r e maxfs)
- increase robustness: what happens if FF cannot meet the demand?
- sensitivity analysis: what happens if K(i) rises?
(- introducing price as a variable and endogenizing non fossil substitutes)
- Add a small climate model to translate cumulative emissions as concentrations
- Hubbert geometric model does apply to coal?
- make a confrontation between historical production and theoretical Hubbert maximum
production to see if different fossil fuels are or not demand driven
$offText


set        t     Time periods (10 years per period)                    /2010
                                                                        2020
                                                                        2030
                                                                        2040
                                                                        2050
                                                                        2060
                                                                        2070
                                                                        2080
                                                                        2090
                                                                        2100/;

set        SSP   Population scenario                              /SSP1*SSP5/;

set        RCP   Radiative forcing scenario                       /BAU,RCP26/;

set        i     Fossil fuels                                  /coal,oil,gas/;

set        tfirst(t);

scalar     tstep   Years per period                                        /10/;


tfirst(t) = yes$(t.val eq 2010);

*multidimensional parameters
$include Datatables

parameters
         K(i)               Maximum recoverable resources by fossil source
                                                                         / coal 10000000000000000000000000000
                                                                            oil 10000000000000000000000000000
                                                                            gas 10000000000000000000000000000 /
         EMFAC(i)           Emission factor by fossil source in MtC02 over TWh
                                                                          / coal 35.7285436
                                                                            oil  25.7992090
                                                                            gas  18.1082399  /
         r(i)               Maximum growth rate of production by fossil source
                                                                          / coal 0.03
                                                                            oil  0.03
                                                                            gas  0.03 / ;
*initial values for the variables

parameters
         P0(i)      Production at time t0 of fossil fuel i in relative terms
                                                                         / coal 0.34
                                                                           oil  0.39
                                                                           gas  0.27  /
         Np0(i)     Cumulative production at time t0 (from 1960)
                                                                         / coal 1171358.577
                                                                           oil  1652597.944
                                                                           gas  882942.0745 / ;

parameter
         D(t)       Demand of fossil fuels through scenarios;

scalars
         CUMEM0     Cumulative emission at period t0 (from 1960)         / 991306.0139 /
         maxfs      Maximum fuel switching velocity in a period          / 0.5 /;

*targets and parameters for confrontations
scalar
         tar        Cumulative emission target by policy scenario        /1000/;

*support parameters for scenario looping and results saving
parameter
         P_res(t,i,SSP,RCP), EM_res(t,SSP,RCP);

*random initialization of support parameters
P_res(t,i,SSP,RCP) = 0;
EM_res(t,SSP,RCP) = 0;


nonnegative variables

         Np(t,i)    Cumulative production of fossil fuel I
         EM(t)      Carbon emissions
         P(t,i)     Production of fossil fuel i (ratio) ;


variable
         CUMEM      Cumulative emissions in 2100. Objective variable ;

equations

         PEQ             Production equation
         PGRRATEEQ       Maximum growth rate costraint for production
         CUMPEQ          Cumulative production equation
         PCOSTR          Geometrical progression constraint on avaliable resources
         EMEQ            Emission equations
         CUMEMEQ         Cumulative emission equations ;


peq(t)$(D(t) gt 0)..                             sum(i,P(t,i)) =E=   1;
pgrrateeq(t+1,i)$(D(t+1) ge maxfs*D(t))..        P(t+1,i)      =G=   maxfs*P(t,i);
cumpeq(t+1,i)..                                  Np(t+1,i)     =E=   Np(t,i) + D(t)*P(t,i);
pcostr(t,i)..                                    D(t)*P(t,i)   =L=   ( r(i) * Np(t,i) * ( 1 - Np(t,i)/K(i) ) ) ;
emeq(t)..                                        EM(t)         =E=   D(t) * sum(i,P(t,i) * EMFAC(i));
CUMEMEQ..                                        CUMEM         =E=   CUMEM0 + sum(t,EM(t));


*set boundaries for stability
P.lo(t,i)=0;
P.up(t,i)=1;


model peakoil /all/;


*running model through scenarios
loop ( (SSP,RCP),

        D(t) = 0;
        D(t) = ( FFIGDP(t,SSP,RCP) * GDP(t,SSP,RCP) )$( FFIGDP(t,SSP,RCP) gt 0);

*set initial values
        P.fx(tfirst,i) = P0(i);
        Np.fx(tfirst,i) = Np0(i);

options nlp=conopt4; solve peakoil minimizing CUMEM using nlp;

         P_res(t,i,SSP,RCP) = P.l(t,i);
         EM_res(t,SSP,RCP) = EM.l(t);

      );

display P_res,EM_res;

execute_unload "results.gdx";
execute '=gdx2xls results.gdx';

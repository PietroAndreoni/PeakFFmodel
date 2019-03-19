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
- introduce a "dummy variable" modelled as a fourth FF with very high emission intensity,
whose production start to zero. This way, if with a very low reserve scenario FF
cannot meet the demand, model is still feasibile, even if results must be manipulated
to remain meaningful.
- IMPORTANT: correct parameter FFIGDP which is now actually CIGDP (even if they are good proxy)
- in the paper: search dependence of oil consumptio vs population and vs GDP to show
that the second is stronger
- GAMS GENERAL KNOWLEDGE: is referring the equation to (t+1) the same as (ord(t) le card(t)-1)?
- Think about first condition in pggrateeq
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

set        i     Fossil fuels                            /coal,oil,gas,dummy/;

set        tfirst(t);

tfirst(t) = yes$(t.val eq 2010);

*multidimensional parameters
$include Datatables

scalar     tstep   Years per period                                        /10/;

parameters
         K(i)               Maximum recoverable resources by fossil source
                                                                          / coal   6908219.983
                                                                            oil    2779569.993
                                                                            gas    2407409.994
                                                                            dummy  10000000000 /

         EMFAC(i)           Emission factor by fossil source in MtC02 over TWh
                                                                          / coal  35.7285436
                                                                            oil   25.7992090
                                                                            gas   18.1082399
                                                                            dummy 50  /

         r(i)               Maximum growth rate of production by fossil source
                                                                          / coal  0.052
                                                                            oil   0.051
                                                                            gas   0.06
                                                                            dummy 1 /

         maxfs(i)           Maximum fuel switching velocity in a period
                                                                          / coal  0.5
                                                                            oil   0.5
                                                                            gas   0.5
                                                                            dummy 0 / ;
*initial values for the variables

parameters
         P0(i)      Production at time t0 of fossil fuel i in relative terms
                                                                         / coal  0.34
                                                                           oil   0.39
                                                                           gas   0.27
                                                                           dummy 0 /
         Np0(i)     Cumulative production at time t0 (from 1960)
                                                                         / coal  1171358.577
                                                                           oil   1652597.944
                                                                           gas   882942.0745
                                                                           dummy 10000000 / ;

parameter
         D(t)       Demand of fossil fuels through scenarios;

scalars
         CUMEM0     Cumulative emission at period t0 (from 1960) relativo to total demand at to
                                                                        / 8.044785625 /;

*targets and parameters for confrontations
scalar
         tar        Cumulative emission target by policy scenario        /1000/;

*support parameters for scenario looping and results saving
parameter
         P_res(t,i,SSP,RCP), EM_res(t,SSP,RCP),D_res(t,SSP,RCP);

*initialization of support parameters
P_res(t,i,SSP,RCP) = 0;
EM_res(t,SSP,RCP) = 0;
D_res(t,SSP,RCP) = 0;


nonnegative variables

         Np(t,i)    Cumulative production of fossil fuel I
         EM(t)      Carbon emissions
         P(t,i)     Production of fossil fuel i (ratio);


variable
         CUMEM      Cumulative emissions in 2100. Objective variable ;

equations

         PEQ             Production equation
         PGRRATEEQ       Maximum growth rate costraint for production
         CUMPEQ          Cumulative production equation
         PCOSTR          Geometrical progression constraint on avaliable resources
         EMEQ            Emission equations
         CUMEMEQ         Cumulative emission equations ;


peq(t)$(D(t) gt 0)..                                             sum(i,P(t,i)) =E=   1;
pgrrateeq(t+1,i)$(D(t) gt 0)..                                   P(t+1,i)      =G=   maxfs(i)*P(t,i);
cumpeq(t+1,i)..                                                  Np(t+1,i)     =E=   Np(t,i) + D(t)*P(t,i);
pcostr(t,i)$(ord(t) gt 1)..                                      D(t)*P(t,i)   =L=   r(i) * Np(t,i) * ( 1 - Np(t,i)/K(i) );
emeq(t)..                                                        EM(t)         =E=   sum(i,P(t,i) * EMFAC(i));
CUMEMEQ..                                                        CUMEM         =E=   CUMEM0 + sum(t,EM(t));


*set boundaries for stability
P.lo(t,i)=0;
P.up(t,i)=1;
Np.lo(t,i)=Np0(i);
Np.up(t,i)=K(i);


model peakoil /all/;


*running model through scenarios
loop ( (SSP,RCP),

         D(t) = 0;
         D(t) = (3.68*CIGDP(t,SSP,RCP) * GDP(t,SSP,RCP)) $ (CIGDP(t,SSP,RCP) gt 0);
*3.68 is a scaling factor based on t0 data that gives the relationship between carbon intensity and FFIofGDP
*the trend is assumed to be the same (hypothesis supported by past data)

*set initial values
         P.fx(tfirst,i) = P0(i);
         Np.fx(tfirst,i) = Np0(i);
*"suggest" production to zero if demand if zero to avoid "free production" at zero demand
*         P.l(t,i)$(D(t) eq 0) =0;
*Force dummy at zero in first periods
*         P.fx(t,"dummy")$(ord(t) le 5) = 0;


options nlp=conopt4; solve peakoil minimizing CUMEM using nlp;

         P_res(t,i,SSP,RCP) = P.l(t,i);
         EM_res(t,SSP,RCP) = EM.l(t);
         D_res(t,SSP,RCP) = D(t);

      );

*post processing for correct results if dummy is actively used

loop ( (t,SSP,RCP),

         if ( P_res(t,"dummy",SSP,RCP) gt 0,

                 EM_res(t,SSP,RCP) = EM_res(t,SSP,RCP) - EMFAC("dummy")*P_res(t,"dummy",SSP,RCP);

             );

      );

execute_unload "results.gdx";
execute_unload "resultsshort.gdx",P_res,EM_res,D_res;
execute '=gdx2xls resultsshort.gdx';

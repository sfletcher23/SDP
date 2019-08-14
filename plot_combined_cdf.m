if false
fig3 = figure;
ax = subplot(4,1,1);
copyobj(fig2.Children(2).Children,ax);
xlim([0.5 8])
xlabel({})
ylabel({'Cumulative', 'density'})
legend('Build', 'No build', 'Flexible')
legend('boxoff')
ax = subplot(4,1,2);
copyobj(fig1.Children(1).Children,ax);
xlim([0.5 8])
ylabel({'Cumulative', 'density'})
xlabel({})
ax = subplot(4,1,3);
copyobj(fig1.Children(2).Children,ax);
xlim([0.5 8])
ylabel({'Cumulative', 'density'})
xlabel({})
ax = subplot(4,1,4);
copyobj(fig1.Children(4).Children,ax);
xlim([0.5 8])
ylabel({'Cumulative', 'density'})
xlabel('30-year Cost with Damages [Bn$]')


end

%%
if true   
clear all; close all
load('/Users/sarahfletcher/Documents/MATLAB/Repository_SDP/08_Jan_2018_14_10_12_sim_discount_rate.mat')
runParam.simNum = 10000;

fig = figure;

% Base: 0% DR, 20$ shortage

costParam.shortage_cost = 20; 
costParam.discount_rate = discount_rate_Output{1,1}{1};
V = discount_rate_Output{1,1}{2};
X1 = discount_rate_Output{1,1}{3};
X2 = discount_rate_Output{1,1}{4};
T_gw_all = discount_rate_Output{1,1}{5};
lowestCostAction = discount_rate_Output{1,1}{5};


useNoInfoPolicy = false;
[ sim ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors );

useNoInfoPolicy = true;
lowestCostAction = 2;
[ simnolearn_build ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors);
lowestCostAction = 0;
[ simnolearn_nobuild ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors);


totalCostFlex = sum(sim.costOverTime,2);
totalCostBuild = sum(simnolearn_build.costOverTime,2);
totalCostNoBuild = sum(simnolearn_nobuild.costOverTime,2);
totalShortageFlex = sum(sim.shortageOverTime,2);
totalShortageBuild = sum(simnolearn_build.shortageOverTime,2);
totalShortageNoBuild = sum(simnolearn_nobuild.shortageOverTime,2);


subplot(4,1,1)
hold on
c = cdfplot(totalCostBuild/1E9);
c.LineWidth = 1.5;
c = cdfplot(totalCostNoBuild/1E9)
c.LineWidth = 1.5;
c = cdfplot(totalCostFlex/1E9)
c.LineWidth = 1.5;
xlim([0.5 12])
xlabel({})
title('Base: 0% DR, $20 Shortage')
legend('Build', 'No build', 'Flexible')

% DR 5%, shortage base 20$
costParam.shortage_cost = 20; 
costParam.discount_rate = discount_rate_Output{3,1}{1};
V = discount_rate_Output{3,1}{2};
X1 = discount_rate_Output{3,1}{3};
X2 = discount_rate_Output{3,1}{4};
T_gw_all = discount_rate_Output{3,1}{5};
lowestCostAction = discount_rate_Output{3,1}{5};


useNoInfoPolicy = false;
[ sim ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors );

useNoInfoPolicy = true;
lowestCostAction = 2;
[ simnolearn_build ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors);
lowestCostAction = 0;
[ simnolearn_nobuild ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors);


totalCostFlex = sum(sim.costOverTime,2);
totalCostBuild = sum(simnolearn_build.costOverTime,2);
totalCostNoBuild = sum(simnolearn_nobuild.costOverTime,2);
totalShortageFlex = sum(sim.shortageOverTime,2);
totalShortageBuild = sum(simnolearn_build.shortageOverTime,2);
totalShortageNoBuild = sum(simnolearn_nobuild.shortageOverTime,2);


subplot(4,1,2)
hold on
c = cdfplot(totalCostBuild/1E9);
c.LineWidth = 1.5;
c = cdfplot(totalCostNoBuild/1E9)
c.LineWidth = 1.5;
c = cdfplot(totalCostFlex/1E9)
c.LineWidth = 1.5;
xlim([0.5 12])
xlabel({})
title('5% Discount Rate')




% DR 0%, shortage base 50$
costParam.discount_rate = 0;
costParam.shortage_cost = shortage_cost_Output{4,1}{1};
V = shortage_cost_Output{4,1}{2};
X1 = shortage_cost_Output{4,1}{3};
X2 = shortage_cost_Output{4,1}{4};
T_gw_all = shortage_cost_Output{4,1}{5};
lowestCostAction = shortage_cost_Output{4,1}{5};


useNoInfoPolicy = false;
[ sim ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors );

useNoInfoPolicy = true;
lowestCostAction = 2;
[ simnolearn_build ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors);
lowestCostAction = 0;
[ simnolearn_nobuild ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors);


totalCostFlex = sum(sim.costOverTime,2);
totalCostBuild = sum(simnolearn_build.costOverTime,2);
totalCostNoBuild = sum(simnolearn_nobuild.costOverTime,2);
totalShortageFlex = sum(sim.shortageOverTime,2);
totalShortageBuild = sum(simnolearn_build.shortageOverTime,2);
totalShortageNoBuild = sum(simnolearn_nobuild.shortageOverTime,2);


subplot(4,1,3)
hold on
c = cdfplot(totalCostBuild/1E9);
c.LineWidth = 1.5;
c = cdfplot(totalCostNoBuild/1E9)
c.LineWidth = 1.5;
c = cdfplot(totalCostFlex/1E9)
c.LineWidth = 1.5;
xlim([0.5 12])
xlabel({})
title('50$ Shortage Cost')



% DR 0%, shortage base 5$
costParam.discount_rate = 0;
costParam.shortage_cost = shortage_cost_Output{1,1}{1};
V = shortage_cost_Output{1,1}{2};
X1 = shortage_cost_Output{1,1}{3};
X2 = shortage_cost_Output{1,1}{4};
T_gw_all = shortage_cost_Output{1,1}{5};
lowestCostAction = shortage_cost_Output{1,1}{5};


useNoInfoPolicy = false;
[ sim ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors );

useNoInfoPolicy = true;
lowestCostAction = 2;
[ simnolearn_build ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors);
lowestCostAction = 0;
[ simnolearn_nobuild ] = sim_sdp_gw( X1, X2, V, T_gw_all, cumTgw, useNoInfoPolicy, lowestCostAction, runParam, gwParam, costParam, water, s_gw, s_expand, exp_vectors);


totalCostFlex = sum(sim.costOverTime,2);
totalCostBuild = sum(simnolearn_build.costOverTime,2);
totalCostNoBuild = sum(simnolearn_nobuild.costOverTime,2);
totalShortageFlex = sum(sim.shortageOverTime,2);
totalShortageBuild = sum(simnolearn_build.shortageOverTime,2);
totalShortageNoBuild = sum(simnolearn_nobuild.shortageOverTime,2);


subplot(4,1,4)
hold on
c = cdfplot(totalCostBuild/1E9);
c.LineWidth = 1.5;
c = cdfplot(totalCostNoBuild/1E9)
c.LineWidth = 1.5;
c = cdfplot(totalCostFlex/1E9)
c.LineWidth = 1.5;
xlabel('30-year Cost with Damages [Bn$]')
xlim([0.5 12])
title('5$ Shortage Cost')





end

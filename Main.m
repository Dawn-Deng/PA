%% Main.m
% 模块化版本：Main / Channel_model / Initialization / AO_model / AO_W / Signal_model / Problem_formulation

clear; clc;

[params, state, channelInfo] = Channel_model();
[state, initInfo] = Initialization(state, params);
[state, aoInfo] = AO_model(state, params);
[state, signalInfo] = Signal_model(state, params);
problemInfo = Problem_formulation(state, params, aoInfo);

fprintf('===== Channel model =====\n');
disp(channelInfo);
disp('Users q_k:');
disp(state.users);
disp('Current X:');
disp(state.X);
disp('Current theta:');
disp(state.theta);
disp('Current phi:');
disp(state.phi);
disp('Current composite channel H:');
disp(state.H);

fprintf('===== Initialization =====\n');
disp('G_pot:');
disp(initInfo.Gpot);
disp('Emax:');
disp(initInfo.Emax);
disp('Candidate pool C:');
disp(initInfo.candidatePool);
disp('Reference y_ref:');
disp(initInfo.referencePositions);
disp('Utility matrix U:');
disp(initInfo.utilityMatrix);
disp('Hungarian selected pairs:');
disp(struct2table(initInfo.matching.selectedPairs));
disp('S^(0):');
disp(initInfo.matching.serviceSet);
disp('X^(0):');
disp(initInfo.X0);
disp('theta^(0):');
disp(initInfo.theta0);
disp('phi^(0):');
disp(initInfo.phi0);

fprintf('===== AO =====\n');
disp('AO sum-rate history:');
disp(aoInfo.sumRateHistory);
disp(['Converged = ', num2str(aoInfo.converged), ', iterations = ', num2str(aoInfo.iterations)]);

fprintf('===== Signal model =====\n');
disp('Current S^(t):');
disp(signalInfo.S);
disp('Current W^(t):');
disp(signalInfo.W);
disp('x_rad:');
disp(signalInfo.xRad);
disp('y:');
disp(signalInfo.y);
disp('SINR:');
disp(signalInfo.sinr);
disp('rate:');
disp(signalInfo.rate);
disp('sumRate:');
disp(signalInfo.sumRate);
disp(struct2table(signalInfo.userMetrics));

fprintf('===== Problem formulation =====\n');
disp(problemInfo);

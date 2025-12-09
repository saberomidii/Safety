clear
clc
close all


% Data for Double Integartor: 
number_of_states=[121,441,961,1681,2601,3721,5041,6561,8281,10201,12321,14641,17161,19881,22801,25921,29241,32761,36481,40401];

di_lp=[1.4200828075408936,1.3708758354187012,1.7657439708709717, ...
           2.560851812362671,3.79272198677063,4.5215301513671875,6.762346982955933,8.056127071380615, ...
           12.915220022201538,16.006544828414917,26.421756982803345,32.52098894119263,41.25309681892395, ...
           48.73053288459778,49.987147092819214,82.24253606796265,92.37376809120178,100.78133201599121,139.29631996154785,132.19421100616455];

di_vi_0=[0.13049101829528809,0.6296248435974121,2.411344051361084, ...
           6.298371076583862,14.967566013336182,12.637799978256226,26.65756106376648,48.93324112892151, ...
           80.97895503044128,136.02607202529907,201.995197057724,263.1986770629883,404.91900086402893, ...
           571.2804019451141,601.072557926178,724.6335070133209,943.6016979217529 ,1224.005635023117 ,1566.9820461273193,1963.8433129787445];

di_vi_1=[0.014407873153686523,0.4863858222961426,2.0311529636383057, ...
           6.1510539054870605,14.476797819137573,12.00294017791748,26.376811027526855,46.0824179649353, ...
           83.37427020072937,136.29552292823792,192.25326085090637,261.0792520046234,403.28091502189636, ...
           551.970685005188,600.2877950668335,688.1570739746094,937.0828731060028,1220.7154140472412 ,1555.3492980003357,2031.1658220291138];

di_vi_2=[0.014101982116699219,0.4842829704284668,2.0293710231781006, ...
           6.160218000411987,14.420847177505493,12.032495975494385,24.916422128677368,46.018481969833374, ...
           83.21190404891968,134.52144289016724,192.71325087547302,260.9727439880371,400.28544998168945, ...
           551.7411789894104,598.8917269706726,686.6963930130005,898.3662049770355,1218.0796990394592,1561.7779350280762,2031.8947021961212];

di_vi_3=[0.014250040054321289,0.48580408096313477,2.03043794631958, ...
           6.170543909072876,14.019774913787842,12.021743059158325,25.046019077301025,46.028796911239624, ...
           79.90062499046326,134.53225898742676,192.28736305236816,260.17567682266235,401.8435158729553, ...
           551.6958830356598,598.5191400051117,686.3625199794769,895.9323570728302,1209.384094953537,1552.0347859859467,2022.2694141864777];


ip_lp=[1.2601981163024902,1.379805088043213,1.7107090950012207, ...
           2.43487286567688,3.2856709957122803,4.490180969238281,5.160366058349609,7.6395790576934814, ...
           9.955899000167847,12.446542978286743,17.322367906570435,21.2745521068573,25.695661067962646, ...
           34.408010959625244,37.44218301773071,54.45881199836731,55.437642097473145,72.29860186576843,94.65880680084229,104.63096594810486];

ip_vi_0=[0.16327619552612305,0.44762301445007324,1.4554669857025146, ...
           3.5769689083099365,7.332987070083618,13.352867841720581,24.332506895065308,41.27987813949585, ...
           64.71582412719727,108.08136796951294,155.5500409603119,218.10163187980652,306.92709398269653, ...
           406.233766078949,555.219484090805,712.5647871494293,918.8732528686523,1129.783457994461,1390.0585930347443,1679.174232006073];

ip_vi_1=[0.06304383277893066,0.32422900199890137,1.07407808303833, ...
           3.1184091567993164,6.862974166870117,13.19712495803833,23.776149034500122,41.23936700820923, ...
           65.0259120464325,108.88135004043579,154.4323980808258,218.6062400341034,308.7973871231079, ...
           409.7326719760895,541.9514539241791,689.2597270011902,888.0747950077057,1103.4443249702454,1365.2883820533752,1664.3920278549194];

ip_vi_2=[0.04392385482788086,0.3144409656524658,1.0645558834075928, ...
           3.1044809818267822,6.836055040359497,13.233392000198364,23.95664691925049,42.08556890487671, ...
           65.03753709793091,103.40324592590332,154.55566596984863,222.56648516654968,311.55576395988464, ...
           393.92195892333984,508.51102685928345,681.8342189788818,884.6576750278473,1100.9138991832733,1353.3363480567932,1644.572942018509];

ip_vi_3=[0.044557809829711914,0.30879783630371094,1.0653541088104248, ...
           2.872252941131592,6.3753721714019775,13.227785110473633,22.624710083007812,42.57619285583496, ...
           65.47860908508301,103.56996297836304,147.8182189464569,221.662761926651,297.5173850059509, ...
           390.6766149997711,504.04807114601135,649.4025690555573,841.6501250267029,1088.8664119243622,1346.4383630752563,1628.8848040103912];


% -------------------------------------------------------------------------
% --- Plotting Double Integrator (DI) ---
h_di = figure('Name', 'Performance Comparison of Planning Methods (Double Integrator)');
ax_di = axes('Parent', h_di);
hold(ax_di, 'on');

% Plot Data
plot(ax_di, number_of_states, di_lp, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'LP');
plot(ax_di, number_of_states, di_vi_0, '-x', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'VI ($\lambda=0$)');
plot(ax_di, number_of_states, di_vi_1, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'VI ($\lambda=0.01$)');
plot(ax_di, number_of_states, di_vi_2, '-d', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'VI ($\lambda=0.015$)');
plot(ax_di, number_of_states, di_vi_3, '-^', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'VI ($\lambda=0.02$)');
xlim([0,40401+2000])
% --- Formatting DI ---
xlabel(ax_di, 'Number of States', 'FontSize', 30, 'Interpreter', 'latex');
ylabel(ax_di, 'Time (seconds)', 'FontSize', 30, 'Interpreter', 'latex');
legend(ax_di, 'Location', 'northwest', 'FontSize', 25, 'Interpreter', 'latex');
grid(ax_di, 'on');
set(ax_di, 'LineWidth', 2, 'TickLabelInterpreter', 'latex', 'FontSize', 25); 

% *** START AXIS ADJUSTMENT DI ***
% Set Y-axis exponent for scientific notation
ax_di.YAxis.Exponent = 3; % Assuming max time is around 2000s (2x10^3) for DI
% *** END AXIS ADJUSTMENT DI ***

% --- Saving the Plot DI ---
output_filename_di = 'di_time_computation.pdf';
fprintf('Saving Double Integrator figure to %s...\n', output_filename_di);
set(ax_di, 'LooseInset', get(ax_di, 'TightInset'));
exportgraphics(h_di, output_filename_di, 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('Figure saved as %s\n', output_filename_di);
hold(ax_di, 'off');

% -------------------------------------------------------------------------
% --- Plotting Inverted Pendulum (IP) ---
h_ip = figure('Name', 'Performance Comparison of Planning Methods (Inverted Pendulum)');
ax_ip = axes('Parent', h_ip);
hold(ax_ip, 'on');

% Plot Data
plot(ax_ip, number_of_states, ip_lp, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'LP');
plot(ax_ip, number_of_states, ip_vi_0, '-x', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'VI ($\lambda=0$)');
plot(ax_ip, number_of_states, ip_vi_1, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'VI ($\lambda=0.1$)');
plot(ax_ip, number_of_states, ip_vi_2, '-d', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'VI ($\lambda=0.2$)');
plot(ax_ip, number_of_states, ip_vi_3, '-^', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'VI ($\lambda=0.3$)');
xlim([0,40401+2000])

% --- Formatting IP ---
xlabel(ax_ip, 'Number of States', 'FontSize', 30, 'Interpreter', 'latex');
ylabel(ax_ip, 'Time (seconds)', 'FontSize', 30, 'Interpreter', 'latex');
legend(ax_ip, 'Location', 'northwest', 'FontSize', 25, 'Interpreter', 'latex');
grid(ax_ip, 'on');
set(ax_ip, 'LineWidth', 2, 'TickLabelInterpreter', 'latex', 'FontSize', 25); 

% *** START AXIS ADJUSTMENT IP ***
% Set Y-axis exponent for scientific notation
ax_ip.YAxis.Exponent = 3; % Assuming max time is around 1700s (1.7x10^3) for IP
% *** END AXIS ADJUSTMENT IP ***

% --- Saving the Plot IP ---
output_filename_ip = 'ip_time_computation.pdf';
fprintf('Saving Inverted Pendulum figure to %s...\n', output_filename_ip);
set(ax_ip, 'LooseInset', get(ax_ip, 'TightInset'));
exportgraphics(h_ip, output_filename_ip, 'ContentType', 'vector', 'BackgroundColor', 'white');
fprintf('Figure saved as %s\n', output_filename_ip);
hold(ax_ip, 'off');
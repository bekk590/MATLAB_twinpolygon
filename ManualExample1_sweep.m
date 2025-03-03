clear
clc
close all

%% Binary tree parameters

%Significant figures should be well stated.
stress = 1.1e9;
h_mbr = 20e-9;
N = 4;
%l0 = (700e-6)*0.4;
N_sol = 30;

%w0 = (700e-9)*sqrt(2);


rl1 = 0.6;
%rl2 =  744/(700*0.4);
rl2 = 1;

values = 1.00:0.1:3.00;
%values = 1.17;

%Since this is before modification, rw1 and rw2 should be the same value
rw1 =  sqrt(2);
rw2=   1;

l0 = (700e-6);
w0 = (700e-9);

l_trans = 50e-6;
l_pad = 1e-9;
w_pad = 2.5e-6;

plot_flag = 1;
plot_op_flag = 1;


%(wc^2)/lc = ???
lc = 10e-6;
wc = 2.1e-6;

pad_trigger = -1;%0: 2_pads, 1:4_pads, -1: 1_pad

[Freqs, Q ,m_eff, S_F, eta, rl2_match, Q_match] = ...
         twin_polygon_sweep(stress, h_mbr, l0, w0, N, ...
                           N_sol,rl1, rl2, rw1, rw2, lc, wc,...
                           l_trans, l_pad, w_pad, values,...
                           plot_flag, plot_op_flag, pad_trigger);

%{
figure
set(gcf, 'color', 'w')
box on
plot(Freqs, Q, 'o')
xlabel('Frequency (Hz)')
ylabel('Q')
ax = gca;
ax.YScale = 'log';
%}

figure
set(gcf, 'color', 'w')
box on
plot(rl2_match, Q_match, 'o')
xlabel('rl2')
ylabel('Q')
ax = gca;
ax.YScale = 'log';

minQ = min(Q_match) * 0.95; 
maxQ = max(Q_match) * 1.05; 
ylim([minQ, maxQ]);    

minParam = min(rl2_match) * 0.995; 
maxParam = max(rl2_match) * 1.005; 
xlim([minParam, maxParam]);  

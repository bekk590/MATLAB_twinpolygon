clear
clc
close all

%% Binary tree parameters
stress = 1.1e9;
h_mbr = 20e-9;
N = 4;
l0 = 700e-6*0.4;
N_sol = 15;

w0 = 200e-9;

rl1 =  1/0.4;
rl2 = 1/0.4;
%rl2 =  744/(700*0.4);

values = 2.59:0.001:2.60;
%values = 2.59;

rw1 =  1/10;
rw2=   1/10;

wc = 2e-6;
lc = 700e-6;

l_trans = 5e-6;
l_pad = 20e-6;
w_pad = 4e-6;

plot_flag = 1;
plot_op_flag = 1;

lc = 350e-6*0.4;
wc = 600e-9;





pad_trigger = 0;%0: 2_pads, 1:4_pads


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

%}
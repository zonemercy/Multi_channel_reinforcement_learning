%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB M-file ZunPhyModel.m
% Authors : Vincent Gauthier and Mohamed-Haykel Zayani
% Emails : {vincent.gauthier, mohamed-haykel.zayani}@it-sudparis.eu
% Address : Laboratory CNRS S.A.M.O.V.A.R. - Dept RS2M
% Telecom Sud Paris * 9 rue C. Fourier * 91011 EVRY CEDEX * FRANCE
% Created : April 28th, 2010
% Updated : May 18th, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brief Description:
% ------------------
% ZunPhyModel function is used to determine the probability of successful
% receive of a packet in fonction of RSSI, Modulation and Coding process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References :
% ------------
% [1] G. Zhou, T. He, S. Krishnamurthy, and J.A. Stankovic, "Models and
% solutions for radio irregularity in wireless sensor networks,"
% ACM Transactions on Sensor Networks, vol. 2, 2006, pp. 221-262.
% [2] M. Zuniga and B. Krishnamachari, "Analyzing the transitional
% region in low power wireless links," Proc. of the Conference on Sensor
% and Ad Hoc Communications and Networks (SECON), 2004.
% [3] M.Z. Zamalloa and B. Krishnamachari, "An analysis of unreliability
% and asymmetry in low-power wireless links," ACM Transactions on Sensor 
% Networks, vol. 3, 2007, pp. 7.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p=ZunPhyModel_lam(sigma, PACKET_LENGTH, DATA_RATE, d1, pw, N0)

%clear all;
addpath('./');

%% Chipcom CC1000 radio and physical layer parameters
% Noise Figure 13dB deviation of 23dB from the standard model 
NOISE_FIGURE = 23;
% Bandwidth (hz)
BW = 30;
% Path Loss Exponent 
PATH_LOSS_EXPONENT = 4;
% Variance (Shadowing model)
SHADOWING_STANDARD_DEVIATION = sigma;
% Definition of the standard Distance
D0 = 1.0;
% Transmision Power
Prdbm = pw; %5
% Avg White Noise (dB)
NOISE = N0; %15
% Wavelength 
lambda = 12.5 * 10^(-2);
% Length in bytes
PREAMBLE_LENGTH = 40;
FRAME_LENGTH = PACKET_LENGTH-PREAMBLE_LENGTH;
% Node Minimum and Maximum Ranges
distmin=max(0,DATA_RATE-20);
distmax=min(19200,DATA_RATE+20);
% Noise floor dBm
NF = CalNoiseFloor(BW, NOISE_FIGURE);
%Covariance Matrix S
S=[6.0 -3.3; -3.3 3.7];

% Cholesky Decomposition is used to generate multivariate random
% variables:
%   covariance matrix COVM = T' x T
%   T is a 2x2 upper triangular 
%       P_T   = P_T + T(1,1) * rn1
%       P_N   = P_N + T(1,2) * rn1 + T(2,2) * rn2
%  where rn1 and rn2 are normal(0,1) random variables.
            
T11 = sqrt(S(1,1));
T12 = S(1,2)/sqrt(S(1,1));
T21 = 0;
T22 = sqrt( (S(1,1)*S(2,2)-S(1,2)^2) / S(1,1) );

noiseFloor = NF + quadl(@(rn1) Calc2(rn1,T11), -2, 2);
outputPower = Prdbm + dblquad(@(rn1,rn2)Calc1(rn1,rn2,T12,T22), -2, 2, -2, 2);


p_temp=dblquad(@(DATA_RATE, ShadowingStd) Calc(d1, distmin, distmax, noiseFloor,DATA_RATE,...
    BW, PATH_LOSS_EXPONENT, lambda, ShadowingStd, outputPower, NOISE, D0, PREAMBLE_LENGTH,FRAME_LENGTH,S),...
    distmin, distmax, -2*SHADOWING_STANDARD_DEVIATION, 2*SHADOWING_STANDARD_DEVIATION);
p=p_temp/((distmax-distmin)*4*SHADOWING_STANDARD_DEVIATION);

% p_temp=quadl(@(ShadowingStd) Calc(d1, distmin, distmax, noiseFloor,DATA_RATE,...
%     BW, PATH_LOSS_EXPONENT, lambda, ShadowingStd, outputPower, NOISE, D0, PREAMBLE_LENGTH,FRAME_LENGTH,S),...
%      -2*SHADOWING_STANDARD_DEVIATION, 2*SHADOWING_STANDARD_DEVIATION);
% p=p_temp/(4*SHADOWING_STANDARD_DEVIATION);


% Function that computes the probability of successful frame receive
function ppr = Calc(dist, distmin, distmax,noiseFloor,DATA_RATE, BW, PATH_LOSS_EXPONENT, lambda, ShadowingStd, outputPower, NOISE, D0, PREAMBLE_LENGTH,FRAME_LENGTH,S)

% PathLoss 
rssi = CalRSSI(dist, PATH_LOSS_EXPONENT, lambda, ShadowingStd, mean(outputPower), NOISE, D0, distmin, distmax);
% SNR
ebno = RSSIToEbNo(rssi,mean(noiseFloor),DATA_RATE,BW,1);
% Packet Error Rate (BPSK Modulation)
per = Perror(ebno, 5);
% Packet Error Rate according with coding gain (NRZ Encoding)
ppr = PacketCoding(per,1,PREAMBLE_LENGTH,FRAME_LENGTH);

function TX_DEVIATION = Calc1(rn1,rn2,T12,T22)
% Output Power Deviation 
TX_DEVIATION = T12 * rn1 + T22 * rn2;

function RX_DEVIATION = Calc2(rn1,T11)
%Noise Floor Deviation
RX_DEVIATION = T11 * rn1;
        


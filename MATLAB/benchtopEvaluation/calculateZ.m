function [Zmag, Zphase_deg, feq] = calculateZ(C, Rseries, tr, Rpar)
% impedance_RparC_seriesR
% Computes magnitude and phase of impedance for:
%   Rseries in series with (Rpar || C)
%
% Inputs:
%   C        - Capacitance (F)
%   Rseries  - Series resistance (Ohm)
%   tr       - Step rise time (s)
%   Rpar     - Parallel resistance (Ohm), default = 560e3
%
% Outputs:
%   Zmag        - |Z| at effective frequency (Ohm)
%   Zphase_deg - phase angle (degrees)
%   feq         - effective frequency from step rise time (Hz)
%
% Assumes current-step excitation.

    if nargin < 4 || isempty(Rpar)
        Rpar = 560e3; % default parallel resistor
    end

    % Effective frequency from step rise time
    feq = 0.35 / tr;
    omega = 2 * pi * feq;

    % Parallel impedance (Rpar || C)
    Zpar = 1 ./ (1./Rpar + 1j * omega * C);

    % Total impedance
    Ztot = Rseries + Zpar;

    % Outputs
    Zmag = abs(Ztot);
    Zphase_deg = angle(Ztot) * 180 / pi;
end
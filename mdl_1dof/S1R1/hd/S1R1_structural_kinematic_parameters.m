% Return Structural Kinematic Parameters of the Robot 
% S1R1
%
% Output:
% v_mdh [1x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [1x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [1x1]
%   Aktives Gelenk = 1, Passiv = 0
% NL [1x1]
%   Anzahl der Starrkörper (inklusive Basis)
% NKP [1x1]
%   Anzahl der Kinematikparameter im Vektor `pkin`
% NQJ [1x1]
%   Anzahl der Minimalkoordinaten der kinematischen Kette
% pkin_names (1x1) cell
%   Namen aller Kinematik-Parameter im Vektor `pkin`

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = S1R1_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [0;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 2;
% Anzahl der Kinematikparameter
NKP = 1;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 1;
% Namen der Kinematikparameter
pkin_names = {'d1'};

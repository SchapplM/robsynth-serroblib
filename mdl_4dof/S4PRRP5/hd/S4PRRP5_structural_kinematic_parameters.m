% Return Structural Kinematic Parameters of the Robot 
% S4PRRP5
%
% Output:
% v_mdh [4x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [4x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [4x1]
%   Aktives Gelenk = 1, Passiv = 0
% NL [1x1]
%   Anzahl der Starrkörper (inklusive Basis)
% NKP [1x1]
%   Anzahl der Kinematikparameter im Vektor `pkin`
% NQJ [1x1]
%   Anzahl der Minimalkoordinaten der kinematischen Kette
% pkin_names (1x6) cell
%   Namen aller Kinematik-Parameter im Vektor `pkin`

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = S4PRRP5_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 2; 3;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [1; 0; 0; 1;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 1; 1; 1;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 5;
% Anzahl der Kinematikparameter
NKP = 6;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 4;
% Namen der Kinematikparameter
pkin_names = {'a2', 'a3', 'a4', 'd2', 'd3', 'theta1'};

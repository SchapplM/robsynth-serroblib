% Return Structural Kinematic Parameters of the Robot 
% S6PRPRRR7
%
% Output:
% v_mdh [6x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [6x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [6x1]
%   Aktives Gelenk = 1, Passiv = 0
% NL [1x1]
%   Anzahl der Starrkörper (inklusive Basis)
% NKP [1x1]
%   Anzahl der Kinematikparameter im Vektor `pkin`
% NQJ [1x1]
%   Anzahl der Minimalkoordinaten der kinematischen Kette
% pkin_names (1x14) cell
%   Namen aller Kinematik-Parameter im Vektor `pkin`

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:05
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = S6PRPRRR7_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 2; 3; 4; 5;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [1; 0; 1; 0; 0; 0;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 1; 1; 1; 1; 1;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 7;
% Anzahl der Kinematikparameter
NKP = 14;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 6;
% Namen der Kinematikparameter
pkin_names = {'a2', 'a3', 'a4', 'a5', 'a6', 'alpha2', 'alpha3', 'alpha4', 'd2', 'd4', 'd5', 'd6', 'theta1', 'theta3'};

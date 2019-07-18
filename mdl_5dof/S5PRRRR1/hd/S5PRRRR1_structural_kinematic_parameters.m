% Return Structural Kinematic Parameters of the Robot 
% S5PRRRR1
%
% Output:
% v_mdh [5x1]
%   Vorgänger-Indizes (0=Basis)
% sigma_mdh [5x1]
%   Dregelenk = 0, Schubgelenk = 1
% mu_mdh [5x1]
%   Aktives Gelenk = 1, Passiv = 0
% NL [1x1]
%   Anzahl der Starrkörper (inklusive Basis)
% NKP [1x1]
%   Anzahl der Kinematikparameter im Vektor `pkin`
% NQJ [1x1]
%   Anzahl der Minimalkoordinaten der kinematischen Kette
% pkin_names (1x2) cell
%   Namen aller Kinematik-Parameter im Vektor `pkin`

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [v_mdh, sigma_mdh, mu_mdh, NL, NKP, NQJ, pkin_names] = S5PRRRR1_structural_kinematic_parameters()

% Aus parameters_mdh_v_matlab.m
t1 = [0; 1; 2; 3; 4;];
v_mdh = uint8(t1);

% Aus parameters_mdh_sigma_matlab.m
t1 = [1; 0; 0; 0; 0;];
sigma_mdh = t1;

% Aus parameters_mdh_mu_matlab.m
t1 = [1; 1; 1; 1; 1;];
mu_mdh = t1;

% Aus Roboterdefinition
% Anzahl der Robotersegmente (inkl Basis)
NL = 6;
% Anzahl der Kinematikparameter
NKP = 2;
% Anzahl der Minimalkoordinaten (für hybride Systeme)
NQJ = 5;
% Namen der Kinematikparameter
pkin_names = {'a2', 'a4'};

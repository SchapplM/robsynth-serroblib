% Inverse Dynamik für komplette Trajektorie für
% S6RPRPRP7
%
% Eingabe:
% Q [NTx6]
%   Trajektorie von Gelenkpositionen (NT Zeitschritte in den Zeilen)
% QD [NTx6]
%   Trajektorie von Gelenkgeschwindigkeiten
% QDD [NTx6]
%   Trajektorie von Gelenkbeschleunigungen
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP7_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx6]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S6RPRPRP7_invdynJ_fixb_mdp_slag_vp_traj(Q, QD, QDD, g, pkin, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,6]),
%$cgargs  coder.newtype('double',[inf,6]),
%$cgargs  coder.newtype('double',[inf,6]),
%$cgargs  zeros(3,1), zeros(9,1), zeros(24,1)}
assert(isreal(Q) && all(size(Q,2) == 6), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp_traj: Q needs to be [NTx6] (double)');
assert(isreal(QD) && all(size(QD,2) == 6), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp_traj: QD needs to be [NTx6] (double)');
assert(isreal(QDD) && all(size(QDD,2) == 6), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp_traj: QDD needs to be [NTx6] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp_traj: Gravity vector g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp_traj: Kinematic parameters pkin have to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPRP7_invdynJ_fixb_mdp_slag_vp_traj: Dynamics parameter vector MDP has to be [24x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(Q));
for k = 1:size(Q,1)
  TAU(k,:) = S6RPRPRP7_invdynJ_fixb_mdp_slag_vp(Q(k,:)', QD(k,:)', QDD(k,:)', g, pkin, MDP);
end

% Inverse Dynamik für komplette Trajektorie für
% S4PRRR2
%
% Eingabe:
% Q [NTx4]
%   Trajektorie von Gelenkpositionen (NT Zeitschritte in den Zeilen)
% QD [NTx4]
%   Trajektorie von Gelenkgeschwindigkeiten
% QDD [NTx4]
%   Trajektorie von Gelenkbeschleunigungen
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR2_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx4]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S4PRRR2_invdynJ_fixb_mdp_slag_vp_traj(Q, QD, QDD, g, pkin, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,4]),
%$cgargs  coder.newtype('double',[inf,4]),
%$cgargs  coder.newtype('double',[inf,4]),
%$cgargs  zeros(3,1), zeros(2,1), zeros(10,1)}
assert(isreal(Q) && all(size(Q,2) == 4), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp_traj: Q needs to be [NTx4] (double)');
assert(isreal(QD) && all(size(QD,2) == 4), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp_traj: QD needs to be [NTx4] (double)');
assert(isreal(QDD) && all(size(QDD,2) == 4), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp_traj: QDD needs to be [NTx4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp_traj: Gravity vector g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp_traj: Kinematic parameters pkin have to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRRR2_invdynJ_fixb_mdp_slag_vp_traj: Dynamics parameter vector MDP has to be [10x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(Q));
for k = 1:size(Q,1)
  TAU(k,:) = S4PRRR2_invdynJ_fixb_mdp_slag_vp(Q(k,:)', QD(k,:)', QDD(k,:)', g, pkin, MDP);
end

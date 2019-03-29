% Inverse Dynamik für komplette Trajektorie für
% S5RRRRR2
%
% Eingabe:
% Q [NTx5]
%   Trajektorie von Gelenkpositionen (NT Zeitschritte in den Zeilen)
% QD [NTx5]
%   Trajektorie von Gelenkgeschwindigkeiten
% QDD [NTx5]
%   Trajektorie von Gelenkbeschleunigungen
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR2_convert_par2_MPV_fixb.m
%
% Ausgabe:
% TAU [NTx5]
%   Time series of inverse Dynamics joint torque

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function TAU = S5RRRRR2_invdynJ_fixb_mdp_slag_vp_traj(Q, QD, QDD, g, pkin, MDP)

%% Coder Information
%#codegen
%$cgargs {coder.newtype('double',[inf,5]),
%$cgargs  coder.newtype('double',[inf,5]),
%$cgargs  coder.newtype('double',[inf,5]),
%$cgargs  zeros(3,1), zeros(2,1), zeros(27,1)}
assert(isreal(Q) && all(size(Q,2) == 5), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp_traj: Q needs to be [NTx5] (double)');
assert(isreal(QD) && all(size(QD,2) == 5), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp_traj: QD needs to be [NTx5] (double)');
assert(isreal(QDD) && all(size(QDD,2) == 5), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp_traj: QDD needs to be [NTx5] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp_traj: Gravity vector g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp_traj: Kinematic parameters pkin have to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR2_invdynJ_fixb_mdp_slag_vp_traj: Dynamics parameter vector MDP has to be [27x1] (double)');

%% Inverse Dynamik für jeden Zeitschritt der Trajektorie berechnen
TAU = NaN(size(Q));
for k = 1:size(Q,1)
  TAU(k,:) = S5RRRRR2_invdynJ_fixb_mdp_slag_vp(Q(k,:)', QD(k,:)', QDD(k,:)', g, pkin, MDP);
end

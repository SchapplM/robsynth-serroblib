% Calculate joint inertia matrix for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPRPP5_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:42
% EndTime: 2019-12-31 18:16:43
% DurationCPUTime: 0.15s
% Computational Cost: add. (136->70), mult. (180->86), div. (0->0), fcn. (96->2), ass. (0->23)
t42 = pkin(3) + pkin(4);
t41 = cos(qJ(3));
t38 = t41 ^ 2;
t40 = sin(qJ(3));
t36 = t40 ^ 2 + t38;
t52 = t40 * qJ(4);
t44 = -pkin(1) - pkin(6);
t51 = qJ(5) + t44;
t50 = -MDP(15) + MDP(20);
t49 = MDP(16) + MDP(19);
t48 = MDP(17) + MDP(21);
t47 = t41 * qJ(4) - qJ(2);
t46 = pkin(3) * MDP(17) + MDP(14);
t45 = qJ(4) ^ 2;
t35 = t41 * pkin(3) + t52;
t34 = t40 * pkin(3) - t47;
t32 = t51 * t41;
t31 = t42 * t41 + t52;
t30 = t51 * t40;
t29 = t36 * t44;
t28 = -t42 * t40 + t47;
t27 = t30 * t40 + t32 * t41;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + t38 * MDP(7) + (t36 * t44 ^ 2 + t34 ^ 2) * MDP(17) + (t28 ^ 2 + t30 ^ 2 + t32 ^ 2) * MDP(21) - 0.2e1 * t29 * MDP(15) + 0.2e1 * t27 * MDP(20) + 0.2e1 * (qJ(2) * MDP(13) - t34 * MDP(16) + t28 * MDP(19)) * t41 + 0.2e1 * (qJ(2) * MDP(12) + t34 * MDP(14) - t28 * MDP(18) - t41 * MDP(8)) * t40; t29 * MDP(17) + t27 * MDP(21) - pkin(1) * MDP(6) + t50 * t36 + MDP(4); t48 * t36 + MDP(6); t41 * MDP(9) - t40 * MDP(10) - t35 * MDP(15) + t32 * MDP(18) + t30 * MDP(19) + t31 * MDP(20) + (t30 * qJ(4) + t32 * t42) * MDP(21) + ((MDP(12) + t46) * t41 + (qJ(4) * MDP(17) - MDP(13) + MDP(16)) * t40) * t44; t35 * MDP(17) + t31 * MDP(21) + (MDP(12) + MDP(14) + MDP(18)) * t41 + (-MDP(13) + t49) * t40; MDP(11) + 0.2e1 * pkin(3) * MDP(14) + (pkin(3) ^ 2 + t45) * MDP(17) + 0.2e1 * t42 * MDP(18) + (t42 ^ 2 + t45) * MDP(21) + 0.2e1 * t49 * qJ(4); -t32 * MDP(21) + (-MDP(17) * t44 - t50) * t41; -t48 * t41; -t42 * MDP(21) - MDP(18) - t46; t48; -t40 * MDP(18) + t41 * MDP(19) + t28 * MDP(21); 0; 0; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

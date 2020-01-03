% Calculate joint inertia matrix for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPRPR10_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:05
% EndTime: 2019-12-31 18:26:06
% DurationCPUTime: 0.12s
% Computational Cost: add. (155->51), mult. (215->73), div. (0->0), fcn. (177->6), ass. (0->29)
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t56 = -pkin(1) - pkin(2);
t41 = qJ(2) * t53 - t55 * t56;
t40 = -pkin(3) - t41;
t42 = qJ(2) * t55 + t53 * t56;
t50 = sin(pkin(8));
t51 = cos(pkin(8));
t33 = t40 * t51 - t42 * t50;
t31 = pkin(4) - t33;
t46 = -pkin(3) * t51 - pkin(4);
t69 = t31 - t46;
t34 = t50 * t40 + t51 * t42;
t68 = MDP(10) * pkin(3);
t67 = t41 * MDP(8);
t66 = t42 * MDP(9);
t52 = sin(qJ(5));
t65 = t52 ^ 2 * MDP(11) + MDP(7);
t54 = cos(qJ(5));
t64 = t54 * MDP(12);
t63 = t54 * MDP(16);
t62 = 0.2e1 * t52 * t64 + t65;
t61 = t55 * MDP(8) - t53 * MDP(9);
t60 = -t52 * MDP(13) - t54 * MDP(14);
t59 = -MDP(17) * t52 + t63;
t58 = -MDP(16) * t52 - MDP(17) * t54;
t39 = t50 * t55 + t51 * t53;
t37 = t50 * t53 - t51 * t55;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t33 ^ 2 + t34 ^ 2) * MDP(10) + 0.2e1 * t59 * t31 + 0.2e1 * t67 + 0.2e1 * t66 + t62; t34 * t39 * MDP(10) - pkin(1) * MDP(6) - MDP(4) + (-MDP(10) * t33 + t59) * t37 - t61; MDP(6) + (t37 ^ 2 + t39 ^ 2) * MDP(10); -t67 - t66 - t69 * t63 + (t33 * t51 + t34 * t50) * t68 + (MDP(17) * t69 - 0.2e1 * t64) * t52 - t65; -t59 * t37 + (-t37 * t51 + t39 * t50) * t68 + t61; (t50 ^ 2 + t51 ^ 2) * MDP(10) * pkin(3) ^ 2 - 0.2e1 * t59 * t46 + t62; 0; 0; 0; MDP(10); t58 * (-pkin(7) + t34) + t60; t58 * t39; t58 * (pkin(3) * t50 + pkin(7)) - t60; t59; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

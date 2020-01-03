% Calculate joint inertia matrix for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RPRR7_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:05
% EndTime: 2019-12-31 16:54:05
% DurationCPUTime: 0.13s
% Computational Cost: add. (142->52), mult. (292->78), div. (0->0), fcn. (284->6), ass. (0->34)
t52 = cos(pkin(7));
t46 = -t52 * pkin(2) - pkin(1);
t73 = 0.2e1 * t46;
t72 = cos(qJ(3));
t71 = pkin(1) * MDP(7);
t51 = sin(pkin(7));
t68 = pkin(5) + qJ(2);
t43 = t68 * t51;
t44 = t68 * t52;
t54 = sin(qJ(3));
t38 = t72 * t43 + t54 * t44;
t42 = t72 * t51 + t54 * t52;
t70 = t38 * t42;
t53 = sin(qJ(4));
t55 = cos(qJ(4));
t69 = t53 * t55;
t66 = t51 * MDP(5);
t65 = t52 * MDP(4);
t41 = t54 * t51 - t72 * t52;
t64 = t41 * MDP(19);
t63 = t42 * MDP(14);
t62 = MDP(16) * t69;
t61 = MDP(17) * t55 - MDP(18) * t53;
t60 = t55 * MDP(20) - t53 * MDP(21);
t59 = -MDP(20) * t53 - MDP(21) * t55;
t58 = MDP(13) + t60;
t57 = t53 * MDP(17) + t55 * MDP(18) + t59 * pkin(6);
t50 = t55 ^ 2;
t49 = t53 ^ 2;
t39 = -t54 * t43 + t72 * t44;
t37 = t41 * pkin(3) - t42 * pkin(6) + t46;
t36 = t53 * t37 + t55 * t39;
t35 = t55 * t37 - t53 * t39;
t1 = [t63 * t73 + MDP(1) + (0.2e1 * t65 - 0.2e1 * t66 + t71) * pkin(1) + (t50 * MDP(15) + MDP(8) - 0.2e1 * t62) * t42 ^ 2 + (MDP(13) * t73 + t64 + 0.2e1 * (-MDP(9) + t61) * t42) * t41 + 0.2e1 * (t35 * t41 + t53 * t70) * MDP(20) + 0.2e1 * (-t36 * t41 + t55 * t70) * MDP(21) + (MDP(7) * qJ(2) + (2 * MDP(6))) * (t51 ^ 2 + t52 ^ 2) * qJ(2); t58 * t41 + t63 - t65 + t66 - t71; MDP(7); -t39 * MDP(14) - t58 * t38 + (-MDP(11) + t57) * t41 + (MDP(10) + MDP(15) * t69 + (-t49 + t50) * MDP(16) + t59 * pkin(3)) * t42; 0; t49 * MDP(15) + 0.2e1 * pkin(3) * t60 + MDP(12) + 0.2e1 * t62; t35 * MDP(20) - t36 * MDP(21) + t61 * t42 + t64; t60; t57; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;

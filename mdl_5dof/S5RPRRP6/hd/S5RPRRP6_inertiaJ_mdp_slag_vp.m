% Calculate joint inertia matrix for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP6_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:17
% EndTime: 2019-12-31 18:43:18
% DurationCPUTime: 0.19s
% Computational Cost: add. (201->88), mult. (373->131), div. (0->0), fcn. (296->6), ass. (0->41)
t82 = 2 * MDP(19);
t59 = sin(pkin(8));
t51 = t59 * pkin(1) + pkin(6);
t61 = sin(qJ(4));
t81 = t51 * t61;
t63 = cos(qJ(4));
t80 = t51 * t63;
t64 = cos(qJ(3));
t79 = t51 * t64;
t78 = t61 * t63;
t77 = -qJ(5) - pkin(7);
t76 = MDP(20) * pkin(4);
t62 = sin(qJ(3));
t75 = qJ(5) * t62;
t74 = MDP(18) * t63;
t54 = -t63 * pkin(4) - pkin(3);
t73 = t54 * MDP(20);
t72 = t61 * MDP(15);
t71 = t63 * t79;
t70 = MDP(13) * t78;
t60 = cos(pkin(8));
t52 = -t60 * pkin(1) - pkin(2);
t69 = -MDP(19) * pkin(4) + MDP(14);
t48 = -t64 * pkin(3) - t62 * pkin(7) + t52;
t46 = t63 * t48;
t42 = -t63 * t75 + t46 + (-pkin(4) - t81) * t64;
t43 = t71 + (t48 - t75) * t61;
t68 = -t42 * t61 + t43 * t63;
t49 = t77 * t61;
t50 = t77 * t63;
t67 = -t49 * t61 - t50 * t63;
t66 = t63 * MDP(17) - t61 * MDP(18);
t58 = t64 ^ 2;
t57 = t63 ^ 2;
t56 = t62 ^ 2;
t55 = t61 ^ 2;
t53 = t57 * t62;
t47 = (pkin(4) * t61 + t51) * t62;
t45 = t61 * t48 + t71;
t44 = -t61 * t79 + t46;
t1 = [MDP(1) - 0.2e1 * t52 * t64 * MDP(10) + t58 * MDP(16) + (t42 ^ 2 + t43 ^ 2 + t47 ^ 2) * MDP(20) + (t59 ^ 2 + t60 ^ 2) * MDP(4) * pkin(1) ^ 2 + (t57 * MDP(12) + MDP(5) - 0.2e1 * t70) * t56 + 0.2e1 * (-t44 * t64 + t56 * t81) * MDP(17) + 0.2e1 * (t45 * t64 + t56 * t80) * MDP(18) + (0.2e1 * t52 * MDP(11) + (-t42 * t63 - t43 * t61) * t82 + 0.2e1 * (-t63 * MDP(14) + MDP(6) + t72) * t64) * t62; (-t47 * t64 + t68 * t62) * MDP(20); MDP(4) + (t58 + (t55 + t57) * t56) * MDP(20); t53 * MDP(13) + t68 * MDP(19) + (t42 * t49 - t43 * t50 + t47 * t54) * MDP(20) + (-t51 * MDP(11) - t61 * MDP(14) - t63 * MDP(15) + MDP(8) + (MDP(17) * t61 + t74) * pkin(7)) * t64 + (MDP(7) - t51 * MDP(10) + MDP(12) * t78 - t55 * MDP(13) + (-pkin(3) * t61 - t80) * MDP(17) + (-pkin(3) * t63 + t81) * MDP(18) + (-t49 * t63 + t50 * t61) * MDP(19)) * t62; t53 * MDP(19) + (t55 * MDP(19) + t67 * MDP(20) - MDP(11)) * t62 + (MDP(10) + t66 - t73) * t64; MDP(9) + t55 * MDP(12) + 0.2e1 * t70 + t67 * t82 + (t49 ^ 2 + t50 ^ 2 + t54 ^ 2) * MDP(20) + 0.2e1 * t66 * pkin(3); t42 * t76 - t64 * MDP(16) + t44 * MDP(17) - t45 * MDP(18) + (t69 * t63 - t72) * t62; (-t74 + (-MDP(17) - t76) * t61) * t62; t49 * t76 + (-MDP(18) * pkin(7) + MDP(15)) * t63 + (-MDP(17) * pkin(7) + t69) * t61; MDP(20) * pkin(4) ^ 2 + MDP(16); t47 * MDP(20); -t64 * MDP(20); t73; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

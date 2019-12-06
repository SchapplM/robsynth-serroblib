% Calculate joint inertia matrix for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRPRR6_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:55
% EndTime: 2019-12-05 15:57:56
% DurationCPUTime: 0.23s
% Computational Cost: add. (210->72), mult. (464->117), div. (0->0), fcn. (502->10), ass. (0->49)
t101 = MDP(8) * qJ(3);
t71 = cos(pkin(10));
t64 = -t71 * pkin(3) - pkin(2);
t100 = 0.2e1 * t64;
t99 = cos(qJ(4));
t98 = pkin(2) * MDP(8);
t69 = sin(pkin(10));
t74 = sin(qJ(4));
t60 = t99 * t69 + t74 * t71;
t73 = sin(qJ(5));
t97 = t60 * t73;
t76 = cos(qJ(5));
t96 = t60 * t76;
t70 = sin(pkin(5));
t75 = sin(qJ(2));
t95 = t70 * t75;
t77 = cos(qJ(2));
t94 = t70 * t77;
t93 = t73 * t76;
t92 = pkin(7) + qJ(3);
t90 = t69 * MDP(6);
t89 = t71 * MDP(5);
t59 = t74 * t69 - t99 * t71;
t88 = t59 * MDP(20);
t87 = t60 * MDP(15);
t86 = MDP(17) * t93;
t84 = MDP(18) * t76 - MDP(19) * t73;
t83 = t76 * MDP(21) - t73 * MDP(22);
t82 = -MDP(21) * t73 - MDP(22) * t76;
t81 = MDP(14) + t83;
t80 = t87 - t89 + t90 - t98;
t79 = t73 * MDP(18) + t76 * MDP(19) + t82 * pkin(8);
t72 = cos(pkin(5));
t68 = t76 ^ 2;
t67 = t73 ^ 2;
t62 = t92 * t71;
t61 = t92 * t69;
t57 = t72 * t69 + t71 * t95;
t56 = -t69 * t95 + t72 * t71;
t55 = -t74 * t61 + t99 * t62;
t54 = t99 * t61 + t74 * t62;
t53 = t59 * pkin(4) - t60 * pkin(8) + t64;
t52 = t74 * t56 + t99 * t57;
t51 = -t99 * t56 + t74 * t57;
t50 = t76 * t52 - t73 * t94;
t49 = -t73 * t52 - t76 * t94;
t48 = t73 * t53 + t76 * t55;
t47 = t76 * t53 - t73 * t55;
t1 = [MDP(1) + (t70 ^ 2 * t77 ^ 2 + t56 ^ 2 + t57 ^ 2) * MDP(8); (t49 * t59 + t51 * t97) * MDP(21) + (-t50 * t59 + t51 * t96) * MDP(22) + (-t75 * MDP(4) + (-t59 * MDP(14) + MDP(3) - t80) * t77) * t70 + (MDP(7) + t101) * (-t56 * t69 + t57 * t71); t87 * t100 + MDP(2) + (0.2e1 * t89 - 0.2e1 * t90 + t98) * pkin(2) + (t68 * MDP(16) + MDP(9) - 0.2e1 * t86) * t60 ^ 2 + (MDP(14) * t100 + t88 + 0.2e1 * (-MDP(10) + t84) * t60) * t59 + 0.2e1 * (t47 * t59 + t54 * t97) * MDP(21) + 0.2e1 * (-t48 * t59 + t54 * t96) * MDP(22) + (0.2e1 * MDP(7) + t101) * (t69 ^ 2 + t71 ^ 2) * qJ(3); -MDP(8) * t94; t81 * t59 + t80; MDP(8); -t52 * MDP(15) - t81 * t51; -t55 * MDP(15) - t81 * t54 + (-MDP(12) + t79) * t59 + (MDP(11) + MDP(16) * t93 + (-t67 + t68) * MDP(17) + t82 * pkin(4)) * t60; 0; t67 * MDP(16) + 0.2e1 * pkin(4) * t83 + MDP(13) + 0.2e1 * t86; t49 * MDP(21) - t50 * MDP(22); t47 * MDP(21) - t48 * MDP(22) + t84 * t60 + t88; t83; t79; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

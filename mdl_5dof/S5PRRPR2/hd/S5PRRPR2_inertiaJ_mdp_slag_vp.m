% Calculate joint inertia matrix for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRRPR2_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:36
% EndTime: 2019-12-05 16:17:36
% DurationCPUTime: 0.17s
% Computational Cost: add. (148->64), mult. (306->89), div. (0->0), fcn. (207->6), ass. (0->47)
t80 = sin(qJ(5));
t82 = cos(qJ(5));
t110 = MDP(17) * t82 - MDP(18) * t80;
t79 = cos(pkin(9));
t109 = 0.2e1 * t79;
t108 = 2 * MDP(10);
t107 = 0.2e1 * MDP(17);
t106 = 0.2e1 * MDP(18);
t83 = cos(qJ(3));
t105 = pkin(2) * t83;
t71 = -pkin(3) - t105;
t104 = pkin(3) - t71;
t81 = sin(qJ(3));
t70 = pkin(2) * t81 + qJ(4);
t103 = t70 * t79;
t78 = sin(pkin(9));
t76 = t78 ^ 2;
t62 = t76 * t70;
t102 = t76 * t82;
t61 = -pkin(4) * t79 - pkin(7) * t78 - pkin(3);
t55 = t61 - t105;
t51 = t82 * t103 + t55 * t80;
t101 = t70 * t102 + t51 * t79;
t95 = qJ(4) * t79;
t54 = t61 * t80 + t82 * t95;
t73 = t76 * qJ(4);
t100 = t54 * t79 + t82 * t73;
t77 = t79 ^ 2;
t99 = t77 * t70 + t62;
t98 = t77 * qJ(4) + t73;
t97 = t76 + t77;
t96 = pkin(3) * MDP(11);
t74 = t78 * MDP(9);
t92 = t71 * MDP(11);
t91 = MDP(15) * t78 * t80;
t66 = t82 * t78 * MDP(14);
t90 = t97 * MDP(11);
t89 = MDP(8) * t109 - 0.2e1 * t74;
t88 = t82 ^ 2 * t76 * MDP(12) - 0.2e1 * t80 * MDP(13) * t102 + t77 * MDP(16) + t91 * t109 - 0.2e1 * t79 * t66 + MDP(5);
t87 = (t83 * MDP(6) - t81 * MDP(7)) * pkin(2);
t86 = -MDP(16) * t79 + t66 - t91;
t85 = t74 + (-MDP(8) - t110) * t79;
t67 = t80 * t73;
t56 = t80 * t62;
t53 = t61 * t82 - t80 * t95;
t50 = -t80 * t103 + t55 * t82;
t1 = [MDP(1) + t90; 0; MDP(2) + t70 ^ 2 * t90 + (-t89 + t92) * t71 + 0.2e1 * t87 + t99 * t108 + (-t50 * t79 + t56) * t107 + t101 * t106 + t88; 0; (t98 + t99) * MDP(10) + (t97 * t70 * qJ(4) - pkin(3) * t71) * MDP(11) + (t56 + t67) * MDP(17) + (t100 + t101) * MDP(18) - t104 * t74 + t87 + (t104 * MDP(8) + (-t50 - t53) * MDP(17)) * t79 + t88; qJ(4) ^ 2 * t90 + (t89 + t96) * pkin(3) + t98 * t108 + (-t53 * t79 + t67) * t107 + t100 * t106 + t88; 0; t85 + t92; t85 - t96; MDP(11); (-MDP(17) * t80 - MDP(18) * t82) * t78; MDP(17) * t50 - MDP(18) * t51 + t86; MDP(17) * t53 - MDP(18) * t54 + t86; t110; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

% Calculate joint inertia matrix for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR6_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:49
% EndTime: 2019-12-05 16:32:50
% DurationCPUTime: 0.32s
% Computational Cost: add. (293->106), mult. (668->169), div. (0->0), fcn. (695->10), ass. (0->59)
t107 = (MDP(15) * qJ(4));
t127 = MDP(14) + t107;
t126 = 2 * MDP(14);
t125 = -2 * MDP(17);
t124 = 2 * MDP(22);
t87 = sin(pkin(10));
t123 = pkin(7) * t87;
t95 = cos(qJ(3));
t122 = pkin(7) * t95;
t92 = sin(qJ(3));
t121 = pkin(8) * t92;
t88 = sin(pkin(5));
t93 = sin(qJ(2));
t120 = t88 * t93;
t96 = cos(qJ(2));
t119 = t88 * t96;
t118 = pkin(8) + qJ(4);
t78 = -t95 * pkin(3) - t92 * qJ(4) - pkin(2);
t89 = cos(pkin(10));
t68 = t89 * t122 + t87 * t78;
t116 = pkin(3) * MDP(15);
t115 = MDP(11) * t92;
t91 = sin(qJ(5));
t94 = cos(qJ(5));
t76 = t94 * t87 + t91 * t89;
t69 = t76 * t92;
t114 = t69 * MDP(19);
t75 = t91 * t87 - t94 * t89;
t70 = t75 * t92;
t113 = t70 * MDP(16);
t112 = t70 * MDP(18);
t111 = t75 * MDP(21);
t110 = t87 * MDP(13);
t109 = t89 * MDP(12);
t108 = t95 * MDP(20);
t104 = t87 * MDP(12) + t89 * MDP(13);
t90 = cos(pkin(5));
t72 = t120 * t95 + t90 * t92;
t62 = -t119 * t89 - t72 * t87;
t63 = -t119 * t87 + t72 * t89;
t103 = (t94 * t62 - t91 * t63) * MDP(21) - (t91 * t62 + t94 * t63) * MDP(22);
t102 = t69 * MDP(21) - t70 * MDP(22);
t101 = -t109 + t110 - t116;
t100 = pkin(7) * MDP(15) + t104;
t79 = t118 * t87;
t80 = t118 * t89;
t99 = t76 * MDP(18) - t75 * MDP(19) + (-t94 * t79 - t91 * t80) * MDP(21) - (-t91 * t79 + t94 * t80) * MDP(22);
t98 = t76 * MDP(22) + t101 + t111;
t86 = t92 ^ 2;
t83 = -t89 * pkin(4) - pkin(3);
t77 = (pkin(4) * t87 + pkin(7)) * t92;
t74 = t89 * t78;
t71 = t120 * t92 - t90 * t95;
t67 = -t87 * t122 + t74;
t64 = -t87 * t121 + t68;
t61 = -t89 * t121 + t74 + (-pkin(4) - t123) * t95;
t58 = t91 * t61 + t94 * t64;
t57 = t94 * t61 - t91 * t64;
t1 = [MDP(1) + (t62 ^ 2 + t63 ^ 2 + t71 ^ 2) * MDP(15); (t62 * t67 + t63 * t68) * MDP(15) + t102 * t71 + (-t62 * MDP(12) + t63 * MDP(13) - t103) * t95 + ((-t62 * t89 - t63 * t87) * MDP(14) + t100 * t71) * t92 + (-t93 * MDP(4) + (MDP(10) * t95 + MDP(3) - t115) * t96) * t88; MDP(2) + t86 * MDP(5) - 0.2e1 * pkin(2) * t115 + (t86 * pkin(7) ^ 2 + t67 ^ 2 + t68 ^ 2) * MDP(15) - (t69 * t125 - t113) * t70 + (0.2e1 * pkin(2) * MDP(10) + 0.2e1 * t92 * MDP(6) + t108 + 0.2e1 * t112 + 0.2e1 * t114) * t95 + 0.2e1 * (t86 * t123 - t67 * t95) * MDP(12) + 0.2e1 * (t86 * pkin(7) * t89 + t68 * t95) * MDP(13) + 0.2e1 * (-t57 * t95 + t77 * t69) * MDP(21) + (t58 * t95 - t77 * t70) * t124 + (-t67 * t89 - t68 * t87) * t92 * t126; -t72 * MDP(11) + (-MDP(10) + t98) * t71 + t127 * (-t62 * t87 + t63 * t89); -t76 * t113 + (-t76 * t69 + t70 * t75) * MDP(17) + (t83 * t69 + t77 * t75) * MDP(21) + (-t83 * t70 + t77 * t76) * MDP(22) + (-pkin(7) * MDP(11) + qJ(4) * t104 + MDP(8) - t99) * t95 + (MDP(7) - t104 * pkin(3) + (-MDP(10) + t101) * pkin(7)) * t92 + t127 * (-t67 * t87 + t68 * t89); 0.2e1 * t83 * t111 + MDP(9) + (0.2e1 * t109 - 0.2e1 * t110 + t116) * pkin(3) + (MDP(16) * t76 + t83 * t124 + t75 * t125) * t76 + (t126 + t107) * (t87 ^ 2 + t89 ^ 2) * qJ(4); t71 * MDP(15); t100 * t92 + t102; t98; MDP(15); t103; t57 * MDP(21) - t58 * MDP(22) - t108 - t112 - t114; t99; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

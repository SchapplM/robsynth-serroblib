% Calculate joint inertia matrix for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRPRRP3_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:56
% EndTime: 2019-03-08 20:07:57
% DurationCPUTime: 0.36s
% Computational Cost: add. (522->122), mult. (1081->180), div. (0->0), fcn. (1207->10), ass. (0->66)
t133 = cos(qJ(4));
t92 = sin(pkin(11));
t94 = cos(pkin(11));
t97 = sin(qJ(4));
t80 = t133 * t92 + t97 * t94;
t136 = 0.2e1 * t80;
t135 = MDP(8) * qJ(3);
t86 = -t94 * pkin(3) - pkin(2);
t134 = 0.2e1 * t86;
t132 = pkin(2) * MDP(8);
t93 = sin(pkin(6));
t98 = sin(qJ(2));
t131 = t93 * t98;
t96 = sin(qJ(5));
t99 = cos(qJ(5));
t130 = t96 * t99;
t128 = pkin(8) + qJ(3);
t81 = t128 * t92;
t82 = t128 * t94;
t75 = t133 * t82 - t97 * t81;
t129 = t99 * t75;
t127 = -qJ(6) - pkin(9);
t90 = t96 ^ 2;
t91 = t99 ^ 2;
t125 = t90 + t91;
t124 = MDP(24) * pkin(5);
t123 = qJ(6) * t80;
t100 = cos(qJ(2));
t122 = t100 * t93;
t121 = t92 * MDP(6);
t120 = t94 * MDP(5);
t119 = MDP(19) * t96;
t95 = cos(pkin(6));
t76 = -t92 * t131 + t95 * t94;
t77 = t94 * t131 + t95 * t92;
t71 = t133 * t77 + t97 * t76;
t68 = -t96 * t122 + t99 * t71;
t118 = t68 * MDP(22);
t79 = -t133 * t94 + t97 * t92;
t117 = t79 * MDP(20);
t116 = t80 * MDP(15);
t87 = -t99 * pkin(5) - pkin(4);
t115 = t87 * MDP(24);
t114 = t96 * MDP(22);
t113 = MDP(17) * t130;
t73 = t79 * pkin(4) - t80 * pkin(9) + t86;
t65 = t99 * t73 - t96 * t75;
t112 = -MDP(23) * pkin(5) + MDP(18);
t111 = MDP(21) + t124;
t63 = t79 * pkin(5) - t99 * t123 + t65;
t64 = t129 + (t73 - t123) * t96;
t110 = t63 * t99 + t64 * t96;
t67 = -t99 * t122 - t96 * t71;
t109 = t67 * t99 + t68 * t96;
t83 = t127 * t96;
t84 = t127 * t99;
t107 = t99 * t83 - t96 * t84;
t106 = t65 * MDP(21) - (t96 * t73 + t129) * MDP(22);
t105 = t99 * MDP(21) - t114;
t104 = MDP(21) * t96 + MDP(22) * t99;
t74 = t133 * t81 + t97 * t82;
t103 = MDP(14) + t105;
t102 = -t120 + t121 - t132;
t70 = -t133 * t76 + t97 * t77;
t69 = t96 * t80 * pkin(5) + t74;
t1 = [MDP(1) + (t93 ^ 2 * t100 ^ 2 + t76 ^ 2 + t77 ^ 2) * MDP(8) + (t67 ^ 2 + t68 ^ 2 + t70 ^ 2) * MDP(24); (t67 * t63 + t68 * t64 + t70 * t69) * MDP(24) + (t67 * MDP(21) - t118) * t79 + (-t109 * MDP(23) + t104 * t70) * t80 + (-t98 * MDP(4) + (-t79 * MDP(14) + MDP(3) - t102 - t116) * t100) * t93 + (MDP(7) + t135) * (-t76 * t92 + t77 * t94); MDP(2) + t116 * t134 + (t63 ^ 2 + t64 ^ 2 + t69 ^ 2) * MDP(24) + (0.2e1 * t120 - 0.2e1 * t121 + t132) * pkin(2) + (t91 * MDP(16) + MDP(9) - 0.2e1 * t113) * t80 ^ 2 + (-t110 * MDP(23) + t104 * t74) * t136 + (MDP(14) * t134 + t117 + (MDP(18) * t99 - MDP(10) - t119) * t136 + 0.2e1 * t106) * t79 + (0.2e1 * MDP(7) + t135) * (t92 ^ 2 + t94 ^ 2) * qJ(3); t109 * MDP(24) - MDP(8) * t122; t110 * MDP(24) + (-t125 * MDP(23) + MDP(15)) * t80 + t103 * t79 + t102; t125 * MDP(24) + MDP(8); -t71 * MDP(15) + (-t67 * t96 + t68 * t99) * MDP(23) + (t67 * t83 - t68 * t84) * MDP(24) + (-t103 + t115) * t70; -t75 * MDP(15) + (-t63 * t96 + t64 * t99) * MDP(23) + (t63 * t83 - t64 * t84 + t69 * t87) * MDP(24) - t103 * t74 + (t96 * MDP(18) + t99 * MDP(19) - t104 * pkin(9) - MDP(12)) * t79 + (MDP(11) + MDP(16) * t130 + (-t90 + t91) * MDP(17) - t107 * MDP(23) - t104 * pkin(4)) * t80; t107 * MDP(24); MDP(13) + t90 * MDP(16) + 0.2e1 * t113 + 0.2e1 * (-t83 * t96 - t84 * t99) * MDP(23) + (t83 ^ 2 + t84 ^ 2 + t87 ^ 2) * MDP(24) + 0.2e1 * t105 * pkin(4); t111 * t67 - t118; t63 * t124 + t117 + (t112 * t99 - t119) * t80 + t106; t111 * t99 - t114; t83 * t124 + (-MDP(22) * pkin(9) + MDP(19)) * t99 + (-MDP(21) * pkin(9) + t112) * t96; MDP(24) * pkin(5) ^ 2 + MDP(20); t70 * MDP(24); t69 * MDP(24); 0; t115; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

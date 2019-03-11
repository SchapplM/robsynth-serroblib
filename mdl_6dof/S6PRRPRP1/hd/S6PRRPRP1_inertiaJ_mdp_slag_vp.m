% Calculate joint inertia matrix for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6PRRPRP1_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:26:45
% EndTime: 2019-03-08 21:26:46
% DurationCPUTime: 0.40s
% Computational Cost: add. (575->133), mult. (1169->206), div. (0->0), fcn. (1302->10), ass. (0->63)
t102 = cos(qJ(3));
t95 = sin(pkin(11));
t97 = cos(pkin(11));
t99 = sin(qJ(3));
t85 = t102 * t95 + t97 * t99;
t134 = 0.2e1 * t85;
t123 = cos(pkin(6));
t100 = sin(qJ(2));
t96 = sin(pkin(6));
t127 = t100 * t96;
t105 = t123 * t102 - t99 * t127;
t80 = t102 * t127 + t123 * t99;
t71 = -t105 * t97 + t80 * t95;
t133 = t71 ^ 2;
t98 = sin(qJ(5));
t132 = t85 * t98;
t131 = -qJ(4) - pkin(8);
t93 = t98 ^ 2;
t101 = cos(qJ(5));
t94 = t101 ^ 2;
t130 = t93 + t94;
t129 = MDP(22) * pkin(5);
t128 = qJ(6) * t85;
t114 = t131 * t99;
t88 = t131 * t102;
t78 = t95 * t114 - t97 * t88;
t126 = t101 * t78;
t103 = cos(qJ(2));
t125 = t103 * t96;
t90 = pkin(3) * t95 + pkin(9);
t124 = qJ(6) + t90;
t122 = MDP(17) * t98;
t84 = -t97 * t102 + t95 * t99;
t121 = MDP(18) * t84;
t73 = t105 * t95 + t97 * t80;
t69 = t101 * t73 - t98 * t125;
t120 = t69 * MDP(20);
t91 = -pkin(3) * t97 - pkin(4);
t87 = -pkin(5) * t101 + t91;
t119 = t87 * MDP(22);
t92 = -pkin(3) * t102 - pkin(2);
t118 = t92 * MDP(13);
t117 = t98 * MDP(20);
t116 = MDP(10) * t102;
t115 = t98 * t101 * MDP(15);
t75 = pkin(4) * t84 - pkin(9) * t85 + t92;
t66 = t101 * t75 - t78 * t98;
t76 = -t97 * t114 - t88 * t95;
t113 = -MDP(21) * pkin(5) + MDP(16);
t112 = MDP(19) + t129;
t111 = -t84 * t90 + t85 * t91;
t64 = pkin(5) * t84 - t101 * t128 + t66;
t65 = t126 + (t75 - t128) * t98;
t110 = t101 * t64 + t65 * t98;
t68 = -t101 * t125 - t73 * t98;
t109 = t101 * t68 + t69 * t98;
t108 = t66 * MDP(19) - (t75 * t98 + t126) * MDP(20);
t107 = t101 * MDP(19) - t117;
t106 = t98 * MDP(19) + t101 * MDP(20) + MDP(12);
t82 = t124 * t101;
t81 = t124 * t98;
t70 = pkin(5) * t132 + t76;
t1 = [MDP(1) + (t103 ^ 2 * t96 ^ 2 + t73 ^ 2 + t133) * MDP(13) + (t68 ^ 2 + t69 ^ 2 + t133) * MDP(22); (t71 * t76 + t73 * t78) * MDP(13) + (t64 * t68 + t65 * t69 + t70 * t71) * MDP(22) + (-t73 * MDP(12) + t68 * MDP(19) - t120) * t84 + (-t109 * MDP(21) + t106 * t71) * t85 + (-t100 * MDP(4) + (-MDP(11) * t99 + MDP(3) + t116 - t118) * t103) * t96; MDP(2) + 0.2e1 * pkin(2) * t116 + (t76 ^ 2 + t78 ^ 2 + t92 ^ 2) * MDP(13) + (t64 ^ 2 + t65 ^ 2 + t70 ^ 2) * MDP(22) + (t94 * MDP(14) - 0.2e1 * t115) * t85 ^ 2 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t99 + 0.2e1 * t102 * MDP(6)) * t99 + (t121 + (MDP(16) * t101 - t122) * t134) * t84 + 0.2e1 * (-t78 * MDP(12) + t108) * t84 + (-t110 * MDP(21) + t106 * t76) * t134; t105 * MDP(10) - t80 * MDP(11) + (t101 * t69 - t68 * t98) * MDP(21) + (-t68 * t81 + t69 * t82) * MDP(22) + (-t107 + t119) * t71 + (-t71 * t97 + t73 * t95) * MDP(13) * pkin(3); t99 * MDP(7) + t102 * MDP(8) + (-t64 * t81 + t65 * t82 + t70 * t87) * MDP(22) + (-t93 + t94) * MDP(15) * t85 + (-t99 * MDP(10) - t102 * MDP(11)) * pkin(8) + (t84 * MDP(16) + t111 * MDP(19) + t76 * MDP(20) + (-t82 * t85 - t64) * MDP(21)) * t98 + (MDP(14) * t132 + t84 * MDP(17) - t76 * MDP(19) + t111 * MDP(20) + (t81 * t85 + t65) * MDP(21)) * t101 + ((-t84 * t95 - t85 * t97) * MDP(12) + (-t76 * t97 + t78 * t95) * MDP(13)) * pkin(3); MDP(9) + t93 * MDP(14) + 0.2e1 * t115 + 0.2e1 * (t101 * t82 + t81 * t98) * MDP(21) + (t81 ^ 2 + t82 ^ 2 + t87 ^ 2) * MDP(22) - 0.2e1 * t107 * t91 + (t95 ^ 2 + t97 ^ 2) * MDP(13) * pkin(3) ^ 2; -MDP(13) * t125 + t109 * MDP(22); -t130 * MDP(21) * t85 + t110 * MDP(22) + t107 * t84 + t118; (-t101 * t81 + t82 * t98) * MDP(22); t130 * MDP(22) + MDP(13); t112 * t68 - t120; t64 * t129 + t121 + (t113 * t101 - t122) * t85 + t108; -t81 * t129 + (-t90 * MDP(20) + MDP(17)) * t101 + (-t90 * MDP(19) + t113) * t98; t112 * t101 - t117; MDP(22) * pkin(5) ^ 2 + MDP(18); t71 * MDP(22); t70 * MDP(22); t119; 0; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

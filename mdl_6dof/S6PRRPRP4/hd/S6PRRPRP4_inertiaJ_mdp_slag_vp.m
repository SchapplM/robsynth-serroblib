% Calculate joint inertia matrix for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPRP4_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:43:59
% EndTime: 2019-03-08 21:44:00
% DurationCPUTime: 0.38s
% Computational Cost: add. (346->130), mult. (663->173), div. (0->0), fcn. (618->8), ass. (0->68)
t140 = pkin(4) + pkin(8);
t101 = cos(qJ(3));
t95 = sin(pkin(6));
t99 = sin(qJ(2));
t136 = t95 * t99;
t96 = cos(pkin(6));
t98 = sin(qJ(3));
t78 = t101 * t136 + t96 * t98;
t139 = t78 ^ 2;
t138 = -0.2e1 * t98;
t103 = -pkin(3) - pkin(9);
t85 = t140 * t98;
t97 = sin(qJ(5));
t137 = t85 * t97;
t86 = t140 * t101;
t135 = MDP(15) * pkin(8);
t134 = MDP(23) * pkin(5);
t133 = MDP(24) * pkin(5);
t132 = pkin(3) * MDP(15);
t100 = cos(qJ(5));
t131 = t100 * t97;
t102 = cos(qJ(2));
t130 = t102 * t95;
t117 = -qJ(4) * t98 - pkin(2);
t84 = -pkin(3) * t101 + t117;
t129 = MDP(15) * t84;
t77 = -t101 * t96 + t98 * t136;
t73 = t100 * t130 - t77 * t97;
t128 = MDP(22) * t73;
t88 = pkin(5) * t97 + qJ(4);
t127 = t88 * MDP(24);
t126 = t97 * MDP(21);
t125 = t97 * MDP(22);
t124 = -qJ(6) + t103;
t123 = MDP(19) * t100;
t122 = qJ(4) * MDP(15);
t121 = t100 * MDP(22);
t120 = pkin(8) ^ 2 * MDP(15);
t119 = -MDP(11) + MDP(14);
t118 = MDP(17) * t131;
t80 = t103 * t101 + t117;
t116 = qJ(6) * t101 - t80;
t115 = MDP(13) - t132;
t114 = MDP(21) + t133;
t113 = MDP(21) * t103 + MDP(18);
t81 = t100 * t85;
t67 = pkin(5) * t98 + t116 * t97 + t81;
t68 = -t116 * t100 + t137;
t112 = t100 * t67 + t68 * t97;
t111 = -MDP(10) + t115;
t110 = (-t80 * t97 + t81) * MDP(21) - (t100 * t80 + t137) * MDP(22);
t109 = (-MDP(22) * t103 - MDP(19)) * t97;
t108 = t121 + t126;
t107 = MDP(21) * t100 - t125;
t106 = MDP(12) + t107;
t105 = t106 + t135;
t94 = t101 ^ 2;
t93 = t100 ^ 2;
t92 = t98 ^ 2;
t91 = t97 ^ 2;
t87 = t91 + t93;
t83 = t124 * t100;
t82 = t124 * t97;
t76 = pkin(5) * t100 * t101 + t86;
t72 = t100 * t77 + t97 * t130;
t71 = t100 * t83 + t82 * t97;
t66 = t100 * t72 - t73 * t97;
t1 = [MDP(1) + (t102 ^ 2 * t95 ^ 2 + t77 ^ 2 + t139) * MDP(15) + (t72 ^ 2 + t73 ^ 2 + t139) * MDP(24); (t67 * t72 - t68 * t73 + t76 * t78) * MDP(24) + (MDP(12) * t77 + MDP(21) * t72 + t128) * t98 + ((t100 * t73 + t72 * t97) * MDP(23) + t106 * t78) * t101 + (t101 * t78 + t77 * t98) * t135 + (-t99 * MDP(4) + (-t129 + MDP(3) + t119 * t98 + (MDP(10) - MDP(13)) * t101) * t102) * t95; MDP(2) + pkin(2) * MDP(11) * t138 + (t67 ^ 2 + t68 ^ 2 + t76 ^ 2) * MDP(24) + (MDP(14) * t138 + t129) * t84 + (t91 * MDP(16) + 0.2e1 * t118 + t120) * t94 + (MDP(20) + MDP(5) + t120) * t92 + 0.2e1 * t110 * t98 + 0.2e1 * (t92 + t94) * MDP(12) * pkin(8) + 0.2e1 * (pkin(2) * MDP(10) + t84 * MDP(13) + (-MDP(18) * t97 + MDP(6) - t123) * t98 + (-t100 * t68 + t67 * t97) * MDP(23) + t107 * t86) * t101; -t66 * MDP(23) + (t72 * t83 - t73 * t82) * MDP(24) + t111 * t77 + (t108 + t119 + t122 + t127) * t78; -t112 * MDP(23) + (t67 * t83 + t68 * t82 + t76 * t88) * MDP(24) + t108 * t86 + (-pkin(3) * MDP(12) + t111 * pkin(8) + t113 * t100 + MDP(7) + t109) * t98 + (MDP(8) - MDP(16) * t131 + (t91 - t93) * MDP(17) + (-t100 * t82 + t83 * t97) * MDP(23) + t119 * pkin(8) + t105 * qJ(4)) * t101; MDP(9) + t93 * MDP(16) - 0.2e1 * t118 - 0.2e1 * t71 * MDP(23) + (t82 ^ 2 + t83 ^ 2 + t88 ^ 2) * MDP(24) + (-0.2e1 * MDP(13) + t132) * pkin(3) + (0.2e1 * MDP(14) + 0.2e1 * t121 + t122 + 0.2e1 * t126) * qJ(4); MDP(15) * t77 + MDP(24) * t66; t112 * MDP(24) + t105 * t98; -MDP(23) * t87 + MDP(24) * t71 + t115; MDP(24) * t87 + MDP(15); t114 * t72 + t128; t67 * t133 + MDP(20) * t98 + (-t123 + (-MDP(18) + t134) * t97) * t101 + t110; t83 * t133 + t109 + (t113 - t134) * t100; t114 * t100 - t125; MDP(24) * pkin(5) ^ 2 + MDP(20); t78 * MDP(24); t76 * MDP(24); t127; 0; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

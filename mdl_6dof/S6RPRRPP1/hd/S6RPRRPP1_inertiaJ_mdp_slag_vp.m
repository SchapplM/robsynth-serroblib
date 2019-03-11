% Calculate joint inertia matrix for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRRPP1_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:29:53
% EndTime: 2019-03-09 04:29:55
% DurationCPUTime: 0.53s
% Computational Cost: add. (897->169), mult. (1618->244), div. (0->0), fcn. (1592->8), ass. (0->64)
t153 = MDP(20) + MDP(24);
t133 = cos(qJ(4));
t164 = -qJ(5) - pkin(8);
t112 = t164 * t133;
t127 = sin(pkin(10));
t129 = cos(pkin(10));
t131 = sin(qJ(4));
t149 = t164 * t131;
t101 = -t129 * t112 + t127 * t149;
t140 = MDP(17) * t131 + MDP(18) * t133;
t99 = -t127 * t112 - t129 * t149;
t169 = t131 * MDP(14) + t133 * MDP(15) - t99 * MDP(21) + t101 * MDP(23) - t140 * pkin(8);
t168 = 2 * MDP(19);
t167 = 0.2e1 * MDP(21);
t166 = 2 * MDP(22);
t165 = 0.2e1 * MDP(23);
t130 = cos(pkin(9));
t121 = -pkin(1) * t130 - pkin(2);
t132 = sin(qJ(3));
t134 = cos(qJ(3));
t108 = -pkin(3) * t134 - pkin(8) * t132 + t121;
t106 = t133 * t108;
t156 = t132 * t133;
t128 = sin(pkin(9));
t118 = pkin(1) * t128 + pkin(7);
t161 = t118 * t131;
t89 = -qJ(5) * t156 + t106 + (-pkin(4) - t161) * t134;
t159 = t118 * t134;
t151 = t133 * t159;
t92 = t151 + (-qJ(5) * t132 + t108) * t131;
t83 = t127 * t89 + t129 * t92;
t109 = t127 * t131 - t129 * t133;
t110 = t127 * t133 + t129 * t131;
t122 = -pkin(4) * t133 - pkin(3);
t94 = pkin(5) * t109 - qJ(6) * t110 + t122;
t163 = MDP(24) * t94;
t160 = t118 * t133;
t158 = t131 * t132;
t157 = t131 * t133;
t107 = pkin(4) * t158 + t132 * t118;
t155 = MDP(21) * t109;
t154 = MDP(23) * t110;
t152 = t101 ^ 2 + t99 ^ 2;
t150 = MDP(13) * t157;
t103 = t110 * t132;
t105 = -t127 * t158 + t129 * t156;
t147 = -t101 * t103 + t105 * t99;
t82 = -t127 * t92 + t129 * t89;
t142 = MDP(14) * t133 - MDP(15) * t131;
t141 = MDP(17) * t133 - MDP(18) * t131;
t139 = -t103 * MDP(21) + t105 * MDP(23);
t137 = MDP(20) * t122 - t154 + t155 + t163;
t126 = t134 ^ 2;
t125 = t133 ^ 2;
t124 = t132 ^ 2;
t123 = t131 ^ 2;
t119 = pkin(4) * t129 + pkin(5);
t116 = pkin(4) * t127 + qJ(6);
t97 = t108 * t131 + t151;
t96 = -t131 * t159 + t106;
t84 = pkin(5) * t103 - qJ(6) * t105 + t107;
t81 = pkin(5) * t134 - t82;
t80 = -qJ(6) * t134 + t83;
t1 = [MDP(1) - 0.2e1 * t121 * t134 * MDP(10) + t126 * MDP(16) + (t107 ^ 2 + t82 ^ 2 + t83 ^ 2) * MDP(20) + (t80 ^ 2 + t81 ^ 2 + t84 ^ 2) * MDP(24) + (t128 ^ 2 + t130 ^ 2) * MDP(4) * pkin(1) ^ 2 + (MDP(12) * t125 + MDP(5) - 0.2e1 * t150) * t124 + 0.2e1 * (t121 * MDP(11) + (MDP(6) - t142) * t134) * t132 + 0.2e1 * (t124 * t161 - t134 * t96) * MDP(17) + 0.2e1 * (t124 * t160 + t134 * t97) * MDP(18) + (-t103 * t83 - t105 * t82) * t168 + (t103 * t84 + t134 * t81) * t167 + (-t103 * t80 + t105 * t81) * t166 + (-t105 * t84 - t134 * t80) * t165; (-t103 * t82 + t105 * t83 - t107 * t134) * MDP(20) + (t103 * t81 + t105 * t80 - t134 * t84) * MDP(24); MDP(4) + t153 * (t103 ^ 2 + t105 ^ 2 + t126); (-t109 * t83 - t110 * t82 + t147) * MDP(19) + (t101 * t83 + t107 * t122 - t82 * t99) * MDP(20) + (t103 * t94 + t109 * t84) * MDP(21) + (-t109 * t80 + t110 * t81 + t147) * MDP(22) + (-t105 * t94 - t110 * t84) * MDP(23) + (t101 * t80 + t81 * t99 + t84 * t94) * MDP(24) + (-t118 * MDP(11) + MDP(8) - t169) * t134 + (MDP(7) - t118 * MDP(10) + MDP(12) * t157 + (-t123 + t125) * MDP(13) + (-pkin(3) * t131 - t160) * MDP(17) + (-pkin(3) * t133 + t161) * MDP(18)) * t132; -t132 * MDP(11) + (MDP(10) - t137 + t141) * t134 + t153 * (t105 * t101 + t103 * t99) + (MDP(19) + MDP(22)) * (t103 * t110 - t105 * t109); MDP(9) + t123 * MDP(12) + 0.2e1 * t150 + (t122 ^ 2 + t152) * MDP(20) + t152 * MDP(24) + (-0.2e1 * t154 + 0.2e1 * t155 + t163) * t94 + 0.2e1 * t141 * pkin(3) + (t168 + t166) * (-t101 * t109 + t110 * t99); t96 * MDP(17) - t97 * MDP(18) + t82 * MDP(21) + (-t103 * t116 - t105 * t119) * MDP(22) + t83 * MDP(23) + (t116 * t80 - t119 * t81) * MDP(24) + t142 * t132 + (-MDP(16) + (-pkin(5) - t119) * MDP(21) + (-qJ(6) - t116) * MDP(23)) * t134 + ((-t103 * t127 - t105 * t129) * MDP(19) + (t127 * t83 + t129 * t82) * MDP(20)) * pkin(4); (-t103 * t119 + t105 * t116) * MDP(24) - t140 * t132 + (-t103 * t129 + t105 * t127) * MDP(20) * pkin(4) + t139; (-t109 * t116 - t110 * t119) * MDP(22) + (t101 * t116 - t119 * t99) * MDP(24) + ((-t109 * t127 - t110 * t129) * MDP(19) + (t101 * t127 - t129 * t99) * MDP(20)) * pkin(4) + t169; MDP(16) + (t116 ^ 2 + t119 ^ 2) * MDP(24) + (t127 ^ 2 + t129 ^ 2) * MDP(20) * pkin(4) ^ 2 + t119 * t167 + t116 * t165; MDP(20) * t107 + MDP(24) * t84 - t139; -t153 * t134; t137; 0; t153; MDP(21) * t134 + MDP(22) * t105 + MDP(24) * t81; t103 * MDP(24); MDP(22) * t110 + MDP(24) * t99; -MDP(24) * t119 - MDP(21); 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

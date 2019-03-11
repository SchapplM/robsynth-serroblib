% Calculate joint inertia matrix for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR3_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:17
% EndTime: 2019-03-09 05:07:19
% DurationCPUTime: 0.70s
% Computational Cost: add. (525->159), mult. (943->222), div. (0->0), fcn. (841->8), ass. (0->72)
t126 = sin(qJ(6));
t129 = cos(qJ(6));
t131 = cos(qJ(3));
t128 = sin(qJ(3));
t130 = cos(qJ(4));
t165 = t128 * t130;
t119 = t131 * pkin(4);
t125 = cos(pkin(10));
t116 = -pkin(1) * t125 - pkin(2);
t101 = -pkin(3) * t131 - pkin(8) * t128 + t116;
t124 = sin(pkin(10));
t115 = pkin(1) * t124 + pkin(7);
t127 = sin(qJ(4));
t166 = t127 * t131;
t161 = -t130 * t101 + t115 * t166;
t90 = t119 + t161;
t87 = pkin(5) * t131 - pkin(9) * t165 + t90;
t168 = t127 * t128;
t169 = t115 * t130;
t93 = t127 * t101 + t131 * t169;
t89 = -qJ(5) * t131 + t93;
t88 = pkin(9) * t168 + t89;
t167 = t127 * t129;
t97 = t126 * t165 - t128 * t167;
t103 = t126 * t127 + t129 * t130;
t98 = t103 * t128;
t181 = -(t126 * t88 - t129 * t87) * MDP(28) - (t126 * t87 + t129 * t88) * MDP(29) + t98 * MDP(25) - t97 * MDP(26);
t180 = t161 * MDP(17) + t93 * MDP(18) + t181;
t170 = qJ(5) * t127;
t174 = pkin(4) + pkin(5);
t102 = t130 * t174 + pkin(3) + t170;
t177 = 0.2e1 * t102;
t175 = -2 * MDP(24);
t173 = pkin(8) - pkin(9);
t172 = pkin(4) * t127;
t171 = pkin(4) * MDP(22);
t120 = t127 ^ 2;
t122 = t130 ^ 2;
t160 = t120 + t122;
t159 = MDP(13) * t130;
t158 = qJ(5) * MDP(21);
t157 = t103 * MDP(28);
t104 = -t126 * t130 + t167;
t156 = t104 * MDP(23);
t155 = (qJ(5) * t126 + t129 * t174) * MDP(28);
t154 = (t129 * qJ(5) - t126 * t174) * MDP(29);
t153 = MDP(17) + MDP(19);
t152 = MDP(18) - MDP(21);
t151 = MDP(27) + MDP(16);
t150 = t173 * t127;
t149 = t160 * pkin(8);
t148 = -pkin(4) * t130 - t170;
t146 = t127 * t90 + t130 * t89;
t142 = t97 * MDP(28) + t98 * MDP(29);
t141 = t130 * MDP(14) - t127 * MDP(15);
t140 = t127 * MDP(19) - t130 * MDP(21);
t139 = t129 * MDP(28) - t126 * MDP(29);
t138 = (-t153 - t171) * t127;
t137 = MDP(19) + t139;
t136 = -MDP(27) - t154 - t155;
t111 = t173 * t130;
t135 = t104 * MDP(25) - t103 * MDP(26) - (t111 * t126 - t129 * t150) * MDP(28) - (t129 * t111 + t126 * t150) * MDP(29);
t134 = -t127 * MDP(14) + t135;
t123 = t131 ^ 2;
t121 = t128 ^ 2;
t117 = t122 * t128;
t114 = pkin(8) * t166;
t113 = qJ(5) * t165;
t110 = -pkin(3) + t148;
t96 = -t113 + (t115 + t172) * t128;
t91 = t113 + (-t127 * t174 - t115) * t128;
t1 = [(t89 ^ 2 + t90 ^ 2 + t96 ^ 2) * MDP(22) + MDP(1) + (t98 * MDP(23) + t175 * t97) * t98 + (t124 ^ 2 + t125 ^ 2) * MDP(4) * pkin(1) ^ 2 + t151 * t123 + 0.2e1 * t142 * t91 + 0.2e1 * (t116 * MDP(11) + (-t127 * t89 + t130 * t90) * MDP(20) + t140 * t96) * t128 + (t122 * MDP(12) - 0.2e1 * t127 * t159 + MDP(5) + 0.2e1 * (t127 * MDP(17) + t130 * MDP(18)) * t115) * t121 + 0.2e1 * (-t116 * MDP(10) + (MDP(6) - t141) * t128 + t90 * MDP(19) - t89 * MDP(21) + t180) * t131; (t128 * t146 - t131 * t96) * MDP(22); MDP(4) + (t121 * t160 + t123) * MDP(22); t117 * MDP(13) + t114 * MDP(17) + (-t130 * t96 + t114) * MDP(19) + t146 * MDP(20) - t127 * t96 * MDP(21) + (pkin(8) * t146 + t110 * t96) * MDP(22) + t98 * t156 + (-t103 * t98 - t104 * t97) * MDP(24) + (t102 * t97 + t103 * t91) * MDP(28) + (t102 * t98 + t104 * t91) * MDP(29) + (MDP(7) - t115 * MDP(10) + t130 * t127 * MDP(12) - t120 * MDP(13) + (-pkin(3) * t127 - t169) * MDP(17) + (-pkin(3) * t130 + t115 * t127) * MDP(18) + t140 * t110) * t128 + (-t115 * MDP(11) + MDP(8) + (pkin(8) * t152 - MDP(15)) * t130 + t134) * t131; t117 * MDP(20) + (t120 * MDP(20) + MDP(22) * t149 - MDP(11)) * t128 + (-t110 * MDP(22) + t104 * MDP(29) - t127 * t152 + t130 * t153 + MDP(10) + t157) * t131; MDP(9) + t120 * MDP(12) + (pkin(8) ^ 2 * t160 + t110 ^ 2) * MDP(22) + t157 * t177 + 0.2e1 * MDP(20) * t149 + (MDP(29) * t177 + t103 * t175 + t156) * t104 + 0.2e1 * (pkin(3) * MDP(17) - t110 * MDP(19)) * t130 + 0.2e1 * (-pkin(3) * MDP(18) - t110 * MDP(21) + t159) * t127; (-0.2e1 * t119 - t161) * MDP(19) + t93 * MDP(21) + (-pkin(4) * t90 + qJ(5) * t89) * MDP(22) + (-MDP(16) + t136 - 0.2e1 * t158) * t131 + (MDP(20) * t148 + t141) * t128 - t180; t113 * MDP(22) + (-t130 * t152 + t138) * t128 + t142; t130 * MDP(15) + (qJ(5) * t130 - t172) * MDP(20) + ((MDP(22) * qJ(5) - t152) * t130 + t138) * pkin(8) - t134; 0.2e1 * pkin(4) * MDP(19) + 0.2e1 * t158 + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(22) + 0.2e1 * t155 + 0.2e1 * t154 + t151; MDP(20) * t165 + MDP(22) * t90 + t131 * t137; MDP(22) * t168; (MDP(22) * pkin(8) + MDP(20)) * t127; -t137 - t171; MDP(22); t131 * MDP(27) + t181; -t142; t135; t136; t139; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

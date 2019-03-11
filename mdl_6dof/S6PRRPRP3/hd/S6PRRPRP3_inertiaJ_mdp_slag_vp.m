% Calculate joint inertia matrix for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP3_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:38:23
% EndTime: 2019-03-08 21:38:25
% DurationCPUTime: 0.60s
% Computational Cost: add. (818->181), mult. (1723->268), div. (0->0), fcn. (1852->10), ass. (0->78)
t158 = (MDP(15) * qJ(4));
t176 = MDP(14) + t158;
t175 = 2 * MDP(14);
t174 = -2 * MDP(17);
t173 = 2 * MDP(22);
t172 = 2 * MDP(23);
t171 = 2 * MDP(24);
t170 = 2 * MDP(25);
t169 = cos(qJ(5));
t132 = cos(qJ(3));
t168 = pkin(5) * t132;
t125 = sin(pkin(11));
t167 = pkin(8) * t125;
t166 = pkin(8) * t132;
t165 = pkin(9) + qJ(4);
t129 = sin(qJ(5));
t130 = sin(qJ(3));
t114 = -pkin(3) * t132 - qJ(4) * t130 - pkin(2);
t127 = cos(pkin(11));
t109 = t127 * t114;
t93 = -pkin(9) * t127 * t130 + t109 + (-pkin(4) - t167) * t132;
t102 = t125 * t114 + t127 * t166;
t162 = t125 * t130;
t98 = -pkin(9) * t162 + t102;
t85 = t129 * t93 + t169 * t98;
t164 = pkin(3) * MDP(15);
t163 = qJ(6) * t132;
t126 = sin(pkin(6));
t131 = sin(qJ(2));
t161 = t126 * t131;
t133 = cos(qJ(2));
t160 = t126 * t133;
t113 = pkin(4) * t162 + t130 * pkin(8);
t112 = t125 * t169 + t129 * t127;
t103 = t112 * t130;
t157 = t103 * MDP(19);
t147 = t169 * t127;
t104 = -t129 * t162 + t130 * t147;
t156 = t104 * MDP(18);
t155 = t112 * MDP(16);
t154 = t125 * MDP(13);
t153 = t127 * MDP(12);
t152 = t130 * MDP(11);
t151 = t132 * MDP(20);
t150 = MDP(15) + MDP(26);
t149 = MDP(21) + MDP(23);
t148 = MDP(22) - MDP(25);
t120 = -pkin(4) * t127 - pkin(3);
t146 = t165 * t125;
t84 = -t129 * t98 + t169 * t93;
t145 = -MDP(26) * pkin(5) - MDP(23);
t144 = t149 * t103;
t143 = -MDP(21) + t145;
t141 = MDP(26) * qJ(6) - t148;
t139 = MDP(12) * t125 + MDP(13) * t127;
t111 = t125 * t129 - t147;
t138 = t112 * MDP(18) - t111 * MDP(19);
t137 = -t153 + t154 - t164;
t136 = pkin(8) * MDP(15) + t139;
t91 = pkin(5) * t111 - qJ(6) * t112 + t120;
t135 = t91 * MDP(26) + t111 * t149 + t112 * t148 + t137;
t128 = cos(pkin(6));
t124 = t130 ^ 2;
t115 = t165 * t127;
t107 = t128 * t130 + t132 * t161;
t106 = -t128 * t132 + t130 * t161;
t105 = t106 ^ 2;
t101 = -t125 * t166 + t109;
t100 = t115 * t169 - t129 * t146;
t99 = t115 * t129 + t146 * t169;
t97 = t107 * t127 - t125 * t160;
t96 = -t107 * t125 - t127 * t160;
t88 = pkin(5) * t103 - qJ(6) * t104 + t113;
t87 = t129 * t96 + t169 * t97;
t86 = t129 * t97 - t169 * t96;
t83 = -t84 + t168;
t82 = -t163 + t85;
t1 = [MDP(1) + (t96 ^ 2 + t97 ^ 2 + t105) * MDP(15) + (t86 ^ 2 + t87 ^ 2 + t105) * MDP(26); (t101 * t96 + t102 * t97) * MDP(15) - t87 * t103 * MDP(24) + (t106 * t88 + t82 * t87 + t83 * t86) * MDP(26) + t106 * t144 + (t86 * MDP(24) + t106 * t148) * t104 + (-t96 * MDP(12) + t97 * MDP(13) + t148 * t87 + t149 * t86) * t132 + ((-t125 * t97 - t127 * t96) * MDP(14) + t136 * t106) * t130 + (-t131 * MDP(4) + (t132 * MDP(10) + MDP(3) - t152) * t133) * t126; MDP(2) + t124 * MDP(5) - 0.2e1 * pkin(2) * t152 + (pkin(8) ^ 2 * t124 + t101 ^ 2 + t102 ^ 2) * MDP(15) + (t82 ^ 2 + t83 ^ 2 + t88 ^ 2) * MDP(26) + (t104 * MDP(16) + t103 * t174) * t104 + (0.2e1 * pkin(2) * MDP(10) + 0.2e1 * t130 * MDP(6) + t151 - 0.2e1 * t156 + 0.2e1 * t157) * t132 + 0.2e1 * (-t101 * t132 + t124 * t167) * MDP(12) + 0.2e1 * (pkin(8) * t124 * t127 + t102 * t132) * MDP(13) + 0.2e1 * (t103 * t113 - t132 * t84) * MDP(21) + (t104 * t113 + t132 * t85) * t173 + (t103 * t88 + t132 * t83) * t172 + (-t103 * t82 + t104 * t83) * t171 + (-t104 * t88 - t132 * t82) * t170 + (-t101 * t127 - t102 * t125) * t130 * t175; -t107 * MDP(11) + (-t111 * t87 + t112 * t86) * MDP(24) + (t100 * t87 + t86 * t99) * MDP(26) + (-MDP(10) + t135) * t106 + t176 * (-t125 * t96 + t127 * t97); t104 * t155 + (-t103 * t112 - t104 * t111) * MDP(17) + (t103 * t120 + t111 * t113) * MDP(21) + (t104 * t120 + t112 * t113) * MDP(22) + (t103 * t91 + t111 * t88) * MDP(23) + (-t100 * t103 + t104 * t99 - t111 * t82 + t112 * t83) * MDP(24) + (-t104 * t91 - t112 * t88) * MDP(25) + (t100 * t82 + t83 * t99 + t88 * t91) * MDP(26) + (-pkin(8) * MDP(11) + qJ(4) * t139 + t100 * t148 + t149 * t99 + MDP(8) - t138) * t132 + (MDP(7) - t139 * pkin(3) + (-MDP(10) + t137) * pkin(8)) * t130 + t176 * (-t101 * t125 + t102 * t127); MDP(9) + (t100 ^ 2 + t91 ^ 2 + t99 ^ 2) * MDP(26) + (0.2e1 * t153 - 0.2e1 * t154 + t164) * pkin(3) + 0.2e1 * (t120 * MDP(21) + t91 * MDP(23) - t100 * MDP(24)) * t111 + (-0.2e1 * t91 * MDP(25) + t111 * t174 + t120 * t173 + t171 * t99 + t155) * t112 + (t175 + t158) * (t125 ^ 2 + t127 ^ 2) * qJ(4); t150 * t106; t88 * MDP(26) + t104 * t148 + t130 * t136 + t144; t135; t150; t141 * t87 + t143 * t86; t156 - t157 - t151 + t84 * MDP(21) - t85 * MDP(22) + (t84 - 0.2e1 * t168) * MDP(23) + (-pkin(5) * t104 - qJ(6) * t103) * MDP(24) + (-0.2e1 * t163 + t85) * MDP(25) + (-pkin(5) * t83 + qJ(6) * t82) * MDP(26); (-pkin(5) * t112 - qJ(6) * t111) * MDP(24) + t143 * t99 + t141 * t100 + t138; 0; MDP(20) + pkin(5) * t172 + qJ(6) * t170 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(26); t86 * MDP(26); t132 * MDP(23) + t104 * MDP(24) + t83 * MDP(26); t112 * MDP(24) + t99 * MDP(26); 0; t145; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

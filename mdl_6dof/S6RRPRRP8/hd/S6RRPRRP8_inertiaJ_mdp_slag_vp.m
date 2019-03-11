% Calculate joint inertia matrix for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP8_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:25:18
% EndTime: 2019-03-09 12:25:22
% DurationCPUTime: 1.16s
% Computational Cost: add. (1664->214), mult. (3251->302), div. (0->0), fcn. (3499->8), ass. (0->83)
t148 = sin(pkin(10));
t179 = pkin(8) + qJ(3);
t134 = t179 * t148;
t149 = cos(pkin(10));
t135 = t179 * t149;
t151 = sin(qJ(4));
t153 = cos(qJ(4));
t119 = -t134 * t153 - t135 * t151;
t120 = -t134 * t151 + t135 * t153;
t130 = t148 * t151 - t149 * t153;
t131 = t148 * t153 + t149 * t151;
t110 = -pkin(9) * t130 + t120;
t150 = sin(qJ(5));
t158 = -pkin(9) * t131 + t119;
t183 = cos(qJ(5));
t100 = t110 * t150 - t158 * t183;
t101 = t110 * t183 + t150 * t158;
t114 = t130 * t183 + t131 * t150;
t115 = -t150 * t130 + t131 * t183;
t165 = MDP(28) - MDP(31);
t166 = MDP(27) + MDP(29);
t162 = t115 * MDP(24) - t114 * MDP(25) - t100 * t166 - t101 * t165;
t198 = t131 * MDP(17) - t130 * MDP(18) + t119 * MDP(20) - t120 * MDP(21) + t162;
t159 = MDP(27) * t183 - t150 * MDP(28);
t195 = t159 * pkin(4);
t152 = sin(qJ(2));
t125 = t131 * t152;
t126 = t130 * t152;
t108 = t125 * t183 - t126 * t150;
t109 = -t150 * t125 - t126 * t183;
t174 = t109 * MDP(24) - t108 * MDP(25);
t154 = cos(qJ(2));
t133 = -pkin(2) * t154 - qJ(3) * t152 - pkin(1);
t129 = t149 * t133;
t181 = pkin(7) * t148;
t116 = -pkin(8) * t149 * t152 + t129 + (-pkin(3) - t181) * t154;
t180 = pkin(7) * t154;
t122 = t148 * t133 + t149 * t180;
t176 = t148 * t152;
t118 = -pkin(8) * t176 + t122;
t104 = t116 * t151 + t118 * t153;
t102 = -pkin(9) * t125 + t104;
t103 = t153 * t116 - t118 * t151;
t97 = -pkin(4) * t154 + pkin(9) * t126 + t103;
t88 = -t150 * t102 + t183 * t97;
t193 = t88 * MDP(27) + t174;
t192 = MDP(14) * qJ(3);
t191 = 2 * MDP(13);
t190 = -2 * MDP(16);
t189 = 0.2e1 * MDP(21);
t188 = -2 * MDP(23);
t187 = 0.2e1 * MDP(28);
t186 = 2 * MDP(29);
t185 = 2 * MDP(30);
t184 = 2 * MDP(31);
t182 = pkin(5) * t154;
t89 = t183 * t102 + t150 * t97;
t178 = pkin(2) * MDP(14);
t177 = qJ(6) * t154;
t132 = pkin(3) * t176 + t152 * pkin(7);
t172 = MDP(11) * t149;
t171 = MDP(12) * t148;
t170 = MDP(15) * t131;
t169 = MDP(20) * t130;
t168 = MDP(22) * t115;
t167 = MDP(19) + MDP(26);
t139 = -pkin(3) * t149 - pkin(2);
t164 = pkin(5) * t186 + MDP(26);
t117 = pkin(4) * t125 + t132;
t124 = pkin(4) * t130 + t139;
t161 = MDP(11) * t148 + MDP(12) * t149;
t160 = -t126 * MDP(17) - t125 * MDP(18);
t157 = t171 - t172 - t178;
t146 = t152 ^ 2;
t142 = t150 * pkin(4);
t140 = pkin(4) * t183 + pkin(5);
t138 = t142 + qJ(6);
t121 = -t148 * t180 + t129;
t99 = pkin(5) * t114 - qJ(6) * t115 + t124;
t90 = pkin(5) * t108 - qJ(6) * t109 + t117;
t87 = -t88 + t182;
t86 = -t177 + t89;
t1 = [t146 * MDP(4) + (pkin(7) ^ 2 * t146 + t121 ^ 2 + t122 ^ 2) * MDP(14) + (t86 ^ 2 + t87 ^ 2 + t90 ^ 2) * MDP(32) - 0.2e1 * pkin(1) * t152 * MDP(10) + MDP(1) + t167 * t154 ^ 2 - (-MDP(15) * t126 + t125 * t190) * t126 + (MDP(22) * t109 + t108 * t188) * t109 + 0.2e1 * (MDP(5) * t152 + MDP(9) * pkin(1) - t160 - t174) * t154 + (t117 * t109 + t89 * t154) * t187 + (t90 * t108 + t87 * t154) * t186 + (-t90 * t109 - t86 * t154) * t184 + 0.2e1 * (-t121 * t154 + t146 * t181) * MDP(11) + 0.2e1 * (t146 * pkin(7) * t149 + t122 * t154) * MDP(12) + 0.2e1 * (-t103 * t154 + t132 * t125) * MDP(20) + (t104 * t154 - t132 * t126) * t189 + 0.2e1 * (t117 * t108 - t88 * t154) * MDP(27) + (-t86 * t108 + t87 * t109) * t185 + (-t121 * t149 - t122 * t148) * t152 * t191; (-t125 * t131 + t126 * t130) * MDP(16) + (t125 * t139 + t130 * t132) * MDP(20) + (-t126 * t139 + t131 * t132) * MDP(21) + (-t108 * t115 - t109 * t114) * MDP(23) + (t108 * t124 + t114 * t117) * MDP(27) + (t109 * t124 + t115 * t117) * MDP(28) + (t108 * t99 + t114 * t90) * MDP(29) + (t100 * t109 - t101 * t108 - t114 * t86 + t115 * t87) * MDP(30) + (-t109 * t99 - t115 * t90) * MDP(31) + (t100 * t87 + t101 * t86 + t90 * t99) * MDP(32) - t126 * t170 + t109 * t168 + (MDP(6) - t161 * pkin(2) + (-MDP(9) + t157) * pkin(7)) * t152 + (-pkin(7) * MDP(10) + qJ(3) * t161 + MDP(7) - t198) * t154 + (MDP(13) + t192) * (-t121 * t148 + t122 * t149); MDP(8) + 0.2e1 * t139 * t169 + (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) * MDP(32) + (-0.2e1 * t171 + 0.2e1 * t172 + t178) * pkin(2) + 0.2e1 * (MDP(27) * t124 + MDP(29) * t99 - MDP(30) * t101) * t114 + (-0.2e1 * MDP(31) * t99 + t100 * t185 + t114 * t188 + t124 * t187 + t168) * t115 + (t130 * t190 + t139 * t189 + t170) * t131 + (t191 + t192) * (t148 ^ 2 + t149 ^ 2) * qJ(3); t125 * MDP(20) - t126 * MDP(21) + MDP(32) * t90 + t165 * t109 + t166 * t108 + (MDP(14) * pkin(7) + t161) * t152; MDP(21) * t131 + MDP(32) * t99 + t114 * t166 + t115 * t165 + t157 + t169; MDP(14) + MDP(32); t103 * MDP(20) - t104 * MDP(21) + t88 * MDP(29) + (-t108 * t138 - t109 * t140) * MDP(30) + (t138 * t86 - t140 * t87) * MDP(32) + ((-pkin(5) - t140) * MDP(29) + (-qJ(6) - t138) * MDP(31) - t195 - t167) * t154 + t160 - t165 * t89 + t193; (-t114 * t138 - t115 * t140) * MDP(30) + (-t100 * t140 + t101 * t138) * MDP(32) + t198; 0; (t138 ^ 2 + t140 ^ 2) * MDP(32) + 0.2e1 * t195 + t140 * t186 + t138 * t184 + t167; -t154 * MDP(26) - t89 * MDP(28) + (t88 - 0.2e1 * t182) * MDP(29) + (-pkin(5) * t109 - qJ(6) * t108) * MDP(30) + (-0.2e1 * t177 + t89) * MDP(31) + (-pkin(5) * t87 + qJ(6) * t86) * MDP(32) + t193; (-pkin(5) * t115 - qJ(6) * t114) * MDP(30) + (-pkin(5) * t100 + qJ(6) * t101) * MDP(32) + t162; 0; (0.2e1 * qJ(6) + t142) * MDP(31) + (pkin(5) * t140 + qJ(6) * t138) * MDP(32) + (MDP(29) * t183 + t159) * pkin(4) + t164; qJ(6) * t184 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32) + t164; t154 * MDP(29) + t109 * MDP(30) + MDP(32) * t87; MDP(30) * t115 + MDP(32) * t100; 0; -MDP(32) * t140 - MDP(29); -MDP(32) * pkin(5) - MDP(29); MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

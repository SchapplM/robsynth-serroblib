% Calculate joint inertia matrix for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR14_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR14_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:36:59
% EndTime: 2019-03-09 11:37:03
% DurationCPUTime: 0.99s
% Computational Cost: add. (870->222), mult. (1846->294), div. (0->0), fcn. (1806->8), ass. (0->92)
t156 = sin(qJ(6));
t159 = cos(qJ(6));
t172 = t159 * MDP(31) - t156 * MDP(32);
t211 = MDP(22) + t172;
t160 = cos(qJ(4));
t210 = 0.2e1 * t160;
t209 = -2 * MDP(23);
t208 = 0.2e1 * MDP(31);
t207 = 0.2e1 * MDP(32);
t163 = -pkin(2) - pkin(9);
t206 = pkin(4) + pkin(10);
t161 = cos(qJ(2));
t205 = pkin(1) * t161;
t204 = pkin(5) - t163;
t203 = (MDP(25) * pkin(4));
t202 = pkin(2) * MDP(14);
t155 = cos(pkin(6));
t157 = sin(qJ(4));
t154 = sin(pkin(6));
t198 = t154 * t161;
t130 = t155 * t160 - t157 * t198;
t201 = t130 * t160;
t136 = t204 * t157;
t200 = t136 * t157;
t158 = sin(qJ(2));
t199 = t154 * t158;
t141 = pkin(8) * t199;
t183 = -pkin(2) - t205;
t117 = pkin(3) * t199 + t141 + (-pkin(9) + t183) * t155;
t181 = -qJ(3) * t158 - pkin(1);
t123 = (t163 * t161 + t181) * t154;
t112 = t157 * t117 + t160 * t123;
t132 = t155 * t158 * pkin(1) + pkin(8) * t198;
t149 = t157 ^ 2;
t152 = t160 ^ 2;
t140 = t149 + t152;
t197 = MDP(31) * t156;
t196 = MDP(32) * t159;
t129 = t155 * t157 + t160 * t198;
t121 = t129 * t156 + t159 * t199;
t195 = t121 * MDP(26);
t126 = t183 * t155 + t141;
t194 = t126 * MDP(14);
t193 = t129 * MDP(16);
t192 = t163 ^ 2 * MDP(25);
t191 = t156 * MDP(26);
t188 = t163 * MDP(25);
t187 = MDP(15) + MDP(30);
t186 = -MDP(21) + MDP(24);
t185 = pkin(4) * t199;
t184 = qJ(5) * t199;
t145 = t155 * qJ(3);
t124 = -t145 - t132;
t182 = t159 * t156 * MDP(27);
t180 = MDP(23) - t203;
t111 = t117 * t160 - t157 * t123;
t179 = (MDP(20) - MDP(23)) * t160;
t122 = pkin(3) * t198 - t124;
t135 = t157 * pkin(4) - qJ(5) * t160 + qJ(3);
t178 = -t132 * MDP(10) + (t155 * t205 - t141) * MDP(9);
t177 = t111 * MDP(20) - t112 * MDP(21);
t176 = t129 * MDP(20) + t130 * MDP(21);
t175 = -t129 * MDP(23) - t130 * MDP(24);
t120 = -t129 * t159 + t156 * t199;
t174 = t121 * MDP(28) - t120 * MDP(29);
t173 = MDP(28) * t156 + MDP(29) * t159;
t171 = t196 + t197;
t109 = -t184 - t112;
t170 = -MDP(16) + t173;
t169 = -qJ(5) * t130 + t122;
t168 = -t172 + t188;
t167 = t186 * t157 + t179;
t106 = pkin(5) * t130 - t206 * t199 - t111;
t108 = t206 * t129 + t169;
t104 = t106 * t159 - t108 * t156;
t105 = t106 * t156 + t108 * t159;
t166 = t104 * MDP(31) - t105 * MDP(32) + t174;
t165 = (-MDP(31) * t206 + MDP(28)) * t159 + (MDP(32) * t206 - MDP(29)) * t156;
t164 = MDP(17) + t165;
t151 = t159 ^ 2;
t148 = t156 ^ 2;
t138 = pkin(4) * t160 + qJ(5) * t157;
t137 = t204 * t160;
t134 = t140 * t163;
t133 = pkin(10) * t157 + t135;
t125 = (-pkin(2) * t161 + t181) * t154;
t116 = t133 * t159 + t137 * t156;
t115 = -t133 * t156 + t137 * t159;
t113 = pkin(4) * t129 + t169;
t110 = -t111 - t185;
t107 = -pkin(5) * t129 - t109;
t1 = [(t109 ^ 2 + t110 ^ 2 + t113 ^ 2) * MDP(25) + (t124 ^ 2 + t125 ^ 2 + t126 ^ 2) * MDP(14) + t155 ^ 2 * MDP(8) + MDP(1) + t187 * t130 ^ 2 + (-0.2e1 * t120 * MDP(27) + t195) * t121 + 0.2e1 * (MDP(17) * t199 + t174 - t193) * t130 + (-t105 * t130 + t107 * t121) * t207 + 0.2e1 * (t109 * t129 + t110 * t130) * MDP(22) + (t104 * t130 + t107 * t120) * t208 + 0.2e1 * (t126 * MDP(12) - t124 * MDP(13) + t178) * t155 + 0.2e1 * t176 * t122 + 0.2e1 * t175 * t113 + ((0.2e1 * t161 * MDP(5) + (MDP(19) + MDP(4)) * t158) * t158 + 0.2e1 * (-MDP(10) * t158 + MDP(9) * t161) * pkin(1)) * t154 ^ 2 + 0.2e1 * ((-t124 * MDP(11) + t125 * MDP(12) + MDP(7) * t155) * t161 + (t126 * MDP(11) - t125 * MDP(13) - t129 * MDP(18) + t110 * MDP(23) - t109 * MDP(24) + MDP(6) * t155 + t177) * t158) * t154; -pkin(2) * t194 + (0.2e1 * t145 + t132) * MDP(13) + (-t116 * t130 - t121 * t136) * MDP(32) + (t115 * t130 - t120 * t136) * MDP(31) + t141 * MDP(12) + (MDP(8) + (-0.2e1 * pkin(2) - t205) * MDP(12)) * t155 + (-t124 * MDP(14) + t176) * qJ(3) + (t113 * MDP(25) + t175) * t135 + (-t193 + t122 * MDP(21) - t113 * MDP(24) + (MDP(22) - t188) * t110 + (-t163 * MDP(22) + t187) * t130 + t166) * t160 + ((-t120 * t156 + t121 * t159) * MDP(27) + (-t129 * t163 + t109) * MDP(22) - t109 * t188 + t122 * MDP(20) - t113 * MDP(23) + t121 * t191 - t172 * t107 + t170 * t130) * t157 + ((qJ(3) * MDP(11) + MDP(7)) * t161 + (-pkin(2) * MDP(11) + t160 * MDP(17) - t157 * MDP(18) + t167 * t163 + MDP(6)) * t158) * t154 + t178; MDP(8) + (-0.2e1 * t160 * MDP(24) + MDP(25) * t135) * t135 + (t187 + t192) * t152 + (t148 * MDP(26) + 0.2e1 * t182 + t192) * t149 + (-0.2e1 * MDP(12) + t202) * pkin(2) + (MDP(14) * qJ(3) + MDP(21) * t210 + 0.2e1 * MDP(13)) * qJ(3) + (0.2e1 * qJ(3) * MDP(20) + t135 * t209 + t170 * t210) * t157 - 0.2e1 * t134 * MDP(22) + (t115 * t160 + t159 * t200) * t208 + (-t116 * t160 - t156 * t200) * t207; t155 * MDP(12) + t194 + (-t129 * t157 - t201) * MDP(22) + (-t109 * t157 - t110 * t160) * MDP(25) + (t120 * t157 - t159 * t201) * MDP(31) + (t121 * t157 + t156 * t201) * MDP(32) + (MDP(11) + t167) * t199; t134 * MDP(25) - t140 * t211 + MDP(12) - t202; MDP(25) * t140 + MDP(14); MDP(19) * t199 + (-t111 - 0.2e1 * t185) * MDP(23) + (0.2e1 * t184 + t112) * MDP(24) + (-pkin(4) * t110 - qJ(5) * t109) * MDP(25) + t159 * t195 + (-t120 * t159 - t121 * t156) * MDP(27) + (qJ(5) * t120 + t107 * t156) * MDP(31) + (qJ(5) * t121 + t107 * t159) * MDP(32) + (-MDP(22) * qJ(5) - MDP(18)) * t129 + (-pkin(4) * MDP(22) + t164) * t130 + t177; -t138 * MDP(22) - t171 * t136 + ((MDP(20) - t180) * t163 + t164) * t160 + (-MDP(18) + t159 * t191 + (-t148 + t151) * MDP(27) + t186 * t163 + t168 * qJ(5)) * t157; MDP(25) * t138 + t179 + (t171 + t186) * t157; -0.2e1 * t182 + t151 * MDP(26) + MDP(19) + (t209 + t203) * pkin(4) + (MDP(25) * qJ(5) + 0.2e1 * MDP(24) + 0.2e1 * t196 + 0.2e1 * t197) * qJ(5); MDP(23) * t199 + t110 * MDP(25) + t130 * t211; (MDP(22) - t168) * t160; -t160 * MDP(25); t180; MDP(25); t130 * MDP(30) + t166; t160 * MDP(30) + t115 * MDP(31) - t116 * MDP(32) + t173 * t157; -t172 * t160; t165; t172; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

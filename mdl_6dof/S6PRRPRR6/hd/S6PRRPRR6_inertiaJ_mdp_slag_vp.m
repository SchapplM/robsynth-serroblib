% Calculate joint inertia matrix for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR6_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:27:41
% EndTime: 2019-03-08 22:27:44
% DurationCPUTime: 1.00s
% Computational Cost: add. (1250->217), mult. (3040->337), div. (0->0), fcn. (3543->14), ass. (0->106)
t163 = sin(pkin(7));
t224 = 0.2e1 * t163;
t223 = (MDP(15) * qJ(4));
t222 = 2 * MDP(14);
t221 = 2 * MDP(21);
t220 = 2 * MDP(28);
t219 = 2 * MDP(29);
t218 = cos(qJ(5));
t170 = sin(qJ(3));
t217 = pkin(2) * t170;
t173 = cos(qJ(3));
t216 = pkin(2) * t173;
t215 = pkin(10) + qJ(4);
t214 = pkin(3) * MDP(15);
t162 = sin(pkin(13));
t165 = cos(pkin(13));
t169 = sin(qJ(5));
t149 = t218 * t162 + t169 * t165;
t168 = sin(qJ(6));
t213 = t149 * t168;
t172 = cos(qJ(6));
t212 = t149 * t172;
t211 = t163 * t170;
t210 = t163 * t173;
t166 = cos(pkin(7));
t209 = t166 * MDP(9);
t174 = cos(qJ(2));
t208 = t166 * t174;
t207 = t170 * MDP(7);
t190 = pkin(9) * t210;
t137 = t190 + (qJ(4) + t217) * t166;
t138 = (-pkin(3) * t173 - qJ(4) * t170 - pkin(2)) * t163;
t125 = t165 * t137 + t162 * t138;
t205 = MDP(20) * t173;
t204 = MDP(22) * t149;
t141 = t162 * t211 - t165 * t166;
t143 = t162 * t166 + t165 * t211;
t127 = -t169 * t141 + t218 * t143;
t122 = t127 * t168 + t172 * t210;
t203 = t122 * MDP(26);
t123 = t127 * t172 - t168 * t210;
t202 = t123 * MDP(23);
t201 = t123 * MDP(25);
t126 = t218 * t141 + t143 * t169;
t200 = t126 * MDP(27);
t199 = t127 * MDP(17);
t198 = t127 * MDP(18);
t197 = t127 * MDP(22);
t151 = t215 * t165;
t187 = t215 * t162;
t133 = t218 * t151 - t169 * t187;
t196 = t133 * MDP(22);
t154 = pkin(9) * t211;
t140 = t154 + (-pkin(3) - t216) * t166;
t195 = t140 * MDP(15);
t148 = t162 * t169 - t218 * t165;
t194 = t148 * MDP(27);
t193 = t162 * MDP(13);
t192 = t165 * MDP(12);
t191 = t172 * MDP(23);
t156 = -pkin(4) * t165 - pkin(3);
t189 = MDP(24) * t168 * t172;
t188 = MDP(19) * t210;
t124 = -t137 * t162 + t165 * t138;
t185 = -t124 * t162 + t125 * t165;
t184 = t172 * MDP(25) - t168 * MDP(26);
t183 = MDP(28) * t172 - MDP(29) * t168;
t182 = t168 * MDP(28) + t172 * MDP(29);
t115 = -pkin(4) * t210 - pkin(10) * t143 + t124;
t118 = -pkin(10) * t141 + t125;
t110 = t218 * t115 - t169 * t118;
t111 = t169 * t115 + t218 * t118;
t181 = -MDP(17) + t184;
t180 = MDP(21) + t183;
t128 = pkin(4) * t141 + t140;
t179 = -t192 + t193 + t204 - t214;
t178 = t168 * MDP(25) + t172 * MDP(26) - t182 * pkin(11);
t109 = -pkin(11) * t210 + t111;
t112 = pkin(5) * t126 - pkin(11) * t127 + t128;
t104 = -t109 * t168 + t112 * t172;
t105 = t109 * t172 + t112 * t168;
t177 = t104 * MDP(28) - t105 * MDP(29) + t200 + t201 - t203;
t176 = -MDP(19) + t178;
t171 = sin(qJ(2));
t167 = cos(pkin(6));
t164 = sin(pkin(6));
t161 = t172 ^ 2;
t160 = t168 ^ 2;
t158 = t163 ^ 2;
t145 = t166 * t217 + t190;
t144 = t166 * t216 - t154;
t142 = -t163 * t164 * t174 + t166 * t167;
t132 = t151 * t169 + t218 * t187;
t131 = t167 * t211 + (t170 * t208 + t171 * t173) * t164;
t130 = -t167 * t210 + (t170 * t171 - t173 * t208) * t164;
t129 = pkin(5) * t148 - pkin(11) * t149 + t156;
t121 = t131 * t165 + t142 * t162;
t120 = -t131 * t162 + t142 * t165;
t117 = t129 * t168 + t133 * t172;
t116 = t129 * t172 - t133 * t168;
t114 = t169 * t120 + t218 * t121;
t113 = -t218 * t120 + t121 * t169;
t108 = pkin(5) * t210 - t110;
t107 = t114 * t172 + t130 * t168;
t106 = -t114 * t168 + t130 * t172;
t1 = [MDP(1) + (t120 ^ 2 + t121 ^ 2 + t130 ^ 2) * MDP(15); (-t130 * t166 - t142 * t210) * MDP(10) + (-t131 * t166 + t142 * t211) * MDP(11) + (-t120 * t210 + t130 * t141) * MDP(12) + (t121 * t210 + t130 * t143) * MDP(13) + (-t120 * t143 - t121 * t141) * MDP(14) + (t120 * t124 + t121 * t125 + t130 * t140) * MDP(15) + (t113 * t210 + t126 * t130) * MDP(21) + (t114 * t210 + t127 * t130) * MDP(22) + (t106 * t126 + t113 * t122) * MDP(28) + (-t107 * t126 + t113 * t123) * MDP(29) + (MDP(3) * t174 - MDP(4) * t171) * t164; (t124 ^ 2 + t125 ^ 2 + t140 ^ 2) * MDP(15) + t127 ^ 2 * MDP(16) + t158 * t170 ^ 2 * MDP(5) + MDP(2) + (t207 * t224 + t209) * t166 + (-0.2e1 * t122 * MDP(24) + t202) * t123 + ((MDP(8) * t166 - t198) * t224 + (0.2e1 * MDP(6) * t170 + t205) * t158) * t173 + (0.2e1 * t188 - 0.2e1 * t199 + t200 + 0.2e1 * t201 - 0.2e1 * t203) * t126 + (-t110 * t210 + t126 * t128) * t221 + 0.2e1 * (t111 * t210 + t127 * t128) * MDP(22) + 0.2e1 * (t144 * t166 + t158 * t216) * MDP(10) + 0.2e1 * (-t145 * t166 - t158 * t217) * MDP(11) + 0.2e1 * (-t124 * t210 + t140 * t141) * MDP(12) + 0.2e1 * (t125 * t210 + t140 * t143) * MDP(13) + (-t124 * t143 - t125 * t141) * t222 + (t104 * t126 + t108 * t122) * t220 + (-t105 * t126 + t108 * t123) * t219; -t131 * MDP(11) + (t106 * t148 + t113 * t213) * MDP(28) + (-t107 * t148 + t113 * t212) * MDP(29) + (MDP(21) * t148 - MDP(10) + t179) * t130 + (MDP(14) + t223) * (-t120 * t162 + t121 * t165); t209 + t144 * MDP(10) - t145 * MDP(11) + (-pkin(3) * t141 - t140 * t165) * MDP(12) + (-pkin(3) * t143 + t140 * t162) * MDP(13) + t185 * MDP(14) - pkin(3) * t195 + (t116 * t126 + t122 * t132) * MDP(28) + (-t117 * t126 + t123 * t132) * MDP(29) + (t126 * MDP(21) + t197) * t156 + (t207 + (t132 * MDP(21) + MDP(8) + t196) * t173) * t163 + ((-t141 * t165 + t143 * t162) * MDP(14) + t185 * MDP(15) + (MDP(12) * t162 + MDP(13) * t165) * t210) * qJ(4) + (t128 * MDP(21) + t177 + t188 - t199) * t148 + (t127 * MDP(16) - MDP(18) * t210 + t128 * MDP(22) + t123 * t191 + (-t122 * t172 - t123 * t168) * MDP(24) + t182 * t108 + t181 * t126) * t149; 0.2e1 * t156 * t204 + MDP(9) + (0.2e1 * t192 - 0.2e1 * t193 + t214) * pkin(3) + (MDP(23) * t161 + MDP(16) - 0.2e1 * t189) * t149 ^ 2 + (0.2e1 * t181 * t149 + t156 * t221 + t194) * t148 + (t116 * t148 + t132 * t213) * t220 + (-t117 * t148 + t132 * t212) * t219 + (t222 + t223) * (t162 ^ 2 + t165 ^ 2) * qJ(4); t130 * MDP(15); t141 * MDP(12) + t143 * MDP(13) + t126 * t180 + t195 + t197; t148 * t180 + t179; MDP(15); -MDP(22) * t114 - t113 * t180; t198 - t163 * t205 + t110 * MDP(21) - t111 * MDP(22) + t168 * t202 + (-t122 * t168 + t123 * t172) * MDP(24) + (-pkin(5) * t122 - t108 * t172) * MDP(28) + (-pkin(5) * t123 + t108 * t168) * MDP(29) + t176 * t126; -t196 - t180 * t132 + t176 * t148 + (MDP(18) + t168 * t191 + (-t160 + t161) * MDP(24) - t182 * pkin(5)) * t149; 0; MDP(23) * t160 + 0.2e1 * pkin(5) * t183 + MDP(20) + 0.2e1 * t189; MDP(28) * t106 - MDP(29) * t107; t177; MDP(28) * t116 - MDP(29) * t117 + t184 * t149 + t194; t183; t178; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

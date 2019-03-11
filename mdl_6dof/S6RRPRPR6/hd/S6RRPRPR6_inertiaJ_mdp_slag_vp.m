% Calculate joint inertia matrix for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR6_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:42:58
% EndTime: 2019-03-09 10:43:01
% DurationCPUTime: 0.94s
% Computational Cost: add. (1298->213), mult. (3070->303), div. (0->0), fcn. (3310->10), ass. (0->102)
t226 = 2 * MDP(19);
t225 = 2 * MDP(21);
t224 = 2 * MDP(22);
t223 = 2 * MDP(29);
t222 = 2 * MDP(30);
t221 = pkin(4) + pkin(10);
t174 = cos(qJ(2));
t220 = pkin(1) * t174;
t165 = sin(pkin(11));
t167 = cos(pkin(11));
t166 = sin(pkin(6));
t213 = t166 * t174;
t171 = sin(qJ(2));
t214 = t166 * t171;
t145 = t165 * t214 - t167 * t213;
t219 = pkin(4) * t145;
t158 = pkin(2) * t165 + pkin(9);
t218 = pkin(5) + t158;
t217 = pkin(8) + qJ(3);
t216 = MDP(23) * pkin(4);
t144 = t145 * qJ(5);
t160 = t166 ^ 2;
t215 = t160 * t171;
t168 = cos(pkin(6));
t212 = t168 * MDP(8);
t169 = sin(qJ(6));
t173 = cos(qJ(4));
t211 = t169 * t173;
t172 = cos(qJ(6));
t210 = t172 * t173;
t155 = t168 * t220;
t140 = pkin(2) * t168 - t214 * t217 + t155;
t193 = pkin(1) * t168 * t171;
t143 = t213 * t217 + t193;
t128 = t140 * t165 + t143 * t167;
t125 = pkin(9) * t168 + t128;
t146 = (t165 * t174 + t167 * t171) * t166;
t153 = (-pkin(2) * t174 - pkin(1)) * t166;
t132 = pkin(3) * t145 - pkin(9) * t146 + t153;
t170 = sin(qJ(4));
t121 = t125 * t173 + t132 * t170;
t162 = t170 ^ 2;
t164 = t173 ^ 2;
t209 = t162 + t164;
t208 = MDP(23) * qJ(5);
t207 = MDP(23) * t158 ^ 2;
t206 = MDP(29) * t169;
t205 = MDP(30) * t172;
t138 = t146 * t170 - t168 * t173;
t130 = -t138 * t172 + t145 * t169;
t204 = t130 * MDP(27);
t131 = t138 * t169 + t145 * t172;
t203 = t131 * MDP(24);
t202 = t131 * MDP(26);
t201 = t138 * MDP(14);
t200 = t145 * MDP(17);
t199 = t158 * MDP(20);
t198 = t158 * MDP(23);
t197 = t169 * MDP(24);
t196 = MDP(13) + MDP(28);
t195 = MDP(18) - MDP(21);
t194 = -MDP(19) + MDP(22);
t192 = t144 + t121;
t159 = -pkin(2) * t167 - pkin(3);
t191 = t172 * t169 * MDP(25);
t190 = MDP(21) - t216;
t120 = -t125 * t170 + t132 * t173;
t127 = t140 * t167 - t143 * t165;
t189 = MDP(18) - t190;
t188 = t194 + t208;
t187 = -MDP(26) * t169 - MDP(27) * t172;
t186 = MDP(29) * t172 - MDP(30) * t169;
t185 = t205 + t206;
t184 = -qJ(5) * t170 + t159;
t124 = -pkin(3) * t168 - t127;
t183 = MDP(14) + t187;
t182 = MDP(20) + t186;
t181 = (MDP(6) * t171 + MDP(7) * t174) * t166;
t139 = t146 * t173 + t168 * t170;
t114 = pkin(5) * t139 - t145 * t221 - t120;
t179 = -qJ(5) * t139 + t124;
t116 = t138 * t221 + t179;
t112 = t114 * t172 - t116 * t169;
t113 = t114 * t169 + t116 * t172;
t180 = t112 * MDP(29) - t113 * MDP(30) - t204;
t178 = (-MDP(29) * t221 + MDP(26)) * t172 + (MDP(30) * t221 - MDP(27)) * t169;
t177 = -pkin(4) * MDP(20) + MDP(15) + t178;
t163 = t172 ^ 2;
t161 = t169 ^ 2;
t152 = t218 * t173;
t151 = t218 * t170;
t150 = -pkin(4) * t173 + t184;
t149 = pkin(8) * t213 + t193;
t148 = -pkin(8) * t214 + t155;
t147 = -t173 * t221 + t184;
t134 = t147 * t172 + t151 * t169;
t133 = -t147 * t169 + t151 * t172;
t126 = t131 * t170;
t119 = pkin(4) * t138 + t179;
t118 = -t120 - t219;
t115 = -pkin(5) * t138 + t192;
t1 = [(t127 ^ 2 + t128 ^ 2 + t153 ^ 2) * MDP(12) + (t118 ^ 2 + t119 ^ 2 + t192 ^ 2) * MDP(23) + MDP(1) + (MDP(4) * t171 + 0.2e1 * MDP(5) * t174) * t215 + (-0.2e1 * t138 * MDP(16) + t200) * t145 + t196 * t139 ^ 2 + (-0.2e1 * MDP(25) * t130 + t203) * t131 + t212 * t168 + 0.2e1 * (t148 * t168 + t160 * t220) * MDP(9) + 0.2e1 * (-pkin(1) * t215 - t149 * t168) * MDP(10) + (t118 * t145 - t119 * t138) * t225 + 0.2e1 * (t120 * t145 + t124 * t138) * MDP(18) + (-t121 * t145 + t124 * t139) * t226 + (-t119 * t139 + t145 * t192) * t224 + 0.2e1 * (-t127 * t146 - t128 * t145) * MDP(11) + 0.2e1 * (t118 * t139 - t138 * t192) * MDP(20) + (t112 * t139 + t115 * t130) * t223 + (-t113 * t139 + t115 * t131) * t222 + 0.2e1 * t181 * t168 + 0.2e1 * (t145 * MDP(15) - t201 + t202 - t204) * t139; t126 * MDP(26) + (t130 * t152 + t133 * t139) * MDP(29) + (t131 * t152 - t134 * t139) * MDP(30) + t212 + t148 * MDP(9) - t149 * MDP(10) + t181 + (MDP(18) * t138 + MDP(19) * t139) * t159 + (-MDP(21) * t138 - MDP(22) * t139 + MDP(23) * t119) * t150 + ((-t145 * t165 - t146 * t167) * MDP(11) + (t127 * t167 + t128 * t165) * MDP(12)) * pkin(2) + (-t201 + t124 * MDP(19) - t119 * MDP(22) + (MDP(20) + t198) * t118 + (t196 + t199) * t139 + (-t158 * t195 + MDP(15)) * t145 + t180) * t170 + (-t124 * MDP(18) + (-t138 * t158 + t192) * MDP(20) + t119 * MDP(21) + t192 * t198 + (t130 * t169 - t131 * t172) * MDP(25) - t131 * t197 + t186 * t115 + (t158 * t194 + MDP(16)) * t145 + t183 * t139) * t173; -0.2e1 * t159 * t173 * MDP(18) + MDP(8) + (t165 ^ 2 + t167 ^ 2) * MDP(12) * pkin(2) ^ 2 + (MDP(23) * t150 + t173 * t225) * t150 + (MDP(24) * t161 + 0.2e1 * t191 + t207) * t164 + (t196 + t207) * t162 + (-0.2e1 * t150 * MDP(22) + t159 * t226 + 0.2e1 * t183 * t173) * t170 + (t133 * t170 + t152 * t210) * t223 + (-t134 * t170 - t152 * t211) * t222 + 0.2e1 * t209 * t199; t153 * MDP(12) + (-t138 * t170 - t139 * t173) * MDP(20) + (-t118 * t173 + t170 * t192) * MDP(23) + (t130 * t170 - t139 * t210) * MDP(29) + (t139 * t211 + t126) * MDP(30) + (t170 * t194 + t173 * t195) * t145; 0; MDP(23) * t209 + MDP(12); t200 + t120 * MDP(18) - t121 * MDP(19) + (-t120 - 0.2e1 * t219) * MDP(21) + (t192 + t144) * MDP(22) + (-pkin(4) * t118 + qJ(5) * t192) * MDP(23) + t172 * t203 + (-t130 * t172 - t131 * t169) * MDP(25) + (qJ(5) * t130 + t115 * t169) * MDP(29) + (qJ(5) * t131 + t115 * t172) * MDP(30) + (-qJ(5) * MDP(20) - MDP(16)) * t138 + t177 * t139; t185 * t152 + t177 * t170 + (MDP(16) - t172 * t197 + (t161 - t163) * MDP(25) + t182 * qJ(5)) * t173 + (-t170 * t189 + t173 * t188) * t158; t189 * t173 + (t185 + t188) * t170; -0.2e1 * t191 + t163 * MDP(24) + MDP(17) + (-(2 * MDP(21)) + t216) * pkin(4) + (t224 + 0.2e1 * t205 + 0.2e1 * t206 + t208) * qJ(5); MDP(21) * t145 + MDP(23) * t118 + t139 * t182; (t182 + t198) * t170; -t173 * MDP(23); t190; MDP(23); t139 * MDP(28) + t180 + t202; t170 * MDP(28) + t133 * MDP(29) - t134 * MDP(30) + t173 * t187; -t186 * t173; t178; t186; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

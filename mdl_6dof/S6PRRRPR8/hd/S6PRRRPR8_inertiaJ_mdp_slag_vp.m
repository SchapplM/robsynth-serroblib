% Calculate joint inertia matrix for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR8_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:53:24
% EndTime: 2019-03-08 23:53:27
% DurationCPUTime: 0.99s
% Computational Cost: add. (882->222), mult. (2143->323), div. (0->0), fcn. (2335->12), ass. (0->101)
t162 = cos(pkin(7));
t165 = sin(qJ(4));
t169 = cos(qJ(4));
t160 = sin(pkin(7));
t166 = sin(qJ(3));
t212 = t160 * t166;
t143 = t162 * t165 + t169 * t212;
t225 = 0.2e1 * t143;
t170 = cos(qJ(3));
t211 = t160 * t170;
t194 = MDP(17) - MDP(20);
t193 = MDP(18) - MDP(21);
t224 = 2 * MDP(17);
t223 = 2 * MDP(19);
t222 = 2 * MDP(20);
t221 = 2 * MDP(21);
t220 = 2 * MDP(28);
t219 = 2 * MDP(29);
t218 = pkin(4) + pkin(11);
t217 = pkin(5) + pkin(10);
t216 = pkin(2) * t166;
t215 = pkin(2) * t170;
t214 = (MDP(22) * pkin(4));
t213 = (MDP(22) * pkin(10));
t210 = t162 * MDP(9);
t171 = cos(qJ(2));
t209 = t162 * t171;
t164 = sin(qJ(6));
t208 = t164 * t169;
t168 = cos(qJ(6));
t207 = t168 * t169;
t192 = pkin(9) * t211;
t138 = t192 + (pkin(10) + t216) * t162;
t139 = (-pkin(3) * t170 - pkin(10) * t166 - pkin(2)) * t160;
t126 = t169 * t138 + t165 * t139;
t206 = MDP(16) * t170;
t205 = MDP(22) * qJ(5);
t204 = MDP(22) * pkin(10) ^ 2;
t142 = -t169 * t162 + t165 * t212;
t132 = t142 * t164 - t168 * t211;
t203 = MDP(23) * t132;
t202 = MDP(23) * t168;
t201 = MDP(28) * t164;
t200 = MDP(29) * t168;
t125 = -t165 * t138 + t139 * t169;
t153 = pkin(4) * t211;
t124 = -t125 + t153;
t199 = t124 * MDP(22);
t131 = t142 * t168 + t164 * t211;
t198 = t131 * MDP(26);
t197 = t142 * MDP(13);
t189 = -qJ(5) * t165 - pkin(3);
t148 = -pkin(4) * t169 + t189;
t196 = t148 * MDP(22);
t195 = MDP(12) + MDP(27);
t191 = qJ(5) * t211;
t190 = t168 * t164 * MDP(24);
t188 = MDP(20) - t214;
t187 = -MDP(17) + t188;
t186 = -t193 + t205;
t184 = -MDP(25) * t164 - MDP(26) * t168;
t183 = MDP(28) * t168 - t164 * MDP(29);
t182 = t200 + t201;
t123 = t191 - t126;
t152 = pkin(9) * t212;
t137 = t152 + (-pkin(3) - t215) * t162;
t179 = MDP(13) + t184;
t178 = MDP(19) + t183;
t177 = -qJ(5) * t143 + t137;
t116 = pkin(5) * t143 + pkin(11) * t211 + t124;
t117 = t218 * t142 + t177;
t112 = t116 * t168 - t117 * t164;
t113 = t116 * t164 + t117 * t168;
t176 = t132 * MDP(25) + t112 * MDP(28) - t113 * MDP(29) + t198;
t175 = (-MDP(28) * t218 + MDP(25)) * t168 + (MDP(29) * t218 - MDP(26)) * t164;
t174 = -pkin(4) * MDP(19) + MDP(14) + t175;
t167 = sin(qJ(2));
t163 = cos(pkin(6));
t161 = sin(pkin(6));
t159 = t169 ^ 2;
t158 = t168 ^ 2;
t157 = t165 ^ 2;
t156 = t164 ^ 2;
t155 = t160 ^ 2;
t151 = t217 * t169;
t150 = t217 * t165;
t146 = -t218 * t169 + t189;
t145 = t162 * t216 + t192;
t144 = t162 * t215 - t152;
t140 = -t160 * t161 * t171 + t162 * t163;
t130 = t163 * t212 + (t166 * t209 + t167 * t170) * t161;
t129 = -t163 * t211 + (t166 * t167 - t170 * t209) * t161;
t128 = t146 * t168 + t150 * t164;
t127 = -t146 * t164 + t150 * t168;
t121 = t130 * t169 + t140 * t165;
t120 = t130 * t165 - t140 * t169;
t119 = pkin(4) * t142 + t177;
t118 = -pkin(5) * t142 - t123;
t115 = t120 * t164 + t129 * t168;
t114 = t120 * t168 - t129 * t164;
t1 = [MDP(1) + (t120 ^ 2 + t121 ^ 2 + t129 ^ 2) * MDP(22); (-t129 * t162 - t140 * t211) * MDP(10) + (-t130 * t162 + t140 * t212) * MDP(11) + (t120 * t143 - t121 * t142) * MDP(19) + (t119 * t129 + t120 * t124 - t121 * t123) * MDP(22) + (t114 * t143 - t121 * t131) * MDP(28) + (-t115 * t143 + t121 * t132) * MDP(29) + (MDP(3) * t171 - MDP(4) * t167) * t161 + t194 * (t120 * t211 + t129 * t142) + t193 * (t121 * t211 + t129 * t143); MDP(2) + (t119 ^ 2 + t123 ^ 2 + t124 ^ 2) * MDP(22) + t155 * t166 ^ 2 * MDP(5) + (0.2e1 * MDP(7) * t212 + t210) * t162 + t195 * t143 ^ 2 + (0.2e1 * t131 * MDP(24) + MDP(25) * t225 + t203) * t132 + (0.2e1 * MDP(6) * t166 + t206) * t155 * t170 + (t112 * t143 - t118 * t131) * t220 + (t123 * t142 + t124 * t143) * t223 + (-t113 * t143 + t118 * t132) * t219 + 0.2e1 * (-t145 * t162 - t155 * t216) * MDP(11) + (-t125 * t211 + t137 * t142) * t224 + (-t119 * t142 - t124 * t211) * t222 + 0.2e1 * (t126 * t211 + t137 * t143) * MDP(18) + (-t119 * t143 + t123 * t211) * t221 + 0.2e1 * (t144 * t162 + t155 * t215) * MDP(10) + (-t197 + t198) * t225 + 0.2e1 * (-t143 * MDP(14) + t142 * MDP(15) + MDP(8) * t162) * t211; -t130 * MDP(11) + (t114 * t165 + t121 * t207) * MDP(28) + (-t115 * t165 - t121 * t208) * MDP(29) + (t193 * t165 - t194 * t169 - MDP(10) + t196) * t129 + (MDP(19) + t213) * (t120 * t165 + t121 * t169); t210 + t144 * MDP(10) - t145 * MDP(11) + (t127 * t143 - t131 * t151) * MDP(28) + (-t128 * t143 + t132 * t151) * MDP(29) + (MDP(7) * t166 + MDP(8) * t170) * t160 + (-t142 * MDP(17) - t143 * MDP(18)) * pkin(3) + (-t142 * MDP(20) - t143 * MDP(21) + MDP(22) * t119) * t148 + (-MDP(14) * t211 - t197 + t137 * MDP(18) + t124 * MDP(19) - t119 * MDP(21) + t195 * t143 + (t143 * MDP(19) + t194 * t211 + t199) * pkin(10) + t176) * t165 + (-MDP(15) * t211 - t137 * MDP(17) - t123 * MDP(19) + t119 * MDP(20) - t164 * t203 + (-t131 * t164 - t132 * t168) * MDP(24) + t183 * t118 + t179 * t143 + (-t142 * MDP(19) - t123 * MDP(22) + t193 * t211) * pkin(10)) * t169; pkin(3) * t169 * t224 + MDP(9) + (t169 * t222 + t196) * t148 + (MDP(23) * t156 + 0.2e1 * t190 + t204) * t159 + (t195 + t204) * t157 + 0.2e1 * (-pkin(3) * MDP(18) - t148 * MDP(21) + t179 * t169) * t165 + (t127 * t165 + t151 * t207) * t220 + (-t128 * t165 - t151 * t208) * t219 + (t157 + t159) * pkin(10) * t223; t187 * t120 + (t182 + t186) * t121; -t160 * t206 + t125 * MDP(17) - t126 * MDP(18) + (-t125 + 0.2e1 * t153) * MDP(20) + (-0.2e1 * t191 + t126) * MDP(21) + (-pkin(4) * t124 - qJ(5) * t123) * MDP(22) + t132 * t202 + (t131 * t168 - t132 * t164) * MDP(24) + (-qJ(5) * t131 + t118 * t164) * MDP(28) + (qJ(5) * t132 + t118 * t168) * MDP(29) + (-qJ(5) * MDP(19) - MDP(15)) * t142 + t174 * t143; t182 * t151 + t174 * t165 + (MDP(15) - t164 * t202 + (t156 - t158) * MDP(24) + t178 * qJ(5)) * t169 + (t187 * t165 + t186 * t169) * pkin(10); -0.2e1 * t190 + t158 * MDP(23) + MDP(16) + (-2 * MDP(20) + t214) * pkin(4) + (t221 + 0.2e1 * t200 + 0.2e1 * t201 + t205) * qJ(5); t120 * MDP(22); -MDP(20) * t211 + t143 * t178 + t199; (t178 + t213) * t165; t188; MDP(22); MDP(28) * t114 - MDP(29) * t115; t143 * MDP(27) + t176; MDP(27) * t165 + MDP(28) * t127 - MDP(29) * t128 + t169 * t184; t175; t183; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

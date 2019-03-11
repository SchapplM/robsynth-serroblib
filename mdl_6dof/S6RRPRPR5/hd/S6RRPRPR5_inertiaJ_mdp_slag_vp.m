% Calculate joint inertia matrix for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR5_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:35:15
% EndTime: 2019-03-09 10:35:20
% DurationCPUTime: 1.39s
% Computational Cost: add. (2271->247), mult. (5397->367), div. (0->0), fcn. (6033->12), ass. (0->121)
t265 = MDP(23) * qJ(5);
t195 = sin(pkin(12));
t198 = cos(pkin(12));
t243 = t195 ^ 2 + t198 ^ 2;
t264 = t243 * t265;
t201 = sin(qJ(6));
t204 = cos(qJ(6));
t177 = t195 * t201 - t198 * t204;
t178 = t195 * t204 + t198 * t201;
t253 = pkin(10) + qJ(5);
t181 = t253 * t195;
t182 = t253 * t198;
t214 = t178 * MDP(26) - t177 * MDP(27) + (-t181 * t204 - t182 * t201) * MDP(29) - (-t181 * t201 + t182 * t204) * MDP(30);
t219 = MDP(20) * t195 + MDP(21) * t198;
t263 = -qJ(5) * t219 - MDP(16) + t214;
t262 = 0.2e1 * MDP(20);
t261 = 0.2e1 * MDP(21);
t260 = 2 * MDP(22);
t259 = -2 * MDP(25);
t258 = 0.2e1 * MDP(29);
t257 = 0.2e1 * MDP(30);
t206 = cos(qJ(2));
t256 = pkin(1) * t206;
t202 = sin(qJ(4));
t255 = pkin(10) * t202;
t254 = pkin(8) + qJ(3);
t252 = pkin(4) * MDP(23);
t200 = cos(pkin(6));
t184 = t200 * t256;
t197 = sin(pkin(6));
t203 = sin(qJ(2));
t246 = t197 * t203;
t162 = t200 * pkin(2) - t246 * t254 + t184;
t224 = t200 * t203 * pkin(1);
t245 = t197 * t206;
t166 = t245 * t254 + t224;
t196 = sin(pkin(11));
t199 = cos(pkin(11));
t148 = t162 * t196 + t166 * t199;
t145 = pkin(9) * t200 + t148;
t167 = t196 * t246 - t199 * t245;
t168 = (t196 * t206 + t199 * t203) * t197;
t179 = (-pkin(2) * t206 - pkin(1)) * t197;
t151 = t167 * pkin(3) - t168 * pkin(9) + t179;
t205 = cos(qJ(4));
t137 = -t145 * t202 + t151 * t205;
t133 = -pkin(4) * t167 - t137;
t251 = t133 * t195;
t250 = t133 * t198;
t187 = pkin(2) * t196 + pkin(9);
t249 = t187 * t195;
t248 = t187 * t205;
t191 = t197 ^ 2;
t247 = t191 * t203;
t244 = t200 * MDP(8);
t138 = t145 * t205 + t151 * t202;
t132 = qJ(5) * t167 + t138;
t147 = t162 * t199 - t166 * t196;
t144 = -pkin(3) * t200 - t147;
t158 = t168 * t202 - t200 * t205;
t159 = t168 * t205 + t200 * t202;
t136 = pkin(4) * t158 - qJ(5) * t159 + t144;
t129 = t132 * t198 + t136 * t195;
t188 = -pkin(2) * t199 - pkin(3);
t176 = -t205 * pkin(4) - qJ(5) * t202 + t188;
t156 = t176 * t195 + t198 * t248;
t242 = MDP(18) * t188;
t241 = MDP(20) * t198;
t240 = MDP(21) * t195;
t239 = MDP(23) * t133;
t238 = MDP(23) * t187;
t170 = t177 * t202;
t237 = MDP(24) * t170;
t236 = MDP(24) * t178;
t235 = MDP(26) * t170;
t169 = t178 * t202;
t234 = MDP(27) * t169;
t233 = MDP(29) * t177;
t149 = t159 * t195 - t167 * t198;
t150 = t159 * t198 + t167 * t195;
t139 = t149 * t204 + t150 * t201;
t232 = t139 * MDP(27);
t140 = -t149 * t201 + t150 * t204;
t231 = t140 * MDP(26);
t230 = t158 * MDP(28);
t229 = t159 * MDP(13);
t228 = t159 * MDP(14);
t227 = t170 * MDP(30);
t226 = t187 * MDP(19);
t225 = t188 * MDP(19);
t128 = -t132 * t195 + t136 * t198;
t222 = t243 * MDP(22);
t221 = -t128 * t195 + t129 * t198;
t172 = t198 * t176;
t155 = -t195 * t248 + t172;
t220 = -t155 * t195 + t156 * t198;
t218 = -t169 * MDP(29) + t227;
t217 = (MDP(6) * t203 + MDP(7) * t206) * t197;
t216 = t240 - t241 - t252;
t153 = -t198 * t255 + t172 + (-pkin(5) - t249) * t205;
t154 = -t195 * t255 + t156;
t141 = t153 * t204 - t154 * t201;
t142 = t153 * t201 + t154 * t204;
t215 = MDP(29) * t141 - MDP(30) * t142 - t235;
t212 = MDP(30) * t178 + t216 + t233;
t211 = MDP(20) * t149 + MDP(21) * t150 + MDP(29) * t139 + MDP(30) * t140 + t239;
t126 = pkin(5) * t158 - pkin(10) * t150 + t128;
t127 = -pkin(10) * t149 + t129;
t124 = t126 * t204 - t127 * t201;
t125 = t126 * t201 + t127 * t204;
t210 = t124 * MDP(29) - t125 * MDP(30) + t230 + t231 - t232;
t209 = (-t149 * t198 + t150 * t195) * MDP(22) + t221 * MDP(23);
t194 = t205 ^ 2;
t193 = t202 ^ 2;
t189 = -pkin(5) * t198 - pkin(4);
t175 = pkin(8) * t245 + t224;
t174 = -pkin(8) * t246 + t184;
t173 = (pkin(5) * t195 + t187) * t202;
t152 = t169 * t158;
t130 = pkin(5) * t149 + t133;
t1 = [(t128 ^ 2 + t129 ^ 2 + t133 ^ 2) * MDP(23) + t167 ^ 2 * MDP(17) + (t147 ^ 2 + t148 ^ 2 + t179 ^ 2) * MDP(12) + MDP(1) + (MDP(4) * t203 + 0.2e1 * MDP(5) * t206) * t247 + (0.2e1 * t167 * MDP(15) + t229) * t159 + (MDP(24) * t140 + t139 * t259) * t140 + (0.2e1 * t217 + t244) * t200 + (-0.2e1 * MDP(16) * t167 - 0.2e1 * t228 + t230 + 0.2e1 * t231 - 0.2e1 * t232) * t158 + 0.2e1 * (t174 * t200 + t191 * t256) * MDP(9) + 0.2e1 * (-pkin(1) * t247 - t175 * t200) * MDP(10) + 0.2e1 * (t137 * t167 + t144 * t158) * MDP(18) + 0.2e1 * (-t138 * t167 + t144 * t159) * MDP(19) + 0.2e1 * (-t147 * t168 - t148 * t167) * MDP(11) + (t128 * t158 + t133 * t149) * t262 + (t124 * t158 + t130 * t139) * t258 + (-t129 * t158 + t133 * t150) * t261 + (-t125 * t158 + t130 * t140) * t257 + (-t128 * t150 - t129 * t149) * t260; t244 + t174 * MDP(9) - t175 * MDP(10) + t159 * t225 + (-t156 * t149 - t155 * t150) * MDP(22) + (t128 * t155 + t129 * t156) * MDP(23) + (t139 * t170 - t140 * t169) * MDP(25) - t152 * MDP(27) + (t130 * t169 + t139 * t173) * MDP(29) + (-t130 * t170 + t140 * t173) * MDP(30) - t140 * t237 + t217 + (MDP(20) * t155 - MDP(21) * t156 + t215 + t242) * t158 + (t228 - t144 * MDP(18) - t128 * MDP(20) + t129 * MDP(21) + (MDP(16) - t226) * t167 - t210) * t205 + ((-t167 * t196 - t168 * t199) * MDP(11) + (t147 * t199 + t148 * t196) * MDP(12)) * pkin(2) + (-t158 * MDP(14) + t144 * MDP(19) + (t149 * t187 + t251) * MDP(20) + (t150 * t187 + t250) * MDP(21) + (-t128 * t198 - t129 * t195) * MDP(22) + t133 * t238 + t229 + (-t187 * MDP(18) + MDP(15)) * t167) * t202; MDP(8) + t193 * MDP(13) + 0.2e1 * t202 * t225 + (t187 ^ 2 * t193 + t155 ^ 2 + t156 ^ 2) * MDP(23) + t194 * MDP(28) + (t196 ^ 2 + t199 ^ 2) * MDP(12) * pkin(2) ^ 2 - (t169 * t259 - t237) * t170 + 0.2e1 * (MDP(14) * t202 + t234 + t235 - t242) * t205 + (-t155 * t205 + t193 * t249) * t262 + (t187 * t193 * t198 + t156 * t205) * t261 + (-t141 * t205 + t169 * t173) * t258 + (t142 * t205 - t170 * t173) * t257 + (-t155 * t198 - t156 * t195) * t202 * t260; t158 * t227 + t179 * MDP(12) - t152 * MDP(29) + (MDP(18) * t167 - t211) * t205 + (-t167 * MDP(19) - t158 * t219 + t209) * t202; (t220 - t248) * t202 * MDP(23); MDP(12) + (t193 * t243 + t194) * MDP(23); t159 * MDP(15) + t167 * MDP(17) + t137 * MDP(18) - t138 * MDP(19) + (-pkin(4) * t149 - t250) * MDP(20) + (-pkin(4) * t150 + t251) * MDP(21) + t221 * MDP(22) - pkin(4) * t239 + t140 * t236 + (-t139 * t178 - t140 * t177) * MDP(25) + (t130 * t177 + t139 * t189) * MDP(29) + (t130 * t178 + t140 * t189) * MDP(30) + t209 * qJ(5) + t263 * t158; -t170 * t236 + (-t169 * t178 + t170 * t177) * MDP(25) + (t169 * t189 + t173 * t177) * MDP(29) + (-t170 * t189 + t173 * t178) * MDP(30) + (-t226 - t263) * t205 + (MDP(15) - t219 * pkin(4) + (-MDP(18) + t216) * t187) * t202 + (MDP(22) + t265) * t220; (-MDP(19) + t222 + t264) * t202 + (MDP(18) - t212) * t205; 0.2e1 * t189 * t233 + MDP(17) + (-0.2e1 * t240 + 0.2e1 * t241 + t252) * pkin(4) + (t177 * t259 + t189 * t257 + t236) * t178 + (0.2e1 * t222 + t264) * qJ(5); t211; (t219 + t238) * t202 - t218; -t205 * MDP(23); t212; MDP(23); t210; -MDP(28) * t205 + t215 - t234; t218; t214; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

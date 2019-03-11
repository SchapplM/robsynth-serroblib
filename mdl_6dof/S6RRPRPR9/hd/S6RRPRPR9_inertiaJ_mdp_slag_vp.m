% Calculate joint inertia matrix for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR9_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:01:13
% EndTime: 2019-03-09 11:01:17
% DurationCPUTime: 1.53s
% Computational Cost: add. (2656->274), mult. (5905->409), div. (0->0), fcn. (6770->12), ass. (0->123)
t203 = sin(pkin(11));
t262 = pkin(9) + qJ(3);
t186 = t262 * t203;
t206 = cos(pkin(11));
t188 = t262 * t206;
t209 = sin(qJ(4));
t266 = cos(qJ(4));
t167 = t266 * t186 + t188 * t209;
t169 = -t209 * t186 + t188 * t266;
t184 = t203 * t266 + t209 * t206;
t274 = t184 * MDP(17) - t167 * MDP(20) - t169 * MDP(21);
t273 = 2 * MDP(13);
t272 = -2 * MDP(16);
t271 = 0.2e1 * MDP(20);
t270 = 0.2e1 * MDP(21);
t269 = 2 * MDP(24);
t268 = -2 * MDP(27);
t267 = 2 * MDP(32);
t265 = cos(qJ(6));
t210 = sin(qJ(2));
t264 = pkin(1) * t210;
t211 = cos(qJ(2));
t263 = pkin(1) * t211;
t261 = pkin(10) + qJ(5);
t260 = pkin(2) * MDP(14);
t259 = pkin(4) * MDP(25);
t202 = sin(pkin(12));
t258 = t184 * t202;
t205 = cos(pkin(12));
t257 = t184 * t205;
t204 = sin(pkin(6));
t256 = t204 * t210;
t255 = t204 * t211;
t207 = cos(pkin(6));
t254 = t207 * MDP(8);
t233 = pkin(8) * t255;
t173 = t233 + (qJ(3) + t264) * t207;
t174 = (-pkin(2) * t211 - qJ(3) * t210 - pkin(1)) * t204;
t154 = -t173 * t203 + t206 * t174;
t177 = t203 * t207 + t206 * t256;
t146 = -pkin(3) * t255 - pkin(9) * t177 + t154;
t155 = t206 * t173 + t203 * t174;
t189 = t203 * t256;
t229 = t206 * t207 - t189;
t150 = pkin(9) * t229 + t155;
t138 = t209 * t146 + t150 * t266;
t135 = -qJ(5) * t255 + t138;
t157 = t177 * t209 - t229 * t266;
t158 = t177 * t266 + t209 * t229;
t191 = pkin(8) * t256;
t196 = -pkin(3) * t206 - pkin(2);
t163 = t189 * pkin(3) + t191 + (t196 - t263) * t207;
t141 = t157 * pkin(4) - t158 * qJ(5) + t163;
t130 = t205 * t135 + t202 * t141;
t182 = t203 * t209 - t206 * t266;
t164 = pkin(4) * t182 - qJ(5) * t184 + t196;
t148 = t202 * t164 + t205 * t169;
t253 = t202 ^ 2 + t205 ^ 2;
t251 = MDP(15) * t184;
t250 = MDP(19) * t211;
t249 = MDP(22) * t205;
t248 = MDP(23) * t202;
t137 = t146 * t266 - t209 * t150;
t136 = pkin(4) * t255 - t137;
t247 = MDP(25) * t136;
t208 = sin(qJ(6));
t218 = -t208 * t202 + t205 * t265;
t160 = t218 * t184;
t246 = MDP(26) * t160;
t183 = t202 * t265 + t208 * t205;
t245 = MDP(26) * t183;
t152 = t158 * t202 + t205 * t255;
t153 = t158 * t205 - t202 * t255;
t143 = -t208 * t152 + t153 * t265;
t244 = MDP(28) * t143;
t243 = MDP(28) * t160;
t142 = t152 * t265 + t153 * t208;
t242 = MDP(29) * t142;
t159 = t183 * t184;
t241 = MDP(29) * t159;
t240 = MDP(30) * t157;
t239 = MDP(30) * t182;
t238 = MDP(31) * t218;
t237 = t158 * MDP(17);
t235 = t203 * MDP(12);
t234 = t206 * MDP(11);
t232 = qJ(3) * t255;
t231 = MDP(6) * t256;
t230 = MDP(18) * t255;
t129 = -t135 * t202 + t205 * t141;
t147 = t205 * t164 - t169 * t202;
t228 = t253 * MDP(25);
t227 = t129 * t205 + t130 * t202;
t226 = -t129 * t202 + t130 * t205;
t225 = t147 * t205 + t148 * t202;
t224 = -t147 * t202 + t148 * t205;
t223 = -t154 * t203 + t155 * t206;
t222 = MDP(22) * t202 + MDP(23) * t205;
t144 = pkin(5) * t182 - pkin(10) * t257 + t147;
t145 = -pkin(10) * t258 + t148;
t132 = t144 * t265 - t208 * t145;
t133 = t208 * t144 + t145 * t265;
t221 = MDP(31) * t132 - MDP(32) * t133;
t220 = MDP(31) * t159 + MDP(32) * t160;
t219 = -MDP(32) * t183 + t238;
t217 = t219 - t248 + t249;
t185 = t261 * t202;
t187 = t261 * t205;
t216 = t183 * MDP(28) + t218 * MDP(29) + (-t185 * t265 - t208 * t187) * MDP(31) - (-t208 * t185 + t187 * t265) * MDP(32);
t215 = MDP(20) + t217;
t214 = -qJ(5) * t222 - MDP(18) + t216;
t199 = t204 ^ 2;
t195 = -pkin(5) * t205 - pkin(4);
t179 = t207 * t264 + t233;
t178 = t207 * t263 - t191;
t176 = t191 + (-pkin(2) - t263) * t207;
t156 = pkin(5) * t258 + t167;
t131 = t152 * pkin(5) + t136;
t128 = -pkin(10) * t152 + t130;
t127 = pkin(5) * t157 - pkin(10) * t153 + t129;
t126 = t208 * t127 + t128 * t265;
t125 = t127 * t265 - t208 * t128;
t1 = [t199 * t210 ^ 2 * MDP(4) + (t129 ^ 2 + t130 ^ 2 + t136 ^ 2) * MDP(25) + t158 ^ 2 * MDP(15) + (t154 ^ 2 + t155 ^ 2 + t176 ^ 2) * MDP(14) + MDP(1) + (0.2e1 * t231 + t254) * t207 + (MDP(26) * t143 + t142 * t268) * t143 + (0.2e1 * (MDP(7) * t207 - t237) * t204 + (0.2e1 * MDP(5) * t210 + t250) * t199) * t211 + (t158 * t272 + 0.2e1 * t230 + t240 - 0.2e1 * t242 + 0.2e1 * t244) * t157 + (-t129 * t153 - t130 * t152) * t269 + 0.2e1 * (t125 * t157 + t131 * t142) * MDP(31) + 0.2e1 * (t129 * t157 + t136 * t152) * MDP(22) + (-t126 * t157 + t131 * t143) * t267 + 0.2e1 * (-t130 * t157 + t136 * t153) * MDP(23) + (-t154 * t177 + t155 * t229) * t273 + 0.2e1 * (-t179 * t207 - t199 * t264) * MDP(10) + (-t137 * t255 + t157 * t163) * t271 + (t138 * t255 + t158 * t163) * t270 + 0.2e1 * (t155 * t255 + t176 * t177) * MDP(12) + 0.2e1 * (t178 * t207 + t199 * t263) * MDP(9) + 0.2e1 * (-t154 * t255 - t176 * t229) * MDP(11); (MDP(7) - t274) * t255 + (-pkin(2) * t176 + qJ(3) * t223) * MDP(14) + (-t147 * t153 - t148 * t152 - t184 * t227) * MDP(24) + (t157 * t196 + t163 * t182) * MDP(20) + (-t142 * t160 - t143 * t159) * MDP(27) + (t129 * t147 + t130 * t148 + t136 * t167) * MDP(25) + t178 * MDP(9) - t179 * MDP(10) + (t125 * t182 + t131 * t159 + t132 * t157 + t142 * t156) * MDP(31) + (-t142 * t182 - t157 * t159) * MDP(29) + (t143 * t182 + t157 * t160) * MDP(28) + (-t126 * t182 + t131 * t160 - t133 * t157 + t143 * t156) * MDP(32) + (-t157 * t184 - t158 * t182) * MDP(16) + t182 * t230 + (t158 * t196 + t163 * t184) * MDP(21) + t143 * t246 + t158 * t251 + t254 + ((t203 * t177 + t206 * t229) * qJ(3) + t223) * MDP(13) + (pkin(2) * t229 - t176 * t206 + t203 * t232) * MDP(11) + (-pkin(2) * t177 + t176 * t203 + t206 * t232) * MDP(12) + (-t130 * t182 + t136 * t257 - t148 * t157 + t153 * t167) * MDP(23) + (t129 * t182 + t136 * t258 + t147 * t157 + t152 * t167) * MDP(22) + t157 * t239 + t231; MDP(8) + (t147 ^ 2 + t148 ^ 2 + t167 ^ 2) * MDP(25) + (t196 * t270 + t251) * t184 + (t159 * t268 + t246) * t160 + (0.2e1 * t234 - 0.2e1 * t235 + t260) * pkin(2) + (t184 * t272 + t196 * t271 + t239 - 0.2e1 * t241 + 0.2e1 * t243) * t182 + 0.2e1 * t220 * t156 + 0.2e1 * (MDP(22) * t147 - MDP(23) * t148 + t221) * t182 + 0.2e1 * (-MDP(24) * t225 + t167 * t222) * t184 + (MDP(14) * qJ(3) + t273) * (t203 ^ 2 + t206 ^ 2) * qJ(3); -t229 * MDP(11) + t177 * MDP(12) + t176 * MDP(14) + t158 * MDP(21) + (-t152 * t202 - t153 * t205) * MDP(24) + t227 * MDP(25) + t215 * t157; -t234 + t235 - t260 + t225 * MDP(25) + (-MDP(24) * t253 + MDP(21)) * t184 + t215 * t182; MDP(14) + t228; t237 - t204 * t250 + t137 * MDP(20) - t138 * MDP(21) + (-pkin(4) * t152 - t136 * t205) * MDP(22) + (-pkin(4) * t153 + t136 * t202) * MDP(23) + t226 * MDP(24) - pkin(4) * t247 + t143 * t245 + (-t142 * t183 + t143 * t218) * MDP(27) + (-t131 * t218 + t142 * t195) * MDP(31) + (t131 * t183 + t143 * t195) * MDP(32) + ((-t152 * t205 + t153 * t202) * MDP(24) + t226 * MDP(25)) * qJ(5) + t214 * t157; (-pkin(4) * t258 - t167 * t205) * MDP(22) + (-pkin(4) * t257 + t167 * t202) * MDP(23) + t224 * MDP(24) + (-pkin(4) * t167 + qJ(5) * t224) * MDP(25) + t160 * t245 + (-t159 * t183 + t160 * t218) * MDP(27) + (-t156 * t218 + t159 * t195) * MDP(31) + (t156 * t183 + t160 * t195) * MDP(32) + t214 * t182 + t274; 0; -0.2e1 * t195 * t238 + MDP(19) + (-0.2e1 * t248 + 0.2e1 * t249 + t259) * pkin(4) + (t195 * t267 - t218 * t268 + t245) * t183 + (qJ(5) * t228 + t253 * t269) * qJ(5); MDP(22) * t152 + MDP(23) * t153 + MDP(31) * t142 + MDP(32) * t143 + t247; MDP(25) * t167 + t184 * t222 + t220; 0; -t217 - t259; MDP(25); MDP(31) * t125 - MDP(32) * t126 + t240 - t242 + t244; t221 + t239 - t241 + t243; t219; t216; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

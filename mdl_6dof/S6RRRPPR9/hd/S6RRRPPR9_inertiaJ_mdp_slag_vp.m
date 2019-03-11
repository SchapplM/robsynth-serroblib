% Calculate joint inertia matrix for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR9_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:18:00
% EndTime: 2019-03-09 16:18:05
% DurationCPUTime: 1.74s
% Computational Cost: add. (1819->298), mult. (4045->424), div. (0->0), fcn. (4322->10), ass. (0->120)
t201 = sin(pkin(11));
t203 = cos(pkin(11));
t231 = MDP(19) - MDP(24);
t232 = MDP(18) + MDP(22);
t271 = t232 * t201 + t231 * t203;
t202 = sin(pkin(6));
t273 = 0.2e1 * t202;
t272 = MDP(21) * pkin(9);
t270 = MDP(21) + MDP(25);
t206 = sin(qJ(3));
t269 = 0.2e1 * t206;
t268 = 2 * MDP(16);
t267 = -2 * MDP(27);
t266 = 2 * MDP(32);
t265 = pkin(4) + pkin(5);
t207 = sin(qJ(2));
t264 = pkin(1) * t207;
t210 = cos(qJ(2));
t263 = pkin(1) * t210;
t209 = cos(qJ(3));
t262 = pkin(9) * t209;
t261 = pkin(10) * t206;
t260 = -pkin(10) + qJ(4);
t259 = pkin(2) * MDP(17);
t258 = pkin(3) * MDP(21);
t257 = pkin(9) * MDP(17);
t208 = cos(qJ(6));
t256 = t201 * t208;
t255 = t202 * t207;
t254 = t202 * t210;
t253 = t203 * t206;
t204 = cos(pkin(6));
t252 = t204 * MDP(8);
t251 = t207 * MDP(6);
t190 = pkin(8) * t255;
t165 = t190 + (-pkin(2) - t263) * t204;
t172 = -t204 * t209 + t206 * t255;
t174 = t204 * t206 + t209 * t255;
t145 = pkin(3) * t172 - qJ(4) * t174 + t165;
t229 = pkin(8) * t254;
t166 = t229 + (pkin(9) + t264) * t204;
t168 = (-pkin(2) * t210 - pkin(9) * t207 - pkin(1)) * t202;
t149 = t166 * t209 + t168 * t206;
t146 = -qJ(4) * t254 + t149;
t136 = t201 * t145 + t203 * t146;
t148 = -t206 * t166 + t209 * t168;
t183 = -pkin(3) * t209 - qJ(4) * t206 - pkin(2);
t164 = t201 * t183 + t203 * t262;
t249 = MDP(15) * t210;
t248 = MDP(18) * t203;
t247 = MDP(19) * t201;
t147 = pkin(3) * t254 - t148;
t246 = MDP(21) * t147;
t245 = MDP(22) * t203;
t244 = MDP(24) * t201;
t228 = qJ(5) * t201 + pkin(3);
t181 = -pkin(4) * t203 - t228;
t243 = MDP(25) * t181;
t205 = sin(qJ(6));
t179 = t201 * t205 + t203 * t208;
t170 = t179 * t206;
t242 = MDP(26) * t170;
t180 = -t203 * t205 + t256;
t241 = MDP(26) * t180;
t240 = MDP(28) * t170;
t169 = t205 * t253 - t206 * t256;
t239 = MDP(29) * t169;
t238 = MDP(30) * t172;
t237 = MDP(30) * t209;
t236 = MDP(31) * t179;
t154 = t174 * t201 + t203 * t254;
t155 = t174 * t203 - t201 * t254;
t140 = -t154 * t208 + t155 * t205;
t235 = t140 * MDP(29);
t141 = t154 * t205 + t155 * t208;
t234 = t141 * MDP(28);
t233 = t174 * MDP(13);
t230 = MDP(20) + MDP(23);
t133 = t172 * qJ(5) + t136;
t227 = qJ(5) * t203 - pkin(9);
t135 = t145 * t203 - t201 * t146;
t189 = t201 * t262;
t163 = t183 * t203 - t189;
t159 = -qJ(5) * t209 + t164;
t134 = -pkin(4) * t172 - t135;
t224 = t133 * t203 + t134 * t201;
t223 = -t135 * t201 + t136 * t203;
t196 = t209 * pkin(4);
t160 = -t163 + t196;
t222 = t159 * t203 + t160 * t201;
t221 = -t163 * t201 + t164 * t203;
t220 = t201 * MDP(18) + MDP(19) * t203;
t219 = t201 * MDP(22) - MDP(24) * t203;
t150 = pkin(5) * t209 + t189 + t196 + (-t183 - t261) * t203;
t152 = t201 * t261 + t159;
t138 = t150 * t208 - t152 * t205;
t139 = t150 * t205 + t152 * t208;
t218 = MDP(31) * t138 - MDP(32) * t139;
t217 = MDP(31) * t169 + MDP(32) * t170;
t216 = MDP(31) * t208 - t205 * MDP(32);
t215 = qJ(5) * t155 - t147;
t214 = MDP(22) + t216;
t184 = t260 * t201;
t185 = t260 * t203;
t213 = MDP(28) * t180 - MDP(29) * t179 + (t184 * t208 - t185 * t205) * MDP(31) - (t184 * t205 + t185 * t208) * MDP(32);
t212 = -MDP(14) - t213;
t198 = t202 ^ 2;
t177 = t204 * t264 + t229;
t176 = t204 * t263 - t190;
t175 = t203 * t265 + t228;
t167 = (pkin(4) * t201 - t227) * t206;
t158 = (-t201 * t265 + t227) * t206;
t151 = t203 * qJ(4) * t154;
t137 = pkin(4) * t154 - t215;
t132 = -t154 * t265 + t215;
t131 = pkin(10) * t154 + t133;
t130 = -pkin(10) * t155 - t172 * t265 - t135;
t129 = t130 * t205 + t131 * t208;
t128 = t130 * t208 - t131 * t205;
t1 = [(t133 ^ 2 + t134 ^ 2 + t137 ^ 2) * MDP(25) + t198 * t207 ^ 2 * MDP(4) + (t135 ^ 2 + t136 ^ 2 + t147 ^ 2) * MDP(21) + t174 ^ 2 * MDP(11) + MDP(1) + (t251 * t273 + t252) * t204 + (MDP(26) * t141 + t140 * t267) * t141 + ((MDP(7) * t204 - t233) * t273 + (0.2e1 * MDP(5) * t207 + t249) * t198) * t210 + (-0.2e1 * t174 * MDP(12) + 0.2e1 * MDP(14) * t254 - 0.2e1 * t234 + 0.2e1 * t235 + t238) * t172 + (-t148 * t254 + t165 * t172) * t268 + 0.2e1 * (-t177 * t204 - t198 * t264) * MDP(10) + 0.2e1 * (t149 * t254 + t165 * t174) * MDP(17) + 0.2e1 * (t176 * t204 + t198 * t263) * MDP(9) + 0.2e1 * (t133 * t172 - t137 * t155) * MDP(24) + 0.2e1 * (t135 * t172 + t147 * t154) * MDP(18) + (t129 * t172 + t132 * t141) * t266 + 0.2e1 * (-t136 * t172 + t147 * t155) * MDP(19) + 0.2e1 * (-t134 * t172 + t137 * t154) * MDP(22) + 0.2e1 * (-t128 * t172 + t132 * t140) * MDP(31) + 0.2e1 * (-t133 * t154 + t134 * t155) * MDP(23) + 0.2e1 * (-t135 * t155 - t136 * t154) * MDP(20); -t172 * t237 + t141 * t242 + (t133 * t159 + t134 * t160 + t137 * t167) * MDP(25) + (-t140 * t170 - t141 * t169) * MDP(27) + t176 * MDP(9) - t177 * MDP(10) + t252 + (t135 * t163 + t136 * t164) * MDP(21) + (-t154 * t159 + t155 * t160) * MDP(23) + (-t154 * t164 - t155 * t163) * MDP(20) + (-t129 * t209 + t132 * t170 + t139 * t172 + t141 * t158) * MDP(32) + (t128 * t209 + t132 * t169 - t138 * t172 + t140 * t158) * MDP(31) + (-t140 * t209 + t169 * t172) * MDP(29) + (t141 * t209 - t170 * t172) * MDP(28) + (-t135 * t209 + t163 * t172) * MDP(18) + (t134 * t209 + t154 * t167 - t160 * t172) * MDP(22) + (-t133 * t209 - t155 * t167 + t159 * t172) * MDP(24) + (t136 * t209 - t164 * t172) * MDP(19) + (-pkin(2) * t172 - t165 * t209) * MDP(16) + (t209 * MDP(12) - t259) * t174 + (t251 + (MDP(7) + (-MDP(14) + t257) * t209) * t210) * t202 + (t174 * MDP(11) - MDP(13) * t254 + (-t133 * t201 + t134 * t203) * MDP(23) + (-t135 * t203 - t136 * t201) * MDP(20) - t172 * MDP(12) + t165 * MDP(17) + t220 * t147 + t219 * t137 + (MDP(16) * t254 + MDP(18) * t154 + MDP(19) * t155 + t246) * pkin(9)) * t206; MDP(8) + (t163 ^ 2 + t164 ^ 2) * MDP(21) + (t159 ^ 2 + t160 ^ 2 + t167 ^ 2) * MDP(25) + (t169 * t267 + t242) * t170 + (MDP(12) * t269 + pkin(2) * t268 + t237 - 0.2e1 * t239 + 0.2e1 * t240) * t209 + 0.2e1 * t217 * t158 + 0.2e1 * (-MDP(18) * t163 + MDP(19) * t164 + MDP(22) * t160 - MDP(24) * t159 + t218) * t209 + ((-t163 * t203 - t164 * t201) * MDP(20) + (-t159 * t201 + t160 * t203) * MDP(23) + t219 * t167) * t269 + (-0.2e1 * t259 + (MDP(11) + (0.2e1 * t220 + t272) * pkin(9)) * t206) * t206; t233 - t202 * t249 + t148 * MDP(16) - t149 * MDP(17) + (-pkin(3) * t154 - t147 * t203) * MDP(18) + (-pkin(3) * t155 + t147 * t201) * MDP(19) + (-t151 + t223) * MDP(20) - pkin(3) * t246 + (-t137 * t203 + t154 * t181) * MDP(22) + (-t151 + t224) * MDP(23) + (-t137 * t201 - t155 * t181) * MDP(24) + t137 * t243 + t141 * t241 + (-t140 * t180 - t141 * t179) * MDP(27) + (t132 * t179 + t140 * t175) * MDP(31) + (t132 * t180 + t141 * t175) * MDP(32) + (t230 * t201 * t155 + t223 * MDP(21) + t224 * MDP(25)) * qJ(4) + (-t271 * qJ(4) + t212) * t172; t221 * MDP(20) + t222 * MDP(23) + t170 * t241 + (-t169 * t180 - t170 * t179) * MDP(27) + (t158 * t179 + t169 * t175) * MDP(31) + (t158 * t180 + t170 * t175) * MDP(32) + (t243 - t244 - t245) * t167 + (-t212 - t257) * t209 + (MDP(13) + t219 * t181 - t220 * pkin(3) + (-MDP(16) + t247 - t248 - t258) * pkin(9)) * t206 + (t221 * MDP(21) + t222 * MDP(25) + t271 * t209) * qJ(4); MDP(15) + 0.2e1 * t175 * t236 + (t243 - 0.2e1 * t244 - 0.2e1 * t245) * t181 + (-0.2e1 * t247 + 0.2e1 * t248 + t258) * pkin(3) + (t175 * t266 + t179 * t267 + t241) * t180 + (qJ(4) * t270 + 0.2e1 * t230) * (t201 ^ 2 + t203 ^ 2) * qJ(4); MDP(25) * t137 - t140 * MDP(31) - t141 * MDP(32) + t154 * t232 + t155 * t231 + t246; MDP(25) * t167 + (t271 + t272) * t206 - t217; -MDP(32) * t180 + t201 * t231 - t203 * t232 - t236 + t243 - t258; t270; t155 * MDP(23) + MDP(25) * t134 - t172 * t214; MDP(23) * t253 + MDP(25) * t160 + t209 * t214; (MDP(25) * qJ(4) + MDP(23)) * t201; 0; MDP(25); t128 * MDP(31) - t129 * MDP(32) + t234 - t235 - t238; t218 + t237 - t239 + t240; t213; 0; t216; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

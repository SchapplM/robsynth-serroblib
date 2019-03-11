% Calculate joint inertia matrix for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR8_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:52:26
% EndTime: 2019-03-09 18:52:30
% DurationCPUTime: 1.42s
% Computational Cost: add. (2462->263), mult. (5484->381), div. (0->0), fcn. (6352->12), ass. (0->131)
t224 = sin(pkin(6));
t298 = 0.2e1 * t224;
t226 = cos(pkin(6));
t229 = sin(qJ(3));
t233 = cos(qJ(3));
t230 = sin(qJ(2));
t272 = t224 * t230;
t193 = t229 * t226 + t233 * t272;
t223 = sin(pkin(12));
t225 = cos(pkin(12));
t212 = t229 * t272;
t252 = t226 * t233 - t212;
t174 = t193 * t223 - t225 * t252;
t258 = MDP(24) + MDP(31);
t297 = t174 * t258;
t175 = t225 * t193 + t223 * t252;
t228 = sin(qJ(5));
t232 = cos(qJ(5));
t234 = cos(qJ(2));
t271 = t224 * t234;
t165 = t175 * t228 + t232 * t271;
t166 = t175 * t232 - t228 * t271;
t227 = sin(qJ(6));
t231 = cos(qJ(6));
t152 = t231 * t165 + t166 * t227;
t153 = -t165 * t227 + t166 * t231;
t296 = t153 * MDP(29) - t152 * MDP(30);
t206 = t223 * t233 + t225 * t229;
t209 = t227 * t232 + t228 * t231;
t176 = t209 * t206;
t171 = t176 * MDP(30);
t208 = t227 * t228 - t231 * t232;
t177 = t208 * t206;
t172 = t177 * MDP(29);
t295 = -t172 - t171;
t216 = pkin(3) * t223 + pkin(10);
t245 = t228 * MDP(25) + t232 * MDP(26);
t279 = pkin(11) + t216;
t201 = t279 * t228;
t202 = t279 * t232;
t251 = t209 * MDP(29) - t208 * MDP(30) + (-t201 * t231 - t202 * t227) * MDP(32) - (-t201 * t227 + t202 * t231) * MDP(33);
t236 = t228 * MDP(22) + t232 * MDP(23) - t216 * t245 + t251;
t246 = MDP(25) * t232 - MDP(26) * t228;
t197 = t208 * MDP(32);
t267 = -t209 * MDP(33) - t197;
t238 = t246 + t267;
t257 = pkin(8) * t271;
t284 = pkin(1) * t230;
t190 = t257 + (pkin(9) + t284) * t226;
t191 = (-pkin(2) * t234 - pkin(9) * t230 - pkin(1)) * t224;
t167 = -t190 * t229 + t233 * t191;
t158 = -pkin(3) * t271 - qJ(4) * t193 + t167;
t168 = t233 * t190 + t229 * t191;
t163 = qJ(4) * t252 + t168;
t147 = t223 * t158 + t225 * t163;
t145 = -pkin(10) * t271 + t147;
t213 = pkin(8) * t272;
t219 = -pkin(3) * t233 - pkin(2);
t283 = pkin(1) * t234;
t181 = t212 * pkin(3) + t213 + (t219 - t283) * t226;
t149 = t174 * pkin(4) - t175 * pkin(10) + t181;
t139 = -t145 * t228 + t232 * t149;
t140 = t145 * t232 + t149 * t228;
t294 = t139 * MDP(25) - t140 * MDP(26);
t293 = -(MDP(16) * t229 + MDP(17) * t233) * pkin(9) + t229 * MDP(13) + t233 * MDP(14);
t292 = 2 * MDP(12);
t291 = 0.2e1 * MDP(16);
t290 = 2 * MDP(18);
t289 = 0.2e1 * MDP(25);
t288 = 0.2e1 * MDP(26);
t287 = -2 * MDP(28);
t286 = 0.2e1 * MDP(32);
t285 = 0.2e1 * MDP(33);
t282 = pkin(5) * t174;
t205 = t223 * t229 - t225 * t233;
t281 = pkin(5) * t205;
t280 = pkin(5) * t231;
t278 = -qJ(4) - pkin(9);
t138 = -pkin(11) * t165 + t140;
t277 = t138 * t231;
t184 = pkin(4) * t205 - pkin(10) * t206 + t219;
t211 = t278 * t233;
t253 = t278 * t229;
t187 = -t225 * t211 + t223 * t253;
t275 = t187 * t232;
t156 = t275 + (-pkin(11) * t206 + t184) * t228;
t276 = t156 * t231;
t274 = t206 * t228;
t273 = t206 * t232;
t270 = t226 * MDP(8);
t269 = t228 * t232;
t268 = t230 * MDP(6);
t266 = MDP(15) * t234;
t263 = t153 * MDP(27);
t262 = t165 * MDP(21);
t261 = t166 * MDP(20);
t260 = t177 * MDP(27);
t259 = t229 * MDP(11);
t256 = t174 * MDP(31) + t296;
t255 = t205 * MDP(31) + t295;
t217 = -pkin(3) * t225 - pkin(4);
t254 = MDP(21) * t269;
t137 = -pkin(11) * t166 + t139 + t282;
t134 = t231 * t137 - t138 * t227;
t159 = t232 * t184 - t187 * t228;
t155 = -pkin(11) * t273 + t159 + t281;
t142 = t231 * t155 - t156 * t227;
t146 = t158 * t225 - t223 * t163;
t185 = -t211 * t223 - t225 * t253;
t144 = pkin(4) * t271 - t146;
t248 = t166 * MDP(22) - t165 * MDP(23);
t160 = t184 * t228 + t275;
t247 = t159 * MDP(25) - t160 * MDP(26);
t135 = t137 * t227 + t277;
t244 = t134 * MDP(32) - t135 * MDP(33);
t143 = t155 * t227 + t276;
t243 = t142 * MDP(32) - t143 * MDP(33);
t242 = (MDP(32) * t231 - MDP(33) * t227) * pkin(5);
t241 = (t232 * MDP(22) - t228 * MDP(23)) * t206;
t240 = -t193 * MDP(13) - MDP(14) * t252;
t239 = t248 + t296;
t222 = t232 ^ 2;
t221 = t228 ^ 2;
t220 = t224 ^ 2;
t210 = -pkin(5) * t232 + t217;
t196 = t226 * t284 + t257;
t195 = t226 * t283 - t213;
t189 = t213 + (-pkin(2) - t283) * t226;
t169 = pkin(5) * t274 + t185;
t141 = pkin(5) * t165 + t144;
t1 = [t220 * t230 ^ 2 * MDP(4) + (t146 ^ 2 + t147 ^ 2 + t181 ^ 2) * MDP(19) + MDP(1) + (t268 * t298 + t270) * t226 + (t193 * MDP(11) + t252 * t292) * t193 + (t261 - 0.2e1 * t262) * t166 + (t152 * t287 + t263) * t153 + ((0.2e1 * t230 * MDP(5) + t266) * t220 + (t226 * MDP(7) + t240) * t298) * t234 + (-t167 * t271 - t189 * t252) * t291 + 0.2e1 * (t168 * t271 + t189 * t193) * MDP(17) + 0.2e1 * (t195 * t226 + t220 * t283) * MDP(9) + 0.2e1 * (-t196 * t226 - t220 * t284) * MDP(10) - t146 * t175 * t290 + (t165 * t289 + t166 * t288) * t144 + (t152 * t286 + t153 * t285) * t141 + (t134 * t286 - t135 * t285 + t139 * t289 - t140 * t288 - t147 * t290 + 0.2e1 * t239 + t297) * t174; t193 * t259 - t153 * t260 + (t193 * t233 + t229 * t252) * MDP(12) + (pkin(2) * t252 - t189 * t233) * MDP(16) + (-pkin(2) * t193 + t189 * t229) * MDP(17) + t270 + (t147 * t187 + t181 * t219) * MDP(19) + (-t141 * t177 + t153 * t169) * MDP(33) + t195 * MDP(9) - t196 * MDP(10) + (t141 * t176 + t152 * t169) * MDP(32) + (t152 * t177 - t153 * t176) * MDP(28) + (t175 * MDP(18) - t146 * MDP(19) + t165 * MDP(25) + t166 * MDP(26)) * t185 + (-t187 * MDP(18) + t243 + t247 + t295) * t174 + (-t147 * MDP(18) + t239 + t244 + t294 + t297) * t205 + (-t146 * MDP(18) + (-t166 * MDP(21) - t174 * MDP(23) + t144 * MDP(25)) * t228 + (t174 * MDP(22) + t144 * MDP(26) + t261 - t262) * t232) * t206 + (t268 + (MDP(7) - t293) * t234) * t224; MDP(8) + pkin(2) * t233 * t291 + (t187 ^ 2 + t219 ^ 2) * MDP(19) + (MDP(20) * t222 - 0.2e1 * t254) * t206 ^ 2 - (t176 * t287 - t260) * t177 + (-0.2e1 * pkin(2) * MDP(17) + t233 * t292 + t259) * t229 + (t176 * t286 - t177 * t285) * t169 + (MDP(19) * t185 + t206 * t290 + t273 * t288 + t274 * t289) * t185 + (t142 * t286 - t143 * t285 + t159 * t289 - t160 * t288 - t187 * t290 + t258 * t205 - 0.2e1 * t171 - 0.2e1 * t172 + 0.2e1 * t241) * t205; -t224 * t266 + t167 * MDP(16) - t168 * MDP(17) + t228 * t261 + (-t165 * t228 + t166 * t232) * MDP(21) + (-t144 * t232 + t165 * t217) * MDP(25) + (t144 * t228 + t166 * t217) * MDP(26) + t209 * t263 + (-t152 * t209 - t153 * t208) * MDP(28) + (t141 * t208 + t152 * t210) * MDP(32) + (t141 * t209 + t153 * t210) * MDP(33) + t236 * t174 + ((-t174 * t223 - t175 * t225) * MDP(18) + (t146 * t225 + t147 * t223) * MDP(19)) * pkin(3) - t240; -t209 * t260 + (-t176 * t209 + t177 * t208) * MDP(28) + (t169 * t208 + t176 * t210) * MDP(32) + (t169 * t209 - t177 * t210) * MDP(33) - t246 * t185 + (MDP(20) * t269 + (-t221 + t222) * MDP(21) + t245 * t217) * t206 + t236 * t205 + ((-t205 * t223 - t206 * t225) * MDP(18) + (-t185 * t225 + t187 * t223) * MDP(19)) * pkin(3) + t293; 0.2e1 * t254 + 0.2e1 * t210 * t197 + t221 * MDP(20) + MDP(15) + (t223 ^ 2 + t225 ^ 2) * MDP(19) * pkin(3) ^ 2 - 0.2e1 * t246 * t217 + (MDP(27) * t209 + t208 * t287 + t210 * t285) * t209; MDP(19) * t181 + t174 * t238; MDP(19) * t219 + t205 * t238; 0; MDP(19); t174 * MDP(24) + (t174 * t280 + t134) * MDP(32) + (-t277 + (-t137 - t282) * t227) * MDP(33) + t248 + t256 + t294; t205 * MDP(24) + (t205 * t280 + t142) * MDP(32) + (-t276 + (-t155 - t281) * t227) * MDP(33) + t241 + t247 + t255; t236; t238; 0.2e1 * t242 + t258; t244 + t256; t243 + t255; t251; t267; MDP(31) + t242; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

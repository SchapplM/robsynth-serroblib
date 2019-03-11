% Calculate joint inertia matrix for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR15_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR15_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRPRR15_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:38:47
% EndTime: 2019-03-09 20:38:55
% DurationCPUTime: 2.70s
% Computational Cost: add. (2726->351), mult. (7067->500), div. (0->0), fcn. (7881->12), ass. (0->145)
t228 = sin(pkin(7));
t234 = sin(qJ(3));
t304 = t228 * t234;
t214 = pkin(10) * t304;
t230 = cos(pkin(7));
t238 = cos(qJ(3));
t312 = pkin(2) * t238;
t266 = -pkin(3) - t312;
t185 = pkin(4) * t304 + t214 + (-pkin(11) + t266) * t230;
t264 = -qJ(4) * t234 - pkin(2);
t314 = pkin(3) + pkin(11);
t193 = (-t314 * t238 + t264) * t228;
t233 = sin(qJ(5));
t237 = cos(qJ(5));
t172 = t185 * t237 - t193 * t233;
t173 = t185 * t233 + t193 * t237;
t303 = t228 * t238;
t205 = t230 * t237 - t233 * t303;
t326 = t205 * MDP(24) + t172 * MDP(27) - t173 * MDP(28);
t199 = (-pkin(3) * t238 + t264) * t228;
t200 = t266 * t230 + t214;
t325 = t200 * MDP(18) - t199 * MDP(20) + t326;
t324 = 0.2e1 * t230;
t231 = cos(pkin(6));
t229 = sin(pkin(6));
t239 = cos(qJ(2));
t299 = t230 * t239;
t267 = t229 * t299;
t235 = sin(qJ(2));
t302 = t229 * t235;
t186 = -t231 * t303 + t234 * t302 - t238 * t267;
t301 = t229 * t239;
t202 = t228 * t301 - t231 * t230;
t323 = -t186 * MDP(12) - t202 * MDP(13);
t232 = sin(qJ(6));
t236 = cos(qJ(6));
t322 = MDP(34) * t232 + MDP(35) * t236;
t209 = pkin(1) * t231 * t235 + pkin(9) * t301;
t183 = (t228 * t231 + t267) * pkin(10) + t209;
t313 = pkin(1) * t239;
t218 = t231 * t313;
t188 = pkin(2) * t231 + t218 + (-pkin(10) * t230 - pkin(9)) * t302;
t197 = (-pkin(10) * t228 * t235 - pkin(2) * t239 - pkin(1)) * t229;
t166 = -t234 * t183 + (t188 * t230 + t197 * t228) * t238;
t320 = 0.2e1 * MDP(20);
t319 = -2 * MDP(25);
t318 = 0.2e1 * MDP(27);
t317 = -2 * MDP(30);
t316 = 0.2e1 * MDP(34);
t315 = 0.2e1 * MDP(35);
t311 = pkin(3) * t202;
t310 = pkin(3) * MDP(21);
t201 = qJ(4) * t202;
t187 = t231 * t304 + (t234 * t299 + t235 * t238) * t229;
t158 = pkin(4) * t187 + t314 * t202 - t166;
t175 = -t188 * t228 + t230 * t197;
t252 = -qJ(4) * t187 + t175;
t159 = t314 * t186 + t252;
t155 = t158 * t237 - t159 * t233;
t153 = -pkin(5) * t187 - t155;
t309 = t153 * t232;
t308 = t153 * t236;
t170 = -pkin(5) * t304 - t172;
t307 = t170 * t232;
t306 = t170 * t236;
t222 = t229 ^ 2;
t305 = t222 * t235;
t300 = t230 * t234;
t298 = t231 * MDP(8);
t297 = t232 * t233;
t296 = t232 * t314;
t295 = t233 * t236;
t294 = t236 * t314;
t208 = pkin(2) * t300 + pkin(10) * t303;
t292 = MDP(28) * t233;
t190 = t205 * t236 + t232 * t304;
t291 = MDP(29) * t190;
t290 = MDP(29) * t236;
t204 = t230 * t233 + t237 * t303;
t289 = MDP(33) * t204;
t286 = qJ(4) * MDP(28);
t177 = t186 * t233 - t202 * t237;
t168 = t177 * t232 - t187 * t236;
t285 = t168 * MDP(32);
t169 = t177 * t236 + t187 * t232;
t284 = t169 * MDP(29);
t283 = t169 * MDP(31);
t176 = -t186 * t237 - t202 * t233;
t282 = t176 * MDP(33);
t281 = t177 * MDP(23);
t280 = t177 * MDP(24);
t278 = t186 * MDP(14);
t189 = t205 * t232 - t236 * t304;
t277 = t189 * MDP(32);
t276 = t190 * MDP(31);
t274 = t202 * MDP(15);
t273 = t205 * MDP(22);
t272 = t205 * MDP(23);
t270 = t205 * MDP(28);
t269 = t314 * MDP(28);
t268 = MDP(11) + MDP(26);
t167 = t238 * t183 + t188 * t300 + t197 * t304;
t220 = t230 * qJ(4);
t198 = -t220 - t208;
t265 = MDP(30) * t232 * t236;
t263 = -MDP(27) * t314 + MDP(24);
t262 = -MDP(25) + t269;
t261 = -t201 + t167;
t192 = pkin(4) * t303 - t198;
t156 = t158 * t233 + t159 * t237;
t206 = t230 * t312 - t214;
t259 = t206 * MDP(16) - t208 * MDP(17);
t258 = -t198 * MDP(18) + t199 * MDP(19);
t256 = MDP(31) * t236 - MDP(32) * t232;
t211 = pkin(5) * t233 - pkin(12) * t237 + qJ(4);
t195 = t236 * t211 + t233 * t296;
t196 = t232 * t211 - t233 * t294;
t255 = MDP(34) * t195 - MDP(35) * t196;
t254 = MDP(34) * t236 - MDP(35) * t232;
t160 = -pkin(4) * t186 + t261;
t251 = MDP(27) * t237 + MDP(18) - t292;
t250 = -MDP(23) + t256;
t249 = (t235 * MDP(6) + t239 * MDP(7)) * t229;
t248 = t155 * MDP(27) - t156 * MDP(28) + t280;
t247 = MDP(27) * qJ(4) + t255;
t244 = t232 * MDP(31) + t236 * MDP(32) - pkin(12) * t322;
t154 = pkin(12) * t187 + t156;
t157 = pkin(5) * t176 - pkin(12) * t177 + t160;
t151 = -t154 * t232 + t157 * t236;
t152 = t154 * t236 + t157 * t232;
t243 = t151 * MDP(34) - t152 * MDP(35) + t282 + t283 - t285;
t171 = pkin(12) * t304 + t173;
t174 = pkin(5) * t204 - pkin(12) * t205 + t192;
t161 = -t171 * t232 + t174 * t236;
t162 = t171 * t236 + t174 * t232;
t242 = t161 * MDP(34) - t162 * MDP(35) + t276 - t277 + t289;
t241 = -MDP(25) + t244;
t227 = t237 ^ 2;
t226 = t236 ^ 2;
t224 = t233 ^ 2;
t223 = t232 ^ 2;
t207 = -pkin(9) * t302 + t218;
t165 = -t166 + t311;
t163 = pkin(3) * t186 + t252;
t1 = [t177 ^ 2 * MDP(22) + (t163 ^ 2 + t165 ^ 2 + t261 ^ 2) * MDP(21) + MDP(1) + (MDP(4) * t235 + 0.2e1 * MDP(5) * t239) * t305 + (t274 + 0.2e1 * t278) * t202 + t268 * t187 ^ 2 + (t168 * t317 + t284) * t169 + t298 * t231 + (t187 * t319 - 0.2e1 * t281 + t282 + 0.2e1 * t283 - 0.2e1 * t285) * t176 + 0.2e1 * (t207 * t231 + t222 * t313) * MDP(9) + 0.2e1 * (-pkin(1) * t305 - t209 * t231) * MDP(10) + (-t163 * t187 - t202 * t261) * t320 + 0.2e1 * (t167 * t202 + t175 * t187) * MDP(17) + 0.2e1 * (-t166 * t202 + t175 * t186) * MDP(16) + 0.2e1 * (-t163 * t186 - t165 * t202) * MDP(19) + 0.2e1 * (-t156 * t187 + t160 * t177) * MDP(28) + 0.2e1 * (t165 * t187 - t186 * t261) * MDP(18) + (t155 * t187 + t160 * t176) * t318 + (t151 * t176 + t153 * t168) * t316 + (-t152 * t176 + t153 * t169) * t315 + 0.2e1 * t249 * t231 + 0.2e1 * (t280 + t323) * t187; (t166 * t230 - t202 * t206) * MDP(16) + (t160 * t205 + t177 * t192) * MDP(28) + (t160 * t204 + t176 * t192) * MDP(27) + (-t167 * t230 + t202 * t208) * MDP(17) + t298 + (-t176 * t205 - t177 * t204) * MDP(23) + t207 * MDP(9) - t209 * MDP(10) + (t151 * t204 + t153 * t189 + t161 * t176 + t168 * t170) * MDP(34) + (-t168 * t204 - t176 * t189) * MDP(32) + (t169 * t204 + t176 * t190) * MDP(31) + (-t152 * t204 + t153 * t190 - t162 * t176 + t169 * t170) * MDP(35) + (t163 * t199 + t165 * t200 - t198 * t261) * MDP(21) + (-t168 * t190 - t169 * t189) * MDP(30) - t230 * t274 + t204 * t282 + t177 * t273 + t190 * t284 + (t165 * t230 - t200 * t202) * MDP(19) + (t198 * t202 + t230 * t261) * MDP(20) + t249 + (-t230 * MDP(14) - t258) * t186 + (MDP(13) * t230 - t204 * MDP(25) + t325) * t187 + ((-t186 * MDP(16) - t187 * MDP(17)) * pkin(2) + (t187 * MDP(12) - t202 * MDP(14) - t175 * MDP(16) + MDP(18) * t261 + t163 * MDP(19)) * t238 + (t175 * MDP(17) + t165 * MDP(18) - t163 * MDP(20) - t176 * MDP(25) + t268 * t187 + t248 + t323) * t234) * t228; t230 ^ 2 * MDP(15) + t205 ^ 2 * MDP(22) + (t198 ^ 2 + t199 ^ 2 + t200 ^ 2) * MDP(21) + MDP(8) + (t189 * t317 + t291) * t190 + (t304 * t319 - 0.2e1 * t272 + 0.2e1 * t276 - 0.2e1 * t277 + t289) * t204 + (-t162 * t204 + t170 * t190) * t315 + (t161 * t204 + t170 * t189) * t316 + (t200 * MDP(19) - t198 * MDP(20) + t259) * t324 + 0.2e1 * (t204 * MDP(27) + t270) * t192 + ((t234 * MDP(13) + t238 * MDP(14)) * t324 + 0.2e1 * t258 * t238 + (t268 * t234 ^ 2 + 0.2e1 * (MDP(16) * t238 - MDP(17) * t234) * pkin(2)) * t228 + 0.2e1 * (MDP(12) * t303 + t325) * t234) * t228; t187 * MDP(13) - t278 - t274 + t166 * MDP(16) - t167 * MDP(17) + (-pkin(3) * t187 - qJ(4) * t186) * MDP(18) + (-t166 + 0.2e1 * t311) * MDP(19) + (t261 - t201) * MDP(20) + (-pkin(3) * t165 + qJ(4) * t261) * MDP(21) + t177 * t286 + t247 * t176 + (t160 * MDP(27) + t262 * t187 + t243 - t281) * t233 + (t177 * MDP(22) + t160 * MDP(28) + t236 * t284 + (-t168 * t236 - t169 * t232) * MDP(30) + (t168 * t314 + t309) * MDP(34) + (t169 * t314 + t308) * MDP(35) + t263 * t187 + t250 * t176) * t237; t214 * MDP(19) + (0.2e1 * t220 + t208) * MDP(20) + (-pkin(3) * t200 - qJ(4) * t198) * MDP(21) + qJ(4) * t270 + (MDP(15) + (-0.2e1 * pkin(3) - t312) * MDP(19)) * t230 + t247 * t204 + (t192 * MDP(27) + t242 - t272) * t233 + ((MDP(18) * qJ(4) + MDP(14)) * t238 + (-pkin(3) * MDP(18) + t262 * t233 + MDP(13)) * t234) * t228 + (t273 + t192 * MDP(28) + t190 * t290 + (-t189 * t236 - t190 * t232) * MDP(30) + (t189 * t314 + t307) * MDP(34) + (t190 * t314 + t306) * MDP(35) + t263 * t304 + t250 * t204) * t237 + t259; t224 * MDP(33) + MDP(15) + (-0.2e1 * MDP(19) + t310) * pkin(3) + (MDP(21) * qJ(4) + t233 * t318 + t320) * qJ(4) + (MDP(29) * t226 + MDP(22) - 0.2e1 * t265) * t227 + (t195 * t233 + t227 * t296) * t316 + (-t196 * t233 + t227 * t294) * t315 + 0.2e1 * (t233 * t250 + t286) * t237; -t202 * MDP(19) + t165 * MDP(21) + (-t168 * t237 - t176 * t297) * MDP(34) + (-t169 * t237 - t176 * t295) * MDP(35) + t251 * t187; t230 * MDP(19) + t200 * MDP(21) + (-t189 * t237 - t204 * t297) * MDP(34) + (-t190 * t237 - t204 * t295) * MDP(35) + t251 * t304; MDP(19) - t310 + t322 * (-t224 - t227); MDP(21); t187 * MDP(26) + t232 * t284 + (-t168 * t232 + t169 * t236) * MDP(30) + (-pkin(5) * t168 - t308) * MDP(34) + (-pkin(5) * t169 + t309) * MDP(35) + t241 * t176 + t248; MDP(26) * t304 + t232 * t291 + (-t189 * t232 + t190 * t236) * MDP(30) + (-pkin(5) * t189 - t306) * MDP(34) + (-pkin(5) * t190 + t307) * MDP(35) + t241 * t204 + t326; (t241 + t269) * t233 + (t232 * t290 + (-t223 + t226) * MDP(30) + (-pkin(5) * t232 - t294) * MDP(34) + (-pkin(5) * t236 + t296) * MDP(35) + t263) * t237; -t292 + (MDP(27) + t254) * t237; MDP(29) * t223 + 0.2e1 * pkin(5) * t254 + MDP(26) + 0.2e1 * t265; t243; t242; MDP(33) * t233 + t256 * t237 + t255; -t322 * t233; t244; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

% Calculate joint inertia matrix for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_inertiaJ_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR12_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:00:49
% EndTime: 2019-03-09 08:00:56
% DurationCPUTime: 2.22s
% Computational Cost: add. (5772->346), mult. (16079->526), div. (0->0), fcn. (19033->16), ass. (0->164)
t243 = sin(pkin(8));
t339 = 0.2e1 * t243;
t250 = sin(qJ(6));
t254 = cos(qJ(6));
t262 = -(MDP(34) * t250 + MDP(35) * t254) * pkin(13) + t250 * MDP(31) + t254 * MDP(32);
t242 = sin(pkin(14));
t245 = sin(pkin(6));
t246 = cos(pkin(14));
t310 = t245 * t246;
t249 = cos(pkin(6));
t330 = pkin(1) * t249;
t222 = qJ(2) * t310 + t242 * t330;
t244 = sin(pkin(7));
t248 = cos(pkin(7));
t309 = t246 * t248;
t271 = t244 * t249 + t245 * t309;
t201 = t271 * pkin(10) + t222;
t210 = (-pkin(10) * t242 * t244 - pkin(2) * t246 - pkin(1)) * t245;
t253 = sin(qJ(3));
t257 = cos(qJ(3));
t234 = t246 * t330;
t315 = t242 * t245;
t203 = pkin(2) * t249 + t234 + (-pkin(10) * t248 - qJ(2)) * t315;
t316 = t203 * t248;
t182 = t257 * t201 + (t210 * t244 + t316) * t253;
t219 = -t244 * t310 + t248 * t249;
t247 = cos(pkin(8));
t231 = t253 * t315;
t267 = -t271 * t257 + t231;
t266 = t267 * t247;
t258 = t219 * t243 - t266;
t173 = t258 * pkin(11) + t182;
t252 = sin(qJ(4));
t256 = cos(qJ(4));
t311 = t244 * t257;
t181 = -t201 * t253 + t210 * t311 + t257 * t316;
t312 = t244 * t253;
t202 = t249 * t312 + (t242 * t257 + t253 * t309) * t245;
t327 = pkin(11) * t202;
t176 = pkin(3) * t219 - t247 * t327 + t181;
t188 = -t244 * t203 + t248 * t210;
t178 = t267 * pkin(3) - t243 * t327 + t188;
t278 = t176 * t247 + t178 * t243;
t167 = -t252 * t173 + t278 * t256;
t338 = 2 * MDP(20);
t337 = 2 * MDP(21);
t336 = 2 * MDP(27);
t335 = 2 * MDP(28);
t334 = -2 * MDP(30);
t333 = 0.2e1 * MDP(34);
t332 = 0.2e1 * MDP(35);
t238 = t245 ^ 2;
t331 = pkin(1) * t238;
t329 = pkin(3) * t252;
t328 = pkin(3) * t256;
t326 = pkin(12) * t250;
t325 = pkin(12) * t254;
t255 = cos(qJ(5));
t324 = pkin(12) * t255;
t323 = MDP(27) * pkin(4);
t322 = pkin(4) * MDP(28);
t168 = t173 * t256 + t278 * t252;
t190 = -t219 * t247 - t267 * t243;
t165 = -pkin(12) * t190 + t168;
t169 = -t176 * t243 + t247 * t178;
t313 = t243 * t256;
t184 = t202 * t252 - t219 * t313 + t256 * t266;
t185 = t202 * t256 + t258 * t252;
t166 = pkin(4) * t184 - pkin(12) * t185 + t169;
t251 = sin(qJ(5));
t161 = -t165 * t251 + t166 * t255;
t159 = -pkin(5) * t184 - t161;
t321 = t159 * t250;
t320 = t159 * t254;
t308 = t247 * t257;
t314 = t243 * t252;
t205 = t248 * t314 + (t252 * t308 + t253 * t256) * t244;
t220 = -t243 * t311 + t247 * t248;
t192 = t205 * t251 - t255 * t220;
t319 = t192 * t251;
t282 = pkin(11) * t313;
t217 = t282 + (pkin(12) + t329) * t247;
t218 = (-pkin(4) * t256 - pkin(12) * t252 - pkin(3)) * t243;
t197 = -t217 * t251 + t218 * t255;
t195 = pkin(5) * t313 - t197;
t318 = t195 * t250;
t317 = t195 * t254;
t307 = MDP(18) * t247;
t306 = MDP(26) * t256;
t224 = t247 * t251 + t255 * t314;
t208 = t224 * t254 - t250 * t313;
t305 = MDP(29) * t208;
t304 = MDP(29) * t250;
t303 = MDP(29) * t254;
t302 = MDP(33) * t255;
t175 = t185 * t255 - t190 * t251;
t170 = t175 * t250 - t184 * t254;
t301 = t170 * MDP(32);
t171 = t175 * t254 + t184 * t250;
t300 = t171 * MDP(31);
t174 = t185 * t251 + t190 * t255;
t299 = t174 * MDP(33);
t298 = t175 * MDP(23);
t297 = t175 * MDP(24);
t296 = t184 * MDP(26);
t295 = t185 * MDP(15);
t294 = t185 * MDP(16);
t293 = t190 * MDP(17);
t292 = t190 * MDP(18);
t207 = t224 * t250 + t254 * t313;
t291 = t207 * MDP(32);
t290 = t208 * MDP(31);
t289 = t219 * MDP(12);
t223 = -t255 * t247 + t251 * t314;
t288 = t223 * MDP(33);
t287 = t224 * MDP(22);
t286 = t224 * MDP(23);
t285 = t224 * MDP(24);
t284 = t247 * MDP(19);
t283 = t252 * MDP(17);
t281 = MDP(30) * t250 * t254;
t280 = pkin(12) * MDP(27) - MDP(24);
t279 = pkin(12) * MDP(28) - MDP(25);
t162 = t165 * t255 + t166 * t251;
t198 = t217 * t255 + t218 * t251;
t235 = pkin(11) * t314;
t225 = t247 * t328 - t235;
t226 = t247 * t329 + t282;
t277 = -t225 * MDP(20) + t226 * MDP(21);
t276 = MDP(31) * t254 - MDP(32) * t250;
t228 = -pkin(5) * t255 - pkin(13) * t251 - pkin(4);
t212 = t228 * t254 - t250 * t324;
t213 = t228 * t250 + t254 * t324;
t274 = MDP(34) * t212 - MDP(35) * t213;
t273 = MDP(34) * t254 - MDP(35) * t250;
t216 = t235 + (-pkin(4) - t328) * t247;
t270 = -MDP(23) + t276;
t269 = t274 - t323;
t268 = t197 * MDP(27) - t198 * MDP(28) + t285;
t164 = pkin(4) * t190 - t167;
t265 = t185 * MDP(17) - t190 * MDP(19) + t167 * MDP(20) - t168 * MDP(21);
t264 = t161 * MDP(27) - t162 * MDP(28) + t296 + t297;
t263 = t267 * MDP(11);
t160 = pkin(13) * t184 + t162;
t163 = pkin(5) * t174 - pkin(13) * t175 + t164;
t157 = -t160 * t250 + t163 * t254;
t158 = t160 * t254 + t163 * t250;
t261 = t157 * MDP(34) - t158 * MDP(35) + t299 + t300 - t301;
t194 = pkin(5) * t223 - pkin(13) * t224 + t216;
t196 = -pkin(13) * t313 + t198;
t179 = t194 * t254 - t196 * t250;
t180 = t194 * t250 + t196 * t254;
t260 = t179 * MDP(34) - t180 * MDP(35) + t288 + t290 - t291;
t259 = -MDP(25) + t262;
t241 = t254 ^ 2;
t240 = t251 ^ 2;
t239 = t250 ^ 2;
t237 = t243 ^ 2;
t221 = -qJ(2) * t315 + t234;
t204 = -t256 * t244 * t308 - t248 * t313 + t252 * t312;
t193 = t205 * t255 + t220 * t251;
t187 = t193 * t254 + t204 * t250;
t186 = -t193 * t250 + t204 * t254;
t1 = [t190 ^ 2 * MDP(19) + MDP(1) + t175 ^ 2 * MDP(22) + (pkin(1) ^ 2 * t238 + t221 ^ 2 + t222 ^ 2) * MDP(7) + (-0.2e1 * t263 + t289) * t219 + (-0.2e1 * t293 + t295) * t185 + (MDP(29) * t171 + t170 * t334) * t171 + (0.2e1 * t219 * MDP(10) + MDP(8) * t202 - 0.2e1 * t267 * MDP(9)) * t202 + (0.2e1 * t292 - 0.2e1 * t294 + t296 + 0.2e1 * t297) * t184 + (-0.2e1 * t184 * MDP(25) - 0.2e1 * t298 + t299 + 0.2e1 * t300 - 0.2e1 * t301) * t174 + (t157 * t174 + t159 * t170) * t333 + (-t158 * t174 + t159 * t171) * t332 + (t161 * t184 + t164 * t174) * t336 + (-t162 * t184 + t164 * t175) * t335 + (t168 * t190 + t169 * t185) * t337 + (-t167 * t190 + t169 * t184) * t338 + 0.2e1 * (-t182 * t219 + t188 * t202) * MDP(14) + 0.2e1 * (-t221 * t242 + t222 * t246) * t245 * MDP(6) + 0.2e1 * (-t222 * t249 - t242 * t331) * MDP(5) + 0.2e1 * (t221 * t249 + t246 * t331) * MDP(4) + 0.2e1 * (t181 * t219 + t188 * t267) * MDP(13); -MDP(4) * t310 + MDP(5) * t315 - t245 * pkin(1) * MDP(7) + (t248 * t231 + (t244 * t219 - t248 * t271) * t257) * MDP(13) + (t202 * t248 - t219 * t312) * MDP(14) + (t184 * t220 + t190 * t204) * MDP(20) + (t185 * t220 + t190 * t205) * MDP(21) + (t174 * t204 - t184 * t192) * MDP(27) + (t175 * t204 - t184 * t193) * MDP(28) + (t170 * t192 + t174 * t186) * MDP(34) + (t171 * t192 - t174 * t187) * MDP(35); MDP(7); t171 * t305 + t174 * t288 + t175 * t287 + t181 * MDP(13) - t182 * MDP(14) + t202 * MDP(10) + (-t170 * t208 - t171 * t207) * MDP(30) + t289 + (t157 * t223 + t159 * t207 + t170 * t195 + t174 * t179) * MDP(34) + (-t170 * t223 - t174 * t207) * MDP(32) + (t171 * t223 + t174 * t208) * MDP(31) + (-t158 * t223 + t159 * t208 + t171 * t195 - t174 * t180) * MDP(35) + (-t174 * t224 - t175 * t223) * MDP(23) + (t164 * t223 + t174 * t216) * MDP(27) + (t164 * t224 + t175 * t216) * MDP(28) - t263 + t277 * t190 + t265 * t247 + (-t223 * MDP(25) + t268 - t307) * t184 + ((-t184 * MDP(20) - t185 * MDP(21)) * pkin(3) + (-t184 * MDP(16) + t169 * MDP(21) - t293 + t295) * t252 + (-t169 * MDP(20) + t174 * MDP(25) - t264 - t292 + t294) * t256) * t243; (-t204 * t247 - t220 * t313) * MDP(20) + (-t205 * t247 + t220 * t314) * MDP(21) + (t192 * t313 + t204 * t223) * MDP(27) + (t193 * t313 + t204 * t224) * MDP(28) + (t186 * t223 + t192 * t207) * MDP(34) + (-t187 * t223 + t192 * t208) * MDP(35) + (MDP(13) * t257 - MDP(14) * t253) * t244; t237 * t252 ^ 2 * MDP(15) + t224 ^ 2 * MDP(22) + MDP(12) + (t283 * t339 + t284) * t247 + (t207 * t334 + t305) * t208 + ((-t285 + t307) * t339 + (0.2e1 * MDP(16) * t252 + t306) * t237) * t256 + (0.2e1 * MDP(25) * t313 - 0.2e1 * t286 + t288 + 0.2e1 * t290 - 0.2e1 * t291) * t223 + (t225 * t247 + t237 * t328) * t338 + (-t226 * t247 - t237 * t329) * t337 + (-t197 * t313 + t216 * t223) * t336 + (t198 * t313 + t216 * t224) * t335 + (t179 * t223 + t195 * t207) * t333 + (-t180 * t223 + t195 * t208) * t332; -t175 * t322 - t184 * MDP(18) + t269 * t174 + (-t164 * MDP(27) - t279 * t184 - t261 + t298) * t255 + (t175 * MDP(22) + t164 * MDP(28) + t171 * t303 + (-t170 * t254 - t171 * t250) * MDP(30) + (pkin(12) * t170 + t321) * MDP(34) + (pkin(12) * t171 + t320) * MDP(35) - t280 * t184 + t270 * t174) * t251 + t265; -t205 * MDP(21) + (-t186 * t255 + t250 * t319) * MDP(34) + (t187 * t255 + t254 * t319) * MDP(35) + (-MDP(27) * t255 + MDP(28) * t251 - MDP(20)) * t204; -t224 * t322 + t284 + (t256 * MDP(18) + t283) * t243 + t269 * t223 + (-t216 * MDP(27) + t279 * t313 - t260 + t286) * t255 + (t287 + t216 * MDP(28) + t208 * t303 + (-t207 * t254 - t208 * t250) * MDP(30) + (pkin(12) * t207 + t318) * MDP(34) + (pkin(12) * t208 + t317) * MDP(35) + t280 * t313 + t270 * t223) * t251 - t277; MDP(19) + (t302 + (2 * t323)) * t255 + (MDP(29) * t241 + MDP(22) - 0.2e1 * t281) * t240 + (-t212 * t255 + t240 * t326) * t333 + (t213 * t255 + t240 * t325) * t332 + 0.2e1 * (-t255 * t270 - t322) * t251; t171 * t304 + (-t170 * t250 + t171 * t254) * MDP(30) + (-pkin(5) * t170 - t320) * MDP(34) + (-pkin(5) * t171 + t321) * MDP(35) + t259 * t174 + t264; -MDP(28) * t193 + (-MDP(27) - t273) * t192; -t243 * t306 + t208 * t304 + (-t207 * t250 + t208 * t254) * MDP(30) + (-pkin(5) * t207 - t317) * MDP(34) + (-pkin(5) * t208 + t318) * MDP(35) + t259 * t223 + t268; (-t262 - t279) * t255 + (t250 * t303 + (-t239 + t241) * MDP(30) + (-pkin(5) * t250 - t325) * MDP(34) + (-pkin(5) * t254 + t326) * MDP(35) - t280) * t251; MDP(29) * t239 + 0.2e1 * pkin(5) * t273 + MDP(26) + 0.2e1 * t281; t261; MDP(34) * t186 - MDP(35) * t187; t260; t276 * t251 + t274 - t302; t262; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

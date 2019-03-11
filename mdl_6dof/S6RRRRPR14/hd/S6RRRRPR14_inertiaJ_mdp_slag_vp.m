% Calculate joint inertia matrix for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR14_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR14_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR14_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:25:48
% EndTime: 2019-03-10 00:25:55
% DurationCPUTime: 2.79s
% Computational Cost: add. (5551->424), mult. (14344->621), div. (0->0), fcn. (16375->14), ass. (0->171)
t267 = sin(pkin(7));
t364 = 0.2e1 * t267;
t363 = (MDP(28) * qJ(5));
t276 = cos(qJ(4));
t362 = t276 * MDP(19);
t266 = sin(pkin(13));
t269 = cos(pkin(13));
t272 = sin(qJ(6));
t351 = cos(qJ(6));
t246 = t266 * t351 + t272 * t269;
t344 = pkin(12) + qJ(5);
t250 = t344 * t266;
t251 = t344 * t269;
t292 = -t272 * t266 + t269 * t351;
t283 = t246 * MDP(31) + t292 * MDP(32) + (-t250 * t351 - t272 * t251) * MDP(34) - (-t272 * t250 + t251 * t351) * MDP(35);
t293 = MDP(25) * t266 + MDP(26) * t269;
t280 = -qJ(5) * t293 - MDP(21) + t283;
t271 = cos(pkin(6));
t275 = sin(qJ(2));
t268 = sin(pkin(6));
t278 = cos(qJ(2));
t335 = t268 * t278;
t242 = pkin(1) * t271 * t275 + pkin(9) * t335;
t270 = cos(pkin(7));
t334 = t270 * t278;
t300 = t268 * t334;
t215 = (t267 * t271 + t300) * pkin(10) + t242;
t274 = sin(qJ(3));
t277 = cos(qJ(3));
t350 = pkin(1) * t278;
t257 = t271 * t350;
t336 = t268 * t275;
t218 = pkin(2) * t271 + t257 + (-pkin(10) * t270 - pkin(9)) * t336;
t226 = (-pkin(10) * t267 * t275 - pkin(2) * t278 - pkin(1)) * t268;
t296 = t218 * t270 + t226 * t267;
t191 = -t274 * t215 + t277 * t296;
t361 = 2 * MDP(16);
t360 = 2 * MDP(17);
t359 = 2 * MDP(23);
t358 = 2 * MDP(24);
t357 = 0.2e1 * MDP(25);
t356 = 0.2e1 * MDP(26);
t355 = 2 * MDP(27);
t354 = -2 * MDP(30);
t353 = 0.2e1 * MDP(34);
t352 = 0.2e1 * MDP(35);
t349 = pkin(2) * t274;
t348 = pkin(2) * t277;
t347 = pkin(11) * t266;
t346 = pkin(11) * t276;
t273 = sin(qJ(4));
t345 = pkin(12) * t273;
t343 = MDP(23) * pkin(3);
t342 = pkin(3) * MDP(24);
t341 = pkin(4) * MDP(28);
t340 = pkin(11) * MDP(24);
t263 = t268 ^ 2;
t339 = t263 * t275;
t338 = t267 * t274;
t337 = t267 * t277;
t333 = t271 * MDP(8);
t201 = -t218 * t267 + t270 * t226;
t216 = -t271 * t337 + t274 * t336 - t277 * t300;
t217 = t271 * t338 + (t274 * t334 + t275 * t277) * t268;
t185 = pkin(3) * t216 - pkin(11) * t217 + t201;
t192 = t215 * t277 + t274 * t296;
t235 = t267 * t335 - t270 * t271;
t188 = -pkin(11) * t235 + t192;
t177 = t185 * t273 + t188 * t276;
t174 = qJ(5) * t216 + t177;
t187 = pkin(3) * t235 - t191;
t206 = t217 * t273 + t235 * t276;
t207 = t217 * t276 - t235 * t273;
t180 = pkin(4) * t206 - qJ(5) * t207 + t187;
t169 = t269 * t174 + t266 * t180;
t254 = pkin(10) * t338;
t230 = t254 + (-pkin(3) - t348) * t270;
t237 = -t276 * t270 + t273 * t338;
t238 = t270 * t273 + t276 * t338;
t205 = pkin(4) * t237 - qJ(5) * t238 + t230;
t301 = pkin(10) * t337;
t231 = t301 + (pkin(11) + t349) * t270;
t232 = (-pkin(3) * t277 - pkin(11) * t274 - pkin(2)) * t267;
t211 = t231 * t276 + t232 * t273;
t208 = -qJ(5) * t337 + t211;
t190 = t266 * t205 + t269 * t208;
t249 = -pkin(4) * t276 - qJ(5) * t273 - pkin(3);
t228 = t266 * t249 + t269 * t346;
t331 = MDP(11) * t217;
t330 = MDP(13) * t274;
t329 = MDP(14) * t270;
t328 = MDP(22) * t277;
t327 = MDP(25) * t269;
t326 = MDP(26) * t266;
t176 = t185 * t276 - t273 * t188;
t175 = -pkin(4) * t216 - t176;
t325 = MDP(28) * t175;
t210 = -t273 * t231 + t232 * t276;
t209 = pkin(4) * t337 - t210;
t324 = MDP(28) * t209;
t219 = t238 * t266 + t269 * t337;
t220 = t238 * t269 - t266 * t337;
t200 = -t272 * t219 + t220 * t351;
t323 = MDP(29) * t200;
t234 = t292 * t273;
t322 = MDP(29) * t234;
t321 = MDP(29) * t246;
t320 = MDP(31) * t200;
t319 = MDP(31) * t234;
t199 = t219 * t351 + t220 * t272;
t318 = MDP(32) * t199;
t233 = t246 * t273;
t317 = MDP(32) * t233;
t316 = MDP(33) * t237;
t315 = MDP(33) * t276;
t314 = MDP(34) * t292;
t194 = t207 * t266 - t216 * t269;
t195 = t207 * t269 + t216 * t266;
t181 = t194 * t351 + t195 * t272;
t313 = t181 * MDP(32);
t182 = -t272 * t194 + t195 * t351;
t312 = t182 * MDP(31);
t311 = t206 * MDP(33);
t310 = t207 * MDP(19);
t309 = t207 * MDP(20);
t308 = t216 * MDP(22);
t307 = t217 * MDP(12);
t306 = t235 * MDP(13);
t305 = t235 * MDP(14);
t304 = t238 * MDP(18);
t303 = t238 * MDP(20);
t302 = t270 * MDP(15);
t299 = -MDP(21) + t340;
t168 = -t174 * t266 + t269 * t180;
t189 = t269 * t205 - t208 * t266;
t298 = -t168 * t266 + t169 * t269;
t297 = -t189 * t266 + t190 * t269;
t239 = t270 * t348 - t254;
t241 = t270 * t349 + t301;
t294 = -t239 * MDP(16) + t241 * MDP(17);
t291 = (MDP(6) * t275 + MDP(7) * t278) * t268;
t290 = t326 - t327 - t341;
t289 = t210 * MDP(23) - t211 * MDP(24) + t303;
t288 = t194 * MDP(25) + t195 * MDP(26) + t325;
t287 = MDP(25) * t219 + MDP(26) * t220 + t324;
t286 = t217 * MDP(13) - t235 * MDP(15) + t191 * MDP(16) - t192 * MDP(17);
t285 = t176 * MDP(23) - t177 * MDP(24) + t308 + t309;
t244 = t269 * t249;
t213 = -t269 * t345 + t244 + (-pkin(5) - t347) * t276;
t221 = -t266 * t345 + t228;
t197 = t213 * t351 - t272 * t221;
t198 = t272 * t213 + t221 * t351;
t284 = MDP(34) * t197 - MDP(35) * t198 - t317 + t319;
t166 = pkin(5) * t206 - pkin(12) * t195 + t168;
t167 = -pkin(12) * t194 + t169;
t164 = t166 * t351 - t272 * t167;
t165 = t272 * t166 + t167 * t351;
t281 = t164 * MDP(34) - t165 * MDP(35) + t311 + t312 - t313;
t265 = t273 ^ 2;
t262 = t267 ^ 2;
t260 = -pkin(5) * t269 - pkin(4);
t247 = (pkin(5) * t266 + pkin(11)) * t273;
t240 = -pkin(9) * t336 + t257;
t227 = -t266 * t346 + t244;
t196 = pkin(5) * t219 + t209;
t184 = -pkin(12) * t219 + t190;
t183 = pkin(5) * t237 - pkin(12) * t220 + t189;
t172 = t272 * t183 + t184 * t351;
t171 = t183 * t351 - t272 * t184;
t170 = pkin(5) * t194 + t175;
t1 = [t235 ^ 2 * MDP(15) + (t168 ^ 2 + t169 ^ 2 + t175 ^ 2) * MDP(28) + MDP(1) + t207 ^ 2 * MDP(18) + (MDP(4) * t275 + 0.2e1 * MDP(5) * t278) * t339 + (-0.2e1 * t306 + t331) * t217 + (MDP(29) * t182 + t181 * t354) * t182 + (0.2e1 * t291 + t333) * t271 + (0.2e1 * t305 - 0.2e1 * t307 + t308 + 0.2e1 * t309) * t216 + (-0.2e1 * t216 * MDP(21) - 0.2e1 * t310 + t311 + 0.2e1 * t312 - 0.2e1 * t313) * t206 + (-t168 * t195 - t169 * t194) * t355 + (t164 * t206 + t170 * t181) * t353 + (t168 * t206 + t175 * t194) * t357 + (-t165 * t206 + t170 * t182) * t352 + (-t169 * t206 + t175 * t195) * t356 + (t176 * t216 + t187 * t206) * t359 + (-t177 * t216 + t187 * t207) * t358 + (t192 * t235 + t201 * t217) * t360 + (-t191 * t235 + t201 * t216) * t361 + 0.2e1 * (-pkin(1) * t339 - t242 * t271) * MDP(10) + 0.2e1 * (t240 * t271 + t263 * t350) * MDP(9); t182 * t323 + t237 * t311 + t207 * t304 + (t187 * t237 + t206 * t230) * MDP(23) + (t187 * t238 + t207 * t230) * MDP(24) + (-t181 * t200 - t182 * t199) * MDP(30) + (t168 * t189 + t169 * t190 + t175 * t209) * MDP(28) + (-t168 * t220 - t169 * t219 - t189 * t195 - t190 * t194) * MDP(27) + (t164 * t237 + t170 * t199 + t171 * t206 + t181 * t196) * MDP(34) + (-t181 * t237 - t199 * t206) * MDP(32) + (t182 * t237 + t200 * t206) * MDP(31) + (t168 * t237 + t175 * t219 + t189 * t206 + t194 * t209) * MDP(25) + (-t165 * t237 + t170 * t200 - t172 * t206 + t182 * t196) * MDP(35) + (-t169 * t237 + t175 * t220 - t190 * t206 + t195 * t209) * MDP(26) + (-t206 * t238 - t207 * t237) * MDP(19) + t240 * MDP(9) - t242 * MDP(10) + t333 + t294 * t235 + t286 * t270 + t291 + (-t237 * MDP(21) + t289 - t329) * t216 + ((-t216 * MDP(16) - t217 * MDP(17)) * pkin(2) + (-t216 * MDP(12) + MDP(17) * t201 - t306 + t331) * t274 + (-t201 * MDP(16) + t206 * MDP(21) - t285 - t305 + t307) * t277) * t267; t262 * t274 ^ 2 * MDP(11) + MDP(8) + (t189 ^ 2 + t190 ^ 2 + t209 ^ 2) * MDP(28) + t238 ^ 2 * MDP(18) + (t330 * t364 + t302) * t270 + (t199 * t354 + t323) * t200 + ((-t303 + t329) * t364 + (0.2e1 * MDP(12) * t274 + t328) * t262) * t277 + (-0.2e1 * MDP(19) * t238 + 0.2e1 * MDP(21) * t337 + t316 - 0.2e1 * t318 + 0.2e1 * t320) * t237 + (-t189 * t220 - t190 * t219) * t355 + (t171 * t237 + t196 * t199) * t353 + (t189 * t237 + t209 * t219) * t357 + (-t172 * t237 + t196 * t200) * t352 + (-t190 * t237 + t209 * t220) * t356 + (-t241 * t270 - t262 * t349) * t360 + (-t210 * t337 + t230 * t237) * t359 + (t211 * t337 + t230 * t238) * t358 + (t239 * t270 + t262 * t348) * t361; -t216 * MDP(14) - t207 * t342 + (-t194 * t228 - t195 * t227) * MDP(27) + (t168 * t227 + t169 * t228) * MDP(28) + t182 * t322 + (-t181 * t234 - t182 * t233) * MDP(30) + (t170 * t233 + t181 * t247) * MDP(34) + (t170 * t234 + t182 * t247) * MDP(35) + (MDP(25) * t227 - MDP(26) * t228 + t284 - t343) * t206 + (-t187 * MDP(23) - t168 * MDP(25) + t169 * MDP(26) - t216 * t299 - t281 + t310) * t276 + (t207 * MDP(18) - t206 * MDP(19) + t216 * MDP(20) + t187 * MDP(24) + (-t168 * t269 - t169 * t266) * MDP(27) + t293 * t175 + (-MDP(23) * t216 + t288) * pkin(11)) * t273 + t286; t302 + (-pkin(3) * t237 - t230 * t276) * MDP(23) + (-t189 * t276 + t227 * t237) * MDP(25) + (t190 * t276 - t228 * t237) * MDP(26) + (-t219 * t228 - t220 * t227) * MDP(27) + (t189 * t227 + t190 * t228) * MDP(28) + t200 * t322 + (-t199 * t234 - t200 * t233) * MDP(30) + (-t200 * t276 + t234 * t237) * MDP(31) + (t199 * t276 - t233 * t237) * MDP(32) - t237 * t315 + (-t171 * t276 + t196 * t233 + t197 * t237 + t199 * t247) * MDP(34) + (t172 * t276 + t196 * t234 - t198 * t237 + t200 * t247) * MDP(35) + (-t342 + t362) * t238 + (t330 + (t299 * t276 + MDP(14)) * t277) * t267 + (t304 - t237 * MDP(19) - MDP(20) * t337 + t230 * MDP(24) + (-t189 * t269 - t190 * t266) * MDP(27) + t293 * t209 + (MDP(23) * t337 + t287) * pkin(11)) * t273 - t294; MDP(15) + t265 * MDP(18) + (pkin(11) ^ 2 * t265 + t227 ^ 2 + t228 ^ 2) * MDP(28) + (t233 * t354 + t322) * t234 + (t315 + 0.2e1 * t317 - 0.2e1 * t319 + (2 * t343)) * t276 + (-t227 * t276 + t265 * t347) * t357 + (pkin(11) * t265 * t269 + t228 * t276) * t356 + (-t197 * t276 + t233 * t247) * t353 + (t198 * t276 + t234 * t247) * t352 + (-(2 * t342) + 0.2e1 * t362 + (-t227 * t269 - t228 * t266) * t355) * t273; (-pkin(4) * t194 - t175 * t269) * MDP(25) + (-pkin(4) * t195 + t175 * t266) * MDP(26) + t298 * MDP(27) - pkin(4) * t325 + t182 * t321 + (-t181 * t246 + t182 * t292) * MDP(30) + (-t170 * t292 + t181 * t260) * MDP(34) + (t170 * t246 + t182 * t260) * MDP(35) + ((-t194 * t269 + t195 * t266) * MDP(27) + t298 * MDP(28)) * qJ(5) + t280 * t206 + t285; -t267 * t328 + (-pkin(4) * t219 - t209 * t269) * MDP(25) + (-pkin(4) * t220 + t209 * t266) * MDP(26) + t297 * MDP(27) - pkin(4) * t324 + t200 * t321 + (-t199 * t246 + t200 * t292) * MDP(30) + (-t196 * t292 + t199 * t260) * MDP(34) + (t196 * t246 + t200 * t260) * MDP(35) + ((-t219 * t269 + t220 * t266) * MDP(27) + t297 * MDP(28)) * qJ(5) + t280 * t237 + t289; t234 * t321 + (-t233 * t246 + t234 * t292) * MDP(30) + (t233 * t260 - t247 * t292) * MDP(34) + (t234 * t260 + t246 * t247) * MDP(35) + (-t280 - t340) * t276 + (MDP(20) - t293 * pkin(4) + (-MDP(23) + t290) * pkin(11)) * t273 + (MDP(27) + t363) * (-t227 * t266 + t228 * t269); -0.2e1 * t260 * t314 + MDP(22) + (-0.2e1 * t326 + 0.2e1 * t327 + t341) * pkin(4) + (t260 * t352 - t292 * t354 + t321) * t246 + (t355 + t363) * (t266 ^ 2 + t269 ^ 2) * qJ(5); t181 * MDP(34) + t182 * MDP(35) + t288; MDP(34) * t199 + MDP(35) * t200 + t287; MDP(34) * t233 + MDP(35) * t234 + (MDP(28) * pkin(11) + t293) * t273; MDP(35) * t246 + t290 - t314; MDP(28); t281; MDP(34) * t171 - MDP(35) * t172 + t316 - t318 + t320; t284 - t315; t283; 0; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

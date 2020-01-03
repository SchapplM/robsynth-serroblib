% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:57
% EndTime: 2019-12-31 17:31:04
% DurationCPUTime: 3.59s
% Computational Cost: add. (1932->320), mult. (5367->473), div. (0->0), fcn. (4025->8), ass. (0->136)
t273 = sin(qJ(3));
t276 = cos(qJ(3));
t271 = cos(pkin(4));
t332 = qJD(1) * t271;
t304 = qJD(2) + t332;
t274 = sin(qJ(2));
t270 = sin(pkin(4));
t333 = qJD(1) * t270;
t317 = t274 * t333;
t359 = -t273 * t317 + t276 * t304;
t277 = cos(qJ(2));
t331 = qJD(1) * t277;
t316 = t270 * t331;
t261 = -qJD(3) + t316;
t225 = qJD(4) - t359;
t267 = t270 ^ 2;
t322 = qJD(1) * qJD(2);
t358 = -0.2e1 * t267 * t322;
t357 = MDP(5) * (t274 ^ 2 - t277 ^ 2);
t342 = t270 * t277;
t321 = pkin(6) * t342;
t353 = pkin(1) * t274;
t239 = t321 + (pkin(7) + t353) * t271;
t240 = (-pkin(2) * t277 - pkin(7) * t274 - pkin(1)) * t270;
t356 = t276 * t239 + t273 * t240;
t320 = pkin(1) * t332;
t244 = pkin(6) * t316 + t274 * t320;
t217 = t304 * pkin(7) + t244;
t224 = qJD(1) * t240;
t197 = t217 * t276 + t224 * t273;
t286 = (pkin(2) * t274 - pkin(7) * t277) * t270;
t243 = qJD(2) * t286;
t235 = qJD(1) * t243;
t343 = t270 * t274;
t265 = pkin(6) * t343;
t352 = pkin(1) * t277;
t245 = (t271 * t352 - t265) * qJD(2);
t236 = qJD(1) * t245;
t279 = -qJD(3) * t197 + t276 * t235 - t273 * t236;
t311 = t270 * t322;
t300 = t274 * t311;
t183 = -pkin(3) * t300 - t279;
t230 = t273 * t304 + t276 * t317;
t355 = (pkin(3) * t230 + t225 * pkin(8)) * t225 + t183;
t354 = -qJD(3) * t356 + t243 * t276 - t245 * t273;
t278 = qJD(1) ^ 2;
t299 = t277 * t311;
t209 = t359 * qJD(3) + t276 * t299;
t272 = sin(qJ(4));
t275 = cos(qJ(4));
t324 = qJD(4) * t275;
t318 = t275 * t209 - t261 * t324 + t272 * t300;
t325 = qJD(4) * t272;
t188 = -t230 * t325 + t318;
t351 = t188 * t272;
t346 = t230 * t272;
t206 = t275 * t261 + t346;
t350 = t206 * t225;
t208 = t230 * t275 - t261 * t272;
t349 = t208 * t225;
t348 = t359 * t261;
t347 = t230 * t261;
t345 = t261 * t276;
t344 = t267 * t278;
t287 = t273 * t299;
t210 = qJD(3) * t230 + t287;
t341 = t272 * t210;
t340 = t273 * t261;
t339 = t275 * t210;
t338 = t276 * t277;
t337 = t244 + t261 * (pkin(3) * t273 - pkin(8) * t276);
t241 = -pkin(6) * t317 + t277 * t320;
t242 = qJD(1) * t286;
t335 = t276 * t241 + t273 * t242;
t330 = qJD(2) * t276;
t329 = qJD(3) * t272;
t328 = qJD(3) * t273;
t327 = qJD(3) * t275;
t326 = qJD(3) * t276;
t323 = t274 * qJD(2);
t319 = t272 * t342;
t314 = t270 * t323;
t313 = qJD(2) * t342;
t283 = -t217 * t328 + t224 * t326 + t273 * t235 + t276 * t236;
t182 = pkin(8) * t300 + t283;
t246 = (t271 * t353 + t321) * qJD(2);
t237 = qJD(1) * t246;
t187 = pkin(3) * t210 - pkin(8) * t209 + t237;
t309 = -t182 * t272 + t275 * t187;
t308 = t209 * t272 - t275 * t300;
t306 = t275 * t225;
t260 = -pkin(3) * t276 - pkin(8) * t273 - pkin(2);
t305 = pkin(8) * t317 - qJD(4) * t260 + t335;
t302 = t267 * t274 * t277 * MDP(4);
t298 = MDP(15) * t317;
t297 = pkin(1) * t358;
t220 = (t272 * t274 + t275 * t338) * t333;
t295 = t275 * t326 - t220;
t294 = t182 * t275 + t187 * t272;
t216 = -t304 * pkin(2) - t241;
t193 = -pkin(3) * t359 - t230 * pkin(8) + t216;
t195 = -pkin(8) * t261 + t197;
t180 = t193 * t275 - t195 * t272;
t181 = t193 * t272 + t195 * t275;
t238 = t265 + (-pkin(2) - t352) * t271;
t247 = -t271 * t276 + t273 * t343;
t248 = t271 * t273 + t276 * t343;
t198 = pkin(3) * t247 - pkin(8) * t248 + t238;
t200 = -pkin(8) * t342 + t356;
t293 = t198 * t275 - t200 * t272;
t292 = t198 * t272 + t200 * t275;
t196 = -t217 * t273 + t224 * t276;
t290 = -t239 * t273 + t240 * t276;
t289 = -t241 * t273 + t242 * t276;
t213 = t248 * t272 + t275 * t342;
t285 = -t225 * t324 - t341;
t284 = -t225 * t325 + t339;
t282 = -t239 * t328 + t240 * t326 + t273 * t243 + t276 * t245;
t281 = pkin(1) * (-t271 * t322 + t344);
t194 = pkin(3) * t261 - t196;
t280 = -pkin(8) * t210 + (t194 + t196) * t225;
t219 = t272 * t276 * t316 - t275 * t317;
t214 = t248 * t275 - t319;
t212 = -qJD(3) * t247 + t276 * t313;
t211 = qJD(3) * t248 + t273 * t313;
t201 = -pkin(3) * t317 - t289;
t199 = pkin(3) * t342 - t290;
t192 = -qJD(4) * t213 + t212 * t275 + t272 * t314;
t191 = -qJD(4) * t319 + t212 * t272 + t248 * t324 - t275 * t314;
t190 = pkin(3) * t211 - pkin(8) * t212 + t246;
t189 = qJD(4) * t208 + t308;
t185 = -pkin(3) * t314 - t354;
t184 = pkin(8) * t314 + t282;
t179 = -qJD(4) * t181 + t309;
t178 = qJD(4) * t180 + t294;
t1 = [0.2e1 * t302 * t322 + t357 * t358 + (-t237 * t271 - t246 * t304 + t274 * t297) * MDP(9) + (-t236 * t271 - t245 * t304 + t277 * t297) * MDP(10) + (t209 * t248 + t212 * t230) * MDP(11) + (-t209 * t247 - t210 * t248 - t211 * t230 + t212 * t359) * MDP(12) + (-t212 * t261 + (-t209 * t277 + (qJD(1) * t248 + t230) * t323) * t270) * MDP(13) + (t211 * t261 + (t210 * t277 + (-qJD(1) * t247 + t359) * t323) * t270) * MDP(14) + (-t261 * t270 - t267 * t331) * MDP(15) * t323 + (-t354 * t261 - t246 * t359 + t238 * t210 + t237 * t247 + t216 * t211 + (-t279 * t277 + (qJD(1) * t290 + t196) * t323) * t270) * MDP(16) + (t282 * t261 + t246 * t230 + t238 * t209 + t237 * t248 + t216 * t212 + (t283 * t277 + (-qJD(1) * t356 - t197) * t323) * t270) * MDP(17) + (t188 * t214 + t192 * t208) * MDP(18) + (-t188 * t213 - t189 * t214 - t191 * t208 - t192 * t206) * MDP(19) + (t188 * t247 + t192 * t225 + t208 * t211 + t210 * t214) * MDP(20) + (-t189 * t247 - t191 * t225 - t206 * t211 - t210 * t213) * MDP(21) + (t210 * t247 + t211 * t225) * MDP(22) + ((-qJD(4) * t292 - t184 * t272 + t190 * t275) * t225 + t293 * t210 + t179 * t247 + t180 * t211 + t185 * t206 + t199 * t189 + t183 * t213 + t194 * t191) * MDP(23) + (-(qJD(4) * t293 + t184 * t275 + t190 * t272) * t225 - t292 * t210 - t178 * t247 - t181 * t211 + t185 * t208 + t199 * t188 + t183 * t214 + t194 * t192) * MDP(24) + (MDP(6) * t313 - MDP(7) * t314) * (qJD(2) + 0.2e1 * t332); t344 * t357 + (-pkin(6) * t299 + t244 * t304 + t274 * t281) * MDP(9) + (pkin(6) * t300 + t241 * t304 + t277 * t281) * MDP(10) + (t209 * t273 - t230 * t345) * MDP(11) + ((t209 - t348) * t276 + (-t210 + t347) * t273) * MDP(12) + (-t261 * t326 + (t261 * t338 + (qJD(2) * t273 - t230) * t274) * t333) * MDP(13) + (t261 * t328 + (-t277 * t340 + (-t359 + t330) * t274) * t333) * MDP(14) + t261 * t298 + (-pkin(2) * t210 - t237 * t276 + t289 * t261 + t244 * t359 + (pkin(7) * t345 + t216 * t273) * qJD(3) + (-t196 * t274 + (-pkin(7) * t323 - t216 * t277) * t273) * t333) * MDP(16) + (-pkin(2) * t209 + t237 * t273 - t335 * t261 - t244 * t230 + (-pkin(7) * t340 + t216 * t276) * qJD(3) + (-t216 * t338 + (-pkin(7) * t330 + t197) * t274) * t333) * MDP(17) + (t188 * t273 * t275 + (-t273 * t325 + t295) * t208) * MDP(18) + (t206 * t220 + t208 * t219 + (-t206 * t275 - t208 * t272) * t326 + (-t351 - t189 * t275 + (t206 * t272 - t208 * t275) * qJD(4)) * t273) * MDP(19) + (-t188 * t276 + t295 * t225 + (-t261 * t208 + t284) * t273) * MDP(20) + (t189 * t276 + (-t272 * t326 + t219) * t225 + (t261 * t206 + t285) * t273) * MDP(21) + (-t210 * t276 - t225 * t340) * MDP(22) + (t260 * t339 - t194 * t219 - t201 * t206 + (t305 * t272 - t337 * t275) * t225 + (t194 * t329 - t179 + (qJD(3) * t206 + t285) * pkin(7)) * t276 + (t194 * t324 + t183 * t272 - t261 * t180 + (t225 * t329 + t189) * pkin(7)) * t273) * MDP(23) + (-t260 * t341 - t194 * t220 - t201 * t208 + (t337 * t272 + t305 * t275) * t225 + (t194 * t327 + t178 + (qJD(3) * t208 - t284) * pkin(7)) * t276 + (-t194 * t325 + t183 * t275 + t261 * t181 + (t225 * t327 + t188) * pkin(7)) * t273) * MDP(24) + (-t302 + (-t277 * MDP(6) + t274 * MDP(7)) * t270 * t271) * t278; -t359 ^ 2 * MDP(12) + (t209 + t348) * MDP(13) + (-t287 - t347) * MDP(14) + qJD(2) * t298 + (-t197 * t261 + t279) * MDP(16) + (-t196 * t261 - t216 * t359 - t283) * MDP(17) + (t208 * t306 + t351) * MDP(18) + ((t188 - t350) * t275 + (-t189 - t349) * t272) * MDP(19) + (t225 * t306 + t341) * MDP(20) + (-t225 ^ 2 * t272 + t339) * MDP(21) + (-pkin(3) * t189 - t197 * t206 + t280 * t272 - t355 * t275) * MDP(23) + (-pkin(3) * t188 - t197 * t208 + t355 * t272 + t280 * t275) * MDP(24) + (-MDP(11) * t359 + t230 * MDP(12) - qJD(3) * MDP(14) - t216 * MDP(16) - t208 * MDP(20) + t206 * MDP(21) - t225 * MDP(22) - t180 * MDP(23) + t181 * MDP(24)) * t230; t208 * t206 * MDP(18) + (-t206 ^ 2 + t208 ^ 2) * MDP(19) + (t318 + t350) * MDP(20) + (-t308 + t349) * MDP(21) + t210 * MDP(22) + (t181 * t225 - t194 * t208 + t309) * MDP(23) + (t180 * t225 + t194 * t206 - t294) * MDP(24) + (-MDP(20) * t346 - MDP(21) * t208 - MDP(23) * t181 - MDP(24) * t180) * qJD(4);];
tauc = t1;

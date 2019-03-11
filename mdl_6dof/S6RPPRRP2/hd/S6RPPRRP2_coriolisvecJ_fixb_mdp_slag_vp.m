% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:01:22
% EndTime: 2019-03-09 02:01:28
% DurationCPUTime: 3.02s
% Computational Cost: add. (3578->335), mult. (8451->433), div. (0->0), fcn. (6065->8), ass. (0->129)
t301 = sin(pkin(9)) * pkin(1) + qJ(3);
t297 = t301 * qJD(1);
t311 = cos(pkin(10));
t305 = t311 * qJD(2);
t309 = sin(pkin(10));
t269 = t305 + (-pkin(7) * qJD(1) - t297) * t309;
t280 = t309 * qJD(2) + t311 * t297;
t353 = qJD(1) * t311;
t270 = pkin(7) * t353 + t280;
t314 = sin(qJ(4));
t316 = cos(qJ(4));
t234 = t314 * t269 + t316 * t270;
t398 = qJD(4) * t234;
t294 = t309 * t316 + t311 * t314;
t397 = t294 * qJD(3);
t315 = cos(qJ(5));
t350 = qJD(5) * t315;
t313 = sin(qJ(5));
t352 = qJD(4) * t313;
t388 = t294 * qJD(1);
t366 = t309 * t314;
t293 = -t316 * t311 + t366;
t389 = t293 * qJD(1);
t243 = t388 * t350 + (qJD(5) - t389) * t352;
t338 = qJD(1) * (t309 ^ 2 + t311 ^ 2);
t396 = MDP(7) * t338;
t286 = -qJD(1) * t366 + t316 * t353;
t283 = qJD(5) - t286;
t320 = t293 * qJD(3);
t390 = t269 * t316 - t314 * t270;
t222 = -qJD(1) * t320 + qJD(4) * t390;
t231 = qJD(4) * pkin(8) + t234;
t295 = -cos(pkin(9)) * pkin(1) - pkin(3) * t311 - pkin(2);
t284 = t295 * qJD(1) + qJD(3);
t241 = -pkin(4) * t286 - pkin(8) * t388 + t284;
t289 = t294 * qJD(4);
t281 = qJD(1) * t289;
t288 = t293 * qJD(4);
t318 = qJD(1) * t288;
t249 = t281 * pkin(4) + pkin(8) * t318;
t351 = qJD(5) * t313;
t324 = t315 * t222 - t231 * t351 + t241 * t350 + t313 * t249;
t381 = qJ(6) * t281;
t203 = qJD(6) * t283 + t324 + t381;
t211 = -t231 * t313 + t241 * t315;
t348 = qJD(6) - t211;
t209 = -pkin(5) * t283 + t348;
t395 = t283 * t209 + t203;
t336 = t313 * t222 + t231 * t350 + t241 * t351 - t315 * t249;
t384 = pkin(5) * t281;
t204 = t336 - t384;
t212 = t231 * t315 + t241 * t313;
t210 = qJ(6) * t283 + t212;
t379 = t210 * t283;
t394 = -t204 + t379;
t349 = t315 * qJD(4);
t271 = t313 * t388 - t349;
t391 = t293 * t243 + t289 * t271;
t242 = -qJD(5) * t349 + t315 * t318 + t351 * t388;
t273 = t315 * t388 + t352;
t359 = -t293 * t242 + t273 * t289;
t253 = pkin(4) * t293 - pkin(8) * t294 + t295;
t382 = pkin(7) + t301;
t290 = t382 * t309;
t291 = t382 * t311;
t257 = -t290 * t314 + t291 * t316;
t358 = t313 * t253 + t315 * t257;
t256 = t290 * t316 + t314 * t291;
t387 = MDP(21) + MDP(23);
t230 = -qJD(4) * pkin(4) - t390;
t215 = pkin(5) * t271 - qJ(6) * t273 + t230;
t383 = pkin(8) * t281;
t386 = t283 * t215 - t383;
t385 = t273 ^ 2;
t223 = qJD(1) * t397 + t398;
t206 = pkin(5) * t243 + qJ(6) * t242 - qJD(6) * t273 + t223;
t380 = t206 * t313;
t378 = t212 * t283;
t377 = t242 * t313;
t375 = t271 * t273;
t374 = t271 * t286;
t373 = t271 * t313;
t337 = t273 * t283;
t372 = t273 * t315;
t371 = t283 * t313;
t370 = t288 * t313;
t369 = t288 * t315;
t367 = t294 * t315;
t276 = t313 * t281;
t278 = t315 * t281;
t334 = pkin(5) * t313 - qJ(6) * t315;
t363 = qJD(6) * t313 - t283 * t334 + t234;
t362 = -t243 * t367 + t271 * t369;
t260 = pkin(4) * t388 - pkin(8) * t286;
t361 = t313 * t260 + t315 * t390;
t360 = -t313 * t243 - t271 * t350;
t357 = t286 * t371 + t278;
t259 = t283 * t369;
t356 = t294 * t278 - t259;
t355 = t283 * t350 + t276;
t345 = qJD(4) * MDP(11);
t343 = MDP(22) - MDP(25);
t342 = pkin(8) * t351;
t341 = t273 * t370;
t340 = t294 * t351;
t335 = pkin(5) * t315 + qJ(6) * t313;
t333 = t203 * t315 + t204 * t313;
t332 = t209 * t315 - t210 * t313;
t330 = (-t297 * t309 + t305) * t309 - t280 * t311;
t328 = t294 * t350 - t370;
t327 = t340 + t369;
t326 = t283 * t230 - t383;
t325 = t215 * t273 + t336;
t237 = -t256 * qJD(4) - t320;
t261 = pkin(4) * t289 + pkin(8) * t288;
t323 = t315 * t237 + t253 * t350 - t257 * t351 + t313 * t261;
t238 = t257 * qJD(4) + t397;
t296 = -pkin(4) - t335;
t239 = pkin(5) * t273 + qJ(6) * t271;
t224 = t334 * t294 + t256;
t219 = t271 * t283 - t242;
t217 = -pkin(5) * t293 - t253 * t315 + t257 * t313;
t216 = qJ(6) * t293 + t358;
t214 = -pkin(5) * t388 - t260 * t315 + t313 * t390;
t213 = qJ(6) * t388 + t361;
t208 = -t334 * t288 + (t335 * qJD(5) - qJD(6) * t315) * t294 + t238;
t207 = -pkin(5) * t289 + t358 * qJD(5) + t237 * t313 - t261 * t315;
t205 = qJ(6) * t289 + qJD(6) * t293 + t323;
t1 = [(-t288 * t388 - t294 * t318) * MDP(9) + (-t294 * t281 - t288 * t286 - t289 * t388 + t293 * t318) * MDP(10) - t288 * t345 - t289 * qJD(4) * MDP(12) + (-qJD(4) * t238 + t281 * t295 + t284 * t289) * MDP(14) + (-t284 * t288 + (-t295 * t389 - t237) * qJD(4)) * MDP(15) + (-t242 * t367 - t327 * t273) * MDP(16) + (t341 + (t377 + (-t372 + t373) * qJD(5)) * t294 + t362) * MDP(17) + (-t283 * t340 + t356 + t359) * MDP(18) + (-t294 * t276 - t328 * t283 - t391) * MDP(19) + (t281 * t293 + t283 * t289) * MDP(20) + (-t336 * t293 + t211 * t289 + t238 * t271 + t256 * t243 + ((-qJD(5) * t257 + t261) * t283 + t253 * t281 + t230 * qJD(5) * t294) * t315 + ((-qJD(5) * t253 - t237) * t283 - t257 * t281 + t223 * t294 - t230 * t288) * t313) * MDP(21) + (-t212 * t289 + t223 * t367 - t327 * t230 + t238 * t273 - t256 * t242 - t358 * t281 - t323 * t283 - t324 * t293) * MDP(22) + (-t204 * t293 - t207 * t283 + t208 * t271 - t209 * t289 + t328 * t215 - t217 * t281 + t224 * t243 + t294 * t380) * MDP(23) + (-t205 * t271 + t207 * t273 - t216 * t243 - t217 * t242 - t332 * t288 + (-t203 * t313 + t204 * t315 + (-t209 * t313 - t210 * t315) * qJD(5)) * t294) * MDP(24) + (t203 * t293 + t205 * t283 - t206 * t367 - t208 * t273 + t210 * t289 + t327 * t215 + t216 * t281 + t224 * t242) * MDP(25) + (t203 * t216 + t204 * t217 + t205 * t210 + t206 * t224 + t207 * t209 + t208 * t215) * MDP(26) + (0.2e1 * t396 + (t301 * t338 - t330) * MDP(8)) * qJD(3); (t259 + t359) * MDP(22) + (-t341 + t362) * MDP(24) + (t356 - t359) * MDP(25) + (t206 * t293 - t209 * t370 - t210 * t369 + t215 * t289) * MDP(26) + (-MDP(14) * t289 + MDP(15) * t288) * qJD(4) + (-MDP(24) * t377 + t333 * MDP(26) + (-t315 * MDP(22) - t313 * t387) * t281 + ((t372 + t373) * MDP(24) + t332 * MDP(26) + (t343 * t313 - t315 * t387) * t283) * qJD(5)) * t294 + t387 * (t283 * t370 + t391); t357 * MDP(21) + t360 * MDP(24) + t355 * MDP(25) + (-t215 * MDP(26) - t271 * t387 - t343 * t273) * t388 + (t281 * MDP(23) + (t242 + t374) * MDP(24) + t394 * MDP(26) + (-t283 * MDP(22) - t286 * MDP(25)) * t283) * t315 + (-t281 * MDP(22) + t395 * MDP(26) + MDP(24) * t337 + (-qJD(5) * MDP(21) - t283 * MDP(23)) * t283) * t313 + (MDP(14) * t388 + t286 * MDP(15)) * qJD(4) + (t330 * MDP(8) + (t294 * MDP(14) - t293 * MDP(15)) * qJD(4) - t396) * qJD(1); -t286 ^ 2 * MDP(10) + (-t286 - t389) * t345 + (-t223 + t398) * MDP(14) + (qJD(3) * t389 - t284 * t286) * MDP(15) + (t315 * t337 - t377) * MDP(16) + ((-t242 + t374) * t315 - t273 * t371 + t360) * MDP(17) + (-t283 * t286 * t315 + t355) * MDP(18) + (-t283 * t351 + t357) * MDP(19) + (-pkin(4) * t243 - t234 * t271 + (-t223 + (-pkin(8) * qJD(5) - t260) * t283) * t315 + (t283 * t390 + t326) * t313) * MDP(21) + (pkin(4) * t242 + t223 * t313 - t234 * t273 + (t342 + t361) * t283 + t326 * t315) * MDP(22) + (-t206 * t315 + t243 * t296 + (-pkin(8) * t350 + t214) * t283 - t363 * t271 + t386 * t313) * MDP(23) + (t213 * t271 - t214 * t273 + ((qJD(5) * t273 - t243) * pkin(8) + t395) * t315 + ((qJD(5) * t271 - t242) * pkin(8) - t394) * t313) * MDP(24) + (-t380 + t242 * t296 + (-t213 - t342) * t283 + t363 * t273 - t386 * t315) * MDP(25) + (t206 * t296 - t209 * t214 - t210 * t213 - t363 * t215 + (t332 * qJD(5) + t333) * pkin(8)) * MDP(26) + (MDP(10) * t388 - t284 * MDP(14) - t273 * MDP(18) + t271 * MDP(19) - t283 * MDP(20) - t211 * MDP(21) + t212 * MDP(22) + t209 * MDP(23) - t210 * MDP(25) - t286 * MDP(9)) * t388; MDP(16) * t375 + (-t271 ^ 2 + t385) * MDP(17) + t219 * MDP(18) + (-t243 + t337) * MDP(19) + t281 * MDP(20) + (-t230 * t273 - t336 + t378) * MDP(21) + (t211 * t283 + t230 * t271 - t324) * MDP(22) + (-t239 * t271 - t325 + t378 + 0.2e1 * t384) * MDP(23) + (pkin(5) * t242 - qJ(6) * t243 + (t210 - t212) * t273 + (t209 - t348) * t271) * MDP(24) + (0.2e1 * t381 - t215 * t271 + t239 * t273 + (0.2e1 * qJD(6) - t211) * t283 + t324) * MDP(25) + (-pkin(5) * t204 + qJ(6) * t203 - t209 * t212 + t348 * t210 - t215 * t239) * MDP(26); (-qJD(4) * t388 + t375) * MDP(23) + t219 * MDP(24) + (-t283 ^ 2 - t385) * MDP(25) + (t325 - t379 - t384) * MDP(26);];
tauc  = t1;

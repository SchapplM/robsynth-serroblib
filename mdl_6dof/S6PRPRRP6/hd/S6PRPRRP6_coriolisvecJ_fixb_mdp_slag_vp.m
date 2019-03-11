% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:21:01
% EndTime: 2019-03-08 20:21:08
% DurationCPUTime: 2.68s
% Computational Cost: add. (2212->352), mult. (4819->477), div. (0->0), fcn. (3148->8), ass. (0->157)
t286 = sin(qJ(4));
t355 = qJD(2) * t286;
t278 = qJD(5) + t355;
t289 = cos(qJ(4));
t282 = t289 ^ 2;
t397 = MDP(9) * (t286 ^ 2 - t282);
t291 = -pkin(2) - pkin(8);
t290 = cos(qJ(2));
t283 = sin(pkin(6));
t359 = qJD(1) * t283;
t332 = t290 * t359;
t312 = qJD(3) - t332;
t257 = qJD(2) * t291 + t312;
t284 = cos(pkin(6));
t287 = sin(qJ(2));
t356 = qJD(2) * t283;
t331 = t287 * t356;
t351 = qJD(4) * t289;
t358 = qJD(1) * t286;
t222 = t257 * t351 + (-qJD(4) * t284 + t331) * t358;
t373 = t284 * t289;
t277 = qJD(1) * t373;
t239 = t286 * t257 + t277;
t230 = qJD(4) * pkin(9) + t239;
t317 = pkin(4) * t289 + pkin(9) * t286;
t261 = qJD(4) * t317 + qJD(3);
t240 = (t261 + t332) * qJD(2);
t270 = pkin(4) * t286 - pkin(9) * t289 + qJ(3);
t333 = t287 * t359;
t248 = qJD(2) * t270 + t333;
t285 = sin(qJ(5));
t288 = cos(qJ(5));
t348 = qJD(5) * t288;
t349 = qJD(5) * t285;
t322 = -t285 * t222 - t230 * t348 + t288 * t240 - t248 * t349;
t342 = qJD(2) * qJD(4);
t326 = t289 * t342;
t323 = pkin(5) * t326;
t203 = -t322 - t323;
t211 = t230 * t288 + t248 * t285;
t208 = qJ(6) * t278 + t211;
t396 = -t208 * t278 + t203;
t371 = t285 * t287;
t395 = -(t286 * t371 - t288 * t290) * t359 - t261 * t288;
t344 = t288 * qJD(4);
t368 = t287 * t288;
t394 = (t285 * t290 + t286 * t368) * t359 - t289 * t291 * t344 - t285 * t261 - t270 * t348;
t369 = t286 * t291;
t362 = t285 * t270 + t288 * t369;
t238 = t257 * t289 - t284 * t358;
t340 = MDP(20) + MDP(22);
t339 = MDP(21) - MDP(24);
t347 = qJD(5) * t289;
t328 = t285 * t347;
t300 = t286 * t344 + t328;
t341 = qJD(4) * qJD(5);
t243 = qJD(2) * t300 - t288 * t341;
t393 = qJD(2) * pkin(2);
t320 = t289 * t331;
t352 = qJD(4) * t286;
t363 = -qJD(4) * t277 - t257 * t352;
t223 = -qJD(1) * t320 - t363;
t372 = t285 * t286;
t274 = t342 * t372;
t353 = qJD(4) * t285;
t354 = qJD(2) * t289;
t266 = t288 * t354 + t353;
t350 = qJD(5) * t266;
t244 = -t274 + t350;
t204 = pkin(5) * t244 + qJ(6) * t243 - qJD(6) * t266 + t223;
t392 = t204 * t285;
t391 = t204 * t288;
t389 = t223 * t285;
t388 = t223 * t288;
t229 = -qJD(4) * pkin(4) - t238;
t387 = t229 * t285;
t386 = t229 * t288;
t385 = t243 * t285;
t384 = t244 * t288;
t264 = t285 * t354 - t344;
t381 = t264 * t278;
t380 = t264 * t285;
t379 = t264 * t288;
t378 = t266 * t285;
t377 = t266 * t288;
t376 = t270 * t288;
t375 = t278 * t288;
t374 = t283 * t290;
t370 = t285 * t291;
t367 = qJ(6) * t351 + (-t291 * t349 + qJD(6)) * t286 - t394;
t324 = -pkin(5) + t370;
t366 = qJD(5) * t362 + t324 * t351 + t395;
t313 = pkin(5) * t285 - qJ(6) * t288;
t365 = qJD(6) * t285 - t278 * t313 + t239;
t268 = t317 * qJD(2);
t364 = t288 * t238 + t285 * t268;
t292 = qJD(4) ^ 2;
t293 = qJD(2) ^ 2;
t360 = -t292 - t293;
t357 = qJD(2) * qJ(3);
t346 = qJD(6) * t278;
t345 = t229 * qJD(5);
t210 = -t230 * t285 + t248 * t288;
t343 = qJD(6) - t210;
t338 = pkin(9) * t278 * t285;
t337 = pkin(9) * t375;
t336 = pkin(9) * t351;
t335 = -t288 * t222 - t285 * t240 - t248 * t348;
t330 = t290 * t356;
t329 = t278 * t348;
t327 = t288 * t347;
t325 = MDP(19) * t351;
t321 = t289 * t333;
t319 = t266 * t333;
t318 = qJ(6) * t326;
t269 = t333 + t357;
t315 = -t269 + t333;
t314 = pkin(5) * t288 + qJ(6) * t285;
t301 = t230 * t349 + t335;
t202 = -t301 + t318 + t346;
t311 = t202 * t288 + t203 * t285;
t207 = -pkin(5) * t278 + t343;
t310 = t207 * t288 - t208 * t285;
t309 = t207 * t285 + t208 * t288;
t308 = -t238 * t285 + t268 * t288;
t307 = qJD(2) * t282 - t278 * t286;
t306 = -t291 + t313;
t212 = pkin(5) * t264 - qJ(6) * t266 + t229;
t305 = -t212 * t286 + t336;
t304 = t229 * t286 - t336;
t256 = -t286 * t374 + t373;
t303 = -t256 * t285 + t283 * t368;
t235 = t256 * t288 + t283 * t371;
t255 = t284 * t286 + t289 * t374;
t302 = t211 * t278 + t322;
t299 = t315 - t357;
t297 = -t285 * t340 - t288 * t339;
t296 = t210 * t278 + t301;
t262 = (qJD(3) + t332) * qJD(2);
t295 = qJD(2) * t312 - t291 * t292 + t262;
t294 = (t377 + t380) * MDP(23) + t310 * MDP(25) + (t285 * t339 - t288 * t340) * t278;
t271 = -pkin(4) - t314;
t263 = t312 - t393;
t249 = t306 * t289;
t247 = t264 * t321;
t242 = t286 * t324 - t376;
t241 = qJ(6) * t286 + t362;
t233 = qJD(4) * t256 - t320;
t232 = -qJD(4) * t255 + t286 * t331;
t231 = pkin(5) * t266 + qJ(6) * t264;
t219 = -t243 + t381;
t216 = (qJD(5) * t314 - qJD(6) * t288) * t289 - t306 * t352;
t215 = -pkin(5) * t354 - t308;
t214 = qJ(6) * t354 + t364;
t206 = qJD(5) * t303 + t232 * t288 + t285 * t330;
t205 = qJD(5) * t235 + t232 * t285 - t288 * t330;
t1 = [(t205 * t266 - t206 * t264 - t235 * t244 + t243 * t303) * MDP(23) + (t202 * t235 - t203 * t303 + t204 * t255 + t205 * t207 + t206 * t208 + t212 * t233) * MDP(25) + (-t233 * MDP(13) - t232 * MDP(14) + (-t235 * t339 + t303 * t340) * t354) * qJD(4) + ((qJD(2) * t269 * MDP(7) + (t286 * MDP(13) + MDP(14) * t289 - MDP(4) + MDP(6)) * t293) * t290 + (t262 * MDP(7) + (-MDP(3) + MDP(5)) * t293 + ((MDP(13) * t289 - MDP(14) * t286) * qJD(4) + (t263 - t332) * MDP(7)) * qJD(2)) * t287) * t283 + t340 * (-t205 * t278 + t233 * t264 + t255 * t244) - t339 * (t206 * t278 - t233 * t266 + t243 * t255); 0.2e1 * qJD(2) * qJD(3) * MDP(6) + (qJ(3) * t262 + qJD(3) * t269 + (-t269 * t290 + (-t263 - t393) * t287) * t359) * MDP(7) + 0.2e1 * t342 * t397 - t292 * t289 * MDP(11) - t299 * t351 * MDP(13) + (t295 * t289 + t299 * t352) * MDP(14) + (-t243 * t288 * t289 - t300 * t266) * MDP(15) + ((t378 + t379) * t352 + (t385 - t384 + (-t377 + t380) * qJD(5)) * t289) * MDP(16) + (-t278 * t328 + (t266 * t289 + t307 * t288) * qJD(4)) * MDP(17) + (-t278 * t327 + (-t264 * t289 - t307 * t285) * qJD(4)) * MDP(18) + (t278 + t355) * t325 + (t247 + (-t270 * t349 - t395) * t278 + (t288 * t345 + t389 - t291 * t244 + (-t278 * t370 + (-t285 * t369 + t376) * qJD(2) + t210) * qJD(4)) * t289) * MDP(20) + (t394 * t278 + (t319 - t285 * t345 + t388 + t291 * t243 + (-t362 * qJD(2) - t211) * qJD(4)) * t289) * MDP(21) + (t216 * t264 + t244 * t249 + t247 - t366 * t278 + (t212 * t348 + t392 + (-qJD(2) * t242 - t207) * qJD(4)) * t289) * MDP(22) + (-t241 * t244 - t242 * t243 + t366 * t266 - t367 * t264 - t310 * t352 + (-t309 * qJD(5) - t202 * t285 + t203 * t288) * t289) * MDP(23) + (-t216 * t266 + t243 * t249 + t367 * t278 + (-t319 + t212 * t349 - t391 + (qJD(2) * t241 + t208) * qJD(4)) * t289) * MDP(24) + (t202 * t241 + t203 * t242 + t204 * t249 + (t216 + t321) * t212 + t367 * t208 + t366 * t207) * MDP(25) + (-0.2e1 * MDP(8) * t326 - t292 * MDP(10) + t295 * MDP(13) - t243 * MDP(17) - t244 * MDP(18) + (-t291 * t329 + (t264 * t291 - t387) * qJD(4) + t322) * MDP(20) + ((t278 * t291 + t230) * t349 + (t266 * t291 - t386) * qJD(4) + t335) * MDP(21) + (-t212 * t353 - t203) * MDP(22) + (t212 * t344 + t202) * MDP(24)) * t286; -t293 * MDP(6) + t340 * t264 * t352 + (MDP(7) * t315 + t294) * qJD(2) + (t360 * MDP(14) - t204 * MDP(25) - t340 * t244 + t339 * t243 + ((t378 - t379) * MDP(23) + t309 * MDP(25) + t297 * t278) * qJD(4)) * t289 + (t360 * MDP(13) + (-t384 - t385) * MDP(23) + t311 * MDP(25) + t294 * qJD(5) + (t212 * MDP(25) + t266 * t339 + t297 * t354) * qJD(4)) * t286; (qJD(4) * t239 + t315 * t354 + t363) * MDP(13) - t315 * t355 * MDP(14) + (t266 * t375 - t385) * MDP(15) + ((-t243 - t381) * t288 + (-t266 * t278 - t244) * t285) * MDP(16) + (t329 + (t286 * t375 + (-t266 + t353) * t289) * qJD(2)) * MDP(17) + (-t278 * t349 + (-t278 * t372 + (t264 + t344) * t289) * qJD(2)) * MDP(18) - t278 * MDP(19) * t354 + (-pkin(4) * t244 - t388 - t308 * t278 - t239 * t264 + (-t337 + t387) * qJD(5) + (-t210 * t289 + t285 * t304) * qJD(2)) * MDP(20) + (pkin(4) * t243 + t389 + t364 * t278 - t239 * t266 + (t338 + t386) * qJD(5) + (t211 * t289 + t288 * t304) * qJD(2)) * MDP(21) + (-t391 + t215 * t278 + t244 * t271 - t365 * t264 + (t212 * t285 - t337) * qJD(5) + (t207 * t289 - t285 * t305) * qJD(2)) * MDP(22) + (t214 * t264 - t215 * t266 + (t202 + t278 * t207 + (-t244 + t350) * pkin(9)) * t288 + ((qJD(5) * t264 - t243) * pkin(9) + t396) * t285) * MDP(23) + (-t392 - t214 * t278 + t243 * t271 + t365 * t266 + (-t212 * t288 - t338) * qJD(5) + (-t208 * t289 + t288 * t305) * qJD(2)) * MDP(24) + (t204 * t271 - t207 * t215 - t208 * t214 - t365 * t212 + (qJD(5) * t310 + t311) * pkin(9)) * MDP(25) + (t289 * t286 * MDP(8) - t397) * t293; t219 * MDP(17) + (-qJD(2) * t327 - t285 * t341 + t274) * MDP(18) + qJD(2) * t325 + t302 * MDP(20) + t296 * MDP(21) + (t302 + 0.2e1 * t323) * MDP(22) + (pkin(5) * t243 - qJ(6) * t244) * MDP(23) + (-t296 + 0.2e1 * t318 + 0.2e1 * t346) * MDP(24) + (-pkin(5) * t203 + qJ(6) * t202 - t207 * t211 + t208 * t343 - t212 * t231) * MDP(25) + (t278 * MDP(18) - t229 * MDP(20) - t212 * MDP(22) + (t208 - t211) * MDP(23) + t231 * MDP(24) + MDP(16) * t266) * t266 + (t266 * MDP(15) + t229 * MDP(21) - t231 * MDP(22) + (t207 - t343) * MDP(23) - t212 * MDP(24) - MDP(16) * t264) * t264; (t264 * t266 - t326) * MDP(22) + t219 * MDP(23) + (-t266 ^ 2 - t278 ^ 2) * MDP(24) + (t212 * t266 + t396) * MDP(25);];
tauc  = t1;

% Calculate vector of inverse dynamics joint torques for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP12_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP12_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:26:03
% EndTime: 2021-01-15 19:26:13
% DurationCPUTime: 3.43s
% Computational Cost: add. (1981->375), mult. (3759->473), div. (0->0), fcn. (2194->6), ass. (0->168)
t298 = cos(qJ(1));
t289 = g(2) * t298;
t295 = sin(qJ(1));
t406 = g(1) * t295;
t414 = -t406 + t289;
t417 = qJDD(2) + t414;
t296 = cos(qJ(4));
t284 = pkin(4) * t296 + pkin(3);
t292 = -qJ(5) - pkin(7);
t294 = sin(qJ(3));
t297 = cos(qJ(3));
t332 = -t284 * t294 - t292 * t297;
t416 = -qJ(2) + t332;
t266 = pkin(3) * t294 - pkin(7) * t297 + qJ(2);
t293 = sin(qJ(4));
t299 = -pkin(1) - pkin(6);
t382 = t294 * t296;
t372 = t293 * t266 + t299 * t382;
t281 = pkin(4) * t293 - t299;
t399 = qJDD(1) * pkin(1);
t415 = t399 - t417;
t365 = qJD(3) * t293;
t368 = qJD(1) * t297;
t261 = t296 * t368 + t365;
t413 = qJD(4) * t261;
t379 = t296 * t298;
t385 = t293 * t295;
t243 = -t294 * t385 + t379;
t381 = t295 * t296;
t383 = t293 * t298;
t245 = t294 * t383 + t381;
t384 = t293 * t297;
t412 = -g(1) * t243 - g(2) * t245 + g(3) * t384;
t357 = t296 * qJD(3);
t341 = t294 * t357;
t359 = qJD(4) * t297;
t311 = t293 * t359 + t341;
t351 = qJDD(1) * t297;
t354 = qJD(3) * qJD(4);
t347 = t293 * qJDD(3) + (t351 + t354) * t296;
t222 = qJD(1) * t311 - t347;
t356 = qJD(1) * qJD(3);
t340 = t294 * t356;
t373 = t296 * qJDD(3) + t293 * t340;
t223 = t293 * t351 - t373 + t413;
t411 = t261 ^ 2;
t339 = t297 * t356;
t352 = qJDD(1) * t294;
t256 = qJDD(4) + t339 + t352;
t410 = pkin(4) * t256;
t259 = t293 * t368 - t357;
t409 = pkin(4) * t259;
t404 = g(3) * t294;
t403 = g(3) * t297;
t301 = qJD(1) ^ 2;
t402 = qJ(2) * t301;
t401 = qJ(5) * t222;
t400 = qJ(5) * t223;
t398 = t222 * t293;
t397 = t223 * t296;
t396 = t256 * t293;
t395 = t256 * t296;
t369 = qJD(1) * t294;
t277 = qJD(4) + t369;
t394 = t259 * t277;
t393 = t259 * t293;
t392 = t259 * t296;
t391 = t261 * t277;
t390 = t261 * t293;
t389 = t261 * t296;
t273 = qJDD(1) * t299 + qJDD(2);
t388 = t273 * t297;
t275 = qJD(1) * t299 + qJD(2);
t387 = t275 * t297;
t265 = t294 * t275;
t380 = t296 * t297;
t242 = t266 * qJD(1);
t249 = qJD(3) * pkin(7) + t265;
t220 = t296 * t242 - t249 * t293;
t216 = -qJ(5) * t261 + t220;
t215 = pkin(4) * t277 + t216;
t378 = -t216 + t215;
t336 = qJD(4) * t292;
t346 = t293 * t369;
t358 = qJD(5) * t296;
t328 = pkin(3) * t297 + pkin(7) * t294;
t264 = t328 * qJD(1);
t375 = t293 * t264 + t275 * t380;
t377 = -qJ(5) * t346 + t293 * t336 + t358 - t375;
t248 = t296 * t264;
t376 = -qJD(5) * t293 + t296 * t336 + t275 * t384 - t248 - (pkin(4) * t297 + qJ(5) * t382) * qJD(1);
t349 = t297 * t289;
t374 = g(3) * t382 + t296 * t349;
t291 = t297 ^ 2;
t371 = t294 ^ 2 - t291;
t300 = qJD(3) ^ 2;
t370 = -t300 - t301;
t367 = qJD(3) * t259;
t366 = qJD(3) * t261;
t364 = qJD(3) * t294;
t363 = qJD(3) * t297;
t362 = qJD(3) * t299;
t361 = qJD(4) * t293;
t360 = qJD(4) * t296;
t355 = qJD(1) * qJD(4);
t353 = qJDD(1) * qJ(2);
t350 = t297 * t406;
t257 = qJD(3) * t328 + qJD(2);
t345 = t297 * t362;
t348 = t293 * t257 + t266 * t360 + t296 * t345;
t250 = -qJD(3) * pkin(3) - t387;
t335 = -qJD(5) - t409;
t229 = t250 - t335;
t344 = t229 * t360;
t343 = t277 * t361;
t342 = t296 * t359;
t338 = -t293 * t299 + pkin(4);
t334 = -qJDD(3) * pkin(3) + t275 * t364;
t331 = qJD(4) * t294 + qJD(1);
t227 = qJD(1) * t257 + qJDD(1) * t266;
t232 = qJDD(3) * pkin(7) + t273 * t294 + t275 * t363;
t330 = t293 * t227 + t296 * t232 + t242 * t360 - t249 * t361;
t231 = t334 - t388;
t329 = -pkin(7) * qJD(4) * t277 - t231;
t327 = g(1) * t245 - g(2) * t243;
t244 = t294 * t381 + t383;
t246 = t294 * t379 - t385;
t326 = -g(1) * t246 - g(2) * t244;
t325 = g(1) * t298 + g(2) * t295;
t221 = t242 * t293 + t249 * t296;
t225 = t296 * t227;
t305 = -qJD(4) * t221 - t232 * t293 + t225;
t208 = -qJD(5) * t261 + t305 + t401 + t410;
t209 = -qJD(5) * t259 + t330 - t400;
t323 = -t208 * t293 + t209 * t296;
t217 = -qJ(5) * t259 + t221;
t322 = t215 * t296 + t217 * t293;
t321 = t215 * t293 - t217 * t296;
t319 = -t349 - t404;
t318 = -t273 - t414;
t315 = t277 * t360 + t396;
t314 = -t343 + t395;
t313 = 0.2e1 * qJ(2) * t356 + qJDD(3) * t299;
t312 = pkin(4) * t223 + qJDD(5) + t334;
t310 = 0.2e1 * qJD(1) * qJD(2) - t325;
t309 = t318 + t402;
t308 = -pkin(7) * t256 + t250 * t277;
t307 = g(1) * t244 - g(2) * t246 + g(3) * t380 - t330;
t306 = t310 + 0.2e1 * t353;
t304 = -t299 * t300 + t306;
t303 = (t389 + t393) * MDP(23) - t322 * MDP(24);
t302 = t305 + t412;
t287 = qJDD(3) * t297;
t271 = t293 * t350;
t269 = t292 * t296;
t268 = t292 * t293;
t258 = t281 * t297;
t255 = t259 ^ 2;
t254 = t296 * t266;
t241 = t296 * t257;
t234 = -pkin(4) * t346 + t265;
t233 = t294 * t362 + (-t293 * t364 + t342) * pkin(4);
t230 = -qJ(5) * t384 + t372;
t226 = -qJ(5) * t380 + t294 * t338 + t254;
t214 = t312 - t388;
t213 = -qJ(5) * t342 + (-qJD(5) * t297 + (qJ(5) * qJD(3) - qJD(4) * t299) * t294) * t293 + t348;
t212 = qJ(5) * t341 + t241 - t372 * qJD(4) + (qJ(5) * t361 + qJD(3) * t338 - t358) * t297;
t1 = [qJDD(1) * MDP(1) - t414 * MDP(2) + t325 * MDP(3) + (-0.2e1 * t399 + t417) * MDP(4) + t306 * MDP(5) + (t415 * pkin(1) + (t310 + t353) * qJ(2)) * MDP(6) + (qJDD(1) * t291 - 0.2e1 * t294 * t339) * MDP(7) + 0.2e1 * (-t294 * t351 + t356 * t371) * MDP(8) + (-t294 * t300 + t287) * MDP(9) + (-qJDD(3) * t294 - t297 * t300) * MDP(10) + (t294 * t304 + t297 * t313) * MDP(12) + (-t294 * t313 + t297 * t304) * MDP(13) + (-t222 * t380 - t261 * t311) * MDP(14) + ((t390 + t392) * t364 + (t398 - t397 + (-t389 + t393) * qJD(4)) * t297) * MDP(15) + ((-t357 * t277 - t222) * t294 + (t314 + t366) * t297) * MDP(16) + ((t277 * t365 - t223) * t294 + (-t315 - t367) * t297) * MDP(17) + (t256 * t294 + t363 * t277) * MDP(18) + (t241 * t277 + t254 * t256 + (t259 * t362 + t225 + (-t277 * t299 - t249) * t360) * t294 + (t220 * qJD(3) - t299 * t223 + t250 * t360) * t297 + ((-qJD(4) * t266 - t345) * t277 + t231 * t297 + (-qJD(3) * t250 - qJD(4) * t242 - t256 * t299 - t232) * t294) * t293 + t326) * MDP(19) + (-t348 * t277 - t372 * t256 + (t299 * t343 + (-t250 * t296 + t261 * t299) * qJD(3) - t330) * t294 + (-t221 * qJD(3) + t299 * t222 + t231 * t296 - t250 * t361) * t297 + t327) * MDP(20) + (t212 * t277 + t223 * t258 + t226 * t256 + t233 * t259 + (-t229 * t365 + t208) * t294 + (qJD(3) * t215 + t214 * t293 + t344) * t297 + t326) * MDP(21) + (-t213 * t277 - t222 * t258 - t230 * t256 + t233 * t261 + (-t229 * t357 - t209) * t294 + (-qJD(3) * t217 + t214 * t296 - t229 * t361) * t297 + t327) * MDP(22) + (-t212 * t261 - t213 * t259 + t222 * t226 - t223 * t230 + t322 * t364 + (qJD(4) * t321 - t208 * t296 - t209 * t293 + t325) * t297) * MDP(23) + (t209 * t230 + t217 * t213 + t208 * t226 + t215 * t212 + t214 * t258 + t229 * t233 - g(1) * (-t281 * t295 - t298 * t416) - g(2) * (t281 * t298 - t295 * t416)) * MDP(24); qJDD(1) * MDP(4) - t301 * MDP(5) + (-t402 - t415) * MDP(6) + t287 * MDP(12) + t414 * MDP(24) + (MDP(19) + MDP(21)) * (-t297 * t223 + (t367 - t396) * t294 + (-t293 * t363 - t296 * t331) * t277) + (MDP(20) + MDP(22)) * (t222 * t297 + (t366 - t395) * t294 + (t293 * t331 - t297 * t357) * t277) + t303 * qJD(1) + (t370 * MDP(13) - t214 * MDP(24) + ((t390 - t392) * MDP(23) - t321 * MDP(24)) * qJD(3)) * t297 + (t370 * MDP(12) - qJDD(3) * MDP(13) + (-t397 - t398) * MDP(23) + (qJD(3) * t229 + t323) * MDP(24) + t303 * qJD(4)) * t294; MDP(9) * t351 - MDP(10) * t352 + qJDD(3) * MDP(11) + (-t297 * t309 + t404) * MDP(12) + (t294 * t309 + t403) * MDP(13) + (t277 * t389 - t398) * MDP(14) + ((-t222 - t394) * t296 + (-t223 - t391) * t293) * MDP(15) + ((-t261 * t297 + t277 * t382) * qJD(1) + t315) * MDP(16) + ((-t277 * t293 * t294 + t259 * t297) * qJD(1) + t314) * MDP(17) - t277 * MDP(18) * t368 + (-t220 * t368 - t259 * t265 - pkin(3) * t223 - t248 * t277 + (t329 - t350) * t296 + (t277 * t387 + t308) * t293 + t374) * MDP(19) + (pkin(3) * t222 + t375 * t277 + t221 * t368 - t261 * t265 + t271 + t308 * t296 + (t319 - t329) * t293) * MDP(20) + (-t215 * t368 - t223 * t284 - t234 * t259 + t256 * t268 + (-t214 - t350) * t296 + t376 * t277 + (t229 * t369 + (t229 + t409) * qJD(4)) * t293 + t374) * MDP(21) + (t344 + t222 * t284 - t234 * t261 + t256 * t269 + t271 - t377 * t277 + (t217 * t297 + t229 * t382) * qJD(1) + (pkin(4) * t413 + t214 + t319) * t293) * MDP(22) + (-t403 + t222 * t268 + t223 * t269 - t376 * t261 - t377 * t259 - t322 * qJD(4) + (-qJD(1) * t322 + t414) * t294 + t323) * MDP(23) + (-t209 * t269 + t208 * t268 - t214 * t284 - g(3) * t332 + (pkin(4) * t361 - t234) * t229 + t377 * t217 + t376 * t215 + t414 * (t284 * t297 - t292 * t294)) * MDP(24) + (MDP(7) * t294 * t297 - MDP(8) * t371) * t301; t261 * t259 * MDP(14) + (-t255 + t411) * MDP(15) + (-t222 + t394) * MDP(16) + (-t223 + t391) * MDP(17) + t256 * MDP(18) + (t221 * t277 - t250 * t261 + t302) * MDP(19) + (t220 * t277 + t250 * t259 + t307) * MDP(20) + (0.2e1 * t410 + t401 + t217 * t277 + (-t229 + t335) * t261 + t302) * MDP(21) + (-pkin(4) * t411 + t400 + t216 * t277 + (qJD(5) + t229) * t259 + t307) * MDP(22) + (pkin(4) * t222 - t259 * t378) * MDP(23) + (t378 * t217 + (-t229 * t261 + t208 + t412) * pkin(4)) * MDP(24); (t293 * t354 - t373 + t391) * MDP(21) + (-t296 * t340 + t347 - t394) * MDP(22) + (-t255 - t411) * MDP(23) + (t215 * t261 + t217 * t259 + t312 - t404) * MDP(24) + ((qJDD(1) * t293 + t296 * t355) * MDP(21) - t293 * MDP(22) * t355 + t318 * MDP(24)) * t297;];
tau = t1;

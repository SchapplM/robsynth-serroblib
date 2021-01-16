% Calculate vector of inverse dynamics joint torques for
% S5RRPRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:13:28
% EndTime: 2021-01-15 21:13:39
% DurationCPUTime: 4.13s
% Computational Cost: add. (1829->348), mult. (4099->461), div. (0->0), fcn. (2792->10), ass. (0->168)
t344 = qJD(2) + qJD(4);
t351 = sin(qJ(4));
t355 = cos(qJ(2));
t453 = cos(qJ(4));
t399 = qJD(1) * t453;
t352 = sin(qJ(2));
t421 = qJD(1) * t352;
t460 = -t351 * t421 + t355 * t399;
t438 = t460 * t344;
t347 = qJ(2) + qJ(4);
t340 = sin(t347);
t353 = sin(qJ(1));
t356 = cos(qJ(1));
t382 = g(1) * t356 + g(2) * t353;
t461 = t382 * t340;
t295 = t351 * t355 + t453 * t352;
t269 = t344 * t295;
t393 = qJDD(1) * t453;
t411 = qJDD(1) * t352;
t253 = t269 * qJD(1) + t351 * t411 - t355 * t393;
t357 = pkin(1) + pkin(2);
t288 = qJD(5) - t460;
t390 = t288 ^ 2;
t459 = pkin(4) * t390;
t446 = pkin(3) + qJ(3);
t303 = t446 * t352;
t296 = qJD(1) * t303;
t349 = qJD(2) * pkin(1);
t280 = qJD(2) * pkin(2) - t296 + t349;
t304 = t446 * t355;
t297 = qJD(1) * t304;
t401 = t453 * t297;
t262 = t351 * t280 + t401;
t394 = qJD(2) * t446;
t286 = -qJD(3) * t352 - t355 * t394;
t348 = qJDD(2) * pkin(1);
t260 = qJDD(2) * pkin(2) + t286 * qJD(1) - qJDD(1) * t303 + t348;
t285 = qJD(3) * t355 - t352 * t394;
t264 = t285 * qJD(1) + qJDD(1) * t304;
t391 = -t453 * t260 + t351 * t264;
t238 = t262 * qJD(4) + t391;
t341 = cos(t347);
t448 = g(3) * t341;
t458 = t238 + t448;
t305 = t357 * t355;
t293 = -qJD(1) * t305 + qJD(3);
t333 = g(3) * t340;
t456 = -t293 * t460 + t382 * t341 + t333;
t370 = -t453 * t303 - t351 * t304;
t247 = t370 * qJD(4) + t453 * t285 + t351 * t286;
t251 = qJDD(5) + t253;
t402 = t453 * t280;
t432 = t351 * t297;
t261 = -t402 + t432;
t369 = -t351 * t352 + t453 * t355;
t268 = t344 * t369;
t276 = -t351 * t303 + t453 * t304;
t278 = -pkin(4) * t295 - t305;
t343 = qJDD(2) + qJDD(4);
t418 = qJD(4) * t351;
t380 = t453 * t264 - t297 * t418;
t398 = qJD(4) * t453;
t362 = t351 * t260 + t280 * t398 + t380;
t236 = t343 * pkin(4) + t362;
t420 = qJD(1) * t355;
t291 = -t351 * t420 - t352 * t399;
t267 = pkin(4) * t291 + t293;
t388 = qJD(5) * t267 + t236;
t455 = t238 * t295 - t276 * t251 + t261 * t268 - (qJD(5) * t278 + t247) * t288 + t388 * t369;
t346 = t355 ^ 2;
t454 = 0.2e1 * t346;
t452 = g(1) * t353;
t449 = g(2) * t356;
t447 = g(3) * t355;
t410 = qJDD(1) * t355;
t252 = t351 * t410 + t352 * t393 + t438;
t350 = sin(qJ(5));
t354 = cos(qJ(5));
t416 = qJD(5) * t354;
t405 = t354 * t252 + t350 * t343 + t344 * t416;
t417 = qJD(5) * t350;
t239 = t291 * t417 + t405;
t445 = t239 * t350;
t444 = t261 * t295;
t272 = -t291 * t350 - t354 * t344;
t443 = t272 * t288;
t442 = t272 * t291;
t275 = -t291 * t354 + t344 * t350;
t441 = t275 * t288;
t440 = t275 * t291;
t439 = t278 * t251;
t437 = t291 * t344;
t435 = t350 * t251;
t434 = t350 * t353;
t433 = t350 * t356;
t431 = t352 * t353;
t430 = t352 * t356;
t359 = qJD(1) ^ 2;
t429 = t352 * t359;
t428 = t353 * t354;
t427 = t353 * t355;
t249 = t354 * t251;
t426 = t354 * t356;
t425 = t355 * t356;
t424 = t355 * t359;
t298 = t357 * t421;
t299 = t357 * qJD(2) * t352;
t345 = t352 ^ 2;
t423 = t345 - t346;
t422 = t345 + t454;
t407 = pkin(1) * t420;
t314 = qJD(3) - t407;
t415 = -qJD(3) + t314;
t414 = qJD(1) * qJD(2);
t413 = qJD(1) * qJD(3);
t412 = qJDD(1) * t346;
t409 = qJDD(2) * t355;
t397 = t352 * t414;
t408 = pkin(1) * t397 + qJDD(3);
t406 = pkin(1) * t410;
t404 = g(1) * t425 + g(2) * t427 + g(3) * t352;
t403 = t357 * t453;
t396 = t355 * t414;
t292 = -t406 + t408;
t358 = qJD(2) ^ 2;
t392 = -qJ(3) * t358 - t292;
t389 = t354 * t288;
t271 = pkin(2) * t397 - qJDD(1) * t305 + t408;
t244 = -pkin(4) * t252 + t271;
t256 = t344 * pkin(4) + t262;
t387 = qJD(5) * t256 - t244;
t318 = t351 * t357 + pkin(4);
t385 = -pkin(4) * t460 + qJD(5) * t318 + t298;
t384 = t352 * t396;
t383 = g(1) * t430 + g(2) * t431 - t447;
t381 = -t449 + t452;
t265 = -t351 * t296 + t401;
t379 = t357 * t418 - t265;
t377 = -t251 * t318 - t261 * t460;
t376 = -t292 + t381;
t243 = t256 * t354 + t267 * t350;
t375 = -t243 * t291 + t261 * t416 + t458 * t350;
t242 = -t256 * t350 + t267 * t354;
t374 = t242 * t291 + t261 * t417 + t354 * t461;
t373 = t249 + (t350 * t460 - t417) * t288;
t266 = -t453 * t296 - t432;
t372 = -t357 * t398 + t266;
t368 = t268 * t354 - t295 * t417;
t367 = -pkin(4) * t251 + (-t288 - t460) * t261;
t366 = -qJ(3) * qJDD(1) + (-qJD(3) - t314) * qJD(1);
t302 = -qJ(3) * t421 + t349;
t364 = -qJD(2) * t302 * t355 - t382;
t363 = t293 * t291 - t391 - t448 + t461;
t323 = t354 * t343;
t240 = qJD(5) * t275 + t252 * t350 - t323;
t361 = ((t239 - t443) * t354 + (-t240 - t441) * t350) * MDP(23) + (t275 * t389 + t445) * MDP(22) + (t373 - t442) * MDP(25) + (t288 * t389 + t435 + t440) * MDP(24) + (t252 - t438) * MDP(17) + (-t253 - t437) * MDP(18) + (t291 ^ 2 - t460 ^ 2) * MDP(16) + t343 * MDP(19) + (MDP(15) * t460 + t288 * MDP(26)) * t291;
t360 = qJ(3) ^ 2;
t329 = g(1) * t427;
t328 = g(2) * t430;
t284 = t341 * t426 + t434;
t283 = -t341 * t433 + t428;
t282 = -t341 * t428 + t433;
t281 = t341 * t434 + t426;
t277 = -t352 * t413 + t348 + (-t396 - t411) * qJ(3);
t255 = -pkin(4) * t268 + t299;
t248 = t276 * qJD(4) + t351 * t285 - t453 * t286;
t241 = t354 * t244;
t1 = [(-g(1) * t281 - g(2) * t283 - t370 * t239 - t243 * t269 + t248 * t275 + (-(-qJD(5) * t276 + t255) * t288 - t439 - t387 * t369 - qJD(5) * t444) * t350 + t455 * t354) * MDP(28) + (-g(1) * t282 - g(2) * t284 - t370 * t240 - t241 * t369 + t242 * t269 + t248 * t272 + (t439 + t255 * t288 + (t256 * t369 - t276 * t288 + t444) * qJD(5)) * t354 + t455 * t350) * MDP(27) + (-t295 * t435 + t240 * t369 - t269 * t272 + (-t268 * t350 - t295 * t416) * t288) * MDP(25) + (-t269 * t344 + t343 * t369) * MDP(18) + (-t251 * t369 + t269 * t288) * MDP(26) + (-t239 * t369 + t295 * t249 + t269 * t275 + t368 * t288) * MDP(24) + (-t352 * t358 + t409) * MDP(7) + (-t247 * t344 - t252 * t305 + t268 * t293 + t271 * t295 - t276 * t343 - t291 * t299 - t381 * t340) * MDP(21) + t381 * MDP(2) + t382 * MDP(3) + (qJDD(1) * t345 + 0.2e1 * t384) * MDP(4) + (t239 * t295 * t354 + t368 * t275) * MDP(22) + (qJDD(2) * t352 + t355 * t358) * MDP(6) + (t268 * t344 + t295 * t343) * MDP(17) + (t252 * t295 - t268 * t291) * MDP(15) + ((-t272 * t354 - t275 * t350) * t268 + (-t445 - t240 * t354 + (t272 * t350 - t275 * t354) * qJD(5)) * t295) * MDP(23) + (t360 * t412 + t376 * t355 * pkin(1) + (t413 * t454 + t364) * qJ(3) + (-t277 * qJ(3) - t302 * qJD(3) + (pkin(1) * t314 - 0.2e1 * t360 * t420) * qJD(2)) * t352) * MDP(14) + (-t248 * t344 - t253 * t305 + t269 * t293 - t271 * t369 - t299 * t460 + t381 * t341 + t343 * t370) * MDP(20) + (t252 * t369 - t253 * t295 + t268 * t460 + t269 * t291) * MDP(16) + (pkin(1) * t412 + t329 + (t392 - t449) * t355 + (-qJ(3) * qJDD(2) + (-0.2e1 * t407 + t415) * qJD(2)) * t352) * MDP(11) + (-qJ(3) * t409 + t328 + (-t392 - t406 - t452) * t352 + (t423 * qJD(1) * pkin(1) + t415 * t355) * qJD(2)) * MDP(12) + 0.2e1 * (t352 * t410 - t423 * t414) * MDP(5) + (-g(2) * t425 + t329) * MDP(9) + (-g(1) * t431 + t328) * MDP(10) + (-t277 * t352 + t422 * t413 + (t422 * qJDD(1) - 0.2e1 * t384) * qJ(3) + t364) * MDP(13) + qJDD(1) * MDP(1); t423 * t359 * MDP(5) + (0.2e1 * t348 + (pkin(1) * t424 + t366) * t352 + t383) * MDP(11) + ((qJ(3) * qJD(1) * t302 + t360 * t429) * t355 + (-t447 + t277 + (-qJD(1) * t314 + t382) * t352) * pkin(1)) * MDP(14) + (-pkin(1) * t345 * t359 + t366 * t355 + t404) * MDP(12) + (-pkin(1) * t411 + (qJ(3) * t429 + (t302 - t349) * qJD(1)) * t355) * MDP(13) + (t266 * t344 + t298 * t291 + (-t343 * t357 - t260) * t351 + (-t344 * t403 - t402) * qJD(4) - t380 + t456) * MDP(21) + MDP(7) * t410 + (t343 * t403 + t265 * t344 + t298 * t460 + (-t401 + (-t344 * t357 - t280) * t351) * qJD(4) + t363) * MDP(20) + t383 * MDP(9) + t404 * MDP(10) + MDP(6) * t411 + t361 + qJDD(2) * MDP(8) - t352 * MDP(4) * t424 + (-t239 * t403 + t377 * t354 - t350 * t461 + t379 * t275 + (t385 * t350 + t372 * t354) * t288 + t375) * MDP(28) + (-t240 * t403 - t458 * t354 + t377 * t350 + t379 * t272 + (t372 * t350 - t385 * t354) * t288 + t374) * MDP(27); (0.2e1 * t397 - t410) * MDP(11) + (0.2e1 * t396 + t411) * MDP(12) + (-t345 - t346) * t359 * MDP(13) + (-qJ(3) * t346 * t359 + t302 * t421 - t376) * MDP(14) + (t253 - t437) * MDP(20) + (t252 + t438) * MDP(21) + (t373 + t442) * MDP(27) + (-t354 * t390 - t435 + t440) * MDP(28); (-t261 * t344 - t362 + t456) * MDP(21) + (-t262 * t272 + t367 * t350 + (-t458 - t459) * t354 + t374) * MDP(27) + (t363 + (-qJD(4) + t344) * t262) * MDP(20) + (-t262 * t275 + t367 * t354 + (-t461 + t459) * t350 + t375) * MDP(28) + t361; t275 * t272 * MDP(22) + (-t272 ^ 2 + t275 ^ 2) * MDP(23) + (t405 + t443) * MDP(24) + (t323 + t441) * MDP(25) + t251 * MDP(26) + (-g(1) * t283 + g(2) * t281 + t243 * t288 - t261 * t275 + t241) * MDP(27) + (g(1) * t284 - g(2) * t282 + t242 * t288 + t261 * t272) * MDP(28) + ((-t236 + t333) * MDP(28) + (MDP(25) * t291 - MDP(27) * t256 - MDP(28) * t267) * qJD(5)) * t354 + (qJD(5) * t291 * MDP(24) + (-qJD(5) * t344 - t252) * MDP(25) + (-t388 + t333) * MDP(27) + t387 * MDP(28)) * t350;];
tau = t1;

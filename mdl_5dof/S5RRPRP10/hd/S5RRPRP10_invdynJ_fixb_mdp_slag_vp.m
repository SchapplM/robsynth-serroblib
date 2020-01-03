% Calculate vector of inverse dynamics joint torques for
% S5RRPRP10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP10_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:14
% EndTime: 2019-12-31 20:11:18
% DurationCPUTime: 3.79s
% Computational Cost: add. (1701->382), mult. (3525->478), div. (0->0), fcn. (1997->6), ass. (0->186)
t348 = sin(qJ(2));
t404 = qJD(1) * qJD(2);
t390 = t348 * t404;
t351 = cos(qJ(2));
t402 = qJDD(1) * t351;
t471 = -t390 + t402;
t459 = pkin(3) + pkin(6);
t347 = sin(qJ(4));
t350 = cos(qJ(4));
t417 = qJD(2) * t347;
t420 = qJD(1) * t351;
t291 = t350 * t420 + t417;
t411 = qJD(4) * t291;
t255 = -t350 * qJDD(2) + t471 * t347 + t411;
t389 = t351 * t404;
t403 = qJDD(1) * t348;
t365 = t389 + t403;
t290 = qJDD(4) + t365;
t391 = t347 * t420;
t415 = qJD(2) * t350;
t293 = -t391 + t415;
t353 = -pkin(2) - pkin(7);
t333 = t348 * qJ(3);
t387 = -pkin(1) - t333;
t363 = t353 * t351 + t387;
t266 = t363 * qJD(1);
t421 = qJD(1) * t348;
t327 = pkin(6) * t421;
t464 = qJD(3) + t327;
t405 = pkin(3) * t421 + t464;
t269 = t353 * qJD(2) + t405;
t249 = t266 * t350 + t269 * t347;
t317 = pkin(2) * t390;
t447 = qJ(3) * t351;
t377 = pkin(7) * t348 - t447;
t413 = qJD(3) * t348;
t360 = t377 * qJD(2) - t413;
t247 = t360 * qJD(1) + t363 * qJDD(1) + t317;
t316 = pkin(6) * t389;
t324 = pkin(6) * t403;
t388 = qJDD(3) + t316 + t324;
t258 = t365 * pkin(3) + t353 * qJDD(2) + t388;
t384 = -t247 * t347 + t350 * t258;
t358 = -t249 * qJD(4) + t384;
t237 = pkin(4) * t290 + qJ(5) * t255 - qJD(5) * t293 + t358;
t244 = -qJ(5) * t291 + t249;
t318 = qJD(4) + t421;
t470 = t318 * t244 + t237;
t256 = -qJD(4) * t391 + qJDD(2) * t347 + (qJD(2) * qJD(4) + t471) * t350;
t409 = qJD(4) * t350;
t397 = -t350 * t247 - t347 * t258 - t269 * t409;
t410 = qJD(4) * t347;
t238 = -qJ(5) * t256 - qJD(5) * t291 - t266 * t410 - t397;
t248 = -t266 * t347 + t350 * t269;
t243 = -qJ(5) * t293 + t248;
t242 = pkin(4) * t318 + t243;
t469 = -t318 * t242 + t238;
t443 = t293 * t318;
t468 = -t256 + t443;
t337 = t351 * pkin(2);
t423 = t337 + t333;
t305 = -pkin(1) - t423;
t287 = -pkin(7) * t351 + t305;
t307 = t459 * t348;
t424 = t350 * t287 + t347 * t307;
t276 = t350 * t290;
t465 = -t318 * t410 + t276;
t349 = sin(qJ(1));
t352 = cos(qJ(1));
t378 = g(1) * t352 + g(2) * t349;
t463 = pkin(4) * t256 + qJDD(5);
t432 = t350 * t352;
t279 = -t347 * t349 + t348 * t432;
t435 = t349 * t350;
t281 = t347 * t352 + t348 * t435;
t433 = t350 * t351;
t462 = -g(1) * t279 - g(2) * t281 + g(3) * t433;
t328 = pkin(6) * t420;
t300 = pkin(3) * t420 + t328;
t343 = qJD(2) * qJ(3);
t277 = t343 + t300;
t461 = t318 * t277 + t353 * t290;
t460 = t293 ^ 2;
t457 = pkin(4) * t347;
t455 = g(1) * t349;
t451 = g(2) * t352;
t450 = g(3) * t348;
t340 = g(3) * t351;
t448 = pkin(6) * qJDD(2);
t446 = qJDD(2) * pkin(2);
t445 = t255 * t350;
t444 = t291 * t318;
t346 = -qJ(5) - pkin(7);
t442 = t346 * t351;
t441 = t347 * t348;
t440 = t347 * t351;
t439 = t348 * t349;
t438 = t348 * t350;
t437 = t348 * t352;
t355 = qJD(1) ^ 2;
t436 = t348 * t355;
t434 = t349 * t351;
t431 = t351 * t352;
t429 = qJ(5) - t353;
t428 = -t243 + t242;
t331 = pkin(2) * t421;
t273 = t377 * qJD(1) + t331;
t427 = t350 * t273 + t347 * t300;
t285 = t350 * t300;
t426 = -qJD(5) * t350 + t429 * t410 + t273 * t347 - t285 - (pkin(4) * t351 - qJ(5) * t441) * qJD(1);
t303 = t429 * t350;
t425 = -qJ(5) * t350 * t421 - qJD(4) * t303 - qJD(5) * t347 - t427;
t308 = t459 * t351;
t344 = t348 ^ 2;
t345 = t351 ^ 2;
t422 = t344 - t345;
t419 = qJD(2) * t291;
t418 = qJD(2) * t293;
t416 = qJD(2) * t348;
t414 = qJD(2) * t351;
t412 = qJD(4) * t266;
t408 = qJD(4) * t351;
t407 = qJD(4) * t353;
t406 = qJD(5) * t351;
t401 = pkin(4) * t440;
t399 = t347 * t437;
t398 = t351 * t436;
t325 = pkin(6) * t402;
t341 = qJDD(2) * qJ(3);
t342 = qJD(2) * qJD(3);
t396 = t325 + t341 + t342;
t395 = -g(1) * t437 - g(2) * t439 + t340;
t323 = pkin(4) * t350 + pkin(3);
t394 = qJD(2) * t459;
t392 = t347 * t408;
t386 = -qJD(2) * pkin(2) + qJD(3);
t385 = qJ(5) * t351 - t287;
t383 = t352 * pkin(1) + pkin(2) * t431 + t349 * pkin(6) + qJ(3) * t437;
t382 = -t324 - t395;
t381 = pkin(3) * t402 + t396;
t380 = qJD(1) * t394;
t354 = qJD(2) ^ 2;
t379 = pkin(6) * t354 + t451;
t304 = t327 + t386;
t306 = -t328 - t343;
t375 = t304 * t351 + t306 * t348;
t374 = t318 ^ 2;
t372 = t387 - t337;
t369 = -0.2e1 * pkin(1) * t404 - t448;
t278 = t372 * qJD(1);
t368 = t278 * t421 + qJDD(3) - t382;
t367 = -t290 * t347 - t318 * t409;
t366 = -qJ(3) * t414 - t413;
t330 = pkin(2) * t416;
t263 = t330 + t360;
t301 = t459 * t414;
t364 = t350 * t263 - t287 * t410 + t347 * t301 + t307 * t409;
t362 = 0.2e1 * qJDD(1) * pkin(1) - t379;
t361 = t448 + (-qJD(1) * t305 - t278) * qJD(2);
t259 = -t348 * t380 + t381;
t359 = -t378 * t351 + t259 - t450;
t257 = t366 * qJD(1) + t372 * qJDD(1) + t317;
t275 = t330 + t366;
t357 = qJD(1) * t275 + qJDD(1) * t305 + t257 + t379;
t267 = pkin(6) * t390 - t396;
t272 = t388 - t446;
t356 = t375 * qJD(2) - t267 * t351 + t272 * t348;
t338 = t352 * pkin(6);
t321 = g(1) * t434;
t315 = qJ(3) * t431;
t313 = qJ(3) * t434;
t302 = t429 * t347;
t299 = t348 * t394;
t297 = -qJ(3) * t420 + t331;
t296 = t350 * t307;
t289 = t291 ^ 2;
t286 = t350 * t301;
t282 = -t347 * t439 + t432;
t280 = t399 + t435;
t260 = pkin(4) * t291 + qJD(5) + t277;
t252 = -qJ(5) * t433 + t424;
t251 = pkin(4) * t348 + t385 * t347 + t296;
t241 = t259 + t463;
t240 = -t350 * t406 + (t348 * t415 + t392) * qJ(5) + t364;
t239 = pkin(4) * t414 + t286 + t385 * t409 + (-qJ(5) * t416 - qJD(4) * t307 - t263 + t406) * t347;
t1 = [qJDD(1) * MDP(1) + (-t451 + t455) * MDP(2) + t378 * MDP(3) + (qJDD(1) * t344 + 0.2e1 * t348 * t389) * MDP(4) + 0.2e1 * (t348 * t402 - t422 * t404) * MDP(5) + (qJDD(2) * t348 + t351 * t354) * MDP(6) + (qJDD(2) * t351 - t348 * t354) * MDP(7) + (t369 * t348 + t362 * t351 + t321) * MDP(9) + (t369 * t351 + (-t362 - t455) * t348) * MDP(10) + ((t344 + t345) * qJDD(1) * pkin(6) + t356 - t378) * MDP(11) + (t361 * t348 + t357 * t351 - t321) * MDP(12) + (t361 * t351 + (-t357 + t455) * t348) * MDP(13) + (t356 * pkin(6) - g(1) * t338 - g(2) * t383 + t257 * t305 + t278 * t275 - t372 * t455) * MDP(14) + (t255 * t440 + (t347 * t416 - t350 * t408) * t293) * MDP(15) + ((-t291 * t347 + t293 * t350) * t416 + (t445 + t256 * t347 + (t291 * t350 + t293 * t347) * qJD(4)) * t351) * MDP(16) + ((t318 * t417 - t255) * t348 + (t367 + t418) * t351) * MDP(17) + ((t318 * t415 - t256) * t348 + (-t419 - t465) * t351) * MDP(18) + (t290 * t348 + t318 * t414) * MDP(19) + ((-t263 * t347 + t286) * t318 + (-t287 * t347 + t296) * t290 + t384 * t348 - t299 * t291 + t308 * t256 + t259 * t433 - g(1) * t282 - g(2) * t280 + (t248 * t351 - t277 * t438) * qJD(2) + (-t249 * t348 - t277 * t440 - t424 * t318) * qJD(4)) * MDP(20) + (-t364 * t318 - t424 * t290 - t299 * t293 - t308 * t255 + g(1) * t281 - g(2) * t279 + ((qJD(2) * t277 + t412) * t347 + t397) * t348 + (-qJD(2) * t249 - t259 * t347 - t277 * t409) * t351) * MDP(21) + (-t239 * t293 - t240 * t291 + t251 * t255 - t252 * t256 + t321 + (-t242 * t347 + t244 * t350) * t416 + (-t451 + t237 * t347 - t238 * t350 + (t242 * t350 + t244 * t347) * qJD(4)) * t351) * MDP(22) + (t238 * t252 + t244 * t240 + t237 * t251 + t242 * t239 + t241 * (pkin(4) * t433 + t308) - g(1) * (t323 * t352 + t338) - g(2) * (pkin(4) * t399 - t346 * t431 + t383) + (-g(1) * (-pkin(4) * t441 + t372 + t442) - g(2) * t323) * t349 + (-pkin(4) * t392 + (-pkin(6) - t323) * t416) * t260) * MDP(23); -MDP(4) * t398 + t422 * t355 * MDP(5) + MDP(6) * t403 + MDP(7) * t402 + qJDD(2) * MDP(8) + (pkin(1) * t436 + t382) * MDP(9) + (t450 - t325 + (pkin(1) * t355 + t378) * t351) * MDP(10) + ((-pkin(2) * t348 + t447) * qJDD(1) + ((-t306 - t343) * t348 + (-t304 + t386) * t351) * qJD(1)) * MDP(11) + (-t297 * t420 + t368 - 0.2e1 * t446) * MDP(12) + (t325 + 0.2e1 * t341 + 0.2e1 * t342 + (qJD(1) * t297 - g(3)) * t348 + (qJD(1) * t278 - t378) * t351) * MDP(13) + (-t267 * qJ(3) - t306 * qJD(3) - t272 * pkin(2) - t278 * t297 - g(1) * (-pkin(2) * t437 + t315) - g(2) * (-pkin(2) * t439 + t313) - g(3) * t423 - t375 * qJD(1) * pkin(6)) * MDP(14) + (-t347 * t443 - t445) * MDP(15) + ((-t256 - t443) * t350 + (t255 + t444) * t347) * MDP(16) + ((-t293 * t351 - t318 * t441) * qJD(1) + t465) * MDP(17) + ((t291 * t351 - t318 * t438) * qJD(1) + t367) * MDP(18) - t318 * MDP(19) * t420 + (-t248 * t420 + qJ(3) * t256 - t285 * t318 + t405 * t291 + t461 * t350 + ((t273 - t407) * t318 + t359) * t347) * MDP(20) + (-qJ(3) * t255 + t427 * t318 + t249 * t420 + t405 * t293 - t461 * t347 + (-t318 * t407 + t359) * t350) * MDP(21) + (-t255 * t303 + t256 * t302 - t425 * t291 - t426 * t293 - t469 * t347 - t470 * t350 - t395) * MDP(22) + (-t238 * t302 - t237 * t303 + t241 * (qJ(3) + t457) - g(1) * (t352 * t401 + t315) - g(2) * (t349 * t401 + t313) - g(3) * (t423 - t442) + (pkin(4) * t409 + t464) * t260 + t425 * t244 + t426 * t242 + (t260 * t323 * qJD(1) - g(3) * t457 + t378 * (pkin(2) - t346)) * t348) * MDP(23); MDP(11) * t403 + (qJDD(2) + t398) * MDP(12) + (-t344 * t355 - t354) * MDP(13) + (qJD(2) * t306 + t316 + t368 - t446) * MDP(14) + (t276 - t419) * MDP(20) - MDP(21) * t418 + (-qJD(2) * t260 + t395) * MDP(23) + ((-t291 * t421 + t255 - t411) * MDP(22) + t470 * MDP(23) - MDP(21) * t374) * t350 + (-MDP(20) * t374 - t290 * MDP(21) + t468 * MDP(22) + t469 * MDP(23)) * t347; t293 * t291 * MDP(15) + (-t289 + t460) * MDP(16) + (t444 - t255) * MDP(17) + t468 * MDP(18) + t290 * MDP(19) + (t249 * t318 - t277 * t293 + t358 + t462) * MDP(20) + (g(1) * t280 - g(2) * t282 + t248 * t318 + t277 * t291 + (t412 - t340) * t347 + t397) * MDP(21) + (pkin(4) * t255 - t428 * t291) * MDP(22) + (t428 * t244 + (-t260 * t293 + t237 + t462) * pkin(4)) * MDP(23); (-t289 - t460) * MDP(22) + (-g(1) * t431 - g(2) * t434 + t242 * t293 + t244 * t291 + t381 + (-g(3) - t380) * t348 + t463) * MDP(23);];
tau = t1;

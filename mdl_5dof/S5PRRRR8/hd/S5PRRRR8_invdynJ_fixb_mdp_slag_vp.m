% Calculate vector of inverse dynamics joint torques for
% S5PRRRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:56
% EndTime: 2019-12-05 17:17:06
% DurationCPUTime: 4.76s
% Computational Cost: add. (2423->362), mult. (5533->506), div. (0->0), fcn. (4412->14), ass. (0->169)
t354 = qJD(3) + qJD(4);
t362 = sin(qJ(4));
t366 = cos(qJ(3));
t469 = cos(qJ(4));
t418 = qJD(2) * t469;
t363 = sin(qJ(3));
t437 = qJD(2) * t363;
t475 = -t362 * t437 + t366 * t418;
t476 = t475 * t354;
t470 = pkin(7) + pkin(8);
t353 = qJDD(3) + qJDD(4);
t364 = sin(qJ(2));
t359 = sin(pkin(5));
t439 = qJD(1) * t359;
t423 = t364 * t439;
t409 = qJD(2) * t470 + t423;
t360 = cos(pkin(5));
t438 = qJD(1) * t360;
t296 = t363 * t438 + t409 * t366;
t428 = t360 * qJDD(1);
t334 = t366 * t428;
t367 = cos(qJ(2));
t432 = qJD(1) * qJD(2);
t306 = qJDD(2) * pkin(7) + (qJDD(1) * t364 + t367 * t432) * t359;
t408 = pkin(8) * qJDD(2) + t306;
t258 = qJDD(3) * pkin(3) - qJD(3) * t296 - t408 * t363 + t334;
t295 = -t409 * t363 + t366 * t438;
t263 = qJD(3) * t295 + t363 * t428 + t408 * t366;
t464 = qJD(3) * pkin(3);
t290 = t295 + t464;
t424 = t469 * t296;
t266 = t362 * t290 + t424;
t471 = t266 * qJD(4) - t469 * t258 + t362 * t263;
t243 = -t353 * pkin(4) + t471;
t463 = cos(pkin(10));
t411 = t463 * t364;
t358 = sin(pkin(10));
t449 = t358 * t367;
t312 = t360 * t411 + t449;
t410 = t463 * t367;
t450 = t358 * t364;
t314 = -t360 * t450 + t410;
t357 = qJ(3) + qJ(4);
t351 = sin(t357);
t352 = cos(t357);
t412 = t359 * t463;
t448 = t359 * t364;
t451 = t358 * t359;
t383 = -g(3) * (-t351 * t448 + t352 * t360) - g(2) * (-t312 * t351 - t352 * t412) - g(1) * (-t314 * t351 + t352 * t451);
t379 = -t243 + t383;
t444 = t362 * t366;
t321 = t469 * t363 + t444;
t292 = t354 * t321;
t413 = qJDD(2) * t469;
t430 = qJDD(2) * t363;
t396 = t362 * t430 - t366 * t413;
t275 = qJD(2) * t292 + t396;
t272 = qJDD(5) + t275;
t350 = -pkin(3) * t366 - pkin(2);
t386 = -t362 * t363 + t469 * t366;
t283 = -pkin(4) * t386 - pkin(9) * t321 + t350;
t311 = -t360 * t410 + t450;
t313 = t360 * t449 + t411;
t400 = g(1) * t313 + g(2) * t311;
t446 = t359 * t367;
t380 = g(3) * t446 - t400;
t377 = t380 * t352;
t474 = t283 * t272 - t377;
t382 = t363 * t464 - t423;
t417 = qJD(4) * t469;
t436 = qJD(4) * t362;
t374 = t362 * t258 + t469 * t263 + t290 * t417 - t296 * t436;
t242 = t353 * pkin(9) + t374;
t445 = t362 * t296;
t265 = t469 * t290 - t445;
t259 = -t354 * pkin(4) - t265;
t422 = t367 * t439;
t309 = t350 * qJD(2) - t422;
t319 = -qJD(2) * t444 - t363 * t418;
t273 = -pkin(4) * t475 + pkin(9) * t319 + t309;
t291 = t354 * t386;
t328 = t470 * t363;
t329 = t470 * t366;
t303 = -t362 * t328 + t469 * t329;
t310 = qJD(5) - t475;
t399 = g(1) * t314 + g(2) * t312;
t425 = qJD(3) * t470;
t322 = t363 * t425;
t323 = t366 * t425;
t387 = -t469 * t328 - t362 * t329;
t442 = -t387 * qJD(4) + t469 * t322 + t362 * t323 + t386 * t422;
t473 = (qJD(5) * t273 + t242) * t386 + t243 * t321 + t259 * t291 + (-qJD(5) * t283 + t442) * t310 - t303 * t272 - g(3) * t448 - t399;
t368 = qJD(3) ^ 2;
t416 = t364 * t432;
t397 = -qJDD(1) * t446 + t359 * t416;
t472 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t368 + (-g(3) * t367 + t416) * t359 - t397 + t400;
t465 = qJD(2) * pkin(2);
t429 = qJDD(2) * t366;
t274 = t362 * t429 + t363 * t413 + t476;
t361 = sin(qJ(5));
t365 = cos(qJ(5));
t434 = qJD(5) * t365;
t426 = t365 * t274 + t361 * t353 + t354 * t434;
t435 = qJD(5) * t361;
t253 = t319 * t435 + t426;
t461 = t253 * t361;
t460 = t259 * t475;
t459 = t272 * t361;
t458 = t272 * t365;
t315 = t360 * t366 - t363 * t448;
t447 = t359 * t366;
t316 = t360 * t363 + t364 * t447;
t279 = t362 * t315 + t469 * t316;
t457 = t279 * t272;
t453 = t319 * t361;
t299 = -t365 * t354 - t453;
t455 = t299 * t310;
t392 = t319 * t365 - t354 * t361;
t454 = t392 * t310;
t452 = t321 * t365;
t443 = qJDD(1) - g(3);
t441 = t303 * qJD(4) - t321 * t422 - t362 * t322 + t469 * t323;
t355 = t363 ^ 2;
t440 = -t366 ^ 2 + t355;
t433 = t259 * qJD(5);
t431 = qJD(2) * qJD(3);
t420 = qJD(2) * t446;
t415 = t363 * t431;
t414 = t366 * t431;
t260 = t354 * pkin(9) + t266;
t393 = t260 * t361 - t273 * t365;
t407 = -t319 * t393 + t361 * t433;
t405 = t274 * t361 - t365 * t353;
t404 = t310 * t365;
t282 = -pkin(4) * t319 - pkin(9) * t475;
t348 = pkin(3) * t362 + pkin(9);
t402 = pkin(3) * t437 + qJD(5) * t348 + t282;
t267 = t362 * t295 + t424;
t401 = pkin(3) * t436 - t267;
t398 = pkin(4) * t292 - pkin(9) * t291 + t382;
t394 = -t348 * t272 - t460;
t250 = t260 * t365 + t273 * t361;
t268 = t469 * t295 - t445;
t391 = -pkin(3) * t417 + t268;
t390 = -g(1) * t358 + t463 * g(2);
t388 = t469 * t315 - t362 * t316;
t385 = t310 * t435 - t458;
t384 = t291 * t365 - t321 * t435;
t381 = -t250 * t319 - t361 * t379 + t365 * t433;
t325 = -t422 - t465;
t378 = -qJD(2) * t325 - t306 + t399;
t284 = pkin(3) * t415 + t350 * qJDD(2) + t397;
t375 = -pkin(7) * qJDD(3) + (t325 + t422 - t465) * qJD(3);
t254 = -t392 * qJD(5) + t405;
t373 = ((t253 - t455) * t365 + (-t254 + t454) * t361) * MDP(20) + (-t392 * t404 + t461) * MDP(19) + (-t310 ^ 2 * t361 - t299 * t319 + t458) * MDP(22) + (t310 * t404 - t319 * t392 + t459) * MDP(21) + (t274 - t476) * MDP(14) + (-t396 + (-qJD(2) * t321 - t319) * t354) * MDP(15) + (t319 ^ 2 - t475 ^ 2) * MDP(13) + t353 * MDP(16) + (MDP(12) * t475 + t310 * MDP(23)) * t319;
t372 = t309 * t319 + t383 - t471;
t286 = t312 * t352 - t351 * t412;
t288 = t314 * t352 + t351 * t451;
t308 = t351 * t360 + t352 * t448;
t371 = g(1) * t288 + g(2) * t286 + g(3) * t308 - t309 * t475 - t374;
t369 = qJD(2) ^ 2;
t349 = -t469 * pkin(3) - pkin(4);
t294 = -t316 * qJD(3) - t363 * t420;
t293 = t315 * qJD(3) + t366 * t420;
t252 = t279 * qJD(4) + t362 * t293 - t469 * t294;
t251 = t388 * qJD(4) + t469 * t293 + t362 * t294;
t248 = pkin(4) * t275 - pkin(9) * t274 + t284;
t247 = t365 * t248;
t1 = [t443 * MDP(1) + (qJD(3) * t294 + qJDD(3) * t315) * MDP(10) + (-qJD(3) * t293 - qJDD(3) * t316) * MDP(11) + (-t252 * t354 + t353 * t388) * MDP(17) + (-t251 * t354 - t279 * t353) * MDP(18) + ((-t251 * t361 - t279 * t434) * t310 - t361 * t457 + t252 * t299 - t388 * t254) * MDP(24) + (-(t251 * t365 - t279 * t435) * t310 - t365 * t457 - t252 * t392 - t388 * t253) * MDP(25) + ((-qJDD(2) * MDP(4) + (-t366 * MDP(10) + t363 * MDP(11) - MDP(3)) * t369 + (-MDP(17) * t475 - MDP(18) * t319 + (MDP(24) * t365 - MDP(25) * t361) * t310) * qJD(2)) * t364 + (qJDD(2) * MDP(3) - t369 * MDP(4) + (-t415 + t429) * MDP(10) + (-t414 - t430) * MDP(11) - t275 * MDP(17) - t274 * MDP(18) + t385 * MDP(24) + (t310 * t434 + t459) * MDP(25)) * t367) * t359; qJDD(2) * MDP(2) + (t443 * t446 + t400) * MDP(3) + (-t443 * t448 + t399) * MDP(4) + (qJDD(2) * t355 + 0.2e1 * t363 * t414) * MDP(5) + 0.2e1 * (t363 * t429 - t440 * t431) * MDP(6) + (qJDD(3) * t363 + t366 * t368) * MDP(7) + (qJDD(3) * t366 - t363 * t368) * MDP(8) + (t375 * t363 + t366 * t472) * MDP(10) + (-t363 * t472 + t375 * t366) * MDP(11) + (t274 * t321 - t291 * t319) * MDP(12) + (t274 * t386 - t275 * t321 + t291 * t475 + t292 * t319) * MDP(13) + (t291 * t354 + t321 * t353) * MDP(14) + (-t292 * t354 + t353 * t386) * MDP(15) + (t275 * t350 - t284 * t386 + t292 * t309 + t353 * t387 - t441 * t354 - t382 * t475 - t377) * MDP(17) + (t274 * t350 + t284 * t321 + t291 * t309 - t303 * t353 - t382 * t319 + t380 * t351 + t442 * t354) * MDP(18) + (t253 * t452 - t384 * t392) * MDP(19) + ((-t299 * t365 + t361 * t392) * t291 + (-t461 - t254 * t365 + (t299 * t361 + t365 * t392) * qJD(5)) * t321) * MDP(20) + (-t253 * t386 + t272 * t452 - t292 * t392 + t310 * t384) * MDP(21) + (-t321 * t459 + t254 * t386 - t292 * t299 + (-t291 * t361 - t321 * t434) * t310) * MDP(22) + (-t272 * t386 + t292 * t310) * MDP(23) + (-t247 * t386 - t393 * t292 - t387 * t254 + t441 * t299 + (t398 * t310 + (t259 * t321 + t260 * t386 - t303 * t310) * qJD(5) + t474) * t365 + t473 * t361) * MDP(24) + (-t250 * t292 - t387 * t253 - t441 * t392 + ((-qJD(5) * t260 + t248) * t386 - t321 * t433 + (qJD(5) * t303 - t398) * t310 - t474) * t361 + t473 * t365) * MDP(25); t373 + MDP(8) * t429 + MDP(7) * t430 + (t349 * t253 + t394 * t365 - t401 * t392 + (t402 * t361 + t391 * t365) * t310 + t381) * MDP(25) + (-g(3) * t315 + t378 * t363 + t390 * t447 + t334) * MDP(10) + (g(3) * t316 + (-t390 * t359 - t428) * t363 + t378 * t366) * MDP(11) + (t268 * t354 + (t319 * t437 - t353 * t362 - t354 * t417) * pkin(3) + t371) * MDP(18) + (t267 * t354 + (t469 * t353 - t354 * t436 + t437 * t475) * pkin(3) + t372) * MDP(17) + (t349 * t254 + t401 * t299 + (t391 * t310 + t394) * t361 + (-t402 * t310 + t379) * t365 + t407) * MDP(24) + qJDD(3) * MDP(9) + (-t363 * t366 * MDP(5) + t440 * MDP(6)) * t369; t373 + (t265 * t354 + t371) * MDP(18) + (t266 * t354 + t372) * MDP(17) + (-pkin(4) * t253 + (t265 * t365 + t282 * t361) * t310 + t266 * t392 - t365 * t460 + t385 * pkin(9) + t381) * MDP(25) + (-pkin(4) * t254 - t266 * t299 + (-pkin(9) * t272 + t265 * t310 - t460) * t361 + ((-pkin(9) * qJD(5) - t282) * t310 + t379) * t365 + t407) * MDP(24); -t392 * t299 * MDP(19) + (-t299 ^ 2 + t392 ^ 2) * MDP(20) + (t426 + t455) * MDP(21) + (-t405 - t454) * MDP(22) + t272 * MDP(23) + (-t361 * t242 + t247 + t250 * t310 + t259 * t392 - g(1) * (-t288 * t361 + t313 * t365) - g(2) * (-t286 * t361 + t311 * t365) - g(3) * (-t308 * t361 - t365 * t446)) * MDP(24) + (-t365 * t242 - t361 * t248 - t393 * t310 + t259 * t299 - g(1) * (-t288 * t365 - t313 * t361) - g(2) * (-t286 * t365 - t311 * t361) - g(3) * (-t308 * t365 + t361 * t446)) * MDP(25) + (MDP(21) * t453 + t392 * MDP(22) - t250 * MDP(24) + t393 * MDP(25)) * qJD(5);];
tau = t1;

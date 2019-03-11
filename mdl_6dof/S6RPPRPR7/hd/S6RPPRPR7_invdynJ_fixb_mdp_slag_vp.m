% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:53
% EndTime: 2019-03-09 01:54:02
% DurationCPUTime: 6.90s
% Computational Cost: add. (4096->460), mult. (8323->583), div. (0->0), fcn. (6164->14), ass. (0->194)
t418 = sin(pkin(9));
t420 = cos(pkin(9));
t424 = sin(qJ(4));
t427 = cos(qJ(4));
t375 = t418 * t427 + t420 * t424;
t364 = t375 * qJD(1);
t362 = qJD(6) + t364;
t480 = qJD(1) * t424;
t468 = t418 * t480;
t479 = qJD(1) * t427;
t470 = t420 * t479;
t366 = -t468 + t470;
t417 = sin(pkin(10));
t419 = cos(pkin(10));
t343 = -t419 * qJD(4) + t366 * t417;
t426 = cos(qJ(6));
t345 = qJD(4) * t417 + t366 * t419;
t423 = sin(qJ(6));
t499 = t345 * t423;
t526 = -t426 * t343 - t499;
t527 = t362 * t526;
t422 = -pkin(1) - qJ(3);
t512 = -qJD(1) * qJD(3) + qJDD(1) * t422;
t377 = qJDD(2) + t512;
t463 = -pkin(7) * qJDD(1) + t377;
t349 = t463 * t418;
t350 = t463 * t420;
t385 = qJD(1) * t422 + qJD(2);
t466 = -pkin(7) * qJD(1) + t385;
t359 = t466 * t418;
t360 = t466 * t420;
t477 = qJD(4) * t427;
t478 = qJD(4) * t424;
t443 = -t424 * t349 + t350 * t427 - t359 * t477 - t360 * t478;
t285 = -qJDD(4) * pkin(4) + qJDD(5) - t443;
t413 = pkin(9) + qJ(4);
t401 = sin(t413);
t403 = cos(t413);
t425 = sin(qJ(1));
t428 = cos(qJ(1));
t518 = g(1) * t425 - g(2) * t428;
t435 = -g(3) * t401 + t403 * t518;
t521 = t285 + t435;
t447 = t343 * t423 - t345 * t426;
t525 = t362 * t447;
t376 = t417 * t426 + t419 * t423;
t368 = t376 * qJD(6);
t484 = t376 * t364 + t368;
t520 = -t424 * t359 + t360 * t427;
t309 = -qJD(4) * pkin(4) + qJD(5) - t520;
t369 = -t418 * t477 - t420 * t478;
t487 = t427 * t420;
t374 = t418 * t424 - t487;
t432 = -t285 * t374 + t309 * t369 + t518;
t439 = t362 * t376;
t481 = t418 ^ 2 + t420 ^ 2;
t523 = t385 * t481;
t459 = g(1) * t428 + g(2) * t425;
t442 = t459 * t401;
t522 = t459 * t403;
t414 = qJDD(1) * qJ(2);
t415 = qJD(1) * qJD(2);
t517 = t414 + t415;
t383 = qJDD(3) + t517;
t519 = t383 - t459;
t516 = t418 * MDP(7) + t420 * MDP(8);
t373 = t417 * t423 - t426 * t419;
t515 = t373 * qJD(6);
t485 = -t373 * t364 - t515;
t469 = t420 * t477;
t384 = qJD(4) * t468;
t513 = -t375 * qJDD(1) + t384;
t335 = qJD(1) * t469 - t513;
t333 = qJDD(6) + t335;
t502 = t333 * t376;
t514 = -t362 * t485 - t502;
t508 = g(3) * t403;
t433 = -t401 * t518 - t508;
t511 = 0.2e1 * t415;
t510 = pkin(8) * t419;
t404 = t418 * pkin(3);
t507 = -pkin(7) + t422;
t506 = pkin(8) + qJ(5);
t505 = pkin(1) * qJDD(1);
t504 = t526 * t366;
t503 = t447 * t366;
t501 = t335 * t417;
t500 = t335 * t419;
t497 = t364 * t417;
t496 = t369 * t417;
t308 = t373 * t333;
t495 = t374 * t417;
t494 = t374 * t419;
t412 = pkin(10) + qJ(6);
t400 = sin(t412);
t493 = t400 * t425;
t492 = t400 * t428;
t402 = cos(t412);
t491 = t402 * t425;
t490 = t402 * t428;
t391 = qJ(2) + t404;
t446 = t349 * t427 + t350 * t424;
t283 = qJDD(4) * qJ(5) + (qJD(5) + t520) * qJD(4) + t446;
t386 = qJDD(1) * t487;
t472 = qJDD(1) * t418;
t456 = -t424 * t472 + t386;
t334 = -qJD(4) * t364 + t456;
t372 = pkin(3) * t472 + t383;
t284 = pkin(4) * t335 - qJ(5) * t334 - qJD(5) * t366 + t372;
t267 = t419 * t283 + t417 * t284;
t370 = -t418 * t478 + t469;
t303 = pkin(4) * t370 - qJ(5) * t369 + qJD(5) * t374 + qJD(2);
t378 = t507 * t418;
t379 = t507 * t420;
t338 = t378 * t424 - t379 * t427;
t310 = -qJD(3) * t375 - qJD(4) * t338;
t278 = t417 * t303 + t419 * t310;
t322 = t427 * t359 + t424 * t360;
t314 = qJD(4) * qJ(5) + t322;
t397 = qJD(1) * qJ(2) + qJD(3);
t380 = qJD(1) * t404 + t397;
t318 = pkin(4) * t364 - qJ(5) * t366 + t380;
t287 = t419 * t314 + t417 * t318;
t332 = pkin(4) * t366 + qJ(5) * t364;
t290 = t417 * t332 + t419 * t520;
t331 = pkin(4) * t375 + qJ(5) * t374 + t391;
t339 = t378 * t427 + t379 * t424;
t295 = t417 * t331 + t419 * t339;
t486 = t369 * qJD(4) - t374 * qJDD(4);
t483 = t428 * pkin(1) + t425 * qJ(2);
t476 = qJD(6) * t423;
t475 = qJD(6) * t426;
t474 = -qJD(5) + t309;
t323 = -t419 * qJDD(4) + t334 * t417;
t324 = qJDD(4) * t417 + t334 * t419;
t471 = -t423 * t323 + t426 * t324 - t343 * t475;
t467 = g(2) * t483;
t465 = t481 * MDP(9);
t464 = t481 * t377;
t266 = -t283 * t417 + t419 * t284;
t262 = pkin(5) * t335 - pkin(8) * t324 + t266;
t265 = -pkin(8) * t323 + t267;
t462 = t426 * t262 - t265 * t423;
t277 = t419 * t303 - t310 * t417;
t286 = -t314 * t417 + t419 * t318;
t461 = t426 * t323 + t423 * t324;
t289 = t419 * t332 - t417 * t520;
t294 = t419 * t331 - t339 * t417;
t460 = qJDD(2) - t505;
t457 = -t362 * t484 - t308;
t455 = t262 * t423 + t265 * t426;
t454 = t266 * t419 + t267 * t417;
t453 = t266 * t375 + t286 * t370;
t452 = -t267 * t375 - t287 * t370;
t272 = pkin(5) * t364 - pkin(8) * t345 + t286;
t274 = -pkin(8) * t343 + t287;
t263 = t272 * t426 - t274 * t423;
t264 = t272 * t423 + t274 * t426;
t282 = pkin(5) * t375 + pkin(8) * t494 + t294;
t288 = pkin(8) * t495 + t295;
t451 = t282 * t426 - t288 * t423;
t450 = t282 * t423 + t288 * t426;
t448 = -t286 * t417 + t287 * t419;
t444 = -qJD(4) * t370 - qJDD(4) * t375;
t382 = t506 * t419;
t441 = pkin(5) * t366 + qJD(5) * t417 + qJD(6) * t382 + t364 * t510 + t289;
t381 = t506 * t417;
t440 = pkin(8) * t497 - qJD(5) * t419 + qJD(6) * t381 + t290;
t438 = t362 * t373;
t268 = -t345 * t476 + t471;
t436 = pkin(4) * t401 - qJ(5) * t403 + t404;
t269 = -qJD(6) * t447 + t461;
t311 = -qJD(3) * t374 + qJD(4) * t339;
t429 = qJD(1) ^ 2;
t421 = -pkin(7) - qJ(3);
t406 = t428 * qJ(2);
t394 = -pkin(5) * t419 - pkin(4);
t363 = t364 ^ 2;
t354 = t401 * t490 - t493;
t353 = t401 * t492 + t491;
t352 = t401 * t491 + t492;
t351 = -t401 * t493 + t490;
t330 = t373 * t374;
t329 = t376 * t374;
t312 = -pkin(5) * t495 + t338;
t297 = -pkin(5) * t497 + t322;
t296 = pkin(5) * t496 + t311;
t293 = pkin(5) * t343 + t309;
t292 = t369 * t376 - t475 * t494 + t476 * t495;
t291 = t368 * t374 - t369 * t373;
t273 = -pkin(8) * t496 + t278;
t271 = pkin(5) * t370 - t369 * t510 + t277;
t270 = pkin(5) * t323 + t285;
t1 = [(t277 * t364 + t294 * t335 + t311 * t343 + t338 * t323 + t417 * t432 - t419 * t442 + t453) * MDP(18) + (qJD(2) * t364 - qJD(4) * t311 - qJDD(4) * t338 + t335 * t391 + t370 * t380 + t372 * t375 - t442) * MDP(16) + t518 * MDP(2) + (qJDD(2) - t518 - 0.2e1 * t505) * MDP(4) + (t518 + t481 * (-t377 - t512)) * MDP(9) + (t268 * t330 - t291 * t447) * MDP(22) + (t268 * t375 + t291 * t362 + t330 * t333 - t370 * t447) * MDP(24) + (-(t271 * t423 + t273 * t426) * t362 - t450 * t333 - t455 * t375 - t264 * t370 - t296 * t447 + t312 * t268 + t270 * t330 + t293 * t291 + g(1) * t353 - g(2) * t351 + (-t263 * t375 - t362 * t451) * qJD(6)) * MDP(28) + (-t460 * pkin(1) - g(1) * (-pkin(1) * t425 + t406) - t467 + (t414 + t511) * qJ(2)) * MDP(6) + (0.2e1 * t414 + t511 - t459) * MDP(5) + t459 * MDP(3) + t444 * MDP(14) + qJDD(1) * MDP(1) + t516 * (t517 + t519) + (qJD(2) * t366 - qJD(4) * t310 - qJDD(4) * t339 + t334 * t391 + t369 * t380 - t372 * t374 - t522) * MDP(17) + (-t277 * t345 - t278 * t343 - t294 * t324 - t295 * t323 + t522 + t454 * t374 + (-t286 * t419 - t287 * t417) * t369) * MDP(20) + (-t278 * t364 - t295 * t335 + t311 * t345 + t338 * t324 + t417 * t442 + t419 * t432 + t452) * MDP(19) + (t383 * qJ(2) + t397 * qJD(2) - g(1) * (t422 * t425 + t406) - g(2) * (qJ(3) * t428 + t483) + t422 * t464 - qJD(3) * t523) * MDP(10) + (-t269 * t375 - t292 * t362 + t329 * t333 + t370 * t526) * MDP(25) + ((t271 * t426 - t273 * t423) * t362 + t451 * t333 + t462 * t375 + t263 * t370 - t296 * t526 + t312 * t269 - t270 * t329 + t293 * t292 - g(1) * t354 - g(2) * t352 + (-t264 * t375 - t362 * t450) * qJD(6)) * MDP(27) + (t268 * t329 - t269 * t330 + t291 * t526 + t292 * t447) * MDP(23) + (-t334 * t374 + t366 * t369) * MDP(11) + (t333 * t375 + t362 * t370) * MDP(26) + (-t334 * t375 + t335 * t374 - t364 * t369 - t366 * t370) * MDP(12) + t486 * MDP(13) + (t267 * t295 + t287 * t278 + t266 * t294 + t286 * t277 + t285 * t338 + t309 * t311 - g(1) * t406 - t467 + (-g(1) * t436 + g(2) * t421) * t428 + (-g(1) * (-pkin(1) + t421) - g(2) * t436) * t425) * MDP(21); (t460 - t518) * MDP(6) + (-qJD(1) * t397 + t464 - t518) * MDP(10) + (-qJD(1) * t364 + t486) * MDP(16) + (-qJD(1) * t366 + t444) * MDP(17) + (-t375 * t501 + t323 * t374 - t343 * t369 + (-qJD(1) * t419 - t370 * t417) * t364) * MDP(18) + (-t375 * t500 + t324 * t374 - t345 * t369 + (qJD(1) * t417 - t370 * t419) * t364) * MDP(19) + ((qJD(1) * t345 - t323 * t375 - t343 * t370) * t419 + (qJD(1) * t343 + t324 * t375 + t345 * t370) * t417) * MDP(20) + ((-qJD(1) * t286 - t452) * t419 + (-qJD(1) * t287 - t453) * t417 - t432) * MDP(21) + (t374 * t269 + t369 * t526 - t370 * t439 + qJD(1) * t438 + (t362 * t515 - t502) * t375) * MDP(27) + (t374 * t268 + t369 * t447 + t370 * t438 + qJD(1) * t439 + (qJD(6) * t439 + t308) * t375) * MDP(28) + (MDP(4) - t465) * qJDD(1) + (-qJ(2) * MDP(6) - MDP(5) - t516) * t429; (qJD(1) * t523 + t519) * MDP(10) - t384 * MDP(16) + t386 * MDP(17) + (-t343 * t366 + t500) * MDP(18) + (-t345 * t366 - t363 * t419 - t501) * MDP(19) + (-t323 * t417 - t324 * t419) * MDP(20) + (-t309 * t366 + t454 - t459) * MDP(21) + (t457 + t504) * MDP(27) + (t503 + t514) * MDP(28) - t429 * t465 + ((-t343 * t419 + t345 * t417) * MDP(20) + t448 * MDP(21) - MDP(18) * t497) * t364 + ((MDP(16) * t424 + MDP(8)) * t420 + (MDP(16) * t427 - MDP(17) * t424 + MDP(7)) * t418) * qJDD(1) + ((t366 + t470) * MDP(16) + (-t418 * t479 - t420 * t480 - t364) * MDP(17)) * qJD(4); t366 * t364 * MDP(11) + (t366 ^ 2 - t363) * MDP(12) + t456 * MDP(13) + ((t366 - t470) * qJD(4) + t513) * MDP(14) + qJDD(4) * MDP(15) + (qJD(4) * t322 - t366 * t380 - t435 + t443) * MDP(16) + (t364 * t380 - t433 - t446) * MDP(17) + (-qJ(5) * t501 - pkin(4) * t323 - t286 * t366 - t322 * t343 + (t417 * t474 - t289) * t364 - t521 * t419) * MDP(18) + (-qJ(5) * t500 - pkin(4) * t324 + t287 * t366 - t322 * t345 + (t419 * t474 + t290) * t364 + t521 * t417) * MDP(19) + (t289 * t345 + t290 * t343 + (-qJ(5) * t323 - qJD(5) * t343 - t286 * t364 + t267) * t419 + (qJ(5) * t324 + qJD(5) * t345 - t287 * t364 - t266) * t417 + t433) * MDP(20) + (-t286 * t289 - t287 * t290 - t309 * t322 + t448 * qJD(5) - t521 * pkin(4) + (-t266 * t417 + t267 * t419 + t433) * qJ(5)) * MDP(21) + (t268 * t376 - t447 * t485) * MDP(22) + (-t268 * t373 - t269 * t376 + t447 * t484 + t485 * t526) * MDP(23) + (t503 - t514) * MDP(24) + (t457 - t504) * MDP(25) - t362 * t366 * MDP(26) + ((-t381 * t426 - t382 * t423) * t333 + t394 * t269 + t270 * t373 - t263 * t366 + t297 * t526 + (t423 * t440 - t426 * t441) * t362 + t484 * t293 - t435 * t402) * MDP(27) + (-(-t381 * t423 + t382 * t426) * t333 + t394 * t268 + t270 * t376 + t264 * t366 + t297 * t447 + (t423 * t441 + t426 * t440) * t362 + t485 * t293 + t435 * t400) * MDP(28); (t345 * t364 + t323) * MDP(18) + (-t343 * t364 + t324) * MDP(19) + (-t343 ^ 2 - t345 ^ 2) * MDP(20) + (t286 * t345 + t287 * t343 + t521) * MDP(21) + (t269 - t525) * MDP(27) + (t268 + t527) * MDP(28); t447 * t526 * MDP(22) + (t447 ^ 2 - t526 ^ 2) * MDP(23) + (t471 - t527) * MDP(24) + (-t461 - t525) * MDP(25) + t333 * MDP(26) + (-g(1) * t351 - g(2) * t353 + t264 * t362 + t293 * t447 + t400 * t508 + t462) * MDP(27) + (g(1) * t352 - g(2) * t354 + t263 * t362 - t293 * t526 + t402 * t508 - t455) * MDP(28) + (-MDP(24) * t499 + MDP(25) * t447 - MDP(27) * t264 - MDP(28) * t263) * qJD(6);];
tau  = t1;

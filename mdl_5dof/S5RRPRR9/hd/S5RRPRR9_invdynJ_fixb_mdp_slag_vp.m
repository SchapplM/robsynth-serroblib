% Calculate vector of inverse dynamics joint torques for
% S5RRPRR9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRPRR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:22:10
% EndTime: 2019-12-31 20:22:17
% DurationCPUTime: 5.64s
% Computational Cost: add. (3724->430), mult. (8739->573), div. (0->0), fcn. (6603->14), ass. (0->202)
t427 = sin(pkin(9));
t432 = sin(qJ(2));
t486 = qJD(1) * t432;
t428 = cos(pkin(9));
t436 = cos(qJ(2));
t499 = t428 * t436;
t382 = qJD(1) * t499 - t427 * t486;
t527 = qJD(4) + qJD(5);
t536 = t382 - t527;
t429 = -qJ(3) - pkin(6);
t470 = qJD(2) * t429;
t381 = -qJD(3) * t432 + t436 * t470;
t404 = t429 * t432;
t342 = qJDD(2) * pkin(2) + qJD(1) * t381 + qJDD(1) * t404;
t380 = qJD(3) * t436 + t432 * t470;
t405 = t429 * t436;
t350 = qJD(1) * t380 - qJDD(1) * t405;
t296 = t342 * t428 - t427 * t350;
t294 = -qJDD(2) * pkin(3) - t296;
t373 = qJD(4) - t382;
t411 = pkin(2) * t427 + pkin(7);
t423 = qJ(2) + pkin(9);
t416 = sin(t423);
t417 = cos(t423);
t433 = sin(qJ(1));
t437 = cos(qJ(1));
t461 = g(1) * t437 + g(2) * t433;
t444 = -g(3) * t417 + t416 * t461;
t535 = -qJD(4) * t411 * t373 - t294 + t444;
t394 = t427 * t436 + t428 * t432;
t384 = t394 * qJD(1);
t431 = sin(qJ(4));
t435 = cos(qJ(4));
t481 = t435 * qJD(2);
t361 = t384 * t431 - t481;
t434 = cos(qJ(5));
t363 = qJD(2) * t431 + t384 * t435;
t430 = sin(qJ(5));
t510 = t363 * t430;
t305 = t434 * t361 + t510;
t371 = qJD(5) + t373;
t534 = t305 * t371;
t454 = t361 * t430 - t434 * t363;
t533 = t371 * t454;
t498 = t430 * t435;
t397 = t431 * t434 + t498;
t489 = t536 * t397;
t485 = qJD(4) * t431;
t509 = t382 * t431;
t532 = t485 - t509;
t415 = pkin(2) * t436 + pkin(1);
t402 = -qJD(1) * t415 + qJD(3);
t314 = -pkin(3) * t382 - pkin(7) * t384 + t402;
t398 = qJD(1) * t404;
t521 = qJD(2) * pkin(2);
t390 = t398 + t521;
t399 = qJD(1) * t405;
t500 = t428 * t399;
t347 = t427 * t390 - t500;
t336 = qJD(2) * pkin(7) + t347;
t289 = t314 * t431 + t336 * t435;
t283 = -pkin(8) * t361 + t289;
t483 = qJD(5) * t430;
t280 = t283 * t483;
t387 = t427 * t399;
t346 = t390 * t428 + t387;
t335 = -qJD(2) * pkin(3) - t346;
t303 = pkin(4) * t361 + t335;
t426 = qJ(4) + qJ(5);
t422 = cos(t426);
t502 = t422 * t433;
t421 = sin(t426);
t503 = t421 * t437;
t365 = -t417 * t502 + t503;
t501 = t422 * t437;
t504 = t421 * t433;
t367 = t417 * t501 + t504;
t525 = g(3) * t416;
t531 = g(1) * t367 - g(2) * t365 + t303 * t305 + t422 * t525 + t280;
t364 = t417 * t504 + t501;
t366 = -t417 * t503 + t502;
t478 = qJDD(1) * t436;
t479 = qJDD(1) * t432;
t455 = -t427 * t479 + t428 * t478;
t348 = -qJD(2) * t384 + t455;
t480 = qJD(1) * qJD(2);
t471 = t436 * t480;
t472 = t432 * t480;
t349 = qJDD(1) * t394 - t427 * t472 + t428 * t471;
t445 = pkin(2) * t472 - qJDD(1) * t415 + qJDD(3);
t290 = -pkin(3) * t348 - pkin(7) * t349 + t445;
t287 = t435 * t290;
t297 = t427 * t342 + t428 * t350;
t295 = qJDD(2) * pkin(7) + t297;
t301 = qJD(4) * t481 + t431 * qJDD(2) + t435 * t349 - t384 * t485;
t383 = t394 * qJD(2);
t343 = qJD(1) * t383 + qJDD(4) - t455;
t268 = pkin(4) * t343 - pkin(8) * t301 - qJD(4) * t289 - t295 * t431 + t287;
t302 = qJD(4) * t363 - t435 * qJDD(2) + t349 * t431;
t484 = qJD(4) * t435;
t448 = t431 * t290 + t435 * t295 + t314 * t484 - t336 * t485;
t269 = -pkin(8) * t302 + t448;
t468 = t434 * t268 - t430 * t269;
t530 = -g(1) * t366 + g(2) * t364 + t303 * t454 + t421 * t525 + t468;
t338 = qJDD(5) + t343;
t529 = t338 * MDP(24) + (-t305 ^ 2 + t454 ^ 2) * MDP(21) - t305 * MDP(20) * t454;
t331 = t397 * t394;
t396 = t430 * t431 - t434 * t435;
t490 = t536 * t396;
t526 = -t338 * t397 - t371 * t490;
t467 = t301 * t430 + t434 * t302;
t274 = -qJD(5) * t454 + t467;
t523 = g(3) * t436;
t522 = pkin(8) + t411;
t288 = t435 * t314 - t336 * t431;
t282 = -pkin(8) * t363 + t288;
t277 = pkin(4) * t373 + t282;
t520 = t277 * t434;
t519 = t283 * t434;
t518 = t301 * t431;
t517 = t305 * t384;
t516 = t454 * t384;
t514 = t361 * t373;
t513 = t361 * t384;
t512 = t363 * t373;
t511 = t363 * t384;
t393 = t427 * t432 - t499;
t386 = t393 * qJD(2);
t508 = t386 * t431;
t507 = t386 * t435;
t506 = t394 * t431;
t505 = t394 * t435;
t497 = t431 * t343;
t496 = t431 * t433;
t495 = t431 * t437;
t494 = t433 * t435;
t330 = t435 * t343;
t360 = t404 * t427 - t405 * t428;
t353 = t435 * t360;
t493 = t435 * t437;
t325 = pkin(2) * t486 + pkin(3) * t384 - pkin(7) * t382;
t352 = t398 * t428 + t387;
t492 = t431 * t325 + t435 * t352;
t345 = pkin(3) * t393 - pkin(7) * t394 - t415;
t491 = t431 * t345 + t353;
t424 = t432 ^ 2;
t488 = -t436 ^ 2 + t424;
t482 = qJD(5) * t434;
t477 = t432 * t521;
t476 = t434 * t301 - t430 * t302 - t361 * t482;
t412 = -pkin(2) * t428 - pkin(3);
t475 = t394 * t485;
t474 = t394 * t484;
t469 = qJD(4) * t522;
t323 = t380 * t427 - t428 * t381;
t351 = t398 * t427 - t500;
t359 = -t428 * t404 - t405 * t427;
t465 = t373 * t435;
t464 = -qJD(4) * t314 - t295;
t463 = qJD(5) * t277 + t269;
t462 = t532 * pkin(4) - t351;
t460 = g(1) * t433 - g(2) * t437;
t459 = -t396 * t338 + t489 * t371;
t458 = -t336 * t484 + t287;
t318 = t435 * t325;
t392 = t522 * t435;
t457 = pkin(4) * t384 + qJD(5) * t392 - t352 * t431 + t318 + (-pkin(8) * t382 + t469) * t435;
t391 = t522 * t431;
t456 = -pkin(8) * t509 + qJD(5) * t391 + t431 * t469 + t492;
t271 = t277 * t430 + t519;
t453 = -t532 * t373 + t330;
t452 = -0.2e1 * pkin(1) * t480 - pkin(6) * qJDD(2);
t451 = t474 - t508;
t450 = -t475 - t507;
t324 = t380 * t428 + t381 * t427;
t326 = pkin(3) * t383 + pkin(7) * t386 + t477;
t447 = t435 * t324 + t431 * t326 + t345 * t484 - t360 * t485;
t273 = -t363 * t483 + t476;
t446 = t335 * t373 - t411 * t343;
t438 = qJD(2) ^ 2;
t441 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t438 + t460;
t439 = qJD(1) ^ 2;
t440 = pkin(1) * t439 - pkin(6) * qJDD(1) + t461;
t403 = -pkin(4) * t435 + t412;
t377 = t417 * t493 + t496;
t376 = -t417 * t495 + t494;
t375 = -t417 * t494 + t495;
t374 = t417 * t496 + t493;
t334 = t435 * t345;
t332 = t396 * t394;
t322 = pkin(4) * t506 + t359;
t319 = t435 * t326;
t293 = pkin(4) * t451 + t323;
t291 = -pkin(8) * t506 + t491;
t284 = pkin(4) * t393 - pkin(8) * t505 - t360 * t431 + t334;
t279 = -t386 * t498 - t430 * t475 - t483 * t506 + (t527 * t505 - t508) * t434;
t278 = -t527 * t331 + t396 * t386;
t276 = pkin(4) * t302 + t294;
t275 = -pkin(8) * t451 + t447;
t272 = pkin(8) * t507 + pkin(4) * t383 - t324 * t431 + t319 + (-t353 + (pkin(8) * t394 - t345) * t431) * qJD(4);
t270 = -t283 * t430 + t520;
t1 = [(-g(1) * t364 - g(2) * t366 - t271 * t383 + t322 * t273 - t276 * t332 + t303 * t278 + t280 * t393 - t293 * t454 + (-(-qJD(5) * t291 + t272) * t371 - t284 * t338 - t268 * t393) * t430 + (-(qJD(5) * t284 + t275) * t371 - t291 * t338 - t463 * t393) * t434) * MDP(26) + (-t273 * t331 + t274 * t332 - t278 * t305 + t279 * t454) * MDP(21) + (-t273 * t332 - t278 * t454) * MDP(20) + (t273 * t393 + t278 * t371 - t332 * t338 - t383 * t454) * MDP(22) + (-t296 * t394 - t297 * t393 + t323 * t384 + t324 * t382 + t346 * t386 - t347 * t383 + t348 * t360 + t349 * t359 - t461) * MDP(11) + (-(-t361 * t435 - t363 * t431) * t386 + (-t518 - t302 * t435 + (t361 * t431 - t363 * t435) * qJD(4)) * t394) * MDP(14) + ((-t360 * t484 + t319) * t373 + t334 * t343 + t458 * t393 + t288 * t383 + t323 * t361 + t359 * t302 + t335 * t474 - g(1) * t375 - g(2) * t377 + ((-qJD(4) * t345 - t324) * t373 - t360 * t343 + t464 * t393 + t294 * t394 - t335 * t386) * t431) * MDP(18) + (t301 * t505 + t363 * t450) * MDP(13) + (-g(1) * t374 - g(2) * t376 - t289 * t383 + t294 * t505 + t359 * t301 + t323 * t363 + t335 * t450 - t343 * t491 - t373 * t447 - t393 * t448) * MDP(19) + (-t302 * t393 - t361 * t383 - t373 * t451 - t394 * t497) * MDP(16) + (t301 * t393 + t330 * t394 + t363 * t383 + t373 * t450) * MDP(15) + 0.2e1 * (t432 * t478 - t480 * t488) * MDP(5) + (t297 * t360 + t347 * t324 - t296 * t359 - t346 * t323 - t445 * t415 + t402 * t477 - g(1) * (-t415 * t433 - t429 * t437) - g(2) * (t415 * t437 - t429 * t433)) * MDP(12) + (qJDD(1) * t424 + 0.2e1 * t432 * t471) * MDP(4) + ((t272 * t434 - t275 * t430) * t371 + (t284 * t434 - t291 * t430) * t338 + t468 * t393 + t270 * t383 + t293 * t305 + t322 * t274 + t276 * t331 + t303 * t279 - g(1) * t365 - g(2) * t367 + ((-t284 * t430 - t291 * t434) * t371 - t271 * t393) * qJD(5)) * MDP(25) + t460 * MDP(2) + t461 * MDP(3) + (t432 * t452 + t436 * t441) * MDP(9) + (-t432 * t441 + t436 * t452) * MDP(10) + (qJDD(2) * t432 + t436 * t438) * MDP(6) + (qJDD(2) * t436 - t432 * t438) * MDP(7) + (t343 * t393 + t373 * t383) * MDP(17) + (-t274 * t393 - t279 * t371 - t305 * t383 - t331 * t338) * MDP(23) + (t338 * t393 + t371 * t383) * MDP(24) + qJDD(1) * MDP(1); MDP(6) * t479 + MDP(7) * t478 + qJDD(2) * MDP(8) + (t432 * t440 - t523) * MDP(9) + (g(3) * t432 + t436 * t440) * MDP(10) + ((t346 - t352) * t382 + (t348 * t427 - t349 * t428) * pkin(2)) * MDP(11) + (t346 * t351 - t347 * t352 + (-t523 + t296 * t428 + t297 * t427 + (-qJD(1) * t402 + t461) * t432) * pkin(2)) * MDP(12) + (t363 * t465 + t518) * MDP(13) + ((t301 - t514) * t435 + (-t302 - t512) * t431) * MDP(14) + (t373 * t465 + t497 - t511) * MDP(15) + (t453 + t513) * MDP(16) + (t412 * t302 - t318 * t373 - t351 * t361 + (t352 * t373 + t446) * t431 + t535 * t435) * MDP(18) + (t412 * t301 - t351 * t363 + t492 * t373 - t535 * t431 + t446 * t435) * MDP(19) + (t273 * t397 - t454 * t490) * MDP(20) + (-t273 * t396 - t274 * t397 - t305 * t490 - t454 * t489) * MDP(21) + (t516 - t526) * MDP(22) + (t459 + t517) * MDP(23) + ((-t391 * t434 - t392 * t430) * t338 + t403 * t274 + t276 * t396 + (t430 * t456 - t434 * t457) * t371 + t462 * t305 - t489 * t303 + t444 * t422) * MDP(25) + (-(-t391 * t430 + t392 * t434) * t338 + t403 * t273 + t276 * t397 + (t430 * t457 + t434 * t456) * t371 - t462 * t454 + t490 * t303 - t444 * t421) * MDP(26) + (-t432 * t436 * MDP(4) + t488 * MDP(5)) * t439 - ((-t347 + t351) * MDP(11) + t373 * MDP(17) + t288 * MDP(18) - t289 * MDP(19) + t371 * MDP(24) + t270 * MDP(25) - t271 * MDP(26)) * t384; (-t382 ^ 2 - t384 ^ 2) * MDP(11) + (t346 * t384 - t347 * t382 + t445 - t460) * MDP(12) + (t453 - t513) * MDP(18) + (-t373 ^ 2 * t435 - t497 - t511) * MDP(19) + (t459 - t517) * MDP(25) + (t516 + t526) * MDP(26); t363 * t361 * MDP(13) + (-t361 ^ 2 + t363 ^ 2) * MDP(14) + (t301 + t514) * MDP(15) + (-t302 + t512) * MDP(16) + t343 * MDP(17) + (-g(1) * t376 + g(2) * t374 + t289 * t373 - t335 * t363 + (t464 + t525) * t431 + t458) * MDP(18) + (g(1) * t377 - g(2) * t375 + t288 * t373 + t335 * t361 + t435 * t525 - t448) * MDP(19) + (t273 + t534) * MDP(22) + (-t274 - t533) * MDP(23) + (-(-t282 * t430 - t519) * t371 - t271 * qJD(5) + (-t305 * t363 + t434 * t338 - t371 * t483) * pkin(4) + t530) * MDP(25) + ((-t283 * t371 - t268) * t430 + (t282 * t371 - t463) * t434 + (-t430 * t338 + t363 * t454 - t371 * t482) * pkin(4) + t531) * MDP(26) + t529; (t476 + t534) * MDP(22) + (-t467 - t533) * MDP(23) + (t271 * t371 + t530) * MDP(25) + (-t430 * t268 - t434 * t269 + t270 * t371 + t531) * MDP(26) + (-MDP(22) * t510 + MDP(23) * t454 - MDP(25) * t271 - MDP(26) * t520) * qJD(5) + t529;];
tau = t1;

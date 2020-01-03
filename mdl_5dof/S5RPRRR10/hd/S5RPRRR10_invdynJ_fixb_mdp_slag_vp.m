% Calculate vector of inverse dynamics joint torques for
% S5RPRRR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:55
% EndTime: 2019-12-31 19:11:04
% DurationCPUTime: 5.85s
% Computational Cost: add. (3522->420), mult. (8363->545), div. (0->0), fcn. (6564->14), ass. (0->193)
t414 = cos(pkin(9));
t421 = cos(qJ(3));
t473 = qJD(1) * t421;
t395 = t414 * t473;
t413 = sin(pkin(9));
t417 = sin(qJ(3));
t474 = qJD(1) * t417;
t460 = t413 * t474;
t371 = t395 - t460;
t517 = qJD(4) + qJD(5);
t531 = t371 - t517;
t467 = qJD(1) * qJD(2);
t510 = pkin(6) + qJ(2);
t514 = t510 * qJDD(1) + t467;
t358 = t514 * t413;
t359 = t514 * t414;
t387 = t510 * t413;
t377 = qJD(1) * t387;
t388 = t510 * t414;
t378 = qJD(1) * t388;
t339 = -t417 * t377 + t421 * t378;
t529 = qJD(3) * t339;
t437 = -t358 * t421 - t417 * t359 - t529;
t289 = -qJDD(3) * pkin(3) - t437;
t363 = qJD(4) - t371;
t411 = pkin(9) + qJ(3);
t402 = sin(t411);
t403 = cos(t411);
t418 = sin(qJ(1));
t422 = cos(qJ(1));
t449 = g(1) * t422 + g(2) * t418;
t428 = -g(3) * t403 + t449 * t402;
t530 = -qJD(4) * pkin(7) * t363 - t289 + t428;
t376 = t413 * t421 + t414 * t417;
t372 = t376 * qJD(1);
t416 = sin(qJ(4));
t420 = cos(qJ(4));
t468 = t420 * qJD(3);
t348 = t372 * t416 - t468;
t419 = cos(qJ(5));
t350 = qJD(3) * t416 + t372 * t420;
t415 = sin(qJ(5));
t498 = t350 * t415;
t299 = t419 * t348 + t498;
t360 = qJD(5) + t363;
t528 = t299 * t360;
t442 = t348 * t415 - t419 * t350;
t527 = t360 * t442;
t486 = t415 * t420;
t380 = t416 * t419 + t486;
t477 = t531 * t380;
t374 = t376 * qJD(3);
t465 = qJDD(1) * t421;
t394 = t414 * t465;
t466 = qJDD(1) * t417;
t526 = -qJD(1) * t374 - t413 * t466 + t394;
t472 = qJD(4) * t416;
t497 = t371 * t416;
t525 = t472 - t497;
t397 = -pkin(2) * t414 - pkin(1);
t384 = t397 * qJD(1) + qJD(2);
t309 = -pkin(3) * t371 - pkin(7) * t372 + t384;
t332 = qJD(3) * pkin(7) + t339;
t283 = t309 * t416 + t332 * t420;
t278 = -pkin(8) * t348 + t283;
t470 = qJD(5) * t415;
t276 = t278 * t470;
t520 = -t377 * t421 - t417 * t378;
t331 = -qJD(3) * pkin(3) - t520;
t297 = pkin(4) * t348 + t331;
t412 = qJ(4) + qJ(5);
t408 = cos(t412);
t489 = t408 * t418;
t407 = sin(t412);
t490 = t407 * t422;
t353 = -t403 * t489 + t490;
t488 = t408 * t422;
t491 = t407 * t418;
t355 = t403 * t488 + t491;
t512 = g(3) * t402;
t524 = g(1) * t355 - g(2) * t353 + t297 * t299 + t408 * t512 + t276;
t352 = t403 * t491 + t488;
t354 = -t403 * t490 + t489;
t462 = qJD(3) * t395 + t413 * t465 + t414 * t466;
t336 = -qJD(3) * t460 + t462;
t383 = t397 * qJDD(1) + qJDD(2);
t291 = -pkin(3) * t526 - pkin(7) * t336 + t383;
t287 = t420 * t291;
t441 = -t358 * t417 + t359 * t421;
t288 = qJDD(3) * pkin(7) + qJD(3) * t520 + t441;
t295 = qJD(4) * t468 + t416 * qJDD(3) + t420 * t336 - t372 * t472;
t330 = qJDD(4) - t526;
t264 = pkin(4) * t330 - pkin(8) * t295 - t283 * qJD(4) - t288 * t416 + t287;
t296 = t350 * qJD(4) - t420 * qJDD(3) + t336 * t416;
t471 = qJD(4) * t420;
t431 = t420 * t288 + t416 * t291 + t309 * t471 - t332 * t472;
t265 = -pkin(8) * t296 + t431;
t457 = t419 * t264 - t415 * t265;
t523 = -g(1) * t354 + g(2) * t352 + t297 * t442 + t407 * t512 + t457;
t327 = qJDD(5) + t330;
t522 = t327 * MDP(26) + (-t299 ^ 2 + t442 ^ 2) * MDP(23) - t299 * MDP(22) * t442;
t322 = t380 * t376;
t518 = qJ(2) * qJDD(1);
t448 = g(1) * t418 - g(2) * t422;
t516 = qJDD(2) - t448;
t379 = t415 * t416 - t419 * t420;
t478 = t531 * t379;
t515 = -t327 * t380 - t478 * t360;
t456 = t295 * t415 + t419 * t296;
t270 = -t442 * qJD(5) + t456;
t513 = pkin(7) + pkin(8);
t509 = qJDD(1) * pkin(1);
t282 = t420 * t309 - t332 * t416;
t277 = -pkin(8) * t350 + t282;
t273 = pkin(4) * t363 + t277;
t508 = t273 * t419;
t507 = t278 * t419;
t506 = t295 * t416;
t505 = t299 * t372;
t504 = t442 * t372;
t502 = t348 * t363;
t501 = t348 * t372;
t500 = t350 * t363;
t499 = t350 * t372;
t375 = t413 * t417 - t421 * t414;
t373 = t375 * qJD(3);
t496 = t373 * t416;
t495 = t373 * t420;
t494 = t376 * t416;
t493 = t376 * t420;
t487 = t414 * MDP(4);
t485 = t416 * t330;
t484 = t416 * t418;
t483 = t416 * t422;
t482 = t418 * t420;
t317 = t420 * t330;
t347 = -t387 * t417 + t388 * t421;
t340 = t420 * t347;
t481 = t420 * t422;
t333 = pkin(3) * t372 - pkin(7) * t371;
t480 = t416 * t333 + t420 * t520;
t335 = pkin(3) * t375 - pkin(7) * t376 + t397;
t479 = t416 * t335 + t340;
t476 = t413 ^ 2 + t414 ^ 2;
t469 = qJD(5) * t419;
t463 = t419 * t295 - t415 * t296 - t348 * t469;
t461 = qJD(4) * t513;
t459 = t376 * t472;
t458 = t376 * t471;
t454 = t363 * t420;
t453 = -qJD(4) * t309 - t288;
t452 = qJD(5) * t273 + t265;
t451 = 0.2e1 * t476;
t450 = t525 * pkin(4) - t339;
t447 = -t379 * t327 + t477 * t360;
t446 = -t332 * t471 + t287;
t320 = t420 * t333;
t390 = t513 * t420;
t445 = pkin(4) * t372 + qJD(5) * t390 - t520 * t416 + t320 + (-pkin(8) * t371 + t461) * t420;
t389 = t513 * t416;
t444 = -pkin(8) * t497 + qJD(5) * t389 + t416 * t461 + t480;
t267 = t273 * t415 + t507;
t439 = -t387 * t421 - t388 * t417;
t438 = t509 - t516;
t436 = -t525 * t363 + t317;
t435 = t458 - t496;
t434 = -t459 - t495;
t433 = -pkin(7) * t330 + t363 * t331;
t310 = -t375 * qJD(2) + t439 * qJD(3);
t334 = pkin(3) * t374 + pkin(7) * t373;
t430 = t420 * t310 + t416 * t334 + t335 * t471 - t347 * t472;
t269 = -t350 * t470 + t463;
t425 = t451 * t467 - t449;
t311 = t376 * qJD(2) + t347 * qJD(3);
t400 = -pkin(4) * t420 - pkin(3);
t367 = t403 * t481 + t484;
t366 = -t403 * t483 + t482;
t365 = -t403 * t482 + t483;
t364 = t403 * t484 + t481;
t325 = t420 * t335;
t323 = t379 * t376;
t321 = t420 * t334;
t312 = pkin(4) * t494 - t439;
t290 = t435 * pkin(4) + t311;
t285 = -pkin(8) * t494 + t479;
t280 = pkin(4) * t375 - pkin(8) * t493 - t347 * t416 + t325;
t275 = -t373 * t486 - t415 * t459 - t470 * t494 + (t517 * t493 - t496) * t419;
t274 = -t517 * t322 + t379 * t373;
t272 = pkin(4) * t296 + t289;
t271 = -t435 * pkin(8) + t430;
t268 = pkin(8) * t495 + pkin(4) * t374 - t310 * t416 + t321 + (-t340 + (pkin(8) * t376 - t335) * t416) * qJD(4);
t266 = -t278 * t415 + t508;
t1 = [(t295 * t493 + t434 * t350) * MDP(15) + (-t296 * t375 - t348 * t374 - t435 * t363 - t376 * t485) * MDP(18) + (t295 * t375 + t376 * t317 + t350 * t374 + t434 * t363) * MDP(17) + ((t268 * t419 - t271 * t415) * t360 + (t280 * t419 - t285 * t415) * t327 + t457 * t375 + t266 * t374 + t290 * t299 + t312 * t270 + t272 * t322 + t297 * t275 - g(1) * t353 - g(2) * t355 + ((-t280 * t415 - t285 * t419) * t360 - t267 * t375) * qJD(5)) * MDP(27) + (-g(1) * t352 - g(2) * t354 - t267 * t374 + t312 * t269 - t272 * t323 + t297 * t274 + t276 * t375 - t290 * t442 + (-(-qJD(5) * t285 + t268) * t360 - t280 * t327 - t264 * t375) * t415 + (-(qJD(5) * t280 + t271) * t360 - t285 * t327 - t452 * t375) * t419) * MDP(28) + (t269 * t375 + t274 * t360 - t323 * t327 - t374 * t442) * MDP(24) + (-t269 * t322 + t270 * t323 - t274 * t299 + t275 * t442) * MDP(23) + (-t269 * t323 - t274 * t442) * MDP(22) + (-g(1) * t364 - g(2) * t366 - t283 * t374 + t289 * t493 - t295 * t439 + t311 * t350 - t479 * t330 + t434 * t331 - t430 * t363 - t431 * t375) * MDP(21) + ((-t347 * t471 + t321) * t363 + t325 * t330 + t446 * t375 + t282 * t374 + t311 * t348 - t439 * t296 + t331 * t458 - g(1) * t365 - g(2) * t367 + ((-qJD(4) * t335 - t310) * t363 - t347 * t330 + t453 * t375 + t289 * t376 - t331 * t373) * t416) * MDP(20) + t448 * MDP(2) + t449 * MDP(3) + (t330 * t375 + t363 * t374) * MDP(19) + (-qJD(3) * t374 - qJDD(3) * t375) * MDP(11) + (t327 * t375 + t360 * t374) * MDP(26) + (-t270 * t375 - t275 * t360 - t299 * t374 - t322 * t327) * MDP(25) + qJDD(1) * MDP(1) + (-qJD(3) * t311 + qJDD(3) * t439 + t374 * t384 + t375 * t383 - t397 * t526 + t448 * t403) * MDP(13) + (-t336 * t375 - t371 * t373 - t372 * t374 + t376 * t526) * MDP(9) + (-MDP(5) * t413 + t487) * (t438 + t509) + (t438 * pkin(1) + (t476 * t518 + t425) * qJ(2)) * MDP(7) + (t451 * t518 + t425) * MDP(6) + (-qJD(3) * t310 - qJDD(3) * t347 + t336 * t397 - t373 * t384 + t376 * t383 - t448 * t402) * MDP(14) + (-qJD(3) * t373 + qJDD(3) * t376) * MDP(10) + (t336 * t376 - t372 * t373) * MDP(8) + (-(-t348 * t420 - t350 * t416) * t373 + (-t506 - t296 * t420 + (t348 * t416 - t350 * t420) * qJD(4)) * t376) * MDP(16); t516 * MDP(7) - t394 * MDP(13) + t462 * MDP(14) + (t436 - t501) * MDP(20) + (-t363 ^ 2 * t420 - t485 - t499) * MDP(21) + (t447 - t505) * MDP(27) + (t504 + t515) * MDP(28) + (-t487 - pkin(1) * MDP(7) + (MDP(13) * t417 + MDP(5)) * t413) * qJDD(1) + ((t413 * t473 + t414 * t474 + t372) * MDP(13) + (t371 - t460) * MDP(14)) * qJD(3) + (-MDP(7) * qJ(2) - MDP(6)) * qJD(1) ^ 2 * t476; -t371 ^ 2 * MDP(9) + ((-t371 - t460) * qJD(3) + t462) * MDP(10) + t526 * MDP(11) + qJDD(3) * MDP(12) + (t428 + t437 + t529) * MDP(13) + (-t371 * t384 + t449 * t403 - t441 + t512) * MDP(14) + (t350 * t454 + t506) * MDP(15) + ((t295 - t502) * t420 + (-t296 - t500) * t416) * MDP(16) + (t363 * t454 + t485 - t499) * MDP(17) + (t436 + t501) * MDP(18) + (-pkin(3) * t296 - t320 * t363 - t339 * t348 + (t363 * t520 + t433) * t416 + t530 * t420) * MDP(20) + (-pkin(3) * t295 - t339 * t350 + t480 * t363 - t530 * t416 + t433 * t420) * MDP(21) + (t269 * t380 - t442 * t478) * MDP(22) + (-t269 * t379 - t270 * t380 - t478 * t299 - t442 * t477) * MDP(23) + (t504 - t515) * MDP(24) + (t447 + t505) * MDP(25) + ((-t389 * t419 - t390 * t415) * t327 + t400 * t270 + t272 * t379 + (t444 * t415 - t445 * t419) * t360 + t450 * t299 - t477 * t297 + t428 * t408) * MDP(27) + (-(-t389 * t415 + t390 * t419) * t327 + t400 * t269 + t272 * t380 + (t445 * t415 + t444 * t419) * t360 - t450 * t442 + t478 * t297 - t428 * t407) * MDP(28) + (qJD(3) * MDP(11) - t384 * MDP(13) - t363 * MDP(19) - t282 * MDP(20) + t283 * MDP(21) - t360 * MDP(26) - t266 * MDP(27) + t267 * MDP(28) - t371 * MDP(8) + MDP(9) * t372) * t372; t350 * t348 * MDP(15) + (-t348 ^ 2 + t350 ^ 2) * MDP(16) + (t295 + t502) * MDP(17) + (-t296 + t500) * MDP(18) + t330 * MDP(19) + (-g(1) * t366 + g(2) * t364 + t283 * t363 - t331 * t350 + (t453 + t512) * t416 + t446) * MDP(20) + (g(1) * t367 - g(2) * t365 + t282 * t363 + t331 * t348 + t420 * t512 - t431) * MDP(21) + (t269 + t528) * MDP(24) + (-t270 - t527) * MDP(25) + (-(-t277 * t415 - t507) * t360 - t267 * qJD(5) + (-t299 * t350 + t327 * t419 - t360 * t470) * pkin(4) + t523) * MDP(27) + ((-t278 * t360 - t264) * t415 + (t277 * t360 - t452) * t419 + (-t327 * t415 + t350 * t442 - t360 * t469) * pkin(4) + t524) * MDP(28) + t522; (t463 + t528) * MDP(24) + (-t456 - t527) * MDP(25) + (t267 * t360 + t523) * MDP(27) + (-t415 * t264 - t419 * t265 + t266 * t360 + t524) * MDP(28) + (-MDP(24) * t498 + t442 * MDP(25) - t267 * MDP(27) - MDP(28) * t508) * qJD(5) + t522;];
tau = t1;

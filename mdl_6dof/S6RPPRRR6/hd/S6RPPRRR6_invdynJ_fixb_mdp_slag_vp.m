% Calculate vector of inverse dynamics joint torques for
% S6RPPRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:31:41
% EndTime: 2019-03-09 02:31:47
% DurationCPUTime: 5.04s
% Computational Cost: add. (2445->466), mult. (4524->606), div. (0->0), fcn. (2956->10), ass. (0->199)
t410 = (qJD(1) * qJD(2));
t411 = (qJ(2) * qJDD(1));
t466 = qJDD(3) + t410 + t411;
t363 = -pkin(7) * qJDD(1) + t466;
t395 = qJ(2) * qJD(1) + qJD(3);
t382 = -pkin(7) * qJD(1) + t395;
t425 = cos(qJ(4));
t421 = sin(qJ(4));
t500 = qJD(4) * t421;
t316 = -qJDD(4) * pkin(4) - t363 * t425 + t382 * t500;
t505 = qJD(1) * t421;
t391 = qJD(5) + t505;
t422 = sin(qJ(1));
t426 = cos(qJ(1));
t453 = g(1) * t426 + g(2) * t422;
t539 = g(3) * t421;
t432 = t425 * t453 - t539;
t550 = qJD(5) * pkin(8) * t391 + t316 + t432;
t420 = sin(qJ(5));
t424 = cos(qJ(5));
t490 = t424 * qJD(4);
t504 = qJD(1) * t425;
t358 = t420 * t504 - t490;
t423 = cos(qJ(6));
t491 = t420 * qJD(4);
t360 = t424 * t504 + t491;
t419 = sin(qJ(6));
t529 = t360 * t419;
t312 = t423 * t358 + t529;
t384 = qJD(6) + t391;
t549 = t312 * t384;
t449 = t358 * t419 - t423 * t360;
t548 = t384 * t449;
t418 = pkin(1) + qJ(3);
t454 = pkin(4) * t421 - pkin(8) * t425;
t369 = t454 + t418;
t328 = qJD(1) * t369 - qJD(2);
t368 = t421 * t382;
t342 = qJD(4) * pkin(8) + t368;
t299 = t328 * t420 + t342 * t424;
t296 = -pkin(9) * t358 + t299;
t494 = qJD(6) * t419;
t292 = t296 * t494;
t527 = t382 * t425;
t343 = -qJD(4) * pkin(4) - t527;
t321 = pkin(5) * t358 + t343;
t415 = qJ(5) + qJ(6);
t400 = sin(t415);
t401 = cos(t415);
t523 = t421 * t422;
t330 = -t400 * t426 - t401 * t523;
t515 = t426 * t401;
t520 = t422 * t400;
t332 = t421 * t515 - t520;
t538 = g(3) * t425;
t547 = g(1) * t332 - g(2) * t330 + t312 * t321 + t401 * t538 + t292;
t329 = t421 * t520 - t515;
t521 = t421 * t426;
t331 = -t400 * t521 - t401 * t422;
t455 = pkin(4) * t425 + pkin(8) * t421;
t356 = qJD(4) * t455 + qJD(3);
t409 = qJDD(1) * qJ(3);
t416 = qJDD(1) * pkin(1);
t481 = qJDD(2) - t416;
t465 = t409 - t481;
t306 = qJD(1) * t356 + qJDD(1) * t454 + t465;
t301 = t424 * t306;
t495 = qJD(5) * t425;
t436 = -t420 * t495 - t421 * t490;
t485 = qJDD(1) * t425;
t307 = qJD(1) * t436 + qJD(5) * t490 + t420 * qJDD(4) + t424 * t485;
t499 = qJD(4) * t425;
t317 = qJDD(4) * pkin(8) + t363 * t421 + t382 * t499;
t489 = qJD(1) * qJD(4);
t468 = t425 * t489;
t486 = qJDD(1) * t421;
t355 = qJDD(5) + t468 + t486;
t282 = pkin(5) * t355 - pkin(9) * t307 - qJD(5) * t299 - t317 * t420 + t301;
t474 = t421 * t491;
t308 = -qJD(1) * t474 + qJD(5) * t360 - t424 * qJDD(4) + t420 * t485;
t496 = qJD(5) * t424;
t479 = -t420 * t306 - t424 * t317 - t328 * t496;
t498 = qJD(5) * t420;
t440 = -t342 * t498 - t479;
t283 = -pkin(9) * t308 + t440;
t464 = t423 * t282 - t419 * t283;
t546 = -g(1) * t331 + g(2) * t329 + t321 * t449 + t400 * t538 + t464;
t350 = qJDD(6) + t355;
t545 = t350 * MDP(28) + (-t312 ^ 2 + t449 ^ 2) * MDP(25) - t312 * MDP(24) * t449;
t362 = t419 * t424 + t420 * t423;
t333 = t362 * t425;
t509 = g(1) * t422 - g(2) * t426;
t543 = -qJD(6) * t424 - t496;
t542 = qJD(1) * t418;
t482 = qJD(5) + qJD(6);
t383 = -qJD(2) + t542;
t417 = -pkin(7) + qJ(2);
t541 = qJD(4) * (qJD(2) + t383 + t542) + qJDD(4) * t417;
t463 = t307 * t419 + t423 * t308;
t287 = -qJD(6) * t449 + t463;
t404 = 2 * t410;
t540 = pkin(8) + pkin(9);
t298 = t424 * t328 - t342 * t420;
t295 = -pkin(9) * t360 + t298;
t291 = pkin(5) * t391 + t295;
t537 = t291 * t423;
t536 = t296 * t423;
t535 = t307 * t420;
t518 = t423 * t424;
t526 = t419 * t420;
t361 = -t518 + t526;
t534 = t350 * t361;
t533 = t350 * t362;
t532 = t355 * t424;
t531 = t358 * t391;
t530 = t360 * t391;
t528 = t360 * t424;
t525 = t420 * t355;
t524 = t420 * t425;
t522 = t421 * t424;
t519 = t422 * t424;
t517 = t424 * t425;
t516 = t424 * t426;
t437 = t543 * t423;
t476 = t420 * t505;
t514 = t419 * t476 + t482 * t526 - t505 * t518 + t437;
t439 = qJD(1) * t362;
t513 = t482 * t362 + t421 * t439;
t365 = t455 * qJD(1);
t512 = t420 * t365 + t382 * t517;
t511 = t420 * t369 + t417 * t522;
t510 = t426 * pkin(1) + t422 * qJ(2);
t414 = t425 ^ 2;
t508 = t421 ^ 2 - t414;
t427 = qJD(4) ^ 2;
t428 = qJD(1) ^ 2;
t507 = -t427 - t428;
t503 = qJD(2) * t421;
t502 = qJD(4) * t358;
t501 = qJD(4) * t360;
t497 = qJD(5) * t421;
t493 = qJD(6) * t423;
t488 = qJD(3) * qJD(1);
t487 = qJDD(1) * t418;
t483 = qJDD(4) * t421;
t478 = t423 * t307 - t419 * t308 - t358 * t493;
t477 = qJD(5) * t540;
t475 = t417 * t499;
t472 = t425 * t490;
t471 = t417 * t497;
t469 = t420 * t494;
t467 = qJDD(2) - t509;
t462 = t391 * t417 + t342;
t461 = -qJD(5) * t328 - t317;
t460 = qJD(6) * t291 + t283;
t459 = qJD(1) + t497;
t458 = t420 * t356 + t369 * t496 + t417 * t472 + t424 * t503;
t457 = -t416 + t467;
t456 = -t368 + (t476 + t498) * pkin(5);
t349 = t424 * t365;
t377 = t540 * t424;
t447 = pkin(5) * t425 + pkin(9) * t522;
t452 = qJD(1) * t447 + qJD(6) * t377 - t382 * t524 + t424 * t477 + t349;
t376 = t540 * t420;
t451 = pkin(9) * t476 + qJD(6) * t376 + t420 * t477 + t512;
t285 = t291 * t419 + t536;
t448 = -t409 + t457;
t445 = t391 * t496 + t525;
t444 = t391 * t498 - t532;
t442 = t404 + (2 * t411) - t453;
t441 = t384 * t361;
t286 = -t360 * t494 + t478;
t438 = qJD(1) * t383 + t453;
t435 = -t424 * t495 + t474;
t434 = -pkin(8) * t355 + t343 * t391;
t433 = -t363 + t438;
t364 = t465 + t488;
t429 = -t417 * t427 + t364 + t487 + t488 + t509;
t403 = t426 * qJ(2);
t398 = qJDD(4) * t425;
t394 = -pkin(5) * t424 - pkin(4);
t354 = (pkin(5) * t420 - t417) * t425;
t352 = t424 * t369;
t347 = -t420 * t422 + t421 * t516;
t346 = -t420 * t521 - t519;
t345 = -t420 * t426 - t421 * t519;
t344 = t420 * t523 - t516;
t338 = t424 * t356;
t334 = t361 * t425;
t322 = -pkin(5) * t435 - qJD(2) * t425 + t417 * t500;
t318 = -pkin(9) * t524 + t511;
t310 = -pkin(9) * t517 + t352 + (-t417 * t420 + pkin(5)) * t421;
t294 = -t425 * t469 + (t482 * t517 - t474) * t423 + t436 * t419;
t293 = -t482 * t333 + t361 * t500;
t290 = pkin(5) * t308 + t316;
t289 = pkin(9) * t435 - t420 * t471 + t458;
t288 = -t424 * t471 + t338 + t447 * qJD(4) + (-t475 - t503 + (pkin(9) * t425 - t369) * qJD(5)) * t420;
t284 = -t296 * t419 + t537;
t1 = [t453 * MDP(3) + (qJDD(3) + t442) * MDP(7) + t442 * MDP(5) + (-t286 * t333 + t287 * t334 - t293 * t312 + t294 * t449) * MDP(25) + (-t286 * t334 - t293 * t449) * MDP(24) + (-t285 * t499 - g(1) * t329 - g(2) * t331 + t354 * t286 - t290 * t334 + t292 * t421 + t321 * t293 - t322 * t449 + (-(-qJD(6) * t318 + t288) * t384 - t310 * t350 - t282 * t421) * t419 + (-(qJD(6) * t310 + t289) * t384 - t318 * t350 - t460 * t421) * t423) * MDP(30) + (t286 * t421 + t293 * t384 - t334 * t350 - t449 * t499) * MDP(26) + qJDD(1) * MDP(1) + ((-t391 * t490 + t307) * t421 + (-t444 + t501) * t425) * MDP(19) + (-g(1) * t345 - g(2) * t347 + t338 * t391 + t352 * t355 + (t417 * t502 - t462 * t496 + t301) * t421 + (-qJD(2) * t358 + qJD(4) * t298 - t308 * t417 + t343 * t496) * t425 + ((-qJD(5) * t369 - t475) * t391 + t316 * t425 + (-qJD(2) * t391 - t343 * qJD(4) - t417 * t355 + t461) * t421) * t420) * MDP(22) + ((t391 * t491 - t308) * t421 + (-t445 - t502) * t425) * MDP(20) + ((t358 * t424 + t360 * t420) * t500 + (-t535 - t308 * t424 + (t358 * t420 - t528) * qJD(5)) * t425) * MDP(18) + (t364 * t418 + t383 * qJD(3) + t466 * qJ(2) + t395 * qJD(2) - g(1) * (-t418 * t422 + t403) - g(2) * (qJ(3) * t426 + t510)) * MDP(9) + (-t448 + t487 + 0.2e1 * t488) * MDP(8) + (-t421 * t427 + t398) * MDP(12) + 0.2e1 * (-t421 * t485 + t508 * t489) * MDP(11) + t509 * MDP(2) + (-t458 * t391 - t511 * t355 - g(1) * t344 - g(2) * t346 + (t462 * t498 + (-t343 * t424 + t360 * t417) * qJD(4) + t479) * t421 + (-qJD(2) * t360 - qJD(4) * t299 - t307 * t417 + t316 * t424 - t343 * t498) * t425) * MDP(23) + (-t425 * t427 - t483) * MDP(13) + (qJDD(1) * t414 - 0.2e1 * t421 * t468) * MDP(10) + (-0.2e1 * t416 + t467) * MDP(4) + (t429 * t421 + t541 * t425) * MDP(15) + (-t541 * t421 + t429 * t425) * MDP(16) + (t307 * t517 + t360 * t436) * MDP(17) + (-t481 * pkin(1) - g(1) * (-pkin(1) * t422 + t403) - g(2) * t510 + (t404 + t411) * qJ(2)) * MDP(6) + (-t287 * t421 - t294 * t384 - t312 * t499 - t333 * t350) * MDP(27) + ((t288 * t423 - t289 * t419) * t384 + (t310 * t423 - t318 * t419) * t350 + t464 * t421 + t284 * t499 + t322 * t312 + t354 * t287 + t290 * t333 + t321 * t294 - g(1) * t330 - g(2) * t332 + ((-t310 * t419 - t318 * t423) * t384 - t285 * t421) * qJD(6)) * MDP(29) + (t355 * t421 + t391 * t499) * MDP(21) + (t350 * t421 + t384 * t499) * MDP(28); t457 * MDP(6) + t448 * MDP(9) + t444 * MDP(22) + t445 * MDP(23) + (-MDP(6) * qJ(2) - MDP(5) - MDP(7)) * t428 + (MDP(29) * t361 + MDP(30) * t362) * t350 + (MDP(29) * t513 - MDP(30) * t514) * t384 + (-MDP(15) * t421 - MDP(16) * t425 + MDP(4) - MDP(8)) * qJDD(1) + ((-qJD(3) - t395) * MDP(9) + (0.2e1 * qJD(4) * MDP(16) + (t420 * MDP(22) + t424 * MDP(23)) * t391) * t421 + (-0.2e1 * qJD(4) * MDP(15) + t358 * MDP(22) + t360 * MDP(23) + t312 * MDP(29) - MDP(30) * t449) * t425) * qJD(1); qJDD(1) * MDP(7) - t428 * MDP(8) + (-t438 + t466) * MDP(9) + (t421 * t507 + t398) * MDP(15) + (t425 * t507 - t483) * MDP(16) + (-t425 * t308 + (t502 - t525) * t421 + (-t424 * t459 - t425 * t491) * t391) * MDP(22) + (-t307 * t425 + (t501 - t532) * t421 + (t420 * t459 - t472) * t391) * MDP(23) + (qJD(1) * t441 + (-qJD(4) * t362 * t384 - t287) * t425 + ((t419 * t498 + t437 + t469) * t384 - t533 + qJD(4) * t312) * t421) * MDP(29) + (t384 * t439 + (qJD(4) * t441 - t286) * t425 + (-(t543 * t419 - t420 * t493 - t423 * t498) * t384 + t534 - qJD(4) * t449) * t421) * MDP(30); MDP(12) * t485 - MDP(13) * t486 + qJDD(4) * MDP(14) + (-t425 * t433 + t539) * MDP(15) + (t421 * t433 + t538) * MDP(16) + (t391 * t528 + t535) * MDP(17) + ((t307 - t531) * t424 + (-t308 - t530) * t420) * MDP(18) + ((-t360 * t425 + t391 * t522) * qJD(1) + t445) * MDP(19) + ((-t391 * t420 * t421 + t358 * t425) * qJD(1) - t444) * MDP(20) + (-t358 * t368 - pkin(4) * t308 - t349 * t391 + (t391 * t527 + t434) * t420 - t550 * t424) * MDP(22) + (-pkin(4) * t307 - t360 * t368 + t512 * t391 + t420 * t550 + t434 * t424) * MDP(23) + (t286 * t362 + t449 * t514) * MDP(24) + (-t286 * t361 - t287 * t362 + t312 * t514 + t449 * t513) * MDP(25) + (-t384 * t514 + t533) * MDP(26) + (-t384 * t513 - t534) * MDP(27) + ((-t376 * t423 - t377 * t419) * t350 + t394 * t287 + t290 * t361 + (t419 * t451 - t423 * t452) * t384 + t513 * t321 + t456 * t312 - t432 * t401) * MDP(29) + (-(-t376 * t419 + t377 * t423) * t350 + t394 * t286 + t290 * t362 + (t419 * t452 + t423 * t451) * t384 - t514 * t321 - t456 * t449 + t432 * t400) * MDP(30) + (-t391 * MDP(21) - MDP(22) * t298 + t299 * MDP(23) + MDP(26) * t449 + t312 * MDP(27) - t384 * MDP(28) - t284 * MDP(29) + t285 * MDP(30)) * t504 + (MDP(10) * t421 * t425 - MDP(11) * t508) * t428; t360 * t358 * MDP(17) + (-t358 ^ 2 + t360 ^ 2) * MDP(18) + (t307 + t531) * MDP(19) + (-t308 + t530) * MDP(20) + t355 * MDP(21) + (-t342 * t496 - g(1) * t346 + g(2) * t344 + t299 * t391 - t343 * t360 + t301 + (t461 + t538) * t420) * MDP(22) + (g(1) * t347 - g(2) * t345 + g(3) * t517 + t298 * t391 + t343 * t358 - t440) * MDP(23) + (t286 + t549) * MDP(26) + (-t287 - t548) * MDP(27) + (-(-t295 * t419 - t536) * t384 - t285 * qJD(6) + (-t312 * t360 + t423 * t350 - t384 * t494) * pkin(5) + t546) * MDP(29) + ((-t296 * t384 - t282) * t419 + (t295 * t384 - t460) * t423 + (-t419 * t350 + t360 * t449 - t384 * t493) * pkin(5) + t547) * MDP(30) + t545; (t478 + t549) * MDP(26) + (-t463 - t548) * MDP(27) + (t285 * t384 + t546) * MDP(29) + (-t419 * t282 - t423 * t283 + t284 * t384 + t547) * MDP(30) + (-MDP(26) * t529 + MDP(27) * t449 - MDP(29) * t285 - MDP(30) * t537) * qJD(6) + t545;];
tau  = t1;

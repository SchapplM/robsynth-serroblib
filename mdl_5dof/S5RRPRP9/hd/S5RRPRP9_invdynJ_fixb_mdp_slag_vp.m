% Calculate vector of inverse dynamics joint torques for
% S5RRPRP9
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP9_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:40
% EndTime: 2019-12-31 20:07:48
% DurationCPUTime: 6.33s
% Computational Cost: add. (3770->475), mult. (8404->593), div. (0->0), fcn. (5789->10), ass. (0->192)
t433 = sin(pkin(8));
t532 = pkin(7) + qJ(3);
t396 = t532 * t433;
t434 = cos(pkin(8));
t397 = t532 * t434;
t436 = sin(qJ(4));
t541 = cos(qJ(4));
t345 = -t436 * t396 + t397 * t541;
t439 = cos(qJ(2));
t424 = t439 * qJDD(1);
t437 = sin(qJ(2));
t498 = qJD(1) * qJD(2);
t461 = t437 * t498 - t424;
t384 = qJDD(4) + t461;
t430 = pkin(8) + qJ(4);
t422 = sin(t430);
t429 = g(3) * t439;
t438 = sin(qJ(1));
t440 = cos(qJ(1));
t479 = g(1) * t440 + g(2) * t438;
t472 = t479 * t437;
t453 = t472 - t429;
t554 = t345 * t384 + t422 * t453;
t506 = qJD(1) * t437;
t489 = t433 * t506;
t503 = qJD(2) * t434;
t381 = -t489 + t503;
t492 = t434 * t506;
t504 = qJD(2) * t433;
t382 = t492 + t504;
t331 = -t541 * t381 + t382 * t436;
t553 = t331 ^ 2;
t466 = -t436 * t381 - t382 * t541;
t542 = t466 ^ 2;
t505 = qJD(1) * t439;
t413 = -qJD(4) + t505;
t552 = t331 * t413;
t551 = t413 * t466;
t493 = t541 * t434;
t464 = -t436 * t433 + t493;
t488 = qJD(4) * t541;
t500 = qJD(4) * t436;
t546 = -t433 * t500 + t434 * t488;
t512 = -t464 * t505 + t546;
t386 = t433 * t541 + t436 * t434;
t373 = t386 * qJD(4);
t455 = t439 * t386;
t511 = -qJD(1) * t455 + t373;
t486 = t439 * t498;
t497 = qJDD(1) * t437;
t462 = t486 + t497;
t476 = pkin(2) * t437 - qJ(3) * t439;
t388 = t476 * qJD(1);
t350 = pkin(6) * t489 + t434 * t388;
t523 = t434 * t439;
t474 = pkin(3) * t437 - pkin(7) * t523;
t323 = qJD(1) * t474 + t350;
t374 = t433 * t388;
t524 = t434 * t437;
t525 = t433 * t439;
t463 = -pkin(6) * t524 - pkin(7) * t525;
t336 = qJD(1) * t463 + t374;
t465 = -t396 * t541 - t436 * t397;
t550 = qJD(3) * t464 + qJD(4) * t465 - t436 * t323 - t541 * t336;
t549 = qJD(3) * t386 + qJD(4) * t345 + t323 * t541 - t436 * t336;
t477 = pkin(2) * t439 + qJ(3) * t437;
t393 = -pkin(1) - t477;
t380 = t434 * t393;
t337 = -pkin(7) * t524 + t380 + (-pkin(6) * t433 - pkin(3)) * t439;
t354 = pkin(6) * t523 + t433 * t393;
t526 = t433 * t437;
t343 = -pkin(7) * t526 + t354;
t548 = t436 * t337 + t541 * t343;
t545 = pkin(6) * t486 + qJDD(3);
t496 = qJDD(2) * t433;
t447 = t434 * t462 + t496;
t510 = t462 * t433;
t475 = qJDD(2) * t434 - t510;
t295 = -qJD(4) * t466 + t436 * t447 - t541 * t475;
t370 = qJD(2) * t476 - qJD(3) * t437;
t502 = qJD(2) * t437;
t494 = pkin(6) * t502;
t341 = t434 * t370 + t433 * t494;
t315 = qJD(2) * t474 + t341;
t357 = t433 * t370;
t324 = qJD(2) * t463 + t357;
t544 = -qJD(4) * t548 + t315 * t541 - t436 * t324;
t543 = -0.2e1 * pkin(1);
t540 = pkin(3) * t433;
t539 = pkin(4) * t384;
t538 = pkin(6) * t381;
t537 = g(1) * t438;
t534 = g(2) * t440;
t533 = g(3) * t437;
t531 = qJ(5) * t384;
t530 = qJDD(2) * pkin(2);
t376 = t393 * qJD(1);
t420 = pkin(6) * t505;
t398 = qJD(2) * qJ(3) + t420;
t338 = t434 * t376 - t398 * t433;
t495 = pkin(3) * t505;
t307 = -pkin(7) * t382 + t338 - t495;
t339 = t433 * t376 + t434 * t398;
t311 = pkin(7) * t381 + t339;
t287 = t436 * t307 + t311 * t541;
t529 = t287 * t413;
t528 = t331 * t466;
t522 = t437 * t440;
t423 = cos(t430);
t521 = t438 * t423;
t520 = t438 * t439;
t519 = t439 * t440;
t518 = t440 * t422;
t377 = t433 * t495 + t420;
t517 = t511 * pkin(4) - t512 * qJ(5) - qJD(5) * t386 - t377;
t516 = -qJ(5) * t506 + t550;
t515 = pkin(4) * t506 + t549;
t328 = qJD(1) * t370 + qJDD(1) * t393;
t359 = -pkin(6) * t461 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t303 = t433 * t328 + t434 * t359;
t418 = pkin(6) * t497;
t509 = -t418 - t429;
t501 = qJD(2) * t439;
t491 = t433 * t501;
t378 = pkin(3) * t491 + pkin(6) * t501;
t389 = pkin(3) * t526 + t437 * pkin(6);
t508 = t440 * pkin(1) + t438 * pkin(6);
t431 = t437 ^ 2;
t507 = -t439 ^ 2 + t431;
t286 = t307 * t541 - t436 * t311;
t499 = qJD(5) - t286;
t417 = pkin(3) * t434 + pkin(2);
t485 = -qJD(2) * pkin(2) + qJD(3);
t302 = t434 * t328 - t433 * t359;
t290 = pkin(3) * t461 - pkin(7) * t447 + t302;
t296 = pkin(7) * t475 + t303;
t484 = -t541 * t290 + t436 * t296 + t307 * t500 + t311 * t488;
t415 = t437 * t537;
t482 = -g(2) * t522 + t415;
t360 = t422 * t520 + t423 * t440;
t362 = t439 * t518 - t521;
t481 = -g(1) * t360 + g(2) * t362;
t361 = t423 * t520 - t518;
t363 = t422 * t438 + t423 * t519;
t480 = g(1) * t361 - g(2) * t363;
t478 = -t534 + t537;
t392 = pkin(6) * t506 + t485;
t473 = pkin(4) * t423 + qJ(5) * t422 + t417;
t368 = t418 - t530 + t545;
t471 = -pkin(6) * qJDD(2) + t498 * t543;
t469 = t337 * t541 - t436 * t343;
t460 = t434 * t497 + t496;
t459 = t436 * t290 + t541 * t296 + t307 * t488 - t311 * t500;
t458 = t436 * t315 + t541 * t324 + t337 * t488 - t343 * t500;
t294 = -t381 * t488 + t382 * t500 - t436 * t475 - t541 * t447;
t442 = qJD(1) ^ 2;
t457 = pkin(1) * t442 + t479;
t349 = -pkin(3) * t381 + t392;
t441 = qJD(2) ^ 2;
t456 = pkin(6) * t441 + qJDD(1) * t543 + t534;
t454 = g(2) * t437 * t521 + t465 * t384 + (g(1) * t522 - t429) * t423;
t452 = -t439 * t479 - t533;
t451 = -t472 - t530;
t449 = g(1) * t362 + g(2) * t360 + t422 * t533 - t484;
t448 = -t368 + t453;
t316 = -pkin(3) * t475 + t368;
t291 = pkin(4) * t331 + qJ(5) * t466 + t349;
t446 = -t291 * t466 + qJDD(5) - t449;
t445 = -g(1) * t363 - g(2) * t361 - t423 * t533 + t459;
t279 = t295 * pkin(4) + t294 * qJ(5) + qJD(5) * t466 + t316;
t427 = t440 * pkin(6);
t367 = t464 * t437;
t366 = t386 * t437;
t353 = -pkin(6) * t525 + t380;
t351 = -pkin(6) * t492 + t374;
t342 = -t434 * t494 + t357;
t329 = -pkin(4) * t464 - qJ(5) * t386 - t417;
t319 = qJD(2) * t455 + t437 * t546;
t318 = t373 * t437 + t436 * t491 - t493 * t501;
t308 = pkin(4) * t366 - qJ(5) * t367 + t389;
t301 = -pkin(4) * t466 + qJ(5) * t331;
t300 = t439 * pkin(4) - t469;
t299 = -qJ(5) * t439 + t548;
t285 = -t413 * qJ(5) + t287;
t284 = t413 * pkin(4) + t499;
t283 = pkin(4) * t319 + qJ(5) * t318 - qJD(5) * t367 + t378;
t282 = -t294 - t552;
t281 = -pkin(4) * t502 - t544;
t280 = qJ(5) * t502 - qJD(5) * t439 + t458;
t278 = qJDD(5) + t484 - t539;
t277 = -qJD(5) * t413 + t459 + t531;
t1 = [(t277 * t299 + t285 * t280 + t279 * t308 + t291 * t283 + t278 * t300 + t284 * t281 - g(1) * (-pkin(4) * t361 - qJ(5) * t360 + t440 * t540 + t427) - g(2) * (pkin(4) * t363 + qJ(5) * t362 + t417 * t519 + t522 * t532 + t508) + (-g(1) * (-t417 * t439 - t437 * t532 - pkin(1)) - g(2) * t540) * t438) * MDP(25) + qJDD(1) * MDP(1) + (qJDD(1) * t431 + 0.2e1 * t437 * t486) * MDP(4) + 0.2e1 * (t424 * t437 - t498 * t507) * MDP(5) + (qJDD(2) * t437 + t439 * t441) * MDP(6) + (qJDD(2) * t439 - t437 * t441) * MDP(7) + (t471 * t437 + (-t456 + t537) * t439) * MDP(9) + (t437 * t456 + t439 * t471 - t415) * MDP(10) + (-t479 * t433 + (-pkin(6) * t475 + t368 * t433 + (t353 * qJD(1) + t338) * qJD(2)) * t437 + (-t341 * qJD(1) - t353 * qJDD(1) - t302 + t478 * t434 + (t392 * t433 - t538) * qJD(2)) * t439) * MDP(11) + (-t479 * t434 + (t368 * t434 + (-qJD(1) * t354 - t339) * qJD(2) + t460 * pkin(6)) * t437 + (t342 * qJD(1) + t354 * qJDD(1) + t303 - t478 * t433 + (t392 * t434 + (t382 + t492) * pkin(6)) * qJD(2)) * t439) * MDP(12) + (t342 * t381 - t354 * t510 - t341 * t382 + (-t353 * qJDD(2) - t303 * t437 - t339 * t501) * t433 + (t354 * qJDD(2) - t302 * t437 - t338 * t501 - t353 * t462) * t434 + t482) * MDP(13) + (t303 * t354 + t339 * t342 + t302 * t353 + t338 * t341 - g(1) * t427 - g(2) * (t440 * t477 + t508) - t393 * t537 + (t368 * t437 + t392 * t501) * pkin(6)) * MDP(14) + (-t294 * t367 + t318 * t466) * MDP(15) + (t294 * t366 - t295 * t367 + t318 * t331 + t319 * t466) * MDP(16) + (t294 * t439 + t318 * t413 + t367 * t384 - t466 * t502) * MDP(17) + (t295 * t439 + t319 * t413 - t331 * t502 - t366 * t384) * MDP(18) + (-t384 * t439 - t413 * t502) * MDP(19) + (t286 * t502 + t389 * t295 + t316 * t366 + t349 * t319 + t378 * t331 + t469 * t384 - t413 * t544 + t484 * t439 + t480) * MDP(20) + (-t287 * t502 - t389 * t294 + t316 * t367 - t349 * t318 - t378 * t466 - t384 * t548 + t413 * t458 + t439 * t459 + t481) * MDP(21) + (t278 * t439 + t279 * t366 + t281 * t413 + t283 * t331 - t284 * t502 + t291 * t319 + t295 * t308 - t300 * t384 + t480) * MDP(22) + (-t277 * t366 + t278 * t367 - t280 * t331 - t281 * t466 - t284 * t318 - t285 * t319 - t294 * t300 - t295 * t299 + t482) * MDP(23) + (-t277 * t439 - t279 * t367 - t280 * t413 + t283 * t466 + t285 * t502 + t291 * t318 + t294 * t308 + t299 * t384 - t481) * MDP(24) + t478 * MDP(2) + t479 * MDP(3); MDP(6) * t497 + MDP(7) * t424 + qJDD(2) * MDP(8) + (t437 * t457 + t509) * MDP(9) + (t533 + (-pkin(6) * qJDD(1) + t457) * t439) * MDP(10) + (t433 * qJ(3) * t424 - pkin(2) * t510 + (t448 + t530) * t434 + ((-qJ(3) * t504 - t338) * t437 + (t538 + t350 + (qJD(3) - t392) * t433) * t439) * qJD(1)) * MDP(11) + (-t476 * t434 * qJDD(1) + (t368 + t451 + t429) * t433 + ((-qJ(3) * t503 + t339) * t437 + (-pkin(6) * t382 - t351 + (-t392 + t485) * t434) * t439) * qJD(1)) * MDP(12) + (t350 * t382 - t351 * t381 + (qJ(3) * t475 + qJD(3) * t381 + t338 * t505 + t303) * t434 + (qJ(3) * t447 + qJD(3) * t382 + t339 * t505 - t302) * t433 + t452) * MDP(13) + (-t392 * t420 - t338 * t350 - t339 * t351 + (-t338 * t433 + t339 * t434) * qJD(3) + t448 * pkin(2) + (-t302 * t433 + t303 * t434 + t452) * qJ(3)) * MDP(14) + (-t294 * t386 - t466 * t512) * MDP(15) + (-t294 * t464 - t295 * t386 - t331 * t512 + t466 * t511) * MDP(16) + (t384 * t386 - t413 * t512 + t466 * t506) * MDP(17) + (t331 * t506 + t384 * t464 + t413 * t511) * MDP(18) + t413 * MDP(19) * t506 + (-t286 * t506 - t417 * t295 - t316 * t464 - t377 * t331 + t511 * t349 + t413 * t549 + t454) * MDP(20) + (t287 * t506 + t417 * t294 + t316 * t386 + t512 * t349 + t377 * t466 + t413 * t550 - t554) * MDP(21) + (-t279 * t464 + t284 * t506 + t291 * t511 + t295 * t329 + t331 * t517 + t413 * t515 + t454) * MDP(22) + (t277 * t464 + t278 * t386 + t284 * t512 - t285 * t511 + t294 * t465 - t295 * t345 - t331 * t516 - t466 * t515 + t452) * MDP(23) + (-t279 * t386 - t285 * t506 - t291 * t512 + t294 * t329 - t413 * t516 + t466 * t517 + t554) * MDP(24) + (t277 * t345 - t278 * t465 + t279 * t329 + t517 * t291 + t516 * t285 + t515 * t284 + (-g(3) * t473 - t479 * t532) * t439 + (-g(3) * t532 + t473 * t479) * t437) * MDP(25) + (-MDP(4) * t437 * t439 + MDP(5) * t507) * t442; (-t382 * t505 - t475) * MDP(11) + ((-t381 + t503) * t505 + t460) * MDP(12) + (-t381 ^ 2 - t382 ^ 2) * MDP(13) + (t338 * t382 - t339 * t381 + t451 - t509 + t545) * MDP(14) + (-t542 - t553) * MDP(23) + (t284 * t466 + t285 * t331 + t279 - t453) * MDP(25) + (-MDP(21) + MDP(24)) * (t294 - t552) + (MDP(20) + MDP(22)) * (t295 + t551); -MDP(15) * t528 + (t542 - t553) * MDP(16) + t282 * MDP(17) + (-t295 + t551) * MDP(18) + t384 * MDP(19) + (t349 * t466 + t449 - t529) * MDP(20) + (-t286 * t413 + t331 * t349 - t445) * MDP(21) + (-t301 * t331 - t446 - t529 + 0.2e1 * t539) * MDP(22) + (pkin(4) * t294 - qJ(5) * t295 - (t285 - t287) * t466 + (t284 - t499) * t331) * MDP(23) + (0.2e1 * t531 - t291 * t331 - t301 * t466 + (-0.2e1 * qJD(5) + t286) * t413 + t445) * MDP(24) + (t277 * qJ(5) - t278 * pkin(4) - t291 * t301 - t284 * t287 - g(1) * (-pkin(4) * t362 + qJ(5) * t363) - g(2) * (-pkin(4) * t360 + qJ(5) * t361) - (-pkin(4) * t422 + qJ(5) * t423) * t533 + t499 * t285) * MDP(25); (-t384 - t528) * MDP(22) + t282 * MDP(23) + (-t413 ^ 2 - t542) * MDP(24) + (t285 * t413 + t446 - t539) * MDP(25);];
tau = t1;

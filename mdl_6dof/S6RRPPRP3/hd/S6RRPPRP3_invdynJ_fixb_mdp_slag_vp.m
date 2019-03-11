% Calculate vector of inverse dynamics joint torques for
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:35:45
% EndTime: 2019-03-09 08:35:51
% DurationCPUTime: 5.70s
% Computational Cost: add. (2630->517), mult. (5106->605), div. (0->0), fcn. (2755->6), ass. (0->228)
t433 = sin(qJ(2));
t521 = qJD(1) * t433;
t398 = pkin(7) * t521;
t358 = -qJ(4) * t521 + t398;
t493 = qJD(3) + t358;
t564 = pkin(2) + pkin(3);
t494 = t564 * qJD(2);
t336 = -t494 + t493;
t573 = t433 * t564;
t436 = cos(qJ(2));
t507 = qJD(1) * qJD(2);
t491 = t436 * t507;
t505 = qJDD(1) * t433;
t572 = t491 + t505;
t520 = qJD(1) * t436;
t571 = qJDD(2) * t564;
t399 = pkin(7) * t520;
t360 = -qJ(4) * t520 + t399;
t423 = qJD(2) * qJ(3);
t345 = -t360 - t423;
t434 = sin(qJ(1));
t417 = g(2) * t434;
t437 = cos(qJ(1));
t570 = g(1) * t437 + t417;
t432 = sin(qJ(5));
t435 = cos(qJ(5));
t534 = t435 * t437;
t541 = t433 * t434;
t348 = t432 * t541 - t534;
t537 = t434 * t435;
t539 = t433 * t437;
t350 = -t432 * t539 - t537;
t569 = -g(1) * t350 + g(2) * t348;
t516 = qJD(2) * t432;
t356 = t435 * t520 + t516;
t492 = t433 * t507;
t504 = qJDD(1) * t436;
t318 = qJD(5) * t356 - t435 * qJDD(2) + (-t492 + t504) * t432;
t568 = -pkin(5) * t318 + qJDD(6);
t372 = qJ(4) * t492;
t515 = qJD(2) * t433;
t502 = pkin(7) * t515;
t463 = qJD(4) * t436 + t502;
t395 = pkin(7) * t504;
t420 = qJDD(2) * qJ(3);
t421 = qJD(2) * qJD(3);
t497 = t395 + t420 + t421;
t313 = qJ(4) * t504 + qJD(1) * t463 - t372 - t497;
t426 = qJDD(2) * pkin(4);
t311 = -t313 + t426;
t378 = qJD(5) + t521;
t422 = -pkin(8) - t564;
t558 = g(3) * t433;
t567 = qJD(5) * t422 * t378 + t436 * t570 - t311 + t558;
t347 = -qJD(1) * pkin(1) - pkin(2) * t520 - qJ(3) * t521;
t406 = t433 * qJ(3);
t413 = t436 * pkin(2);
t525 = t413 + t406;
t364 = -pkin(1) - t525;
t556 = pkin(7) * qJDD(2);
t566 = (qJD(1) * t364 + t347) * qJD(2) - t556;
t565 = t356 ^ 2;
t562 = pkin(5) * t432;
t418 = g(1) * t434;
t559 = g(2) * t437;
t416 = g(3) * t436;
t412 = t436 * pkin(3);
t557 = pkin(7) - qJ(4);
t393 = t435 * pkin(5) + pkin(4);
t555 = qJ(3) * t436;
t427 = qJDD(1) * pkin(1);
t554 = qJDD(2) * pkin(2);
t510 = t435 * qJD(2);
t513 = qJD(5) * t432;
t453 = t433 * t510 + t436 * t513;
t317 = qJD(1) * t453 - qJD(5) * t510 - t432 * qJDD(2) - t435 * t504;
t553 = t317 * t432;
t552 = t318 * t435;
t354 = qJDD(5) + t572;
t551 = t354 * t435;
t355 = t432 * t520 - t510;
t550 = t355 * t378;
t549 = t355 * t432;
t548 = t355 * t435;
t547 = t356 * t378;
t546 = t356 * t432;
t545 = t356 * t435;
t544 = t432 * t354;
t543 = t432 * t436;
t542 = t432 * t437;
t540 = t433 * t435;
t440 = qJD(1) ^ 2;
t538 = t433 * t440;
t536 = t434 * t436;
t367 = t557 * t433;
t357 = t435 * t367;
t535 = t435 * t436;
t533 = t436 * t437;
t532 = qJ(6) - t422;
t333 = pkin(3) * t520 + qJD(4) - t347;
t477 = pkin(4) * t433 + pkin(8) * t436;
t323 = qJD(1) * t477 + t333;
t331 = qJD(2) * t422 + t493;
t302 = t435 * t323 - t331 * t432;
t299 = qJ(6) * t356 + t302;
t295 = pkin(5) * t378 + t299;
t531 = -t299 + t295;
t484 = qJD(5) * t532;
t390 = qJ(3) * t520;
t459 = pkin(4) * t436 + t422 * t433;
t326 = qJD(1) * t459 + t390;
t528 = t432 * t326 + t435 * t360;
t530 = -qJD(6) * t435 - t528 + (qJ(6) * t521 + t484) * t432;
t325 = t435 * t326;
t465 = pkin(5) * t436 - qJ(6) * t540;
t529 = -qJD(1) * t465 + t435 * t484 - t325 + (qJD(6) + t360) * t432;
t496 = t412 + t525;
t352 = pkin(1) + t496;
t330 = t477 + t352;
t527 = t432 * t330 + t357;
t404 = t433 * qJD(3);
t514 = qJD(2) * t436;
t526 = qJ(3) * t514 + t404;
t424 = t433 ^ 2;
t425 = t436 ^ 2;
t523 = t424 - t425;
t522 = t424 + t425;
t519 = qJD(2) * t345;
t518 = qJD(2) * t355;
t517 = qJD(2) * t356;
t512 = qJD(5) * t435;
t511 = qJD(6) * t436;
t508 = -qJD(4) - t333;
t506 = qJD(1) * qJD(4);
t430 = -qJ(6) - pkin(8);
t503 = -t430 + t564;
t501 = qJ(6) * t535;
t500 = t436 * t538;
t449 = t459 * qJD(2);
t322 = t449 + t526;
t343 = -qJD(4) * t433 + t514 * t557;
t499 = t432 * t322 + t330 * t512 + t435 * t343;
t498 = t395 + 0.2e1 * t420 + 0.2e1 * t421;
t495 = -g(1) * t539 - g(2) * t541 + t416;
t377 = pkin(7) * t491;
t394 = pkin(7) * t505;
t490 = qJDD(3) + t377 + t394;
t489 = -pkin(1) - t406;
t488 = t418 - t559;
t487 = -qJD(2) * pkin(2) + qJD(3);
t414 = t437 * pkin(7);
t486 = -qJ(4) * t437 + t414;
t485 = g(1) * t503;
t483 = qJD(1) * t352 + t333;
t443 = -qJ(4) * t572 - t433 * t506 + t490;
t310 = qJDD(2) * t422 + t443;
t481 = -qJD(5) * t323 - t310;
t480 = t437 * pkin(1) + pkin(2) * t533 + t434 * pkin(7) + qJ(3) * t539;
t479 = t394 + t495;
t478 = t433 * t494;
t439 = qJD(2) ^ 2;
t476 = pkin(7) * t439 + t559;
t473 = pkin(2) * t504 + qJ(3) * t572 + qJD(1) * t404 + t427;
t454 = pkin(3) * t504 + qJDD(4) + t473;
t301 = qJD(1) * t449 + qJDD(1) * t477 + t454;
t298 = t435 * t301;
t474 = -t331 * t512 + t298;
t472 = pkin(3) * t533 + t480;
t303 = t323 * t432 + t331 * t435;
t300 = qJ(6) * t355 + t303;
t471 = -t295 * t435 - t300 * t432;
t470 = -t295 * t432 + t300 * t435;
t363 = t398 + t487;
t365 = t399 + t423;
t469 = t363 * t436 - t365 * t433;
t468 = -qJDD(3) - t479;
t467 = MDP(24) * t435 - MDP(25) * t432;
t466 = t489 - t413;
t464 = -0.2e1 * pkin(1) * t507 - t556;
t462 = t377 - t468;
t461 = -t378 * t512 - t544;
t460 = t378 * t513 - t551;
t339 = qJD(2) * pkin(4) - t345;
t457 = -t432 * t301 - t435 * t310 - t323 * t512 + t331 * t513;
t455 = -t476 + 0.2e1 * t427;
t452 = -qJ(4) * qJDD(1) - t570;
t291 = pkin(5) * t354 - qJ(6) * t317 - qJD(5) * t303 + qJD(6) * t356 - t310 * t432 + t298;
t292 = qJ(6) * t318 + qJD(6) * t355 - t457;
t451 = t291 * t435 + t292 * t432 - t559;
t334 = -pkin(7) * t492 + t497;
t448 = -t339 * t378 - t422 * t354;
t307 = -qJD(1) * t478 + t454;
t332 = -t478 + t526;
t447 = -qJD(1) * t332 - qJDD(1) * t352 - t307 + t559;
t446 = t462 - t554;
t320 = pkin(2) * t492 - t473;
t344 = pkin(2) * t515 - t526;
t445 = -qJD(1) * t344 - qJDD(1) * t364 - t320 - t476;
t340 = t490 - t554;
t444 = qJD(2) * t469 + t334 * t436 + t340 * t433;
t442 = (-t546 + t548) * MDP(26) + t470 * MDP(27) + (-t432 * MDP(24) - t435 * MDP(25)) * t378;
t441 = (-t545 - t549) * MDP(26) + t471 * MDP(27) - t467 * t378;
t431 = qJ(3) + pkin(4);
t407 = t436 * qJ(4);
t391 = qJ(4) * t515;
t383 = g(1) * t536;
t382 = g(1) * t541;
t376 = qJ(3) * t533;
t374 = qJ(3) * t536;
t368 = pkin(7) * t436 - t407;
t362 = t532 * t435;
t361 = t532 * t432;
t359 = pkin(2) * t521 - t390;
t353 = t355 ^ 2;
t351 = -t432 * t434 + t433 * t534;
t349 = -t433 * t537 - t542;
t342 = -t391 + t463;
t341 = -t521 * t564 + t390;
t329 = t435 * t330;
t321 = -pkin(5) * t355 + qJD(6) + t339;
t316 = t435 * t322;
t312 = t443 - t571;
t309 = qJ(6) * t543 + t527;
t306 = pkin(5) * t433 - t367 * t432 + t329 + t501;
t296 = t311 + t568;
t294 = qJD(5) * t501 + (-qJ(6) * t515 - qJD(5) * t367 + t511) * t432 + t499;
t293 = t435 * t511 - t343 * t432 + t316 + t465 * qJD(2) + (-t357 + (-qJ(6) * t436 - t330) * t432) * qJD(5);
t1 = [(t292 * t309 + t300 * t294 + t291 * t306 + t295 * t293 - t296 * t407 + t321 * (t515 * t562 + t391 - t502) - g(1) * (-pkin(5) * t542 + t486) - g(2) * (t393 * t539 + t472) + (t296 * (pkin(7) - t562) + t321 * (-pkin(5) * t512 - qJD(4)) + t430 * t559) * t436 + (-g(1) * (-t393 * t433 + t489) - g(2) * (-qJ(4) - t562) + t436 * t485) * t434) * MDP(27) + (t433 * t566 + t445 * t436 + t383) * MDP(11) + (t445 * t433 - t436 * t566 + t382) * MDP(13) + ((-qJD(2) * t336 - qJDD(1) * t368 + t313 + (-qJD(2) * t367 + t342) * qJD(1)) * t436 + (-t519 - qJDD(1) * t367 - t312 + (qJD(2) * t368 - t343) * qJD(1)) * t433 + t570) * MDP(17) + (pkin(7) * qJDD(1) * t522 + t444 - t570) * MDP(12) + t570 * MDP(3) + (qJDD(1) * t424 + 0.2e1 * t433 * t491) * MDP(4) + t488 * MDP(2) + (-(-t367 * t513 + t499) * t378 - t527 * t354 + t342 * t356 + t368 * t317 - g(1) * t348 - g(2) * t350 + (t339 * t510 + t457) * t433 + (-qJD(2) * t303 - t311 * t435 + t339 * t513) * t436) * MDP(25) + 0.2e1 * (t433 * t504 - t507 * t523) * MDP(5) + (t354 * t433 + t378 * t514) * MDP(23) + (t293 * t356 + t294 * t355 - t306 * t317 + t309 * t318 + t383 + t471 * t515 + (qJD(5) * t470 + t451) * t436) * MDP(26) + ((t378 * t510 + t317) * t433 + (t460 - t517) * t436) * MDP(21) + ((-t378 * t516 + t318) * t433 + (-t461 + t518) * t436) * MDP(22) + ((-t367 * t512 + t316) * t378 + t329 * t354 + t474 * t433 + t342 * t355 - t368 * t318 - g(1) * t349 - g(2) * t351 + (qJD(2) * t302 - t339 * t512) * t436 + ((-qJD(5) * t330 - t343) * t378 - t367 * t354 - t311 * t436 + (qJD(2) * t339 + t481) * t433) * t432) * MDP(24) + (qJDD(2) * t367 - t383 + (t433 * t483 + t343) * qJD(2) + t447 * t436) * MDP(16) + (qJDD(2) * t368 + t382 + (t436 * t483 - t342) * qJD(2) - t447 * t433) * MDP(15) + (qJDD(2) * t433 + t436 * t439) * MDP(6) + (qJDD(2) * t436 - t433 * t439) * MDP(7) + (pkin(7) * t444 - g(1) * t414 - g(2) * t480 + t320 * t364 + t347 * t344 - t418 * t466) * MDP(14) + (t312 * t367 + t336 * t343 - t313 * t368 + t345 * t342 + t307 * t352 + t333 * t332 - g(1) * t486 - g(2) * t472 + (-g(1) * (t466 - t412) + g(2) * qJ(4)) * t434) * MDP(18) + ((t546 + t548) * t515 + (t553 - t552 + (-t545 + t549) * qJD(5)) * t436) * MDP(20) + (t433 * t464 + t436 * t455 + t383) * MDP(9) + (-t433 * t455 + t436 * t464 - t382) * MDP(10) + qJDD(1) * MDP(1) + (-t317 * t535 - t356 * t453) * MDP(19); -MDP(4) * t500 - t378 * MDP(23) * t520 + t523 * MDP(5) * t440 + MDP(6) * t505 + MDP(7) * t504 + qJDD(2) * MDP(8) + (pkin(1) * t538 - t479) * MDP(9) + (t558 - t395 + (pkin(1) * t440 + t570) * t436) * MDP(10) + (0.2e1 * t554 + (-t347 * t433 + t359 * t436) * qJD(1) + t468) * MDP(11) + ((-pkin(2) * t433 + t555) * qJDD(1) + ((t365 - t423) * t433 + (-t363 + t487) * t436) * qJD(1)) * MDP(12) + ((qJD(1) * t359 - g(3)) * t433 + (qJD(1) * t347 - t570) * t436 + t498) * MDP(13) + (t334 * qJ(3) + t365 * qJD(3) - t340 * pkin(2) - t347 * t359 - g(1) * (-pkin(2) * t539 + t376) - g(2) * (-pkin(2) * t541 + t374) - g(3) * t525 - t469 * qJD(1) * pkin(7)) * MDP(14) + (qJD(2) * t358 + t372 + (-g(3) + (-pkin(7) * qJD(2) - t341) * qJD(1)) * t433 + (qJD(1) * t508 + t452) * t436 + t498) * MDP(15) + (-qJ(4) * t505 - qJD(2) * t360 - 0.2e1 * t571 + ((-qJ(4) * qJD(2) + t341) * t436 + t508 * t433) * qJD(1) + t462) * MDP(16) + (-t555 + t573) * qJDD(1) * MDP(17) + (-g(1) * t376 - g(2) * t374 - g(3) * t496 - t313 * qJ(3) - t312 * t564 - t333 * t341 - t336 * t360 - t345 * t493 + t570 * t573) * MDP(18) + (t378 * t545 - t553) * MDP(19) + ((-t317 - t550) * t435 + (-t318 - t547) * t432) * MDP(20) + ((t356 * t436 - t378 * t540) * qJD(1) + t461) * MDP(21) + ((t378 * t432 * t433 - t355 * t436) * qJD(1) + t460) * MDP(22) + (-t302 * t520 - t431 * t318 - t325 * t378 - t493 * t355 + (t360 * t378 + t448) * t432 - t567 * t435) * MDP(24) + (t303 * t520 + t431 * t317 - t356 * t493 + t528 * t378 + t432 * t567 + t448 * t435) * MDP(25) + (-t317 * t361 - t318 * t362 + t529 * t356 + t530 * t355 + (t295 * t378 - t292) * t435 + (t300 * t378 + t291) * t432 - t495) * MDP(26) + (-t292 * t362 + t291 * t361 + t296 * (qJ(3) + t393) - g(1) * (t393 * t533 + t376) - g(2) * (t393 * t536 + t374) - g(3) * (-t430 * t436 + t496) + (-pkin(5) * t513 + t493) * t321 + t530 * t300 + t529 * t295 + (-qJD(1) * t321 * t562 - g(3) * t393 + t417 * t503 + t437 * t485) * t433) * MDP(27); (-qJD(2) * t365 + t446) * MDP(14) + (-qJDD(2) * pkin(3) - qJ(4) * t491 + t446 + t519) * MDP(18) + (t518 - t544) * MDP(24) + (t517 - t551) * MDP(25) + (t552 + t553) * MDP(26) + (-t321 * qJD(2) - t291 * t432 + t292 * t435 + t495) * MDP(27) + (MDP(13) + MDP(15)) * (-t424 * t440 - t439) + (-MDP(11) + MDP(16)) * (qJDD(2) + t500) + t441 * qJD(5) + ((-MDP(18) * qJ(4) + MDP(12) - MDP(17)) * qJDD(1) + (t347 * MDP(14) + MDP(18) * t508 + t441) * qJD(1)) * t433; (t454 + t488) * MDP(18) + (-t317 * t435 + t318 * t432) * MDP(26) + (t418 + t451) * MDP(27) - t522 * MDP(17) * t440 + t467 * t354 + (MDP(15) * t433 - MDP(16) * t436) * qJDD(1) + t442 * qJD(5) + ((0.2e1 * qJD(2) * MDP(15) - t345 * MDP(18) - t355 * MDP(24) - t356 * MDP(25) + t321 * MDP(27)) * t436 + (t336 * MDP(18) + (-MDP(18) * t564 + 0.2e1 * MDP(16)) * qJD(2) + t442) * t433) * qJD(1); t356 * t355 * MDP(19) + (-t353 + t565) * MDP(20) + (t317 - t550) * MDP(21) + (t318 - t547) * MDP(22) + t354 * MDP(23) + (t303 * t378 + t339 * t356 + (t481 - t416) * t432 + t474 + t569) * MDP(24) + (g(1) * t351 - g(2) * t349 - g(3) * t535 + t302 * t378 - t339 * t355 + t457) * MDP(25) + (-pkin(5) * t317 + t355 * t531) * MDP(26) + (t531 * t300 + (-g(3) * t543 + t321 * t356 + t291 + t569) * pkin(5)) * MDP(27); (-t353 - t565) * MDP(26) + (-t295 * t356 - t300 * t355 + t334 + t372 + t426 - t558 + (t452 - t506) * t436 + t568) * MDP(27);];
tau  = t1;

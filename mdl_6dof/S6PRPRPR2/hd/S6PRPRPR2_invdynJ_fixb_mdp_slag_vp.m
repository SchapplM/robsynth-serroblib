% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:33:06
% EndTime: 2019-03-08 19:33:13
% DurationCPUTime: 6.66s
% Computational Cost: add. (3257->485), mult. (7686->684), div. (0->0), fcn. (6475->16), ass. (0->211)
t474 = cos(pkin(11));
t479 = sin(qJ(2));
t472 = sin(pkin(6));
t550 = qJD(1) * t472;
t532 = t479 * t550;
t436 = t474 * t532;
t470 = sin(pkin(11));
t482 = cos(qJ(2));
t531 = t482 * t550;
t391 = t470 * t531 + t436;
t478 = sin(qJ(4));
t481 = cos(qJ(4));
t517 = pkin(4) * t478 - qJ(5) * t481;
t409 = qJD(4) * t517 - qJD(5) * t478;
t596 = -t391 + t409;
t547 = qJD(2) * t481;
t453 = -qJD(6) + t547;
t469 = sin(pkin(12));
t473 = cos(pkin(12));
t541 = t473 * qJD(4);
t548 = qJD(2) * t478;
t418 = t469 * t548 - t541;
t546 = qJD(4) * t469;
t420 = t473 * t548 + t546;
t477 = sin(qJ(6));
t480 = cos(qJ(6));
t507 = t418 * t477 - t420 * t480;
t595 = t453 * t507;
t549 = qJD(2) * t472;
t527 = qJD(1) * t549;
t539 = qJDD(1) * t472;
t594 = t479 * t539 + t482 * t527;
t435 = t470 * t532;
t394 = t474 * t531 - t435;
t457 = pkin(2) * t470 + pkin(8);
t545 = qJD(4) * t478;
t529 = t457 * t545;
t574 = t469 * t481;
t556 = t394 * t574 + t469 * t529 + t473 * t596;
t432 = qJD(2) * pkin(2) + t531;
t382 = t470 * t432 + t436;
t380 = qJD(2) * pkin(8) + t382;
t476 = cos(pkin(6));
t452 = qJD(1) * t476 + qJD(3);
t593 = -t478 * t380 + t452 * t481;
t566 = t473 * t481;
t592 = -t394 * t566 + t469 * t596;
t471 = sin(pkin(10));
t475 = cos(pkin(10));
t562 = t476 * t482;
t591 = -t471 * t562 - t475 * t479;
t590 = -qJDD(4) * pkin(4) + qJDD(5);
t561 = t482 * t474;
t424 = t470 * t479 - t561;
t499 = t476 * t424;
t506 = t470 * t482 + t474 * t479;
t359 = -t471 * t506 - t475 * t499;
t362 = t471 * t499 - t475 * t506;
t570 = t472 * t479;
t401 = t470 * t570 - t472 * t561;
t589 = -g(1) * t362 - g(2) * t359 + g(3) * t401;
t563 = t476 * t479;
t552 = -t470 * t562 - t474 * t563;
t363 = -t475 * t424 + t471 * t552;
t358 = t471 * t424 + t475 * t552;
t588 = pkin(2) * t474;
t584 = pkin(9) + qJ(5);
t580 = t420 * t477;
t368 = t480 * t418 + t580;
t582 = t368 * t453;
t581 = t394 * t478;
t450 = t476 * qJDD(1) + qJDD(3);
t579 = t450 * t478;
t466 = pkin(12) + qJ(6);
t463 = sin(t466);
t577 = t463 * t481;
t464 = cos(t466);
t576 = t464 * t481;
t575 = t469 * t477;
t572 = t471 * t479;
t571 = t472 * t478;
t569 = t472 * t481;
t568 = t472 * t482;
t567 = t473 * t478;
t560 = qJDD(1) - g(3);
t449 = t482 * t539;
t397 = qJDD(2) * pkin(2) - t479 * t527 + t449;
t346 = t470 * t397 + t474 * t594;
t340 = qJDD(2) * pkin(8) + t346;
t309 = qJDD(4) * qJ(5) + t340 * t481 + t579 + (qJD(5) + t593) * qJD(4);
t345 = t397 * t474 - t470 * t594;
t504 = pkin(4) * t481 + qJ(5) * t478 + pkin(3);
t321 = qJD(2) * t409 - qJDD(2) * t504 - t345;
t306 = t473 * t309 + t469 * t321;
t505 = pkin(5) * t478 - pkin(9) * t566;
t559 = -qJD(4) * t505 - t556;
t558 = (-pkin(9) * t574 - t457 * t567) * qJD(4) + t592;
t352 = t481 * t380 + t478 * t452;
t348 = qJD(4) * qJ(5) + t352;
t381 = t432 * t474 - t435;
t357 = -qJD(2) * t504 - t381;
t314 = t473 * t348 + t469 * t357;
t425 = t469 * t480 + t473 * t477;
t498 = t425 * t481;
t542 = qJD(6) * t480;
t543 = qJD(6) * t478;
t356 = qJD(4) * t498 + t542 * t567 - t543 * t575;
t403 = t425 * t478;
t536 = t481 * qJDD(2);
t540 = qJD(2) * qJD(4);
t495 = t478 * t540 - t536;
t422 = qJDD(6) + t495;
t557 = t356 * t453 - t403 * t422;
t428 = t517 * qJD(2);
t325 = t469 * t428 + t473 * t593;
t555 = -t473 * t529 + t592;
t423 = -t480 * t473 + t575;
t497 = t423 * t481;
t554 = qJD(2) * t497 - t423 * qJD(6);
t553 = -qJD(2) * t498 + t425 * qJD(6);
t416 = -t504 - t588;
t376 = t469 * t416 + t457 * t566;
t467 = t478 ^ 2;
t551 = -t481 ^ 2 + t467;
t544 = qJD(4) * t481;
t538 = qJDD(2) * t478;
t537 = qJDD(4) * t469;
t534 = t475 * t562;
t460 = t473 * qJDD(4);
t526 = t481 * t540;
t496 = t526 + t538;
t387 = t469 * t496 - t460;
t388 = t473 * t496 + t537;
t533 = -t477 * t387 + t480 * t388 - t418 * t542;
t530 = t469 * t547;
t528 = qJ(5) * t536;
t525 = t469 * t538;
t524 = t473 * t538;
t522 = pkin(5) * t469 + t457;
t305 = -t309 * t469 + t473 * t321;
t301 = pkin(5) * t495 - pkin(9) * t388 + t305;
t304 = -pkin(9) * t387 + t306;
t521 = t480 * t301 - t477 * t304;
t313 = -t348 * t469 + t473 * t357;
t324 = t473 * t428 - t469 * t593;
t520 = t480 * t387 + t477 * t388;
t519 = t522 * t544 - t581;
t516 = t477 * t301 + t480 * t304;
t515 = -t305 * t469 + t306 * t473;
t311 = -pkin(5) * t547 - pkin(9) * t420 + t313;
t312 = -pkin(9) * t418 + t314;
t302 = t311 * t480 - t312 * t477;
t303 = t311 * t477 + t312 * t480;
t514 = -t313 * t469 + t314 * t473;
t402 = t506 * t472;
t378 = t402 * t481 + t476 * t478;
t329 = -t378 * t469 + t401 * t473;
t330 = t378 * t473 + t401 * t469;
t513 = t329 * t480 - t330 * t477;
t512 = t329 * t477 + t330 * t480;
t406 = t473 * t416;
t354 = -pkin(9) * t567 + t406 + (-t457 * t469 - pkin(5)) * t481;
t365 = -pkin(9) * t469 * t478 + t376;
t511 = t354 * t480 - t365 * t477;
t510 = t354 * t477 + t365 * t480;
t355 = -qJD(4) * t497 - t425 * t543;
t404 = t423 * t478;
t509 = t355 * t453 + t404 * t422;
t377 = t402 * t478 - t476 * t481;
t503 = t478 * t340 + t380 * t544 - t450 * t481 + t452 * t545;
t442 = t584 * t473;
t502 = qJD(2) * t505 + qJD(5) * t469 + qJD(6) * t442 + t324;
t441 = t584 * t469;
t501 = pkin(9) * t530 + qJD(5) * t473 - qJD(6) * t441 - t325;
t500 = -g(3) * t476 + (-g(1) * t471 + g(2) * t475) * t472;
t315 = -qJD(6) * t580 + t533;
t494 = g(1) * (t363 * t478 - t471 * t569) + g(2) * (-t358 * t478 + t475 * t569) + g(3) * t377;
t335 = -t358 * t481 - t475 * t571;
t337 = t363 * t481 + t471 * t571;
t493 = g(1) * t337 + g(2) * t335 + g(3) * t378;
t492 = -g(1) * t363 + g(2) * t358 - g(3) * t402;
t344 = -qJD(4) * pkin(4) + qJD(5) - t593;
t310 = t503 + t590;
t490 = -t310 + t494;
t489 = -qJ(5) * t545 + (qJD(5) - t344) * t481;
t379 = -qJD(2) * pkin(3) - t381;
t458 = -pkin(3) - t588;
t488 = -qJDD(4) * t457 + (qJD(2) * t458 + t379 + t394) * qJD(4);
t316 = -qJD(6) * t507 + t520;
t487 = -g(1) * t591 - g(3) * t568;
t486 = t494 - t503;
t483 = qJD(4) ^ 2;
t485 = -qJD(2) * t391 + t457 * t483 - t345 - t589 + (-pkin(3) + t458) * qJDD(2);
t484 = qJD(2) ^ 2;
t459 = -pkin(5) * t473 - pkin(4);
t443 = pkin(2) * t534;
t439 = qJDD(4) * t481 - t478 * t483;
t438 = qJDD(4) * t478 + t481 * t483;
t408 = t522 * t478;
t393 = t424 * t549;
t392 = qJD(2) * t402;
t375 = -t457 * t574 + t406;
t364 = t507 * t545;
t338 = pkin(5) * t530 + t352;
t328 = -qJD(4) * t377 - t393 * t481;
t327 = qJD(4) * t378 - t393 * t478;
t326 = pkin(5) * t418 + t344;
t318 = t328 * t473 + t392 * t469;
t317 = -t328 * t469 + t392 * t473;
t307 = pkin(5) * t387 + t310;
t1 = [t560 * MDP(1) + (-t345 * t401 + t346 * t402 - t381 * t392 - t382 * t393 + t450 * t476 - g(3)) * MDP(5) + (-qJD(4) * t327 - qJDD(4) * t377 - t401 * t536) * MDP(11) + (-qJD(4) * t328 - qJDD(4) * t378 + t401 * t538) * MDP(12) + (t327 * t418 - t329 * t536 + t377 * t387) * MDP(13) + (t327 * t420 + t330 * t536 + t377 * t388) * MDP(14) + (-t317 * t420 - t318 * t418 - t329 * t388 - t330 * t387) * MDP(15) + (t305 * t329 + t306 * t330 + t310 * t377 + t313 * t317 + t314 * t318 + t327 * t344 - g(3)) * MDP(16) + (-(-qJD(6) * t512 + t317 * t480 - t318 * t477) * t453 + t513 * t422 + t327 * t368 + t377 * t316) * MDP(22) + ((qJD(6) * t513 + t317 * t477 + t318 * t480) * t453 - t512 * t422 - t327 * t507 + t377 * t315) * MDP(23) + ((qJDD(2) * t482 - t479 * t484) * MDP(3) + (-qJDD(2) * t479 - t482 * t484) * MDP(4)) * t472 + ((-t392 * t481 + t401 * t545) * MDP(11) + (t392 * t478 + t401 * t544) * MDP(12) + (-t317 * t481 + t329 * t545) * MDP(13) + (t318 * t481 - t330 * t545) * MDP(14)) * qJD(2); qJDD(2) * MDP(2) + (t449 - g(2) * (t534 - t572) + t487) * MDP(3) + (-g(1) * (t471 * t563 - t475 * t482) - g(2) * (-t471 * t482 - t475 * t563) - t560 * t570) * MDP(4) + (-g(2) * t443 + t381 * t391 - t382 * t394 + (g(2) * t572 + t345 * t474 + t346 * t470 + t487) * pkin(2)) * MDP(5) + (qJDD(2) * t467 + 0.2e1 * t478 * t526) * MDP(6) + 0.2e1 * (t478 * t536 - t540 * t551) * MDP(7) + t438 * MDP(8) + t439 * MDP(9) + (t478 * t488 - t481 * t485) * MDP(11) + (t478 * t485 + t481 * t488) * MDP(12) + (t492 * t469 + (t310 * t469 + t457 * t387 - t394 * t418 + (qJD(2) * t375 + t313) * qJD(4)) * t478 + (-t375 * qJDD(2) - t305 + (t344 * t469 + t418 * t457) * qJD(4) - t556 * qJD(2) + t589 * t473) * t481) * MDP(13) + (t492 * t473 + (t310 * t473 + t457 * t388 - t394 * t420 + (-qJD(2) * t376 - t314) * qJD(4)) * t478 + (t376 * qJDD(2) + t306 + (t344 * t473 + t420 * t457) * qJD(4) + t555 * qJD(2) - t589 * t469) * t481) * MDP(14) + (-t375 * t388 - t376 * t387 - t556 * t420 - t555 * t418 + (-t313 * t473 - t314 * t469) * t544 + (-t305 * t473 - t306 * t469 + t589) * t478) * MDP(15) + (t306 * t376 + t305 * t375 + t310 * t478 * t457 - g(1) * (pkin(2) * t591 + pkin(8) * t363) - g(2) * (-pkin(2) * t572 - pkin(8) * t358 + t443) - g(3) * (pkin(2) * t568 + pkin(8) * t402) + (t457 * t544 - t581) * t344 + t555 * t314 + t556 * t313 + t589 * t504) * MDP(16) + (-t315 * t404 - t355 * t507) * MDP(17) + (-t315 * t403 + t316 * t404 - t355 * t368 + t356 * t507) * MDP(18) + (-t315 * t481 - t364 - t509) * MDP(19) + (t316 * t481 - t368 * t545 + t557) * MDP(20) + (-t422 * t481 - t453 * t545) * MDP(21) + (t511 * t422 - t521 * t481 + t302 * t545 + t408 * t316 + t307 * t403 + t326 * t356 - g(1) * (t362 * t576 + t363 * t463) - g(2) * (-t358 * t463 + t359 * t576) - g(3) * (-t401 * t576 + t402 * t463) + (t477 * t558 + t480 * t559) * t453 + t519 * t368 + (t303 * t481 + t453 * t510) * qJD(6)) * MDP(22) + (-t510 * t422 + t516 * t481 - t303 * t545 + t408 * t315 - t307 * t404 + t326 * t355 - g(1) * (-t362 * t577 + t363 * t464) - g(2) * (-t358 * t464 - t359 * t577) - g(3) * (t401 * t577 + t402 * t464) + (-t477 * t559 + t480 * t558) * t453 - t519 * t507 + (t302 * t481 + t453 * t511) * qJD(6)) * MDP(23); (t500 + t450) * MDP(5) + t439 * MDP(11) - t438 * MDP(12) + t500 * MDP(16) + t557 * MDP(22) + (-t364 + t509) * MDP(23) + ((-t387 * t473 + t388 * t469) * MDP(15) + t515 * MDP(16)) * t478 + ((-t387 + t525) * MDP(13) + (-t388 + t524) * MDP(14) - t310 * MDP(16) - t316 * MDP(22) - t315 * MDP(23)) * t481 + ((MDP(13) * t418 + MDP(14) * t420 + MDP(16) * t344 + MDP(22) * t368) * t478 + ((-t418 * t473 + t420 * t469) * MDP(15) + t514 * MDP(16)) * t481 + (-MDP(13) * t469 - MDP(14) * t473) * qJD(2) * t551) * qJD(4); MDP(8) * t538 + MDP(9) * t536 + qJDD(4) * MDP(10) + (qJD(4) * t352 - t379 * t548 + t486) * MDP(11) + (-t579 + (-qJD(2) * t379 - t340) * t481 + t493) * MDP(12) + (t469 * t528 - pkin(4) * t387 - t352 * t418 + t490 * t473 + (-t313 * t478 + t324 * t481 + t469 * t489) * qJD(2)) * MDP(13) + (t473 * t528 - pkin(4) * t388 - t352 * t420 - t490 * t469 + (t314 * t478 - t325 * t481 + t473 * t489) * qJD(2)) * MDP(14) + (t324 * t420 + t325 * t418 + (-qJ(5) * t387 - qJD(5) * t418 + t313 * t547 + t306) * t473 + (qJ(5) * t388 + qJD(5) * t420 + t314 * t547 - t305) * t469 - t493) * MDP(15) + (-t313 * t324 - t314 * t325 - t344 * t352 + t514 * qJD(5) + t490 * pkin(4) + (-t493 + t515) * qJ(5)) * MDP(16) + (t315 * t425 - t507 * t554) * MDP(17) + (-t315 * t423 - t316 * t425 - t368 * t554 + t507 * t553) * MDP(18) + (t422 * t425 - t453 * t554 + t507 * t548) * MDP(19) + (t368 * t548 - t422 * t423 + t453 * t553) * MDP(20) + t453 * MDP(21) * t548 + ((-t441 * t480 - t442 * t477) * t422 + t459 * t316 + t307 * t423 - t302 * t548 - t338 * t368 + (t477 * t501 + t480 * t502) * t453 + t553 * t326 + t494 * t464) * MDP(22) + (-(-t441 * t477 + t442 * t480) * t422 + t459 * t315 + t307 * t425 + t303 * t548 + t338 * t507 + (-t477 * t502 + t480 * t501) * t453 + t554 * t326 - t494 * t463) * MDP(23) + (-MDP(6) * t478 * t481 + MDP(7) * t551) * t484; (t525 - t460 + (-t420 + t546) * t547) * MDP(13) + (t524 + t537 + (t418 + t541) * t547) * MDP(14) + (-t418 ^ 2 - t420 ^ 2) * MDP(15) + (t313 * t420 + t314 * t418 - t486 + t590) * MDP(16) + (t316 + t595) * MDP(22) + (t315 + t582) * MDP(23); -t507 * t368 * MDP(17) + (-t368 ^ 2 + t507 ^ 2) * MDP(18) + (t533 - t582) * MDP(19) + (-t520 + t595) * MDP(20) + t422 * MDP(21) + (-t303 * t453 + t326 * t507 - g(1) * (-t337 * t463 - t362 * t464) - g(2) * (-t335 * t463 - t359 * t464) - g(3) * (-t378 * t463 + t401 * t464) + t521) * MDP(22) + (-t302 * t453 + t326 * t368 - g(1) * (-t337 * t464 + t362 * t463) - g(2) * (-t335 * t464 + t359 * t463) - g(3) * (-t378 * t464 - t401 * t463) - t516) * MDP(23) + (-MDP(19) * t580 + MDP(20) * t507 - MDP(22) * t303 - MDP(23) * t302) * qJD(6);];
tau  = t1;

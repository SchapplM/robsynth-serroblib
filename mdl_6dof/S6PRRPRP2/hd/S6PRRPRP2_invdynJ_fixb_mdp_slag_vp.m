% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:32:29
% EndTime: 2019-03-08 21:32:39
% DurationCPUTime: 7.84s
% Computational Cost: add. (5593->547), mult. (12955->721), div. (0->0), fcn. (10207->14), ass. (0->241)
t512 = sin(qJ(3));
t515 = cos(qJ(3));
t633 = qJ(4) + pkin(8);
t566 = qJD(3) * t633;
t455 = qJD(4) * t515 - t512 * t566;
t456 = -qJD(4) * t512 - t515 * t566;
t505 = sin(pkin(11));
t630 = cos(pkin(11));
t393 = t455 * t630 + t505 * t456;
t564 = t630 * t515;
t542 = -t505 * t512 + t564;
t516 = cos(qJ(2));
t507 = sin(pkin(6));
t593 = qJD(1) * t507;
t573 = t516 * t593;
t429 = t542 * t573;
t649 = t393 - t429;
t487 = qJD(2) * t564;
t591 = qJD(2) * t512;
t462 = -t505 * t591 + t487;
t454 = qJD(5) - t462;
t471 = t505 * t515 + t512 * t630;
t601 = t455 * t505 - t630 * t456 - t471 * t573;
t463 = t471 * qJD(3);
t466 = t542 * qJD(3);
t631 = qJD(3) * pkin(3);
t579 = t512 * t631;
t395 = pkin(4) * t463 - pkin(9) * t466 + t579;
t496 = t515 * pkin(3) + pkin(2);
t408 = -pkin(4) * t542 - pkin(9) * t471 - t496;
t477 = t633 * t512;
t478 = t633 * t515;
t431 = -t505 * t477 + t478 * t630;
t511 = sin(qJ(5));
t514 = cos(qJ(5));
t513 = sin(qJ(2));
t574 = t513 * t593;
t587 = qJD(5) * t514;
t588 = qJD(5) * t511;
t648 = -t408 * t587 + t431 * t588 - t649 * t514 + (-t395 + t574) * t511;
t598 = t511 * t408 + t514 * t431;
t615 = t507 * t513;
t577 = t512 * t615;
t509 = cos(pkin(6));
t611 = t509 * t515;
t467 = -t577 + t611;
t614 = t507 * t515;
t468 = t509 * t512 + t513 * t614;
t400 = t505 * t467 + t468 * t630;
t613 = t507 * t516;
t484 = t511 * t613;
t384 = t400 * t514 - t484;
t506 = sin(pkin(10));
t508 = cos(pkin(10));
t610 = t509 * t516;
t458 = t506 * t513 - t508 * t610;
t460 = t506 * t610 + t508 * t513;
t557 = g(1) * t460 + g(2) * t458;
t647 = -g(3) * t613 + t557;
t646 = MDP(19) + MDP(21);
t645 = -MDP(20) + MDP(23);
t584 = qJD(2) * qJD(3);
t569 = t512 * t584;
t644 = pkin(3) * t569 + qJDD(4);
t612 = t509 * t513;
t459 = t506 * t516 + t508 * t612;
t461 = -t506 * t612 + t508 * t516;
t502 = qJ(3) + pkin(11);
t497 = sin(t502);
t498 = cos(t502);
t616 = t507 * t508;
t617 = t506 * t507;
t536 = -g(3) * (-t497 * t615 + t498 * t509) - g(2) * (-t459 * t497 - t498 * t616) - g(1) * (-t461 * t497 + t498 * t617);
t560 = t633 * qJD(2) + t574;
t592 = qJD(1) * t509;
t427 = t512 * t592 + t515 * t560;
t417 = t505 * t427;
t426 = -t512 * t560 + t515 * t592;
t421 = t426 + t631;
t368 = t421 * t630 - t417;
t363 = -qJD(3) * pkin(4) - t368;
t464 = t471 * qJD(2);
t432 = -t514 * qJD(3) + t464 * t511;
t434 = qJD(3) * t511 + t464 * t514;
t346 = t432 * pkin(5) - t434 * qJ(6) + t363;
t582 = qJDD(2) * t512;
t551 = qJDD(2) * t564 - t505 * t582;
t406 = qJD(2) * t463 + qJDD(5) - t551;
t639 = pkin(3) * t505;
t493 = pkin(9) + t639;
t620 = t493 * t406;
t643 = t346 * t454 - t620;
t567 = qJDD(1) * t613;
t585 = qJD(1) * qJD(2);
t570 = t513 * t585;
t552 = t507 * t570 - t567;
t628 = qJDD(2) * pkin(2);
t441 = t552 - t628;
t517 = qJD(3) ^ 2;
t642 = -pkin(8) * t517 + t507 * (-g(3) * t516 + t570) - t441 + t557 + t628;
t641 = qJD(3) * t487 + qJDD(2) * t471 - t505 * t569;
t640 = t434 ^ 2;
t638 = pkin(3) * t512;
t637 = pkin(5) * t406;
t632 = qJD(2) * pkin(2);
t629 = qJ(6) * t406;
t565 = t630 * t427;
t369 = t505 * t421 + t565;
t364 = qJD(3) * pkin(9) + t369;
t453 = -qJD(2) * t496 + qJD(4) - t573;
t382 = -pkin(4) * t462 - pkin(9) * t464 + t453;
t344 = t364 * t514 + t382 * t511;
t627 = t344 * t454;
t583 = qJD(3) * qJD(5);
t365 = -t511 * qJDD(3) + t464 * t588 + (-t641 - t583) * t514;
t626 = t365 * t511;
t624 = t432 * t462;
t623 = t434 * t432;
t561 = t434 * t454;
t622 = t454 * t511;
t621 = t471 * t514;
t619 = t498 * t511;
t618 = t498 * t514;
t609 = t633 * t513;
t401 = t511 * t406;
t402 = t514 * t406;
t608 = t514 * t516;
t607 = qJDD(1) - g(3);
t606 = qJ(6) * t463 - qJD(6) * t542 - t648;
t396 = t429 * t511 - t514 * t574;
t605 = -pkin(5) * t463 + qJD(5) * t598 + t393 * t511 - t395 * t514 - t396;
t553 = pkin(5) * t511 - qJ(6) * t514;
t554 = t514 * pkin(5) + t511 * qJ(6);
t604 = t553 * t466 + (qJD(5) * t554 - qJD(6) * t514) * t471 + t601;
t580 = t509 * qJDD(1);
t486 = t515 * t580;
t442 = qJDD(2) * pkin(8) + (qJDD(1) * t513 + t516 * t585) * t507;
t529 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t592 + t442;
t549 = t560 * qJD(3);
t359 = qJDD(3) * pkin(3) - t512 * t529 - t515 * t549 + t486;
t360 = (-t549 + t580) * t512 + t529 * t515;
t336 = t505 * t359 + t630 * t360;
t519 = -t514 * qJDD(3) + t511 * t641;
t366 = t434 * qJD(5) + t519;
t603 = -t511 * t366 - t432 * t587;
t373 = t426 * t630 - t417;
t394 = pkin(3) * t591 + pkin(4) * t464 - pkin(9) * t462;
t602 = t514 * t373 + t511 * t394;
t600 = t454 * t587 + t401;
t599 = t462 * t622 + t402;
t597 = -t458 * t496 + t459 * t633;
t596 = -t460 * t496 + t461 * t633;
t371 = t426 * t505 + t565;
t595 = -qJD(6) * t511 + t454 * t553 - t371;
t503 = t512 ^ 2;
t594 = -t515 ^ 2 + t503;
t590 = qJD(2) * t513;
t589 = qJD(5) * t493;
t343 = -t364 * t511 + t382 * t514;
t586 = qJD(6) - t343;
t581 = qJDD(2) * t515;
t578 = t508 * t614;
t576 = t507 * t608;
t575 = t630 * pkin(3);
t572 = t507 * t590;
t571 = qJD(2) * t613;
t568 = t516 * t584;
t334 = qJDD(3) * pkin(9) + t336;
t410 = -qJD(3) * t464 + t551;
t524 = -pkin(3) * t581 + t441 + t644;
t351 = -t410 * pkin(4) - pkin(9) * t641 + t524;
t541 = t514 * t334 + t511 * t351 - t364 * t588 + t382 * t587;
t326 = qJD(6) * t454 + t541 + t629;
t337 = -pkin(5) * t454 + t586;
t563 = -t337 * t462 + t326;
t559 = t511 * t334 - t514 * t351 + t364 * t587 + t382 * t588;
t327 = qJDD(6) + t559 - t637;
t338 = qJ(6) * t454 + t344;
t562 = t338 * t462 + t327;
t430 = t630 * t477 + t478 * t505;
t558 = pkin(4) * t498 + pkin(9) * t497;
t556 = g(1) * t461 + g(2) * t459;
t555 = g(1) * t506 - g(2) * t508;
t335 = t359 * t630 - t505 * t360;
t550 = t337 * t514 - t338 * t511;
t518 = qJD(2) ^ 2;
t548 = qJDD(2) * t516 - t513 * t518;
t547 = pkin(4) + t554;
t383 = t400 * t511 + t576;
t544 = t466 * t511 + t471 * t587;
t543 = -t466 * t514 + t471 * t588;
t539 = t363 * t454 - t620;
t388 = -t458 * t619 - t459 * t514;
t390 = -t460 * t619 - t461 * t514;
t436 = t484 * t498 - t514 * t615;
t538 = g(1) * t390 + g(2) * t388 + g(3) * t436;
t389 = -t458 * t618 + t459 * t511;
t391 = -t460 * t618 + t461 * t511;
t437 = (t498 * t608 + t511 * t513) * t507;
t537 = -g(1) * t391 - g(2) * t389 - g(3) * t437;
t333 = -qJDD(3) * pkin(4) - t335;
t533 = -g(3) * t615 - t556;
t475 = -t573 - t632;
t532 = -qJD(2) * t475 - t442 + t556;
t528 = -t454 * t589 + t536;
t412 = t459 * t498 - t497 * t616;
t377 = t412 * t511 - t458 * t514;
t414 = t461 * t498 + t497 * t617;
t379 = t414 * t511 - t460 * t514;
t445 = t497 * t509 + t498 * t615;
t415 = t445 * t511 + t576;
t527 = g(1) * t379 + g(2) * t377 + g(3) * t415 - t559;
t328 = t366 * pkin(5) + t365 * qJ(6) - t434 * qJD(6) + t333;
t526 = -t328 + t528;
t525 = -pkin(8) * qJDD(3) + (t475 + t573 - t632) * qJD(3);
t378 = t412 * t514 + t458 * t511;
t380 = t414 * t514 + t460 * t511;
t416 = t445 * t514 - t484;
t523 = -g(1) * t380 - g(2) * t378 - g(3) * t416 + t541;
t522 = t346 * t434 + qJDD(6) - t527;
t494 = -t575 - pkin(4);
t489 = pkin(3) * t611;
t476 = t506 * pkin(3) * t614;
t473 = t496 * t613;
t469 = -t575 - t547;
t425 = -qJD(3) * t468 - t512 * t571;
t424 = qJD(3) * t467 + t515 * t571;
t409 = -qJDD(2) * t496 + t552 + t644;
t399 = -t467 * t630 + t468 * t505;
t381 = pkin(5) * t434 + qJ(6) * t432;
t374 = t471 * t553 + t430;
t372 = t424 * t630 + t505 * t425;
t370 = t424 * t505 - t425 * t630;
t355 = pkin(5) * t542 - t408 * t514 + t431 * t511;
t354 = -qJ(6) * t542 + t598;
t347 = t432 * t454 - t365;
t345 = -pkin(5) * t464 + t373 * t511 - t394 * t514;
t342 = qJ(6) * t464 + t602;
t341 = qJD(5) * t384 + t372 * t511 - t514 * t572;
t340 = -qJD(5) * t383 + t372 * t514 + t511 * t572;
t1 = [t607 * MDP(1) + (qJD(3) * t425 + qJDD(3) * t467) * MDP(10) + (-qJD(3) * t424 - qJDD(3) * t468) * MDP(11) + (t370 * t464 + t372 * t462 + t399 * t641 + t400 * t410) * MDP(12) + (-t335 * t399 + t336 * t400 - t368 * t370 + t369 * t372 - g(3)) * MDP(13) + (-t340 * t432 + t341 * t434 - t365 * t383 - t366 * t384) * MDP(22) + (t326 * t384 + t327 * t383 + t328 * t399 + t337 * t341 + t338 * t340 + t346 * t370 - g(3)) * MDP(24) + t646 * (-t341 * t454 + t399 * t366 + t370 * t432 - t383 * t406) + t645 * (t340 * t454 + t365 * t399 - t370 * t434 + t384 * t406) + (t548 * MDP(3) + (-qJDD(2) * t513 - t516 * t518) * MDP(4) + (-t512 * t568 + t515 * t548) * MDP(10) + (-t512 * t548 - t515 * t568) * MDP(11) + (-t409 * t516 + t453 * t590) * MDP(13)) * t507; qJDD(2) * MDP(2) + (t647 + t567) * MDP(3) + (-t607 * t615 + t556) * MDP(4) + (qJDD(2) * t503 + 0.2e1 * t515 * t569) * MDP(5) + 0.2e1 * (t512 * t581 - t584 * t594) * MDP(6) + (qJDD(3) * t512 + t515 * t517) * MDP(7) + (qJDD(3) * t515 - t512 * t517) * MDP(8) + (t525 * t512 + t515 * t642) * MDP(10) + (-t512 * t642 + t525 * t515) * MDP(11) + (-t335 * t471 + t336 * t542 - t368 * t466 - t369 * t463 + t431 * t410 + t430 * t641 + t649 * t462 + t464 * t601 + t533) * MDP(12) + (t336 * t431 - t335 * t430 - t409 * t496 - g(1) * t596 - g(2) * t597 - g(3) * (t507 * t609 + t473) + (-t574 + t579) * t453 + t649 * t369 - t601 * t368) * MDP(13) + (-t365 * t621 - t434 * t543) * MDP(14) + ((-t432 * t514 - t434 * t511) * t466 + (t626 - t366 * t514 + (t432 * t511 - t434 * t514) * qJD(5)) * t471) * MDP(15) + (t365 * t542 + t402 * t471 + t434 * t463 - t454 * t543) * MDP(16) + (t366 * t542 - t401 * t471 - t432 * t463 - t454 * t544) * MDP(17) + (-t406 * t542 + t454 * t463) * MDP(18) + (t559 * t542 + t343 * t463 + t430 * t366 + t396 * t454 + t601 * t432 + ((-qJD(5) * t431 + t395) * t454 + t408 * t406 + t363 * qJD(5) * t471) * t514 + ((-qJD(5) * t408 - t393) * t454 - t431 * t406 + t333 * t471 + t363 * t466) * t511 + t537) * MDP(19) + (t333 * t621 - t344 * t463 - t543 * t363 - t430 * t365 - t598 * t406 + t601 * t434 + t454 * t648 + t541 * t542 + t538) * MDP(20) + (t328 * t471 * t511 + t327 * t542 - t337 * t463 + t346 * t544 - t355 * t406 + t366 * t374 + t432 * t604 - t454 * t605 + t537) * MDP(21) + (-t354 * t366 - t355 * t365 + t550 * t466 + t605 * t434 - t606 * t432 + t647 * t497 + (-t326 * t511 + t327 * t514 + (-t337 * t511 - t338 * t514) * qJD(5)) * t471) * MDP(22) + (-t326 * t542 - t328 * t621 + t338 * t463 + t346 * t543 + t354 * t406 + t365 * t374 - t434 * t604 + t454 * t606 - t538) * MDP(23) + (t326 * t354 + t328 * t374 + t327 * t355 - g(1) * (pkin(5) * t391 + qJ(6) * t390 - t460 * t558 + t596) - g(2) * (pkin(5) * t389 + qJ(6) * t388 - t458 * t558 + t597) + t604 * t346 + t606 * t338 + t605 * t337 + (-pkin(5) * t437 - qJ(6) * t436 - t473 - (t516 * t558 + t609) * t507) * g(3)) * MDP(24); MDP(7) * t582 + MDP(8) * t581 + qJDD(3) * MDP(9) + (-g(3) * t467 + t512 * t532 - t555 * t614 + t486) * MDP(10) + (g(3) * t468 + (t507 * t555 - t580) * t512 + t532 * t515) * MDP(11) + (t410 * t639 - t641 * t575 - (-t369 + t371) * t464 + (t368 - t373) * t462) * MDP(12) + (-g(1) * t476 - g(3) * t489 + t368 * t371 - t369 * t373 + (g(2) * t578 + t335 * t630 + t336 * t505 + (-qJD(2) * t453 - t533) * t512) * pkin(3)) * MDP(13) + (t514 * t561 - t626) * MDP(14) + ((-t365 + t624) * t514 - t434 * t622 + t603) * MDP(15) + (-t454 * t462 * t514 - t434 * t464 + t600) * MDP(16) + (t432 * t464 - t454 * t588 + t599) * MDP(17) - t454 * t464 * MDP(18) + (-t343 * t464 + t494 * t366 - t371 * t432 + (t373 * t454 + t539) * t511 + (-t333 + (-t394 - t589) * t454 + t536) * t514) * MDP(19) + (-t494 * t365 + t602 * t454 + t344 * t464 - t371 * t434 + t539 * t514 + (t333 - t528) * t511) * MDP(20) + (t337 * t464 + t345 * t454 + t366 * t469 + t595 * t432 + t511 * t643 + t526 * t514) * MDP(21) + (-g(1) * t414 - g(2) * t412 - g(3) * t445 + t342 * t432 - t345 * t434 + (-t366 * t493 + (t434 * t493 + t337) * qJD(5) + t563) * t514 + (-t365 * t493 + (t432 * t493 - t338) * qJD(5) + t562) * t511) * MDP(22) + (-t338 * t464 - t342 * t454 + t365 * t469 - t595 * t434 + t526 * t511 - t514 * t643) * MDP(23) + (t328 * t469 - t338 * t342 - t337 * t345 - g(1) * (pkin(9) * t414 - t461 * t638 + t476) - g(2) * (-pkin(3) * t578 + pkin(9) * t412 - t459 * t638) - g(3) * (-pkin(3) * t577 + pkin(9) * t445 + t489) + t595 * t346 + (qJD(5) * t550 + t326 * t514 + t327 * t511) * t493 + t536 * t547) * MDP(24) + (-MDP(5) * t512 * t515 + MDP(6) * t594) * t518; -t462 ^ 2 * MDP(12) + (-t369 * t462 + t524 - t647) * MDP(13) + t599 * MDP(19) + t603 * MDP(22) + t600 * MDP(23) - t647 * MDP(24) + (-MDP(12) * t464 + t368 * MDP(13) - t346 * MDP(24) - t432 * t646 + t434 * t645) * t464 + (t406 * MDP(21) + (t365 + t624) * MDP(22) + (qJD(5) * t338 - t562) * MDP(24) + (-MDP(20) * t454 - t462 * MDP(23)) * t454) * t514 + (-t406 * MDP(20) + (qJD(5) * t337 + t563) * MDP(24) + MDP(22) * t561 + (-qJD(5) * MDP(19) - t454 * MDP(21)) * t454) * t511; MDP(14) * t623 + (-t432 ^ 2 + t640) * MDP(15) + t347 * MDP(16) + (-t464 * t587 - t511 * t583 - t519 + t561) * MDP(17) + t406 * MDP(18) + (-t363 * t434 + t527 + t627) * MDP(19) + (t343 * t454 + t363 * t432 - t523) * MDP(20) + (-t381 * t432 - t522 + t627 + 0.2e1 * t637) * MDP(21) + (pkin(5) * t365 - qJ(6) * t366 + (t338 - t344) * t434 + (t337 - t586) * t432) * MDP(22) + (0.2e1 * t629 - t346 * t432 + t381 * t434 + (0.2e1 * qJD(6) - t343) * t454 + t523) * MDP(23) + (t326 * qJ(6) - t327 * pkin(5) - t346 * t381 - t337 * t344 - g(1) * (-pkin(5) * t379 + qJ(6) * t380) - g(2) * (-pkin(5) * t377 + qJ(6) * t378) - g(3) * (-pkin(5) * t415 + qJ(6) * t416) + t586 * t338) * MDP(24); (-qJDD(5) + t410 + t623) * MDP(21) + t347 * MDP(22) + (-t454 ^ 2 - t640) * MDP(23) + (-t338 * t454 + t522 - t637) * MDP(24);];
tau  = t1;

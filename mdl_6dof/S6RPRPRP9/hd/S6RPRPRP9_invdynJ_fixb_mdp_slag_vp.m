% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPRP9_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:29:41
% EndTime: 2019-03-09 03:29:52
% DurationCPUTime: 7.64s
% Computational Cost: add. (5114->563), mult. (9942->695), div. (0->0), fcn. (6582->10), ass. (0->226)
t485 = sin(pkin(9));
t618 = pkin(8) + qJ(4);
t452 = t618 * t485;
t486 = cos(pkin(9));
t453 = t618 * t486;
t488 = sin(qJ(5));
t491 = cos(qJ(5));
t397 = -t452 * t488 + t453 * t491;
t492 = cos(qJ(3));
t570 = qJD(1) * qJD(3);
t549 = t492 * t570;
t489 = sin(qJ(3));
t566 = qJDD(1) * t489;
t513 = t549 + t566;
t439 = qJDD(5) + t513;
t482 = pkin(9) + qJ(5);
t474 = sin(t482);
t620 = g(3) * t489;
t493 = cos(qJ(1));
t481 = g(2) * t493;
t490 = sin(qJ(1));
t621 = g(1) * t490;
t639 = t481 - t621;
t626 = -t492 * t639 - t620;
t642 = t397 * t439 - t474 * t626;
t579 = qJD(1) * t492;
t551 = t485 * t579;
t577 = qJD(3) * t486;
t436 = -t551 + t577;
t558 = t486 * t579;
t578 = qJD(3) * t485;
t437 = t558 + t578;
t382 = -t491 * t436 + t437 * t488;
t641 = t382 ^ 2;
t580 = qJD(1) * t489;
t467 = qJD(5) + t580;
t640 = t382 * t467;
t559 = t485 * t580;
t598 = t491 * t486;
t572 = qJD(5) * t491;
t573 = qJD(5) * t488;
t630 = -t485 * t573 + t486 * t572;
t586 = -t488 * t559 + t580 * t598 + t630;
t441 = t485 * t491 + t486 * t488;
t425 = t441 * qJD(1);
t427 = t441 * qJD(5);
t585 = t489 * t425 + t427;
t638 = -MDP(24) + MDP(27);
t520 = t436 * t488 + t437 * t491;
t624 = t520 ^ 2;
t637 = t467 * t520;
t636 = t639 * t486;
t530 = pkin(3) * t492 + qJ(4) * t489;
t443 = t530 * qJD(1);
t494 = -pkin(1) - pkin(7);
t461 = qJD(1) * t494 + qJD(2);
t606 = t485 * t492;
t393 = t486 * t443 - t461 * t606;
t563 = pkin(8) * t486 * t489;
t361 = (pkin(4) * t492 + t563) * qJD(1) + t393;
t605 = t486 * t492;
t394 = t485 * t443 + t461 * t605;
t373 = pkin(8) * t559 + t394;
t440 = t485 * t488 - t598;
t519 = -t452 * t491 - t453 * t488;
t635 = qJD(4) * t440 - qJD(5) * t519 + t488 * t361 + t491 * t373;
t634 = -qJD(4) * t441 - qJD(5) * t397 - t361 * t491 + t373 * t488;
t529 = pkin(3) * t489 - qJ(4) * t492;
t446 = qJ(2) + t529;
t434 = t486 * t446;
t548 = -t485 * t494 + pkin(4);
t380 = -pkin(8) * t605 + t489 * t548 + t434;
t600 = t489 * t494;
t402 = t485 * t446 + t486 * t600;
t392 = -pkin(8) * t606 + t402;
t633 = t488 * t380 + t491 * t392;
t599 = t490 * t492;
t544 = -g(1) * t599 + t620;
t459 = qJDD(1) * t494 + qJDD(2);
t546 = -t459 - t481;
t631 = t546 * t492 - t544;
t629 = MDP(23) + MDP(25);
t565 = qJDD(1) * t492;
t538 = qJDD(3) * t486 - t485 * t565;
t550 = t489 * t570;
t509 = t485 * t550 + t538;
t539 = t485 * qJDD(3) - t486 * t550;
t510 = t565 * t486 + t539;
t341 = qJD(5) * t520 + t488 * t510 - t491 * t509;
t419 = qJD(3) * t530 - qJD(4) * t492 + qJD(2);
t404 = t486 * t419;
t358 = t404 + (t492 * t548 + t563) * qJD(3);
t574 = qJD(3) * t494;
t554 = t492 * t574;
t389 = t485 * t419 + t486 * t554;
t576 = qJD(3) * t489;
t557 = t485 * t576;
t370 = pkin(8) * t557 + t389;
t625 = -qJD(5) * t633 + t358 * t491 - t370 * t488;
t623 = pkin(4) * t485;
t622 = pkin(5) * t439;
t619 = g(3) * t492;
t617 = pkin(1) * qJDD(1);
t496 = qJD(1) ^ 2;
t616 = qJ(2) * t496;
t615 = qJ(6) * t439;
t430 = t446 * qJD(1);
t445 = t489 * t461;
t431 = qJD(3) * qJ(4) + t445;
t374 = t486 * t430 - t431 * t485;
t350 = pkin(4) * t580 - pkin(8) * t437 + t374;
t375 = t485 * t430 + t486 * t431;
t354 = pkin(8) * t436 + t375;
t330 = t350 * t488 + t354 * t491;
t614 = t330 * t467;
t613 = t374 * t486;
t612 = t382 * t520;
t527 = -qJDD(3) * pkin(3) + t461 * t576 + qJDD(4);
t597 = t492 * t459;
t395 = t527 - t597;
t611 = t395 * t485;
t610 = t395 * t492;
t608 = t461 * t492;
t475 = cos(t482);
t607 = t475 * t490;
t603 = t618 * t492;
t602 = t489 * t490;
t601 = t489 * t493;
t596 = t492 * t493;
t595 = t493 * t475;
t495 = qJD(3) ^ 2;
t594 = t494 * t495;
t407 = -pkin(4) * t559 + t445;
t593 = pkin(5) * t585 - qJ(6) * t586 - qJD(6) * t441 - t407;
t592 = -qJ(6) * t579 - t635;
t591 = pkin(5) * t579 - t634;
t418 = t440 * t492;
t589 = -qJD(3) * t418 - t427 * t489 - t425;
t416 = t441 * t492;
t417 = t440 * t489;
t588 = -t440 * qJD(1) + qJD(3) * t416 - qJD(5) * t417;
t371 = qJD(1) * t419 + qJDD(1) * t446;
t390 = qJDD(3) * qJ(4) + t459 * t489 + (qJD(4) + t608) * qJD(3);
t343 = t485 * t371 + t486 * t390;
t584 = g(1) * t596 + g(2) * t599;
t583 = t493 * pkin(1) + t490 * qJ(2);
t483 = t489 ^ 2;
t484 = t492 ^ 2;
t582 = t483 - t484;
t581 = -t495 - t496;
t575 = qJD(3) * t492;
t329 = t350 * t491 - t354 * t488;
t571 = qJD(6) - t329;
t569 = qJDD(1) * qJ(2);
t568 = qJDD(1) * t485;
t567 = qJDD(1) * t486;
t564 = qJDD(3) * t489;
t562 = g(1) * t607;
t561 = 0.2e1 * qJD(1) * qJD(2);
t560 = t493 * pkin(7) + t583;
t473 = pkin(4) * t486 + pkin(3);
t556 = t488 * t576;
t555 = t491 * t576;
t545 = -g(2) * t601 + t619;
t543 = qJD(3) * pkin(3) - qJD(4);
t342 = t486 * t371 - t485 * t390;
t435 = pkin(4) * t606 - t492 * t494;
t542 = qJD(1) * t402 + t375;
t541 = qJDD(2) - t617;
t401 = -t485 * t600 + t434;
t540 = -t401 * qJDD(1) - t342;
t333 = pkin(4) * t513 - pkin(8) * t510 + t342;
t336 = pkin(8) * t509 + t343;
t537 = -t491 * t333 + t488 * t336 + t350 * t573 + t354 * t572;
t408 = t474 * t602 - t595;
t410 = t474 * t601 + t607;
t536 = g(1) * t410 + g(2) * t408;
t409 = t474 * t493 + t475 * t602;
t411 = -t474 * t490 + t489 * t595;
t535 = -g(1) * t411 - g(2) * t409;
t534 = g(1) * t493 + g(2) * t490;
t532 = g(2) * t492 * t595 + t439 * t519 + t475 * t620;
t528 = t616 + t621;
t526 = -t342 * t485 + t343 * t486;
t523 = -t374 * t485 + t375 * t486;
t522 = t380 * t491 - t392 * t488;
t423 = -pkin(4) * t557 + t489 * t574;
t518 = t561 + 0.2e1 * t569;
t517 = pkin(5) * t475 + qJ(6) * t474 + t473;
t516 = t639 * t485;
t422 = -t543 - t608;
t515 = -g(1) * t602 - t545;
t514 = 0.2e1 * qJ(2) * t570 + qJDD(3) * t494;
t512 = t488 * t333 + t491 * t336 + t350 * t572 - t354 * t573;
t511 = t488 * t358 + t491 * t370 + t380 * t572 - t392 * t573;
t340 = -t436 * t572 + t437 * t573 - t488 * t509 - t491 * t510;
t391 = -pkin(4) * t436 + t422;
t505 = g(1) * t408 - g(2) * t410 + t474 * t619 - t537;
t503 = t518 - t534;
t337 = pkin(5) * t382 - qJ(6) * t520 + t391;
t501 = t337 * t520 + qJDD(6) - t505;
t500 = -g(1) * t409 + g(2) * t411 - t475 * t619 + t512;
t499 = -pkin(4) * t509 + t527;
t497 = t341 * pkin(5) + t340 * qJ(6) - qJD(6) * t520 + t499;
t478 = t493 * qJ(2);
t476 = qJDD(3) * t492;
t415 = t441 * t489;
t388 = -t485 * t554 + t404;
t378 = pkin(5) * t440 - qJ(6) * t441 - t473;
t369 = -t485 * t555 - t486 * t556 + t492 * t630;
t367 = t427 * t492 - t485 * t556 + t486 * t555;
t353 = pkin(5) * t416 + qJ(6) * t418 + t435;
t352 = t499 - t597;
t346 = pkin(5) * t520 + qJ(6) * t382;
t345 = -pkin(5) * t489 - t522;
t344 = qJ(6) * t489 + t633;
t328 = pkin(5) * t369 + qJ(6) * t367 + qJD(6) * t418 + t423;
t327 = qJ(6) * t467 + t330;
t326 = -pkin(5) * t467 + t571;
t325 = -t340 + t640;
t324 = -pkin(5) * t575 - t625;
t323 = qJ(6) * t575 + qJD(6) * t489 + t511;
t322 = t497 - t597;
t321 = qJDD(6) + t537 - t622;
t320 = qJD(6) * t467 + t512 + t615;
t1 = [(-t340 * t489 - t367 * t467 - t418 * t439 + t520 * t575) * MDP(20) + (t320 * t489 + t322 * t418 + t323 * t467 + t327 * t575 - t328 * t520 + t337 * t367 + t340 * t353 + t344 * t439 - t536) * MDP(27) + (-t320 * t416 - t321 * t418 - t323 * t382 + t324 * t520 - t326 * t367 - t327 * t369 - t340 * t345 - t341 * t344 + t584) * MDP(26) + (t340 * t418 - t367 * t520) * MDP(18) + (t340 * t416 + t341 * t418 + t367 * t382 - t369 * t520) * MDP(19) + (t343 * t402 + t375 * t389 + t342 * t401 + t374 * t388 - g(1) * (pkin(3) * t601 - qJ(4) * t596 + t478) - g(2) * t560 + (t422 * t576 - t610) * t494 + (-g(1) * t494 - g(2) * t529) * t490) * MDP(17) + (t389 * t436 + t402 * t538 - t388 * t437 - t401 * t539 + (-t343 * t485 + t486 * t540) * t492 + (t485 * t542 + t613) * t576 + t584) * MDP(16) + (t514 * t492 + (t503 - t594) * t489) * MDP(12) + (-t514 * t489 + (t518 - t594) * t492 - t584) * MDP(13) + (-t341 * t489 - t369 * t467 - t382 * t575 - t416 * t439) * MDP(21) + (-t321 * t489 + t322 * t416 - t324 * t467 - t326 * t575 + t328 * t382 + t337 * t369 + t341 * t353 - t345 * t439 + t535) * MDP(25) + (t439 * t489 + t467 * t575) * MDP(22) + (-t492 * t495 - t564) * MDP(10) + t503 * MDP(5) + 0.2e1 * (-t489 * t565 + t570 * t582) * MDP(8) + (-t541 * pkin(1) - g(1) * (-pkin(1) * t490 + t478) - g(2) * t583 + (t561 + t569) * qJ(2)) * MDP(6) + (-t636 + (-t494 * t539 + (-t494 * t565 + t395) * t486 - t542 * qJD(3)) * t492 + (-t389 * qJD(1) - t402 * qJDD(1) - t343 + t534 * t485 + (-t422 * t486 + t437 * t494) * qJD(3)) * t489) * MDP(15) + (t329 * t575 + t435 * t341 + t352 * t416 + t391 * t369 + t423 * t382 + t522 * t439 + t467 * t625 - t537 * t489 + t535) * MDP(23) + (t320 * t344 + t327 * t323 + t322 * t353 + t337 * t328 + t321 * t345 + t326 * t324 - g(1) * (pkin(5) * t411 + qJ(6) * t410 + t473 * t601 - t596 * t618 + t478) - g(2) * (pkin(5) * t409 + qJ(6) * t408 + t493 * t623 + t560) + (-g(1) * (-t623 + t494) - g(2) * (t473 * t489 - t603)) * t490) * MDP(28) + (-t330 * t575 - t435 * t340 - t352 * t418 - t391 * t367 + t423 * t520 - t439 * t633 - t467 * t511 - t489 * t512 + t536) * MDP(24) + t534 * MDP(3) + (qJDD(1) * t484 - 0.2e1 * t489 * t549) * MDP(7) + qJDD(1) * MDP(1) + (-t489 * t495 + t476) * MDP(9) + (qJDD(2) + t639 - 0.2e1 * t617) * MDP(4) - t639 * MDP(2) + (-t516 + (t494 * t538 + t611 + (qJD(1) * t401 + t374) * qJD(3)) * t492 + (t388 * qJD(1) - t534 * t486 + (-t436 * t494 + (t494 * t579 - t422) * t485) * qJD(3) - t540) * t489) * MDP(14); qJDD(1) * MDP(4) - t496 * MDP(5) + (t481 - t528 + t541) * MDP(6) + (t489 * t581 + t476) * MDP(12) + (t492 * t581 - t564) * MDP(13) + (-t483 * t568 + t492 * t538 + (-t486 * t496 + (-t436 - t551) * qJD(3)) * t489) * MDP(14) + (-t483 * t567 - t492 * t510 + (t485 * t496 + (t437 - 0.2e1 * t558) * qJD(3)) * t489) * MDP(15) + ((t486 * t538 + ((t550 + t565) * t486 + t539) * t485) * t489 + (qJD(1) * t486 + t485 * t575) * t437 + (-qJD(1) * t485 + t486 * t575) * t436) * MDP(16) + (-t610 + t526 * t489 + (-t375 * t485 - t613) * qJD(1) + (t422 * t489 + t492 * t523) * qJD(3) + t639) * MDP(17) + (-t340 * t415 + t341 * t417 - t382 * t589 + t520 * t588) * MDP(26) + (-t320 * t417 + t321 * t415 - t322 * t492 + t326 * t588 + t327 * t589 + t337 * t576 + t639) * MDP(28) + t629 * (-t341 * t492 + t382 * t576 - t415 * t439 - t467 * t588) + t638 * (-t340 * t492 - t417 * t439 + t467 * t589 - t520 * t576); MDP(9) * t565 - MDP(10) * t566 + qJDD(3) * MDP(11) + ((-t546 - t616) * t492 + t544) * MDP(12) + ((-t459 + t528) * t489 + t545) * MDP(13) + (pkin(3) * t538 - t395 * t486 + (t636 + (-qJ(4) * t578 - t374) * qJD(1)) * t492 + (-qJ(4) * t568 + g(3) * t486 + t461 * t436 + (-t393 + (t422 + t543) * t485) * qJD(1)) * t489) * MDP(14) + (-pkin(3) * t539 + t611 + (-pkin(3) * t567 - t516 + (-qJ(4) * t577 + t375) * qJD(1)) * t492 + (-qJ(4) * t567 - g(3) * t485 - t461 * t437 + (t394 + (-qJD(4) + t422) * t486) * qJD(1)) * t489) * MDP(15) + (t393 * t437 - t394 * t436 + (qJD(4) * t436 - t374 * t580 + t343) * t486 + (qJD(4) * t437 - t375 * t580 - t342) * t485 + (t485 * t510 + t486 * t509) * qJ(4) + t515) * MDP(16) + (-t422 * t445 - t374 * t393 - t375 * t394 + t523 * qJD(4) + (-t395 - t626) * pkin(3) + (t489 * t639 + t526 - t619) * qJ(4)) * MDP(17) + (-t340 * t441 + t520 * t586) * MDP(18) + (t340 * t440 - t341 * t441 - t382 * t586 - t520 * t585) * MDP(19) + (t439 * t441 + t467 * t586 - t520 * t579) * MDP(20) + (t382 * t579 - t439 * t440 - t467 * t585) * MDP(21) - t467 * MDP(22) * t579 + (-t473 * t341 + t352 * t440 - t407 * t382 + (-qJD(1) * t329 - t562) * t492 + t634 * t467 + t585 * t391 + t532) * MDP(23) + (t330 * t579 + t473 * t340 + t352 * t441 + t586 * t391 - t407 * t520 + t467 * t635 - t642) * MDP(24) + (t322 * t440 + t341 * t378 + (qJD(1) * t326 - t562) * t492 - t591 * t467 + t593 * t382 + t585 * t337 + t532) * MDP(25) + (-t320 * t440 + t321 * t441 + t326 * t586 - t327 * t585 + t340 * t519 - t341 * t397 - t382 * t592 + t520 * t591 + t515) * MDP(26) + (-t322 * t441 - t327 * t579 - t337 * t586 + t340 * t378 + t467 * t592 - t520 * t593 + t642) * MDP(27) + (-g(3) * t603 + t320 * t397 - t321 * t519 + t322 * t378 + t591 * t326 + t592 * t327 + t593 * t337 + t517 * t620 + t639 * (t489 * t618 + t492 * t517)) * MDP(28) + (MDP(7) * t489 * t492 - MDP(8) * t582) * t496; ((t437 - t578) * t580 - t538) * MDP(14) + (t436 * t580 + t510) * MDP(15) + (-t436 ^ 2 - t437 ^ 2) * MDP(16) + (t374 * t437 - t375 * t436 + t527 + t631) * MDP(17) + (-t624 - t641) * MDP(26) + (-t326 * t520 + t327 * t382 + t497 + t631) * MDP(28) + t638 * (t340 + t640) + t629 * (t341 + t637); MDP(18) * t612 + (t624 - t641) * MDP(19) + t325 * MDP(20) + (-t341 + t637) * MDP(21) + t439 * MDP(22) + (-t391 * t520 + t505 + t614) * MDP(23) + (t329 * t467 + t382 * t391 - t500) * MDP(24) + (-t346 * t382 - t501 + t614 + 0.2e1 * t622) * MDP(25) + (pkin(5) * t340 - qJ(6) * t341 + (t327 - t330) * t520 + (t326 - t571) * t382) * MDP(26) + (0.2e1 * t615 - t337 * t382 + t346 * t520 + (0.2e1 * qJD(6) - t329) * t467 + t500) * MDP(27) + (t320 * qJ(6) - t321 * pkin(5) - t337 * t346 - t326 * t330 - g(1) * (-pkin(5) * t408 + qJ(6) * t409) - g(2) * (pkin(5) * t410 - qJ(6) * t411) - (-pkin(5) * t474 + qJ(6) * t475) * t619 + t571 * t327) * MDP(28); (-t439 + t612) * MDP(25) + t325 * MDP(26) + (-t467 ^ 2 - t624) * MDP(27) + (-t327 * t467 + t501 - t622) * MDP(28);];
tau  = t1;

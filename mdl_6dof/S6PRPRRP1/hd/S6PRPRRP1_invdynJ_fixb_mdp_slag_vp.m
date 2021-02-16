% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:24:34
% EndTime: 2021-01-16 01:24:45
% DurationCPUTime: 6.84s
% Computational Cost: add. (3790->503), mult. (7993->665), div. (0->0), fcn. (6245->18), ass. (0->236)
t504 = cos(pkin(11));
t510 = sin(qJ(2));
t503 = sin(pkin(6));
t583 = qJD(1) * t503;
t559 = t510 * t583;
t462 = t504 * t559;
t501 = sin(pkin(11));
t513 = cos(qJ(2));
t558 = t513 * t583;
t413 = t501 * t558 + t462;
t509 = sin(qJ(4));
t512 = cos(qJ(4));
t541 = pkin(4) * t509 - pkin(9) * t512;
t455 = t541 * qJD(4);
t653 = -t413 + t455;
t461 = t501 * t559;
t415 = t504 * t558 - t461;
t534 = -pkin(4) * t512 - pkin(9) * t509 - pkin(3);
t638 = pkin(2) * t504;
t441 = t534 - t638;
t508 = sin(qJ(5));
t511 = cos(qJ(5));
t572 = qJD(5) * t511;
t597 = t511 * t512;
t652 = -t415 * t597 + t441 * t572 + t653 * t508;
t582 = qJD(2) * t503;
t552 = qJD(1) * t582;
t567 = qJDD(1) * t503;
t651 = t510 * t567 + t513 * t552;
t576 = qJD(4) * t509;
t599 = t508 * t512;
t485 = pkin(2) * t501 + pkin(8);
t609 = t485 * t508;
t650 = t415 * t599 + t653 * t511 + t576 * t609;
t569 = qJD(2) * qJD(4);
t551 = t512 * t569;
t566 = qJDD(2) * t509;
t649 = qJD(4) * qJD(5) + t551 + t566;
t498 = qJ(2) + pkin(11);
t493 = cos(t498);
t502 = sin(pkin(10));
t472 = t502 * t493;
t492 = sin(t498);
t505 = cos(pkin(10));
t506 = cos(pkin(6));
t603 = t505 * t506;
t421 = t492 * t603 + t472;
t605 = t503 * t512;
t391 = t421 * t509 + t505 * t605;
t473 = t505 * t493;
t607 = t502 * t506;
t418 = t492 * t607 - t473;
t393 = t418 * t509 + t502 * t605;
t484 = t506 * t512;
t606 = t503 * t509;
t426 = t492 * t606 - t484;
t525 = -g(1) * t393 + g(2) * t391 + g(3) * t426;
t573 = qJD(5) * t509;
t550 = qJD(2) * t573;
t531 = (-qJDD(4) + t550) * t511;
t580 = qJD(2) * t512;
t389 = ((qJD(5) + t580) * qJD(4) + t566) * t508 + t531;
t474 = t513 * t567;
t417 = qJDD(2) * pkin(2) - t510 * t552 + t474;
t372 = t501 * t417 + t651 * t504;
t370 = qJDD(2) * pkin(8) + t372;
t457 = qJD(2) * pkin(2) + t558;
t410 = t501 * t457 + t462;
t403 = qJD(2) * pkin(8) + t410;
t475 = t506 * qJDD(1) + qJDD(3);
t478 = qJD(1) * t506 + qJD(3);
t575 = qJD(4) * t512;
t533 = t509 * t370 + t403 * t575 - t512 * t475 + t478 * t576;
t337 = -qJDD(4) * pkin(4) + t533;
t331 = t389 * pkin(5) + qJDD(6) + t337;
t648 = -t331 + t525;
t453 = t485 * t597;
t532 = pkin(5) * t509 - qJ(6) * t597;
t571 = qJD(6) * t511;
t624 = qJ(6) * t509;
t647 = -t509 * t571 + t532 * qJD(4) + (-t453 + (-t441 + t624) * t508) * qJD(5) + t650;
t578 = qJD(4) * t485;
t598 = t509 * t511;
t593 = (-qJ(6) * qJD(5) - t578) * t598 + (-qJD(6) * t509 + (-qJ(6) * qJD(4) - qJD(5) * t485) * t512) * t508 + t652;
t646 = -t509 * t403 + t478 * t512;
t424 = (t501 * t513 + t504 * t510) * t503;
t645 = -t424 * t509 + t484;
t390 = t418 * t512 - t502 * t606;
t392 = t421 * t512 - t505 * t606;
t490 = pkin(6) + t498;
t482 = cos(t490);
t491 = pkin(6) - t498;
t483 = cos(t491);
t443 = t482 + t483;
t608 = t502 * t492;
t639 = -t505 / 0.2e1;
t405 = t443 * t639 + t608;
t604 = t505 * t492;
t640 = t502 / 0.2e1;
t407 = t443 * t640 + t604;
t602 = t506 * t509;
t427 = t492 * t605 + t602;
t480 = sin(t490);
t481 = sin(t491);
t438 = -t481 / 0.2e1 - t480 / 0.2e1;
t644 = -g(3) * (-t427 * t508 + t438 * t511) - g(2) * (-t392 * t508 + t405 * t511) - g(1) * (t390 * t508 + t407 * t511);
t565 = t512 * qJDD(2);
t446 = t509 * t569 + qJDD(5) - t565;
t479 = -qJD(5) + t580;
t570 = t511 * qJD(4);
t555 = t512 * t570;
t524 = -t508 * t573 + t555;
t643 = t446 * t598 - t524 * t479;
t577 = qJD(4) * t508;
t581 = qJD(2) * t509;
t449 = t511 * t581 + t577;
t641 = t449 ^ 2;
t637 = pkin(5) * t446;
t447 = t508 * t581 - t570;
t636 = pkin(5) * t447;
t635 = pkin(5) * t508;
t628 = g(3) * t503;
t627 = qJ(6) + pkin(9);
t543 = t508 * qJDD(4) + t649 * t511;
t388 = t508 * t550 - t543;
t626 = qJ(6) * t388;
t625 = qJ(6) * t389;
t623 = t331 * t508;
t622 = t337 * t508;
t378 = t512 * t403 + t509 * t478;
t375 = qJD(4) * pkin(9) + t378;
t409 = t457 * t504 - t461;
t382 = t534 * qJD(2) - t409;
t343 = t375 * t511 + t382 * t508;
t339 = -qJ(6) * t447 + t343;
t621 = t339 * t479;
t442 = -t480 + t481;
t406 = t442 * t639 + t472;
t620 = t406 * t508;
t408 = t442 * t640 + t473;
t619 = t408 * t508;
t618 = t415 * t447;
t617 = t415 * t449;
t439 = t483 / 0.2e1 - t482 / 0.2e1;
t615 = t439 * t508;
t614 = t447 * t479;
t613 = t449 * t479;
t612 = t475 * t509;
t610 = t479 * t511;
t601 = t506 * t510;
t600 = t506 * t513;
t596 = qJDD(1) - g(3);
t342 = -t375 * t508 + t511 * t382;
t338 = -qJ(6) * t449 + t342;
t335 = -pkin(5) * t479 + t338;
t595 = -t338 + t335;
t454 = t541 * qJD(2);
t592 = -t508 * t454 - t511 * t646;
t591 = -t389 * t598 - t447 * t555;
t548 = qJD(5) * t627;
t557 = t508 * t580;
t589 = qJ(6) * t557 - t508 * t548 + t571 + t592;
t436 = t511 * t454;
t588 = -t532 * qJD(2) - t511 * t548 - t436 + (-qJD(6) + t646) * t508;
t586 = t508 * t441 + t453;
t499 = t509 ^ 2;
t585 = -t512 ^ 2 + t499;
t584 = MDP(22) * t508;
t579 = qJD(4) * t447;
t574 = qJD(5) * t508;
t564 = MDP(18) + MDP(20);
t563 = MDP(19) + MDP(21);
t561 = t493 * t605;
t560 = t525 * t508;
t556 = t479 * t577;
t554 = t479 * t574;
t553 = t509 * t572;
t547 = -qJD(6) - t636;
t546 = t388 * t512 + t449 * t576;
t336 = qJDD(4) * pkin(9) + qJD(4) * t646 + t370 * t512 + t612;
t371 = t417 * t504 - t651 * t501;
t349 = qJD(2) * t455 + t534 * qJDD(2) - t371;
t544 = t511 * t336 + t508 * t349 - t375 * t574 + t382 * t572;
t539 = -t335 * t511 - t339 * t508;
t538 = t335 * t508 - t339 * t511;
t401 = t424 * t512 + t602;
t535 = t501 * t510 - t504 * t513;
t423 = t535 * t503;
t363 = t401 * t511 + t423 * t508;
t362 = -t401 * t508 + t423 * t511;
t488 = pkin(5) * t511 + pkin(4);
t536 = t488 * t512 + t509 * t627;
t374 = -qJD(4) * pkin(4) - t646;
t530 = -t446 * t508 + t479 * t572;
t529 = -g(3) * t506 + (-g(1) * t502 + g(2) * t505) * t503;
t419 = t493 * t607 + t604;
t420 = t493 * t603 - t608;
t528 = -g(1) * (t408 * t511 + t419 * t599) - g(2) * (t406 * t511 - t420 * t599) - g(3) * (t439 * t511 - t508 * t561);
t527 = -g(1) * (-t419 * t597 + t619) - g(2) * (t420 * t597 + t620) - g(3) * (t511 * t561 + t615);
t526 = -g(1) * t390 + g(2) * t392 + g(3) * t427;
t523 = g(1) * t419 - g(2) * t420 - t493 * t628;
t522 = -pkin(9) * t446 - t479 * t374;
t402 = -qJD(2) * pkin(3) - t409;
t486 = -pkin(3) - t638;
t521 = -qJDD(4) * t485 + (qJD(2) * t486 + t402 + t415) * qJD(4);
t520 = -g(1) * (t390 * t511 - t407 * t508) - g(2) * (-t392 * t511 - t405 * t508) - g(3) * (-t427 * t511 - t438 * t508) - t544;
t347 = t511 * t349;
t519 = -t343 * qJD(5) - t336 * t508 + t347;
t514 = qJD(4) ^ 2;
t518 = -qJD(2) * t413 + t485 * t514 - t371 - t523 + (-pkin(3) + t486) * qJDD(2);
t517 = t519 + t644;
t516 = -g(1) * (-t502 * t600 - t505 * t510) - g(2) * (-t502 * t510 + t505 * t600) - t513 * t628;
t515 = qJD(2) ^ 2;
t468 = t627 * t511;
t467 = t627 * t508;
t466 = -pkin(3) * t501 + pkin(8) * t504;
t465 = qJDD(4) * t512 - t509 * t514;
t464 = qJDD(4) * t509 + t512 * t514;
t460 = pkin(3) * t504 + pkin(8) * t501 + pkin(2);
t445 = t447 ^ 2;
t433 = (t485 + t635) * t509;
t429 = t511 * t441;
t416 = t535 * t582;
t414 = qJD(2) * t424;
t404 = t485 * t575 + (t508 * t575 + t553) * pkin(5);
t387 = -t508 * t624 + t586;
t381 = -qJ(6) * t598 + t429 + (-pkin(5) - t609) * t512;
t367 = pkin(5) * t557 + t378;
t361 = t401 * qJD(4) - t416 * t509;
t360 = t645 * qJD(4) - t416 * t512;
t355 = t374 - t547;
t333 = -t363 * qJD(5) - t360 * t508 + t414 * t511;
t332 = t362 * qJD(5) + t360 * t511 + t414 * t508;
t330 = -qJD(6) * t447 + t544 - t625;
t327 = -qJD(6) * t449 + t519 + t626 + t637;
t1 = [t596 * MDP(1) + (-t371 * t423 + t372 * t424 - t409 * t414 - t410 * t416 + t475 * t506 - g(3)) * MDP(5) + (-qJD(4) * t361 + qJDD(4) * t645 - t423 * t565) * MDP(11) + (-qJD(4) * t360 - qJDD(4) * t401 + t423 * t566) * MDP(12) + (-t332 * t447 - t333 * t449 + t362 * t388 - t363 * t389) * MDP(22) + (t327 * t362 + t330 * t363 - t331 * t645 + t332 * t339 + t333 * t335 + t355 * t361 - g(3)) * MDP(23) + t563 * (t332 * t479 + t361 * t449 - t363 * t446 + t388 * t645) + t564 * (-t333 * t479 + t361 * t447 + t362 * t446 - t389 * t645) + ((-t414 * t512 + t423 * t576) * MDP(11) + (t414 * t509 + t423 * t575) * MDP(12)) * qJD(2) + ((qJDD(2) * t513 - t510 * t515) * MDP(3) + (-qJDD(2) * t510 - t513 * t515) * MDP(4)) * t503; qJDD(2) * MDP(2) + (t474 + t516) * MDP(3) + (-g(1) * (t502 * t601 - t505 * t513) - g(2) * (-t502 * t513 - t505 * t601) - t596 * t510 * t503) * MDP(4) + (t409 * t413 - t410 * t415 + (t371 * t504 + t372 * t501 + t516) * pkin(2)) * MDP(5) + (qJDD(2) * t499 + 0.2e1 * t509 * t551) * MDP(6) + 0.2e1 * (t509 * t565 - t585 * t569) * MDP(7) + t464 * MDP(8) + t465 * MDP(9) + (t521 * t509 - t518 * t512) * MDP(11) + (t518 * t509 + t521 * t512) * MDP(12) + (-t388 * t598 + t524 * t449) * MDP(13) + (-t449 * t553 + (-t449 * t575 + (qJD(5) * t447 + t388) * t509) * t508 + t591) * MDP(14) + (t546 + t643) * MDP(15) + ((t389 + t556) * t512 + (t530 - t579) * t509) * MDP(16) + (-t446 * t512 - t479 * t576) * MDP(17) + (t429 * t446 + (t441 * t574 - t650) * t479 + (t447 * t578 - t347 + (t479 * t485 + t375) * t572 + (qJD(4) * t374 + qJD(5) * t382 - t446 * t485 + t336) * t508) * t512 + (qJD(4) * t342 + t374 * t572 + t389 * t485 - t618 + t622) * t509 + t527) * MDP(18) + (-t586 * t446 + t652 * t479 + (-t485 * t554 + (t374 * t511 + t449 * t485) * qJD(4) + t544) * t512 + (-t374 * t574 + t337 * t511 - t485 * t388 - t617 + (-t485 * t610 - t343) * qJD(4)) * t509 + t528) * MDP(19) + (t381 * t446 + t389 * t433 + t404 * t447 + (t355 * t577 - t327) * t512 - t647 * t479 + (qJD(4) * t335 + t355 * t572 - t618 + t623) * t509 + t527) * MDP(20) + (-t387 * t446 - t388 * t433 + t404 * t449 + (t355 * t570 + t330) * t512 + t593 * t479 + (-qJD(4) * t339 + t331 * t511 - t355 * t574 - t617) * t509 + t528) * MDP(21) + (t381 * t388 - t387 * t389 - t647 * t449 - t593 * t447 + t539 * t575 + (t538 * qJD(5) - t327 * t511 - t330 * t508 + t523) * t509) * MDP(22) + (t330 * t387 + t327 * t381 + t331 * t433 - g(1) * (pkin(5) * t619 - (t460 * t505 + t466 * t607) * t510 + (-t460 * t607 + t466 * t505) * t513 - t536 * t419) - g(2) * (pkin(5) * t620 - (t460 * t502 - t466 * t603) * t510 + (t460 * t603 + t466 * t502) * t513 + t536 * t420) - g(3) * (pkin(5) * t615 + (t460 * t513 + t466 * t510 + t536 * t493) * t503) + (-t509 * t415 + t404) * t355 + t593 * t339 + t647 * t335) * MDP(23); (t529 + t475) * MDP(5) + t465 * MDP(11) - t464 * MDP(12) + t591 * MDP(22) + t529 * MDP(23) + (-t331 * MDP(23) + (-t538 * MDP(23) + t449 * t584) * qJD(4)) * t512 + t563 * (t546 - t643) + t564 * ((-t389 + t556) * t512 + (t530 + t579) * t509) + (-t388 * t584 + (qJD(4) * t355 - t327 * t508 + t330 * t511) * MDP(23) + ((t447 * t508 + t449 * t511) * MDP(22) + t539 * MDP(23)) * qJD(5)) * t509; MDP(8) * t566 + MDP(9) * t565 + qJDD(4) * MDP(10) + (qJD(4) * t378 - t402 * t581 + t525 - t533) * MDP(11) + (-t612 + (-qJD(2) * t402 - t370) * t512 + t526) * MDP(12) + (-t388 * t508 - t449 * t610) * MDP(13) + ((-t388 + t614) * t511 + (-t389 + t613) * t508) * MDP(14) + ((-t449 * t509 + t479 * t597) * qJD(2) - t530) * MDP(15) + (t554 + t446 * t511 + (t447 * t509 - t479 * t599) * qJD(2)) * MDP(16) + t479 * MDP(17) * t581 + (-t342 * t581 - pkin(4) * t389 - t378 * t447 + t436 * t479 + (-t479 * t646 + t522) * t508 + (pkin(9) * qJD(5) * t479 - t337 + t525) * t511) * MDP(18) + (t343 * t581 + pkin(4) * t388 + t622 - t378 * t449 + (-pkin(9) * t574 + t592) * t479 + t522 * t511 - t560) * MDP(19) + (-t335 * t581 - t367 * t447 - t389 * t488 - t446 * t467 - t588 * t479 + (-t355 * t580 + (t355 + t636) * qJD(5)) * t508 + t648 * t511) * MDP(20) + (t623 - t367 * t449 + t388 * t488 - t446 * t468 + t589 * t479 + (t355 * t511 + t449 * t635) * qJD(5) + (t339 * t509 - t355 * t597) * qJD(2) - t560) * MDP(21) + (-t388 * t467 - t389 * t468 - t588 * t449 - t589 * t447 + (t479 * t335 + t330) * t511 + (-t327 + t621) * t508 - t526) * MDP(22) + (t330 * t468 - t327 * t467 - t331 * t488 - g(1) * (-t390 * t627 + t393 * t488) - g(2) * (-t391 * t488 + t392 * t627) - g(3) * (-t426 * t488 + t427 * t627) + (pkin(5) * t574 - t367) * t355 + t589 * t339 + t588 * t335) * MDP(23) + (-t509 * t512 * MDP(6) + t585 * MDP(7)) * t515; t449 * t447 * MDP(13) + (-t445 + t641) * MDP(14) + (-t388 - t614) * MDP(15) + (-t389 - t613) * MDP(16) + t446 * MDP(17) + (-t343 * t479 - t374 * t449 + t517) * MDP(18) + (-t342 * t479 + t374 * t447 + t520) * MDP(19) + (0.2e1 * t637 + t626 - t621 + (-t355 + t547) * t449 + t517) * MDP(20) + (-pkin(5) * t641 + t625 - t338 * t479 + (qJD(6) + t355) * t447 + t520) * MDP(21) + (pkin(5) * t388 - t595 * t447) * MDP(22) + (t595 * t339 + (-t355 * t449 + t327 + t644) * pkin(5)) * MDP(23); (t531 - t613) * MDP(20) + (t543 + t614) * MDP(21) + (-t445 - t641) * MDP(22) + (t335 * t449 + t339 * t447 - t648) * MDP(23) + (t649 * MDP(20) - MDP(21) * t550) * t508;];
tau = t1;

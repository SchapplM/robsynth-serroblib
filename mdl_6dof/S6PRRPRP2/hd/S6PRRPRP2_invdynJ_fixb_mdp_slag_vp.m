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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:57
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
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
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:55:07
% EndTime: 2021-01-16 02:55:28
% DurationCPUTime: 10.20s
% Computational Cost: add. (5821->587), mult. (13485->752), div. (0->0), fcn. (10635->14), ass. (0->249)
t528 = sin(qJ(3));
t531 = cos(qJ(3));
t526 = qJ(4) + pkin(8);
t596 = qJD(3) * t526;
t470 = qJD(4) * t531 - t528 * t596;
t522 = sin(pkin(11));
t553 = -qJD(4) * t528 - t531 * t596;
t660 = cos(pkin(11));
t409 = t470 * t660 + t522 * t553;
t590 = t660 * t531;
t644 = t522 * t528;
t563 = t590 - t644;
t532 = cos(qJ(2));
t524 = sin(pkin(6));
t623 = qJD(1) * t524;
t603 = t532 * t623;
t445 = t563 * t603;
t678 = t409 - t445;
t503 = qJD(2) * t590;
t621 = qJD(2) * t528;
t476 = t522 * t621 - t503;
t469 = qJD(5) + t476;
t525 = cos(pkin(6));
t529 = sin(qJ(2));
t640 = t524 * t529;
t482 = t525 * t531 - t528 * t640;
t604 = t529 * t623;
t586 = qJD(2) * t526 + t604;
t622 = qJD(1) * t525;
t442 = -t528 * t586 + t531 * t622;
t443 = t528 * t622 + t531 * t586;
t592 = t660 * t443;
t387 = t442 * t522 + t592;
t527 = sin(qJ(5));
t530 = cos(qJ(5));
t578 = pkin(5) * t527 - qJ(6) * t530;
t631 = qJD(6) * t527 - t469 * t578 + t387;
t591 = t660 * t528;
t486 = t522 * t531 + t591;
t628 = t470 * t522 - t486 * t603 - t660 * t553;
t478 = t486 * qJD(3);
t481 = t563 * qJD(3);
t662 = qJD(3) * pkin(3);
t608 = t528 * t662;
t411 = pkin(4) * t478 - pkin(9) * t481 + t608;
t513 = pkin(3) * t531 + pkin(2);
t424 = -pkin(4) * t563 - pkin(9) * t486 - t513;
t494 = t526 * t531;
t447 = t494 * t660 - t526 * t644;
t617 = qJD(5) * t530;
t618 = qJD(5) * t527;
t677 = -t424 * t617 + t447 * t618 - t678 * t530 + (-t411 + t604) * t527;
t627 = t527 * t424 + t530 * t447;
t639 = t524 * t531;
t483 = t525 * t528 + t529 * t639;
t416 = t522 * t482 + t483 * t660;
t638 = t524 * t532;
t500 = t527 * t638;
t400 = t416 * t530 - t500;
t614 = qJD(2) * qJD(3);
t599 = t528 * t614;
t548 = qJDD(2) * t486 - t522 * t599;
t676 = qJD(3) * t503 + t548;
t661 = cos(pkin(10));
t593 = t661 * t532;
t523 = sin(pkin(10));
t642 = t523 * t529;
t473 = -t525 * t593 + t642;
t594 = t661 * t529;
t641 = t523 * t532;
t475 = t525 * t641 + t594;
t580 = g(1) * t475 + g(2) * t473;
t551 = -g(3) * t638 + t580;
t675 = t604 - t608;
t479 = t486 * qJD(2);
t674 = qJD(3) * t479;
t673 = pkin(3) * t599 + qJDD(4);
t434 = t522 * t443;
t438 = t442 + t662;
t384 = t438 * t660 - t434;
t379 = -qJD(3) * pkin(4) - t384;
t448 = -t530 * qJD(3) + t479 * t527;
t450 = qJD(3) * t527 + t479 * t530;
t361 = t448 * pkin(5) - t450 * qJ(6) + t379;
t612 = qJDD(2) * t528;
t576 = -qJDD(2) * t590 + t522 * t612;
t426 = qJD(2) * t478 + t576;
t422 = qJDD(5) + t426;
t668 = pkin(3) * t522;
t510 = pkin(9) + t668;
t647 = t510 * t422;
t672 = t361 * t469 - t647;
t597 = qJDD(1) * t638;
t615 = qJD(1) * qJD(2);
t600 = t529 * t615;
t577 = t524 * t600 - t597;
t658 = qJDD(2) * pkin(2);
t457 = t577 - t658;
t533 = qJD(3) ^ 2;
t671 = -pkin(8) * t533 + t524 * (-g(3) * t532 + t600) - t457 + t580 + t658;
t670 = t450 ^ 2;
t669 = t469 ^ 2;
t667 = pkin(3) * t528;
t666 = pkin(5) * t422;
t663 = qJD(2) * pkin(2);
t659 = qJ(6) * t422;
t385 = t522 * t438 + t592;
t380 = qJD(3) * pkin(9) + t385;
t468 = -qJD(2) * t513 + qJD(4) - t603;
t398 = pkin(4) * t476 - pkin(9) * t479 + t468;
t359 = t380 * t530 + t398 * t527;
t353 = qJ(6) * t469 + t359;
t657 = t353 * t469;
t656 = t359 * t469;
t613 = qJD(3) * qJD(5);
t381 = -t527 * qJDD(3) + t479 * t618 + (-t676 - t613) * t530;
t655 = t381 * t527;
t653 = t448 * t476;
t652 = t448 * t479;
t651 = t450 * t448;
t588 = t450 * t469;
t650 = t450 * t479;
t649 = t469 * t527;
t648 = t486 * t530;
t519 = qJ(3) + pkin(11);
t515 = cos(t519);
t646 = t515 * t527;
t645 = t515 * t530;
t643 = t523 * t524;
t417 = t527 * t422;
t418 = t530 * t422;
t636 = t530 * t532;
t635 = qJDD(1) - g(3);
t634 = qJ(6) * t478 - qJD(6) * t563 - t677;
t412 = t445 * t527 - t530 * t604;
t633 = -pkin(5) * t478 + qJD(5) * t627 + t409 * t527 - t411 * t530 - t412;
t579 = t530 * pkin(5) + t527 * qJ(6);
t632 = t578 * t481 + (qJD(5) * t579 - qJD(6) * t530) * t486 + t628;
t610 = t525 * qJDD(1);
t502 = t531 * t610;
t458 = qJDD(2) * pkin(8) + (qJDD(1) * t529 + t532 * t615) * t524;
t545 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t622 + t458;
t574 = t586 * qJD(3);
t375 = qJDD(3) * pkin(3) - t528 * t545 - t531 * t574 + t502;
t376 = (-t574 + t610) * t528 + t545 * t531;
t350 = t660 * t375 - t522 * t376;
t351 = t522 * t375 + t660 * t376;
t535 = -t530 * qJDD(3) + t527 * t676;
t382 = t450 * qJD(5) + t535;
t630 = -t527 * t382 - t448 * t617;
t389 = t442 * t660 - t434;
t609 = pkin(3) * t621;
t410 = pkin(4) * t479 + pkin(9) * t476 + t609;
t629 = t530 * t389 + t527 * t410;
t472 = t525 * t642 - t593;
t626 = -t472 * t526 - t475 * t513;
t625 = t513 * t638 + t526 * t640;
t520 = t528 ^ 2;
t624 = -t531 ^ 2 + t520;
t620 = qJD(2) * t529;
t619 = qJD(5) * t510;
t358 = -t380 * t527 + t398 * t530;
t616 = qJD(6) - t358;
t611 = qJDD(2) * t531;
t606 = t524 * t636;
t605 = t660 * pkin(3);
t602 = t524 * t620;
t601 = qJD(2) * t638;
t598 = t531 * t614;
t595 = t524 * t661;
t474 = t525 * t594 + t641;
t589 = -t473 * t513 + t474 * t526;
t514 = sin(t519);
t427 = t472 * t515 - t514 * t643;
t446 = t494 * t522 + t526 * t591;
t349 = qJDD(3) * pkin(9) + t351;
t366 = -pkin(3) * t611 + t426 * pkin(4) - pkin(9) * t676 + t457 + t673;
t585 = t527 * t349 - t530 * t366 + t380 * t617 + t398 * t618;
t584 = t528 * t601;
t583 = t523 * pkin(3) * t639 + t472 * t667;
t582 = pkin(4) * t515 + pkin(9) * t514;
t581 = g(1) * t472 - g(2) * t474;
t348 = -qJDD(3) * pkin(4) - t350;
t575 = t482 * pkin(3);
t534 = qJD(2) ^ 2;
t573 = qJDD(2) * t532 - t529 * t534;
t572 = pkin(4) + t579;
t570 = t417 + (t476 * t530 + t617) * t469;
t568 = -t469 * t618 - t476 * t649 + t418;
t567 = -g(1) * t523 + g(2) * t661;
t399 = t416 * t527 + t606;
t565 = t481 * t527 + t486 * t617;
t564 = -t481 * t530 + t486 * t618;
t562 = t530 * t349 + t527 * t366 - t380 * t618 + t398 * t617;
t560 = t379 * t469 - t647;
t404 = -t473 * t646 - t474 * t530;
t406 = t472 * t530 - t475 * t646;
t452 = t500 * t515 - t530 * t640;
t559 = g(1) * t406 + g(2) * t404 + g(3) * t452;
t405 = -t473 * t645 + t474 * t527;
t407 = -t472 * t527 - t475 * t645;
t453 = (t515 * t636 + t527 * t529) * t524;
t558 = -g(1) * t407 - g(2) * t405 - g(3) * t453;
t429 = t474 * t515 - t514 * t595;
t461 = t514 * t525 + t515 * t640;
t557 = g(1) * t427 - g(2) * t429 - g(3) * t461;
t428 = -t474 * t514 - t515 * t595;
t430 = t472 * t514 + t515 * t643;
t460 = -t514 * t640 + t515 * t525;
t556 = -g(1) * t430 - g(2) * t428 - g(3) * t460;
t555 = t573 * t524;
t552 = t483 * qJD(3);
t492 = -t603 - t663;
t550 = -qJD(2) * t492 - t458 - t581;
t549 = (-t474 * t528 - t531 * t595) * pkin(3);
t544 = -t469 * t619 + t556;
t393 = t429 * t527 - t473 * t530;
t395 = -t427 * t527 - t475 * t530;
t432 = t461 * t527 + t606;
t543 = g(1) * t395 + g(2) * t393 + g(3) * t432 - t585;
t343 = pkin(5) * t382 + qJ(6) * t381 - qJD(6) * t450 + t348;
t542 = -t343 + t544;
t541 = -pkin(8) * qJDD(3) + (t492 + t603 - t663) * qJD(3);
t425 = -qJDD(2) * t513 + t577 + t673;
t540 = -t552 - t584;
t394 = t429 * t530 + t473 * t527;
t396 = -t427 * t530 + t475 * t527;
t433 = t461 * t530 - t500;
t539 = -g(1) * t396 - g(2) * t394 - g(3) * t433 + t562;
t538 = t361 * t450 + qJDD(6) - t543;
t511 = -t605 - pkin(4);
t484 = -t605 - t572;
t441 = qJD(3) * t482 + t531 * t601;
t415 = -t482 * t660 + t483 * t522;
t397 = pkin(5) * t450 + qJ(6) * t448;
t390 = t486 * t578 + t446;
t388 = t441 * t660 + t522 * t540;
t386 = t441 * t522 - t540 * t660;
t370 = pkin(5) * t563 - t424 * t530 + t447 * t527;
t369 = -qJ(6) * t563 + t627;
t362 = t448 * t469 - t381;
t360 = -pkin(5) * t479 + t389 * t527 - t410 * t530;
t357 = qJ(6) * t479 + t629;
t356 = qJD(5) * t400 + t388 * t527 - t530 * t602;
t355 = -qJD(5) * t399 + t388 * t530 + t527 * t602;
t352 = -pkin(5) * t469 + t616;
t342 = qJDD(6) + t585 - t666;
t341 = qJD(6) * t469 + t562 + t659;
t1 = [t635 * MDP(1) + MDP(3) * t555 + (-qJDD(2) * t529 - t532 * t534) * t524 * MDP(4) + (t482 * qJDD(3) + t531 * t555 + (-t552 - 0.2e1 * t584) * qJD(3)) * MDP(10) + (-qJD(3) * t441 - qJDD(3) * t483 + (-t528 * t573 - t532 * t598) * t524) * MDP(11) + (-qJD(3) * t386 - qJDD(3) * t415 + (-t426 * t532 + t476 * t620) * t524) * MDP(12) + (-t388 * qJD(3) - t416 * qJDD(3) + (t479 * t620 - t532 * t676) * t524) * MDP(13) + (t386 * t479 - t388 * t476 + t415 * t676 - t416 * t426) * MDP(14) + (-t350 * t415 + t351 * t416 - t384 * t386 + t385 * t388 - g(3) + (-t425 * t532 + t468 * t620) * t524) * MDP(15) + (-t355 * t448 + t356 * t450 - t381 * t399 - t382 * t400) * MDP(24) + (t341 * t400 + t342 * t399 + t343 * t415 + t352 * t356 + t353 * t355 + t361 * t386 - g(3)) * MDP(26) + (MDP(21) + MDP(23)) * (-t356 * t469 + t415 * t382 + t386 * t448 - t399 * t422) + (-MDP(22) + MDP(25)) * (t355 * t469 + t381 * t415 - t386 * t450 + t400 * t422); (-g(1) * t626 - g(2) * t589 - g(3) * t625 - t350 * t446 + t351 * t447 - t628 * t384 + t385 * t678 - t425 * t513 - t675 * t468) * MDP(15) + (-qJD(3) * t678 - t447 * qJDD(3) + t425 * t486 + t468 * t481 - t513 * t676) * MDP(13) + (-g(3) * t640 - t350 * t486 + t351 * t563 - t384 * t481 - t385 * t478 - t447 * t426 + t446 * t676 - t476 * t678 + t581) * MDP(14) - (MDP(13) * t675 - MDP(14) * t628) * t479 + (t541 * t528 + t531 * t671) * MDP(10) + (-t528 * t671 + t541 * t531) * MDP(11) + (-t369 * t382 - t370 * t381 + (t352 * t530 - t353 * t527) * t481 + t633 * t450 - t634 * t448 + (-t341 * t527 + t342 * t530 + (-t352 * t527 - t353 * t530) * qJD(5)) * t486) * MDP(24) + (t551 + t597) * MDP(3) + ((-MDP(13) + MDP(24)) * t514 + t515 * MDP(12)) * t551 + (qJDD(3) * t528 + t531 * t533) * MDP(7) + (qJDD(3) * t531 - t528 * t533) * MDP(8) + (t348 * t648 - t359 * t478 - t564 * t379 - t446 * t381 - t627 * t422 + t628 * t450 + t469 * t677 + t562 * t563 + t559) * MDP(22) + (-t422 * t563 + t469 * t478) * MDP(20) + (-t341 * t563 - t343 * t648 + t353 * t478 + t361 * t564 + t369 * t422 + t381 * t390 - t450 * t632 + t469 * t634 - t559) * MDP(25) + (t382 * t563 - t417 * t486 - t448 * t478 - t469 * t565) * MDP(19) + (t381 * t563 + t418 * t486 + t450 * t478 - t469 * t564) * MDP(18) + (t343 * t486 * t527 + t342 * t563 - t352 * t478 + t361 * t565 - t370 * t422 + t382 * t390 + t448 * t632 - t469 * t633 + t558) * MDP(23) + (t585 * t563 + t358 * t478 + t446 * t382 + t412 * t469 + t628 * t448 + ((-qJD(5) * t447 + t411) * t469 + t424 * t422 + t379 * qJD(5) * t486) * t530 + ((-qJD(5) * t424 - t409) * t469 - t447 * t422 + t348 * t486 + t379 * t481) * t527 + t558) * MDP(21) + (-t476 * t604 - qJDD(3) * t446 - t425 * t563 - t426 * t513 + t468 * t478 + (t476 * t667 - t628) * qJD(3)) * MDP(12) + qJDD(2) * MDP(2) + ((-t448 * t530 - t450 * t527) * t481 + (t655 - t382 * t530 + (t448 * t527 - t450 * t530) * qJD(5)) * t486) * MDP(17) + (-t381 * t648 - t450 * t564) * MDP(16) + (-t635 * t640 - t581) * MDP(4) + (t341 * t369 + t343 * t390 + t342 * t370 - g(1) * (pkin(5) * t407 + qJ(6) * t406 - t475 * t582 + t626) - g(2) * (pkin(5) * t405 + qJ(6) * t404 - t473 * t582 + t589) - g(3) * (pkin(5) * t453 + qJ(6) * t452 + t582 * t638 + t625) + t632 * t361 + t634 * t353 + t633 * t352) * MDP(26) + 0.2e1 * (t528 * t611 - t614 * t624) * MDP(6) + (qJDD(2) * t520 + 0.2e1 * t528 * t598) * MDP(5); MDP(7) * t612 + MDP(8) * t611 + qJDD(3) * MDP(9) + (-g(3) * t482 + t528 * t550 + t567 * t639 + t502) * MDP(10) + (g(3) * t483 + (-t524 * t567 - t610) * t528 + t550 * t531) * MDP(11) + (t387 * qJD(3) - t468 * t479 + (qJDD(3) * t660 - t476 * t621) * pkin(3) + t556 + t350) * MDP(12) + (qJD(3) * t389 + t468 * t476 + (-qJDD(3) * t522 - t479 * t621) * pkin(3) - t557 - t351) * MDP(13) + (-t426 * t668 - t676 * t605 - (-t385 + t387) * t479 + (-t384 + t389) * t476) * MDP(14) + (-g(1) * t583 - g(2) * t549 - g(3) * t575 + t350 * t605 + t351 * t668 + t384 * t387 - t385 * t389 - t468 * t609) * MDP(15) + (t530 * t588 - t655) * MDP(16) + ((-t381 - t653) * t530 - t450 * t649 + t630) * MDP(17) + (t570 - t650) * MDP(18) + (t568 + t652) * MDP(19) - t469 * t479 * MDP(20) + (-t358 * t479 + t511 * t382 - t387 * t448 + (t389 * t469 + t560) * t527 + (-t348 + (-t410 - t619) * t469 + t556) * t530) * MDP(21) + (-t511 * t381 + t629 * t469 + t359 * t479 - t387 * t450 + t560 * t530 + (t348 - t544) * t527) * MDP(22) + (t352 * t479 + t360 * t469 + t382 * t484 - t631 * t448 + t527 * t672 + t542 * t530) * MDP(23) + (t357 * t448 - t360 * t450 + (t352 * t476 - t382 * t510 + t341 + (t450 * t510 + t352) * qJD(5)) * t530 + (-t353 * t476 - t381 * t510 + t342 + (t448 * t510 - t353) * qJD(5)) * t527 + t557) * MDP(24) + (-t353 * t479 - t357 * t469 + t381 * t484 + t631 * t450 + t542 * t527 - t530 * t672) * MDP(25) + (t343 * t484 - t353 * t357 - t352 * t360 - g(1) * (-pkin(9) * t427 + t430 * t572 + t583) - g(2) * (t429 * pkin(9) + t428 * t572 + t549) - g(3) * (pkin(9) * t461 + t460 * t572 + t575) - t631 * t361 + (t341 * t530 + t342 * t527 + t352 * t617 - t353 * t618) * t510) * MDP(26) + (-MDP(5) * t528 * t531 + MDP(6) * t624) * t534; (t576 + 0.2e1 * t674) * MDP(12) + ((t503 - t476) * qJD(3) + t548) * MDP(13) + (-t476 ^ 2 - t479 ^ 2) * MDP(14) + (t384 * t479 + t385 * t476 + t425 - t551) * MDP(15) + (t568 - t652) * MDP(21) + (-t530 * t669 - t417 - t650) * MDP(22) + (-t469 * t649 + t418 - t652) * MDP(23) + ((t381 - t653) * t530 + t527 * t588 + t630) * MDP(24) + (t570 + t650) * MDP(25) + (-t361 * t479 + (-t342 + t657) * t530 + (t352 * t469 + t341) * t527 - t551) * MDP(26); MDP(16) * t651 + (-t448 ^ 2 + t670) * MDP(17) + t362 * MDP(18) + (-t479 * t617 - t527 * t613 - t535 + t588) * MDP(19) + t422 * MDP(20) + (-t379 * t450 + t543 + t656) * MDP(21) + (t358 * t469 + t379 * t448 - t539) * MDP(22) + (-t397 * t448 - t538 + t656 + 0.2e1 * t666) * MDP(23) + (pkin(5) * t381 - qJ(6) * t382 + (t353 - t359) * t450 + (t352 - t616) * t448) * MDP(24) + (0.2e1 * t659 - t361 * t448 + t397 * t450 + (0.2e1 * qJD(6) - t358) * t469 + t539) * MDP(25) + (t341 * qJ(6) - t342 * pkin(5) - t361 * t397 - t352 * t359 - g(1) * (-pkin(5) * t395 + qJ(6) * t396) - g(2) * (-pkin(5) * t393 + qJ(6) * t394) - g(3) * (-pkin(5) * t432 + qJ(6) * t433) + t616 * t353) * MDP(26); (-qJDD(5) - t576 + t651 - t674) * MDP(23) + t362 * MDP(24) + (-t669 - t670) * MDP(25) + (t538 - t657 - t666) * MDP(26);];
tau = t1;

% Calculate vector of inverse dynamics joint torques for
% S5RRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:03:08
% EndTime: 2021-01-15 22:03:30
% DurationCPUTime: 8.58s
% Computational Cost: add. (5486->583), mult. (15497->803), div. (0->0), fcn. (12610->18), ass. (0->282)
t553 = sin(pkin(10));
t554 = sin(pkin(5));
t562 = cos(qJ(2));
t704 = cos(pkin(10));
t619 = t704 * t562;
t600 = qJD(1) * t619;
t558 = sin(qJ(2));
t659 = qJD(1) * t558;
t631 = t554 * t659;
t481 = t553 * t631 - t554 * t600;
t478 = qJD(4) + t481;
t580 = t553 * t562 + t558 * t704;
t492 = t580 * t554;
t484 = qJD(1) * t492;
t561 = cos(qJ(4));
t555 = cos(pkin(5));
t660 = qJD(1) * t555;
t607 = qJD(2) + t660;
t520 = t561 * t607;
t557 = sin(qJ(4));
t453 = t484 * t557 - t520;
t452 = qJD(5) + t453;
t559 = sin(qJ(1));
t673 = t559 * t562;
t563 = cos(qJ(1));
t677 = t558 * t563;
t506 = -t555 * t673 - t677;
t667 = t562 * t563;
t678 = t558 * t559;
t583 = t555 * t667 - t678;
t724 = -g(1) * t506 - g(2) * t583;
t572 = qJD(2) * t580;
t601 = t554 * t619;
t648 = qJDD(1) * t558;
t448 = -qJDD(1) * t601 + (qJD(1) * t572 + t553 * t648) * t554;
t447 = qJDD(4) + t448;
t709 = pkin(1) * t562;
t534 = t555 * t709;
t528 = qJD(1) * t534;
t705 = pkin(7) + qJ(3);
t626 = t705 * t558;
t603 = t554 * t626;
t575 = t555 * pkin(2) - t603;
t456 = qJD(2) * pkin(2) + qJD(1) * t575 + t528;
t682 = t555 * t558;
t533 = pkin(1) * t682;
t627 = t554 * t705;
t470 = (t562 * t627 + t533) * qJD(1);
t620 = t704 * t470;
t410 = t553 * t456 + t620;
t403 = pkin(8) * t607 + t410;
t544 = pkin(2) * t562 + pkin(1);
t512 = t544 * t554;
t503 = -qJD(1) * t512 + qJD(3);
t420 = pkin(3) * t481 - pkin(8) * t484 + t503;
t383 = t403 * t561 + t420 * t557;
t647 = qJDD(1) * t562;
t643 = pkin(1) * t647;
t527 = t555 * t643;
t649 = qJDD(1) * t555;
t531 = qJDD(2) + t649;
t644 = pkin(1) * qJD(2) * t555;
t604 = qJD(1) * t644;
t621 = qJD(2) * t705;
t657 = qJD(3) * t558;
t407 = -t558 * t604 + pkin(2) * t531 + t527 + (-qJDD(1) * t626 + (-t562 * t621 - t657) * qJD(1)) * t554;
t574 = qJD(3) * t562 - t558 * t621;
t623 = t554 * t647;
t632 = pkin(7) * t623 + qJDD(1) * t533 + t562 * t604;
t417 = (qJ(3) * t647 + qJD(1) * t574) * t554 + t632;
t381 = t553 * t407 + t704 * t417;
t378 = pkin(8) * t531 + t381;
t658 = qJD(2) * t558;
t630 = t554 * t658;
t602 = qJD(1) * t630;
t449 = -t553 * t602 + (qJD(2) * t600 + qJDD(1) * t580) * t554;
t596 = t544 * qJDD(1);
t646 = pkin(2) * t602 + qJDD(3);
t468 = -t554 * t596 + t646;
t393 = pkin(3) * t448 - pkin(8) * t449 + t468;
t615 = t378 * t557 - t561 * t393;
t715 = t383 * qJD(4) + t615;
t364 = -pkin(4) * t447 + t715;
t550 = qJ(2) + pkin(10);
t547 = sin(t550);
t666 = t563 * t547;
t548 = cos(t550);
t675 = t559 * t548;
t490 = t555 * t666 + t675;
t570 = -t561 * t484 - t557 * t607;
t597 = g(1) * t559 - g(2) * t563;
t685 = t554 * t561;
t539 = t555 * t561;
t687 = t554 * t557;
t707 = g(3) * (-t547 * t687 + t539);
t530 = t563 * t548;
t676 = t559 * t547;
t718 = t555 * t676 - t530;
t723 = t557 * (g(1) * t718 - g(2) * t490) + (-pkin(4) * t570 + pkin(9) * t452) * t452 + t597 * t685 + t364 + t707;
t546 = pkin(5) - t550;
t722 = sin(t546) / 0.2e1;
t560 = cos(qJ(5));
t556 = sin(qJ(5));
t693 = t570 * t556;
t414 = -t560 * t478 - t693;
t721 = t414 * t478;
t670 = t560 * t561;
t441 = -t481 * t670 + t484 * t556;
t652 = qJD(5) * t556;
t653 = qJD(4) * t561;
t571 = -t557 * t652 + t560 * t653 - t441;
t720 = t571 * t452;
t467 = t534 + t575;
t684 = t554 * t562;
t662 = pkin(7) * t684 + t533;
t479 = qJ(3) * t684 + t662;
t433 = t553 * t467 + t704 * t479;
t422 = pkin(8) * t555 + t433;
t686 = t554 * t558;
t491 = t553 * t686 - t601;
t439 = pkin(3) * t491 - pkin(8) * t492 - t512;
t719 = t561 * t422 + t557 * t439;
t717 = -t492 * t557 + t539;
t716 = (-t684 * t705 - t533) * qJD(2) - t554 * t657;
t398 = -qJD(4) * t570 + t557 * t449 - t561 * t531;
t549 = t554 ^ 2;
t712 = 0.2e1 * t549;
t564 = qJD(1) ^ 2;
t545 = pkin(5) + t550;
t535 = sin(t545);
t711 = -t535 / 0.2e1 + t722;
t710 = pkin(1) * t549;
t706 = g(3) * t562;
t654 = qJD(4) * t557;
t397 = qJD(4) * t520 + t561 * t449 - t484 * t654 + t557 * t531;
t651 = qJD(5) * t560;
t633 = t560 * t397 + t556 * t447 + t478 * t651;
t371 = t570 * t652 + t633;
t703 = t371 * t556;
t416 = t478 * t556 - t560 * t570;
t614 = t397 * t556 - t560 * t447;
t372 = qJD(5) * t416 + t614;
t702 = t372 * t561;
t396 = qJDD(5) + t398;
t701 = t396 * t556;
t700 = t414 * t452;
t699 = t416 * t452;
t698 = t416 * t481;
t610 = t452 * t560;
t697 = t453 * t478;
t696 = t453 * t484;
t695 = t570 * t478;
t694 = t570 * t484;
t692 = t478 * t557;
t690 = t531 * MDP(8);
t689 = t548 * t554;
t688 = t549 * t564;
t461 = t553 * t470;
t683 = t555 * t557;
t681 = t556 * t559;
t680 = t556 * t563;
t679 = t557 * t447;
t674 = t559 * t561;
t672 = t560 * t396;
t671 = t560 * t559;
t669 = t560 * t563;
t668 = t561 * t563;
t469 = -qJD(1) * t603 + t528;
t423 = t469 * t553 + t620;
t665 = t423 - t478 * (pkin(4) * t557 - pkin(9) * t561);
t380 = t704 * t407 - t553 * t417;
t424 = t469 * t704 - t461;
t435 = pkin(2) * t631 + pkin(3) * t484 + pkin(8) * t481;
t663 = t561 * t424 + t557 * t435;
t551 = t558 ^ 2;
t661 = -t562 ^ 2 + t551;
t382 = -t403 * t557 + t420 * t561;
t374 = -pkin(4) * t478 - t382;
t656 = qJD(4) * t374;
t540 = pkin(2) * t553 + pkin(8);
t655 = qJD(4) * t540;
t650 = qJD(1) * qJD(2);
t645 = g(3) * t547 * t554;
t641 = t556 * t689;
t640 = t560 * t689;
t639 = t558 * t688;
t638 = t559 * t687;
t637 = t556 * t674;
t636 = t556 * t668;
t635 = t559 * t670;
t634 = t560 * t668;
t524 = t563 * t687;
t628 = t554 * t555 * t564;
t625 = t562 * t650;
t624 = t554 * t648;
t578 = -t561 * t378 - t557 * t393 + t403 * t654 - t420 * t653;
t363 = pkin(9) * t447 - t578;
t377 = -pkin(3) * t531 - t380;
t366 = pkin(4) * t398 - pkin(9) * t397 + t377;
t617 = -t556 * t363 + t560 * t366;
t616 = -t371 * t561 + t416 * t654;
t529 = t562 * t644;
t457 = t554 * t574 + t529;
t404 = t457 * t553 - t704 * t716;
t612 = t490 * t561 - t524;
t611 = t478 * t561;
t541 = -pkin(2) * t704 - pkin(3);
t508 = -t561 * pkin(4) - t557 * pkin(9) + t541;
t609 = pkin(9) * t484 - qJD(5) * t508 + t663;
t608 = pkin(2) * t630;
t606 = qJD(2) + 0.2e1 * t660;
t605 = t531 + t649;
t598 = g(1) * t563 + g(2) * t559;
t409 = t456 * t704 - t461;
t432 = t467 * t704 - t553 * t479;
t595 = t560 * t363 + t556 * t366;
t375 = pkin(9) * t478 + t383;
t402 = -pkin(3) * t607 - t409;
t379 = t453 * pkin(4) + pkin(9) * t570 + t402;
t368 = t375 * t560 + t379 * t556;
t594 = t375 * t556 - t379 * t560;
t389 = pkin(9) * t491 + t719;
t421 = -t555 * pkin(3) - t432;
t465 = t492 * t561 + t683;
t392 = -pkin(4) * t717 - t465 * pkin(9) + t421;
t593 = t389 * t560 + t392 * t556;
t592 = -t389 * t556 + t392 * t560;
t405 = t457 * t704 + t553 * t716;
t483 = t554 * t572;
t486 = (-t553 * t558 + t619) * t554 * qJD(2);
t436 = pkin(3) * t483 - pkin(8) * t486 + t608;
t591 = -t405 * t557 + t436 * t561;
t589 = -t422 * t557 + t439 * t561;
t438 = t465 * t560 + t491 * t556;
t437 = t465 * t556 - t560 * t491;
t537 = cos(t545);
t538 = cos(t546);
t511 = t537 + t538;
t588 = t676 - t563 * t511 / 0.2e1;
t587 = t666 + t559 * t511 / 0.2e1;
t586 = t561 * t447 - t478 * t654 - t481 * t692;
t585 = -t561 * t718 + t638;
t440 = -t481 * t556 * t561 - t560 * t484;
t582 = (-t556 * t653 + t440) * t452;
t581 = -t452 * t651 - t701;
t577 = t561 * t405 - t422 * t654 + t557 * t436 + t439 * t653;
t576 = t402 * t478 - t540 * t447;
t573 = t662 * t555;
t569 = g(1) * (t555 * t675 + t666) - g(2) * (t530 * t555 - t676) - g(3) * t689 - t377;
t568 = -pkin(9) * t396 + (t374 + t382) * t452;
t494 = t555 * t635 - t680;
t502 = t555 * t681 + t634;
t567 = -t494 * t547 + t502 * t548 + t560 * t638;
t495 = t555 * t636 - t671;
t500 = -t555 * t669 - t637;
t566 = -t495 * t547 + t500 * t548 + t524 * t556;
t509 = pkin(2) * t682 - t627;
t507 = -t555 * t678 + t667;
t505 = -t555 * t677 - t673;
t501 = t555 * t671 - t636;
t499 = t555 * t680 - t635;
t498 = t547 * t685 + t683;
t496 = t555 * t634 + t681;
t493 = t555 * t637 + t669;
t473 = t559 * t711 + t530;
t472 = t563 * t711 - t675;
t459 = t554 * t674 + t557 * t718;
t458 = t490 * t557 + t554 * t668;
t431 = qJD(4) * t465 + t486 * t557;
t430 = qJD(4) * t717 + t486 * t561;
t426 = t493 * t547 + t501 * t548 - t556 * t638;
t425 = -t496 * t547 + t499 * t548 + t524 * t560;
t388 = -pkin(4) * t491 - t589;
t387 = qJD(5) * t438 + t430 * t556 - t560 * t483;
t386 = -qJD(5) * t437 + t430 * t560 + t483 * t556;
t384 = -pkin(4) * t484 + t424 * t557 - t435 * t561;
t373 = pkin(4) * t431 - pkin(9) * t430 + t404;
t370 = -pkin(4) * t483 + qJD(4) * t719 - t591;
t369 = pkin(9) * t483 + t577;
t362 = -t368 * qJD(5) + t617;
t361 = -t594 * qJD(5) + t595;
t1 = [(-t380 * t492 - t381 * t491 + t404 * t484 - t405 * t481 - t409 * t486 - t410 * t483 - t432 * t449 - t433 * t448) * MDP(13) + (-(qJD(5) * t592 + t369 * t560 + t373 * t556) * t452 - t593 * t396 + t361 * t717 - t368 * t431 + t370 * t416 + t388 * t371 + t364 * t438 + t374 * t386 + g(1) * t566 - g(2) * t426) * MDP(28) + ((-qJD(5) * t593 - t369 * t556 + t373 * t560) * t452 + t592 * t396 - t362 * t717 - t594 * t431 + t370 * t414 + t388 * t372 + t364 * t437 + t374 * t387 - g(1) * t425 - g(2) * t567) * MDP(27) + (t397 * t717 - t398 * t465 - t430 * t453 + t431 * t570) * MDP(16) + (t591 * t478 + t589 * t447 - t615 * t491 + t382 * t483 + t404 * t453 + t421 * t398 - t377 * t717 + t402 * t431 + g(1) * t612 - g(2) * t585 + (-t383 * t491 - t478 * t719) * qJD(4)) * MDP(20) + (-t398 * t491 - t431 * t478 + t447 * t717 - t453 * t483) * MDP(18) + (t372 * t717 - t387 * t452 - t396 * t437 - t414 * t431) * MDP(25) + (-t371 * t717 + t386 * t452 + t396 * t438 + t416 * t431) * MDP(24) + (-t396 * t717 + t431 * t452) * MDP(26) + t597 * MDP(2) + t598 * MDP(3) + (t397 * t465 - t430 * t570) * MDP(15) + (t397 * t491 + t430 * t478 + t447 * t465 - t483 * t570) * MDP(17) + (-(-pkin(7) * t630 + t529) * t607 - t662 * t531 - (-pkin(7) * t602 + t632) * t555 + g(1) * t583 - g(2) * t506 + 0.2e1 * (-t625 - t648) * t710) * MDP(10) + (-g(1) * t458 - g(2) * t459 + t377 * t465 - t383 * t483 + t421 * t397 + t402 * t430 - t404 * t570 - t447 * t719 - t478 * t577 + t491 * t578) * MDP(21) + (t643 * t712 + (-pkin(7) * t686 + t534) * t531 + (-pkin(7) * t624 + t527) * t555 - g(1) * t505 - g(2) * t507 - t662 * qJD(2) ^ 2 + 0.2e1 * (-t558 * t710 - t573) * t650) * MDP(9) + (t558 * t647 - t650 * t661) * MDP(5) * t712 + t555 * t690 + (t447 * t491 + t478 * t483) * MDP(19) + (t371 * t438 + t386 * t416) * MDP(22) + (-t371 * t437 - t372 * t438 - t386 * t414 - t387 * t416) * MDP(23) + (-g(1) * t472 - g(2) * t473 + t380 * t555 - t404 * t607 + t432 * t531 - t512 * t448 + t468 * t491 + t481 * t608 + t503 * t483) * MDP(11) + (t381 * t433 + t410 * t405 + t380 * t432 - t409 * t404 - t468 * t512 + t503 * t608 - g(1) * (-t509 * t563 - t544 * t559) - g(2) * (-t509 * t559 + t544 * t563)) * MDP(14) + (-g(1) * t588 + g(2) * t587 - t381 * t555 - t405 * t607 - t433 * t531 - t512 * t449 + t468 * t492 + t484 * t608 + t503 * t486) * MDP(12) + (qJDD(1) * t551 + 0.2e1 * t558 * t625) * t549 * MDP(4) + qJDD(1) * MDP(1) + (-t598 * MDP(13) + (t562 * t605 - t606 * t658) * MDP(7) + (qJD(2) * t562 * t606 + t558 * t605) * MDP(6)) * t554; (-t562 * t628 + t624) * MDP(6) + (t558 * t628 + t623) * MDP(7) + (t508 * t672 - t384 * t414 - t374 * t440 - g(1) * (-t494 * t548 - t502 * t547) - g(2) * (t496 * t548 + t499 * t547) - t556 * t645 + (t556 * t609 - t560 * t665) * t452 + (-g(3) * t640 + t556 * t656 - t362 + (qJD(4) * t414 + t581) * t540) * t561 + (t374 * t651 + t364 * t556 - t594 * t481 + t540 * t372 + (t452 * t540 * t556 - t594) * qJD(4)) * t557) * MDP(27) + (t397 * t557 - t570 * t611) * MDP(15) + (t541 * t397 + t663 * t478 + t383 * t484 + t423 * t570 + t576 * t561 + (t478 * t655 - t569) * t557) * MDP(21) - t478 * t484 * MDP(19) + t661 * MDP(5) * t688 + (-t508 * t701 - t384 * t416 - t374 * t441 - g(1) * (t493 * t548 - t501 * t547) - g(2) * (-t495 * t548 - t500 * t547) - t560 * t645 + (t556 * t665 + t560 * t609) * t452 + (g(3) * t641 + t560 * t656 + t361 + (qJD(4) * t416 + t452 * t652 - t672) * t540) * t561 + (-t374 * t652 + t364 * t560 - t368 * t481 + t540 * t371 + (t540 * t610 - t368) * qJD(4)) * t557) * MDP(28) - t562 * MDP(4) * t639 + (t423 * t607 - t503 * t484 + g(1) * t587 + g(2) * t588 - g(3) * (t722 + t535 / 0.2e1) + (-t481 * t631 + t531 * t704) * pkin(2) + t380) * MDP(11) + (t414 * t441 + t416 * t440 + (-t414 * t560 - t416 * t556) * t653 + (-t703 - t372 * t560 + (t414 * t556 - t416 * t560) * qJD(5)) * t557) * MDP(23) + (t702 + t582 + (t581 - t721) * t557) * MDP(25) + (t371 * t557 * t560 + t416 * t571) * MDP(22) + (t688 * t709 + (-pkin(7) * t631 + t528) * t660 + g(1) * t507 - g(2) * t505 + g(3) * t686 + t528 * qJD(2) - t632) * MDP(10) + ((t672 + t698) * t557 + t720 + t616) * MDP(24) + t690 + (-t396 * t561 + t452 * t692) * MDP(26) + ((t410 - t423) * t484 + (-t409 + t424) * t481 + (-t448 * t553 - t449 * t704) * pkin(2)) * MDP(13) + (-t382 * t484 + t541 * t398 - t423 * t453 + (t424 * t478 + t576) * t557 + ((-t435 - t655) * t478 + t569) * t561) * MDP(20) + (t424 * t607 + t503 * t481 + g(1) * t473 - g(2) * t472 - g(3) * (-t538 / 0.2e1 + t537 / 0.2e1) + (-t484 * t631 - t531 * t553) * pkin(2) - t381) * MDP(12) + ((t397 - t697) * t561 + (-t398 + t695) * t557) * MDP(16) + (t586 + t696) * MDP(18) + (t478 * t611 + t679 + t694) * MDP(17) + (pkin(1) * t639 + t527 + (-pkin(7) * t648 - t706) * t554 + t564 * t573 + t724) * MDP(9) + (t409 * t423 - t410 * t424 + (t381 * t553 + t380 * t704 + (-t503 * t659 - t706) * t554 + t724) * pkin(2)) * MDP(14); (t484 * t607 + t448) * MDP(11) + (-t481 * t607 + t449) * MDP(12) + (-t481 ^ 2 - t484 ^ 2) * MDP(13) + (-g(3) * t555 + t409 * t484 + t410 * t481 + (-t596 - t597) * t554 + t646) * MDP(14) + (t586 - t696) * MDP(20) + (-t478 ^ 2 * t561 - t679 + t694) * MDP(21) + (-t702 + t582 + (t581 + t721) * t557) * MDP(27) + ((-t672 + t698) * t557 - t720 + t616) * MDP(28); -t453 ^ 2 * MDP(16) + (t397 + t697) * MDP(17) + (-t398 - t695) * MDP(18) + t447 * MDP(19) + (-g(1) * t459 + g(2) * t458 + t383 * t478 - t707 - t715) * MDP(20) + (g(1) * t585 + g(2) * t612 + g(3) * t498 + t382 * t478 + t402 * t453 + t578) * MDP(21) + (t416 * t610 + t703) * MDP(22) + ((t371 - t700) * t560 + (-t372 - t699) * t556) * MDP(23) + (t452 * t610 + t701) * MDP(24) + (-t452 ^ 2 * t556 + t672) * MDP(25) + (-pkin(4) * t372 - t383 * t414 + t568 * t556 - t560 * t723) * MDP(27) + (-pkin(4) * t371 - t383 * t416 + t556 * t723 + t568 * t560) * MDP(28) - (MDP(15) * t453 - MDP(16) * t570 - t402 * MDP(20) - t416 * MDP(24) + t414 * MDP(25) - t452 * MDP(26) + MDP(27) * t594 + t368 * MDP(28)) * t570; t416 * t414 * MDP(22) + (-t414 ^ 2 + t416 ^ 2) * MDP(23) + (t633 + t700) * MDP(24) + (-t614 + t699) * MDP(25) + t396 * MDP(26) + (t368 * t452 - t374 * t416 - g(1) * t426 - g(2) * t566 - g(3) * (-t498 * t556 - t640) + t617) * MDP(27) + (-t594 * t452 + t374 * t414 + g(1) * t567 - g(2) * t425 - g(3) * (-t498 * t560 + t641) - t595) * MDP(28) + (MDP(24) * t693 - MDP(25) * t416 - MDP(27) * t368 + MDP(28) * t594) * qJD(5);];
tau = t1;

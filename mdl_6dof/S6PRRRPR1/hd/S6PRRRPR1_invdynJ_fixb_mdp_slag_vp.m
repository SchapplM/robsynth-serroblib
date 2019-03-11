% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:03:36
% EndTime: 2019-03-08 23:03:46
% DurationCPUTime: 6.72s
% Computational Cost: add. (6498->505), mult. (15531->686), div. (0->0), fcn. (12484->18), ass. (0->238)
t574 = sin(qJ(4));
t575 = sin(qJ(3));
t578 = cos(qJ(4));
t579 = cos(qJ(3));
t513 = t574 * t579 + t575 * t578;
t702 = pkin(8) + pkin(9);
t645 = qJD(3) * t702;
t518 = t575 * t645;
t519 = t579 * t645;
t569 = sin(pkin(6));
t580 = cos(qJ(2));
t673 = t569 * t580;
t643 = qJD(1) * t673;
t527 = t702 * t575;
t528 = t702 * t579;
t663 = -t527 * t574 + t528 * t578;
t720 = -qJD(4) * t663 + t513 * t643 + t574 * t518 - t519 * t578;
t512 = t574 * t575 - t578 * t579;
t656 = qJD(4) * t578;
t657 = qJD(4) * t574;
t719 = -t512 * t643 + t518 * t578 + t519 * t574 + t527 * t656 + t528 * t657;
t563 = qJD(3) + qJD(4);
t577 = cos(qJ(6));
t573 = sin(qJ(6));
t658 = qJD(2) * t579;
t638 = t578 * t658;
t659 = qJD(2) * t575;
t640 = t574 * t659;
t505 = -t638 + t640;
t507 = -t574 * t658 - t578 * t659;
t567 = sin(pkin(12));
t570 = cos(pkin(12));
t605 = -t505 * t567 - t507 * t570;
t679 = t605 * t573;
t438 = -t563 * t577 + t679;
t624 = -t505 * t570 + t507 * t567;
t709 = qJD(6) - t624;
t718 = t438 * t709;
t440 = t563 * t573 + t577 * t605;
t717 = t440 * t709;
t474 = t563 * t513;
t716 = -qJ(5) * t474 - qJD(5) * t512 - t719;
t473 = t563 * t512;
t715 = qJ(5) * t473 - qJD(5) * t513 + t720;
t651 = qJD(2) * qJD(3);
t635 = t579 * t651;
t650 = qJDD(2) * t575;
t714 = t635 + t650;
t622 = t709 * t577;
t649 = qJDD(2) * t579;
t435 = qJD(4) * t638 - t563 * t640 + t574 * t649 + t578 * t714;
t612 = t574 * t650 - t578 * t649;
t436 = qJD(2) * t474 + t612;
t402 = -t435 * t567 - t436 * t570;
t400 = qJDD(6) - t402;
t684 = t400 * t573;
t713 = -t622 * t709 - t684;
t561 = qJDD(3) + qJDD(4);
t576 = sin(qJ(2));
t660 = qJD(1) * t576;
t644 = t569 * t660;
t632 = qJD(2) * t702 + t644;
t572 = cos(pkin(6));
t661 = qJD(1) * t572;
t478 = t575 * t661 + t579 * t632;
t648 = t572 * qJDD(1);
t533 = t579 * t648;
t652 = qJD(1) * qJD(2);
t485 = qJDD(2) * pkin(8) + (qJDD(1) * t576 + t580 * t652) * t569;
t630 = pkin(9) * qJDD(2) + t485;
t426 = qJDD(3) * pkin(3) - qJD(3) * t478 - t575 * t630 + t533;
t477 = -t575 * t632 + t579 * t661;
t429 = qJD(3) * t477 + t575 * t648 + t579 * t630;
t470 = t578 * t478;
t690 = qJD(3) * pkin(3);
t471 = t477 + t690;
t606 = -t471 * t574 - t470;
t590 = qJD(4) * t606 + t426 * t578 - t574 * t429;
t377 = pkin(4) * t561 - qJ(5) * t435 + qJD(5) * t507 + t590;
t703 = -(qJD(4) * t471 + t429) * t578 - t574 * t426 + t478 * t657;
t379 = -qJ(5) * t436 - qJD(5) * t505 - t703;
t367 = t377 * t570 - t379 * t567;
t365 = -pkin(5) * t561 - t367;
t568 = sin(pkin(11));
t571 = cos(pkin(11));
t671 = t572 * t576;
t499 = t568 * t580 + t571 * t671;
t501 = -t568 * t671 + t571 * t580;
t566 = qJ(3) + qJ(4);
t557 = pkin(12) + t566;
t546 = sin(t557);
t547 = cos(t557);
t675 = t569 * t576;
t676 = t569 * t571;
t677 = t568 * t569;
t712 = t365 + g(3) * (-t546 * t675 + t547 * t572) + g(2) * (-t499 * t546 - t547 * t676) + g(1) * (-t501 * t546 + t547 * t677);
t560 = t579 * pkin(3);
t693 = pkin(2) + t560;
t496 = t507 * qJ(5);
t468 = t574 * t478;
t626 = t471 * t578 - t468;
t421 = t496 + t626;
t398 = t577 * t400;
t655 = qJD(6) * t573;
t710 = t655 * t709 - t398;
t708 = pkin(4) * t436 + qJDD(5);
t558 = sin(t566);
t559 = cos(t566);
t707 = -g(3) * (-t558 * t675 + t559 * t572) - g(2) * (-t499 * t558 - t559 * t676) - g(1) * (-t501 * t558 + t559 * t677);
t368 = t377 * t567 + t379 * t570;
t366 = pkin(10) * t561 + t368;
t411 = pkin(4) * t563 + t421;
t689 = qJ(5) * t505;
t422 = -t606 - t689;
t682 = t422 * t567;
t388 = t411 * t570 - t682;
t386 = -pkin(5) * t563 - t388;
t497 = -qJD(2) * t693 - t643;
t457 = pkin(4) * t505 + qJD(5) + t497;
t397 = -pkin(5) * t624 - pkin(10) * t605 + t457;
t466 = t512 * t570 + t513 * t567;
t467 = -t512 * t567 + t513 * t570;
t617 = pkin(4) * t512 - t693;
t412 = pkin(5) * t466 - pkin(10) * t467 + t617;
t446 = -qJ(5) * t512 + t663;
t623 = -t527 * t578 - t528 * t574;
t601 = -qJ(5) * t513 + t623;
t416 = t446 * t570 + t567 * t601;
t431 = -t473 * t570 - t474 * t567;
t615 = g(1) * t501 + g(2) * t499;
t594 = -g(3) * t675 - t615;
t667 = t567 * t715 + t570 * t716;
t706 = -(qJD(6) * t397 + t366) * t466 + t365 * t467 + t386 * t431 + (-qJD(6) * t412 - t667) * t709 - t416 * t400 + t594;
t670 = t572 * t580;
t498 = t568 * t576 - t571 * t670;
t500 = t568 * t670 + t571 * t576;
t616 = g(1) * t500 + g(2) * t498;
t595 = g(3) * t673 - t616;
t705 = t412 * t400 - t595 * t547;
t637 = t576 * t652;
t529 = t569 * t637;
t581 = qJD(3) ^ 2;
t634 = qJDD(1) * t673;
t704 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t581 + t569 * (-g(3) * t580 + t637) - t529 + t616 + t634;
t700 = pkin(4) * t507;
t692 = pkin(3) * qJD(4);
t691 = qJD(2) * pkin(2);
t687 = t386 * t624;
t686 = t386 * t467;
t403 = t435 * t570 - t436 * t567;
t654 = qJD(6) * t577;
t647 = t403 * t577 + t561 * t573 + t563 * t654;
t393 = -t605 * t655 + t647;
t685 = t393 * t573;
t681 = t438 * t605;
t680 = t440 * t605;
t678 = t567 * t574;
t674 = t569 * t579;
t413 = t570 * t422;
t672 = t570 * t574;
t669 = qJDD(1) - g(3);
t668 = t567 * t716 - t570 * t715;
t389 = t411 * t567 + t413;
t666 = t477 * t578 - t468;
t424 = t496 + t666;
t625 = -t477 * t574 - t470;
t602 = t625 + t689;
t665 = -t424 * t567 + t570 * t602 + (t567 * t578 + t672) * t692;
t664 = -t570 * t424 - t567 * t602 + (t570 * t578 - t678) * t692;
t553 = pkin(3) * t578 + pkin(4);
t495 = pkin(3) * t672 + t553 * t567;
t523 = pkin(4) * t559 + t560;
t564 = t575 ^ 2;
t662 = -t579 ^ 2 + t564;
t556 = t575 * t690;
t641 = qJD(2) * t673;
t636 = t575 * t651;
t633 = pkin(4) * t474 + t556;
t631 = t569 * t669;
t387 = pkin(10) * t563 + t389;
t608 = t387 * t573 - t397 * t577;
t629 = t386 * t655 + t605 * t608;
t628 = t403 * t573 - t561 * t577;
t409 = pkin(5) * t605 - pkin(10) * t624 - t700;
t491 = pkin(10) + t495;
t555 = pkin(3) * t659;
t619 = qJD(6) * t491 + t409 + t555;
t548 = pkin(4) * t567 + pkin(10);
t618 = qJD(6) * t548 + t409;
t614 = g(1) * t568 - g(2) * t571;
t430 = -t473 * t567 + t474 * t570;
t613 = pkin(5) * t430 - pkin(10) * t431 + t633 - t644;
t610 = -t400 * t491 - t687;
t609 = -t400 * t548 - t687;
t374 = t387 * t577 + t397 * t573;
t607 = t388 * t624 + t389 * t605;
t503 = t572 * t579 - t575 * t675;
t504 = t572 * t575 + t576 * t674;
t447 = t503 * t578 - t504 * t574;
t448 = t503 * t574 + t504 * t578;
t604 = t573 * t624 * t709 - t710;
t494 = -pkin(3) * t678 + t553 * t570;
t600 = t431 * t577 - t467 * t655;
t598 = -t644 + t556;
t597 = pkin(3) * t636 - qJDD(2) * t693 + t529;
t596 = t374 * t605 + t386 * t654 + t573 * t712;
t521 = -t643 - t691;
t592 = -qJD(2) * t521 - t485 + t615;
t460 = t597 - t634;
t588 = -pkin(8) * qJDD(3) + (t521 + t643 - t691) * qJD(3);
t417 = t460 + t708;
t394 = qJD(6) * t440 + t628;
t587 = -t507 * t505 * MDP(12) - t709 * t605 * MDP(25) + ((t393 - t718) * t577 + (-t394 - t717) * t573) * MDP(22) + (t604 + t681) * MDP(24) + (-t680 - t713) * MDP(23) + (t440 * t622 + t685) * MDP(21) + (t505 * t563 + t435) * MDP(14) + (-t612 + (-qJD(2) * t513 - t507) * t563) * MDP(15) + (-t505 ^ 2 + t507 ^ 2) * MDP(13) + t561 * MDP(16);
t585 = -g(1) * (-t501 * t559 - t558 * t677) - g(2) * (-t499 * t559 + t558 * t676) - g(3) * (-t558 * t572 - t559 * t675) + t497 * t505 + t703;
t584 = t497 * t507 + t590 + t707;
t582 = qJD(2) ^ 2;
t562 = -qJ(5) - t702;
t549 = -pkin(4) * t570 - pkin(5);
t522 = -pkin(3) * t575 - pkin(4) * t558;
t517 = pkin(2) + t523;
t490 = -pkin(5) - t494;
t483 = t546 * t572 + t547 * t675;
t476 = -qJD(3) * t504 - t575 * t641;
t475 = qJD(3) * t503 + t579 * t641;
t456 = t501 * t547 + t546 * t677;
t454 = t499 * t547 - t546 * t676;
t419 = t447 * t567 + t448 * t570;
t418 = -t447 * t570 + t448 * t567;
t415 = t446 * t567 - t570 * t601;
t406 = -qJD(4) * t448 - t475 * t574 + t476 * t578;
t405 = qJD(4) * t447 + t475 * t578 + t476 * t574;
t391 = t421 * t570 - t682;
t390 = t421 * t567 + t413;
t381 = t405 * t570 + t406 * t567;
t380 = t405 * t567 - t406 * t570;
t375 = -pkin(5) * t402 - pkin(10) * t403 + t417;
t372 = t577 * t375;
t1 = [t669 * MDP(1) + (qJD(3) * t476 + qJDD(3) * t503) * MDP(10) + (-qJD(3) * t475 - qJDD(3) * t504) * MDP(11) + (t406 * t563 + t447 * t561) * MDP(17) + (-t405 * t563 - t448 * t561) * MDP(18) + (t380 * t605 + t381 * t624 + t402 * t419 + t403 * t418) * MDP(19) + (-t367 * t418 + t368 * t419 - t380 * t388 + t381 * t389 - g(3)) * MDP(20) + ((-t381 * t573 - t419 * t654) * t709 - t419 * t684 + t380 * t438 + t418 * t394) * MDP(26) + (-(t381 * t577 - t419 * t655) * t709 - t419 * t398 + t380 * t440 + t418 * t393) * MDP(27) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t579 + MDP(11) * t575 - MDP(3)) * t582 + (MDP(17) * t505 - MDP(18) * t507 + MDP(20) * t457 + (MDP(26) * t577 - MDP(27) * t573) * t709) * qJD(2)) * t576 + (qJDD(2) * MDP(3) - t582 * MDP(4) + (-t636 + t649) * MDP(10) - t714 * MDP(11) - t436 * MDP(17) - t435 * MDP(18) - t417 * MDP(20) + t710 * MDP(26) + (t654 * t709 + t684) * MDP(27)) * t580) * t569; qJDD(2) * MDP(2) + (t669 * t673 + t616) * MDP(3) + (-t576 * t631 + t615) * MDP(4) + (qJDD(2) * t564 + 0.2e1 * t575 * t635) * MDP(5) + 0.2e1 * (t575 * t649 - t651 * t662) * MDP(6) + (qJDD(3) * t575 + t579 * t581) * MDP(7) + (qJDD(3) * t579 - t575 * t581) * MDP(8) + (t588 * t575 + t579 * t704) * MDP(10) + (-t575 * t704 + t588 * t579) * MDP(11) + (t435 * t513 + t473 * t507) * MDP(12) + (-t435 * t512 - t436 * t513 + t473 * t505 + t474 * t507) * MDP(13) + (-t473 * t563 + t513 * t561) * MDP(14) + (-t474 * t563 - t512 * t561) * MDP(15) + (-t693 * t436 + t460 * t512 + t497 * t474 + t598 * t505 - t595 * t559 + t623 * t561 + t563 * t720) * MDP(17) + (-t693 * t435 + t460 * t513 - t497 * t473 - t598 * t507 + t595 * t558 - t663 * t561 + t563 * t719) * MDP(18) + (-t367 * t467 - t368 * t466 - t388 * t431 - t389 * t430 + t402 * t416 + t403 * t415 + t605 * t668 + t624 * t667 + t594) * MDP(19) + (t368 * t416 - t367 * t415 + t417 * t617 + t457 * t633 - g(1) * (-t500 * t517 - t501 * t562) - g(2) * (-t498 * t517 - t499 * t562) + t667 * t389 - t668 * t388 + (-t457 * t660 - g(3) * (t517 * t580 - t562 * t576)) * t569) * MDP(20) + (t393 * t467 * t577 + t440 * t600) * MDP(21) + ((-t438 * t577 - t440 * t573) * t431 + (-t685 - t394 * t577 + (t438 * t573 - t440 * t577) * qJD(6)) * t467) * MDP(22) + (t393 * t466 + t398 * t467 + t430 * t440 + t600 * t709) * MDP(23) + (-t467 * t684 - t394 * t466 - t430 * t438 + (-t431 * t573 - t467 * t654) * t709) * MDP(24) + (t400 * t466 + t430 * t709) * MDP(25) + (t372 * t466 - t608 * t430 + t415 * t394 + t668 * t438 + (t613 * t709 + (-t387 * t466 - t416 * t709 + t686) * qJD(6) + t705) * t577 + t706 * t573) * MDP(26) + (-t374 * t430 + t415 * t393 + t668 * t440 + (-(-qJD(6) * t387 + t375) * t466 - qJD(6) * t686 + (qJD(6) * t416 - t613) * t709 - t705) * t573 + t706 * t577) * MDP(27); t587 + (g(3) * t504 + (t569 * t614 - t648) * t575 + t592 * t579) * MDP(11) + (t402 * t495 - t403 * t494 + t605 * t665 + t624 * t664 + t607) * MDP(19) + MDP(8) * t649 + (-g(3) * t503 + t575 * t592 - t614 * t674 + t533) * MDP(10) + (t490 * t394 + t665 * t438 + (-t664 * t709 + t610) * t573 + (-t619 * t709 - t712) * t577 + t629) * MDP(26) + (t368 * t495 + t367 * t494 - t457 * (t555 - t700) - g(1) * (t501 * t522 + t523 * t677) - g(2) * (t499 * t522 - t523 * t676) - g(3) * (t522 * t675 + t523 * t572) + t664 * t389 - t665 * t388) * MDP(20) + (t490 * t393 + t610 * t577 + t665 * t440 + (t573 * t619 - t577 * t664) * t709 + t596) * MDP(27) + (-t625 * t563 + (-t505 * t659 + t561 * t578 - t563 * t657) * pkin(3) + t584) * MDP(17) + (t666 * t563 + (t507 * t659 - t561 * t574 - t563 * t656) * pkin(3) + t585) * MDP(18) + qJDD(3) * MDP(9) + MDP(7) * t650 + (-MDP(5) * t575 * t579 + MDP(6) * t662) * t582; t587 + (t388 * t390 - t389 * t391 + (t367 * t570 + t368 * t567 + t457 * t507 + t707) * pkin(4)) * MDP(20) + (-t390 * t440 + t549 * t393 + t609 * t577 + (t391 * t577 + t573 * t618) * t709 + t596) * MDP(27) + (-t390 * t438 + t549 * t394 + (t391 * t709 + t609) * t573 + (-t618 * t709 - t712) * t577 + t629) * MDP(26) + (-t390 * t605 - t391 * t624 + (t402 * t567 - t403 * t570) * pkin(4) + t607) * MDP(19) + (-t563 * t606 + t584) * MDP(17) + (t563 * t626 + t585) * MDP(18); (-t605 ^ 2 - t624 ^ 2) * MDP(19) + (t388 * t605 - t389 * t624 - t580 * t631 + t597 - t616 + t708) * MDP(20) + (t604 - t681) * MDP(26) + (-t680 + t713) * MDP(27); t440 * t438 * MDP(21) + (-t438 ^ 2 + t440 ^ 2) * MDP(22) + (t647 + t718) * MDP(23) + (-t628 + t717) * MDP(24) + t400 * MDP(25) + (-t573 * t366 + t372 + t374 * t709 - t386 * t440 - g(1) * (-t456 * t573 + t500 * t577) - g(2) * (-t454 * t573 + t498 * t577) - g(3) * (-t483 * t573 - t577 * t673)) * MDP(26) + (-t577 * t366 - t573 * t375 - t608 * t709 + t386 * t438 - g(1) * (-t456 * t577 - t500 * t573) - g(2) * (-t454 * t577 - t498 * t573) - g(3) * (-t483 * t577 + t573 * t673)) * MDP(27) + (-MDP(23) * t679 - MDP(24) * t440 - MDP(26) * t374 + MDP(27) * t608) * qJD(6);];
tau  = t1;

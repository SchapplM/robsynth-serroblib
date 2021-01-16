% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:13:08
% EndTime: 2021-01-16 03:13:23
% DurationCPUTime: 9.48s
% Computational Cost: add. (3473->567), mult. (7558->740), div. (0->0), fcn. (5385->10), ass. (0->273)
t560 = sin(qJ(5));
t563 = cos(qJ(5));
t677 = qJD(3) * t560;
t564 = cos(qJ(3));
t681 = qJD(2) * t564;
t501 = t563 * t681 + t677;
t672 = qJD(5) * t501;
t561 = sin(qJ(3));
t664 = qJD(2) * qJD(3);
t638 = t561 * t664;
t660 = qJDD(2) * t564;
t766 = t638 - t660;
t431 = -t563 * qJDD(3) - t766 * t560 + t672;
t682 = qJD(2) * t561;
t537 = qJD(5) + t682;
t719 = t501 * t537;
t775 = t431 + t719;
t557 = sin(pkin(6));
t562 = sin(qJ(2));
t565 = cos(qJ(2));
t665 = qJD(1) * qJD(2);
t640 = t565 * t665;
t559 = cos(pkin(6));
t678 = qJD(3) * t559;
t726 = qJDD(2) * pkin(8);
t774 = t726 + (qJDD(1) * t562 + t640) * t557 + qJD(1) * t678;
t566 = pkin(4) + pkin(8);
t642 = t560 * t681;
t675 = qJD(3) * t563;
t503 = -t642 + t675;
t567 = -pkin(3) - pkin(9);
t686 = qJD(1) * t557;
t650 = t562 * t686;
t510 = qJD(2) * pkin(8) + t650;
t685 = qJD(1) * t559;
t464 = t561 * t510 - t564 * t685;
t765 = qJD(4) + t464;
t666 = pkin(4) * t682 + t765;
t433 = qJD(3) * t567 + t666;
t731 = qJ(4) * t561;
t497 = t564 * t567 - pkin(2) - t731;
t649 = t565 * t686;
t454 = qJD(2) * t497 - t649;
t407 = t433 * t560 + t454 * t563;
t663 = qJDD(1) * t559;
t636 = t564 * t663;
t674 = qJD(3) * t564;
t746 = t510 * t674 + t561 * t774;
t582 = qJDD(4) - t636 + t746;
t637 = t564 * t664;
t661 = qJDD(2) * t561;
t592 = t637 + t661;
t402 = pkin(4) * t592 + qJDD(3) * t567 + t582;
t730 = qJ(4) * t564;
t615 = pkin(9) * t561 - t730;
t673 = qJD(4) * t561;
t581 = qJD(3) * t615 - t673;
t641 = t562 * t665;
t713 = t557 * t565;
t614 = -qJDD(1) * t713 + t557 * t641;
t603 = pkin(3) * t638 + t614;
t410 = qJD(2) * t581 + qJDD(2) * t497 + t603;
t631 = t563 * t402 - t410 * t560;
t577 = -t407 * qJD(5) + t631;
t729 = qJ(6) * t431;
t500 = qJDD(5) + t592;
t741 = pkin(5) * t500;
t390 = -qJD(6) * t503 + t577 + t729 + t741;
t399 = -qJ(6) * t501 + t407;
t724 = t399 * t537;
t773 = t390 + t724;
t670 = qJD(5) * t563;
t671 = qJD(5) * t560;
t627 = t560 * t402 + t563 * t410 + t433 * t670 - t454 * t671;
t432 = -qJD(5) * t642 + qJDD(3) * t560 + (qJD(3) * qJD(5) - t766) * t563;
t728 = qJ(6) * t432;
t391 = -qJD(6) * t501 + t627 - t728;
t406 = t563 * t433 - t454 * t560;
t398 = -qJ(6) * t503 + t406;
t397 = pkin(5) * t537 + t398;
t772 = -t397 * t537 + t391;
t718 = t503 * t537;
t771 = -t432 + t718;
t770 = t432 + t718;
t700 = t563 * t565;
t470 = (-t560 * t562 + t561 * t700) * t557;
t509 = t566 * t674;
t769 = -qJD(1) * t470 + t563 * t509;
t676 = qJD(3) * t561;
t543 = pkin(3) * t676;
t467 = t543 + t581;
t703 = t562 * t563;
t704 = t561 * t565;
t471 = (t560 * t704 + t703) * t557;
t518 = t566 * t561;
t768 = qJD(1) * t471 - t563 * t467 + t497 * t671 - t560 * t509 - t518 * t670;
t556 = sin(pkin(10));
t710 = t559 * t562;
t558 = cos(pkin(10));
t759 = t558 * t565;
t482 = t556 * t710 - t759;
t763 = t556 * t565;
t484 = t558 * t710 + t763;
t714 = t557 * t564;
t706 = t561 * t562;
t656 = t557 * t706;
t490 = -t559 * t564 + t656;
t736 = g(3) * t490;
t767 = g(2) * (t484 * t561 + t558 * t714) + t736 - g(1) * (t482 * t561 + t556 * t714);
t657 = MDP(12) * qJDD(2);
t540 = pkin(5) * t560 + qJ(4);
t756 = qJ(6) - t567;
t611 = t540 * t561 + t564 * t756;
t762 = t557 * t611;
t618 = pkin(3) * t564 + t731;
t761 = t557 * t618;
t760 = t558 * t562;
t740 = pkin(5) * t563;
t635 = t566 + t740;
t758 = t562 * t635;
t757 = t565 * (pkin(2) + t611);
t632 = qJ(6) * t564 - t497;
t667 = qJD(6) * t564;
t695 = pkin(5) * t674 + t632 * t670 + (-qJ(6) * t676 - qJD(5) * t518 - t467 + t667) * t560 + t769;
t669 = qJD(5) * t564;
t643 = t560 * t669;
t755 = -t563 * t667 + (t675 * t561 + t643) * qJ(6) - t768;
t754 = t397 - t398;
t513 = t756 * t563;
t465 = t564 * t510 + t561 * t685;
t453 = pkin(4) * t681 + t465;
t544 = pkin(3) * t682;
t478 = qJD(2) * t615 + t544;
t692 = t560 * t453 + t563 * t478;
t693 = qJ(6) * t563 * t682 + qJD(5) * t513 + qJD(6) * t560 + t692;
t437 = t563 * t453;
t708 = t560 * t561;
t691 = -qJD(6) * t563 + t671 * t756 + t478 * t560 - t437 - (pkin(5) * t564 - qJ(6) * t708) * qJD(2);
t690 = t563 * t497 + t560 * t518;
t595 = -qJ(4) * t674 - t673;
t753 = t543 + t595 - t650;
t651 = -pkin(4) - t740;
t689 = pkin(5) * t670 - t651 * t682 + t765;
t553 = qJD(3) * qJ(4);
t458 = -t553 - t465;
t457 = -qJD(3) * pkin(3) + t765;
t481 = t563 * t500;
t752 = -t537 * t671 + t481;
t551 = qJDD(3) * qJ(4);
t552 = qJD(3) * qJD(4);
t626 = t510 * t676 - t561 * t663 - t564 * t774;
t411 = -t551 - t552 + t626;
t725 = qJDD(3) * pkin(3);
t412 = t582 - t725;
t751 = -t411 * t564 + t412 * t561;
t702 = t562 * t564;
t653 = t559 * t702;
t715 = t557 * t561;
t492 = t653 - t715;
t699 = t564 * t565;
t491 = t557 * t702 + t559 * t561;
t735 = g(3) * t491;
t750 = -t735 - g(2) * (t492 * t558 + t556 * t699);
t749 = -t465 * qJD(3) + t746 - t767;
t748 = -MDP(10) + MDP(13);
t747 = MDP(11) - MDP(14);
t435 = t553 + t453;
t745 = t435 * t537 + t567 * t500;
t709 = t559 * t565;
t485 = t556 * t709 + t760;
t600 = -t556 * t562 + t558 * t709;
t621 = g(1) * t485 - g(2) * t600;
t568 = qJD(3) ^ 2;
t739 = pkin(8) * t568;
t744 = 0.2e1 * qJDD(2) * pkin(2) + t557 * (-g(3) * t565 + t641) - t614 + t621 - t739;
t514 = pkin(2) + t618;
t662 = qJDD(2) * t514;
t416 = qJD(2) * t595 + t603 - t662;
t585 = -g(3) * t713 + t621;
t743 = -qJD(2) * t753 - t416 + t585 + t662 - t739;
t742 = t503 ^ 2;
t734 = g(3) * t557;
t733 = qJD(2) * pkin(2);
t732 = pkin(8) * qJDD(3);
t721 = t431 * t563;
t720 = t500 * t560;
t717 = pkin(2) * t760;
t716 = t556 * t559;
t712 = t558 * t559;
t711 = t559 * t560;
t707 = t560 * t565;
t705 = t561 * t563;
t701 = t563 * t564;
t696 = qJDD(1) - g(3);
t519 = t566 * t564;
t554 = t561 ^ 2;
t555 = t564 ^ 2;
t688 = t554 - t555;
t687 = t554 + t555;
t684 = qJD(2) * t514;
t683 = qJD(2) * t557;
t680 = qJD(3) * t501;
t679 = qJD(3) * t503;
t668 = qJD(5) * t567;
t659 = MDP(21) + MDP(23);
t658 = MDP(22) + MDP(24);
t655 = t557 * t701;
t654 = t559 * t705;
t569 = qJD(2) ^ 2;
t652 = t561 * t564 * t569;
t648 = t562 * t683;
t647 = t565 * t683;
t634 = -pkin(5) * t501 - qJD(6);
t421 = t435 - t634;
t646 = t421 * t671;
t645 = t421 * t670;
t489 = t559 * t706 + t714;
t630 = t489 * t558 + t556 * t704;
t629 = -t489 * t556 + t558 * t704;
t628 = pkin(8) - t651;
t625 = t501 * t649;
t624 = t503 * t649;
t623 = t564 * t649;
t622 = -g(1) * t482 + g(2) * t484;
t617 = pkin(3) * t561 - t730;
t616 = pkin(8) * t562 + t514 * t565;
t612 = t540 * t564 - t561 * t756;
t610 = t537 ^ 2;
t451 = -t490 * t560 + t557 * t700;
t450 = t490 * t563 + t557 * t707;
t599 = -t537 * t670 - t720;
t598 = t562 * t617;
t597 = t565 * t617;
t596 = t618 * t559;
t594 = t562 * t612;
t593 = t565 * t612;
t468 = t600 * t561;
t469 = t485 * t561;
t589 = -g(1) * (-t469 * t560 - t482 * t563) - g(2) * (t468 * t560 + t484 * t563) - g(3) * t471;
t588 = -g(1) * (-t469 * t563 + t482 * t560) - g(2) * (t468 * t563 - t484 * t560) - g(3) * t470;
t587 = -g(1) * t629 - g(2) * t630 - t736;
t526 = t558 * t699;
t586 = -g(1) * (-t492 * t556 + t526) + t750;
t403 = -t766 * pkin(4) - t411;
t394 = t432 * pkin(5) + qJDD(6) + t403;
t584 = t394 + t586;
t583 = t403 + t586;
t578 = -g(1) * (-t485 * t563 - t560 * t629) - g(2) * (-t560 * t630 + t563 * t600) - g(3) * t451 - t627;
t511 = -t649 - t733;
t576 = -t732 + (t511 + t649 - t733) * qJD(3);
t466 = -t649 - t684;
t575 = t732 + (-t466 - t649 + t684) * qJD(3);
t574 = -t562 * t611 + t565 * t635;
t516 = t556 * t715;
t573 = g(1) * (t482 * t564 - t516) - g(2) * (t484 * t564 - t558 * t715) + qJD(3) * t464 - t626 - t735;
t570 = -g(1) * (-t485 * t560 + t563 * t629) - g(2) * (t560 * t600 + t563 * t630) - g(3) * t450 + t577;
t546 = t556 * pkin(2);
t534 = pkin(2) * t712;
t512 = t756 * t560;
t508 = t566 * t676;
t507 = -qJ(4) * t681 + t544;
t506 = t563 * t518;
t499 = t501 ^ 2;
t488 = pkin(5) * t701 + t519;
t456 = t466 * t682;
t455 = -pkin(5) * t643 - t628 * t676;
t449 = qJD(3) * t491 + t561 * t647;
t448 = -qJD(3) * t656 + (t647 + t678) * t564;
t429 = -qJ(6) * t701 + t690;
t426 = pkin(5) * t561 + t560 * t632 + t506;
t405 = qJD(5) * t450 + t449 * t560 + t563 * t648;
t404 = qJD(5) * t451 + t449 * t563 - t560 * t648;
t1 = [t696 * MDP(1) + (-t411 * t491 + t412 * t490 - t448 * t458 + t449 * t457 - g(3)) * MDP(15) + (-t404 * t503 - t405 * t501 + t431 * t450 + t432 * t451) * MDP(25) + (t390 * t450 - t391 * t451 + t394 * t491 + t397 * t404 + t399 * t405 + t421 * t448 - g(3)) * MDP(26) + t658 * (-t405 * t537 - t431 * t491 + t448 * t503 + t451 * t500) + t659 * (t404 * t537 + t432 * t491 + t448 * t501 + t450 * t500) + (t490 * t561 + t491 * t564) * t657 + (t448 * t564 + t449 * t561 + (t490 * t564 - t491 * t561) * qJD(3)) * MDP(12) * qJD(2) + ((qJD(2) * t466 * MDP(15) - qJDD(2) * MDP(4) + (t561 * t747 + t564 * t748 - MDP(3)) * t569) * t562 + (-t416 * MDP(15) + qJDD(2) * MDP(3) - t569 * MDP(4) - t747 * t592 + t748 * t766) * t565) * t557 - t747 * (qJD(3) * t448 + qJDD(3) * t491) + t748 * (qJD(3) * t449 + qJDD(3) * t490); qJDD(2) * MDP(2) + (t696 * t713 + t621) * MDP(3) + (-t557 * t562 * t696 + t622) * MDP(4) + (qJDD(2) * t554 + 0.2e1 * t561 * t637) * MDP(5) + 0.2e1 * (t561 * t660 - t664 * t688) * MDP(6) + (qJDD(3) * t561 + t564 * t568) * MDP(7) + (qJDD(3) * t564 - t561 * t568) * MDP(8) + (t576 * t561 + t564 * t744) * MDP(10) + (-t561 * t744 + t576 * t564) * MDP(11) + ((t457 * t564 + t458 * t561) * qJD(3) + t687 * t726 + (-g(3) * t562 - t640 * t687) * t557 - t622 + t751) * MDP(12) + (t575 * t561 - t564 * t743) * MDP(13) + (t561 * t743 + t575 * t564) * MDP(14) + (-t416 * t514 + t458 * t623 - t457 * t561 * t649 - g(1) * (-t616 * t716 - t618 * t760 - t717) - g(2) * (-(t556 * t618 + t546) * t562 + (t558 * t596 + t534) * t565) - t616 * t734 + t753 * t466 + (t458 * t676 + t457 * t674 - g(1) * t759 - g(2) * (t562 * t712 + t763) + t751) * pkin(8)) * MDP(15) + (t431 * t560 * t564 + (t560 * t676 - t563 * t669) * t503) * MDP(16) + ((-t501 * t560 + t503 * t563) * t676 + (t721 + t432 * t560 + (t501 * t563 + t503 * t560) * qJD(5)) * t564) * MDP(17) + ((t537 * t677 - t431) * t561 + (t599 + t679) * t564) * MDP(18) + ((t537 * t675 - t432) * t561 + (-t680 - t752) * t564) * MDP(19) + (t500 * t561 + t537 * t674) * MDP(20) + ((-t497 * t560 + t506) * t500 - t508 * t501 + t519 * t432 + (-t435 * t675 + t631) * t561 + (-t467 * t560 + t769) * t537 + (-t407 * t561 - t537 * t690) * qJD(5) + (qJD(3) * t406 + t403 * t563 - t435 * t671 - t625) * t564 + t589) * MDP(21) + (-t690 * t500 - t508 * t503 - t519 * t431 + (t435 * t677 - t627) * t561 + t768 * t537 + (-qJD(3) * t407 - t403 * t560 - t435 * t670 - t624) * t564 + t588) * MDP(22) + (t426 * t500 + t432 * t488 + t455 * t501 + (-t421 * t675 + t390) * t561 + t695 * t537 + (qJD(3) * t397 + t394 * t563 - t625 - t646) * t564 + t589) * MDP(23) + (-t429 * t500 - t431 * t488 + t455 * t503 + (t421 * t677 - t391) * t561 - t755 * t537 + (-qJD(3) * t399 - t394 * t560 - t624 - t645) * t564 + t588) * MDP(24) + (t426 * t431 - t429 * t432 - t695 * t503 - t755 * t501 + (-t397 * t560 + t399 * t563) * t676 + (t390 * t560 - t391 * t563 + (t397 * t563 + t399 * t560) * qJD(5) + t585) * t564) * MDP(25) + (t391 * t429 + t390 * t426 + t394 * t488 - g(1) * (-t717 + t574 * t558 + (-t757 - t758) * t716) - g(2) * (t534 * t565 - t546 * t562 + (t565 * t611 + t758) * t712 + t574 * t556) - (t562 * t628 + t757) * t734 + (t455 - t623) * t421 + t755 * t399 + t695 * t397) * MDP(26); -MDP(5) * t652 + t688 * t569 * MDP(6) + MDP(7) * t661 + MDP(8) * t660 + qJDD(3) * MDP(9) + (-t511 * t682 + t636 - t749) * MDP(10) + (-t511 * t681 - t573) * MDP(11) - t617 * t657 + (-0.2e1 * t725 + qJDD(4) + t456 + (-qJD(2) * t507 - t663) * t564 + t749) * MDP(13) + (0.2e1 * t551 + 0.2e1 * t552 + (t466 * t564 + t507 * t561) * qJD(2) + t573) * MDP(14) + (-t411 * qJ(4) - t412 * pkin(3) - t466 * t507 - t457 * t465 - g(1) * (-t558 * t597 + (t617 * t710 + t761) * t556) - g(2) * (-t556 * t597 + (-t559 * t598 - t761) * t558) - g(3) * (-t557 * t598 + t596) - t765 * t458) * MDP(15) + (-t560 * t718 - t721) * MDP(16) + (t560 * t775 - t770 * t563) * MDP(17) + ((-t503 * t564 - t537 * t708) * qJD(2) + t752) * MDP(18) + ((t501 * t564 - t537 * t705) * qJD(2) + t599) * MDP(19) - t537 * MDP(20) * t681 + (-t406 * t681 + qJ(4) * t432 - t437 * t537 + t666 * t501 + t745 * t563 + ((t478 - t668) * t537 + t583) * t560) * MDP(21) + (-qJ(4) * t431 + t692 * t537 + t407 * t681 + t666 * t503 - t745 * t560 + (-t537 * t668 + t583) * t563) * MDP(22) + (t645 + t432 * t540 - t500 * t513 + t691 * t537 + t689 * t501 + (-t397 * t564 + t421 * t705) * qJD(2) + t584 * t560) * MDP(23) + (-t646 - t431 * t540 + t500 * t512 + t693 * t537 + t689 * t503 + (t399 * t564 - t421 * t708) * qJD(2) + t584 * t563) * MDP(24) + (-t431 * t513 + t432 * t512 + t693 * t501 - t691 * t503 - t772 * t560 - t773 * t563 + t767) * MDP(25) + (-t391 * t512 - t390 * t513 + t394 * t540 - g(1) * (t558 * t593 + (-t612 * t710 + t762) * t556) - g(2) * (t556 * t593 + (t559 * t594 - t762) * t558) - g(3) * (t557 * t594 + t559 * t611) + t689 * t421 - t693 * t399 + t691 * t397) * MDP(26); t561 * t657 + (qJDD(3) + t652) * MDP(13) + (-t554 * t569 - t568) * MDP(14) + (qJD(3) * t458 + t412 + t456 + t587) * MDP(15) + (-qJD(3) * t421 + t587) * MDP(26) + t658 * (-t563 * t610 - t679 - t720) + t659 * (-t560 * t610 + t481 - t680) + ((-t501 * t682 + t431 - t672) * MDP(25) + t773 * MDP(26)) * t563 + (t771 * MDP(25) + t772 * MDP(26)) * t560; t503 * t501 * MDP(16) + (-t499 + t742) * MDP(17) + (-t431 + t719) * MDP(18) + t771 * MDP(19) + t500 * MDP(20) + (t407 * t537 - t435 * t503 + t570) * MDP(21) + (t406 * t537 + t435 * t501 + t578) * MDP(22) + (0.2e1 * t741 + t729 + t724 + (-t421 + t634) * t503 + t570) * MDP(23) + (-pkin(5) * t742 + t728 + t398 * t537 + (qJD(6) + t421) * t501 + t578) * MDP(24) + (pkin(5) * t431 - t501 * t754) * MDP(25) + (t754 * t399 + (t390 - t421 * t503 - g(1) * ((-t556 * t711 + t558 * t705) * t565 + (-t556 * t654 - t558 * t560) * t562 - t556 * t655) - g(2) * ((t556 * t705 + t558 * t711) * t565 + (-t556 * t560 + t558 * t654) * t562 + t558 * t655) - g(3) * (-t559 * t701 + (t561 * t703 + t707) * t557)) * pkin(5)) * MDP(26); t770 * MDP(23) - t775 * MDP(24) + (-t499 - t742) * MDP(25) + (-g(1) * (-t556 * t653 + t516 + t526) + t397 * t503 + t399 * t501 + t394 + t750) * MDP(26);];
tau = t1;

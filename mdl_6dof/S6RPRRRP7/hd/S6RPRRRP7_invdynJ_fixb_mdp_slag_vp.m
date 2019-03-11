% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:21:17
% EndTime: 2019-03-09 06:21:33
% DurationCPUTime: 11.58s
% Computational Cost: add. (10819->603), mult. (25114->728), div. (0->0), fcn. (19416->14), ass. (0->254)
t617 = sin(qJ(4));
t620 = cos(qJ(4));
t613 = sin(pkin(10));
t614 = cos(pkin(10));
t618 = sin(qJ(3));
t764 = cos(qJ(3));
t571 = t764 * t613 + t618 * t614;
t776 = t571 * qJD(1);
t529 = qJD(3) * t620 - t617 * t776;
t530 = qJD(3) * t617 + t620 * t776;
t616 = sin(qJ(5));
t763 = cos(qJ(5));
t473 = -t763 * t529 + t530 * t616;
t649 = t616 * t529 + t530 * t763;
t807 = t473 * t649;
t691 = t764 * t614;
t652 = -t618 * t613 + t691;
t802 = t652 * qJD(1);
t512 = pkin(3) * t776 - pkin(8) * t802;
t498 = t620 * t512;
t622 = -pkin(9) - pkin(8);
t695 = qJD(4) * t622;
t741 = t802 * t620;
t754 = pkin(7) + qJ(2);
t583 = t754 * t613;
t572 = qJD(1) * t583;
t584 = t754 * t614;
t573 = qJD(1) * t584;
t782 = -t764 * t572 - t618 * t573;
t806 = pkin(4) * t776 - pkin(9) * t741 - t617 * t782 - t620 * t695 + t498;
t715 = t617 * t512 + t620 * t782;
t742 = t802 * t617;
t805 = -pkin(9) * t742 - t617 * t695 + t715;
t564 = t571 * qJD(3);
t701 = qJDD(1) * t613;
t662 = -qJDD(1) * t691 + t618 * t701;
t515 = qJD(1) * t564 + t662;
t509 = qJDD(4) + t515;
t506 = qJDD(5) + t509;
t552 = qJD(4) - t802;
t548 = qJD(5) + t552;
t729 = t616 * t617;
t647 = t763 * t620 - t729;
t689 = t763 * t617;
t575 = t616 * t620 + t689;
t774 = qJD(4) + qJD(5);
t521 = t774 * t575;
t712 = -t575 * t802 + t521;
t665 = t647 * t506 - t548 * t712;
t681 = t763 * qJD(5);
t780 = t763 * qJD(4) + t681;
t713 = t620 * t780 - t647 * t802 - t729 * t774;
t664 = t575 * t506 + t548 * t713;
t707 = qJD(4) * t617;
t804 = t707 - t742;
t638 = t571 * qJDD(1);
t626 = (qJD(3) * (qJD(4) + t802) + t638) * t620 + (-qJD(4) * t776 + qJDD(3)) * t617;
t706 = qJD(4) * t620;
t563 = t652 * qJD(3);
t627 = qJD(1) * t563 + t638;
t797 = qJD(4) * qJD(3) + t627;
t696 = t617 * t797 + t776 * t706;
t658 = qJDD(3) * t620 - t696;
t705 = qJD(5) * t616;
t427 = -t529 * t681 + t530 * t705 - t616 * t658 - t763 * t626;
t413 = t473 * t548 - t427;
t428 = t529 * t705 + t530 * t681 + t616 * t626 - t763 * t658;
t765 = t649 ^ 2;
t803 = t506 * MDP(26) + (t548 * t649 - t428) * MDP(25) + MDP(22) * t807 + (-t473 ^ 2 + t765) * MDP(23) + t413 * MDP(24);
t510 = -qJD(3) * pkin(3) - t782;
t469 = -t529 * pkin(4) + t510;
t431 = t473 * pkin(5) - qJ(6) * t649 + t469;
t801 = t431 * t473;
t800 = t469 * t473;
t753 = qJDD(1) * pkin(1);
t603 = qJDD(2) - t753;
t619 = sin(qJ(1));
t621 = cos(qJ(1));
t781 = g(1) * t619 - g(2) * t621;
t661 = -t603 + t781;
t796 = t563 * t617 + t571 * t706;
t611 = pkin(10) + qJ(3);
t604 = sin(t611);
t733 = t604 * t621;
t734 = t604 * t619;
t795 = -g(1) * t733 - g(2) * t734;
t585 = t622 * t617;
t586 = t622 * t620;
t528 = t616 * t585 - t586 * t763;
t612 = qJ(4) + qJ(5);
t606 = sin(t612);
t668 = g(1) * t621 + g(2) * t619;
t605 = cos(t611);
t756 = g(3) * t605;
t636 = t604 * t668 - t756;
t794 = t528 * t506 + t606 * t636;
t793 = qJD(3) * t776;
t443 = pkin(5) * t649 + qJ(6) * t473;
t791 = -t662 - t793;
t747 = t509 * t617;
t598 = pkin(2) * t614 + pkin(1);
t514 = -pkin(3) * t652 - pkin(8) * t571 - t598;
t504 = t620 * t514;
t526 = -t618 * t583 + t584 * t764;
t737 = t571 * t620;
t450 = -pkin(4) * t652 - pkin(9) * t737 - t526 * t617 + t504;
t519 = t620 * t526;
t714 = t617 * t514 + t519;
t738 = t571 * t617;
t457 = -pkin(9) * t738 + t714;
t787 = t616 * t450 + t763 * t457;
t648 = t585 * t763 + t616 * t586;
t786 = -qJD(5) * t648 + t616 * t806 + t763 * t805;
t785 = -qJD(5) * t528 + t616 * t805 - t763 * t806;
t518 = -t618 * t572 + t764 * t573;
t672 = pkin(4) * t804 - t518;
t494 = t506 * qJ(6);
t539 = t548 * qJD(6);
t784 = t494 + t539;
t495 = t620 * t509;
t783 = t552 * t707 - t495;
t686 = t571 * t707;
t739 = t563 * t620;
t779 = -t686 + t739;
t778 = -t764 * t583 - t618 * t584;
t722 = t620 * t621;
t728 = t617 * t619;
t553 = t605 * t728 + t722;
t723 = t619 * t620;
t727 = t617 * t621;
t555 = -t605 * t727 + t723;
t777 = -g(1) * t555 + g(2) * t553;
t775 = qJ(2) * qJDD(1);
t500 = t506 * pkin(5);
t773 = t500 - qJDD(6);
t607 = cos(t612);
t731 = t607 * t621;
t732 = t606 * t619;
t535 = t605 * t732 + t731;
t721 = t621 * t606;
t724 = t619 * t607;
t537 = t605 * t721 - t724;
t703 = qJD(1) * qJD(2);
t767 = qJDD(1) * t754 + t703;
t545 = t767 * t613;
t546 = t767 * t614;
t654 = -t618 * t545 + t546 * t764;
t461 = qJDD(3) * pkin(8) + qJD(3) * t782 + t654;
t580 = -qJD(1) * t598 + qJD(2);
t486 = -pkin(3) * t802 - pkin(8) * t776 + t580;
t700 = qJDD(1) * t614;
t465 = -pkin(2) * t700 + t515 * pkin(3) - pkin(8) * t627 + t603;
t511 = qJD(3) * pkin(8) + t518;
t663 = t620 * t465 - t511 * t706;
t409 = t509 * pkin(4) - pkin(9) * t626 - t617 * t461 - t486 * t707 + t663;
t644 = t620 * t461 + t617 * t465 + t486 * t706 - t511 * t707;
t412 = pkin(9) * t658 + t644;
t454 = t620 * t486 - t511 * t617;
t445 = -pkin(9) * t530 + t454;
t437 = pkin(4) * t552 + t445;
t455 = t486 * t617 + t511 * t620;
t446 = pkin(9) * t529 + t455;
t675 = -t763 * t409 + t616 * t412 + t437 * t705 + t446 * t681;
t736 = t604 * t606;
t634 = g(1) * t537 + g(2) * t535 + g(3) * t736 - t675;
t629 = t431 * t649 - t634 - t773;
t772 = -t469 * t649 + t634;
t771 = -t552 ^ 2 * t620 - t747;
t757 = g(3) * t604;
t770 = t668 * t605 + t757;
t769 = -t427 * t647 - t649 * t712;
t487 = t652 * qJD(2) + qJD(3) * t778;
t513 = pkin(3) * t564 - pkin(8) * t563;
t499 = t620 * t513;
t424 = -pkin(9) * t739 + pkin(4) * t564 - t487 * t617 + t499 + (-t519 + (pkin(9) * t571 - t514) * t617) * qJD(4);
t643 = t620 * t487 + t617 * t513 + t514 * t706 - t526 * t707;
t430 = -pkin(9) * t796 + t643;
t768 = -qJD(5) * t787 + t424 * t763 - t616 * t430;
t755 = g(3) * t617;
t752 = qJDD(3) * pkin(3);
t690 = t763 * t446;
t417 = t616 * t437 + t690;
t751 = t417 * t548;
t749 = t473 * t776;
t748 = t649 * t776;
t745 = t529 * t776;
t744 = t530 * t552;
t743 = t530 * t776;
t735 = t604 * t607;
t730 = t616 * t446;
t718 = pkin(5) * t712 - qJ(6) * t713 - qJD(6) * t575 + t672;
t717 = -qJ(6) * t776 - t786;
t716 = pkin(5) * t776 - t785;
t419 = t445 * t763 - t730;
t711 = pkin(4) * t681 + qJD(6) - t419;
t710 = t613 ^ 2 + t614 ^ 2;
t708 = qJD(3) * t618;
t416 = t437 * t763 - t730;
t704 = qJD(6) - t416;
t602 = pkin(4) * t620 + pkin(3);
t684 = t510 * t706;
t683 = qJD(3) * t764;
t680 = pkin(4) * t617 + t754;
t678 = t710 * qJD(1) ^ 2;
t677 = -qJD(4) * t486 - t461;
t676 = t616 * t409 + t763 * t412 + t437 * t681 - t446 * t705;
t488 = qJD(2) * t571 - t583 * t708 + t584 * t683;
t673 = 0.2e1 * t710;
t418 = t616 * t445 + t690;
t671 = pkin(4) * t705 - t418;
t670 = -g(1) * t535 + g(2) * t537;
t536 = t605 * t724 - t721;
t538 = t605 * t731 + t732;
t669 = g(1) * t536 - g(2) * t538;
t666 = -t575 * t428 - t473 * t713;
t489 = pkin(4) * t738 - t778;
t660 = t552 * t742 - t783;
t659 = pkin(5) * t607 + qJ(6) * t606 + t602;
t463 = pkin(4) * t796 + t488;
t657 = t602 * t605 - t604 * t622 + t598;
t655 = t450 * t763 - t616 * t457;
t645 = -t545 * t764 - t618 * t546 + t572 * t708 - t573 * t683;
t642 = t616 * t424 + t763 * t430 + t450 * t681 - t457 * t705;
t641 = t648 * t506 + (-t756 - t795) * t607;
t462 = -t645 - t752;
t635 = g(1) * t538 + g(2) * t536 + g(3) * t735 - t676;
t632 = t673 * t703 - t668;
t630 = t416 * t548 + t635;
t628 = -g(1) * (-t537 * pkin(5) + qJ(6) * t538) - g(2) * (-t535 * pkin(5) + qJ(6) * t536) - g(3) * (-pkin(5) * t736 + qJ(6) * t735);
t432 = -pkin(4) * t658 + t462;
t625 = t626 * t620;
t601 = -pkin(4) * t763 - pkin(5);
t597 = pkin(4) * t616 + qJ(6);
t578 = -qJDD(1) * t598 + qJDD(2);
t556 = t605 * t722 + t728;
t554 = -t605 * t723 + t727;
t516 = -pkin(5) * t647 - qJ(6) * t575 - t602;
t502 = t647 * t571;
t501 = t575 * t571;
t440 = t501 * pkin(5) - t502 * qJ(6) + t489;
t439 = t563 * t689 - t616 * t686 - t705 * t738 + (t563 * t616 + t571 * t780) * t620;
t438 = t521 * t571 - t563 * t647;
t434 = pkin(4) * t530 + t443;
t426 = pkin(5) * t652 - t655;
t425 = -qJ(6) * t652 + t787;
t415 = t548 * qJ(6) + t417;
t414 = -t548 * pkin(5) + t704;
t406 = pkin(5) * t439 + qJ(6) * t438 - qJD(6) * t502 + t463;
t405 = t428 * pkin(5) + t427 * qJ(6) - qJD(6) * t649 + t432;
t404 = -t564 * pkin(5) - t768;
t403 = qJ(6) * t564 - qJD(6) * t652 + t642;
t402 = t675 - t773;
t401 = t676 + t784;
t1 = [((-t526 * t706 + t499) * t552 + t504 * t509 - t663 * t652 + t454 * t564 - t488 * t529 + t778 * t658 + t571 * t684 - g(1) * t554 - g(2) * t556 + ((-qJD(4) * t514 - t487) * t552 - t526 * t509 - t677 * t652 + t462 * t571 + t510 * t563) * t617) * MDP(20) + (-g(1) * t553 - g(2) * t555 - t455 * t564 + t462 * t737 + t488 * t530 - t509 * t714 + t510 * t779 - t552 * t643 - t626 * t778 + t644 * t652) * MDP(21) + (-qJD(3) * t488 + qJDD(3) * t778 - t515 * t598 + t564 * t580 - t578 * t652 + t605 * t781) * MDP(13) + (t563 * t776 + t571 * t627) * MDP(8) + (t416 * t564 + t489 * t428 + t432 * t501 + t469 * t439 + t463 * t473 + t655 * t506 + t548 * t768 + t652 * t675 + t669) * MDP(27) + (t427 * t652 - t438 * t548 + t502 * t506 + t564 * t649) * MDP(24) + (-t401 * t652 + t403 * t548 - t405 * t502 - t406 * t649 + t415 * t564 + t425 * t506 + t427 * t440 + t431 * t438 - t670) * MDP(31) + (-t509 * t652 + t552 * t564) * MDP(19) + (-t506 * t652 + t548 * t564) * MDP(26) + (t428 * t652 - t439 * t548 - t473 * t564 - t501 * t506) * MDP(25) + (-qJD(3) * t564 + qJDD(3) * t652) * MDP(11) + (t402 * t652 - t404 * t548 + t405 * t501 + t406 * t473 - t414 * t564 - t426 * t506 + t428 * t440 + t431 * t439 + t669) * MDP(29) + (t495 * t571 + t530 * t564 + t552 * t779 - t626 * t652) * MDP(17) + (-t401 * t501 + t402 * t502 - t403 * t473 + t404 * t649 - t414 * t438 - t415 * t439 - t425 * t428 - t426 * t427 + t604 * t781) * MDP(30) + t781 * MDP(2) + (-t417 * t564 - t489 * t427 + t432 * t502 - t469 * t438 + t463 * t649 - t506 * t787 - t548 * t642 + t652 * t676 + t670) * MDP(28) + (-t571 * t515 + t563 * t802 - t564 * t776 + t627 * t652) * MDP(9) + (-t509 * t738 + t529 * t564 - t552 * t796 - t652 * t658) * MDP(18) + (t529 * t779 - t530 * t796 - t626 * t738 + t658 * t737) * MDP(16) + (t614 * MDP(4) - t613 * MDP(5)) * (t661 + t753) + (t427 * t501 - t428 * t502 + t438 * t473 - t439 * t649) * MDP(23) + (-t427 * t502 - t438 * t649) * MDP(22) + (t530 * t779 + t571 * t625) * MDP(15) + (pkin(1) * t661 + (t710 * t775 + t632) * qJ(2)) * MDP(7) + (t673 * t775 + t632) * MDP(6) + qJDD(1) * MDP(1) + (qJD(3) * t563 + qJDD(3) * t571) * MDP(10) + (-g(1) * t734 + g(2) * t733 - t487 * qJD(3) - t526 * qJDD(3) + t580 * t563 + t578 * t571 - t598 * t627) * MDP(14) + t668 * MDP(3) + (t401 * t425 + t415 * t403 + t405 * t440 + t431 * t406 + t402 * t426 + t414 * t404 - g(1) * (-pkin(5) * t536 - qJ(6) * t535) - g(2) * (pkin(5) * t538 + qJ(6) * t537) + (-g(1) * t680 - g(2) * t657) * t621 + (g(1) * t657 - g(2) * t680) * t619) * MDP(32); -MDP(4) * t700 + MDP(5) * t701 - MDP(6) * t678 + (-qJ(2) * t678 - t661) * MDP(7) + (t662 + 0.2e1 * t793) * MDP(13) + (0.2e1 * qJD(3) * t802 + t638) * MDP(14) + (t660 + t745) * MDP(20) + (-t743 + t771) * MDP(21) + (t666 - t769) * MDP(30) + (t401 * t575 - t402 * t647 + t414 * t712 + t415 * t713 - t431 * t776 - t781) * MDP(32) + (-MDP(28) + MDP(31)) * (t664 + t748) + (MDP(27) + MDP(29)) * (t665 - t749); (t602 * t427 + t432 * t575 + t713 * t469 + t548 * t786 + t672 * t649 - t794) * MDP(28) + (-t405 * t575 + t427 * t516 - t431 * t713 + t548 * t717 - t649 * t718 + t794) * MDP(31) + t638 * MDP(10) + (-pkin(3) * t696 + t518 * t529 - pkin(8) * t747 + (-t462 + t636 + t752) * t620 + (-t498 + (t510 + t782) * t617 - pkin(8) * t706) * t552) * MDP(20) + (t617 * t626 + t744 * t620) * MDP(15) + (t617 * t658 + t625 - t804 * t530 + (t706 - t741) * t529) * MDP(16) + t791 * MDP(11) + (t518 * qJD(3) + t636 + t645) * MDP(13) + (-t405 * t647 + t428 * t516 + t431 * t712 + t473 * t718 - t548 * t716 + t641) * MDP(29) + (-t427 * t575 + t649 * t713) * MDP(22) + (-pkin(3) * t626 - t510 * t741 - t518 * t530 + t552 * t715 + t605 * t755 + t684 + t783 * pkin(8) + (t462 + t795) * t617) * MDP(21) + (-t602 * t428 - t432 * t647 + t712 * t469 + t672 * t473 + t548 * t785 + t641) * MDP(27) + (t401 * t528 - t402 * t648 + t405 * t516 + t718 * t431 + t717 * t415 + t716 * t414 + (-g(3) * t659 + t622 * t668) * t605 + (g(3) * t622 + t659 * t668) * t604) * MDP(32) + (t666 + t769) * MDP(23) + (t401 * t647 + t402 * t575 + t414 * t713 - t415 * t712 + t427 * t648 - t428 * t528 - t473 * t717 + t649 * t716 - t770) * MDP(30) + (-t580 * t802 - t654 + t770) * MDP(14) + (-t743 - t771) * MDP(17) + qJDD(3) * MDP(12) - t802 ^ 2 * MDP(9) + (t664 - t748) * MDP(24) + (t665 + t749) * MDP(25) + (t660 - t745) * MDP(18) + (MDP(11) * qJD(3) - MDP(13) * t580 - MDP(19) * t552 - MDP(20) * t454 + MDP(21) * t455 - MDP(26) * t548 - MDP(27) * t416 + MDP(28) * t417 + MDP(29) * t414 - MDP(31) * t415 - MDP(8) * t802 + MDP(9) * t776) * t776; -t530 * t529 * MDP(15) + (-t529 ^ 2 + t530 ^ 2) * MDP(16) + (t617 * qJDD(3) - t529 * t552 + t620 * t797 - t776 * t707) * MDP(17) + (t658 + t744) * MDP(18) + t509 * MDP(19) + (t455 * t552 - t510 * t530 + (t677 + t757) * t617 + t663 + t777) * MDP(20) + (g(1) * t556 - g(2) * t554 + t454 * t552 - t510 * t529 + t620 * t757 - t644) * MDP(21) + (t418 * t548 + (-t473 * t530 + t506 * t763 - t548 * t705) * pkin(4) + t772) * MDP(27) + (t419 * t548 + t800 + (-t506 * t616 - t530 * t649 - t548 * t681) * pkin(4) + t635) * MDP(28) + (-t434 * t473 - t506 * t601 - t548 * t671 - t629) * MDP(29) + (-t427 * t601 - t428 * t597 + (t415 + t671) * t649 + (t414 - t711) * t473) * MDP(30) + (t434 * t649 + t506 * t597 + t548 * t711 - t635 + t784 - t801) * MDP(31) + (t401 * t597 + t402 * t601 - t431 * t434 - t414 * t418 + t711 * t415 + (t414 * t705 + t604 * t755 + t777) * pkin(4) + t628) * MDP(32) + t803; (t751 + t772) * MDP(27) + (t630 + t800) * MDP(28) + (-t443 * t473 + t500 - t629 + t751) * MDP(29) + (pkin(5) * t427 - qJ(6) * t428 + (t415 - t417) * t649 + (t414 - t704) * t473) * MDP(30) + (t443 * t649 + 0.2e1 * t494 + 0.2e1 * t539 - t630 - t801) * MDP(31) + (-t402 * pkin(5) + t401 * qJ(6) - t414 * t417 + t415 * t704 - t431 * t443 + t628) * MDP(32) + t803; (-qJDD(4) - qJDD(5) + t791 + t807) * MDP(29) + t413 * MDP(30) + (-t548 ^ 2 - t765) * MDP(31) + (-t415 * t548 + t629) * MDP(32);];
tau  = t1;

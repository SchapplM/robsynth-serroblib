% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:04:26
% EndTime: 2019-03-09 18:04:42
% DurationCPUTime: 11.26s
% Computational Cost: add. (13379->531), mult. (32929->678), div. (0->0), fcn. (26066->16), ass. (0->258)
t671 = cos(qJ(5));
t672 = cos(qJ(3));
t673 = cos(qJ(2));
t753 = qJD(1) * t673;
t734 = t672 * t753;
t667 = sin(qJ(3));
t668 = sin(qJ(2));
t754 = qJD(1) * t668;
t735 = t667 * t754;
t586 = -t734 + t735;
t588 = -t667 * t753 - t672 * t754;
t663 = sin(pkin(11));
t664 = cos(pkin(11));
t703 = t586 * t664 - t588 * t663;
t551 = t671 * t703;
t555 = t586 * t663 + t588 * t664;
t666 = sin(qJ(5));
t512 = -t555 * t666 + t551;
t670 = cos(qJ(6));
t748 = qJD(6) * t670;
t826 = t512 * t670 + t748;
t659 = qJD(2) + qJD(3);
t742 = qJDD(1) * t673;
t744 = qJD(1) * qJD(2);
t732 = t673 * t744;
t743 = qJDD(1) * t668;
t806 = t732 + t743;
t539 = qJD(3) * t734 - t659 * t735 + t667 * t742 + t672 * t806;
t600 = t667 * t673 + t668 * t672;
t563 = t659 * t600;
t710 = t667 * t743 - t672 * t742;
t540 = qJD(1) * t563 + t710;
t490 = -t539 * t663 - t540 * t664;
t491 = t539 * t664 - t540 * t663;
t750 = qJD(5) * t666;
t448 = -qJD(5) * t551 + t666 * t490 + t671 * t491 + t555 * t750;
t657 = qJDD(2) + qJDD(3);
t652 = qJDD(5) + t657;
t653 = qJD(5) + t659;
t665 = sin(qJ(6));
t738 = t670 * t448 + t665 * t652 + t653 * t748;
t749 = qJD(6) * t665;
t804 = -t671 * t555 - t666 * t703;
t431 = -t749 * t804 + t738;
t430 = t431 * t670;
t504 = t653 * t665 + t670 * t804;
t629 = t670 * t652;
t432 = qJD(6) * t504 + t448 * t665 - t629;
t502 = -t670 * t653 + t665 * t804;
t825 = -t665 * t432 - t502 * t826 + t430;
t429 = t431 * t665;
t449 = qJD(5) * t804 - t671 * t490 + t666 * t491;
t447 = qJDD(6) + t449;
t444 = t665 * t447;
t777 = t512 * t653;
t779 = t804 * t653;
t781 = t504 * t804;
t805 = qJD(6) + t512;
t824 = t652 * MDP(24) + (-t449 + t779) * MDP(23) - t512 ^ 2 * MDP(21) + (MDP(20) * t512 + MDP(21) * t804 - MDP(31) * t805) * t804 + (t448 + t777) * MDP(22) + (t504 * t826 + t429) * MDP(27) + (t805 * t826 + t444 - t781) * MDP(29);
t583 = t588 * qJ(4);
t797 = pkin(7) + pkin(8);
t620 = t797 * t673;
t607 = qJD(1) * t620;
t589 = t667 * t607;
t619 = t797 * t668;
t605 = qJD(1) * t619;
t787 = qJD(2) * pkin(2);
t595 = -t605 + t787;
t724 = t672 * t595 - t589;
t537 = t583 + t724;
t528 = pkin(3) * t659 + t537;
t593 = t672 * t607;
t702 = -t595 * t667 - t593;
t786 = qJ(4) * t586;
t538 = -t702 - t786;
t529 = t663 * t538;
t483 = t664 * t528 - t529;
t552 = t555 * pkin(9);
t474 = pkin(4) * t659 + t483 + t552;
t774 = t664 * t538;
t484 = t663 * t528 + t774;
t794 = pkin(9) * t703;
t476 = t484 - t794;
t442 = t474 * t671 - t476 * t666;
t440 = -pkin(5) * t653 - t442;
t785 = t440 * t512;
t656 = t673 * pkin(2);
t789 = pkin(1) + t656;
t618 = t789 * qJD(1);
t565 = t586 * pkin(3) + qJD(4) - t618;
t523 = pkin(4) * t703 + t565;
t662 = qJ(2) + qJ(3);
t645 = pkin(11) + qJ(5) + t662;
t635 = sin(t645);
t621 = g(3) * t635;
t636 = cos(t645);
t669 = sin(qJ(1));
t674 = cos(qJ(1));
t712 = g(1) * t674 + g(2) * t669;
t564 = qJDD(2) * pkin(2) - t797 * t806;
t733 = t668 * t744;
t566 = t797 * (-t733 + t742);
t688 = qJD(3) * t702 + t672 * t564 - t667 * t566;
t469 = pkin(3) * t657 - qJ(4) * t539 + qJD(4) * t588 + t688;
t752 = qJD(3) * t667;
t800 = (qJD(3) * t595 + t566) * t672 + t667 * t564 - t607 * t752;
t473 = -qJ(4) * t540 - qJD(4) * t586 + t800;
t435 = t664 * t469 - t473 * t663;
t426 = pkin(4) * t657 - pkin(9) * t491 + t435;
t436 = t663 * t469 + t664 * t473;
t427 = pkin(9) * t490 + t436;
t801 = (qJD(5) * t474 + t427) * t671 + t666 * t426 - t476 * t750;
t683 = t512 * t523 + t636 * t712 + t621 - t801;
t821 = t712 * t635;
t817 = pkin(5) * t804;
t782 = t502 * t804;
t723 = t605 * t667 - t593;
t542 = t723 + t786;
t760 = -t672 * t605 - t589;
t543 = t583 + t760;
t773 = t664 * t667;
t788 = pkin(2) * qJD(3);
t761 = -t664 * t542 + t543 * t663 + (-t663 * t672 - t773) * t788;
t775 = t663 * t667;
t802 = -t663 * t542 - t664 * t543 + (t664 * t672 - t775) * t788;
t654 = sin(t662);
t655 = cos(t662);
t816 = -g(3) * t655 + t654 * t712;
t443 = t474 * t666 + t476 * t671;
t441 = pkin(10) * t653 + t443;
t458 = t512 * pkin(5) - pkin(10) * t804 + t523;
t422 = -t441 * t665 + t458 * t670;
t815 = -t422 * t804 + t440 * t749 + t670 * t821;
t423 = t441 * t670 + t458 * t665;
t799 = qJD(5) * t443 - t671 * t426 + t666 * t427;
t414 = -pkin(5) * t652 + t799;
t791 = g(3) * t636;
t810 = t414 + t791;
t814 = t423 * t804 + t440 * t748 + t665 * t810;
t679 = -t523 * t804 - t791 - t799 + t821;
t809 = t794 - t761;
t808 = t552 - t802;
t803 = t805 * (pkin(10) * t805 + t817);
t648 = pkin(2) * t672 + pkin(3);
t581 = -pkin(2) * t775 + t664 * t648;
t575 = pkin(4) + t581;
t582 = pkin(2) * t773 + t648 * t663;
t764 = t666 * t575 + t671 * t582;
t759 = -t667 * t619 + t672 * t620;
t642 = pkin(3) * t664 + pkin(4);
t795 = pkin(3) * t663;
t758 = t666 * t642 + t671 * t795;
t599 = t667 * t668 - t672 * t673;
t736 = qJD(2) * t797;
t606 = t668 * t736;
t608 = t673 * t736;
t751 = qJD(3) * t672;
t694 = -t672 * t606 - t667 * t608 - t619 * t751 - t620 * t752;
t496 = -qJ(4) * t563 - qJD(4) * t599 + t694;
t562 = t659 * t599;
t687 = -qJD(3) * t759 + t606 * t667 - t672 * t608;
t497 = qJ(4) * t562 - qJD(4) * t600 + t687;
t463 = -t496 * t663 + t664 * t497;
t522 = -t562 * t664 - t563 * t663;
t454 = -pkin(9) * t522 + t463;
t464 = t664 * t496 + t663 * t497;
t521 = t562 * t663 - t563 * t664;
t455 = pkin(9) * t521 + t464;
t722 = -t672 * t619 - t620 * t667;
t553 = -qJ(4) * t600 + t722;
t554 = -qJ(4) * t599 + t759;
t506 = t664 * t553 - t554 * t663;
t560 = -t599 * t663 + t600 * t664;
t481 = -pkin(9) * t560 + t506;
t507 = t663 * t553 + t664 * t554;
t559 = -t599 * t664 - t600 * t663;
t482 = pkin(9) * t559 + t507;
t707 = t481 * t671 - t482 * t666;
t415 = qJD(5) * t707 + t454 * t666 + t455 * t671;
t457 = t481 * t666 + t482 * t671;
t705 = t671 * t559 - t560 * t666;
t461 = qJD(5) * t705 + t521 * t666 + t522 * t671;
t517 = t559 * t666 + t560 * t671;
t713 = pkin(3) * t599 - t789;
t533 = -pkin(4) * t559 + t713;
t465 = -pkin(5) * t705 - pkin(10) * t517 + t533;
t413 = pkin(10) * t652 + t801;
t718 = qJD(6) * t458 + t413;
t798 = t414 * t517 + t440 * t461 - t457 * t447 - (qJD(6) * t465 + t415) * t805 + t705 * t718;
t796 = pkin(3) * t588;
t784 = t440 * t517;
t783 = t465 * t447;
t780 = t504 * t665;
t772 = t665 * t669;
t771 = t665 * t674;
t770 = t669 * t670;
t445 = t670 * t447;
t769 = t670 * t674;
t704 = t575 * t671 - t582 * t666;
t766 = -qJD(5) * t704 + t666 * t809 + t671 * t808;
t765 = t764 * qJD(5) - t666 * t808 + t671 * t809;
t489 = t664 * t537 - t529;
t488 = -t537 * t663 - t774;
t477 = t488 + t794;
t478 = t552 + t489;
t700 = t642 * t671 - t666 * t795;
t763 = -t700 * qJD(5) + t477 * t666 + t478 * t671;
t762 = qJD(5) * t758 + t477 * t671 - t478 * t666;
t757 = pkin(3) * t655 + t656;
t660 = t668 ^ 2;
t756 = -t673 ^ 2 + t660;
t651 = t668 * t787;
t731 = pkin(3) * t563 + t651;
t721 = t805 * t665;
t584 = pkin(2) * t733 - qJDD(1) * t789;
t689 = pkin(3) * t540 + qJDD(4) + t584;
t466 = -pkin(4) * t490 + t689;
t418 = pkin(5) * t449 - pkin(10) * t448 + t466;
t716 = qJD(6) * t441 - t418;
t527 = -pkin(4) * t555 - t796;
t460 = pkin(10) * t512 + t527 + t817;
t545 = pkin(10) + t764;
t650 = pkin(2) * t754;
t715 = qJD(6) * t545 + t460 + t650;
t577 = pkin(10) + t758;
t714 = qJD(6) * t577 + t460;
t711 = g(1) * t669 - g(2) * t674;
t709 = -t447 * t545 + t785;
t708 = -t447 * t577 + t785;
t706 = -t483 * t703 - t484 * t555;
t495 = -pkin(4) * t521 + t731;
t701 = t445 - (t512 * t665 + t749) * t805;
t698 = -0.2e1 * pkin(1) * t744 - pkin(7) * qJDD(2);
t697 = t461 * t670 - t517 * t749;
t693 = -pkin(10) * t447 + t442 * t805 + t785;
t675 = qJD(2) ^ 2;
t691 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t675 + t711;
t676 = qJD(1) ^ 2;
t690 = pkin(1) * t676 - pkin(7) * qJDD(1) + t712;
t686 = -t670 * t810 + t815;
t684 = -t665 * t821 + t814;
t682 = g(3) * t654 - t618 * t586 + t655 * t712 - t800;
t678 = -t618 * t588 + t688 + t816;
t677 = -t588 * t586 * MDP(11) + (-t780 * t805 + t825) * MDP(28) + (t701 + t782) * MDP(30) + (t586 * t659 + t539) * MDP(13) + (-t710 + (-qJD(1) * t600 - t588) * t659) * MDP(14) + (-t586 ^ 2 + t588 ^ 2) * MDP(12) + t657 * MDP(15) + t824;
t658 = -qJ(4) - t797;
t604 = pkin(1) + t757;
t576 = -pkin(5) - t700;
t574 = t636 * t769 + t772;
t573 = -t636 * t771 + t770;
t572 = -t636 * t770 + t771;
t571 = t636 * t772 + t769;
t544 = -pkin(5) - t704;
t524 = t527 + t650;
t462 = qJD(5) * t517 - t671 * t521 + t522 * t666;
t424 = pkin(5) * t462 - pkin(10) * t461 + t495;
t417 = t670 * t418;
t416 = qJD(5) * t457 - t454 * t671 + t455 * t666;
t1 = [(-t416 * t653 + t449 * t533 + t462 * t523 - t466 * t705 + t495 * t512 + t636 * t711 + t652 * t707) * MDP(25) + (-t462 * t653 + t652 * t705) * MDP(23) + (-t435 * t560 + t436 * t559 + t463 * t555 - t464 * t703 - t483 * t522 + t484 * t521 + t507 * t490 - t506 * t491 - t712) * MDP(18) + (-t562 * t659 + t600 * t657) * MDP(13) + (t539 * t600 + t562 * t588) * MDP(11) + (-t539 * t599 - t540 * t600 + t562 * t586 + t563 * t588) * MDP(12) + (t448 * t705 - t449 * t517 - t461 * t512 - t462 * t804) * MDP(21) + (-t415 * t653 + t448 * t533 - t457 * t652 + t461 * t523 + t466 * t517 + t495 * t804 - t635 * t711) * MDP(26) + (t448 * t517 + t461 * t804) * MDP(20) + (-t540 * t789 - t618 * t563 + t584 * t599 + t586 * t651 + t655 * t711 + t657 * t722 + t659 * t687) * MDP(16) + (-t539 * t789 + t618 * t562 + t584 * t600 - t588 * t651 - t654 * t711 - t657 * t759 - t659 * t694) * MDP(17) + (-g(1) * t571 - g(2) * t573 + t416 * t504 - t423 * t462 - t707 * t431 + (-(-qJD(6) * t457 + t424) * t805 - t783 - t716 * t705 - qJD(6) * t784) * t665 + t798 * t670) * MDP(33) + (-g(1) * t572 - g(2) * t574 + t416 * t502 - t417 * t705 + t422 * t462 - t707 * t432 + (t424 * t805 + t783 + (t441 * t705 - t457 * t805 + t784) * qJD(6)) * t670 + t798 * t665) * MDP(32) + (-t431 * t705 + t445 * t517 + t462 * t504 + t697 * t805) * MDP(29) + (-t517 * t444 + t432 * t705 - t462 * t502 + (-t461 * t665 - t517 * t748) * t805) * MDP(30) + (-t447 * t705 + t462 * t805) * MDP(31) + (t430 * t517 + t504 * t697) * MDP(27) + ((-t502 * t670 - t780) * t461 + (-t429 - t432 * t670 + (t502 * t665 - t504 * t670) * qJD(6)) * t517) * MDP(28) + 0.2e1 * (t668 * t742 - t744 * t756) * MDP(5) + (qJDD(1) * t660 + 0.2e1 * t668 * t732) * MDP(4) + t711 * MDP(2) + t712 * MDP(3) + (t461 * t653 + t517 * t652) * MDP(22) + (-t563 * t659 - t599 * t657) * MDP(14) + (qJDD(2) * t668 + t673 * t675) * MDP(6) + (qJDD(2) * t673 - t668 * t675) * MDP(7) + (t436 * t507 + t484 * t464 + t435 * t506 + t483 * t463 + t689 * t713 + t565 * t731 - g(1) * (-t604 * t669 - t658 * t674) - g(2) * (t604 * t674 - t658 * t669)) * MDP(19) + (t668 * t698 + t673 * t691) * MDP(9) + (-t668 * t691 + t673 * t698) * MDP(10) + qJDD(1) * MDP(1); (t544 * t431 + t709 * t670 + t765 * t504 + (t665 * t715 + t670 * t766) * t805 + t684) * MDP(33) + t677 + (-g(3) * t673 + t668 * t690) * MDP(9) + (g(3) * t668 + t673 * t690) * MDP(10) + (-t524 * t512 + t652 * t704 - t653 * t765 + t679) * MDP(25) + (t436 * t582 + t435 * t581 - t565 * (t650 - t796) - g(3) * t757 - t712 * (-pkin(2) * t668 - pkin(3) * t654) + t802 * t484 + t761 * t483) * MDP(19) + (t760 * t659 + (t588 * t754 - t667 * t657 - t659 * t751) * pkin(2) + t682) * MDP(17) + (-t723 * t659 + (-t586 * t754 + t657 * t672 - t659 * t752) * pkin(2) + t678) * MDP(16) + (-t524 * t804 - t652 * t764 + t653 * t766 + t683) * MDP(26) + (t582 * t490 - t581 * t491 + t555 * t761 - t703 * t802 + t706) * MDP(18) + MDP(7) * t742 + MDP(6) * t743 + qJDD(2) * MDP(8) + (t544 * t432 + t709 * t665 + t765 * t502 + (t665 * t766 - t670 * t715) * t805 + t686) * MDP(32) + (-MDP(4) * t668 * t673 + MDP(5) * t756) * t676; (t576 * t431 + t708 * t670 + t762 * t504 + (t665 * t714 + t670 * t763) * t805 + t684) * MDP(33) + t677 + (-t527 * t804 - t652 * t758 + t653 * t763 + t683) * MDP(26) + (-t483 * t488 - t484 * t489 + (t435 * t664 + t436 * t663 + t565 * t588 + t816) * pkin(3)) * MDP(19) + (t576 * t432 + t708 * t665 + t762 * t502 + (t665 * t763 - t670 * t714) * t805 + t686) * MDP(32) + (t489 * t703 - t488 * t555 + (t490 * t663 - t491 * t664) * pkin(3) + t706) * MDP(18) + (-t527 * t512 + t652 * t700 - t653 * t762 + t679) * MDP(25) + (-t659 * t702 + t678) * MDP(16) + (t659 * t724 + t682) * MDP(17); (-t555 ^ 2 - t703 ^ 2) * MDP(18) + (-t483 * t555 + t484 * t703 + t689 - t711) * MDP(19) + (t449 + t779) * MDP(25) + (t448 - t777) * MDP(26) + (t701 - t782) * MDP(32) + (-t670 * t805 ^ 2 - t444 - t781) * MDP(33); (t443 * t653 + t679) * MDP(25) + (t442 * t653 + t683) * MDP(26) + (-t504 * t721 + t825) * MDP(28) + (-t721 * t805 + t445 + t782) * MDP(30) + (-pkin(5) * t432 - t443 * t502 + t693 * t665 + (-t810 - t803) * t670 + t815) * MDP(32) + (-pkin(5) * t431 - t443 * t504 + t693 * t670 + (-t821 + t803) * t665 + t814) * MDP(33) + t824; t504 * t502 * MDP(27) + (-t502 ^ 2 + t504 ^ 2) * MDP(28) + (t502 * t805 + t738) * MDP(29) + (t504 * t805 + t629) * MDP(30) + t447 * MDP(31) + (-g(1) * t573 + g(2) * t571 + t423 * t805 - t440 * t504 + t417) * MDP(32) + (g(1) * t574 - g(2) * t572 + t422 * t805 + t440 * t502) * MDP(33) + ((-t413 + t621) * MDP(33) + (-MDP(30) * t804 - MDP(32) * t441 - MDP(33) * t458) * qJD(6)) * t670 + (-qJD(6) * t804 * MDP(29) + (-qJD(6) * t653 - t448) * MDP(30) + (-t718 + t621) * MDP(32) + t716 * MDP(33)) * t665;];
tau  = t1;

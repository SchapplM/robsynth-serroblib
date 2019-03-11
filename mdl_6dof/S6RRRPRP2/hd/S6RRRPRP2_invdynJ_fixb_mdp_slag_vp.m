% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:37:37
% EndTime: 2019-03-09 16:37:50
% DurationCPUTime: 9.88s
% Computational Cost: add. (12955->628), mult. (30710->771), div. (0->0), fcn. (22955->14), ass. (0->290)
t657 = cos(qJ(5));
t750 = qJD(5) * t657;
t654 = sin(qJ(3));
t655 = sin(qJ(2));
t658 = cos(qJ(3));
t659 = cos(qJ(2));
t588 = t654 * t655 - t658 * t659;
t571 = t588 * qJD(1);
t755 = qJD(1) * t655;
t774 = t654 * t659;
t573 = -qJD(1) * t774 - t658 * t755;
t652 = sin(pkin(10));
t804 = cos(pkin(10));
t722 = -t804 * t571 + t573 * t652;
t786 = t722 * t657;
t836 = t750 - t786;
t653 = sin(qJ(5));
t751 = qJD(5) * t653;
t787 = t722 * t653;
t835 = t751 - t787;
t646 = t659 * pkin(2);
t807 = pkin(1) + t646;
t712 = pkin(3) * t588 - t807;
t606 = qJD(1) * t807;
t834 = qJDD(1) * t807;
t566 = t573 * qJ(4);
t821 = pkin(7) + pkin(8);
t608 = t821 * t659;
t596 = qJD(1) * t608;
t574 = t654 * t596;
t607 = t821 * t655;
t594 = qJD(1) * t607;
t805 = qJD(2) * pkin(2);
t582 = -t594 + t805;
t721 = t658 * t582 - t574;
t522 = t566 + t721;
t743 = qJD(2) + qJD(3);
t515 = pkin(3) * t743 + t522;
t578 = t658 * t596;
t698 = -t582 * t654 - t578;
t803 = qJ(4) * t571;
t523 = -t698 - t803;
t726 = t804 * t523;
t480 = t652 * t515 + t726;
t478 = pkin(9) * t743 + t480;
t550 = pkin(3) * t571 + qJD(4) - t606;
t690 = -t652 * t571 - t573 * t804;
t492 = -pkin(4) * t722 - pkin(9) * t690 + t550;
t453 = t478 * t657 + t492 * t653;
t824 = qJD(5) - t722;
t445 = qJ(6) * t824 + t453;
t833 = t445 * t824;
t527 = t653 * t743 + t657 * t690;
t719 = t527 * t824;
t589 = t655 * t658 + t774;
t548 = t743 * t589;
t529 = t594 * t654 - t578 + t803;
t760 = -t658 * t594 - t574;
t530 = t566 + t760;
t498 = t652 * t529 + t530 * t804;
t778 = t652 * t654;
t806 = pkin(2) * qJD(3);
t563 = (t658 * t804 - t778) * t806;
t762 = t563 - t498;
t832 = t762 * t653;
t716 = qJD(1) * t743;
t704 = t655 * t716;
t744 = qJDD(1) * t659;
t747 = qJD(1) * qJD(2);
t732 = t659 * t747;
t745 = qJDD(1) * t655;
t689 = -t732 - t745;
t746 = qJD(1) * qJD(3);
t825 = -t659 * t746 + t689;
t713 = t654 * t744 - t658 * t825;
t677 = -t654 * t704 + t713;
t777 = t652 * t658;
t663 = t804 * t677 + (-t704 + t744) * t777 + t825 * t778;
t831 = qJD(5) * t743 + t663;
t651 = qJ(2) + qJ(3);
t643 = pkin(10) + t651;
t628 = sin(t643);
t660 = cos(qJ(1));
t783 = t628 * t660;
t656 = sin(qJ(1));
t813 = g(2) * t656;
t830 = g(1) * t783 + t628 * t813;
t706 = t657 * pkin(5) + t653 * qJ(6);
t709 = g(1) * t660 + t813;
t549 = qJDD(2) * pkin(2) + t689 * t821;
t546 = t658 * t549;
t733 = t655 * t747;
t688 = -t733 + t744;
t552 = t821 * t688;
t718 = qJD(3) * t582 + t552;
t753 = qJD(3) * t658;
t734 = t596 * t753;
t741 = qJDD(2) + qJDD(3);
t467 = -t734 + t546 - t713 * qJ(4) + t573 * qJD(4) + t741 * pkin(3) + (qJ(4) * t704 - t718) * t654;
t754 = qJD(3) * t654;
t569 = t596 * t754;
t472 = -t571 * qJD(4) - t569 + (qJ(4) * t825 + t549) * t654 + ((-t655 * t746 + t688) * qJ(4) + t718) * t658;
t435 = t652 * t467 + t804 * t472;
t433 = pkin(9) * t741 + t435;
t668 = t589 * t716;
t684 = t588 * qJDD(1);
t664 = -t684 - t668;
t490 = -t652 * t677 + t804 * t664;
t626 = pkin(2) * t733;
t742 = qJDD(4) + t626;
t801 = qJDD(1) * pkin(1);
t442 = -pkin(2) * t744 - pkin(3) * t664 - t490 * pkin(4) - pkin(9) * t663 + t742 - t801;
t686 = t657 * t433 + t653 * t442 - t478 * t751 + t492 * t750;
t487 = qJDD(5) - t490;
t802 = qJ(6) * t487;
t421 = qJD(6) * t824 + t686 + t802;
t714 = t653 * t433 - t657 * t442 + t478 * t750 + t492 * t751;
t817 = pkin(5) * t487;
t423 = qJDD(6) + t714 - t817;
t829 = t421 * t657 + t423 * t653;
t544 = t588 * t804 + t589 * t652;
t545 = -t652 * t588 + t589 * t804;
t506 = pkin(4) * t544 - pkin(9) * t545 + t712;
t720 = -t658 * t607 - t608 * t654;
t536 = -qJ(4) * t589 + t720;
t758 = -t654 * t607 + t658 * t608;
t537 = -qJ(4) * t588 + t758;
t508 = t652 * t536 + t537 * t804;
t765 = t653 * t506 + t657 * t508;
t725 = t804 * t654;
t763 = t804 * t529 - t530 * t652 + (t725 + t777) * t806;
t828 = t835 * pkin(5) - qJ(6) * t836 - t653 * qJD(6);
t629 = cos(t643);
t827 = -t629 * pkin(4) - t628 * pkin(9);
t826 = g(1) * t656 - g(2) * t660;
t823 = pkin(3) * t668 + qJDD(1) * t712 + t742;
t822 = t527 ^ 2;
t820 = pkin(3) * t573;
t818 = pkin(3) * t652;
t816 = pkin(5) * t690;
t812 = g(3) * t628;
t811 = g(3) * t629;
t645 = cos(t651);
t810 = g(3) * t645;
t809 = g(3) * t653;
t434 = t467 * t804 - t652 * t472;
t432 = -pkin(4) * t741 - t434;
t470 = -t653 * t741 - t657 * t831 + t690 * t751;
t471 = t653 * t831 - t657 * t741 + t690 * t750;
t424 = t471 * pkin(5) + t470 * qJ(6) - t527 * qJD(6) + t432;
t800 = t424 * t653;
t799 = t453 * t824;
t516 = t652 * t523;
t479 = t515 * t804 - t516;
t477 = -pkin(4) * t743 - t479;
t525 = t653 * t690 - t657 * t743;
t455 = t525 * pkin(5) - t527 * qJ(6) + t477;
t798 = t455 * t722;
t797 = t470 * t653;
t796 = t471 * t657;
t795 = t477 * t722;
t638 = pkin(2) * t658 + pkin(3);
t565 = pkin(2) * t725 + t652 * t638;
t560 = pkin(9) + t565;
t794 = t487 * t560;
t630 = pkin(9) + t818;
t793 = t487 * t630;
t792 = t525 * t722;
t791 = t525 * t653;
t790 = t527 * t525;
t789 = t527 * t653;
t788 = t527 * t657;
t785 = t545 * t657;
t784 = t563 * t657;
t782 = t629 * t660;
t644 = sin(t651);
t781 = t644 * t656;
t780 = t644 * t660;
t648 = -qJ(4) - t821;
t779 = t648 * t660;
t481 = t653 * t487;
t775 = t653 * t656;
t773 = t656 * t657;
t482 = t657 * t487;
t772 = t657 * t660;
t771 = t660 * t653;
t770 = -t653 * t471 - t525 * t750;
t769 = t750 * t824 + t481;
t768 = t787 * t824 + t482;
t489 = t522 * t804 - t516;
t502 = pkin(4) * t690 - pkin(9) * t722 - t820;
t767 = t657 * t489 + t653 * t502;
t640 = pkin(2) * t755;
t499 = t502 + t640;
t766 = t657 * t498 + t653 * t499;
t764 = t828 + t763;
t488 = t522 * t652 + t726;
t761 = -t488 + t828;
t759 = t830 * t657;
t633 = pkin(3) * t645;
t757 = t633 + t646;
t649 = t655 ^ 2;
t756 = -t659 ^ 2 + t649;
t752 = qJD(5) * t630;
t452 = -t478 * t653 + t492 * t657;
t748 = qJD(6) - t452;
t642 = t655 * t805;
t738 = t629 * t809 - t653 * t830;
t736 = qJD(2) * t821;
t735 = t804 * pkin(3);
t730 = pkin(3) * t548 + t642;
t728 = -t424 - t811;
t727 = -t432 - t811;
t534 = t690 * qJ(6);
t450 = t534 + t766;
t724 = -t450 + t784;
t723 = t499 * t657 + t816 + t832;
t595 = t655 * t736;
t597 = t659 * t736;
t687 = -t658 * t595 - t654 * t597 - t607 * t753 - t608 * t754;
t495 = -qJ(4) * t548 - qJD(4) * t588 + t687;
t547 = t743 * t588;
t676 = -qJD(3) * t758 + t595 * t654 - t658 * t597;
t496 = qJ(4) * t547 - qJD(4) * t589 + t676;
t459 = t495 * t652 - t804 * t496;
t507 = -t804 * t536 + t537 * t652;
t717 = t743 * t655;
t555 = t629 * t775 + t772;
t557 = t629 * t771 - t773;
t711 = -g(1) * t555 + g(2) * t557;
t556 = t629 * t773 - t771;
t558 = t629 * t772 + t775;
t710 = g(1) * t556 - g(2) * t558;
t705 = pkin(5) * t653 - qJ(6) * t657;
t564 = -pkin(2) * t778 + t638 * t804;
t703 = t629 * t706 + t633 - t827;
t444 = -pkin(5) * t824 + t748;
t702 = t444 * t657 - t445 * t653;
t701 = -t794 - t795;
t700 = -t793 - t795;
t699 = t479 * t722 + t480 * t690;
t697 = pkin(4) + t706;
t696 = t444 * t690 + t455 * t751 + t759;
t695 = -t452 * t690 + t477 * t751 + t759;
t559 = -pkin(4) - t564;
t694 = -0.2e1 * pkin(1) * t747 - pkin(7) * qJDD(2);
t513 = -t547 * t804 - t652 * t548;
t693 = t513 * t653 + t545 * t750;
t692 = -t513 * t657 + t545 * t751;
t691 = -t560 * t751 + t784;
t460 = t495 * t804 + t652 * t496;
t512 = -t547 * t652 + t548 * t804;
t465 = pkin(4) * t512 - pkin(9) * t513 + t730;
t685 = t657 * t460 + t653 * t465 + t506 * t750 - t508 * t751;
t683 = t432 * t653 + t453 * t690 + t477 * t750 + t738;
t682 = g(1) * t780 + g(2) * t781 - t606 * t573 + t546 - t810;
t681 = -t445 * t690 + t455 * t786 - t738 - t800;
t661 = qJD(2) ^ 2;
t680 = -pkin(7) * t661 + 0.2e1 * t801 + t826;
t662 = qJD(1) ^ 2;
t679 = pkin(1) * t662 - pkin(7) * qJDD(1) + t709;
t678 = g(1) * t557 + g(2) * t555 + t628 * t809 - t714;
t674 = qJD(5) * t702 + t829;
t673 = -t797 - t796 + (t788 + t791) * qJD(5);
t672 = t455 * t527 + qJDD(6) - t678;
t671 = -g(1) * t558 - g(2) * t556 - t657 * t812 + t686;
t670 = g(3) * t644 - t654 * t549 - t606 * t571 + t645 * t709 - t718 * t658 + t569;
t669 = -t573 * t571 * MDP(11) - t824 * t690 * MDP(24) + ((-t470 + t792) * t657 - t824 * t789 + t770) * MDP(21) + (t525 * t690 - t751 * t824 + t768) * MDP(23) + (-t527 * t690 - t786 * t824 + t769) * MDP(22) + (t657 * t719 - t797) * MDP(20) + (t571 * t743 + t677) * MDP(13) + (-t573 * t743 + t664) * MDP(14) + (-t571 ^ 2 + t573 ^ 2) * MDP(12) + t741 * MDP(15);
t666 = t709 * t697 * t628;
t665 = t444 * t836 - t835 * t445 - t629 * t709 - t812 + t829;
t631 = -t735 - pkin(4);
t602 = pkin(9) * t782;
t600 = t656 * t629 * pkin(9);
t598 = -pkin(2) * t655 - pkin(3) * t644;
t593 = pkin(1) + t757;
t581 = -t735 - t697;
t580 = t660 * t593;
t568 = t626 - t834;
t551 = t559 - t706;
t500 = pkin(5) * t527 + qJ(6) * t525;
t473 = t545 * t705 + t507;
t458 = -pkin(5) * t544 - t506 * t657 + t508 * t653;
t457 = qJ(6) * t544 + t765;
t449 = t489 * t653 - t502 * t657 - t816;
t448 = t534 + t767;
t447 = t525 * t824 - t470;
t428 = t705 * t513 + (qJD(5) * t706 - qJD(6) * t657) * t545 + t459;
t426 = -pkin(5) * t512 + qJD(5) * t765 + t460 * t653 - t465 * t657;
t425 = qJ(6) * t512 + qJD(6) * t544 + t685;
t1 = [(t421 * t457 + t445 * t425 + t424 * t473 + t455 * t428 + t423 * t458 + t444 * t426 - g(1) * (-pkin(5) * t556 - qJ(6) * t555 - t779) - g(2) * (pkin(4) * t782 + pkin(5) * t558 + pkin(9) * t783 + qJ(6) * t557 + t580) + (-g(1) * (-t593 + t827) + g(2) * t648) * t656) * MDP(30) + (-t807 * t713 + t568 * t589 + t606 * t547 - t687 * t743 - t758 * t741 + (t654 * t716 * t807 - t573 * t805) * t655 - t826 * t644) * MDP(17) + (-t434 * t545 - t435 * t544 + t459 * t690 + t460 * t722 - t479 * t513 - t480 * t512 + t508 * t490 + t507 * t663 - t709) * MDP(18) + (t435 * t508 + t480 * t460 - t434 * t507 - t479 * t459 + t823 * t712 + t550 * t730 - g(1) * (-t593 * t656 - t779) - g(2) * (-t648 * t656 + t580)) * MDP(19) + qJDD(1) * MDP(1) + (qJDD(2) * t655 + t659 * t661) * MDP(6) + (qJDD(2) * t659 - t655 * t661) * MDP(7) + (t571 * t642 + t645 * t826 + t676 * t743 + t720 * t741 + (t568 - t834) * t588 - 0.2e1 * t606 * t548) * MDP(16) + (t547 * t571 - t713 * t588 + t573 * t548 - t589 * t684 + (-t589 * t658 * t717 + (-t548 * t659 + t588 * t717) * t654) * qJD(1)) * MDP(12) + (-t425 * t525 + t426 * t527 - t457 * t471 - t458 * t470 + t826 * t628 + t702 * t513 + (-t421 * t653 + t423 * t657 + (-t444 * t653 - t445 * t657) * qJD(5)) * t545) * MDP(28) + t826 * MDP(2) + (t573 * t547 + t589 * t677) * MDP(11) + (-t547 * t743 + t589 * t741) * MDP(13) + (-t423 * t544 - t426 * t824 + t428 * t525 - t444 * t512 + t455 * t693 - t458 * t487 + t471 * t473 + t545 * t800 + t710) * MDP(27) + (t421 * t544 - t424 * t785 + t425 * t824 - t428 * t527 + t445 * t512 + t455 * t692 + t457 * t487 + t470 * t473 - t711) * MDP(29) + (t432 * t785 - t453 * t512 + t459 * t527 - t507 * t470 - t477 * t692 - t487 * t765 - t544 * t686 - t685 * t824 + t711) * MDP(26) + (-t471 * t544 - t481 * t545 - t512 * t525 - t693 * t824) * MDP(23) + (-t470 * t544 + t482 * t545 + t512 * t527 - t692 * t824) * MDP(22) + (t487 * t544 + t512 * t824) * MDP(24) + (-t714 * t544 + t452 * t512 + t459 * t525 + t507 * t471 + ((-qJD(5) * t508 + t465) * t824 + t506 * t487 + t477 * qJD(5) * t545) * t657 + ((-qJD(5) * t506 - t460) * t824 - t508 * t487 + t432 * t545 + t477 * t513) * t653 + t710) * MDP(25) + ((-t525 * t657 - t789) * t513 + (t797 - t796 + (-t788 + t791) * qJD(5)) * t545) * MDP(21) + (-t470 * t785 - t527 * t692) * MDP(20) + 0.2e1 * (t655 * t744 - t747 * t756) * MDP(5) + (-t548 * t743 - t588 * t741) * MDP(14) + (t655 * t694 + t659 * t680) * MDP(9) + (-t655 * t680 + t659 * t694) * MDP(10) + t709 * MDP(3) + (qJDD(1) * t649 + 0.2e1 * t655 * t732) * MDP(4); (g(3) * t655 + t659 * t679) * MDP(10) + (-g(3) * t659 + t655 * t679) * MDP(9) + (t559 * t471 + t727 * t657 + t701 * t653 + t763 * t525 + ((-qJD(5) * t560 - t499) * t657 - t832) * t824 + t695) * MDP(25) + (-t525 * t724 + t527 * t723 + t560 * t673 + t665) * MDP(28) + (t424 * t551 - g(1) * (t598 * t660 + t602) - g(2) * (t598 * t656 + t600) - g(3) * (t646 + t703) + t764 * t455 + t724 * t445 + t723 * t444 + t666 + t674 * t560) * MDP(30) + (t565 * t490 - t564 * t663 + t690 * t763 + t722 * t762 + t699) * MDP(18) + (-t559 * t470 + t701 * t657 + t763 * t527 + (-t691 + t766) * t824 + t683) * MDP(26) + (t471 * t551 + t728 * t657 + (-t794 - t798) * t653 + t764 * t525 + (-t560 * t750 - t723) * t824 + t696) * MDP(27) + qJDD(2) * MDP(8) + (t760 * t743 + (t573 * t755 - t654 * t741 - t743 * t753) * pkin(2) + t670) * MDP(17) + (-t734 + t578 * t743 + (-t594 * t743 - t718) * t654 + (-t571 * t755 + t658 * t741 - t743 * t754) * pkin(2) + t682) * MDP(16) + t669 + (t435 * t565 + t434 * t564 - t550 * (t640 - t820) - g(3) * t757 - t709 * t598 + t762 * t480 - t763 * t479) * MDP(19) + (t470 * t551 + (-qJD(5) * t455 + t794) * t657 - t764 * t527 + (-t450 + t691) * t824 + t681) * MDP(29) + MDP(7) * t744 + MDP(6) * t745 + (-MDP(4) * t655 * t659 + MDP(5) * t756) * t662; (t479 * t488 - t480 * t489 + (t434 * t804 + t435 * t652 + t550 * t573 + t644 * t709 - t810) * pkin(3)) * MDP(19) + (t630 * t482 - t448 * t824 + t470 * t581 - t761 * t527 + (-t630 * t653 * t824 - t455 * t657) * qJD(5) + t681) * MDP(29) + (t448 * t525 - t449 * t527 + t630 * t673 + t665) * MDP(28) + (t424 * t581 - t445 * t448 - t444 * t449 - g(1) * (-pkin(3) * t780 + t602) - g(2) * (-pkin(3) * t781 + t600) - g(3) * t703 + t761 * t455 + t674 * t630 + t666) * MDP(30) + (-t488 * t690 - t489 * t722 + t490 * t818 - t663 * t735 + t699) * MDP(18) + (t631 * t471 - t488 * t525 + (t489 * t824 + t700) * t653 + ((-t502 - t752) * t824 + t727) * t657 + t695) * MDP(25) + (t721 * t743 + t670) * MDP(17) + t669 + (t449 * t824 + t471 * t581 + (-t793 - t798) * t653 + t761 * t525 + (-t752 * t824 + t728) * t657 + t696) * MDP(27) + (-t631 * t470 - t488 * t527 + t700 * t657 + (t630 * t751 + t767) * t824 + t683) * MDP(26) + (-qJD(2) * t698 - t654 * t552 + t682) * MDP(16); -t722 ^ 2 * MDP(18) + t768 * MDP(25) + t770 * MDP(28) + t769 * MDP(29) - t826 * MDP(30) + (-t690 * MDP(18) - t455 * MDP(30) + (-MDP(26) + MDP(29)) * t527 + (-MDP(25) - MDP(27)) * t525) * t690 + (t487 * MDP(27) + (t470 + t792) * MDP(28) + (-t423 + t833) * MDP(30) + (-MDP(26) * t824 - MDP(29) * t722) * t824) * t657 + (-t487 * MDP(26) + (t444 * t824 + t421) * MDP(30) + MDP(28) * t719 + (-qJD(5) * MDP(25) - MDP(27) * t824) * t824) * t653 + (t479 * t690 - t480 * t722 + t823 - t826) * MDP(19); MDP(20) * t790 + (-t525 ^ 2 + t822) * MDP(21) + t447 * MDP(22) + (-t471 + t719) * MDP(23) + t487 * MDP(24) + (-t477 * t527 + t678 + t799) * MDP(25) + (t452 * t824 + t477 * t525 - t671) * MDP(26) + (-t500 * t525 - t672 + t799 + 0.2e1 * t817) * MDP(27) + (pkin(5) * t470 - qJ(6) * t471 + (t445 - t453) * t527 + (t444 - t748) * t525) * MDP(28) + (0.2e1 * t802 - t455 * t525 + t500 * t527 + (0.2e1 * qJD(6) - t452) * t824 + t671) * MDP(29) + (t421 * qJ(6) - t423 * pkin(5) - t455 * t500 - t444 * t453 - g(1) * (-pkin(5) * t557 + qJ(6) * t558) - g(2) * (-pkin(5) * t555 + qJ(6) * t556) + t705 * t812 + t748 * t445) * MDP(30); (-t487 + t790) * MDP(27) + t447 * MDP(28) + (-t824 ^ 2 - t822) * MDP(29) + (t672 - t817 - t833) * MDP(30);];
tau  = t1;

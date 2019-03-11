% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:23:27
% EndTime: 2019-03-10 01:23:45
% DurationCPUTime: 12.28s
% Computational Cost: add. (11003->569), mult. (26667->753), div. (0->0), fcn. (19277->8), ass. (0->251)
t658 = cos(qJ(3));
t656 = sin(qJ(2));
t659 = cos(qJ(2));
t771 = t658 * t659;
t690 = pkin(3) * t656 - pkin(9) * t771;
t791 = pkin(8) + pkin(9);
t725 = qJD(3) * t791;
t696 = pkin(2) * t656 - pkin(8) * t659;
t615 = t696 * qJD(1);
t655 = sin(qJ(3));
t742 = qJD(1) * t656;
t723 = t655 * t742;
t748 = pkin(7) * t723 + t658 * t615;
t831 = qJD(1) * t690 + t658 * t725 + t748;
t595 = t655 * t615;
t773 = t656 * t658;
t774 = t655 * t659;
t830 = -t595 - (-pkin(7) * t773 - pkin(9) * t774) * qJD(1) - t655 * t725;
t729 = qJD(1) * qJD(2);
t716 = t659 * t729;
t736 = qJD(3) * t655;
t721 = t656 * t736;
t728 = qJD(2) * qJD(3);
t570 = -qJD(1) * t721 + (t716 + t728) * t658;
t735 = qJD(3) * t658;
t719 = t656 * t735;
t737 = qJD(2) * t659;
t722 = t655 * t737;
t678 = t719 + t722;
t571 = qJD(1) * t678 + t655 * t728;
t730 = t658 * qJD(2);
t608 = t723 - t730;
t739 = qJD(2) * t655;
t610 = t658 * t742 + t739;
t654 = sin(qJ(4));
t657 = cos(qJ(4));
t732 = qJD(4) * t657;
t734 = qJD(4) * t654;
t487 = t657 * t570 - t654 * t571 - t608 * t732 - t610 * t734;
t653 = sin(qJ(5));
t689 = t654 * t570 + t657 * t571 - t608 * t734 + t610 * t732;
t694 = -t608 * t654 + t657 * t610;
t695 = -t608 * t657 - t654 * t610;
t790 = cos(qJ(5));
t718 = qJD(5) * t790;
t731 = qJD(5) * t653;
t440 = -t790 * t487 + t653 * t689 + t694 * t731 - t695 * t718;
t795 = t653 * t695 + t694 * t790;
t441 = qJD(5) * t795 + t653 * t487 + t790 * t689;
t498 = t653 * t694 - t695 * t790;
t495 = t498 ^ 2;
t741 = qJD(1) * t659;
t638 = -qJD(3) + t741;
t629 = -qJD(4) + t638;
t621 = -qJD(5) + t629;
t717 = t656 * t729;
t792 = t795 ^ 2;
t829 = MDP(29) * t717 + (-t621 * t795 - t441) * MDP(28) + t498 * MDP(25) * t795 + (-t498 * t621 - t440) * MDP(27) + (-t495 + t792) * MDP(26);
t828 = MDP(22) * t717 + (t694 ^ 2 - t695 ^ 2) * MDP(19) + (t629 * t695 + t487) * MDP(20) + (-t629 * t694 - t689) * MDP(21) - t695 * MDP(18) * t694 + t829;
t692 = t654 * t655 - t657 * t658;
t800 = qJD(4) + qJD(3);
t817 = t800 * t692;
t753 = t692 * t741 - t817;
t625 = t791 * t655;
t626 = t791 * t658;
t749 = -t654 * t625 + t657 * t626;
t827 = qJD(4) * t749 + t830 * t654 + t831 * t657;
t611 = t654 * t658 + t655 * t657;
t577 = t611 * t741;
t681 = t611 * qJD(3);
t666 = -qJD(4) * t611 - t681;
t825 = t577 + t666;
t824 = -t625 * t732 - t831 * t654 + t830 * t657;
t818 = t498 * qJ(6);
t823 = pkin(4) * t742 + pkin(10) * t753 + t827;
t777 = t654 * t626;
t787 = pkin(10) * t611;
t822 = (-t777 - t787) * qJD(4) + t824 + (t577 - t681) * pkin(10);
t646 = pkin(7) * t741;
t789 = pkin(3) * t655;
t604 = t741 * t789 + t646;
t804 = pkin(3) * t736 - t604;
t620 = -pkin(2) * t659 - pkin(8) * t656 - pkin(1);
t601 = t620 * qJD(1);
t624 = qJD(2) * pkin(8) + t646;
t561 = t658 * t601 - t624 * t655;
t524 = -pkin(9) * t610 + t561;
t516 = -pkin(3) * t638 + t524;
t776 = t655 * t601;
t562 = t624 * t658 + t776;
t525 = -pkin(9) * t608 + t562;
t519 = t654 * t525;
t467 = t657 * t516 - t519;
t812 = pkin(10) * t694;
t454 = t467 - t812;
t449 = -pkin(4) * t629 + t454;
t521 = t657 * t525;
t468 = t654 * t516 + t521;
t811 = pkin(10) * t695;
t455 = t468 + t811;
t451 = t653 * t455;
t432 = t790 * t449 - t451;
t808 = qJ(6) * t795;
t426 = t432 - t808;
t425 = -pkin(5) * t621 + t426;
t453 = t790 * t455;
t433 = t653 * t449 + t453;
t427 = t433 - t818;
t821 = -t425 * t498 + t427 * t795;
t623 = -qJD(2) * pkin(2) + pkin(7) * t742;
t573 = pkin(3) * t608 + t623;
t514 = -pkin(4) * t695 + t573;
t618 = t696 * qJD(2);
t602 = qJD(1) * t618;
t700 = pkin(7) * t717;
t752 = -t658 * t602 - t655 * t700;
t674 = -qJD(3) * t562 - t752;
t477 = pkin(3) * t717 - pkin(9) * t570 + t674;
t688 = t601 * t735 + t655 * t602 - t624 * t736;
t671 = -t658 * t700 + t688;
t492 = -pkin(9) * t571 + t671;
t672 = -qJD(4) * t468 + t657 * t477 - t654 * t492;
t421 = pkin(4) * t717 - pkin(10) * t487 + t672;
t698 = -t654 * t477 - t657 * t492 - t516 * t732 + t525 * t734;
t424 = -pkin(10) * t689 - t698;
t699 = -t653 * t421 - t790 * t424 - t449 * t718 + t455 * t731;
t815 = t498 * t514 + t699;
t677 = t790 * t692;
t762 = qJD(5) * t677 + t611 * t731 - t653 * t825 - t753 * t790;
t558 = t611 * t790 - t653 * t692;
t761 = qJD(5) * t558 + t653 * t753 - t790 * t825;
t816 = -pkin(4) * t825 + t804;
t669 = -qJD(5) * t433 + t790 * t421 - t653 * t424;
t799 = -t514 * t795 + t669;
t813 = -0.2e1 * t729;
t810 = MDP(4) * t656;
t651 = t656 ^ 2;
t809 = MDP(5) * (-t659 ^ 2 + t651);
t459 = pkin(5) * t498 + qJD(6) + t514;
t807 = t459 * t795;
t584 = t692 * t656;
t607 = t658 * t620;
t788 = pkin(7) * t655;
t560 = -pkin(9) * t773 + t607 + (-pkin(3) - t788) * t659;
t640 = pkin(7) * t771;
t746 = t655 * t620 + t640;
t775 = t655 * t656;
t566 = -pkin(9) * t775 + t746;
t704 = t657 * t560 - t566 * t654;
t490 = -pkin(4) * t659 + pkin(10) * t584 + t704;
t685 = t611 * t656;
t754 = t654 * t560 + t657 * t566;
t493 = -pkin(10) * t685 + t754;
t760 = t653 * t490 + t790 * t493;
t703 = -t657 * t625 - t777;
t534 = t703 - t787;
t535 = -pkin(10) * t692 + t749;
t758 = t653 * t534 + t790 * t535;
t706 = -t524 * t654 - t521;
t460 = t706 - t811;
t759 = t657 * t524 - t519;
t461 = t759 - t812;
t643 = pkin(3) * t657 + pkin(4);
t778 = t653 * t654;
t806 = t643 * t718 + (-t654 * t731 + (t657 * t790 - t778) * qJD(4)) * pkin(3) - t653 * t460 - t790 * t461;
t724 = t790 * t654;
t805 = t643 * t731 - (-t654 * t718 + (-t653 * t657 - t724) * qJD(4)) * pkin(3) + t790 * t460 - t461 * t653;
t803 = qJD(5) * t758 + t653 * t822 + t790 * t823;
t740 = qJD(2) * t611;
t676 = t659 * t740;
t662 = t656 * t817 - t676;
t802 = t534 * t718 - t535 * t731 - t653 * t823 + t790 * t822;
t801 = t432 * MDP(30);
t798 = -t573 * t694 + t672;
t797 = -t573 * t695 + t698;
t786 = t487 * t611;
t549 = t571 * pkin(3) + pkin(7) * t716;
t785 = t549 * t611;
t784 = t570 * t655;
t783 = t608 * t638;
t782 = t610 * t638;
t781 = t623 * t655;
t780 = t623 * t658;
t779 = t638 * t658;
t660 = qJD(2) ^ 2;
t772 = t656 * t660;
t770 = t659 * t660;
t661 = qJD(1) ^ 2;
t769 = t659 * t661;
t557 = t611 * t653 + t677;
t768 = -qJ(6) * t761 - qJD(6) * t557 + t802;
t767 = -pkin(5) * t742 + qJ(6) * t762 - t558 * qJD(6) - t803;
t766 = t425 - t426;
t765 = t790 * t454 - t451;
t757 = t806 + t808;
t756 = -t805 - t818;
t751 = t655 * t618 + t620 * t735;
t738 = qJD(2) * t656;
t750 = t658 * t618 + t738 * t788;
t619 = pkin(3) * t775 + t656 * pkin(7);
t733 = qJD(4) * t655;
t506 = t690 * qJD(2) + (-t640 + (pkin(9) * t656 - t620) * t655) * qJD(3) + t750;
t720 = t659 * t736;
t511 = -t678 * pkin(9) + (-t656 * t730 - t720) * pkin(7) + t751;
t727 = t654 * t506 + t657 * t511 + t560 * t732;
t574 = pkin(3) * t678 + pkin(7) * t737;
t644 = -t658 * pkin(3) - pkin(2);
t714 = MDP(15) * t738;
t713 = pkin(1) * t813;
t712 = -t454 * t653 - t453;
t707 = t790 * t490 - t493 * t653;
t705 = t790 * t534 - t535 * t653;
t702 = t608 + t730;
t701 = -t610 + t739;
t697 = -pkin(3) * t778 + t790 * t643;
t518 = pkin(3) * t610 + pkin(4) * t694;
t691 = qJD(1) * t651 - t638 * t659;
t686 = t573 * t611;
t680 = t692 * qJD(2);
t512 = t656 * t666 - t659 * t680;
t673 = -qJD(4) * t754 + t657 * t506 - t511 * t654;
t436 = pkin(4) * t738 - pkin(10) * t512 + t673;
t679 = t659 * t730 - t721;
t438 = (-t773 * t800 - t722) * pkin(10) * t657 + (-qJD(4) * t566 + (t656 * t733 - t679) * pkin(10)) * t654 + t727;
t684 = t653 * t436 + t790 * t438 + t490 * t718 - t493 * t731;
t675 = t790 * t685;
t579 = pkin(4) * t692 + t644;
t563 = pkin(4) * t685 + t619;
t462 = pkin(4) * t689 + t549;
t528 = -t584 * t790 - t653 * t685;
t668 = -qJD(5) * t760 + t790 * t436 - t653 * t438;
t423 = t441 * pkin(5) + t462;
t664 = t611 * t800 - t577;
t494 = -pkin(4) * t662 + t574;
t642 = pkin(4) * t790 + pkin(5);
t590 = pkin(3) * t724 + t653 * t643;
t585 = pkin(5) + t697;
t527 = -t584 * t653 + t675;
t464 = -qJ(6) * t557 + t758;
t463 = -qJ(6) * t558 + t705;
t446 = qJD(5) * t528 + t653 * t512 - t662 * t790;
t445 = qJD(5) * t675 - t512 * t790 - t584 * t731 - t653 * t662;
t444 = -qJ(6) * t527 + t760;
t443 = -pkin(5) * t659 - qJ(6) * t528 + t707;
t429 = t765 - t808;
t428 = t712 + t818;
t416 = -qJ(6) * t446 - qJD(6) * t527 + t684;
t415 = pkin(5) * t738 + t445 * qJ(6) - t528 * qJD(6) + t668;
t414 = -qJ(6) * t441 - qJD(6) * t498 - t699;
t413 = pkin(5) * t717 + t440 * qJ(6) - qJD(6) * t795 + t669;
t1 = [(t638 * t719 + t571 * t659 + (-t608 * t656 - t655 * t691) * qJD(2)) * MDP(14) + (t638 * t721 - t570 * t659 + (t610 * t656 + t658 * t691) * qJD(2)) * MDP(13) + (t512 * t695 + t584 * t689 - t694 * t676 + (-t786 + t694 * (-t657 * t735 - t658 * t732 + (t733 + t736) * t654)) * t656) * MDP(19) + ((-qJD(1) * t584 + t694) * MDP(20) + (qJD(1) * t528 + t795) * MDP(27) + (-qJD(1) * t527 - t498) * MDP(28) + (-qJD(1) * t754 - t468) * MDP(24) + t801 + (-qJD(1) * t760 - t433) * MDP(31) + (-t629 - t741) * MDP(22) + (-t621 - t741) * MDP(29)) * t738 + (-t563 * t440 - t514 * t445 + t462 * t528 + t494 * t795 + t684 * t621 - t699 * t659) * MDP(31) + (t440 * t527 - t441 * t528 + t445 * t498 - t446 * t795) * MDP(26) + (-t440 * t528 - t445 * t795) * MDP(25) + (-t413 * t528 - t414 * t527 - t415 * t795 - t416 * t498 + t425 * t445 - t427 * t446 + t440 * t443 - t441 * t444) * MDP(32) + (-t673 * t629 - t574 * t695 + t619 * t689 + (qJD(2) * t686 - t672) * t659 + (t785 - t573 * t817 + (qJD(1) * t704 + t467) * qJD(2)) * t656) * MDP(23) + ((t629 * t740 + t689) * t659 + (-t817 * t629 + (-qJD(1) * t685 + t695) * qJD(2)) * t656) * MDP(21) + (t563 * t441 + t514 * t446 + t462 * t527 + t494 * t498 - t621 * t668 - t659 * t669 + t707 * t717) * MDP(30) + 0.2e1 * t716 * t810 + (-pkin(7) * t770 + t656 * t713) * MDP(9) + ((-t608 * t658 - t610 * t655) * t737 + (-t784 - t571 * t658 + (t608 * t655 - t610 * t658) * qJD(3)) * t656) * MDP(12) + (-t487 * t584 + t512 * t694) * MDP(18) + ((-t566 * t734 + t727) * t629 - t698 * t659 + t574 * t694 + t619 * t487 - t549 * t584 + t573 * t512) * MDP(24) + (-t487 * t659 - t512 * t629) * MDP(20) + (t441 * t659 + t446 * t621) * MDP(28) + (t440 * t659 + t445 * t621) * MDP(27) + ((-pkin(7) * t720 + t751) * t638 + t688 * t659 + (pkin(7) * t570 - t623 * t736) * t656 + ((pkin(7) * t610 + t780) * t659 + (-pkin(7) * t779 - qJD(1) * t746 - t562) * t656) * qJD(2)) * MDP(17) + (-(-t620 * t736 + t750) * t638 + (t623 * t735 + pkin(7) * t571 + (t607 * qJD(1) + t561) * qJD(2)) * t656 + ((pkin(7) * t608 + t781) * qJD(2) + (t776 + (pkin(7) * t638 + t624) * t658) * qJD(3) + t752) * t659) * MDP(16) + (-t638 - t741) * t714 - MDP(7) * t772 + (pkin(7) * t772 + t659 * t713) * MDP(10) + (t414 * t444 + t427 * t416 + t413 * t443 + t425 * t415 + t423 * (t527 * pkin(5) + t563) + t459 * (t446 * pkin(5) + t494)) * MDP(33) + (t570 * t773 + t610 * t679) * MDP(11) + t809 * t813 + MDP(6) * t770; (t579 * t441 + t462 * t557 + t498 * t816 + t514 * t761 + t705 * t717) * MDP(30) + (t414 * t464 + t413 * t463 + t423 * (t557 * pkin(5) + t579) + (pkin(4) * t664 + pkin(5) * t761 + t804) * t459 + t768 * t427 + t767 * t425) * MDP(33) + (-t413 * t558 - t414 * t557 + t425 * t762 - t427 * t761 + t440 * t463 - t441 * t464 - t498 * t768 - t767 * t795) * MDP(32) + ((t570 + t783) * t658 + (-t571 + t782) * t655) * MDP(12) + (-t610 * t779 + t784) * MDP(11) + (t644 * t487 + t753 * t573 + t694 * t804 + t785) * MDP(24) + (t644 * t689 + t549 * t692 + t604 * t695 - t573 * t577 + t686 * qJD(4) + (-t695 * t789 + t686) * qJD(3)) * MDP(23) + (-t487 * t692 - t611 * t689 + t694 * t825 + t753 * t695) * MDP(19) + (t694 * t753 + t786) * MDP(18) + (-pkin(2) * t570 - t595 * t638 + (-pkin(8) * t638 * t655 + t780) * qJD(3) + (-t623 * t771 + (-pkin(8) * t730 + t562) * t656 + (t638 * t773 + t659 * t701) * pkin(7)) * qJD(1)) * MDP(17) + (-pkin(2) * t571 + t748 * t638 + (pkin(8) * t779 + t781) * qJD(3) + ((-pkin(8) * t739 - t561) * t656 + (-pkin(7) * t702 - t781) * t659) * qJD(1)) * MDP(16) + (-t638 * t735 + (t638 * t771 + t656 * t701) * qJD(1)) * MDP(13) + (t440 * t557 - t441 * t558 + t498 * t762 - t761 * t795) * MDP(26) + (-t440 * t558 - t762 * t795) * MDP(25) + (-t579 * t440 + t462 * t558 - t762 * t514 + t816 * t795) * MDP(31) + (t638 * t736 + (-t638 * t774 + t656 * t702) * qJD(1)) * MDP(14) - t769 * t810 + t661 * t809 + ((-t626 * t734 + t824) * MDP(24) + t827 * MDP(23) - t753 * MDP(20) + t664 * MDP(21)) * t629 + (MDP(9) * t656 * t661 + MDP(10) * t769) * pkin(1) + (t762 * MDP(27) + t761 * MDP(28) + MDP(30) * t803 + MDP(31) * t802) * t621 + (t638 * MDP(15) + t629 * MDP(22) + t621 * MDP(29) - t801 + (-qJD(2) * t749 + t468) * MDP(24) + (qJD(2) * t703 - t467) * MDP(23) + (-t694 + t740) * MDP(20) + (-t680 - t695) * MDP(21) + (-qJD(2) * t557 + t498) * MDP(28) + (qJD(2) * t558 - t795) * MDP(27) + (-qJD(2) * t758 + t433) * MDP(31)) * t742; (-t518 * t498 + t621 * t805 + t697 * t717 + t799) * MDP(30) + t610 * t608 * MDP(11) + (-t518 * t795 - t590 * t717 + t621 * t806 + t815) * MDP(31) + (t440 * t585 - t441 * t590 - t498 * t757 - t756 * t795 + t821) * MDP(32) + (t706 * t629 + (t610 * t695 + t629 * t734 + t657 * t717) * pkin(3) + t798) * MDP(23) + (-t759 * t629 + (-t610 * t694 + t629 * t732 - t654 * t717) * pkin(3) + t797) * MDP(24) + (-t571 - t782) * MDP(14) + (t570 - t783) * MDP(13) + qJD(1) * t714 + (-t561 * t638 + t608 * t623 - t671) * MDP(17) + (-t562 * t638 - t610 * t623 + t674) * MDP(16) + (t414 * t590 + t413 * t585 - t459 * (pkin(5) * t795 + t518) + t757 * t427 + t756 * t425) * MDP(33) + (-t608 ^ 2 + t610 ^ 2) * MDP(12) + t828; (-t468 * t629 + t798) * MDP(23) + (-t467 * t629 + t797) * MDP(24) + (t712 * t621 + (-t498 * t694 + t621 * t731 + t717 * t790) * pkin(4) + t799) * MDP(30) + (-t765 * t621 + (t621 * t718 - t653 * t717 - t694 * t795) * pkin(4) + t815) * MDP(31) + (t428 * t795 + t429 * t498 + t642 * t440 + (-t441 * t653 + (-t498 * t790 + t653 * t795) * qJD(5)) * pkin(4) + t821) * MDP(32) + (-pkin(5) * t807 + t413 * t642 - t425 * t428 - t427 * t429 + (t414 * t653 - t459 * t694 + (-t425 * t653 + t427 * t790) * qJD(5)) * pkin(4)) * MDP(33) + t828; (-t433 * t621 + t799) * MDP(30) + (-t432 * t621 + t815) * MDP(31) + (pkin(5) * t440 - t498 * t766) * MDP(32) + (t766 * t427 + (t413 - t807) * pkin(5)) * MDP(33) + t829; (-t495 - t792) * MDP(32) + (t425 * t795 + t427 * t498 + t423) * MDP(33);];
tauc  = t1;

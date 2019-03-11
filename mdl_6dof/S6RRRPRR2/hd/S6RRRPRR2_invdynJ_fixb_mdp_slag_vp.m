% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR2
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
%   see S6RRRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:08:58
% EndTime: 2019-03-09 18:09:13
% DurationCPUTime: 10.33s
% Computational Cost: add. (12109->607), mult. (28851->781), div. (0->0), fcn. (22647->18), ass. (0->300)
t730 = qJ(2) + qJ(3);
t717 = pkin(11) + t730;
t704 = sin(t717);
t736 = sin(qJ(1));
t741 = cos(qJ(1));
t786 = g(1) * t741 + g(2) * t736;
t915 = t786 * t704;
t734 = sin(qJ(3));
t735 = sin(qJ(2));
t739 = cos(qJ(3));
t740 = cos(qJ(2));
t666 = t734 * t735 - t739 * t740;
t648 = t666 * qJD(1);
t831 = qJD(1) * t735;
t846 = t734 * t740;
t650 = -qJD(1) * t846 - t739 * t831;
t731 = sin(pkin(11));
t873 = cos(pkin(11));
t797 = -t873 * t648 + t650 * t731;
t891 = qJD(5) + qJD(6);
t914 = t797 - t891;
t766 = -t731 * t648 - t650 * t873;
t888 = pkin(3) * t650;
t558 = pkin(4) * t766 - pkin(9) * t797 - t888;
t715 = pkin(2) * t831;
t555 = t558 + t715;
t733 = sin(qJ(5));
t738 = cos(qJ(5));
t889 = pkin(7) + pkin(8);
t686 = t889 * t740;
t678 = qJD(1) * t686;
t655 = t739 * t678;
t685 = t889 * t735;
t676 = qJD(1) * t685;
t795 = t676 * t734 - t655;
t872 = qJ(4) * t648;
t595 = t795 + t872;
t642 = t650 * qJ(4);
t651 = t734 * t678;
t836 = -t739 * t676 - t651;
t596 = t642 + t836;
t851 = t731 * t734;
t875 = pkin(2) * qJD(3);
t837 = -t731 * t595 - t596 * t873 + (t739 * t873 - t851) * t875;
t913 = -t738 * t555 - t733 * t837;
t828 = qJD(5) * t733;
t859 = t797 * t733;
t912 = t828 - t859;
t723 = t740 * pkin(2);
t878 = pkin(1) + t723;
t789 = pkin(3) * t666 - t878;
t719 = sin(t730);
t721 = cos(t730);
t911 = -g(3) * t721 + t719 * t786;
t684 = qJD(1) * t878;
t910 = qJDD(1) * t878;
t726 = qJD(2) + qJD(3);
t592 = -t738 * t726 + t733 * t766;
t737 = cos(qJ(6));
t594 = t726 * t733 + t738 * t766;
t732 = sin(qJ(6));
t860 = t594 * t732;
t547 = t737 * t592 + t860;
t892 = qJD(5) - t797;
t599 = qJD(6) + t892;
t909 = t547 * t599;
t908 = t592 * t892;
t775 = t592 * t732 - t737 * t594;
t907 = t599 * t775;
t668 = t735 * t739 + t846;
t617 = t726 * t668;
t820 = qJDD(1) * t740;
t823 = qJD(1) * qJD(2);
t808 = t740 * t823;
t821 = qJDD(1) * t735;
t764 = -t808 - t821;
t822 = qJD(1) * qJD(3);
t901 = -t740 * t822 + t764;
t588 = (-t726 * t831 + t820) * t734 - t901 * t739;
t724 = qJDD(2) + qJDD(3);
t618 = qJDD(2) * pkin(2) + t764 * t889;
t809 = t735 * t823;
t763 = -t809 + t820;
t622 = t889 * t763;
t874 = qJD(2) * pkin(2);
t657 = -t676 + t874;
t774 = -t657 * t734 - t655;
t755 = qJD(3) * t774 + t739 * t618 - t734 * t622;
t512 = pkin(3) * t724 - qJ(4) * t588 + qJD(4) * t650 + t755;
t830 = qJD(3) * t734;
t646 = t678 * t830;
t792 = -qJD(3) * t657 - t622;
t518 = -t648 * qJD(4) - t646 + (qJ(4) * t901 + t618) * t734 + ((-t735 * t822 + t763) * qJ(4) - t792) * t739;
t477 = t512 * t873 - t731 * t518;
t475 = -t724 * pkin(4) - t477;
t705 = cos(t717);
t904 = g(3) * t705 + t475;
t850 = t732 * t738;
t667 = t733 * t737 + t850;
t903 = t914 * t667;
t665 = t732 * t733 - t737 * t738;
t902 = t914 * t665;
t793 = t892 * t738;
t749 = qJD(1) * t617;
t745 = -t666 * qJDD(1) - t749;
t539 = -t588 * t731 + t873 * t745;
t536 = qJDD(5) - t539;
t849 = t733 * t536;
t899 = -t892 * t793 - t849;
t796 = t739 * t657 - t651;
t586 = t642 + t796;
t577 = pkin(3) * t726 + t586;
t587 = -t774 - t872;
t802 = t873 * t587;
t528 = t731 * t577 + t802;
t526 = pkin(9) * t726 + t528;
t619 = pkin(3) * t648 + qJD(4) - t684;
t542 = -pkin(4) * t797 - pkin(9) * t766 + t619;
t498 = t526 * t738 + t542 * t733;
t493 = -pkin(10) * t592 + t498;
t826 = qJD(6) * t732;
t490 = t493 * t826;
t578 = t731 * t587;
t527 = t577 * t873 - t578;
t525 = -t726 * pkin(4) - t527;
t517 = t592 * pkin(5) + t525;
t729 = qJ(5) + qJ(6);
t720 = cos(t729);
t853 = t720 * t736;
t718 = sin(t729);
t854 = t718 * t741;
t625 = -t705 * t853 + t854;
t852 = t720 * t741;
t855 = t718 * t736;
t627 = t705 * t852 + t855;
t881 = g(3) * t720;
t898 = g(1) * t627 - g(2) * t625 + t517 * t547 + t704 * t881 + t490;
t624 = t705 * t855 + t852;
t626 = -t705 * t854 + t853;
t478 = t731 * t512 + t873 * t518;
t476 = pkin(9) * t724 + t478;
t540 = t873 * t588 + t731 * t745;
t702 = pkin(2) * t809;
t744 = pkin(3) * t749 + qJDD(1) * t789 + qJDD(4) + t702;
t491 = -t539 * pkin(4) - t540 * pkin(9) + t744;
t489 = t738 * t491;
t827 = qJD(5) * t738;
t515 = t738 * t540 + t733 * t724 + t726 * t827 - t766 * t828;
t458 = pkin(5) * t536 - pkin(10) * t515 - qJD(5) * t498 - t476 * t733 + t489;
t798 = t540 * t733 - t738 * t724;
t516 = qJD(5) * t594 + t798;
t761 = t738 * t476 + t733 * t491 - t526 * t828 + t542 * t827;
t459 = -pkin(10) * t516 + t761;
t800 = t737 * t458 - t732 * t459;
t882 = g(3) * t718;
t897 = -g(1) * t626 + g(2) * t624 + t517 * t775 + t704 * t882 + t800;
t534 = qJDD(6) + t536;
t896 = t534 * MDP(31) + (-t547 ^ 2 + t775 ^ 2) * MDP(28) - t547 * MDP(27) * t775;
t612 = -t731 * t666 + t668 * t873;
t569 = t667 * t612;
t801 = t873 * t734;
t838 = t873 * t595 - t596 * t731 + (t731 * t739 + t801) * t875;
t835 = -t734 * t685 + t739 * t686;
t894 = t912 * pkin(5);
t893 = t733 * t555 - t738 * t837;
t890 = t534 * t667 + t599 * t902;
t799 = t515 * t732 + t737 * t516;
t472 = -qJD(6) * t775 + t799;
t884 = g(3) * t704;
t879 = t738 * pkin(5);
t722 = t738 * pkin(10);
t713 = pkin(2) * t739 + pkin(3);
t641 = pkin(2) * t801 + t731 * t713;
t636 = pkin(9) + t641;
t877 = -pkin(10) - t636;
t706 = pkin(3) * t731 + pkin(9);
t876 = -pkin(10) - t706;
t497 = -t526 * t733 + t738 * t542;
t492 = -pkin(10) * t594 + t497;
t482 = pkin(5) * t892 + t492;
t871 = t482 * t737;
t870 = t493 * t737;
t869 = t515 * t733;
t868 = t525 * t797;
t866 = t547 * t766;
t865 = t775 * t766;
t616 = t726 * t666;
t575 = -t616 * t873 - t731 * t617;
t864 = t575 * t733;
t863 = t575 * t738;
t862 = t592 * t766;
t861 = t594 * t766;
t858 = t612 * t733;
t857 = t612 * t738;
t848 = t733 * t736;
t847 = t733 * t741;
t845 = t736 * t738;
t530 = t738 * t536;
t794 = -t739 * t685 - t686 * t734;
t604 = -qJ(4) * t668 + t794;
t605 = -qJ(4) * t666 + t835;
t567 = t731 * t604 + t605 * t873;
t562 = t738 * t567;
t844 = t738 * t741;
t538 = t586 * t873 - t578;
t843 = t738 * t538 + t733 * t558;
t611 = t666 * t873 + t668 * t731;
t565 = pkin(4) * t611 - pkin(9) * t612 + t789;
t841 = t733 * t565 + t562;
t839 = t894 + t838;
t834 = pkin(3) * t721 + t723;
t727 = t735 ^ 2;
t833 = -t740 ^ 2 + t727;
t829 = qJD(3) * t739;
t825 = qJD(6) * t737;
t818 = pkin(10) * t859;
t716 = t735 * t874;
t814 = t737 * t515 - t732 * t516 - t592 * t825;
t812 = qJD(2) * t889;
t811 = t612 * t828;
t810 = qJD(5) * t706 * t892;
t523 = t525 * t827;
t806 = pkin(3) * t617 + t716;
t804 = qJD(5) * t877;
t803 = qJD(5) * t876;
t677 = t735 * t812;
t679 = t740 * t812;
t762 = -t739 * t677 - t734 * t679 - t685 * t829 - t686 * t830;
t551 = -qJ(4) * t617 - qJD(4) * t666 + t762;
t754 = -qJD(3) * t835 + t677 * t734 - t739 * t679;
t552 = qJ(4) * t616 - qJD(4) * t668 + t754;
t501 = t551 * t731 - t873 * t552;
t537 = t586 * t731 + t802;
t566 = -t873 * t604 + t605 * t731;
t791 = -qJD(5) * t542 - t476;
t790 = qJD(6) * t482 + t459;
t788 = pkin(5) * t766 - t722 * t797;
t707 = -pkin(3) * t873 - pkin(4);
t787 = -t537 + t894;
t785 = g(1) * t736 - g(2) * t741;
t784 = -t665 * t534 + t599 * t903;
t783 = -t526 * t827 + t489;
t557 = t738 * t558;
t661 = t706 * t738 + t722;
t782 = qJD(6) * t661 - t538 * t733 - t738 * t803 + t557 + t788;
t621 = t636 * t738 + t722;
t781 = qJD(6) * t621 - t738 * t804 + t788 - t913;
t660 = t876 * t733;
t780 = -qJD(6) * t660 - t733 * t803 - t818 + t843;
t620 = t877 * t733;
t779 = -qJD(6) * t620 - t733 * t804 - t818 + t893;
t640 = -pkin(2) * t851 + t713 * t873;
t465 = t482 * t732 + t870;
t778 = -t536 * t636 - t868;
t777 = -t536 * t706 - t868;
t776 = t527 * t797 + t528 * t766;
t773 = t498 * t766 + t733 * t904 + t523;
t772 = -t497 * t766 + t525 * t828 + t738 * t915;
t635 = -pkin(4) - t640;
t771 = -t892 * t912 + t530;
t769 = -0.2e1 * pkin(1) * t823 - pkin(7) * qJDD(2);
t768 = t612 * t827 + t864;
t767 = -t811 + t863;
t502 = t551 * t873 + t731 * t552;
t574 = -t616 * t731 + t617 * t873;
t508 = pkin(4) * t574 - pkin(9) * t575 + t806;
t760 = t738 * t502 + t733 * t508 + t565 * t827 - t567 * t828;
t471 = -t594 * t826 + t814;
t742 = qJD(2) ^ 2;
t757 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t742 + t785;
t743 = qJD(1) ^ 2;
t756 = pkin(1) * t743 - pkin(7) * qJDD(1) + t786;
t464 = -t493 * t732 + t871;
t468 = t516 * pkin(5) + t475;
t752 = -t464 * t766 + t468 * t665 - t517 * t903 - t705 * t881 + t720 * t915;
t751 = g(3) * t719 - t734 * t618 - t684 * t648 + t721 * t786 + t792 * t739 + t646;
t750 = t465 * t766 + t468 * t667 + t517 * t902 + t705 * t882 - t718 * t915;
t747 = -t684 * t650 + t755 + t911;
t746 = -t650 * t648 * MDP(11) + (-t471 * t665 - t472 * t667 - t547 * t902 - t775 * t903) * MDP(28) + (t471 * t667 - t775 * t902) * MDP(27) + ((t515 - t908) * t738 + (-t594 * t892 - t516) * t733) * MDP(21) + (t865 + t890) * MDP(29) + (t784 + t866) * MDP(30) + (t771 + t862) * MDP(23) + (-t861 - t899) * MDP(22) + (t594 * t793 + t869) * MDP(20) + (t648 * t726 + t588) * MDP(13) + (-t650 * t726 + t745) * MDP(14) + (-t648 ^ 2 + t650 ^ 2) * MDP(12) + t724 * MDP(15) + (-MDP(24) * t892 - MDP(31) * t599) * t766;
t725 = -qJ(4) - t889;
t681 = t707 - t879;
t675 = pkin(1) + t834;
t645 = t702 - t910;
t634 = t705 * t844 + t848;
t633 = -t705 * t847 + t845;
t632 = -t705 * t845 + t847;
t631 = t705 * t848 + t844;
t628 = t635 - t879;
t570 = t665 * t612;
t561 = t738 * t565;
t524 = pkin(5) * t858 + t566;
t506 = t738 * t508;
t499 = -pkin(10) * t858 + t841;
t496 = pkin(5) * t611 - pkin(10) * t857 - t567 * t733 + t561;
t485 = t575 * t850 - t732 * t811 - t826 * t858 + (t857 * t891 + t864) * t737;
t484 = -t569 * t891 - t665 * t575;
t481 = pkin(5) * t768 + t501;
t462 = -pkin(10) * t768 + t760;
t461 = -pkin(10) * t863 + pkin(5) * t574 - t502 * t733 + t506 + (-t562 + (pkin(10) * t612 - t565) * t733) * qJD(5);
t1 = [(-t471 * t569 + t472 * t570 - t484 * t547 + t485 * t775) * MDP(28) + (-t471 * t570 - t484 * t775) * MDP(27) + (t471 * t611 + t484 * t599 - t534 * t570 - t574 * t775) * MDP(29) + (-g(1) * t624 - g(2) * t626 - t465 * t574 - t468 * t570 + t524 * t471 - t481 * t775 + t517 * t484 + t490 * t611 + (-(-qJD(6) * t499 + t461) * t599 - t496 * t534 - t458 * t611) * t732 + (-(qJD(6) * t496 + t462) * t599 - t499 * t534 - t790 * t611) * t737) * MDP(33) + (t536 * t611 + t574 * t892) * MDP(24) + (-g(1) * t631 - g(2) * t633 + t475 * t857 - t498 * t574 + t501 * t594 + t566 * t515 + t525 * t767 - t536 * t841 - t611 * t761 - t760 * t892) * MDP(26) + (-t516 * t611 - t574 * t592 - t612 * t849 - t768 * t892) * MDP(23) + (t515 * t611 + t530 * t612 + t574 * t594 + t767 * t892) * MDP(22) + ((-t567 * t827 + t506) * t892 + t561 * t536 + t783 * t611 + t497 * t574 + t501 * t592 + t566 * t516 + t612 * t523 - g(1) * t632 - g(2) * t634 + ((-qJD(5) * t565 - t502) * t892 - t567 * t536 + t791 * t611 + t475 * t612 + t525 * t575) * t733) * MDP(25) + (-t616 * t726 + t668 * t724) * MDP(13) + (-t588 * t666 + t616 * t648 + t650 * t617 + t668 * t745) * MDP(12) + (t588 * t668 + t616 * t650) * MDP(11) + (-t477 * t612 - t478 * t611 + t501 * t766 + t502 * t797 - t527 * t575 - t528 * t574 + t539 * t567 + t540 * t566 - t786) * MDP(18) + (-t617 * t726 - t666 * t724) * MDP(14) + (t648 * t716 + t721 * t785 + t724 * t794 + t726 * t754 + (t645 - t910) * t666 - 0.2e1 * t684 * t617) * MDP(16) + (qJDD(2) * t735 + t740 * t742) * MDP(6) + (qJDD(2) * t740 - t735 * t742) * MDP(7) + (t534 * t611 + t574 * t599) * MDP(31) + (-t472 * t611 - t485 * t599 - t534 * t569 - t547 * t574) * MDP(30) + (-t588 * t878 + t684 * t616 + t645 * t668 - t650 * t716 - t719 * t785 - t724 * t835 - t726 * t762) * MDP(17) + ((-t592 * t738 - t594 * t733) * t575 + (-t869 - t516 * t738 + (t592 * t733 - t594 * t738) * qJD(5)) * t612) * MDP(21) + (t515 * t857 + t594 * t767) * MDP(20) + 0.2e1 * (t735 * t820 - t823 * t833) * MDP(5) + (qJDD(1) * t727 + 0.2e1 * t735 * t808) * MDP(4) + (t478 * t567 + t528 * t502 - t477 * t566 - t527 * t501 + t744 * t789 + t619 * t806 - g(1) * (-t675 * t736 - t725 * t741) - g(2) * (t675 * t741 - t725 * t736)) * MDP(19) + ((t461 * t737 - t462 * t732) * t599 + (t496 * t737 - t499 * t732) * t534 + t800 * t611 + t464 * t574 + t481 * t547 + t524 * t472 + t468 * t569 + t517 * t485 - g(1) * t625 - g(2) * t627 + ((-t496 * t732 - t499 * t737) * t599 - t465 * t611) * qJD(6)) * MDP(32) + qJDD(1) * MDP(1) + t785 * MDP(2) + t786 * MDP(3) + (t735 * t769 + t740 * t757) * MDP(9) + (-t735 * t757 + t740 * t769) * MDP(10); (t635 * t515 + t778 * t738 - t733 * t915 + t838 * t594 + (t636 * t828 + t893) * t892 + t773) * MDP(26) + (t836 * t726 + (t650 * t831 - t724 * t734 - t726 * t829) * pkin(2) + t751) * MDP(17) + (-t795 * t726 + (-t648 * t831 + t724 * t739 - t726 * t830) * pkin(2) + t747) * MDP(16) + (t478 * t641 + t477 * t640 - t619 * (t715 - t888) - g(3) * t834 - t786 * (-pkin(2) * t735 - pkin(3) * t719) + t837 * t528 - t838 * t527) * MDP(19) + (t539 * t641 - t540 * t640 + t766 * t838 + t797 * t837 + t776) * MDP(18) + (g(3) * t735 + t740 * t756) * MDP(10) + (-(t620 * t732 + t621 * t737) * t534 + t628 * t471 + (t732 * t781 + t737 * t779) * t599 - t839 * t775 + t750) * MDP(33) + ((t620 * t737 - t621 * t732) * t534 + t628 * t472 + (t732 * t779 - t737 * t781) * t599 + t839 * t547 + t752) * MDP(32) + (t635 * t516 - t904 * t738 + t778 * t733 + t838 * t592 + (-t636 * t827 + t913) * t892 + t772) * MDP(25) + MDP(7) * t820 + MDP(6) * t821 + qJDD(2) * MDP(8) + t746 + (-g(3) * t740 + t735 * t756) * MDP(9) + (-MDP(4) * t735 * t740 + MDP(5) * t833) * t743; (t726 * t796 + t751) * MDP(17) + ((t660 * t737 - t661 * t732) * t534 + t681 * t472 + (t732 * t780 - t737 * t782) * t599 + t787 * t547 + t752) * MDP(32) + (t707 * t515 + t843 * t892 - t537 * t594 + t777 * t738 + (-t915 + t810) * t733 + t773) * MDP(26) + (-(t660 * t732 + t661 * t737) * t534 + t681 * t471 + (t732 * t782 + t737 * t780) * t599 - t787 * t775 + t750) * MDP(33) + (t707 * t516 - t537 * t592 - t557 * t892 + (t538 * t892 + t777) * t733 + (-t904 - t810) * t738 + t772) * MDP(25) + (-t726 * t774 + t747) * MDP(16) + (t527 * t537 - t528 * t538 + (t477 * t873 + t478 * t731 + t619 * t650 + t911) * pkin(3)) * MDP(19) + t746 + (-t537 * t766 - t538 * t797 + (t539 * t731 - t540 * t873) * pkin(3) + t776) * MDP(18); (-t766 ^ 2 - t797 ^ 2) * MDP(18) + (t771 - t862) * MDP(25) + (-t861 + t899) * MDP(26) + (t784 - t866) * MDP(32) + (t865 - t890) * MDP(33) + (t527 * t766 - t528 * t797 + t744 - t785) * MDP(19); t594 * t592 * MDP(20) + (-t592 ^ 2 + t594 ^ 2) * MDP(21) + (t515 + t908) * MDP(22) + (-t798 + (-qJD(5) + t892) * t594) * MDP(23) + t536 * MDP(24) + (-g(1) * t633 + g(2) * t631 + t498 * t892 - t525 * t594 + (t791 + t884) * t733 + t783) * MDP(25) + (g(1) * t634 - g(2) * t632 + t497 * t892 + t525 * t592 + t738 * t884 - t761) * MDP(26) + (t471 + t909) * MDP(29) + (-t472 - t907) * MDP(30) + (-(-t492 * t732 - t870) * t599 - t465 * qJD(6) + (t534 * t737 - t547 * t594 - t599 * t826) * pkin(5) + t897) * MDP(32) + ((-t493 * t599 - t458) * t732 + (t492 * t599 - t790) * t737 + (-t534 * t732 + t594 * t775 - t599 * t825) * pkin(5) + t898) * MDP(33) + t896; (t814 + t909) * MDP(29) + (-t799 - t907) * MDP(30) + (t465 * t599 + t897) * MDP(32) + (-t732 * t458 - t737 * t459 + t464 * t599 + t898) * MDP(33) + (-MDP(29) * t860 + MDP(30) * t775 - MDP(32) * t465 - MDP(33) * t871) * qJD(6) + t896;];
tau  = t1;

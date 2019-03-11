% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:19:32
% EndTime: 2019-03-09 13:19:51
% DurationCPUTime: 13.83s
% Computational Cost: add. (11544->604), mult. (27487->772), div. (0->0), fcn. (22252->16), ass. (0->282)
t733 = sin(pkin(11));
t734 = cos(pkin(11));
t739 = sin(qJ(2));
t743 = cos(qJ(2));
t679 = -t733 * t739 + t734 * t743;
t667 = t679 * qJD(1);
t680 = t733 * t743 + t734 * t739;
t669 = t680 * qJD(1);
t738 = sin(qJ(4));
t877 = cos(qJ(4));
t616 = t877 * t667 - t669 * t738;
t881 = qJD(5) + qJD(6);
t912 = t616 - t881;
t736 = sin(qJ(6));
t737 = sin(qJ(5));
t741 = cos(qJ(6));
t742 = cos(qJ(5));
t686 = t736 * t737 - t741 * t742;
t908 = t912 * t686;
t810 = qJD(1) * qJD(2);
t795 = t743 * t810;
t796 = t739 * t810;
t625 = qJDD(1) * t680 - t733 * t796 + t734 * t795;
t668 = t680 * qJD(2);
t756 = qJD(1) * t668;
t747 = qJDD(1) * t679 - t756;
t765 = -t738 * t667 - t669 * t877;
t750 = qJD(4) * t765 - t738 * t625 + t877 * t747;
t534 = qJDD(5) - t750;
t531 = qJDD(6) + t534;
t811 = -qJD(5) + t616;
t603 = qJD(6) - t811;
t838 = t736 * t742;
t687 = t737 * t741 + t838;
t899 = t687 * t531 + t603 * t908;
t911 = t912 * t687;
t815 = qJD(5) * t742;
t891 = t616 * t742;
t910 = t815 - t891;
t816 = qJD(5) * t737;
t892 = t616 * t737;
t909 = t816 - t892;
t724 = qJ(2) + pkin(11) + qJ(4);
t715 = sin(t724);
t740 = sin(qJ(1));
t744 = cos(qJ(1));
t781 = g(1) * t744 + g(2) * t740;
t907 = t781 * t715;
t906 = pkin(10) * t892;
t570 = -pkin(4) * t765 - pkin(9) * t616;
t639 = t739 * qJD(1) * pkin(2) + pkin(3) * t669;
t555 = t570 + t639;
t864 = qJ(3) + pkin(7);
t701 = t864 * t739;
t688 = qJD(1) * t701;
t702 = t864 * t743;
t689 = qJD(1) * t702;
t839 = t734 * t689;
t629 = t688 * t733 - t839;
t874 = pkin(8) * t667;
t600 = t629 - t874;
t672 = t733 * t689;
t630 = -t734 * t688 - t672;
t873 = pkin(8) * t669;
t601 = t630 - t873;
t717 = pkin(2) * t734 + pkin(3);
t876 = pkin(2) * t733;
t778 = t877 * t717 - t738 * t876;
t824 = -t778 * qJD(4) + t738 * t600 + t601 * t877;
t905 = -t742 * t555 + t737 * t824;
t904 = t909 * pkin(5);
t903 = -pkin(5) * t765 - pkin(10) * t891;
t728 = qJDD(2) + qJDD(4);
t793 = qJD(2) * t864;
t665 = -qJD(3) * t739 - t743 * t793;
t622 = qJDD(2) * pkin(2) + qJD(1) * t665 - qJDD(1) * t701;
t664 = qJD(3) * t743 - t739 * t793;
t628 = qJD(1) * t664 + qJDD(1) * t702;
t573 = t734 * t622 - t628 * t733;
t538 = qJDD(2) * pkin(3) - pkin(8) * t625 + t573;
t574 = t733 * t622 + t734 * t628;
t543 = pkin(8) * t747 + t574;
t863 = qJD(2) * pkin(2);
t676 = -t688 + t863;
t623 = t734 * t676 - t672;
t587 = qJD(2) * pkin(3) + t623 - t873;
t624 = t733 * t676 + t839;
t594 = t624 + t874;
t797 = qJD(4) * t877;
t817 = qJD(4) * t738;
t783 = -t877 * t538 + t738 * t543 + t587 * t817 + t594 * t797;
t479 = -pkin(4) * t728 + t783;
t729 = qJD(2) + qJD(4);
t599 = t729 * t737 - t742 * t765;
t535 = t877 * t625 + t667 * t797 - t669 * t817 + t738 * t747;
t789 = t535 * t737 - t742 * t728;
t510 = qJD(5) * t599 + t789;
t470 = pkin(5) * t510 + t479;
t539 = t587 * t877 - t738 * t594;
t528 = -t729 * pkin(4) - t539;
t597 = -t742 * t729 - t737 * t765;
t514 = t597 * pkin(5) + t528;
t716 = cos(t724);
t732 = qJ(5) + qJ(6);
t725 = sin(t732);
t869 = g(3) * t725;
t902 = t470 * t687 + t514 * t908 + t716 * t869 - t725 * t907;
t726 = cos(t732);
t868 = g(3) * t726;
t901 = t470 * t686 - t514 * t911 - t716 * t868 + t726 * t907;
t509 = t742 * t535 + t737 * t728 + t729 * t815 + t765 * t816;
t900 = t509 * t742 - t737 * t510 - t597 * t910;
t779 = -t686 * t531 + t603 * t911;
t526 = t737 * t534;
t898 = -t811 * t910 + t526;
t813 = qJD(6) * t741;
t803 = t741 * t509 - t736 * t510 - t597 * t813;
t814 = qJD(6) * t736;
t474 = -t599 * t814 + t803;
t770 = t597 * t736 - t741 * t599;
t790 = t509 * t736 + t741 * t510;
t475 = -qJD(6) * t770 + t790;
t507 = t509 * t737;
t854 = t599 * t736;
t549 = t741 * t597 + t854;
t897 = t728 * MDP(17) + (t599 * t910 + t507) * MDP(20) + t474 * t687 * MDP(27) + (-t474 * t686 - t687 * t475 - t549 * t908) * MDP(28) + (-t908 * MDP(27) - MDP(28) * t911) * t770;
t895 = t528 * t616;
t894 = t549 * t603;
t893 = t603 * t770;
t870 = g(3) * t716;
t887 = t479 + t870;
t540 = t738 * t587 + t877 * t594;
t529 = pkin(9) * t729 + t540;
t720 = t743 * pkin(2) + pkin(1);
t695 = -qJD(1) * t720 + qJD(3);
t633 = -pkin(3) * t667 + t695;
t545 = -pkin(4) * t616 + pkin(9) * t765 + t633;
t496 = t529 * t742 + t545 * t737;
t491 = -pkin(10) * t597 + t496;
t487 = t491 * t814;
t841 = t726 * t740;
t842 = t725 * t744;
t642 = -t716 * t841 + t842;
t840 = t726 * t744;
t843 = t725 * t740;
t644 = t716 * t840 + t843;
t885 = g(1) * t644 - g(2) * t642 + t514 * t549 + t715 * t868 + t487;
t641 = t716 * t843 + t840;
t643 = -t716 * t842 + t841;
t749 = t738 * t538 + t543 * t877 + t587 * t797 - t594 * t817;
t478 = t728 * pkin(9) + t749;
t647 = -pkin(3) * t679 - t720;
t807 = pkin(2) * t796 + qJDD(3);
t591 = pkin(3) * t756 + qJDD(1) * t647 + t807;
t488 = -pkin(4) * t750 - t535 * pkin(9) + t591;
t486 = t742 * t488;
t462 = pkin(5) * t534 - pkin(10) * t509 - qJD(5) * t496 - t478 * t737 + t486;
t760 = t742 * t478 + t737 * t488 - t529 * t816 + t545 * t815;
t463 = -pkin(10) * t510 + t760;
t791 = t741 * t462 - t736 * t463;
t884 = -g(1) * t643 + g(2) * t641 + t514 * t770 + t715 * t869 + t791;
t883 = t531 * MDP(31) + (-t549 ^ 2 + t770 ^ 2) * MDP(28) - t549 * MDP(27) * t770;
t627 = t738 * t679 + t680 * t877;
t571 = t687 * t627;
t822 = t738 * t717 + t877 * t876;
t823 = t822 * qJD(4) + t600 * t877 - t738 * t601;
t882 = t737 * t555 + t742 * t824;
t880 = MDP(13) * t765 - MDP(14) * t616 - t633 * MDP(19);
t495 = -t529 * t737 + t742 * t545;
t490 = -pkin(10) * t599 + t495;
t481 = -pkin(5) * t811 + t490;
t862 = t481 * t741;
t466 = -t491 * t736 + t862;
t861 = t491 * t741;
t467 = t481 * t736 + t861;
t879 = -MDP(14) * t765 - MDP(18) * t633 + MDP(24) * t811 - MDP(25) * t495 + MDP(26) * t496 - MDP(31) * t603 - MDP(32) * t466 + MDP(33) * t467;
t878 = -pkin(9) - pkin(10);
t708 = g(3) * t715;
t867 = g(3) * t743;
t866 = t742 * pkin(5);
t663 = pkin(9) + t822;
t865 = -pkin(10) - t663;
t860 = t549 * t765;
t859 = t770 * t765;
t671 = t679 * qJD(2);
t764 = t679 * t877 - t680 * t738;
t575 = qJD(4) * t764 - t738 * t668 + t671 * t877;
t858 = t575 * t737;
t857 = t575 * t742;
t856 = t597 * t765;
t855 = t599 * t765;
t853 = t599 * t737;
t852 = t765 * t729;
t849 = t616 * t729;
t846 = t627 * t737;
t845 = t627 * t742;
t837 = t737 * t740;
t836 = t737 * t744;
t835 = t740 * t742;
t527 = t742 * t534;
t634 = -t734 * t701 - t702 * t733;
t606 = -pkin(8) * t680 + t634;
t635 = -t733 * t701 + t734 * t702;
t607 = pkin(8) * t679 + t635;
t566 = t738 * t606 + t607 * t877;
t559 = t742 * t566;
t834 = t742 * t744;
t828 = t742 * t539 + t737 * t570;
t560 = -pkin(4) * t764 - pkin(9) * t627 + t647;
t826 = t737 * t560 + t559;
t825 = t823 + t904;
t605 = t734 * t664 + t733 * t665;
t730 = t739 ^ 2;
t821 = -t743 ^ 2 + t730;
t809 = qJDD(1) * t739;
t808 = qJDD(1) * t743;
t723 = t739 * t863;
t806 = qJD(5) * pkin(9) * t811;
t524 = t528 * t815;
t804 = t737 * t887 + t524;
t801 = t528 * t816 + t742 * t907;
t800 = qJD(5) * t878;
t798 = t627 * t816;
t640 = pkin(3) * t668 + t723;
t792 = qJD(5) * t865;
t604 = -t664 * t733 + t734 * t665;
t787 = t811 * t737;
t786 = -qJD(5) * t545 - t478;
t785 = qJD(6) * t481 + t463;
t782 = -t540 + t904;
t780 = g(1) * t740 - g(2) * t744;
t777 = -t529 * t815 + t486;
t727 = t742 * pkin(10);
t637 = t663 * t742 + t727;
t776 = qJD(6) * t637 - t742 * t792 + t903 - t905;
t568 = t742 * t570;
t705 = pkin(9) * t742 + t727;
t775 = qJD(6) * t705 - t539 * t737 - t742 * t800 + t568 + t903;
t636 = t865 * t737;
t774 = -qJD(6) * t636 - t737 * t792 + t882 - t906;
t704 = t878 * t737;
t773 = -qJD(6) * t704 - t737 * t800 + t828 - t906;
t772 = -pkin(9) * t534 - t895;
t771 = -t534 * t663 - t895;
t662 = -pkin(4) - t778;
t769 = t811 * t909 + t527;
t767 = -0.2e1 * pkin(1) * t810 - pkin(7) * qJDD(2);
t766 = t606 * t877 - t738 * t607;
t763 = t627 * t815 + t858;
t762 = -t798 + t857;
t580 = -pkin(8) * t671 + t604;
t581 = -pkin(8) * t668 + t605;
t499 = qJD(4) * t766 + t738 * t580 + t581 * t877;
t576 = qJD(4) * t627 + t668 * t877 + t738 * t671;
t513 = pkin(4) * t576 - pkin(9) * t575 + t640;
t759 = t742 * t499 + t737 * t513 + t560 * t815 - t566 * t816;
t757 = -qJDD(1) * t720 + t807;
t755 = -t783 - t870 + t907;
t745 = qJD(2) ^ 2;
t753 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t745 + t780;
t746 = qJD(1) ^ 2;
t752 = pkin(1) * t746 - pkin(7) * qJDD(1) + t781;
t500 = qJD(4) * t566 - t580 * t877 + t738 * t581;
t748 = t716 * t781 + t708 - t749;
t721 = -pkin(4) - t866;
t661 = t716 * t834 + t837;
t660 = -t716 * t836 + t835;
t659 = -t716 * t835 + t836;
t658 = t716 * t837 + t834;
t645 = t662 - t866;
t572 = t686 * t627;
t558 = t742 * t560;
t522 = pkin(5) * t846 - t766;
t512 = t742 * t513;
t497 = -pkin(10) * t846 + t826;
t493 = -pkin(5) * t764 - pkin(10) * t845 - t566 * t737 + t558;
t483 = t575 * t838 - t736 * t798 - t814 * t846 + (t845 * t881 + t858) * t741;
t482 = -t571 * t881 - t686 * t575;
t480 = pkin(5) * t763 + t500;
t465 = -pkin(10) * t763 + t759;
t464 = -pkin(10) * t857 + pkin(5) * t576 - t499 * t737 + t512 + (-t559 + (pkin(10) * t627 - t560) * t737) * qJD(5);
t1 = [(-t500 * t729 + t576 * t633 - t591 * t764 - t616 * t640 - t647 * t750 + t716 * t780 + t728 * t766) * MDP(18) + (t535 * t764 + t575 * t616 + t576 * t765 + t627 * t750) * MDP(14) + (-t573 * t680 + t574 * t679 - t604 * t669 + t605 * t667 - t623 * t671 - t624 * t668 - t634 * t625 + t635 * t747 - t781) * MDP(11) + (-t474 * t571 + t475 * t572 - t482 * t549 + t483 * t770) * MDP(28) + (-t474 * t572 - t482 * t770) * MDP(27) + (qJDD(2) * t739 + t743 * t745) * MDP(6) + (qJDD(2) * t743 - t739 * t745) * MDP(7) + (t535 * t627 - t575 * t765) * MDP(13) + (-t499 * t729 + t535 * t647 - t566 * t728 + t575 * t633 + t591 * t627 - t640 * t765 - t715 * t780) * MDP(19) + (-t576 * t729 + t728 * t764) * MDP(16) + (t475 * t764 - t483 * t603 - t531 * t571 - t549 * t576) * MDP(30) + (-t531 * t764 + t576 * t603) * MDP(31) + ((t464 * t741 - t465 * t736) * t603 + (t493 * t741 - t497 * t736) * t531 - t791 * t764 + t466 * t576 + t480 * t549 + t522 * t475 + t470 * t571 + t514 * t483 - g(1) * t642 - g(2) * t644 + ((-t493 * t736 - t497 * t741) * t603 + t467 * t764) * qJD(6)) * MDP(32) + (-t474 * t764 + t482 * t603 - t531 * t572 - t576 * t770) * MDP(29) + (-g(1) * t641 - g(2) * t643 - t467 * t576 - t470 * t572 + t522 * t474 - t480 * t770 + t514 * t482 - t487 * t764 + (-(-qJD(6) * t497 + t464) * t603 - t493 * t531 + t462 * t764) * t736 + (-(qJD(6) * t493 + t465) * t603 - t497 * t531 + t785 * t764) * t741) * MDP(33) + (-t534 * t764 - t576 * t811) * MDP(24) + (t510 * t764 - t526 * t627 - t576 * t597 + t763 * t811) * MDP(23) + (-t509 * t764 + t527 * t627 + t576 * t599 - t762 * t811) * MDP(22) + (-g(1) * t658 - g(2) * t660 + t479 * t845 - t496 * t576 + t500 * t599 - t509 * t766 + t528 * t762 - t534 * t826 + t759 * t811 + t760 * t764) * MDP(26) + (-(-t566 * t815 + t512) * t811 + t558 * t534 - t777 * t764 + t495 * t576 + t500 * t597 - t766 * t510 + t627 * t524 - g(1) * t659 - g(2) * t661 + (-(-qJD(5) * t560 - t499) * t811 - t566 * t534 - t786 * t764 + t479 * t627 + t528 * t575) * t737) * MDP(25) + (t574 * t635 + t624 * t605 + t573 * t634 + t623 * t604 - t757 * t720 + t695 * t723 - g(1) * (-t720 * t740 + t744 * t864) - g(2) * (t720 * t744 + t740 * t864)) * MDP(12) + (t575 * t729 + t627 * t728) * MDP(15) + ((-t597 * t742 - t853) * t575 + (-t507 - t510 * t742 + (t597 * t737 - t599 * t742) * qJD(5)) * t627) * MDP(21) + (t509 * t845 + t599 * t762) * MDP(20) + 0.2e1 * (t739 * t808 - t810 * t821) * MDP(5) + (qJDD(1) * t730 + 0.2e1 * t739 * t795) * MDP(4) + t780 * MDP(2) + t781 * MDP(3) + (t739 * t767 + t743 * t753) * MDP(9) + (-t739 * t753 + t743 * t767) * MDP(10) + qJDD(1) * MDP(1); ((t636 * t741 - t637 * t736) * t531 + t645 * t475 + (t736 * t774 - t741 * t776) * t603 + t825 * t549 + t901) * MDP(32) + (-t859 + t899) * MDP(29) + (-(t636 * t736 + t637 * t741) * t531 + t645 * t474 + (t736 * t776 + t741 * t774) * t603 - t825 * t770 + t902) * MDP(33) + (t855 + t898) * MDP(22) + (t811 * t853 + t900) * MDP(21) + (t662 * t510 - t887 * t742 + t771 * t737 + t823 * t597 - (-t663 * t815 + t905) * t811 + t801) * MDP(25) + ((t624 + t629) * t669 + (-t630 + t623) * t667 + (-t734 * t625 + ((-t795 - t809) * t733 + (-t796 + t808) * t734) * t733) * pkin(2)) * MDP(11) + (t750 - t852) * MDP(16) + (t769 - t856) * MDP(23) + (t779 - t860) * MDP(30) + (-MDP(4) * t739 * t743 + MDP(5) * t821) * t746 + (t662 * t509 + t771 * t742 - t737 * t907 + t823 * t599 - (t663 * t816 + t882) * t811 + t804) * MDP(26) - t879 * t765 + (t616 * t639 + t728 * t778 - t729 * t823 + t755) * MDP(18) + (t639 * t765 - t728 * t822 + t729 * t824 + t748) * MDP(19) + t880 * t616 + t897 + (-t623 * t629 - t624 * t630 + (-t867 + t573 * t734 + t574 * t733 + (-qJD(1) * t695 + t781) * t739) * pkin(2)) * MDP(12) + (t739 * t752 - t867) * MDP(9) + (t535 - t849) * MDP(15) + (g(3) * t739 + t743 * t752) * MDP(10) + MDP(7) * t808 + MDP(6) * t809 + qJDD(2) * MDP(8); (-t667 ^ 2 - t669 ^ 2) * MDP(11) + (t623 * t669 - t624 * t667 + t757 - t780) * MDP(12) + (-t750 - t852) * MDP(18) + (t535 + t849) * MDP(19) + (t769 + t856) * MDP(25) + (-t742 * t811 ^ 2 - t526 + t855) * MDP(26) + (t779 + t860) * MDP(32) + (-t859 - t899) * MDP(33); t535 * MDP(15) + t750 * MDP(16) + (t540 * t729 + t755) * MDP(18) + (t539 * t729 + t748) * MDP(19) + (t599 * t787 + t900) * MDP(21) + t898 * MDP(22) + (-t787 * t811 + t527) * MDP(23) + (-pkin(4) * t510 - t540 * t597 + t568 * t811 + (-t539 * t811 + t772) * t737 + (-t887 + t806) * t742 + t801) * MDP(25) + (-pkin(4) * t509 - t828 * t811 - t540 * t599 + t772 * t742 + (-t907 - t806) * t737 + t804) * MDP(26) + t899 * MDP(29) + t779 * MDP(30) + ((t704 * t741 - t705 * t736) * t531 + t721 * t475 + (t736 * t773 - t741 * t775) * t603 + t782 * t549 + t901) * MDP(32) + (-(t704 * t736 + t705 * t741) * t531 + t721 * t474 + (t736 * t775 + t741 * t773) * t603 - t782 * t770 + t902) * MDP(33) - (t729 * MDP(15) - t880) * t616 - (MDP(16) * t729 - MDP(22) * t599 + MDP(23) * t597 + MDP(29) * t770 + MDP(30) * t549 + t879) * t765 + t897; t599 * t597 * MDP(20) + (-t597 ^ 2 + t599 ^ 2) * MDP(21) + (-t597 * t811 + t509) * MDP(22) + (-t789 + (-qJD(5) - t811) * t599) * MDP(23) + t534 * MDP(24) + (-g(1) * t660 + g(2) * t658 - t496 * t811 - t528 * t599 + (t786 + t708) * t737 + t777) * MDP(25) + (g(1) * t661 - g(2) * t659 - t495 * t811 + t528 * t597 + t708 * t742 - t760) * MDP(26) + (t474 + t894) * MDP(29) + (-t475 - t893) * MDP(30) + (-(-t490 * t736 - t861) * t603 - t467 * qJD(6) + (t531 * t741 - t549 * t599 - t603 * t814) * pkin(5) + t884) * MDP(32) + ((-t491 * t603 - t462) * t736 + (t490 * t603 - t785) * t741 + (-t531 * t736 + t599 * t770 - t603 * t813) * pkin(5) + t885) * MDP(33) + t883; (t803 + t894) * MDP(29) + (-t790 - t893) * MDP(30) + (t467 * t603 + t884) * MDP(32) + (-t736 * t462 - t741 * t463 + t466 * t603 + t885) * MDP(33) + (-MDP(29) * t854 + MDP(30) * t770 - MDP(32) * t467 - MDP(33) * t862) * qJD(6) + t883;];
tau  = t1;

% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:47:09
% EndTime: 2019-03-09 20:47:24
% DurationCPUTime: 12.20s
% Computational Cost: add. (14108->659), mult. (31189->799), div. (0->0), fcn. (22516->14), ass. (0->292)
t742 = qJD(2) + qJD(3);
t751 = sin(qJ(3));
t755 = cos(qJ(2));
t894 = cos(qJ(3));
t824 = qJD(1) * t894;
t752 = sin(qJ(2));
t844 = qJD(1) * t752;
t914 = -t751 * t844 + t755 * t824;
t915 = t742 * t914;
t862 = t751 * t755;
t680 = t752 * t894 + t862;
t620 = t742 * t680;
t816 = qJDD(1) * t894;
t838 = qJDD(1) * t752;
t579 = qJD(1) * t620 + t751 * t838 - t755 * t816;
t577 = qJDD(4) + t579;
t747 = sin(pkin(10));
t748 = cos(pkin(10));
t750 = sin(qJ(4));
t754 = cos(qJ(4));
t901 = -t747 * t750 + t748 * t754;
t598 = t901 * t680;
t677 = t747 * t754 + t748 * t750;
t658 = t677 * qJD(4);
t911 = -t677 * t914 + t658;
t841 = qJD(4) * t754;
t842 = qJD(4) * t750;
t910 = -t747 * t842 + t748 * t841 - t901 * t914;
t662 = -qJD(1) * t862 - t752 * t824;
t629 = -t662 * t750 - t754 * t742;
t795 = t662 * t754 - t742 * t750;
t568 = t748 * t629 - t747 * t795;
t653 = qJD(4) - t914;
t913 = t568 * t653;
t741 = qJDD(2) + qJDD(3);
t839 = qJD(1) * qJD(2);
t820 = t755 * t839;
t895 = pkin(8) + pkin(7);
t624 = qJDD(2) * pkin(2) + t895 * (-t820 - t838);
t821 = t752 * t839;
t837 = qJDD(1) * t755;
t628 = t895 * (-t821 + t837);
t694 = t895 * t752;
t682 = qJD(1) * t694;
t886 = qJD(2) * pkin(2);
t667 = -t682 + t886;
t695 = t895 * t755;
t684 = qJD(1) * t695;
t823 = qJD(3) * t894;
t843 = qJD(3) * t751;
t808 = -t894 * t624 + t751 * t628 + t667 * t843 + t684 * t823;
t542 = -pkin(3) * t741 + t808;
t746 = qJ(2) + qJ(3);
t739 = cos(t746);
t726 = g(3) * t739;
t912 = t542 + t726;
t876 = t914 * t750;
t905 = (t842 - t876) * pkin(4);
t738 = sin(t746);
t756 = cos(qJ(1));
t870 = t738 * t756;
t753 = sin(qJ(1));
t871 = t738 * t753;
t909 = g(1) * t870 + g(2) * t871;
t783 = -t751 * t752 + t755 * t894;
t619 = t742 * t783;
t825 = t680 * t841;
t908 = t619 * t750 + t825;
t796 = -t629 * t747 - t748 * t795;
t907 = t796 ^ 2;
t608 = -pkin(3) * t662 - pkin(9) * t914;
t592 = pkin(2) * t844 + t608;
t588 = t754 * t592;
t665 = t751 * t684;
t617 = -t682 * t894 - t665;
t740 = t754 * qJ(5);
t803 = -t662 * pkin(4) - t740 * t914;
t535 = -t617 * t750 + t588 + t803;
t833 = qJ(5) * t876;
t847 = t750 * t592 + t754 * t617;
t547 = -t833 + t847;
t737 = t754 * qJD(5);
t809 = pkin(2) * t823;
t799 = t754 * t809;
t729 = pkin(2) * t751 + pkin(9);
t857 = -qJ(5) - t729;
t812 = qJD(4) * t857;
t618 = t750 * t812 + t737 + t799;
t763 = (-t809 - qJD(5)) * t750 + t754 * t812;
t856 = (-t535 + t763) * t748 + (t547 - t618) * t747;
t906 = pkin(5) * t911 - qJ(6) * t910 - qJD(6) * t677 + t905;
t600 = t754 * t608;
t612 = t667 * t894 - t665;
t538 = -t612 * t750 + t600 + t803;
t848 = t750 * t608 + t754 * t612;
t549 = -t833 + t848;
t749 = -qJ(5) - pkin(9);
t817 = qJD(4) * t749;
t651 = t750 * t817 + t737;
t777 = -t750 * qJD(5) + t754 * t817;
t851 = (-t538 + t777) * t748 + (t549 - t651) * t747;
t666 = t894 * t684;
t616 = -t751 * t682 + t666;
t806 = pkin(2) * t843 - t616;
t743 = qJ(4) + pkin(10);
t736 = cos(t743);
t872 = t736 * t739;
t735 = sin(t743);
t873 = t735 * t739;
t904 = pkin(5) * t872 + qJ(6) * t873;
t902 = -t894 * t694 - t751 * t695;
t805 = g(1) * t756 + g(2) * t753;
t900 = t726 - t909;
t578 = t751 * t837 + t752 * t816 + t915;
t551 = t754 * t578 + t662 * t842 + t750 * t741 + t742 * t841;
t552 = -qJD(4) * t795 + t578 * t750 - t754 * t741;
t514 = t551 * t747 + t748 * t552;
t515 = t551 * t748 - t552 * t747;
t591 = t748 * t651 + t747 * t777;
t692 = pkin(9) * t754 + t740;
t822 = t749 * t750;
t626 = t692 * t747 - t748 * t822;
t627 = t748 * t692 + t747 * t822;
t899 = -t627 * t514 + t515 * t626 - t591 * t568;
t564 = t748 * t618 + t747 * t763;
t668 = t729 * t754 + t740;
t815 = t857 * t750;
t604 = t668 * t747 - t748 * t815;
t605 = t748 * t668 + t747 * t815;
t898 = -t605 * t514 + t515 * t604 - t564 * t568;
t594 = -t742 * pkin(3) - t612;
t562 = t629 * pkin(4) + qJD(5) + t594;
t511 = t568 * pkin(5) - qJ(6) * t796 + t562;
t868 = t739 * t753;
t634 = t735 * t868 + t736 * t756;
t858 = t756 * t735;
t861 = t753 * t736;
t636 = t739 * t858 - t861;
t725 = g(3) * t738;
t892 = pkin(2) * t755;
t732 = pkin(1) + t892;
t655 = pkin(2) * t821 - qJDD(1) * t732;
t530 = pkin(3) * t579 - pkin(9) * t578 + t655;
t526 = t754 * t530;
t769 = t751 * t624 + t628 * t894 + t667 * t823 - t684 * t843;
t541 = t741 * pkin(9) + t769;
t693 = t732 * qJD(1);
t589 = -pkin(3) * t914 + pkin(9) * t662 - t693;
t613 = t751 * t667 + t666;
t595 = pkin(9) * t742 + t613;
t554 = t589 * t750 + t595 * t754;
t476 = pkin(4) * t577 - qJ(5) * t551 - qJD(4) * t554 + qJD(5) * t795 - t541 * t750 + t526;
t779 = t750 * t530 + t754 * t541 + t589 * t841 - t595 * t842;
t479 = -qJ(5) * t552 - qJD(5) * t629 + t779;
t465 = t748 * t476 - t747 * t479;
t819 = -qJDD(6) + t465;
t897 = g(1) * t636 + g(2) * t634 - t511 * t796 + t725 * t735 + t819;
t896 = t653 ^ 2;
t893 = pkin(2) * t752;
t891 = pkin(5) * t577;
t890 = pkin(5) * t662;
t887 = t754 * pkin(4);
t533 = -qJ(5) * t629 + t554;
t528 = t748 * t533;
t553 = t754 * t589 - t595 * t750;
t532 = qJ(5) * t795 + t553;
t494 = t532 * t747 + t528;
t885 = t494 * t796;
t882 = t533 * t747;
t881 = t551 * t750;
t880 = t594 * t914;
t878 = t629 * t653;
t877 = t795 * t653;
t875 = t680 * t750;
t874 = t680 * t754;
t730 = pkin(3) + t887;
t690 = t739 * t730;
t869 = t739 * t749;
t867 = t739 * t756;
t864 = t750 * t753;
t863 = t750 * t756;
t860 = t753 * t754;
t633 = -t751 * t694 + t695 * t894;
t625 = t754 * t633;
t859 = t754 * t756;
t466 = t747 * t476 + t748 * t479;
t835 = t752 * t886;
t561 = pkin(3) * t620 - pkin(9) * t619 + t835;
t558 = t754 * t561;
t828 = qJD(2) * t895;
t683 = t752 * t828;
t685 = t755 * t828;
t572 = qJD(3) * t902 - t894 * t683 - t751 * t685;
t611 = -pkin(3) * t783 - pkin(9) * t680 - t732;
t793 = -qJ(5) * t619 - qJD(5) * t680;
t483 = pkin(4) * t620 - t572 * t750 + t558 + t793 * t754 + (-t625 + (qJ(5) * t680 - t611) * t750) * qJD(4);
t830 = t750 * t561 + t754 * t572 + t611 * t841;
t491 = -qJ(5) * t825 + (-qJD(4) * t633 + t793) * t750 + t830;
t473 = t747 * t483 + t748 * t491;
t523 = pkin(4) * t653 + t532;
t493 = t747 * t523 + t528;
t502 = t747 * t535 + t748 * t547;
t504 = t747 * t538 + t748 * t549;
t602 = t754 * t611;
t548 = -pkin(4) * t783 - t633 * t750 - t680 * t740 + t602;
t846 = t750 * t611 + t625;
t556 = -qJ(5) * t875 + t846;
t513 = t747 * t548 + t748 * t556;
t855 = -t890 - t856;
t646 = t662 * qJ(6);
t496 = -t646 + t502;
t854 = t564 - t496;
t853 = t806 + t906;
t852 = -t613 + t906;
t850 = -t890 - t851;
t498 = -t646 + t504;
t849 = t591 - t498;
t744 = t752 ^ 2;
t845 = -t755 ^ 2 + t744;
t495 = t532 * t748 - t882;
t840 = qJD(6) - t495;
t836 = t894 * pkin(2);
t834 = qJD(4) * pkin(9) * t653;
t832 = t739 * t863;
t831 = t577 * qJ(6) + t466;
t829 = g(1) * t867 + g(2) * t868 + t725;
t586 = t594 * t841;
t813 = -t738 * t749 + t690;
t811 = t653 * t754;
t810 = -qJD(4) * t589 - t541;
t731 = -t836 - pkin(3);
t807 = -g(1) * t871 + g(2) * t870;
t804 = g(1) * t753 - g(2) * t756;
t802 = -t595 * t841 + t526;
t801 = -pkin(9) * t577 - t880;
t472 = t483 * t748 - t491 * t747;
t492 = t523 * t748 - t882;
t512 = t548 * t748 - t556 * t747;
t797 = -t577 * t729 - t880;
t794 = t730 * t738 + t869;
t792 = -t554 * t662 + t750 * t912 + t586;
t791 = t553 * t662 + t594 * t842 + t754 * t909;
t790 = pkin(4) * t875 - t902;
t789 = t813 + t892;
t787 = pkin(5) * t736 + qJ(6) * t735 + t730;
t786 = t805 * t738;
t785 = t805 * t739;
t784 = -0.2e1 * pkin(1) * t839 - pkin(7) * qJDD(2);
t647 = t739 * t864 + t859;
t782 = t619 * t754 - t680 * t842;
t573 = -t751 * t683 + t685 * t894 - t694 * t843 + t695 * t823;
t607 = -pkin(5) * t901 - t677 * qJ(6) - t730;
t776 = pkin(4) * t908 + t573;
t775 = pkin(4) * t864 + t730 * t867 + t756 * t732 - t749 * t870 + t753 * t895;
t758 = qJD(2) ^ 2;
t774 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t758 + t804;
t759 = qJD(1) ^ 2;
t773 = pkin(1) * t759 - pkin(7) * qJDD(1) + t805;
t772 = -t693 * t662 - t808 - t900;
t510 = pkin(4) * t552 + qJDD(5) + t542;
t771 = t756 * t895 + t749 * t871 + pkin(4) * t863 + (-t732 - t690) * t753;
t460 = qJD(6) * t653 + t831;
t462 = -t819 - t891;
t486 = -pkin(5) * t653 + qJD(6) - t492;
t487 = qJ(6) * t653 + t493;
t768 = t460 * t901 + t462 * t677 + t486 * t910 - t487 * t911 - t829;
t767 = -t465 * t677 + t466 * t901 - t492 * t910 - t493 * t911 - t829;
t469 = pkin(5) * t514 - qJ(6) * t515 - qJD(6) * t796 + t510;
t766 = -g(3) * t872 - t469 * t901 - t486 * t662 + t511 * t911 + t736 * t909;
t765 = -g(3) * t873 - t469 * t677 + t487 * t662 - t511 * t910 + t735 * t909;
t762 = ((t551 - t878) * t754 + (-t552 + t877) * t750) * MDP(19) + (-t795 * t811 + t881) * MDP(18) + (t577 * t754 - t629 * t662 - t750 * t896) * MDP(21) + (t577 * t750 + t653 * t811 - t662 * t795) * MDP(20) + (t578 - t915) * MDP(13) + (-t662 * t742 - t579) * MDP(14) + (t662 ^ 2 - t914 ^ 2) * MDP(12) + t741 * MDP(15) + (MDP(11) * t914 + MDP(22) * t653) * t662;
t761 = t693 * t914 - t769 + t829;
t724 = -pkin(4) * t748 - pkin(5);
t722 = pkin(4) * t747 + qJ(6);
t719 = pkin(4) * t860;
t650 = t739 * t859 + t864;
t649 = -t832 + t860;
t648 = -t739 * t860 + t863;
t637 = t735 * t753 + t736 * t867;
t635 = t739 * t861 - t858;
t597 = t677 * t680;
t593 = -t836 + t607;
t537 = -t619 * t901 + t658 * t680;
t536 = -qJD(4) * t598 - t619 * t677;
t527 = t597 * pkin(5) - t598 * qJ(6) + t790;
t519 = -pkin(4) * t795 + pkin(5) * t796 + qJ(6) * t568;
t509 = pkin(5) * t783 - t512;
t508 = -qJ(6) * t783 + t513;
t480 = -t536 * pkin(5) + t537 * qJ(6) - t598 * qJD(6) + t776;
t471 = -pkin(5) * t620 - t472;
t470 = qJ(6) * t620 - qJD(6) * t783 + t473;
t1 = [(t578 * t783 - t579 * t680 + t619 * t914 + t620 * t662) * MDP(12) + (-t573 * t742 - t579 * t732 - t620 * t693 - t655 * t783 + t739 * t804 + t741 * t902 - t835 * t914) * MDP(16) + (-t465 * t598 - t466 * t597 - t472 * t796 - t473 * t568 + t492 * t537 + t493 * t536 - t512 * t515 - t513 * t514 - t807) * MDP(25) + (-t460 * t597 + t462 * t598 - t470 * t568 + t471 * t796 - t486 * t537 + t487 * t536 - t508 * t514 + t509 * t515 - t807) * MDP(28) + (t752 * t784 + t755 * t774) * MDP(9) + (-t752 * t774 + t755 * t784) * MDP(10) + t804 * MDP(2) + (t460 * t508 + t487 * t470 + t469 * t527 + t511 * t480 + t462 * t509 + t486 * t471 - g(1) * (-pkin(5) * t635 - qJ(6) * t634 + t771) - g(2) * (pkin(5) * t637 + qJ(6) * t636 + t775)) * MDP(30) + (-t572 * t742 - t578 * t732 - t619 * t693 - t633 * t741 + t655 * t680 - t662 * t835 + t807) * MDP(17) + (t552 * t783 - t577 * t875 - t620 * t629 - t653 * t908) * MDP(21) + qJDD(1) * MDP(1) + ((-t629 * t754 + t750 * t795) * t619 + (-t881 - t552 * t754 + (t629 * t750 + t754 * t795) * qJD(4)) * t680) * MDP(19) + (t551 * t874 - t782 * t795) * MDP(18) + (g(1) * t634 - g(2) * t636 - t460 * t783 - t469 * t598 + t470 * t653 - t480 * t796 + t487 * t620 + t508 * t577 + t511 * t537 - t515 * t527) * MDP(29) + (-t620 * t742 + t741 * t783) * MDP(14) + (g(1) * t635 - g(2) * t637 + t462 * t783 + t469 * t597 - t471 * t653 + t480 * t568 - t486 * t620 - t509 * t577 - t511 * t536 + t514 * t527) * MDP(27) + (-t577 * t783 + t620 * t653) * MDP(22) + (-t551 * t783 + t577 * t874 - t620 * t795 + t653 * t782) * MDP(20) + (qJDD(2) * t752 + t755 * t758) * MDP(6) + (qJDD(2) * t755 - t752 * t758) * MDP(7) + (-g(1) * t771 - g(2) * t775 + t465 * t512 + t466 * t513 + t492 * t472 + t493 * t473 + t510 * t790 + t562 * t776) * MDP(26) + (t619 * t742 + t680 * t741) * MDP(13) + (t578 * t680 - t619 * t662) * MDP(11) + (-(-t633 * t842 + t830) * t653 - t846 * t577 + t779 * t783 - t554 * t620 - t573 * t795 - t902 * t551 + t542 * t874 - g(1) * t647 - g(2) * t649 + t782 * t594) * MDP(24) + ((-t633 * t841 + t558) * t653 + t602 * t577 - t802 * t783 + t553 * t620 + t573 * t629 - t902 * t552 + t680 * t586 - g(1) * t648 - g(2) * t650 + ((-qJD(4) * t611 - t572) * t653 - t633 * t577 - t810 * t783 + t542 * t680 + t594 * t619) * t750) * MDP(23) + t805 * MDP(3) + (qJDD(1) * t744 + 0.2e1 * t752 * t820) * MDP(4) + 0.2e1 * (t752 * t837 - t839 * t845) * MDP(5); (t731 * t551 + t797 * t754 - t750 * t786 - t806 * t795 + (t729 * t842 - t799 + t847) * t653 + t792) * MDP(24) + t762 + (t616 * t742 + (t741 * t894 - t742 * t843 + t844 * t914) * pkin(2) + t772) * MDP(16) + (t514 * t593 + t568 * t853 - t577 * t604 - t653 * t855 + t766) * MDP(27) + (-t515 * t593 + t577 * t605 + t653 * t854 - t796 * t853 + t765) * MDP(29) + (t496 * t568 + t796 * t855 + t768 + t898) * MDP(28) + (t502 * t568 - t796 * t856 + t767 + t898) * MDP(25) + (t466 * t605 - t465 * t604 + t510 * (t731 - t887) - g(3) * t789 + (t905 + t806) * t562 + (t564 - t502) * t493 + t856 * t492 + t805 * (t794 + t893)) * MDP(26) + (-g(3) * t755 + t752 * t773) * MDP(9) + (g(3) * t752 + t755 * t773) * MDP(10) + MDP(7) * t837 + (t731 * t552 - t912 * t754 + t797 * t750 + t806 * t629 + (-t729 * t841 - t588 + (-t809 + t617) * t750) * t653 + t791) * MDP(23) + (t460 * t605 + t469 * t593 + t462 * t604 - g(3) * (t789 + t904) + t853 * t511 + t854 * t487 + t855 * t486 + t805 * (t738 * t787 + t869 + t893)) * MDP(30) + (t617 * t742 + (t662 * t844 - t741 * t751 - t742 * t823) * pkin(2) + t761) * MDP(17) + qJDD(2) * MDP(8) + MDP(6) * t838 + (-MDP(4) * t752 * t755 + MDP(5) * t845) * t759; t762 + (t460 * t627 + t469 * t607 + t462 * t626 - g(3) * (t690 + t904) + t749 * t785 + t852 * t511 + t849 * t487 + t850 * t486 + (g(3) * t749 + t787 * t805) * t738) * MDP(30) + (-pkin(3) * t551 + t848 * t653 + t613 * t795 + t801 * t754 + (-t786 + t834) * t750 + t792) * MDP(24) + (-t515 * t607 + t577 * t627 + t653 * t849 - t796 * t852 + t765) * MDP(29) + (t504 * t568 - t796 * t851 + t767 + t899) * MDP(25) + (t498 * t568 + t796 * t850 + t768 + t899) * MDP(28) + (t466 * t627 - t465 * t626 - t510 * t730 - g(3) * t813 + (-t613 + t905) * t562 + (t591 - t504) * t493 + t851 * t492 + t805 * t794) * MDP(26) + (-pkin(3) * t552 - t600 * t653 - t613 * t629 + (t612 * t653 + t801) * t750 + (-t912 - t834) * t754 + t791) * MDP(23) + (t514 * t607 + t568 * t852 - t577 * t626 - t653 * t850 + t766) * MDP(27) + (t612 * t742 + t761) * MDP(17) + (t613 * t742 + t772) * MDP(16); -t795 * t629 * MDP(18) + (-t629 ^ 2 + t795 ^ 2) * MDP(19) + (t551 + t878) * MDP(20) + (-t552 - t877) * MDP(21) + t577 * MDP(22) + (-g(1) * t649 + g(2) * t647 + t554 * t653 + t594 * t795 + (t810 + t725) * t750 + t802) * MDP(23) + (g(1) * t650 - g(2) * t648 + t553 * t653 + t594 * t629 + t725 * t754 - t779) * MDP(24) + (t493 * t796 - t885 + (-t514 * t747 - t515 * t748) * pkin(4) + (-t492 + t495) * t568) * MDP(25) + (-g(1) * t719 + t492 * t494 - t493 * t495 + (g(2) * t859 + t465 * t748 + t466 * t747 + t562 * t795 + (t785 + t725) * t750) * pkin(4)) * MDP(26) + (t494 * t653 - t519 * t568 + (pkin(5) - t724) * t577 + t897) * MDP(27) + (t487 * t796 - t514 * t722 + t515 * t724 - t885 + (t486 - t840) * t568) * MDP(28) + (-t736 * t725 - g(1) * t637 - g(2) * t635 - t511 * t568 + t519 * t796 + t577 * t722 + (0.2e1 * qJD(6) - t495) * t653 + t831) * MDP(29) + (t460 * t722 + t462 * t724 - t511 * t519 - t486 * t494 - g(1) * (-pkin(4) * t832 - pkin(5) * t636 + qJ(6) * t637 + t719) - g(2) * (-pkin(4) * t647 - pkin(5) * t634 + qJ(6) * t635) + t840 * t487 - (-pkin(4) * t750 - pkin(5) * t735 + qJ(6) * t736) * t725) * MDP(30); (t492 * t796 + t493 * t568 + t510 + t900) * MDP(26) + (t653 * t796 + t514) * MDP(27) + (-t515 + t913) * MDP(29) + (-t486 * t796 + t487 * t568 + t469 + t900) * MDP(30) + (MDP(25) + MDP(28)) * (-t568 ^ 2 - t907); (t515 + t913) * MDP(28) + (-t896 - t907) * MDP(29) + (-t487 * t653 - t891 - t897) * MDP(30) + (t568 * t796 - t577) * MDP(27);];
tau  = t1;

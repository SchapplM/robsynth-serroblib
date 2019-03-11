% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:23:31
% EndTime: 2019-03-10 01:23:54
% DurationCPUTime: 16.81s
% Computational Cost: add. (13093->690), mult. (29287->888), div. (0->0), fcn. (21375->14), ass. (0->296)
t769 = cos(qJ(3));
t766 = sin(qJ(2));
t770 = cos(qJ(2));
t897 = t769 * t770;
t812 = pkin(3) * t766 - pkin(9) * t897;
t924 = pkin(8) + pkin(9);
t850 = qJD(3) * t924;
t823 = pkin(2) * t766 - pkin(8) * t770;
t709 = t823 * qJD(1);
t765 = sin(qJ(3));
t870 = qJD(1) * t766;
t848 = t765 * t870;
t875 = pkin(7) * t848 + t769 * t709;
t966 = qJD(1) * t812 + t769 * t850 + t875;
t686 = t765 * t709;
t899 = t766 * t769;
t900 = t765 * t770;
t965 = -t686 - (-pkin(7) * t899 - pkin(9) * t900) * qJD(1) - t765 * t850;
t855 = qJD(1) * qJD(2);
t840 = t770 * t855;
t854 = qJDD(1) * t766;
t856 = t769 * qJD(2);
t862 = qJD(3) * t766;
t949 = -qJD(1) * t862 + qJDD(2);
t612 = qJD(3) * t856 + (t840 + t854) * t769 + t949 * t765;
t869 = qJD(1) * t770;
t613 = t765 * (qJD(2) * (qJD(3) + t869) + t854) - t949 * t769;
t702 = t848 - t856;
t866 = qJD(2) * t765;
t704 = t769 * t870 + t866;
t764 = sin(qJ(4));
t768 = cos(qJ(4));
t858 = qJD(4) * t768;
t860 = qJD(4) * t764;
t542 = t768 * t612 - t764 * t613 - t702 * t858 - t704 * t860;
t763 = sin(qJ(5));
t811 = t764 * t612 + t768 * t613 - t702 * t860 + t704 * t858;
t815 = -t702 * t764 + t768 * t704;
t816 = -t702 * t768 - t764 * t704;
t923 = cos(qJ(5));
t842 = qJD(5) * t923;
t857 = qJD(5) * t763;
t501 = -t923 * t542 + t763 * t811 + t815 * t857 - t816 * t842;
t928 = t763 * t816 + t815 * t923;
t502 = qJD(5) * t928 + t763 * t542 + t923 * t811;
t567 = t763 * t815 - t816 * t923;
t564 = t567 ^ 2;
t752 = t770 * qJDD(1);
t935 = -t766 * t855 + t752;
t699 = qJDD(3) - t935;
t692 = qJDD(4) + t699;
t677 = qJDD(5) + t692;
t735 = -qJD(3) + t869;
t726 = -qJD(4) + t735;
t717 = -qJD(5) + t726;
t925 = t928 ^ 2;
t964 = t677 * MDP(29) + (-t717 * t928 - t502) * MDP(28) + t567 * MDP(25) * t928 + (-t567 * t717 - t501) * MDP(27) + (-t564 + t925) * MDP(26);
t813 = t764 * t765 - t768 * t769;
t934 = qJD(4) + qJD(3);
t641 = t934 * t813;
t802 = t813 * t770;
t960 = -qJD(1) * t802 + t641;
t963 = t692 * MDP(22) + (t815 ^ 2 - t816 ^ 2) * MDP(19) + (t726 * t816 + t542) * MDP(20) + (-t726 * t815 - t811) * MDP(21) - t816 * MDP(18) * t815 + t964;
t962 = t965 * t764 + t768 * t966;
t705 = t764 * t769 + t765 * t768;
t653 = t705 * t869;
t797 = t705 * qJD(3);
t780 = -qJD(4) * t705 - t797;
t956 = t653 + t780;
t720 = t924 * t765;
t959 = -t720 * t858 - t764 * t966 + t965 * t768;
t951 = t567 * qJ(6);
t721 = t924 * t769;
t876 = -t764 * t720 + t768 * t721;
t958 = pkin(4) * t870 - pkin(10) * t960 + qJD(4) * t876 + t962;
t902 = t764 * t721;
t920 = pkin(10) * t705;
t957 = (-t902 - t920) * qJD(4) + t959 + (t653 - t797) * pkin(10);
t745 = pkin(7) * t854;
t675 = -qJDD(2) * pkin(2) + pkin(7) * t840 + t745;
t767 = sin(qJ(1));
t771 = cos(qJ(1));
t822 = g(1) * t771 + g(2) * t767;
t914 = g(3) * t770;
t791 = t766 * t822 - t914;
t955 = qJD(3) * pkin(8) * t735 - t675 + t791;
t714 = -pkin(2) * t770 - pkin(8) * t766 - pkin(1);
t693 = t714 * qJD(1);
t747 = pkin(7) * t869;
t719 = qJD(2) * pkin(8) + t747;
t637 = t769 * t693 - t719 * t765;
t594 = -pkin(9) * t704 + t637;
t586 = -pkin(3) * t735 + t594;
t638 = t693 * t765 + t719 * t769;
t595 = -pkin(9) * t702 + t638;
t589 = t764 * t595;
t545 = t768 * t586 - t589;
t946 = pkin(10) * t815;
t526 = t545 - t946;
t519 = -pkin(4) * t726 + t526;
t591 = t768 * t595;
t546 = t586 * t764 + t591;
t945 = pkin(10) * t816;
t527 = t546 + t945;
t521 = t763 * t527;
t503 = t923 * t519 - t521;
t944 = qJ(6) * t928;
t495 = t503 - t944;
t494 = -pkin(5) * t717 + t495;
t523 = t923 * t527;
t504 = t763 * t519 + t523;
t496 = t504 - t951;
t954 = -t494 * t567 + t496 * t928;
t718 = -qJD(2) * pkin(2) + pkin(7) * t870;
t645 = pkin(3) * t702 + t718;
t583 = -pkin(4) * t816 + t645;
t762 = qJ(3) + qJ(4);
t757 = qJ(5) + t762;
t740 = sin(t757);
t741 = cos(t757);
t898 = t767 * t770;
t649 = t740 * t771 - t741 * t898;
t896 = t770 * t771;
t651 = t740 * t767 + t741 * t896;
t712 = t823 * qJD(2);
t642 = qJD(1) * t712 + qJDD(1) * t714;
t632 = t769 * t642;
t674 = pkin(7) * t935 + qJDD(2) * pkin(8);
t528 = pkin(3) * t699 - pkin(9) * t612 - qJD(3) * t638 - t674 * t765 + t632;
t861 = qJD(3) * t769;
t863 = qJD(3) * t765;
t801 = t765 * t642 + t769 * t674 + t693 * t861 - t719 * t863;
t536 = -pkin(9) * t613 + t801;
t785 = -qJD(4) * t546 + t768 * t528 - t764 * t536;
t488 = pkin(4) * t692 - pkin(10) * t542 + t785;
t825 = t764 * t528 + t768 * t536 + t586 * t858 - t595 * t860;
t490 = -pkin(10) * t811 + t825;
t826 = -t763 * t488 - t923 * t490 - t519 * t842 + t527 * t857;
t915 = g(3) * t766;
t948 = g(1) * t651 - g(2) * t649 + t567 * t583 + t741 * t915 + t826;
t794 = t923 * t813;
t889 = qJD(5) * t794 + t705 * t857 - t763 * t956 + t923 * t960;
t634 = t705 * t923 - t763 * t813;
t888 = qJD(5) * t634 - t763 * t960 - t923 * t956;
t922 = pkin(3) * t765;
t695 = t869 * t922 + t747;
t938 = pkin(3) * t863 - t695;
t864 = qJD(2) * t770;
t846 = t765 * t864;
t950 = t766 * t861 + t846;
t783 = -qJD(5) * t504 + t923 * t488 - t763 * t490;
t904 = t741 * t771;
t648 = t740 * t898 + t904;
t905 = t741 * t767;
t650 = -t740 * t896 + t905;
t933 = -g(1) * t650 + g(2) * t648 + t740 * t915;
t932 = -t583 * t928 + t783 + t933;
t532 = pkin(5) * t567 + qJD(6) + t583;
t943 = t532 * t928;
t667 = t813 * t766;
t701 = t769 * t714;
t921 = pkin(7) * t765;
t636 = -pkin(9) * t899 + t701 + (-pkin(3) - t921) * t770;
t737 = pkin(7) * t897;
t874 = t765 * t714 + t737;
t901 = t765 * t766;
t643 = -pkin(9) * t901 + t874;
t830 = t768 * t636 - t643 * t764;
t560 = -pkin(4) * t770 + pkin(10) * t667 + t830;
t803 = t705 * t766;
t881 = t764 * t636 + t768 * t643;
t562 = -pkin(10) * t803 + t881;
t887 = t763 * t560 + t923 * t562;
t829 = -t768 * t720 - t902;
t607 = t829 - t920;
t608 = -pkin(10) * t813 + t876;
t885 = t763 * t607 + t923 * t608;
t942 = -pkin(4) * t956 + t938;
t832 = -t594 * t764 - t591;
t533 = t832 - t945;
t886 = t768 * t594 - t589;
t534 = t886 - t946;
t743 = pkin(3) * t768 + pkin(4);
t903 = t763 * t764;
t941 = t743 * t842 + (-t764 * t857 + (t768 * t923 - t903) * qJD(4)) * pkin(3) - t763 * t533 - t923 * t534;
t849 = t923 * t764;
t940 = t743 * t857 - (-t764 * t842 + (-t763 * t768 - t849) * qJD(4)) * pkin(3) + t923 * t533 - t534 * t763;
t937 = qJD(5) * t885 + t763 * t957 + t923 * t958;
t798 = qJD(2) * t705;
t793 = t770 * t798;
t774 = t934 * t667 - t793;
t936 = t607 * t842 - t608 * t857 - t763 * t958 + t923 * t957;
t754 = sin(t762);
t755 = cos(t762);
t659 = t754 * t898 + t755 * t771;
t661 = -t754 * t896 + t755 * t767;
t931 = -g(1) * t661 + g(2) * t659 - t645 * t815 + t754 * t915 + t785;
t660 = t754 * t771 - t755 * t898;
t662 = t754 * t767 + t755 * t896;
t930 = g(1) * t662 - g(2) * t660 - t645 * t816 + t755 * t915 - t825;
t759 = t769 * pkin(3);
t744 = -pkin(2) - t759;
t696 = -pkin(4) * t754 - pkin(5) * t740;
t664 = -t696 + t922;
t913 = pkin(7) + t664;
t912 = t542 * t705;
t585 = t613 * pkin(3) + t675;
t911 = t585 * t705;
t910 = t612 * t765;
t909 = t692 * t705;
t908 = t702 * t735;
t907 = t704 * t735;
t906 = t704 * t769;
t633 = t705 * t763 + t794;
t895 = -qJ(6) * t888 - qJD(6) * t633 + t936;
t894 = -pkin(5) * t870 + qJ(6) * t889 - t634 * qJD(6) - t937;
t893 = -t495 + t494;
t892 = t923 * t526 - t521;
t884 = t941 + t944;
t883 = -t940 - t951;
t878 = t765 * t712 + t714 * t861;
t865 = qJD(2) * t766;
t877 = t769 * t712 + t865 * t921;
t697 = pkin(4) * t755 + pkin(5) * t741;
t713 = pkin(3) * t901 + t766 * pkin(7);
t760 = t766 ^ 2;
t873 = -t770 ^ 2 + t760;
t868 = qJD(2) * t702;
t867 = qJD(2) * t704;
t859 = qJD(4) * t765;
t575 = t812 * qJD(2) + (-t737 + (pkin(9) * t766 - t714) * t765) * qJD(3) + t877;
t580 = -t950 * pkin(9) + (-t766 * t856 - t770 * t863) * pkin(7) + t878;
t852 = t764 * t575 + t768 * t580 + t636 * t858;
t665 = t759 + t697;
t646 = pkin(3) * t950 + pkin(7) * t864;
t847 = t735 * t856;
t845 = t735 * t863;
t844 = t735 * t861;
t838 = -t526 * t763 - t523;
t833 = t923 * t560 - t562 * t763;
t831 = t923 * t607 - t608 * t763;
t828 = -qJD(3) * t693 - t674;
t824 = -pkin(3) * t903 + t923 * t743;
t588 = pkin(3) * t704 + pkin(4) * t815;
t821 = g(1) * t767 - g(2) * t771;
t820 = t719 * t861 - t632;
t819 = -pkin(8) * t699 + qJD(3) * t718;
t658 = pkin(2) + t665;
t756 = -qJ(6) - pkin(10) - t924;
t818 = t658 * t770 - t756 * t766;
t809 = pkin(1) + t818;
t808 = -0.2e1 * pkin(1) * t855 - pkin(7) * qJDD(2);
t806 = t699 * t765 - t844;
t805 = t699 * t769 + t845;
t804 = t645 * t705;
t581 = -qJD(2) * t802 + t766 * t780;
t786 = -qJD(4) * t881 + t768 * t575 - t580 * t764;
t507 = pkin(4) * t865 - pkin(10) * t581 + t786;
t795 = -t765 * t862 + t770 * t856;
t509 = (-t899 * t934 - t846) * pkin(10) * t768 + (-qJD(4) * t643 + (t766 * t859 - t795) * pkin(10)) * t764 + t852;
t800 = t763 * t507 + t923 * t509 + t560 * t842 - t562 * t857;
t773 = qJD(1) ^ 2;
t796 = pkin(1) * t773 + t822;
t655 = pkin(4) * t813 + t744;
t792 = t923 * t803;
t772 = qJD(2) ^ 2;
t789 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t772 + t821;
t639 = pkin(4) * t803 + t713;
t598 = -t667 * t923 - t763 * t803;
t782 = -qJD(5) * t887 + t923 * t507 - t763 * t509;
t516 = pkin(4) * t811 + t585;
t779 = -t768 * t861 - t769 * t858 + (t859 + t863) * t764;
t778 = t705 * t934 - t653;
t491 = t502 * pkin(5) + qJDD(6) + t516;
t563 = -pkin(4) * t774 + t646;
t742 = pkin(4) * t923 + pkin(5);
t684 = t765 * t767 + t769 * t896;
t683 = -t765 * t896 + t767 * t769;
t682 = t765 * t771 - t767 * t897;
t681 = t765 * t898 + t769 * t771;
t676 = pkin(3) * t849 + t763 * t743;
t668 = pkin(5) + t824;
t597 = -t667 * t763 + t792;
t538 = -qJ(6) * t633 + t885;
t537 = -qJ(6) * t634 + t831;
t515 = qJD(5) * t598 + t763 * t581 - t774 * t923;
t514 = qJD(5) * t792 - t581 * t923 - t667 * t857 - t763 * t774;
t513 = -qJ(6) * t597 + t887;
t512 = -pkin(5) * t770 - qJ(6) * t598 + t833;
t498 = t892 - t944;
t497 = t838 + t951;
t485 = -qJ(6) * t515 - qJD(6) * t597 + t800;
t484 = pkin(5) * t865 + t514 * qJ(6) - t598 * qJD(6) + t782;
t483 = -qJ(6) * t502 - qJD(6) * t567 - t826;
t482 = t677 * pkin(5) + t501 * qJ(6) - qJD(6) * t928 + t783;
t1 = [t821 * MDP(2) + t822 * MDP(3) + ((t726 * t798 + t811) * t770 + (qJD(2) * t816 - t726 * t779 - t909) * t766) * MDP(21) + (qJDD(1) * t760 + 0.2e1 * t766 * t840) * MDP(4) + (-t542 * t667 + t581 * t815) * MDP(18) + (-t542 * t770 - t581 * t726 - t667 * t692 + t815 * t865) * MDP(20) + ((-t643 * t860 + t852) * t726 - t881 * t692 + t825 * t770 - t546 * t865 + t646 * t815 + t713 * t542 - t585 * t667 + t645 * t581 - g(1) * t659 - g(2) * t661) * MDP(24) + (t766 * t808 + t770 * t789) * MDP(9) + (-t766 * t789 + t770 * t808) * MDP(10) + (-t786 * t726 + t830 * t692 - t646 * t816 + t713 * t811 - g(1) * t660 - g(2) * t662 + (qJD(2) * t804 - t785) * t770 + (t545 * qJD(2) - t645 * t641 + t911) * t766) * MDP(23) + (t581 * t816 + t667 * t811 - t815 * t793 + (t779 * t815 - t912) * t766) * MDP(19) + (-t482 * t598 - t483 * t597 - t484 * t928 - t485 * t567 + t494 * t514 - t496 * t515 + t501 * t512 - t502 * t513 + t766 * t821) * MDP(32) + (t501 * t597 - t502 * t598 + t514 * t567 - t515 * t928) * MDP(26) + (-t501 * t598 - t514 * t928) * MDP(25) + (t501 * t770 + t514 * t717 + t598 * t677 + t865 * t928) * MDP(27) + (-g(1) * t648 - g(2) * t650 - t639 * t501 - t504 * t865 - t583 * t514 + t516 * t598 + t563 * t928 - t677 * t887 + t717 * t800 - t770 * t826) * MDP(31) + (qJDD(2) * t766 + t770 * t772) * MDP(6) + (qJDD(2) * t770 - t766 * t772) * MDP(7) + (-t699 * t770 - t735 * t865) * MDP(15) + (-t692 * t770 - t726 * t865) * MDP(22) + (-t677 * t770 - t717 * t865) * MDP(29) + (t502 * t770 + t515 * t717 - t567 * t865 - t597 * t677) * MDP(28) + (-g(1) * t649 - g(2) * t651 + t639 * t502 + t503 * t865 + t583 * t515 + t516 * t597 + t563 * t567 + t677 * t833 - t717 * t782 - t770 * t783) * MDP(30) + ((-t612 - t847) * t770 + (t805 + t867) * t766) * MDP(13) + ((t735 * t866 + t613) * t770 + (-t806 - t868) * t766) * MDP(14) + 0.2e1 * (t752 * t766 - t855 * t873) * MDP(5) + (-(-t714 * t863 + t877) * t735 + t701 * t699 - g(1) * t682 - g(2) * t684 + ((t844 + t868) * pkin(7) + (-pkin(7) * t699 + qJD(2) * t718 - t828) * t765 + t820) * t770 + (pkin(7) * t613 + qJD(2) * t637 + t675 * t765 + t718 * t861) * t766) * MDP(16) + (t878 * t735 - t874 * t699 - g(1) * t681 - g(2) * t683 + (t718 * t856 + (-t845 + t867) * pkin(7) + t801) * t770 + (-t718 * t863 - t638 * qJD(2) + t675 * t769 + (t612 - t847) * pkin(7)) * t766) * MDP(17) + (t612 * t899 + t704 * t795) * MDP(11) + ((-t702 * t769 - t704 * t765) * t864 + (-t910 - t613 * t769 + (t702 * t765 - t906) * qJD(3)) * t766) * MDP(12) + (t483 * t513 + t496 * t485 + t482 * t512 + t494 * t484 + t491 * (t597 * pkin(5) + t639) + t532 * (t515 * pkin(5) + t563) - g(1) * (-t767 * t809 + t771 * t913) - g(2) * (t767 * t913 + t771 * t809)) * MDP(33) + qJDD(1) * MDP(1); (t634 * t677 + t717 * t889) * MDP(27) + (-t633 * t677 + t717 * t888) * MDP(28) + (-MDP(4) * t766 * t770 + MDP(5) * t873) * t773 + (-t655 * t501 + t516 * t634 - t889 * t583 - t885 * t677 + t717 * t936 - t740 * t791 + t928 * t942) * MDP(31) + (-pkin(2) * t613 + t875 * t735 + t819 * t765 + (-t637 * t766 + (-pkin(7) * t702 - t718 * t765) * t770) * qJD(1) + t955 * t769) * MDP(16) + (-pkin(2) * t612 - t686 * t735 + t819 * t769 + (-t718 * t897 + t638 * t766 + (-t704 * t770 + t735 * t899) * pkin(7)) * qJD(1) - t955 * t765) * MDP(17) + (t501 * t633 - t502 * t634 + t567 * t889 - t888 * t928) * MDP(26) + (-t501 * t634 - t889 * t928) * MDP(25) + (-t482 * t634 - t483 * t633 + t494 * t889 - t496 * t888 + t501 * t537 - t502 * t538 - t567 * t895 - t770 * t822 - t894 * t928 - t915) * MDP(32) + (t735 * MDP(15) - MDP(20) * t815 - MDP(21) * t816 + t726 * MDP(22) - t545 * MDP(23) + t546 * MDP(24) - MDP(27) * t928 + t567 * MDP(28) + t717 * MDP(29) - t503 * MDP(30) + t504 * MDP(31)) * t870 + (-t692 * t813 + t726 * t778) * MDP(21) + (t829 * t692 + t744 * t811 + t585 * t813 + t695 * t816 - t645 * t653 + t962 * t726 + (t726 * t876 + t804) * qJD(4) + (-t816 * t922 + t804) * qJD(3) + t791 * t755) * MDP(23) + MDP(7) * t752 + MDP(6) * t854 + (t655 * t502 + t516 * t633 + t677 * t831 - t741 * t914 + (g(1) * t904 + g(2) * t905) * t766 + t937 * t717 + t888 * t583 + t942 * t567) * MDP(30) + (-t876 * t692 + t744 * t542 + t911 + (-t721 * t860 + t959) * t726 - t960 * t645 + t938 * t815 - t791 * t754) * MDP(24) + (-t815 * t960 + t912) * MDP(18) + (-t542 * t813 - t705 * t811 + t815 * t956 - t816 * t960) * MDP(19) + (t726 * t960 + t909) * MDP(20) + (t483 * t538 + t482 * t537 + t491 * (t633 * pkin(5) + t655) - g(3) * t818 + (pkin(4) * t778 + pkin(5) * t888 + t938) * t532 + t895 * t496 + t894 * t494 + t822 * (t658 * t766 + t756 * t770)) * MDP(33) + ((-t704 * t766 + t735 * t897) * qJD(1) + t806) * MDP(13) + ((t702 * t766 - t735 * t900) * qJD(1) + t805) * MDP(14) + ((t612 + t908) * t769 + (-t613 + t907) * t765) * MDP(12) + (-t735 * t906 + t910) * MDP(11) + (t766 * t796 - t745 - t914) * MDP(9) + (t915 + (-pkin(7) * qJDD(1) + t796) * t770) * MDP(10) + qJDD(2) * MDP(8); (-t886 * t726 + (-t764 * t692 - t704 * t815 + t726 * t858) * pkin(3) + t930) * MDP(24) + (-t588 * t928 - t676 * t677 + t717 * t941 + t948) * MDP(31) + (t501 * t668 - t502 * t676 - t567 * t884 - t883 * t928 + t954) * MDP(32) + (t832 * t726 + (t768 * t692 + t704 * t816 + t726 * t860) * pkin(3) + t931) * MDP(23) + (-t588 * t567 + t824 * t677 + t717 * t940 + t932) * MDP(30) + t704 * t702 * MDP(11) + (-t702 ^ 2 + t704 ^ 2) * MDP(12) + t699 * MDP(15) + (g(1) * t684 - g(2) * t682 + g(3) * t899 - t637 * t735 + t702 * t718 - t801) * MDP(17) + (-t613 - t907) * MDP(14) + (t612 - t908) * MDP(13) + (t483 * t676 + t482 * t668 - t532 * (pkin(5) * t928 + t588) - g(1) * (-t664 * t896 + t665 * t767) - g(2) * (-t664 * t898 - t665 * t771) + t664 * t915 + t884 * t496 + t883 * t494) * MDP(33) + (-g(1) * t683 + g(2) * t681 - t638 * t735 - t704 * t718 + (t828 + t915) * t765 - t820) * MDP(16) + t963; (-t546 * t726 + t931) * MDP(23) + (-t545 * t726 + t930) * MDP(24) + (t838 * t717 + (-t567 * t815 + t677 * t923 + t717 * t857) * pkin(4) + t932) * MDP(30) + (-t892 * t717 + (-t763 * t677 + t717 * t842 - t815 * t928) * pkin(4) + t948) * MDP(31) + (t497 * t928 + t498 * t567 + t742 * t501 + (-t502 * t763 + (-t567 * t923 + t763 * t928) * qJD(5)) * pkin(4) + t954) * MDP(32) + (t482 * t742 - t496 * t498 - t494 * t497 - pkin(5) * t943 - g(1) * (t696 * t896 + t697 * t767) - g(2) * (t696 * t898 - t697 * t771) - t696 * t915 + (t483 * t763 - t532 * t815 + (-t494 * t763 + t496 * t923) * qJD(5)) * pkin(4)) * MDP(33) + t963; (-t504 * t717 + t932) * MDP(30) + (-t503 * t717 + t948) * MDP(31) + (pkin(5) * t501 - t567 * t893) * MDP(32) + (t893 * t496 + (t482 + t933 - t943) * pkin(5)) * MDP(33) + t964; (-t564 - t925) * MDP(32) + (t494 * t928 + t496 * t567 + t491 - t791) * MDP(33);];
tau  = t1;

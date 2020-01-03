% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPRPR3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:44
% EndTime: 2020-01-03 11:36:57
% DurationCPUTime: 13.63s
% Computational Cost: add. (51409->499), mult. (76171->682), div. (0->0), fcn. (45228->10), ass. (0->346)
t904 = qJD(1) + qJD(3);
t900 = t904 ^ 2;
t913 = sin(qJ(3));
t901 = qJDD(1) + qJDD(3);
t916 = cos(qJ(3));
t971 = t916 * t901;
t863 = t913 * t900 - t971;
t907 = g(1) - qJDD(2);
t1014 = pkin(6) * t863 - t913 * t907;
t977 = t913 * t901;
t860 = t916 * t900 + t977;
t909 = sin(pkin(8));
t911 = cos(pkin(8));
t808 = t909 * t860 + t911 * t863;
t830 = pkin(6) * t860 - t916 * t907;
t1021 = qJ(2) * t808 + t1014 * t911 + t909 * t830;
t804 = t911 * t860 - t909 * t863;
t734 = qJ(2) * t804 - t1014 * t909 + t911 * t830;
t914 = sin(qJ(1));
t917 = cos(qJ(1));
t750 = t804 * t917 - t808 * t914;
t1033 = -pkin(5) * t750 + t1021 * t914 - t917 * t734;
t1012 = t804 * t914 + t808 * t917;
t1032 = pkin(5) * t1012 + t1021 * t917 + t914 * t734;
t886 = t917 * g(2) + t914 * g(3);
t870 = qJDD(1) * pkin(1) - t886;
t885 = t914 * g(2) - t917 * g(3);
t918 = qJD(1) ^ 2;
t871 = -t918 * pkin(1) - t885;
t815 = t909 * t870 + t911 * t871;
t811 = -t918 * pkin(2) + t815;
t937 = t911 * t870 - t909 * t871;
t926 = qJDD(1) * pkin(2) + t937;
t753 = t913 * t811 - t916 * t926;
t754 = t916 * t811 + t913 * t926;
t956 = t913 * t753 + t916 * t754;
t689 = t916 * t753 - t913 * t754;
t984 = t909 * t689;
t1018 = t911 * t956 + t984;
t982 = t911 * t689;
t638 = -t909 * t956 + t982;
t1031 = t1018 * t917 + t638 * t914;
t1030 = t1018 * t914 - t638 * t917;
t908 = sin(pkin(9));
t903 = t908 ^ 2;
t910 = cos(pkin(9));
t920 = t910 ^ 2;
t849 = (t903 + t920) * t910 * t900;
t960 = t910 * t971;
t818 = -t913 * t849 + t960;
t820 = t916 * t849 + t910 * t977;
t762 = t911 * t818 - t909 * t820;
t764 = t909 * t818 + t911 * t820;
t1023 = t762 * t917 - t764 * t914;
t1022 = t762 * t914 + t764 * t917;
t955 = t911 * t815 - t909 * t937;
t758 = -t909 * t815 - t911 * t937;
t991 = t758 * t917;
t1017 = t914 * t955 - t991;
t992 = t758 * t914;
t1016 = -t917 * t955 - t992;
t874 = t909 * qJDD(1) + t911 * t918;
t841 = qJ(2) * t874 - t911 * t907;
t875 = t911 * qJDD(1) - t909 * t918;
t938 = -qJ(2) * t875 - t909 * t907;
t999 = -t917 * t874 - t914 * t875;
t1015 = pkin(5) * t999 - t917 * t841 + t914 * t938;
t1000 = t914 * t874 - t917 * t875;
t1013 = pkin(5) * t1000 + t914 * t841 + t917 * t938;
t985 = t908 * t910;
t825 = t860 * t985;
t878 = t900 * t985;
t826 = -t913 * t878 + t908 * t960;
t770 = t911 * t825 + t909 * t826;
t773 = t909 * t825 - t911 * t826;
t1010 = t770 * t917 - t773 * t914;
t1009 = t770 * t914 + t773 * t917;
t912 = sin(qJ(5));
t915 = cos(qJ(5));
t833 = (qJD(5) * t904 * t915 + t901 * t912) * t908;
t882 = t910 * t904 - qJD(5);
t966 = t882 * t904 * t908;
t846 = t915 * t966;
t793 = t846 - t833;
t1006 = t912 * t793;
t892 = t910 * t907;
t744 = -t900 * pkin(3) + t901 * qJ(4) + t754;
t997 = 2 * qJD(4);
t952 = t904 * t997 + t744;
t719 = t908 * t952 + t892;
t986 = t908 * t907;
t720 = t910 * t952 - t986;
t672 = t908 * t719 + t910 * t720;
t988 = t904 * t912;
t963 = t908 * t988;
t873 = qJD(5) * t963;
t965 = t882 * t988;
t990 = t901 * t915;
t795 = t908 * (t965 + t990) - t873;
t879 = t882 ^ 2;
t996 = pkin(1) * t907;
t995 = pkin(4) * t908;
t994 = pkin(4) * t910;
t989 = t903 * t900;
t739 = -t901 * pkin(3) - t900 * qJ(4) + qJDD(4) + t753;
t735 = t908 * t739;
t983 = t910 * t901;
t880 = -qJDD(5) + t983;
t987 = t908 * t880;
t889 = t908 * t901;
t736 = t910 * t739;
t944 = -pkin(7) * t908 - t994;
t935 = t744 + (t944 * t904 + t997) * t904;
t702 = t908 * t935 + t892;
t981 = t912 * t702;
t964 = t900 * t912 * t915;
t868 = t903 * t964;
t831 = -t868 + t880;
t980 = t912 * t831;
t832 = -t868 - t880;
t979 = t912 * t832;
t978 = t913 * t739;
t975 = t915 * t702;
t974 = t915 * t831;
t973 = t915 * t832;
t972 = t916 * t739;
t967 = -pkin(3) * t739 + qJ(4) * t672;
t906 = t915 ^ 2;
t962 = t906 * t989;
t961 = t908 * t983;
t959 = pkin(3) * t983 - qJ(4) * t849 - t736;
t642 = t913 * t672 - t972;
t958 = pkin(2) * t642 + t967;
t876 = -t914 * qJDD(1) - t917 * t918;
t957 = pkin(5) * t876 + t917 * g(1);
t703 = t910 * t935 - t986;
t723 = t901 * t944 + t739;
t657 = t912 * t703 - t915 * t723;
t658 = t915 * t703 + t912 * t723;
t627 = t912 * t657 + t915 * t658;
t953 = -t914 * t885 - t886 * t917;
t902 = t908 * t903;
t951 = t902 * t964;
t950 = t882 * t963;
t822 = -t962 - t879;
t775 = t915 * t822 + t980;
t651 = -pkin(4) * t775 + t658;
t679 = -pkin(7) * t775 + t975;
t776 = -t912 * t822 + t974;
t722 = t910 * t776 + t795 * t908;
t949 = -pkin(3) * t775 + qJ(4) * t722 + t910 * t651 + t908 * t679;
t905 = t912 ^ 2;
t881 = t905 * t989;
t836 = -t881 - t879;
t782 = t912 * t836 + t973;
t653 = -pkin(4) * t782 + t657;
t681 = -pkin(7) * t782 + t981;
t785 = t915 * t836 - t979;
t730 = t910 * t785 - t908 * t793;
t948 = -pkin(3) * t782 + qJ(4) * t730 + t910 * t653 + t908 * t681;
t888 = t903 * t901;
t890 = t920 * t901;
t856 = t890 + t888;
t891 = t920 * t900;
t858 = t891 + t989;
t947 = pkin(3) * t858 + qJ(4) * t856 + t672;
t946 = pkin(2) * t818 + t959;
t945 = -pkin(4) * t702 + pkin(7) * t627;
t943 = -pkin(2) * t863 - t753;
t942 = t910 * t868;
t683 = t913 * t722 - t916 * t775;
t941 = pkin(2) * t683 + t949;
t694 = t913 * t730 - t916 * t782;
t940 = pkin(2) * t694 + t948;
t800 = t913 * t856 + t916 * t858;
t939 = pkin(2) * t800 + t947;
t626 = -t915 * t657 + t912 * t658;
t671 = t910 * t719 - t908 * t720;
t936 = t917 * t885 - t914 * t886;
t848 = (t908 * t920 + t902) * t900;
t934 = -pkin(3) * t889 + qJ(4) * t848 + t735;
t816 = t913 * t848 - t908 * t971;
t931 = pkin(2) * t816 + t934;
t792 = t846 + t833;
t794 = -t873 + (-t965 + t990) * t908;
t740 = -t912 * t792 - t915 * t794;
t619 = -pkin(7) * t740 - t626;
t742 = -t915 * t792 + t912 * t794;
t842 = t881 + t962;
t709 = t910 * t742 - t908 * t842;
t930 = qJ(4) * t709 + t908 * t619 + (-pkin(3) - t994) * t740;
t929 = -pkin(4) * t795 + pkin(7) * t776 + t981;
t928 = pkin(4) * t793 + pkin(7) * t785 - t975;
t674 = t913 * t709 - t916 * t740;
t927 = pkin(2) * t674 + t930;
t925 = pkin(4) * t842 + pkin(7) * t742 + t627;
t613 = t910 * t627 + t908 * t702;
t924 = qJ(4) * t613 + (-pkin(3) + t944) * t626;
t599 = t913 * t613 - t916 * t626;
t923 = pkin(2) * t599 + t924;
t922 = -pkin(2) * t860 - t754;
t877 = t917 * qJDD(1) - t914 * t918;
t872 = 0.2e1 * t961;
t866 = t910 * t880;
t859 = -t891 + t989;
t857 = t890 - t888;
t852 = pkin(5) * t877 + t914 * g(1);
t843 = -t881 + t962;
t837 = t879 - t962;
t835 = t881 - t879;
t834 = t915 * t889 - t873;
t819 = t916 * t848 + t908 * t977;
t803 = t916 * t857 + t913 * t859;
t802 = t916 * t856 - t913 * t858;
t801 = t913 * t857 - t916 * t859;
t798 = (-t905 - t906) * t966;
t791 = pkin(1) * t875 + t937;
t790 = -pkin(1) * t874 - t815;
t789 = t912 * t834 - t906 * t966;
t788 = -t915 * t833 - t905 * t966;
t787 = (t834 + t950) * t915;
t784 = t915 * t835 + t980;
t783 = -t912 * t837 + t973;
t781 = t912 * t835 - t974;
t780 = t915 * t837 + t979;
t778 = -t913 * t798 - t916 * t987;
t777 = t916 * t798 - t913 * t987;
t769 = t910 * t787 + t951;
t768 = -t1006 * t910 - t951;
t767 = t908 * t787 - t942;
t766 = -t1006 * t908 + t942;
t763 = -t909 * t816 + t911 * t819;
t760 = t911 * t816 + t909 * t819;
t755 = pkin(1) * t758;
t752 = qJ(2) * t955 + t996;
t748 = -t909 * t801 + t911 * t803;
t747 = -t909 * t800 + t911 * t802;
t746 = t911 * t801 + t909 * t803;
t745 = t911 * t800 + t909 * t802;
t743 = t915 * t793 - t912 * t795;
t741 = t915 * t795 + t1006;
t729 = t910 * t784 - t908 * t792;
t728 = t910 * t783 + t908 * t794;
t727 = t908 * t783 - t910 * t794;
t726 = t908 * t784 + t910 * t792;
t725 = t908 * t785 + t910 * t793;
t721 = t908 * t776 - t795 * t910;
t714 = t916 * t769 + t913 * t789;
t713 = t916 * t768 + t913 * t788;
t712 = t913 * t769 - t916 * t789;
t711 = t913 * t768 - t916 * t788;
t710 = t910 * t743 + t908 * t843;
t708 = t908 * t743 - t910 * t843;
t707 = t908 * t742 + t910 * t842;
t705 = -t909 * t777 + t911 * t778;
t704 = t911 * t777 + t909 * t778;
t701 = -pkin(1) * t808 + t943;
t700 = -pkin(1) * t804 + t922;
t699 = t914 * t760 - t917 * t763;
t698 = t917 * t760 + t914 * t763;
t697 = t916 * t730 + t913 * t782;
t696 = t916 * t729 + t913 * t781;
t695 = t916 * t728 + t913 * t780;
t693 = t913 * t729 - t916 * t781;
t692 = t913 * t728 - t916 * t780;
t686 = pkin(2) * t689;
t685 = pkin(2) * t907 + pkin(6) * t956;
t684 = t916 * t722 + t913 * t775;
t677 = t916 * t710 + t913 * t741;
t676 = t916 * t709 + t913 * t740;
t675 = t913 * t710 - t916 * t741;
t668 = pkin(1) * t762 + t946;
t667 = pkin(1) * t760 + t931;
t666 = -t909 * t712 + t911 * t714;
t665 = -t909 * t711 + t911 * t713;
t664 = t911 * t712 + t909 * t714;
t663 = t911 * t711 + t909 * t713;
t662 = -pkin(6) * t816 - t913 * t720 + t910 * t972;
t661 = -pkin(6) * t818 - t913 * t719 + t908 * t972;
t660 = pkin(6) * t819 + t916 * t720 + t910 * t978;
t659 = -pkin(6) * t820 + t916 * t719 + t908 * t978;
t655 = -pkin(6) * t800 + t916 * t671;
t654 = pkin(6) * t802 + t913 * t671;
t649 = -t909 * t694 + t911 * t697;
t648 = -t909 * t693 + t911 * t696;
t647 = -t909 * t692 + t911 * t695;
t646 = t911 * t694 + t909 * t697;
t645 = t911 * t693 + t909 * t696;
t644 = t911 * t692 + t909 * t695;
t643 = t916 * t672 + t978;
t640 = -pkin(3) * t725 - t928;
t635 = -pkin(3) * t721 - t929;
t634 = -t909 * t683 + t911 * t684;
t633 = t911 * t683 + t909 * t684;
t632 = pkin(1) * t745 + t939;
t631 = -t909 * t675 + t911 * t677;
t630 = -t909 * t674 + t911 * t676;
t629 = t911 * t675 + t909 * t677;
t628 = t911 * t674 + t909 * t676;
t625 = -pkin(1) * t638 - t686;
t623 = -qJ(2) * t760 - t909 * t660 + t911 * t662;
t622 = -qJ(2) * t762 - t909 * t659 + t911 * t661;
t621 = qJ(2) * t763 + t911 * t660 + t909 * t662;
t620 = -qJ(2) * t764 + t911 * t659 + t909 * t661;
t618 = -qJ(4) * t725 - t908 * t653 + t910 * t681;
t616 = -qJ(4) * t721 - t908 * t651 + t910 * t679;
t615 = -qJ(2) * t745 - t909 * t654 + t911 * t655;
t614 = qJ(2) * t747 + t911 * t654 + t909 * t655;
t612 = t908 * t627 - t910 * t702;
t610 = -t909 * t642 + t911 * t643;
t609 = t911 * t642 + t909 * t643;
t608 = pkin(6) * t982 + qJ(2) * t638 - t909 * t685;
t607 = pkin(6) * t984 + qJ(2) * t1018 + t911 * t685 + t996;
t606 = -pkin(3) * t707 - t925;
t605 = -pkin(6) * t642 - (pkin(3) * t913 - qJ(4) * t916) * t671;
t604 = -qJ(4) * t707 + t910 * t619 + t740 * t995;
t603 = pkin(6) * t643 - (-pkin(3) * t916 - qJ(4) * t913 - pkin(2)) * t671;
t602 = -pkin(6) * t694 + t916 * t618 - t913 * t640;
t601 = pkin(1) * t646 + t940;
t600 = t916 * t613 + t913 * t626;
t597 = -pkin(6) * t683 + t916 * t616 - t913 * t635;
t596 = -pkin(2) * t725 + pkin(6) * t697 + t913 * t618 + t916 * t640;
t595 = pkin(1) * t633 + t941;
t594 = -pkin(2) * t721 + pkin(6) * t684 + t913 * t616 + t916 * t635;
t593 = -pkin(3) * t612 - t945;
t592 = pkin(1) * t609 + t958;
t591 = pkin(1) * t628 + t927;
t590 = -qJ(4) * t612 + (-pkin(7) * t910 + t995) * t626;
t589 = -pkin(6) * t674 + t916 * t604 - t913 * t606;
t588 = -pkin(2) * t707 + pkin(6) * t676 + t913 * t604 + t916 * t606;
t587 = -t909 * t599 + t911 * t600;
t586 = t911 * t599 + t909 * t600;
t585 = -qJ(2) * t609 - t909 * t603 + t911 * t605;
t584 = -qJ(2) * t646 - t909 * t596 + t911 * t602;
t583 = pkin(1) * t671 + qJ(2) * t610 + t911 * t603 + t909 * t605;
t582 = -pkin(1) * t725 + qJ(2) * t649 + t911 * t596 + t909 * t602;
t581 = -qJ(2) * t633 - t909 * t594 + t911 * t597;
t580 = -pkin(1) * t721 + qJ(2) * t634 + t911 * t594 + t909 * t597;
t579 = -pkin(6) * t599 + t916 * t590 - t913 * t593;
t578 = -qJ(2) * t628 - t909 * t588 + t911 * t589;
t577 = -pkin(1) * t707 + qJ(2) * t630 + t911 * t588 + t909 * t589;
t576 = -pkin(2) * t612 + pkin(6) * t600 + t913 * t590 + t916 * t593;
t575 = pkin(1) * t586 + t923;
t574 = -qJ(2) * t586 - t909 * t576 + t911 * t579;
t573 = -pkin(1) * t612 + qJ(2) * t587 + t911 * t576 + t909 * t579;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t886, t885, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t791, t790, 0, -t755, 0, 0, 0, 0, 0, t901, t701, t700, 0, t625, t888, t872, 0, t890, 0, 0, t668, t667, t632, t592, t767, t708, t727, t766, t726, t866, t601, t595, t591, t575; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t876, 0, t877, 0, t957, -t852, -t936, -pkin(5) * t936, 0, 0, -t999, 0, -t1000, 0, t1015, t1013, -t1016, -pkin(5) * t1016 + qJ(2) * t992 + t917 * t752, 0, 0, t750, 0, -t1012, 0, t1033, t1032, t1031, pkin(5) * t1031 + t917 * t607 + t914 * t608, t1010, t746 * t917 + t748 * t914, t698, -t1010, -t1023, 0, -pkin(5) * t1022 + t917 * t620 + t914 * t622, -pkin(5) * t699 + t621 * t917 + t623 * t914, t914 * t615 + t917 * t614 - pkin(5) * (t745 * t914 - t747 * t917), t914 * t585 + t917 * t583 - pkin(5) * (t609 * t914 - t610 * t917), t664 * t917 + t666 * t914, t629 * t917 + t631 * t914, t644 * t917 + t647 * t914, t663 * t917 + t665 * t914, t645 * t917 + t648 * t914, t704 * t917 + t705 * t914, t914 * t584 + t917 * t582 - pkin(5) * (t646 * t914 - t649 * t917), t914 * t581 + t917 * t580 - pkin(5) * (t633 * t914 - t634 * t917), t914 * t578 + t917 * t577 - pkin(5) * (t628 * t914 - t630 * t917), t914 * t574 + t917 * t573 - pkin(5) * (t586 * t914 - t587 * t917); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t877, 0, -t876, 0, t852, t957, t953, pkin(5) * t953, 0, 0, t1000, 0, -t999, 0, -t1013, t1015, t1017, pkin(5) * t1017 - qJ(2) * t991 + t914 * t752, 0, 0, t1012, 0, t750, 0, -t1032, t1033, t1030, pkin(5) * t1030 + t914 * t607 - t917 * t608, t1009, t746 * t914 - t748 * t917, t699, -t1009, -t1022, 0, pkin(5) * t1023 + t914 * t620 - t917 * t622, pkin(5) * t698 + t621 * t914 - t623 * t917, -t917 * t615 + t914 * t614 + pkin(5) * (t745 * t917 + t747 * t914), -t917 * t585 + t914 * t583 + pkin(5) * (t609 * t917 + t610 * t914), t664 * t914 - t666 * t917, t629 * t914 - t631 * t917, t644 * t914 - t647 * t917, t663 * t914 - t665 * t917, t645 * t914 - t648 * t917, t704 * t914 - t705 * t917, -t917 * t584 + t914 * t582 + pkin(5) * (t646 * t917 + t649 * t914), -t917 * t581 + t914 * t580 + pkin(5) * (t633 * t917 + t634 * t914), -t917 * t578 + t914 * t577 + pkin(5) * (t628 * t917 + t630 * t914), -t917 * t574 + t914 * t573 + pkin(5) * (t586 * t917 + t587 * t914); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t918, 0, 0, -g(1), t886, 0, 0, 0, t875, 0, -t874, 0, t938, t841, t758, qJ(2) * t758, 0, 0, -t808, 0, -t804, 0, t1021, t734, t638, t608, -t773, t748, t763, t773, t764, 0, t622, t623, t615, t585, t666, t631, t647, t665, t648, t705, t584, t581, t578, t574; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t918, 0, qJDD(1), 0, g(1), 0, -t885, 0, 0, 0, t874, 0, t875, 0, -t841, t938, t955, t752, 0, 0, t804, 0, -t808, 0, -t734, t1021, t1018, t607, t770, t746, t760, -t770, -t762, 0, t620, t621, t614, t583, t664, t629, t644, t663, t645, t704, t582, t580, t577, t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t886, t885, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t791, t790, 0, -t755, 0, 0, 0, 0, 0, t901, t701, t700, 0, t625, t888, t872, 0, t890, 0, 0, t668, t667, t632, t592, t767, t708, t727, t766, t726, t866, t601, t595, t591, t575; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t918, 0, 0, -t907, -t937, 0, 0, 0, -t863, 0, -t860, 0, t1014, t830, t689, pkin(6) * t689, t826, t803, t819, -t826, t820, 0, t661, t662, t655, t605, t714, t677, t695, t713, t696, t778, t602, t597, t589, t579; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t918, 0, qJDD(1), 0, t907, 0, t815, 0, 0, 0, t860, 0, -t863, 0, -t830, t1014, t956, t685, t825, t801, t816, -t825, -t818, 0, t659, t660, t654, t603, t712, t675, t692, t711, t693, t777, t596, t594, t588, t576; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t937, -t815, 0, 0, 0, 0, 0, 0, 0, t901, t943, t922, 0, -t686, t888, t872, 0, t890, 0, 0, t946, t931, t939, t958, t767, t708, t727, t766, t726, t866, t940, t941, t927, t923; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t901, 0, -t900, 0, 0, -t907, t753, 0, t961, t857, t848, -t961, t849, 0, t735, t736, t671, qJ(4) * t671, t769, t710, t728, t768, t729, -t987, t618, t616, t604, t590; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t900, 0, t901, 0, t907, 0, t754, 0, t878, -t859, -t889, -t878, -t983, 0, t719, t720, 0, pkin(3) * t671, -t789, -t741, -t780, -t788, -t781, t798, t640, t635, t606, t593; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t901, -t753, -t754, 0, 0, t888, t872, 0, t890, 0, 0, t959, t934, t947, t967, t767, t708, t727, t766, t726, t866, t948, t949, t930, t924; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t889, t983, t878, 0, t891, 0, 0, t739, t719, 0, t787, t743, t783, -t1006, t784, 0, t681, t679, t619, -pkin(7) * t626; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t889, -t989, t983, -t878, 0, -t739, 0, t720, 0, -t868, -t843, -t794, t868, t792, t880, t653, t651, -pkin(4) * t740, -pkin(4) * t626; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t878, t859, t889, t878, t983, 0, -t719, -t720, 0, 0, t789, t741, t780, t788, t781, -t798, t928, t929, t925, t945; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t834, t793, t832, -t950, t835, t950, 0, t702, t657, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t846, t795, t837, -t833, -t831, t846, -t702, 0, t658, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t868, t843, t794, -t868, -t792, -t880, -t657, -t658, 0, 0;];
m_new_reg = t1;
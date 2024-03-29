% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPRPR2
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPRPR2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:18
% EndTime: 2022-01-23 09:19:30
% DurationCPUTime: 11.75s
% Computational Cost: add. (59264->491), mult. (88678->688), div. (0->0), fcn. (55375->10), ass. (0->335)
t933 = cos(qJ(3));
t980 = qJD(1) + qJD(3);
t979 = t980 ^ 2;
t922 = qJDD(1) + qJDD(3);
t930 = sin(qJ(3));
t984 = t930 * t922;
t889 = t933 * t979 + t984;
t983 = t933 * t922;
t892 = t930 * t979 - t983;
t926 = sin(pkin(8));
t928 = cos(pkin(8));
t832 = t889 * t928 - t892 * t926;
t836 = t889 * t926 + t892 * t928;
t931 = sin(qJ(1));
t934 = cos(qJ(1));
t1019 = t832 * t934 - t836 * t931;
t924 = g(3) - qJDD(2);
t1021 = pkin(6) * t892 - t924 * t930;
t860 = pkin(6) * t889 - t924 * t933;
t1030 = qJ(2) * t836 + t1021 * t928 + t860 * t926;
t762 = qJ(2) * t832 - t1021 * t926 + t860 * t928;
t1042 = pkin(5) * t1019 - t1030 * t931 + t934 * t762;
t1031 = t832 * t931 + t836 * t934;
t1041 = pkin(5) * t1031 + t1030 * t934 + t931 * t762;
t907 = g(1) * t934 + g(2) * t931;
t936 = qJD(1) ^ 2;
t945 = pkin(1) * t936 + t907;
t906 = g(1) * t931 - t934 * g(2);
t949 = qJDD(1) * pkin(1) + t906;
t843 = t926 * t949 - t928 * t945;
t839 = -t936 * pkin(2) + t843;
t942 = t926 * t945 + t928 * t949;
t941 = qJDD(1) * pkin(2) + t942;
t780 = t839 * t930 - t933 * t941;
t781 = t933 * t839 + t930 * t941;
t973 = t780 * t930 + t933 * t781;
t727 = t780 * t933 - t781 * t930;
t997 = t727 * t926;
t1025 = t928 * t973 + t997;
t996 = t727 * t928;
t664 = t926 * t973 - t996;
t1039 = t1025 * t934 - t664 * t931;
t1038 = t1025 * t931 + t664 * t934;
t927 = cos(pkin(9));
t938 = t927 ^ 2;
t916 = t938 * t979;
t925 = sin(pkin(9));
t923 = t925 ^ 2;
t968 = t927 * t979;
t879 = t916 * t927 + t923 * t968;
t846 = -t879 * t930 + t927 * t983;
t848 = t879 * t933 + t927 * t984;
t797 = t846 * t928 - t848 * t926;
t800 = t846 * t926 + t848 * t928;
t1033 = t797 * t934 - t800 * t931;
t1032 = t797 * t931 + t800 * t934;
t972 = t928 * t843 - t926 * t942;
t791 = -t843 * t926 - t928 * t942;
t992 = t791 * t934;
t1027 = t931 * t972 - t992;
t993 = t791 * t931;
t1026 = t934 * t972 + t993;
t932 = cos(qJ(5));
t966 = t932 * t980;
t929 = sin(qJ(5));
t967 = t929 * t980;
t872 = t925 * t967 - t927 * t966;
t874 = t925 * t966 + t927 * t967;
t825 = t874 * t872;
t1006 = qJDD(5) - t825;
t1024 = t1006 * t929;
t1023 = t1006 * t932;
t950 = -t979 * pkin(3) + t922 * qJ(4) + 0.2e1 * qJD(4) * t980 + t781;
t1022 = t922 * pkin(7) + t950;
t898 = t926 * qJDD(1) + t928 * t936;
t899 = t928 * qJDD(1) - t926 * t936;
t1007 = t934 * t898 + t899 * t931;
t871 = qJ(2) * t898 - t924 * t928;
t953 = -qJ(2) * t899 - t924 * t926;
t1020 = pkin(5) * t1007 + t934 * t871 - t931 * t953;
t850 = -t898 * t931 + t934 * t899;
t1018 = -pkin(5) * t850 + t931 * t871 + t934 * t953;
t852 = t889 * t927 * t925;
t902 = t925 * t968;
t915 = t927 * t922;
t976 = t925 * t915;
t853 = -t902 * t930 + t933 * t976;
t806 = t852 * t928 + t853 * t926;
t809 = t852 * t926 - t853 * t928;
t1017 = t806 * t934 - t809 * t931;
t1016 = t806 * t931 + t809 * t934;
t917 = t927 * t924;
t755 = t925 * t950 + t917;
t985 = t925 * t924;
t756 = t927 * t950 - t985;
t699 = t925 * t755 + t927 * t756;
t969 = t923 * t979;
t887 = t969 + t916;
t866 = t872 ^ 2;
t867 = t874 ^ 2;
t1004 = pkin(1) * t924;
t745 = -t917 + (pkin(4) * t968 - t1022) * t925;
t746 = -pkin(4) * t916 + t1022 * t927 - t985;
t678 = -t932 * t745 + t746 * t929;
t679 = t745 * t929 + t746 * t932;
t644 = -t678 * t932 + t679 * t929;
t1003 = pkin(4) * t644;
t913 = t925 * t922;
t864 = t929 * t913 - t932 * t915;
t865 = (t932 * t925 + t929 * t927) * t922;
t765 = -t864 * t929 - t865 * t932;
t1002 = pkin(4) * t765;
t999 = t644 * t925;
t998 = t644 * t927;
t770 = -t922 * pkin(3) - t979 * qJ(4) + qJDD(4) + t780;
t750 = -pkin(4) * t915 - t887 * pkin(7) + t770;
t995 = t750 * t929;
t994 = t750 * t932;
t818 = qJDD(5) + t825;
t991 = t818 * t929;
t990 = t818 * t932;
t763 = t925 * t770;
t764 = t927 * t770;
t982 = -pkin(3) * t770 + qJ(4) * t699;
t863 = t872 * qJD(5);
t981 = t874 * qJD(5);
t978 = t930 * t825;
t977 = t933 * t825;
t975 = pkin(3) * t915 - qJ(4) * t879 - t764;
t668 = t699 * t930 - t770 * t933;
t974 = pkin(2) * t668 + t982;
t645 = t678 * t929 + t932 * t679;
t970 = -t906 * t931 - t934 * t907;
t767 = -t864 * t932 + t865 * t929;
t805 = -t866 - t867;
t624 = -pkin(4) * t805 + pkin(7) * t767 + t645;
t635 = -pkin(7) * t765 - t644;
t716 = -t765 * t925 + t767 * t927;
t964 = -pkin(3) * t805 + qJ(4) * t716 + t927 * t624 + t925 * t635;
t935 = qJD(5) ^ 2;
t816 = -t935 - t866;
t758 = t816 * t932 - t1024;
t820 = t864 + 0.2e1 * t981;
t680 = -pkin(4) * t820 + pkin(7) * t758 - t994;
t757 = t816 * t929 + t1023;
t701 = -pkin(7) * t757 + t995;
t713 = -t757 * t925 + t758 * t927;
t963 = -pkin(3) * t820 + qJ(4) * t713 + t927 * t680 + t925 * t701;
t856 = -t867 - t935;
t787 = -t856 * t929 - t990;
t822 = -0.2e1 * t863 + t865;
t693 = -pkin(4) * t822 + pkin(7) * t787 + t995;
t784 = t856 * t932 - t991;
t718 = -pkin(7) * t784 + t994;
t735 = -t784 * t925 + t787 * t927;
t962 = -pkin(3) * t822 + qJ(4) * t735 + t927 * t693 + t925 * t718;
t912 = t923 * t922;
t914 = t938 * t922;
t885 = t914 + t912;
t961 = pkin(3) * t887 + qJ(4) * t885 + t699;
t960 = pkin(2) * t846 + t975;
t901 = qJDD(1) * t934 - t931 * t936;
t959 = -pkin(5) * t901 - g(3) * t931;
t958 = -pkin(2) * t892 - t780;
t674 = t716 * t930 - t805 * t933;
t957 = pkin(2) * t674 + t964;
t682 = t713 * t930 - t820 * t933;
t956 = pkin(2) * t682 + t963;
t703 = t735 * t930 - t822 * t933;
t955 = pkin(2) * t703 + t962;
t828 = t885 * t930 + t887 * t933;
t954 = pkin(2) * t828 + t961;
t698 = t755 * t927 - t756 * t925;
t952 = t906 * t934 - t907 * t931;
t878 = t887 * t925;
t951 = -pkin(3) * t913 + qJ(4) * t878 + t763;
t948 = pkin(4) * t757 - t678;
t844 = t878 * t930 - t925 * t983;
t947 = pkin(2) * t844 + t951;
t620 = t645 * t927 - t999;
t631 = -pkin(4) * t750 + pkin(7) * t645;
t946 = -pkin(3) * t750 - pkin(7) * t999 + qJ(4) * t620 + t927 * t631;
t616 = t620 * t930 - t750 * t933;
t944 = pkin(2) * t616 + t946;
t943 = pkin(4) * t784 - t679;
t940 = -pkin(2) * t889 - t781;
t900 = qJDD(1) * t931 + t934 * t936;
t897 = 0.2e1 * t976;
t888 = t969 - t916;
t886 = t914 - t912;
t881 = -pkin(5) * t900 + g(3) * t934;
t855 = -t867 + t935;
t854 = t866 - t935;
t847 = t878 * t933 + t925 * t984;
t831 = t886 * t933 + t888 * t930;
t830 = t885 * t933 - t887 * t930;
t829 = t886 * t930 - t888 * t933;
t824 = t867 - t866;
t823 = -t863 + t865;
t821 = -t864 - t981;
t813 = pkin(1) * t899 + t942;
t812 = -pkin(1) * t898 - t843;
t811 = (-t872 * t932 + t874 * t929) * qJD(5);
t810 = (-t872 * t929 - t874 * t932) * qJD(5);
t802 = t823 * t932 - t929 * t981;
t799 = -t844 * t926 + t847 * t928;
t798 = t823 * t929 + t932 * t981;
t795 = t844 * t928 + t847 * t926;
t794 = -t821 * t929 + t932 * t863;
t793 = t821 * t932 + t929 * t863;
t788 = pkin(1) * t791;
t786 = -t855 * t929 + t1023;
t785 = t854 * t932 - t991;
t783 = t855 * t932 + t1024;
t782 = t854 * t929 + t990;
t779 = qJ(2) * t972 + t1004;
t775 = -t829 * t926 + t831 * t928;
t774 = -t828 * t926 + t830 * t928;
t773 = t829 * t928 + t831 * t926;
t772 = t828 * t928 + t830 * t926;
t768 = -t820 * t932 - t822 * t929;
t766 = -t820 * t929 + t822 * t932;
t754 = -t810 * t925 + t811 * t927;
t753 = t810 * t927 + t811 * t925;
t748 = qJDD(5) * t930 + t754 * t933;
t747 = -qJDD(5) * t933 + t754 * t930;
t743 = -pkin(1) * t836 + t958;
t742 = -pkin(1) * t832 + t940;
t741 = -t795 * t931 + t799 * t934;
t740 = t795 * t934 + t799 * t931;
t739 = -t798 * t925 + t802 * t927;
t738 = -t793 * t925 + t794 * t927;
t737 = t798 * t927 + t802 * t925;
t736 = t793 * t927 + t794 * t925;
t734 = -t783 * t925 + t786 * t927;
t733 = -t782 * t925 + t785 * t927;
t732 = t784 * t927 + t787 * t925;
t731 = t783 * t927 + t786 * t925;
t730 = t782 * t927 + t785 * t925;
t724 = pkin(2) * t727;
t723 = pkin(2) * t924 + pkin(6) * t973;
t722 = t734 * t933 + t865 * t930;
t721 = t733 * t933 - t864 * t930;
t720 = t734 * t930 - t865 * t933;
t719 = t733 * t930 + t864 * t933;
t717 = -t766 * t925 + t768 * t927;
t715 = t766 * t927 + t768 * t925;
t714 = t765 * t927 + t767 * t925;
t712 = t757 * t927 + t758 * t925;
t708 = t739 * t933 + t978;
t707 = t738 * t933 - t978;
t706 = t739 * t930 - t977;
t705 = t738 * t930 + t977;
t704 = t735 * t933 + t822 * t930;
t695 = pkin(1) * t797 + t960;
t694 = pkin(1) * t795 + t947;
t692 = -pkin(6) * t844 - t756 * t930 + t933 * t764;
t691 = -pkin(6) * t846 - t755 * t930 + t933 * t763;
t690 = pkin(6) * t847 + t756 * t933 + t930 * t764;
t689 = -pkin(6) * t848 + t755 * t933 + t930 * t763;
t687 = -t747 * t926 + t748 * t928;
t686 = t747 * t928 + t748 * t926;
t685 = t717 * t933 + t824 * t930;
t684 = t717 * t930 - t824 * t933;
t683 = t713 * t933 + t820 * t930;
t675 = t716 * t933 + t805 * t930;
t672 = -pkin(6) * t828 + t698 * t933;
t671 = pkin(6) * t830 + t698 * t930;
t670 = -pkin(3) * t714 - t1002;
t669 = t699 * t933 + t770 * t930;
t662 = -t720 * t926 + t722 * t928;
t661 = -t719 * t926 + t721 * t928;
t660 = t720 * t928 + t722 * t926;
t659 = t719 * t928 + t721 * t926;
t658 = pkin(1) * t772 + t954;
t657 = -t706 * t926 + t708 * t928;
t656 = -t705 * t926 + t707 * t928;
t655 = t706 * t928 + t708 * t926;
t654 = t705 * t928 + t707 * t926;
t653 = -t703 * t926 + t704 * t928;
t652 = t703 * t928 + t704 * t926;
t651 = -pkin(3) * t732 - t943;
t650 = -t684 * t926 + t685 * t928;
t649 = t684 * t928 + t685 * t926;
t648 = pkin(1) * t664 - t724;
t647 = -t682 * t926 + t683 * t928;
t646 = t682 * t928 + t683 * t926;
t643 = -t674 * t926 + t675 * t928;
t642 = t674 * t928 + t675 * t926;
t641 = -pkin(3) * t712 - t948;
t640 = -qJ(2) * t795 - t690 * t926 + t692 * t928;
t639 = -qJ(2) * t797 - t689 * t926 + t691 * t928;
t638 = qJ(2) * t799 + t690 * t928 + t692 * t926;
t637 = -qJ(2) * t800 + t689 * t928 + t691 * t926;
t636 = -qJ(4) * t732 - t693 * t925 + t718 * t927;
t633 = -qJ(2) * t772 - t671 * t926 + t672 * t928;
t632 = qJ(2) * t774 + t671 * t928 + t672 * t926;
t629 = -t668 * t926 + t669 * t928;
t628 = t668 * t928 + t669 * t926;
t627 = pkin(6) * t996 - qJ(2) * t664 - t723 * t926;
t626 = -qJ(4) * t712 - t680 * t925 + t701 * t927;
t625 = pkin(6) * t997 + qJ(2) * t1025 + t723 * t928 + t1004;
t622 = -pkin(6) * t668 - (pkin(3) * t930 - qJ(4) * t933) * t698;
t621 = pkin(6) * t669 - (-pkin(3) * t933 - qJ(4) * t930 - pkin(2)) * t698;
t619 = t645 * t925 + t998;
t617 = t620 * t933 + t750 * t930;
t614 = pkin(1) * t652 + t955;
t613 = pkin(1) * t628 + t974;
t612 = -pkin(6) * t703 + t636 * t933 - t651 * t930;
t611 = -pkin(2) * t732 + pkin(6) * t704 + t636 * t930 + t651 * t933;
t610 = pkin(1) * t646 + t956;
t609 = -pkin(6) * t682 + t626 * t933 - t641 * t930;
t608 = -qJ(4) * t714 - t624 * t925 + t635 * t927;
t607 = -pkin(3) * t619 - t1003;
t606 = -pkin(2) * t712 + pkin(6) * t683 + t626 * t930 + t641 * t933;
t605 = -pkin(6) * t674 + t608 * t933 - t670 * t930;
t604 = -pkin(7) * t998 - qJ(4) * t619 - t631 * t925;
t603 = -pkin(2) * t714 + pkin(6) * t675 + t608 * t930 + t670 * t933;
t602 = -t616 * t926 + t617 * t928;
t601 = t616 * t928 + t617 * t926;
t600 = pkin(1) * t642 + t957;
t599 = -qJ(2) * t628 - t621 * t926 + t622 * t928;
t598 = pkin(1) * t698 + qJ(2) * t629 + t621 * t928 + t622 * t926;
t597 = -qJ(2) * t652 - t611 * t926 + t612 * t928;
t596 = -pkin(1) * t732 + qJ(2) * t653 + t611 * t928 + t612 * t926;
t595 = -qJ(2) * t646 - t606 * t926 + t609 * t928;
t594 = -pkin(1) * t712 + qJ(2) * t647 + t606 * t928 + t609 * t926;
t593 = -qJ(2) * t642 - t603 * t926 + t605 * t928;
t592 = -pkin(6) * t616 + t604 * t933 - t607 * t930;
t591 = -pkin(1) * t714 + qJ(2) * t643 + t603 * t928 + t605 * t926;
t590 = pkin(1) * t601 + t944;
t589 = -pkin(2) * t619 + pkin(6) * t617 + t604 * t930 + t607 * t933;
t588 = -qJ(2) * t601 - t589 * t926 + t592 * t928;
t587 = -pkin(1) * t619 + qJ(2) * t602 + t589 * t928 + t592 * t926;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t901, 0, -t900, 0, t959, -t881, -t952, -pkin(5) * t952, 0, 0, t850, 0, -t1007, 0, t1018, t1020, -t1027, -pkin(5) * t1027 + qJ(2) * t992 - t931 * t779, 0, 0, -t1031, 0, -t1019, 0, t1041, t1042, -t1038, -pkin(5) * t1038 - t931 * t625 + t934 * t627, -t1016, -t773 * t931 + t775 * t934, t741, t1016, t1032, 0, -pkin(5) * t1033 - t931 * t637 + t934 * t639, -pkin(5) * t740 - t638 * t931 + t640 * t934, t934 * t633 - t931 * t632 - pkin(5) * (t772 * t934 + t774 * t931), t934 * t599 - t931 * t598 - pkin(5) * (t628 * t934 + t629 * t931), -t655 * t931 + t657 * t934, -t649 * t931 + t650 * t934, -t660 * t931 + t662 * t934, -t654 * t931 + t656 * t934, -t659 * t931 + t661 * t934, -t686 * t931 + t687 * t934, t934 * t595 - t931 * t594 - pkin(5) * (t646 * t934 + t647 * t931), t934 * t597 - t931 * t596 - pkin(5) * (t652 * t934 + t653 * t931), t934 * t593 - t931 * t591 - pkin(5) * (t642 * t934 + t643 * t931), t934 * t588 - t931 * t587 - pkin(5) * (t601 * t934 + t602 * t931); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t900, 0, t901, 0, t881, t959, t970, pkin(5) * t970, 0, 0, t1007, 0, t850, 0, -t1020, t1018, t1026, pkin(5) * t1026 + qJ(2) * t993 + t934 * t779, 0, 0, t1019, 0, -t1031, 0, -t1042, t1041, t1039, pkin(5) * t1039 + t934 * t625 + t931 * t627, t1017, t773 * t934 + t775 * t931, t740, -t1017, -t1033, 0, -pkin(5) * t1032 + t934 * t637 + t931 * t639, pkin(5) * t741 + t638 * t934 + t640 * t931, t931 * t633 + t934 * t632 + pkin(5) * (-t772 * t931 + t774 * t934), t931 * t599 + t934 * t598 + pkin(5) * (-t628 * t931 + t629 * t934), t655 * t934 + t657 * t931, t649 * t934 + t650 * t931, t660 * t934 + t662 * t931, t654 * t934 + t656 * t931, t659 * t934 + t661 * t931, t686 * t934 + t687 * t931, t931 * t595 + t934 * t594 + pkin(5) * (-t646 * t931 + t647 * t934), t931 * t597 + t934 * t596 + pkin(5) * (-t652 * t931 + t653 * t934), t931 * t593 + t934 * t591 + pkin(5) * (-t642 * t931 + t643 * t934), t931 * t588 + t934 * t587 + pkin(5) * (-t601 * t931 + t602 * t934); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t906, t907, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t813, t812, 0, -t788, 0, 0, 0, 0, 0, t922, t743, t742, 0, t648, t912, t897, 0, t914, 0, 0, t695, t694, t658, t613, t737, t715, t731, t736, t730, t753, t610, t614, t600, t590; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t936, 0, 0, -g(3), -t906, 0, 0, 0, t899, 0, -t898, 0, t953, t871, t791, qJ(2) * t791, 0, 0, -t836, 0, -t832, 0, t1030, t762, -t664, t627, -t809, t775, t799, t809, t800, 0, t639, t640, t633, t599, t657, t650, t662, t656, t661, t687, t595, t597, t593, t588; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t936, 0, qJDD(1), 0, g(3), 0, -t907, 0, 0, 0, t898, 0, t899, 0, -t871, t953, t972, t779, 0, 0, t832, 0, -t836, 0, -t762, t1030, t1025, t625, t806, t773, t795, -t806, -t797, 0, t637, t638, t632, t598, t655, t649, t660, t654, t659, t686, t594, t596, t591, t587; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t906, t907, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t813, t812, 0, -t788, 0, 0, 0, 0, 0, t922, t743, t742, 0, t648, t912, t897, 0, t914, 0, 0, t695, t694, t658, t613, t737, t715, t731, t736, t730, t753, t610, t614, t600, t590; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t936, 0, 0, -t924, -t942, 0, 0, 0, -t892, 0, -t889, 0, t1021, t860, t727, pkin(6) * t727, t853, t831, t847, -t853, t848, 0, t691, t692, t672, t622, t708, t685, t722, t707, t721, t748, t609, t612, t605, t592; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t936, 0, qJDD(1), 0, t924, 0, t843, 0, 0, 0, t889, 0, -t892, 0, -t860, t1021, t973, t723, t852, t829, t844, -t852, -t846, 0, t689, t690, t671, t621, t706, t684, t720, t705, t719, t747, t606, t611, t603, t589; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t942, -t843, 0, 0, 0, 0, 0, 0, 0, t922, t958, t940, 0, -t724, t912, t897, 0, t914, 0, 0, t960, t947, t954, t974, t737, t715, t731, t736, t730, t753, t956, t955, t957, t944; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t922, 0, -t979, 0, 0, -t924, t780, 0, t976, t886, t878, -t976, t879, 0, t763, t764, t698, qJ(4) * t698, t739, t717, t734, t738, t733, t754, t626, t636, t608, t604; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t979, 0, t922, 0, t924, 0, t781, 0, t902, -t888, -t913, -t902, -t915, 0, t755, t756, 0, pkin(3) * t698, -t825, -t824, -t865, t825, t864, -qJDD(5), t641, t651, t670, t607; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t922, -t780, -t781, 0, 0, t912, t897, 0, t914, 0, 0, t975, t951, t961, t982, t737, t715, t731, t736, t730, t753, t963, t962, t964, t946; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t913, t915, t902, 0, t916, 0, 0, t770, t755, 0, t802, t768, t786, t794, t785, t811, t701, t718, t635, -pkin(7) * t644; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t913, -t969, t915, -t902, 0, -t770, 0, t756, 0, t798, t766, t783, t793, t782, t810, t680, t693, t624, t631; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t902, t888, t913, t902, t915, 0, -t755, -t756, 0, 0, t825, t824, t865, -t825, -t864, qJDD(5), t948, t943, t1002, t1003; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t823, -t820, t1006, t863, t854, -t863, 0, t750, t678, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t981, t822, t855, t821, t818, -t981, -t750, 0, t679, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t825, t824, t865, -t825, -t864, qJDD(5), -t678, -t679, 0, 0;];
m_new_reg = t1;

% Calculate inertial parameters regressor of coriolis matrix for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRRR8_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:53
% EndTime: 2019-12-31 22:26:26
% DurationCPUTime: 23.32s
% Computational Cost: add. (25338->756), mult. (50077->1004), div. (0->0), fcn. (54959->8), ass. (0->595)
t618 = sin(qJ(4));
t968 = cos(qJ(5));
t773 = t968 * t618;
t617 = sin(qJ(5));
t621 = cos(qJ(4));
t863 = t617 * t621;
t582 = t773 + t863;
t978 = t582 / 0.2e1;
t772 = t968 * t621;
t864 = t617 * t618;
t1007 = t772 - t864;
t1066 = -t1007 / 0.2e1;
t622 = cos(qJ(2));
t998 = -pkin(7) - pkin(6);
t597 = t998 * t622;
t969 = cos(qJ(3));
t586 = t969 * t597;
t619 = sin(qJ(3));
t620 = sin(qJ(2));
t858 = t619 * t620;
t725 = t998 * t858;
t1012 = -t586 + t725;
t1028 = t1012 * t618;
t780 = t969 * t620;
t857 = t619 * t622;
t583 = t780 + t857;
t612 = -t622 * pkin(2) - pkin(1);
t580 = -t622 * t969 + t858;
t964 = t580 * pkin(3);
t689 = t612 + t964;
t997 = -pkin(9) - pkin(8);
t239 = -t1028 + (t583 * t997 + t689) * t621;
t963 = t580 * pkin(4);
t629 = t239 + t963;
t217 = t968 * t629;
t1027 = t1012 * t621;
t960 = t583 * pkin(8);
t656 = t689 - t960;
t278 = t618 * t656 + t1027;
t885 = t583 * t618;
t240 = -pkin(9) * t885 + t278;
t867 = t617 * t240;
t144 = -t217 + t867;
t777 = t968 * t239;
t159 = t777 - t867;
t1063 = t144 + t159;
t626 = t617 * t629;
t776 = t968 * t240;
t145 = t776 + t626;
t868 = t617 * t239;
t158 = -t776 - t868;
t1056 = t145 + t158;
t1065 = t580 / 0.2e1;
t795 = t969 * pkin(2);
t610 = -t795 - pkin(3);
t957 = t621 * pkin(4);
t594 = t610 - t957;
t1064 = -t594 / 0.2e1;
t611 = -pkin(3) - t957;
t974 = -t611 / 0.2e1;
t596 = t997 * t621;
t735 = -t617 * t596 - t997 * t773;
t1035 = t735 / 0.2e1;
t958 = t619 * pkin(2);
t609 = pkin(8) + t958;
t956 = pkin(9) + t609;
t567 = t956 * t621;
t741 = t956 * t618;
t736 = t617 * t567 + t968 * t741;
t993 = t736 / 0.2e1;
t1062 = t1056 * t978 + t1063 * t1066;
t801 = qJD(2) + qJD(3);
t800 = qJD(4) + qJD(5);
t413 = t1007 * t583;
t492 = -t596 * t968 + t864 * t997;
t1061 = t492 * t1065 + t413 * t974;
t437 = t567 * t968 - t617 * t741;
t1060 = t413 * t1064 + t437 * t1065;
t410 = t582 * t583;
t1047 = t580 * t1035 + t410 * t974;
t1048 = t410 * t1064 + t580 * t993;
t1037 = -t437 / 0.2e1;
t1036 = -t492 / 0.2e1;
t985 = -t735 / 0.2e1;
t754 = t1028 / 0.2e1;
t734 = -t619 * t597 - t998 * t780;
t1031 = t734 * t621;
t761 = t1031 / 0.2e1;
t995 = -t413 / 0.2e1;
t207 = -t1007 * t995 - t410 * t978;
t1057 = t800 * t207;
t441 = t618 * t580;
t1024 = -pkin(4) * t441 + t1012;
t377 = pkin(4) * t885 + t734;
t1054 = t1024 * t377;
t1053 = t1024 * t582;
t1026 = t207 * qJD(1);
t732 = t801 * t582;
t1051 = -t1007 * t732 - t1026;
t888 = t1007 * t582;
t1050 = t801 * t888 + t1026;
t1049 = t1024 * t1007;
t916 = t377 * t582;
t323 = t916 / 0.2e1;
t1032 = t734 * t618;
t961 = t583 * pkin(3);
t962 = t580 * pkin(8);
t473 = t961 + t962;
t453 = t621 * t473;
t312 = t453 + t1032;
t444 = t621 * t580;
t566 = t583 * pkin(4);
t713 = pkin(9) * t444 + t566;
t223 = t312 + t713;
t778 = t968 * t223;
t452 = t618 * t473;
t313 = -t1031 + t452;
t797 = pkin(9) * t441;
t251 = t797 + t313;
t865 = t617 * t251;
t659 = -t865 / 0.2e1 + t778 / 0.2e1;
t76 = t323 + t659 - t1061;
t967 = pkin(2) * t620;
t448 = t473 + t967;
t425 = t621 * t448;
t284 = t425 + t1032;
t222 = t284 + t713;
t779 = t968 * t222;
t424 = t618 * t448;
t285 = -t1031 + t424;
t244 = t797 + t285;
t866 = t617 * t244;
t660 = -t866 / 0.2e1 + t779 / 0.2e1;
t66 = t323 + t660 - t1060;
t775 = t968 * t244;
t716 = -t775 / 0.2e1;
t870 = t617 * t222;
t662 = -t870 / 0.2e1 + t716;
t1046 = t662 - t1048;
t774 = t968 * t251;
t715 = -t774 / 0.2e1;
t869 = t617 * t223;
t661 = -t869 / 0.2e1 + t715;
t1045 = t661 - t1047;
t1044 = t660 + t1060;
t1043 = t659 + t1061;
t643 = t621 * t656;
t277 = -t643 + t1028;
t1042 = (-t277 + t1028) * t583;
t1041 = (-t278 + t1027) * t583;
t409 = t582 * t580;
t1040 = t1024 * t410 - t144 * t583 - t377 * t409;
t528 = t617 * t441;
t412 = -t580 * t772 + t528;
t1039 = t1024 * t413 - t145 * t583 + t377 * t412;
t1038 = -0.2e1 * t583;
t965 = t1012 * pkin(3);
t568 = t780 / 0.2e1 + t857 / 0.2e1;
t834 = qJD(1) * t583;
t766 = t580 * t834;
t366 = t568 * qJD(4) + t766;
t898 = t1012 * t734;
t1016 = t801 * t583;
t1025 = t580 * t1016;
t838 = qJD(1) * t413;
t1023 = t207 * t801 - t410 * t838;
t1022 = t801 * t734;
t1021 = t800 * t736;
t1020 = t800 * t735;
t1019 = t800 * t437;
t1018 = t800 * t492;
t1017 = 0.2e1 * t583;
t785 = t958 / 0.2e1;
t471 = t801 * t580;
t781 = t969 * t580;
t884 = t583 * t619;
t886 = t583 * t609;
t887 = t580 * t610;
t1006 = t964 / 0.2e1 - t887 / 0.2e1 - t886 / 0.2e1 + (t884 / 0.2e1 - t781 / 0.2e1) * pkin(2);
t316 = qJD(5) * t568 + t366;
t577 = t583 ^ 2;
t799 = -t580 ^ 2 + t577;
t883 = t583 * t621;
t759 = t883 / 0.2e1;
t727 = pkin(4) * t759;
t796 = t968 / 0.2e1;
t862 = t618 * t410;
t1002 = pkin(4) * t862 / 0.2e1 - t1007 * t727 + t796 * t566;
t173 = -t1007 * t410 - t582 * t413;
t197 = t410 ^ 2 - t413 ^ 2;
t41 = qJD(1) * t197 + t173 * t801;
t399 = t1007 ^ 2 - t582 ^ 2;
t98 = t173 * qJD(1) + t399 * t801;
t615 = t618 ^ 2;
t616 = t621 ^ 2;
t440 = (t615 / 0.2e1 - t616 / 0.2e1) * t583;
t859 = t618 * t621;
t769 = qJD(1) * t859;
t224 = t440 * t801 + t577 * t769;
t601 = t616 - t615;
t402 = t769 * t1038 + t601 * t801;
t999 = pkin(8) / 0.2e1;
t156 = t778 - t865;
t996 = -t156 / 0.2e1;
t990 = t437 / 0.2e1;
t987 = t734 / 0.2e1;
t983 = t492 / 0.2e1;
t720 = t969 * t968;
t782 = t621 * t969;
t522 = (-t617 * t782 - t618 * t720) * pkin(2);
t982 = -t522 / 0.2e1;
t783 = t618 * t969;
t523 = (-t617 * t783 + t621 * t720) * pkin(2);
t981 = t523 / 0.2e1;
t980 = t528 / 0.2e1;
t979 = -t580 / 0.2e1;
t717 = t586 / 0.2e1;
t977 = t594 / 0.2e1;
t976 = -t609 / 0.2e1;
t975 = t610 / 0.2e1;
t973 = -t617 / 0.2e1;
t972 = t617 / 0.2e1;
t971 = -t618 / 0.2e1;
t970 = -t621 / 0.2e1;
t966 = pkin(4) * t618;
t959 = t617 * pkin(4);
t745 = t983 + t1036;
t746 = t985 + t1035;
t748 = t990 + t1037;
t750 = -t736 / 0.2e1 + t993;
t62 = (-t745 - t748) * t582 - (-t746 - t750) * t1007;
t955 = t62 * qJD(4);
t954 = pkin(2) * qJD(3);
t953 = pkin(3) * qJD(3);
t952 = pkin(4) * qJD(4);
t951 = qJD(2) * pkin(2);
t51 = (-t963 + t754 - t643 / 0.2e1 + pkin(9) * t759 + t239 / 0.2e1) * t617;
t948 = qJD(1) * t51;
t798 = pkin(4) * t883;
t918 = t377 * t413;
t70 = t158 * t580 + t410 * t798 + t918;
t947 = qJD(1) * t70;
t919 = t377 * t410;
t71 = -t159 * t580 + t413 * t798 - t919;
t946 = qJD(1) * t71;
t87 = t144 * t580 - t919;
t945 = qJD(1) * t87;
t88 = -t145 * t580 + t918;
t944 = qJD(1) * t88;
t148 = t779 - t866;
t941 = t148 * t582;
t149 = t775 + t870;
t940 = t149 * t1007;
t939 = t156 * t582;
t157 = t774 + t869;
t938 = t157 * t1007;
t699 = t144 * t412 + t145 * t409;
t18 = -t148 * t413 - t149 * t410 + t699;
t935 = t18 * qJD(1);
t21 = -t156 * t413 - t157 * t410 + t699;
t934 = t21 * qJD(1);
t22 = -t1056 * t413 - t1063 * t410;
t933 = t22 * qJD(1);
t25 = -t144 * t148 + t145 * t149 + t1054;
t932 = t25 * qJD(1);
t26 = -t144 * t156 + t145 * t157 + t1054;
t931 = t26 * qJD(1);
t27 = -t144 * t158 + t145 * t159 + t377 * t798;
t930 = t27 * qJD(1);
t929 = t284 * t618;
t928 = t285 * t621;
t30 = t148 * t580 + t1040;
t927 = t30 * qJD(1);
t31 = -t149 * t580 + t1039;
t926 = t31 * qJD(1);
t925 = t312 * t618;
t924 = t313 * t621;
t32 = t156 * t580 + t1040;
t923 = t32 * qJD(1);
t33 = -t157 * t580 + t1039;
t922 = t33 * qJD(1);
t917 = t377 * t1007;
t911 = t736 * t412;
t909 = t736 * t583;
t907 = t437 * t409;
t904 = t437 * t583;
t901 = t735 * t412;
t899 = t735 * t583;
t893 = t492 * t409;
t890 = t492 * t583;
t788 = -t963 / 0.2e1;
t675 = -t217 / 0.2e1 + t968 * t788;
t53 = t777 / 0.2e1 + t675;
t889 = t53 * qJD(1);
t882 = t594 * t409;
t880 = t594 * t412;
t878 = t594 * t1007;
t877 = t594 * t582;
t876 = t611 * t409;
t874 = t611 * t412;
t872 = t611 * t1007;
t871 = t611 * t582;
t674 = (t277 * t621 - t278 * t618) * t580;
t67 = (t284 * t621 + t285 * t618) * t583 + t674;
t855 = t67 * qJD(1);
t72 = (t312 * t621 + t313 * t618) * t583 + t674;
t854 = t72 * qJD(1);
t77 = -t277 * t284 + t278 * t285 + t898;
t853 = t77 * qJD(1);
t78 = -t277 * t312 + t278 * t313 + t898;
t852 = t78 * qJD(1);
t96 = t1042 + (t284 - t1032) * t580;
t851 = t96 * qJD(1);
t97 = t1041 + (-t285 - t1031) * t580;
t850 = t97 * qJD(1);
t848 = t800 * t173;
t845 = (-t1007 * t968 - t582 * t617) * t952;
t724 = -t795 / 0.2e1;
t686 = -t610 / 0.2e1 + t724;
t627 = (t785 + t976 + t999) * t583 + (-pkin(3) / 0.2e1 + t686) * t580;
t152 = t618 * t627;
t844 = qJD(1) * t152;
t175 = t409 * t413 - t410 * t412;
t843 = qJD(1) * t175;
t193 = -t277 * t580 + t734 * t885;
t842 = qJD(1) * t193;
t194 = -t278 * t580 + t734 * t883;
t841 = qJD(1) * t194;
t839 = qJD(1) * t410;
t422 = t580 * t967 + t583 * t612;
t837 = qJD(1) * t422;
t423 = -t580 * t612 + t583 * t967;
t836 = qJD(1) * t423;
t835 = qJD(1) * t580;
t833 = qJD(1) * t612;
t832 = qJD(1) * t622;
t831 = qJD(2) * t594;
t830 = qJD(2) * t610;
t829 = qJD(3) * t611;
t828 = qJD(3) * t612;
t827 = qJD(4) * t618;
t826 = qJD(4) * t621;
t825 = qJD(5) * t594;
t824 = qJD(5) * t611;
t102 = t1042 + (t312 - t1032) * t580;
t823 = t102 * qJD(1);
t103 = t1041 + (-t313 - t1031) * t580;
t822 = t103 * qJD(1);
t199 = t612 * t967;
t819 = t199 * qJD(1);
t210 = t409 * t580 - t410 * t583;
t816 = t210 * qJD(1);
t211 = t412 * t580 + t413 * t583;
t815 = t211 * qJD(1);
t658 = t863 / 0.2e1 + t773 / 0.2e1;
t286 = (t978 + t658) * t580;
t275 = t286 * qJD(1);
t714 = -t772 / 0.2e1;
t287 = t980 + (t714 + t1066) * t580;
t276 = t287 * qJD(1);
t364 = t799 * t618;
t814 = t364 * qJD(1);
t365 = t799 * t621;
t813 = t365 * qJD(1);
t812 = t799 * qJD(1);
t811 = t440 * qJD(1);
t810 = t440 * qJD(4);
t809 = t441 * qJD(1);
t430 = t444 * qJD(1);
t451 = t601 * t577;
t808 = t451 * qJD(1);
t486 = t717 - t586 / 0.2e1;
t807 = t486 * qJD(1);
t806 = t568 * qJD(1);
t602 = -t620 ^ 2 + t622 ^ 2;
t804 = t602 * qJD(1);
t803 = t620 * qJD(2);
t802 = t622 * qJD(2);
t794 = pkin(1) * t620 * qJD(1);
t793 = pkin(1) * t832;
t792 = t619 * t951;
t791 = t619 * t954;
t790 = t966 / 0.2e1;
t787 = -t566 / 0.2e1;
t786 = pkin(8) * t971;
t322 = t917 / 0.2e1;
t784 = t413 * t790 + t582 * t727 + t322;
t771 = t580 * t833;
t770 = t583 * t833;
t768 = qJD(4) * t580 * t583;
t605 = t618 * t826;
t765 = t620 * t802;
t764 = -t917 / 0.2e1;
t763 = -t916 / 0.2e1;
t762 = t377 * t971;
t760 = -t883 / 0.2e1;
t752 = -t424 / 0.2e1 + t761;
t751 = t993 + t985;
t749 = t1037 + t983;
t747 = t452 / 0.2e1 - t1031 / 0.2e1;
t743 = t977 + t974;
t742 = t611 / 0.2e1 + t977;
t740 = t969 * qJD(2);
t739 = t969 * qJD(3);
t738 = t968 * qJD(4);
t737 = t968 * qJD(5);
t730 = t801 * t621;
t729 = t800 * t580;
t728 = t800 * t582;
t726 = -qJD(4) - t835;
t723 = t577 * t605;
t719 = t783 / 0.2e1;
t718 = -t782 / 0.2e1;
t712 = t801 * t958;
t711 = -t960 + t964;
t708 = t801 * t859;
t706 = t618 * t730;
t704 = -qJD(5) + t726;
t165 = t594 * t966;
t637 = t1063 * t1037 + t1056 * t993;
t665 = t148 * t796 + t149 * t972;
t3 = (t594 * t760 + t665 + t762) * pkin(4) + t637;
t703 = -t3 * qJD(1) + t165 * qJD(2);
t191 = t437 * t523 - t522 * t736 + t594 * t958;
t624 = t1024 * t977 + t144 * t982 + t145 * t981 + t157 * t990 + t377 * t785 + t736 * t996;
t642 = t1024 * t974 + t1035 * t148 + t1036 * t149;
t2 = t624 + t642;
t702 = t2 * qJD(1) + t191 * qJD(2);
t283 = t1007 * t523 - t522 * t582;
t682 = -t410 * t981 + t413 * t982;
t8 = (t996 + t148 / 0.2e1) * t582 - (-t157 / 0.2e1 + t149 / 0.2e1) * t1007 + t751 * t412 - t749 * t409 + t682;
t701 = -qJD(1) * t8 - qJD(2) * t283;
t696 = t928 - t929;
t695 = t924 - t925;
t694 = -t886 - t887;
t683 = t924 / 0.2e1 - t925 / 0.2e1;
t623 = t683 * t609 + (t278 * t782 / 0.2e1 + t277 * t719 + t619 * t987) * pkin(2) + t1012 * t975;
t684 = -t928 / 0.2e1 + t929 / 0.2e1;
t35 = t965 / 0.2e1 + t684 * pkin(8) + t623;
t657 = (t615 + t616) * t969;
t431 = (t609 * t657 + t610 * t619) * pkin(2);
t693 = -t35 * qJD(1) - t431 * qJD(2);
t631 = (t413 * t971 + (t582 * t970 + t973) * t583) * pkin(4) + t764;
t37 = t631 + t1046;
t558 = t582 * t966;
t386 = t558 + t878;
t692 = qJD(1) * t37 - qJD(2) * t386;
t628 = (-t862 / 0.2e1 + (-t1007 * t970 + t796) * t583) * pkin(4) + t763;
t38 = t628 + t1044;
t557 = t1007 * t966;
t385 = -t557 + t877;
t691 = qJD(1) * t38 - qJD(2) * t385;
t565 = t657 * pkin(2);
t90 = (t313 / 0.2e1 - t285 / 0.2e1) * t621 + (-t312 / 0.2e1 + t284 / 0.2e1) * t618;
t690 = -qJD(1) * t90 - qJD(2) * t565;
t688 = t583 * t726;
t687 = -t720 / 0.2e1;
t685 = t961 / 0.2e1 + t962 / 0.2e1;
t681 = t580 * t976 + t583 * t975;
t64 = t764 + t1046;
t680 = qJD(1) * t64 - t1007 * t831;
t63 = t763 + t1044;
t679 = qJD(1) * t63 - t582 * t831;
t678 = t621 * t688;
t638 = -t618 * t681 + t761;
t183 = t638 - t752;
t677 = -qJD(1) * t183 - t621 * t830;
t653 = t681 * t621;
t185 = -t425 / 0.2e1 + t653 + (t987 - t734 / 0.2e1) * t618;
t676 = -qJD(1) * t185 - t618 * t830;
t673 = pkin(3) / 0.2e1 + t686;
t672 = t522 * t1065 + t410 * t785;
t671 = t413 * t785 + t523 * t979;
t56 = -t409 * t743 - t583 * t751 + t672;
t670 = -qJD(1) * t56 + t1007 * t792;
t58 = t412 * t743 + t583 * t749 + t671;
t669 = -qJD(1) * t58 - t582 * t792;
t154 = t621 * t627;
t668 = -qJD(1) * t154 - t618 * t792;
t667 = t685 * t621;
t666 = t706 * t1017;
t664 = t156 * t796 + t157 * t972;
t663 = t522 * t796 + t523 * t972;
t648 = (-t409 * t973 - t968 * t412 / 0.2e1) * pkin(4);
t625 = t648 + t1062;
t9 = t410 * t750 + t413 * t748 + t625;
t655 = qJD(1) * t9 - qJD(3) * t62;
t198 = t611 * t966;
t636 = t1056 * t1035 + t1063 * t1036;
t5 = (t611 * t760 + t664 + t762) * pkin(4) + t636;
t68 = (-t618 * t742 + t663) * pkin(4);
t654 = -t5 * qJD(1) - t68 * qJD(2) + t198 * qJD(3);
t11 = t410 * t746 + t413 * t745 + t625;
t652 = qJD(1) * t11 - qJD(2) * t62;
t634 = (t617 * t718 + t618 * t687) * pkin(2);
t298 = -t582 * t742 + t634;
t245 = t557 + t298;
t420 = -t557 + t871;
t44 = t628 + t1043;
t651 = qJD(1) * t44 + qJD(2) * t245 - qJD(3) * t420;
t633 = (t617 * t719 + t621 * t687) * pkin(2);
t299 = -t1007 * t742 + t633;
t246 = -t558 + t299;
t421 = t558 + t872;
t43 = t631 + t1045;
t650 = qJD(1) * t43 + qJD(2) * t246 - qJD(3) * t421;
t639 = t618 * t685 + t761;
t187 = t639 + t747;
t509 = t673 * t621;
t647 = -qJD(1) * t187 + qJD(2) * t509 + t621 * t953;
t189 = -t453 / 0.2e1 - t667;
t508 = t673 * t618;
t646 = -qJD(1) * t189 + qJD(2) * t508 + t618 * t953;
t73 = t763 + t1043;
t645 = qJD(1) * t73 + qJD(2) * t298 - t582 * t829;
t74 = t764 + t1045;
t644 = qJD(1) * t74 + qJD(2) * t299 - t1007 * t829;
t632 = t648 - t1062;
t300 = t872 / 0.2e1 + t878 / 0.2e1 + t633;
t301 = t871 / 0.2e1 + t877 / 0.2e1 + t634;
t606 = t620 * t832;
t600 = t618 * t791;
t595 = t601 * qJD(4);
t559 = t565 * qJD(3);
t547 = t582 * t791;
t546 = t1007 * t791;
t511 = pkin(2) * t718 + pkin(3) * t970 + t621 * t975;
t510 = pkin(3) * t971 + (t724 + t975) * t618;
t449 = t801 * t568;
t406 = 0.2e1 * t717 - t725;
t382 = t430 + t826;
t381 = -t809 - t827;
t358 = t377 * t790;
t340 = t706 - t811;
t339 = -t708 + t811;
t329 = t1007 * t728;
t327 = t800 * t888;
t326 = 0.2e1 * t618 * t678;
t289 = t580 * t658 + t582 * t979;
t288 = -t1007 * t979 + t580 * t714 + t980;
t280 = t616 * t766 - t810;
t279 = t615 * t766 + t810;
t274 = t283 * qJD(3);
t248 = t558 + t300;
t247 = -t557 + t301;
t238 = qJD(4) * t444 - t813;
t237 = -qJD(4) * t441 + t814;
t236 = t800 * t399;
t219 = -t728 - t275;
t218 = t1007 * t800 - t276;
t201 = -t810 + (-t616 * t834 - t708) * t580;
t200 = t810 + (-t615 * t834 + t706) * t580;
t196 = t1016 * t618 + t813;
t195 = t583 * t730 - t814;
t192 = (-qJD(4) + t835) * t859 * t1017 - t601 * t471;
t190 = t1032 + t453 / 0.2e1 - t667;
t188 = t639 - t747;
t186 = t1032 + t425 / 0.2e1 + t653;
t184 = t638 + t752;
t155 = pkin(8) * t760 + t1006 * t621 + 0.2e1 * t754;
t153 = t1006 * t618 + t583 * t786 - t1027;
t147 = t286 * t801 + t413 * t835;
t146 = t287 * t801 + t410 * t835;
t101 = -t287 * t800 - t815;
t100 = -t286 * t800 - t816;
t93 = t289 * t801 + t413 * t704;
t92 = t288 * t801 + t410 * t704;
t89 = t683 - t684;
t86 = -t412 * t838 + t1057;
t81 = t409 * t839 - t1057;
t80 = t1016 * t582 + t288 * t800 + t815;
t79 = t1007 * t1016 + t289 * t800 + t816;
t75 = t322 + t661 + t1047;
t69 = pkin(4) * t663 + (t594 + t611) * t790;
t65 = t322 + t662 + t1048;
t57 = -t904 / 0.2e1 + t880 / 0.2e1 + t1053 - t890 / 0.2e1 + t874 / 0.2e1 + t671;
t55 = -t909 / 0.2e1 - t882 / 0.2e1 - t1049 - t899 / 0.2e1 - t876 / 0.2e1 + t672;
t54 = t867 - t777 / 0.2e1 + t675;
t52 = t617 * t788 - t776 - t626 / 0.2e1 - t868 / 0.2e1;
t48 = (t732 + t838) * t412 + t1057;
t47 = -(-t1007 * t801 + t839) * t409 - t1057;
t46 = t715 + (t787 - t223 / 0.2e1) * t617 + t784 + t1047;
t45 = t1002 + t76;
t40 = t716 + (t787 - t222 / 0.2e1) * t617 + t784 + t1048;
t39 = t1002 + t66;
t36 = -t843 + t848;
t34 = t928 * t999 + t284 * t786 - t965 / 0.2e1 + t623;
t29 = t843 + t801 * (t1007 * t412 + t409 * t582) + t848;
t12 = t413 * t1036 - t492 * t995 + t632;
t10 = t413 * t1037 - t437 * t995 + t632;
t7 = t907 / 0.2e1 + t938 / 0.2e1 + t911 / 0.2e1 - t939 / 0.2e1 + t893 / 0.2e1 + t940 / 0.2e1 + t901 / 0.2e1 - t941 / 0.2e1 + t682;
t6 = pkin(4) * t664 + t611 * t727 + t358 - t636;
t4 = pkin(4) * t665 + t594 * t727 + t358 - t637;
t1 = t624 - t642;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t765, t602 * qJD(2), 0, -t765, 0, 0, -pkin(1) * t803, -pkin(1) * t802, 0, 0, -t1025, -t801 * t799, 0, t1025, 0, 0, qJD(2) * t422 + t583 * t828, qJD(2) * t423 - t580 * t828, 0, qJD(2) * t199, -t1025 * t616 - t723, -t451 * qJD(4) + t580 * t666, t365 * t801 - t618 * t768, -t1025 * t615 + t723, -t364 * t801 - t621 * t768, t1025, qJD(2) * t96 + qJD(3) * t102 + qJD(4) * t194, qJD(2) * t97 + qJD(3) * t103 - qJD(4) * t193, -qJD(2) * t67 - qJD(3) * t72, qJD(2) * t77 + qJD(3) * t78, (-t410 * t800 + t412 * t801) * t413, t175 * t801 + t197 * t800, t211 * t801 - t410 * t729, (-t409 * t801 + t413 * t800) * t410, t210 * t801 - t413 * t729, t1025, qJD(2) * t30 + qJD(3) * t32 + qJD(4) * t70 + qJD(5) * t88, qJD(2) * t31 + qJD(3) * t33 + qJD(4) * t71 + qJD(5) * t87, qJD(2) * t18 + qJD(3) * t21 + qJD(4) * t22, qJD(2) * t25 + qJD(3) * t26 + qJD(4) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t606, t804, t802, -t606, -t803, 0, -pkin(6) * t802 - t794, pkin(6) * t803 - t793, 0, 0, -t766, -t812, -t471, t766, -t1016, 0, -qJD(2) * t1012 + qJD(3) * t406 + t837, t1022 + t836, (t781 - t884) * t951, t819 + (-t1012 * t969 - t619 * t734) * t951, t201, t192, t196, t200, t195, t366, t851 + (t618 * t694 - t1027) * qJD(2) + t153 * qJD(3) + t186 * qJD(4), t850 + (t621 * t694 + t1028) * qJD(2) + t155 * qJD(3) + t184 * qJD(4), qJD(2) * t696 + t89 * qJD(3) - t855, t853 + (t1012 * t610 + t609 * t696) * qJD(2) + t34 * qJD(3), t48, t29, t80, t47, t79, t316, t927 + (-t882 - t909 - t1049) * qJD(2) + t55 * qJD(3) + t39 * qJD(4) + t66 * qJD(5), t926 + (t880 - t904 + t1053) * qJD(2) + t57 * qJD(3) + t40 * qJD(4) + t65 * qJD(5), t935 + (t907 + t911 + t940 - t941) * qJD(2) + t7 * qJD(3) + t10 * qJD(4), t932 + (t1024 * t594 - t148 * t736 + t149 * t437) * qJD(2) + t1 * qJD(3) + t4 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t766, -t812, -t471, t766, -t1016, 0, qJD(2) * t406 - qJD(3) * t1012 + t770, t1022 - t771, 0, 0, t201, t192, t196, t200, t195, t366, t823 + t153 * qJD(2) + (t618 * t711 - t1027) * qJD(3) + t190 * qJD(4), t822 + t155 * qJD(2) + (t621 * t711 + t1028) * qJD(3) + t188 * qJD(4), t89 * qJD(2) + qJD(3) * t695 - t854, t852 + t34 * qJD(2) + (pkin(8) * t695 - t965) * qJD(3), t48, t29, t80, t47, t79, t316, t923 + t55 * qJD(2) + (-t876 - t899 - t1049) * qJD(3) + t45 * qJD(4) + t76 * qJD(5), t922 + t57 * qJD(2) + (t874 - t890 + t1053) * qJD(3) + t46 * qJD(4) + t75 * qJD(5), t934 + t7 * qJD(2) + (t893 + t901 + t938 - t939) * qJD(3) + t12 * qJD(4), t931 + t1 * qJD(2) + (t1024 * t611 - t156 * t735 + t157 * t492) * qJD(3) + t6 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224, t708 * t1038 - t808, t618 * t688, t224, t678, t449, qJD(2) * t186 + qJD(3) * t190 - qJD(4) * t278 + t841, qJD(2) * t184 + qJD(3) * t188 + qJD(4) * t277 - t842, 0, 0, t1023, t41, t92, -t1023, t93, t449, qJD(2) * t39 + qJD(3) * t45 + qJD(4) * t158 + qJD(5) * t52 + t947, qJD(2) * t40 + qJD(3) * t46 - qJD(4) * t159 + qJD(5) * t54 + t946, t933 + t10 * qJD(2) + t12 * qJD(3) + (t410 * t968 - t413 * t617) * t952, t930 + t4 * qJD(2) + t6 * qJD(3) + (t158 * t968 + t159 * t617) * t952; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1023, t41, t92, -t1023, t93, t449, qJD(2) * t66 + qJD(3) * t76 + qJD(4) * t52 - qJD(5) * t145 + t944, qJD(2) * t65 + qJD(3) * t75 + qJD(4) * t54 + qJD(5) * t144 + t945, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t606, -t804, 0, t606, 0, 0, t794, t793, 0, 0, t766, t812, 0, -t766, 0, 0, qJD(3) * t486 - t837, -t836, 0, -t819, t280, t326, t238, t279, t237, -t366, qJD(3) * t152 + qJD(4) * t185 - t851, qJD(3) * t154 + qJD(4) * t183 - t850, qJD(3) * t90 + t855, qJD(3) * t35 - t853, t86, t36, t101, t81, t100, -t316, qJD(3) * t56 - qJD(4) * t38 - qJD(5) * t63 - t927, qJD(3) * t58 - qJD(4) * t37 - qJD(5) * t64 - t926, qJD(3) * t8 - qJD(4) * t9 - t935, qJD(3) * t2 - qJD(4) * t3 - t932; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t791, -pkin(2) * t739, 0, 0, t605, t595, 0, -t605, 0, 0, t610 * t827 - t621 * t791, t610 * t826 + t600, t559, t431 * qJD(3), t329, t236, 0, -t327, 0, 0, qJD(4) * t385 + t582 * t825 - t546, qJD(4) * t386 + t1007 * t825 + t547, t274, qJD(3) * t191 + qJD(4) * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t712 + t807, (-t740 - t739) * pkin(2), 0, 0, t605, t595, 0, -t605, 0, 0, qJD(4) * t510 - t621 * t712 + t844, qJD(4) * t511 + t600 - t668, t559 - t690, (-pkin(3) * t619 + pkin(8) * t657) * t954 - t693, t329, t236, 0, -t327, 0, 0, qJD(4) * t247 + qJD(5) * t301 - t546 - t670, qJD(4) * t248 + qJD(5) * t300 + t547 - t669, t274 - t701 + t955, (t492 * t523 - t522 * t735 + t611 * t958) * qJD(3) + t69 * qJD(4) + t702; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, t402, t382, t339, t381, -t806, qJD(3) * t510 - t609 * t826 - t676, qJD(3) * t511 + t609 * t827 - t677, 0, 0, t1050, t98, t218, t1051, t219, -t806, qJD(3) * t247 - t1019 - t691, qJD(3) * t248 + t1021 - t692, -t655 + t845, t69 * qJD(3) + (-t437 * t968 - t617 * t736) * t952 + t703; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1050, t98, t218, t1051, t219, -t806, qJD(3) * t301 - t1019 - t679, qJD(3) * t300 + t1021 - t680, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t766, t812, 0, -t766, 0, 0, -qJD(2) * t486 - t770, t771, 0, 0, t280, t326, t238, t279, t237, -t366, -qJD(2) * t152 + qJD(4) * t189 - t823, -qJD(2) * t154 + qJD(4) * t187 - t822, -qJD(2) * t90 + t854, -qJD(2) * t35 - t852, t86, t36, t101, t81, t100, -t316, -qJD(2) * t56 - qJD(4) * t44 - qJD(5) * t73 - t923, -qJD(2) * t58 - qJD(4) * t43 - qJD(5) * t74 - t922, -qJD(2) * t8 - qJD(4) * t11 - t934, -qJD(2) * t2 - qJD(4) * t5 - t931; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t792 - t807, pkin(2) * t740, 0, 0, t605, t595, 0, -t605, 0, 0, -qJD(4) * t508 + t621 * t792 - t844, -qJD(4) * t509 + t668, t690, t693, t329, t236, 0, -t327, 0, 0, -qJD(4) * t245 - qJD(5) * t298 + t670, -qJD(4) * t246 - qJD(5) * t299 + t669, t701 + t955, -qJD(4) * t68 - t702; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t605, t595, 0, -t605, 0, 0, -pkin(3) * t827, -pkin(3) * t826, 0, 0, t329, t236, 0, -t327, 0, 0, qJD(4) * t420 + t582 * t824, qJD(4) * t421 + t1007 * t824, 0, qJD(4) * t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, t402, t382, t339, t381, -t806, -pkin(8) * t826 - t646, pkin(8) * t827 - t647, 0, 0, t1050, t98, t218, t1051, t219, -t806, -t1018 - t651, t1020 - t650, -t652 + t845, (-t492 * t968 - t617 * t735) * t952 + t654; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1050, t98, t218, t1051, t219, -t806, -t1018 - t645, t1020 - t644, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, t666 + t808, -t444 * t801 + t618 * t766, -t224, t441 * t801 + t621 * t766, t449, -qJD(2) * t185 - qJD(3) * t189 - t841, -qJD(2) * t183 - qJD(3) * t187 + t842, 0, 0, -t1023, -t41, t146, t1023, t147, t449, qJD(2) * t38 + qJD(3) * t44 + qJD(5) * t51 - t947, qJD(2) * t37 + qJD(3) * t43 + qJD(5) * t53 - t946, qJD(2) * t9 + qJD(3) * t11 - t933, qJD(2) * t3 + qJD(3) * t5 - t930; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t339, -t402, -t430, t340, t809, t806, qJD(3) * t508 + t676, qJD(3) * t509 + t677, 0, 0, t1051, -t98, t276, t1050, t275, t806, qJD(3) * t245 + t691, qJD(3) * t246 + t692, t655, qJD(3) * t68 - t703; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t339, -t402, -t430, t340, t809, t806, t646, t647, 0, 0, t1051, -t98, t276, t1050, t275, t806, t651, t650, t652, -t654; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t959, -pkin(4) * t737, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t800 * t959 + t948, t889 + (-t738 - t737) * pkin(4), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1023, -t41, t146, t1023, t147, t449, qJD(2) * t63 + qJD(3) * t73 - qJD(4) * t51 - t944, qJD(2) * t64 + qJD(3) * t74 - qJD(4) * t53 - t945, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1051, -t98, t276, t1050, t275, t806, qJD(3) * t298 + t679, qJD(3) * t299 + t680, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1051, -t98, t276, t1050, t275, t806, t645, t644, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t617 * t952 - t948, pkin(4) * t738 - t889, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t13;

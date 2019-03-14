% Calculate inertial parameters regressor of coriolis matrix for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRRPRR4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:36
% EndTime: 2019-03-08 22:15:04
% DurationCPUTime: 23.95s
% Computational Cost: add. (12344->724), mult. (27125->997), div. (0->0), fcn. (29802->10), ass. (0->545)
t810 = qJD(3) - qJD(5);
t623 = cos(qJ(3));
t998 = pkin(8) - pkin(9);
t586 = t998 * t623;
t618 = sin(qJ(5));
t622 = cos(qJ(5));
t619 = sin(qJ(3));
t752 = t998 * t619;
t451 = -t586 * t622 - t618 * t752;
t621 = cos(qJ(6));
t1015 = t451 * t621;
t617 = sin(qJ(6));
t925 = qJ(4) * t619;
t970 = pkin(3) + pkin(4);
t542 = t623 * t970 + pkin(2) + t925;
t604 = t618 * t623;
t876 = t622 * t619;
t563 = -t604 + t876;
t881 = t618 * t619;
t560 = t622 * t623 + t881;
t934 = t560 * pkin(5);
t724 = -t563 * pkin(10) + t934;
t631 = t542 + t724;
t171 = t617 * t631 - t1015;
t1026 = (t171 + t1015) * t563;
t1016 = t451 * t617;
t170 = -t621 * t631 - t1016;
t1025 = (t170 + t1016) * t563;
t616 = sin(pkin(6));
t620 = sin(qJ(2));
t890 = t616 * t620;
t926 = cos(pkin(6));
t534 = t619 * t890 - t623 * t926;
t535 = t619 * t926 + t623 * t890;
t215 = t534 * t622 - t535 * t618;
t747 = t534 * t618 + t535 * t622;
t323 = t747 * t621;
t624 = cos(qJ(2));
t889 = t616 * t624;
t873 = t617 * t889 + t323;
t749 = t873 * t621;
t912 = t747 * t617;
t265 = -t621 * t889 + t912;
t918 = t265 * t617;
t1017 = (t747 - t749 - t918) * t215;
t917 = t1017 * qJD(1);
t961 = -t1015 / 0.2e1;
t802 = t621 * t890;
t457 = t560 * t889;
t898 = t457 * t617;
t367 = -t802 - t898;
t897 = t457 * t621;
t368 = -t617 * t890 + t897;
t939 = -t621 / 0.2e1;
t944 = t617 / 0.2e1;
t872 = t367 * t944 + t368 * t939;
t652 = t265 * t939 + t873 * t944;
t988 = t652 * t560;
t52 = t988 + t872;
t1024 = qJD(2) * t52;
t960 = -t451 / 0.2e1;
t1007 = t215 * t960;
t51 = t988 - t872;
t1023 = t51 * qJD(2);
t1022 = t52 * qJD(1);
t882 = t618 * t586;
t989 = t622 * t752;
t994 = -t989 + t882;
t1002 = t994 * t617;
t403 = t1002 / 0.2e1;
t404 = -t1002 / 0.2e1;
t1021 = t403 + t404;
t1001 = t994 * t621;
t406 = t1001 / 0.2e1;
t407 = -t1001 / 0.2e1;
t1020 = t406 + t407;
t935 = t451 * pkin(5);
t1008 = t215 * t617;
t776 = t1008 / 0.2e1;
t1014 = t451 * t994;
t614 = t621 ^ 2;
t945 = t614 / 0.2e1;
t612 = t617 ^ 2;
t946 = t612 / 0.2e1;
t753 = t945 + t946;
t1003 = t753 * t215;
t1013 = t1003 * t618;
t848 = qJD(3) * t215;
t1012 = qJD(5) * t215 - t848;
t551 = t563 ^ 2;
t341 = -t560 ^ 2 + t551;
t1005 = t341 * t621;
t1011 = qJD(2) * t1005;
t1006 = t341 * t617;
t1010 = qJD(2) * t1006;
t1009 = t810 * t451;
t942 = t618 / 0.2e1;
t310 = t215 * t942;
t866 = t612 + t614;
t748 = t866 * t215;
t405 = t994 * t942;
t999 = t341 * qJD(2);
t997 = -t170 / 0.2e1;
t964 = t747 / 0.2e1;
t954 = t989 / 0.2e1;
t984 = t810 * t563;
t996 = t560 * t984;
t571 = qJ(4) * t618 + t622 * t970;
t566 = pkin(5) + t571;
t913 = t747 * t566;
t613 = t619 ^ 2;
t615 = t623 ^ 2;
t995 = t613 + t615;
t993 = t810 * t747;
t992 = 0.2e1 * t563;
t592 = t614 - t612;
t991 = t592 * t810;
t694 = t170 * t621 - t171 * t617;
t987 = t694 * t560;
t572 = qJ(4) * t622 - t618 * t970;
t567 = -pkin(10) + t572;
t986 = t753 * t567;
t985 = t810 * t560;
t726 = t753 * pkin(10);
t947 = t572 / 0.2e1;
t953 = -t566 / 0.2e1;
t971 = -pkin(5) / 0.2e1;
t111 = (t571 * t753 + t953 + t971) * t618 + (t726 + t947 - t986) * t622;
t932 = t563 * pkin(5);
t933 = t560 * pkin(10);
t423 = t932 + t933;
t887 = t617 * t423;
t213 = t887 - t1001;
t920 = t213 * t621;
t879 = t621 * t423;
t212 = t879 + t1002;
t921 = t212 * t617;
t672 = t920 / 0.2e1 - t921 / 0.2e1;
t902 = t451 * t622;
t772 = t902 / 0.2e1;
t936 = t622 / 0.2e1;
t763 = t617 * t936;
t924 = t171 * t621;
t783 = t924 / 0.2e1;
t800 = t170 * t763 + t622 * t783 + t405;
t15 = t618 * t672 + t772 + t800;
t911 = t747 * t622;
t774 = -t911 / 0.2e1;
t719 = t749 / 0.2e1;
t799 = t265 * t763 + t622 * t719 - t310;
t24 = t774 + t799 + t1013;
t455 = (-0.1e1 + t866) * t622 * t618;
t446 = t455 * qJD(4);
t983 = -t24 * qJD(1) - t15 * qJD(2) + t111 * qJD(3) - t446;
t610 = t623 * qJ(4);
t548 = -t619 * t970 + t610;
t328 = -t423 + t548;
t888 = t617 * t328;
t173 = t1001 + t888;
t922 = t173 * t621;
t880 = t621 * t328;
t172 = t880 - t1002;
t923 = t172 * t617;
t673 = t922 / 0.2e1 - t923 / 0.2e1;
t948 = t571 / 0.2e1;
t755 = t948 + t953;
t951 = -t567 / 0.2e1;
t756 = t951 + t947;
t981 = t560 * t755 + t563 * t756 - t934 / 0.2e1;
t562 = t866 * t622;
t819 = t562 * qJD(3);
t980 = -qJD(5) * t562 + t819;
t543 = t876 / 0.2e1 - t604 / 0.2e1;
t856 = qJD(2) * t563;
t790 = t560 * t856;
t979 = qJD(6) * t543 + t790;
t383 = t621 * t560;
t353 = t617 * t383;
t754 = t946 - t614 / 0.2e1;
t376 = t754 * t563;
t826 = t376 * qJD(6);
t978 = -qJD(5) * t353 - t826;
t977 = -t614 * t790 + t826;
t976 = t612 * t790 + t826;
t975 = -qJD(3) * t353 + t826;
t269 = -t560 * t622 + t563 * t618;
t587 = t616 ^ 2 * t620 * t624;
t638 = -t587 + (t534 * t619 + t535 * t623) * t889;
t974 = qJD(2) * t638;
t972 = t638 * qJD(1);
t852 = qJD(2) * t621;
t796 = t617 * t852;
t132 = t376 * t810 - t551 * t796;
t969 = t172 / 0.2e1;
t968 = -t173 / 0.2e1;
t967 = -t213 / 0.2e1;
t966 = -t265 / 0.2e1;
t965 = t265 / 0.2e1;
t963 = -t747 / 0.2e1;
t962 = t215 / 0.2e1;
t959 = -t994 / 0.2e1;
t958 = t994 / 0.2e1;
t957 = t451 / 0.2e1;
t955 = t563 / 0.2e1;
t952 = t566 / 0.2e1;
t950 = t567 / 0.2e1;
t949 = -t571 / 0.2e1;
t943 = -t618 / 0.2e1;
t941 = -t619 / 0.2e1;
t940 = t619 / 0.2e1;
t938 = t621 / 0.2e1;
t937 = -t622 / 0.2e1;
t931 = t618 * pkin(5);
t930 = t619 * pkin(3);
t929 = t15 * qJD(5);
t758 = -t215 / 0.2e1 + t962;
t723 = t758 * t560;
t759 = t963 + t964;
t55 = t563 * t759 - t723;
t928 = t55 * qJD(3);
t927 = t24 * qJD(5);
t919 = t265 * t560;
t801 = t619 * t889;
t569 = t622 * t801;
t456 = t604 * t889 - t569;
t914 = t215 * t456;
t146 = t215 * t621;
t905 = t994 * t456;
t901 = t456 * t617;
t900 = t456 * t621;
t899 = t456 * t622;
t46 = -t265 * t367 + t368 * t873 - t914;
t896 = t46 * qJD(1);
t894 = t55 * qJD(2);
t893 = t566 * t618;
t892 = t571 * t618;
t891 = t572 * t622;
t379 = t617 * t563;
t884 = t617 * t571;
t384 = t621 * t563;
t96 = t457 * t747 - t587 - t914;
t874 = t96 * qJD(1);
t764 = t889 / 0.2e1;
t727 = t619 * t764;
t732 = t623 * t764;
t871 = t618 * t727 + t622 * t732;
t765 = -t889 / 0.2e1;
t733 = t623 * t765;
t870 = t622 * t733 + t765 * t881;
t869 = t569 / 0.2e1 + t618 * t733;
t868 = -t569 / 0.2e1 + t618 * t732;
t867 = t995 * pkin(8) * t889;
t664 = t560 * t942 + t563 * t936;
t649 = t940 + t664;
t279 = t649 * t617;
t865 = qJD(2) * t279;
t282 = t649 * t621;
t864 = qJD(2) * t282;
t859 = qJD(2) * t376;
t391 = t592 * t551;
t858 = qJD(2) * t391;
t857 = qJD(2) * t560;
t855 = qJD(2) * t617;
t854 = qJD(2) * t619;
t853 = qJD(2) * t620;
t851 = qJD(2) * t623;
t850 = qJD(3) * qJ(4);
t847 = qJD(3) * t617;
t846 = qJD(3) * t621;
t319 = t911 / 0.2e1;
t150 = t319 + t774;
t845 = qJD(4) * t150;
t410 = -t902 / 0.2e1;
t217 = t410 + t772;
t844 = qJD(4) * t217;
t843 = qJD(4) * t562;
t842 = qJD(4) * t618;
t841 = qJD(4) * t619;
t840 = qJD(4) * t622;
t838 = qJD(5) * t542;
t836 = qJD(5) * t617;
t835 = qJD(5) * t621;
t834 = qJD(6) * t617;
t833 = qJD(6) * t621;
t832 = qJD(6) * t622;
t363 = t866 * t619 * t563;
t827 = t363 * qJD(2);
t377 = t617 * t560;
t825 = t377 * qJD(2);
t824 = t383 * qJD(2);
t709 = -pkin(3) * t623 - t925;
t573 = -pkin(2) + t709;
t581 = -t610 + t930;
t441 = t573 * t623 + t581 * t619;
t823 = t441 * qJD(2);
t442 = -t573 * t619 + t581 * t623;
t822 = t442 * qJD(2);
t821 = t543 * qJD(2);
t818 = t592 * qJD(6);
t593 = t615 - t613;
t817 = t593 * qJD(2);
t816 = t593 * qJD(3);
t815 = t613 * qJD(2);
t814 = t618 * qJD(3);
t813 = t619 * qJD(3);
t812 = t622 * qJD(3);
t609 = t623 * qJD(3);
t811 = t623 * qJD(4);
t807 = pkin(2) * t854;
t806 = pkin(2) * t851;
t805 = pkin(8) * t813;
t804 = pkin(8) * t609;
t798 = t612 * t856;
t797 = t614 * t856;
t795 = t617 * t846;
t794 = t560 * t841;
t793 = t617 * t835;
t792 = t560 * t834;
t791 = t560 * t833;
t789 = t573 * t854;
t788 = t616 * t853;
t787 = qJD(2) * t889;
t599 = t617 * t833;
t786 = t560 * t854;
t785 = t563 * t854;
t784 = t560 * t852;
t779 = t747 * t939;
t778 = t747 * t937;
t775 = t146 / 0.2e1;
t773 = t451 * t937;
t771 = -t899 / 0.2e1;
t770 = -t898 / 0.2e1;
t769 = -t897 / 0.2e1;
t766 = t890 / 0.2e1;
t761 = t384 / 0.2e1;
t757 = t960 + t957;
t751 = t873 * t560;
t750 = t873 * t563;
t412 = t866 * t571;
t745 = t810 * t617;
t743 = t810 * t621;
t578 = t810 * t622;
t742 = qJD(6) + t857;
t738 = t617 * t786;
t737 = t619 * t784;
t736 = t551 * t599;
t735 = t563 * t599;
t734 = t563 * t796;
t731 = t560 * t765;
t730 = t560 * t764;
t729 = t563 * t765;
t728 = t563 * t764;
t725 = t971 + t755;
t718 = -0.2e1 * t734;
t717 = 0.2e1 * t734;
t712 = t621 * t745;
t710 = t617 * t743;
t671 = t963 * t994 + t1007;
t627 = t172 * t965 + t873 * t968 - t671;
t630 = t456 * t952 - t567 * t872;
t674 = -t924 / 0.2e1 + t617 * t997;
t1 = -t215 * t674 + t627 + t630;
t18 = -t170 * t172 + t171 * t173 + t1014;
t708 = -t1 * qJD(1) + t18 * qJD(2);
t707 = t560 * t734;
t13 = (t172 * t621 + t173 * t617) * t563 - t987;
t706 = -qJD(2) * t13 - t1022;
t625 = t548 * t764 + t747 * t958 - t1007 + t671;
t668 = t456 * t948 + t457 * t947;
t16 = -t625 + t668;
t98 = t542 * t548;
t705 = -qJD(1) * t16 + qJD(2) * t98;
t19 = (t212 * t621 + t213 * t617) * t563 + t987;
t704 = -qJD(2) * t19 + t1022;
t27 = (t967 + t968) * t621 + (t212 / 0.2e1 + t969) * t617;
t703 = qJD(2) * t27;
t428 = -t900 / 0.2e1;
t645 = -t747 * t955 - t723;
t628 = t265 * t955 + t617 * t645;
t35 = t428 + t628;
t42 = t1025 + (t172 + t1002) * t560;
t702 = t35 * qJD(1) + t42 * qJD(2);
t425 = t901 / 0.2e1;
t626 = t645 * t621 + t750 / 0.2e1;
t36 = t425 + t626;
t43 = t1026 + (-t173 + t1001) * t560;
t701 = t36 * qJD(1) + t43 * qJD(2);
t427 = t900 / 0.2e1;
t675 = (t966 + t912 / 0.2e1) * t563;
t38 = t427 + t675;
t44 = -t1025 + (t212 - t1002) * t560;
t700 = t38 * qJD(1) + t44 * qJD(2);
t426 = -t901 / 0.2e1;
t642 = -t750 / 0.2e1 + t747 * t761;
t41 = t426 + t642;
t45 = -t1026 + (-t213 - t1001) * t560;
t699 = t41 * qJD(1) + t45 * qJD(2);
t47 = t771 + (t265 * t940 + t368 * t942) * t621 + (t367 * t943 + t873 * t941) * t617;
t73 = t694 * t619;
t698 = qJD(1) * t47 + qJD(2) * t73;
t697 = t55 * qJD(1);
t665 = -t215 * t955 + t766;
t84 = t769 - t919 / 0.2e1 + t665 * t617;
t99 = t170 * t560 - t379 * t994;
t696 = qJD(1) * t84 - qJD(2) * t99;
t695 = t24 * qJD(4) - t917;
t693 = -t922 + t923;
t692 = t920 - t921;
t691 = t560 * t566 + t563 * t567;
t100 = -t171 * t560 + t384 * t994;
t83 = t770 + t751 / 0.2e1 - t665 * t621;
t690 = qJD(1) * t83 - qJD(2) * t100;
t138 = t759 * t617;
t629 = (-pkin(10) / 0.2e1 + t756) * t563 + (pkin(5) / 0.2e1 + t755) * t560;
t62 = t617 * t757 + t621 * t629;
t689 = -qJD(1) * t138 + qJD(2) * t62;
t139 = -t779 - t323 / 0.2e1;
t59 = t961 + t1015 / 0.2e1 + t629 * t617;
t688 = -qJD(1) * t139 + qJD(2) * t59;
t140 = t758 * t617;
t666 = t560 * t950 + t563 * t953;
t650 = t328 / 0.2e1 + t666;
t77 = -t621 * t650 + t1021;
t687 = -qJD(1) * t140 - qJD(2) * t77;
t141 = t758 * t621;
t75 = t617 * t650 + t1020;
t686 = -qJD(1) * t141 - qJD(2) * t75;
t685 = t742 * t621;
t684 = qJD(1) * t150 + qJD(2) * t217;
t248 = -t542 * t563 + t548 * t560;
t287 = t728 + t868;
t683 = qJD(1) * t287 - qJD(2) * t248;
t249 = t542 * t560 + t548 * t563;
t289 = t731 + t871;
t682 = qJD(1) * t289 - qJD(2) * t249;
t444 = t954 - t989 / 0.2e1;
t681 = qJD(2) * t444;
t680 = qJD(3) * t566 - t840;
t679 = qJD(3) * t572 + t842;
t678 = qJD(5) * t572 + t842;
t677 = t933 / 0.2e1 + t932 / 0.2e1;
t676 = t610 / 0.2e1 - t930 / 0.2e1;
t669 = -t451 * t952 + t947 * t994;
t667 = t899 / 0.2e1 + t457 * t943;
t163 = (t581 / 0.2e1 + t676) * t889;
t663 = -t573 * t581 * qJD(2) + t163 * qJD(1);
t209 = t727 + t667;
t662 = -qJD(1) * t209 - t542 * t854;
t288 = t729 + t869;
t661 = qJD(1) * t288 - t542 * t856;
t290 = t730 + t870;
t660 = qJD(1) * t290 + t542 * t857;
t658 = t423 / 0.2e1 + t677;
t657 = t712 * t992;
t656 = -qJD(3) * t534 + t623 * t787;
t655 = qJD(3) * t535 + t619 * t787;
t169 = -t412 * t567 + t566 * t572;
t6 = t935 / 0.2e1 + (pkin(10) * t968 + t171 * t949 + t213 * t950) * t621 + (pkin(10) * t969 + t170 * t949 + t212 * t951) * t617 + t669;
t641 = pkin(5) * t964 - pkin(10) * t1003;
t651 = t719 + t918 / 0.2e1;
t9 = -t913 / 0.2e1 + t651 * t571 - (-t572 / 0.2e1 + t986) * t215 + t641;
t654 = -t9 * qJD(1) + t6 * qJD(2) + t169 * qJD(3);
t20 = -t170 * t212 + t171 * t213 - t1014;
t632 = t212 * t965 + t747 * t959 + t873 * t967;
t633 = -pkin(10) * t872 + t456 * t971;
t3 = -(t957 - t674) * t215 + t632 + t633;
t653 = -t3 * qJD(1) + t20 * qJD(2) + t15 * qJD(4);
t11 = (t960 + t674) * t622 + (t959 + t673) * t618;
t23 = (t964 - t651) * t622 + (t962 - t1003) * t618;
t266 = t562 * t567 + t893;
t648 = qJD(1) * t23 + qJD(2) * t11 - qJD(3) * t266;
t118 = t757 * t622 + (t958 + t959) * t618;
t413 = t891 + t892;
t67 = t618 * t758 - t622 * t759;
t646 = -qJD(1) * t67 - qJD(2) * t118 + qJD(3) * t413;
t598 = t617 * t814;
t644 = t618 * t836 - t621 * t832 - t598;
t602 = t621 * t814;
t643 = t617 * t832 + t618 * t835 - t602;
t107 = t617 * t658 + t1020;
t370 = t725 * t621;
t640 = pkin(5) * t835 - qJD(2) * t107 + qJD(3) * t370;
t109 = -t621 * t658 + t1021;
t369 = t725 * t617;
t639 = pkin(5) * t836 - qJD(2) * t109 + qJD(3) * t369;
t635 = qJD(3) * t709 + t811;
t601 = t619 * t609;
t600 = t619 * t851;
t577 = t810 * t618;
t475 = t995 * t787;
t466 = -0.2e1 * t735;
t465 = 0.2e1 * t735;
t463 = (-t620 * t851 - t624 * t813) * t616;
t462 = (-t609 * t624 + t619 * t853) * t616;
t389 = t810 * t543;
t372 = 0.2e1 * t566 * t938;
t371 = t884 / 0.2e1 + (pkin(5) + t566) * t944;
t346 = -t882 + 0.2e1 * t954;
t344 = t718 - t991;
t343 = t717 + t991;
t294 = t729 + t868;
t293 = t728 + t869;
t292 = t731 + t870;
t291 = t730 + t871;
t281 = t619 * t938 - t621 * t664;
t280 = (t664 + t941) * t617;
t268 = -t712 - t859;
t267 = t710 + t859;
t227 = t269 * t621;
t226 = t269 * t617;
t214 = t217 * qJD(5);
t208 = t727 - t667;
t187 = (-t612 / 0.2e1 + t945 - t754) * t560;
t164 = t581 * t765 + t676 * t889;
t149 = t323 / 0.2e1 - t779;
t148 = -t912 / 0.2e1 - t747 * t944;
t147 = -t215 * t939 + t775;
t144 = 0.2e1 * t776;
t137 = t150 * qJD(5);
t119 = 0.2e1 * t405 + t410 + t773;
t112 = -t891 / 0.2e1 + t893 / 0.2e1 - t931 / 0.2e1 + t622 * t726 + t866 * (-t892 / 0.2e1 + t567 * t936);
t110 = 0.2e1 * t403 + t879 / 0.2e1 - t677 * t621;
t108 = 0.2e1 * t406 - t887 / 0.2e1 + t677 * t617;
t86 = -t751 / 0.2e1 - t215 * t761 + t770 - t802 / 0.2e1;
t85 = t919 / 0.2e1 + t563 * t776 + t769 + t617 * t766;
t78 = 0.2e1 * t404 + t880 / 0.2e1 - t666 * t621;
t76 = 0.2e1 * t407 - t888 / 0.2e1 + t666 * t617;
t70 = -0.2e1 * t310 + t319 - t778;
t69 = -0.2e1 * t1003;
t61 = pkin(10) * t761 + t621 * t981 + t1016;
t60 = 0.2e1 * t961 + pkin(10) * t379 / 0.2e1 + t981 * t617;
t48 = -t618 * t872 + t619 * t652 + t771;
t40 = t425 + t642;
t39 = t428 + t675;
t37 = t426 + t626;
t34 = t427 + t628;
t26 = -t672 + t673;
t25 = -t778 + t799 - t1013;
t17 = t625 + t668;
t12 = t618 * t673 + t773 + t800;
t10 = t749 * t949 + t884 * t966 + t913 / 0.2e1 - t215 * t947 + t641 - t951 * t748;
t5 = -t935 / 0.2e1 + t674 * t571 + t672 * t567 + t669 + t673 * pkin(10);
t4 = t170 * t776 + t171 * t775 - t1007 - t632 + t633;
t2 = t1008 * t997 - t215 * t783 - t627 + t630;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t974, 0, 0, 0, 0, 0, 0, 0, 0, 0, t974, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t46 + t1017 * t810; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t788, -t787, 0, 0, 0, 0, 0, 0, 0, 0, t463, t462, t475, t972 + (-pkin(2) * t890 + t867) * qJD(2), 0, 0, 0, 0, 0, 0, t463, t475, -t462, t972 + (t573 * t890 + t867) * qJD(2) + t164 * qJD(3) + qJD(4) * t801, 0, 0, 0, 0, 0, 0, qJD(3) * t294 + qJD(5) * t293 - t560 * t788, qJD(3) * t291 + qJD(5) * t292 - t563 * t788 (t456 * t563 - t457 * t560) * qJD(2) + t928, t874 + (-t451 * t457 - t542 * t890 + t905) * qJD(2) + t17 * qJD(3) + t208 * qJD(4), 0, 0, 0, 0, 0, 0 (t367 * t560 + t379 * t456) * qJD(2) + t34 * qJD(3) + t39 * qJD(5) + t86 * qJD(6) (-t368 * t560 + t384 * t456) * qJD(2) + t37 * qJD(3) + t40 * qJD(5) + t85 * qJD(6) (-t367 * t621 - t368 * t617) * t856 - t810 * t51, t896 + (-t170 * t367 + t171 * t368 + t905) * qJD(2) + t2 * qJD(3) + t48 * qJD(4) + t4 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t655, -t656, 0, 0, 0, 0, 0, 0, 0, 0, -t655, 0, t656, t164 * qJD(2) + (-pkin(3) * t535 - qJ(4) * t534) * qJD(3) + t535 * qJD(4), 0, 0, 0, 0, 0, 0, qJD(2) * t294 - t993, qJD(2) * t291 + t1012, t894, t17 * qJD(2) + (-t215 * t572 - t571 * t747) * qJD(3) + t70 * qJD(4), 0, 0, 0, 0, 0, 0, qJD(2) * t34 + qJD(5) * t149 + qJD(6) * t144 - t747 * t846, qJD(2) * t37 + qJD(5) * t148 + qJD(6) * t147 + t747 * t847, t69 * qJD(5) + t848 * t866 - t1023, t917 + t2 * qJD(2) + (-t567 * t748 - t913) * qJD(3) + t25 * qJD(4) + t10 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t655, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t208 + qJD(3) * t70 + t137, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t48 + qJD(3) * t25 + t927; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t293 + t993, qJD(2) * t292 - t1012, 0, t845, 0, 0, 0, 0, 0, 0, qJD(2) * t39 + qJD(3) * t149 - qJD(6) * t1008 - t747 * t835, qJD(2) * t40 + qJD(3) * t148 - qJD(6) * t146 + t747 * t836, t69 * qJD(3) + qJD(5) * t748 + t1023, t4 * qJD(2) + t10 * qJD(3) + (-pkin(5) * t747 + pkin(10) * t748) * qJD(5) + t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t86 + qJD(3) * t144 - qJD(5) * t1008 - qJD(6) * t873, qJD(2) * t85 + qJD(3) * t147 - qJD(5) * t146 + qJD(6) * t265, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t972, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t163 - t972, 0, 0, 0, 0, 0, 0, -qJD(3) * t287 - qJD(5) * t288, -qJD(3) * t289 - qJD(5) * t290, t928, -qJD(3) * t16 + qJD(4) * t209 - t874, 0, 0, 0, 0, 0, 0, qJD(3) * t35 + qJD(5) * t38 - qJD(6) * t83, qJD(3) * t36 + qJD(5) * t41 - qJD(6) * t84, -t810 * t52, -qJD(3) * t1 - qJD(4) * t47 - qJD(5) * t3 - t896; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t601, t816, 0, -t601, 0, 0, -pkin(2) * t813, -pkin(2) * t609, 0, 0, t601, 0, -t816, 0, 0, -t601, -qJD(3) * t442 + t619 * t811, 0, -qJD(3) * t441 + qJD(4) * t613 (qJD(3) * t581 - t841) * t573, t996, t810 * t341, 0, -t996, 0, 0, qJD(3) * t248 + t563 * t838 + t794, qJD(3) * t249 - t560 * t838 + t563 * t841, 0, qJD(3) * t98 + t542 * t841, t614 * t996 - t736, -qJD(6) * t391 - t560 * t657, -t1005 * t810 - t563 * t792, t612 * t996 + t736, t1006 * t810 - t563 * t791, -t996, qJD(3) * t42 + qJD(5) * t44 + qJD(6) * t100 + t621 * t794, qJD(3) * t43 + qJD(5) * t45 + qJD(6) * t99 - t617 * t794, -qJD(3) * t13 - qJD(4) * t363 - qJD(5) * t19, qJD(3) * t18 - qJD(4) * t73 + qJD(5) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t600, t817, t609, -t600, -t813, 0, -t804 - t807, t805 - t806, 0, 0, t600, t609, -t817, 0, t813, -t600, -t804 - t822, t635, -t805 - t823, pkin(8) * t635 - t663, t790, t999, -t985, -t790, -t984, 0, t1009 - t683, qJD(3) * t994 + qJD(5) * t346 - t682 (t571 * t560 + t572 * t563) * qJD(3) + t269 * qJD(4) + t697 (t451 * t571 + t572 * t994) * qJD(3) + t119 * qJD(4) + t705 (-t795 + t797) * t560 - t978, t187 * qJD(5) + t465 + (-qJD(3) * t592 + t718) * t560, -qJD(5) * t379 + t563 * t847 - t1011 (t795 + t798) * t560 + t978, -qJD(5) * t384 + t563 * t846 + t1010, -t979 (t617 * t691 + t1015) * qJD(3) + t226 * qJD(4) + t60 * qJD(5) + t78 * qJD(6) + t702 (t621 * t691 - t1016) * qJD(3) + t227 * qJD(4) + t61 * qJD(5) + t76 * qJD(6) + t701, qJD(3) * t693 + qJD(5) * t26 + t706 (t451 * t566 - t567 * t693) * qJD(3) + t12 * qJD(4) + t5 * qJD(5) + t708; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t600, t609, t815, -t789 + t804, 0, 0, 0, 0, 0, 0, t786, t785, qJD(3) * t269, qJD(3) * t119 + t214 - t662, 0, 0, 0, 0, 0, 0, qJD(3) * t226 + qJD(6) * t281 + t737, qJD(3) * t227 + qJD(6) * t280 - t738, -t827, qJD(3) * t12 - t698 + t929; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t790, -t999, t985, t790, t984, 0, -t1009 - t661, qJD(3) * t346 + qJD(5) * t994 - t660, 0, t844 (-t793 - t797) * t560 - t975, t187 * qJD(3) + t466 + (-qJD(5) * t592 + t717) * t560, -qJD(3) * t379 + t563 * t836 + t1011 (t793 - t798) * t560 + t975, -qJD(3) * t384 + t563 * t835 - t1010, t979, t60 * qJD(3) + (t617 * t724 + t1015) * qJD(5) + t110 * qJD(6) + t700, t61 * qJD(3) + (t621 * t724 - t1016) * qJD(5) + t108 * qJD(6) + t699, qJD(3) * t26 + qJD(5) * t692 + t704, t5 * qJD(3) + (pkin(10) * t692 + t935) * qJD(5) + t653; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t710 * t992 - t858, -t742 * t379, -t132, -t563 * t685, -t389, qJD(3) * t78 + qJD(4) * t281 + qJD(5) * t110 - qJD(6) * t171 - t690, qJD(3) * t76 + qJD(4) * t280 + qJD(5) * t108 + qJD(6) * t170 - t696, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t163, 0, 0, 0, 0, 0, 0, qJD(2) * t287, qJD(2) * t289, -t894, qJD(2) * t16 - qJD(4) * t67, 0, 0, 0, 0, 0, 0, -qJD(2) * t35 - qJD(5) * t139 + qJD(6) * t140, -qJD(2) * t36 - qJD(5) * t138 + qJD(6) * t141, t1024, qJD(2) * t1 - qJD(4) * t23 - qJD(5) * t9 - t917; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t600, -t817, 0, t600, 0, 0, t807, t806, 0, 0, -t600, 0, t817, 0, 0, t600, t822, 0, t823, t663, -t790, -t999, 0, t790, 0, 0, t683, qJD(5) * t444 + t682, -t697, -qJD(4) * t118 - t705, t977, t465 + 0.2e1 * t707, -t791 + t1011, -t976, t792 - t1010, t979, qJD(5) * t59 + qJD(6) * t77 - t702, qJD(5) * t62 + qJD(6) * t75 - t701, qJD(5) * t27 - t706, -qJD(4) * t11 + qJD(5) * t6 - t708; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), 0, 0, 0, 0, 0, 0, t678, -qJD(5) * t571 + t840, 0, qJD(4) * t413, t599, t818, 0, -t599, 0, 0, -t566 * t834 + t621 * t678, -t566 * t833 - t617 * t678, qJD(5) * t412 - t843, qJD(4) * t266 + qJD(5) * t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t850, 0, 0, 0, 0, 0, 0, t814, t812, 0, t646, 0, 0, 0, 0, 0, 0, t602, -t598, -t819, qJD(5) * t112 + t446 - t648; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t810 * t572, -t571 * t810 + t681, 0, 0, -t599, -t818, 0, t599, 0, 0, qJD(6) * t371 + t572 * t743 + t688, qJD(6) * t372 - t572 * t745 + t689, t412 * t810 + t703, t112 * qJD(4) + (-pkin(5) * t572 - pkin(10) * t412) * qJD(5) + t654; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t267, t343, -t685, t268, t742 * t617, t821, qJD(5) * t371 - t566 * t847 - t567 * t833 - t687, qJD(5) * t372 - t566 * t846 + t567 * t834 - t686, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t209 + qJD(3) * t67 + t137, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t47 + qJD(3) * t23 + t927; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t600, 0, -t815, t789, 0, 0, 0, 0, 0, 0, -t786, -t785, 0, qJD(3) * t118 + t214 + t662, 0, 0, 0, 0, 0, 0, -qJD(6) * t282 - t737, qJD(6) * t279 + t738, t827, qJD(3) * t11 + t698 + t929; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t850, 0, 0, 0, 0, 0, 0, -t577, -t578, 0, -t646, 0, 0, 0, 0, 0, 0, t643, -t644, t980, -qJD(5) * t111 + t648; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t455 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t577, t578, 0, t684, 0, 0, 0, 0, 0, 0, -t643, t644, -t980 (pkin(10) * t562 - t931) * qJD(5) - t983; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t578 * t617 - t618 * t833 - t864, t578 * t621 + t618 * t834 + t865, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t288, qJD(2) * t290, 0, -t845, 0, 0, 0, 0, 0, 0, -qJD(2) * t38 + qJD(3) * t139, -qJD(2) * t41 + qJD(3) * t138, -t1024, qJD(2) * t3 + qJD(3) * t9 - t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t790, t999, 0, -t790, 0, 0, t661, -qJD(3) * t444 + t660, 0, -t844, -t977, t466 - 0.2e1 * t707, qJD(6) * t383 - t1011, t976, -qJD(6) * t377 + t1010, -t979, -qJD(3) * t59 + qJD(6) * t109 - t700, -qJD(3) * t62 + qJD(6) * t107 - t699, -qJD(3) * t27 - t704, -qJD(3) * t6 - t653; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t679, qJD(3) * t571 - t681 - t840, 0, 0, -t599, -t818, 0, t599, 0, 0, -qJD(6) * t369 - t621 * t679 - t688, -qJD(6) * t370 + t617 * t679 - t689, -qJD(3) * t412 - t703 + t843, qJD(4) * t111 - t654; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t814, -t812, 0, -t684, 0, 0, 0, 0, 0, 0, -t602, t598, t819, t983; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t599, t818, 0, -t599, 0, 0, -pkin(5) * t834, -pkin(5) * t833, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, t344, t824 + t833, t267, -t825 - t834, -t821, -pkin(10) * t833 - t639, pkin(10) * t834 - t640, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t83 - qJD(3) * t140, qJD(2) * t84 - qJD(3) * t141, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t132, -t657 + t858, -qJD(5) * t383 + (t563 * t855 + t846) * t560, t132, qJD(5) * t377 + (t563 * t852 - t847) * t560, -t389, -qJD(3) * t77 + qJD(4) * t282 - qJD(5) * t109 + t690, -qJD(3) * t75 - qJD(4) * t279 - qJD(5) * t107 + t696, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, t344, t784, t267, -t560 * t855, -t821, qJD(5) * t369 + t617 * t680 + t687, qJD(5) * t370 + t621 * t680 + t686, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t617 * t812 + t864, -t621 * t812 - t865, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t267, t343, -t824, t268, t825, t821, t639, t640, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t7;
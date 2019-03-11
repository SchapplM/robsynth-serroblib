% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRR10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_invdynJ_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 06:17:10
% EndTime: 2019-03-10 06:19:51
% DurationCPUTime: 94.25s
% Computational Cost: add. (143123->1395), mult. (420672->1926), div. (0->0), fcn. (362042->18), ass. (0->583)
t608 = cos(qJ(2));
t912 = cos(pkin(6));
t856 = pkin(1) * t912;
t593 = t608 * t856;
t579 = qJD(1) * t593;
t605 = sin(qJ(2));
t601 = sin(pkin(6));
t911 = cos(pkin(7));
t736 = t601 * (-pkin(11) * t911 - pkin(10));
t718 = t605 * t736;
t464 = qJD(1) * t718 + t579;
t592 = t605 * t856;
t667 = t608 * t736 - t592;
t465 = t667 * qJD(1);
t909 = sin(pkin(7));
t821 = t608 * t909;
t709 = pkin(2) * t605 - pkin(11) * t821;
t872 = qJD(1) * t601;
t505 = t709 * t872;
t935 = cos(qJ(3));
t798 = t911 * t935;
t749 = qJD(3) * t798;
t932 = sin(qJ(3));
t794 = t909 * t932;
t797 = t911 * t932;
t991 = pkin(2) * t749 - t935 * t464 - t465 * t797 - t505 * t794;
t682 = -t605 * t797 + t608 * t935;
t500 = t682 * t601;
t491 = qJD(1) * t500;
t795 = t909 * t935;
t748 = qJD(3) * t795;
t707 = t748 - t491;
t910 = cos(pkin(8));
t777 = t910 * t909;
t715 = t932 * t777;
t766 = pkin(11) * t794;
t657 = -pkin(12) * t715 - t766;
t684 = t605 * t798 + t608 * t932;
t490 = t684 * t872;
t908 = sin(pkin(8));
t776 = t909 * t908;
t893 = t601 * t605;
t728 = t776 * t893;
t665 = qJD(1) * t728 - t490 * t910;
t990 = pkin(12) * t665 - qJD(3) * t657 - t991;
t305 = -t464 * t932 + t465 * t798 + t505 * t795;
t539 = pkin(2) * t797 + pkin(11) * t795;
t714 = t935 * t777;
t792 = t909 * t872;
t754 = t605 * t792;
t855 = pkin(12) * t910;
t989 = pkin(3) * t754 - t491 * t855 + t305 - (-pkin(12) * t714 - t539) * qJD(3);
t366 = -t465 * t909 + t911 * t505;
t710 = t935 * t776;
t854 = pkin(12) * t908;
t988 = t490 * pkin(3) - t491 * t854 + t366 - (pkin(3) * t794 - pkin(12) * t710) * qJD(3);
t604 = sin(qJ(4));
t662 = t911 * t908 + t714;
t934 = cos(qJ(4));
t452 = t604 * t794 - t662 * t934;
t883 = qJD(4) * t452 - t707 * t934 + (qJD(3) * t715 + t665) * t604;
t455 = t604 * t662 + t794 * t934;
t713 = t934 * t776;
t692 = t713 * t893;
t796 = t910 * t934;
t882 = qJD(1) * t692 - t490 * t796 - t491 * t604 + (t604 * t795 + t715 * t934) * qJD(3) + t455 * qJD(4);
t729 = t777 * t893;
t402 = -qJD(1) * t729 - t490 * t908;
t712 = t932 * t776;
t687 = qJD(3) * t712 + t402;
t936 = cos(qJ(1));
t800 = t912 * t936;
t933 = sin(qJ(1));
t532 = t605 * t933 - t608 * t800;
t533 = t605 * t800 + t608 * t933;
t759 = t601 * t794;
t381 = -t532 * t797 + t533 * t935 - t936 * t759;
t760 = t601 * t795;
t382 = t532 * t798 + t533 * t932 + t760 * t936;
t828 = t601 * t911;
t680 = t532 * t909 - t936 * t828;
t941 = t910 * t382 - t680 * t908;
t236 = -t381 * t934 + t604 * t941;
t315 = t382 * t908 + t680 * t910;
t603 = sin(qJ(5));
t607 = cos(qJ(5));
t163 = t236 * t607 - t315 * t603;
t602 = sin(qJ(6));
t987 = t163 * t602;
t606 = cos(qJ(6));
t986 = t163 * t606;
t443 = pkin(12) * t662 + t539;
t591 = pkin(2) * t798;
t462 = pkin(3) * t911 + t591 + t657;
t851 = t909 * pkin(2);
t497 = -pkin(3) * t795 - pkin(12) * t712 - t851;
t793 = t908 * t934;
t746 = qJD(4) * t793;
t750 = qJD(4) * t796;
t825 = t604 * t908;
t826 = t604 * t910;
t870 = qJD(4) * t604;
t951 = t443 * t870 - t462 * t750 - t497 * t746 + t988 * t825 + t989 * t826 + t934 * t990;
t881 = t908 * t989 - t910 * t988;
t985 = -pkin(13) * t687 + t951;
t984 = t882 * pkin(4) + pkin(13) * t883 + t881;
t892 = t601 * t608;
t875 = pkin(10) * t892 + t592;
t516 = t875 * qJD(1);
t817 = t912 * qJD(1);
t763 = t817 + qJD(2);
t697 = t763 * t909;
t803 = t608 * t828;
t410 = t516 + (qJD(1) * t803 + t697) * pkin(11);
t672 = pkin(2) * t912 + t718;
t413 = qJD(2) * pkin(2) + qJD(1) * t672 + t579;
t824 = t605 * t909;
t702 = -pkin(2) * t608 - pkin(11) * t824 - pkin(1);
t484 = t702 * t872;
t262 = -t410 * t932 + t413 * t798 + t484 * t795;
t853 = t935 * t605;
t681 = t608 * t797 + t853;
t671 = t681 * t601;
t677 = t932 * t697;
t416 = qJD(1) * t671 + t677;
t211 = -t416 * t855 + t262;
t653 = -t410 * t935 - t413 * t797 - t484 * t794;
t757 = t608 * t798;
t732 = t601 * t757;
t809 = t932 * t872;
t966 = -t605 * t809 + t935 * t697;
t415 = qJD(1) * t732 + t966;
t819 = t910 * t415;
t212 = -pkin(12) * t819 + t653;
t408 = t908 * t415;
t322 = t416 * pkin(3) - pkin(12) * t408;
t535 = pkin(3) * t796 - pkin(12) * t825;
t878 = t535 * qJD(4) - t934 * t211 - t212 * t826 - t322 * t825;
t321 = t415 * t934 - t416 * t826;
t705 = t746 - t321;
t560 = t608 * t792;
t686 = t763 * t911 - t560;
t675 = -qJD(3) - t686;
t658 = t908 * t675;
t647 = t819 - t658;
t285 = t934 * t416 + t604 * t647;
t812 = t912 * qJDD(1);
t752 = t812 + qJDD(2);
t695 = t752 * t909;
t862 = qJD(1) * qJD(2);
t785 = t911 * t862;
t745 = t605 * t785;
t965 = -t601 * t745 + t695;
t631 = qJD(3) * t966 + t932 * t965;
t802 = qJDD(1) * t853;
t938 = qJD(1) * (t935 * qJD(2) + t749) + qJDD(1) * t797;
t618 = -t631 - (t608 * t938 + t802) * t601;
t751 = qJD(3) * t797;
t861 = qJDD(1) * t605;
t841 = t601 * t861;
t847 = qJD(3) * t935;
t849 = t605 * t872;
t871 = qJD(2) * t608;
t651 = -t608 * t751 * t872 - qJD(3) * t677 - t809 * t871 - t932 * t841 - t847 * t849 + t935 * t965;
t630 = -qJDD(1) * t732 - t651;
t805 = t601 * t821;
t670 = -qJDD(1) * t805 + t752 * t911 + qJDD(3);
t806 = t601 * t824;
t753 = qJD(2) * t806;
t646 = qJD(1) * t753 + t670;
t968 = t630 * t910 - t646 * t908;
t612 = -t604 * t618 + t934 * t968;
t131 = t285 * qJD(4) + t612;
t610 = qJDD(5) + t131;
t835 = t416 * t908;
t983 = pkin(13) * t835 - t878;
t157 = -t212 * t908 + t910 * t322;
t320 = t415 * t604 + t416 * t796;
t982 = -t320 * pkin(4) + t321 * pkin(13) - t157 + (pkin(4) * t825 - pkin(13) * t793) * qJD(4);
t789 = qJD(4) * t825;
t791 = qJD(4) * t826;
t845 = qJD(4) * t934;
t950 = -t443 * t845 - t462 * t791 - t497 * t789 + t604 * t990 - t988 * t793 - t989 * t796;
t527 = -t910 * t911 + t710;
t867 = qJD(5) * t607;
t868 = qJD(5) * t603;
t887 = t455 * t868 + t527 * t867 - t603 * t687 + t607 * t883;
t376 = t455 * t607 - t527 * t603;
t886 = qJD(5) * t376 - t603 * t883 - t607 * t687;
t162 = t236 * t603 + t315 * t607;
t353 = -t462 * t908 + t910 * t497;
t258 = t452 * pkin(4) - t455 * pkin(13) + t353;
t299 = t934 * t443 + t462 * t826 + t497 * t825;
t267 = -pkin(13) * t527 + t299;
t958 = t258 * t867 - t267 * t868 + t984 * t603 - t607 * t985;
t952 = -pkin(4) * t687 - t950;
t530 = t603 * t825 - t607 * t910;
t827 = t603 * t908;
t880 = qJD(5) * t530 + t416 * t827 - t607 * t705;
t822 = t607 * t908;
t531 = t603 * t910 + t604 * t822;
t879 = qJD(5) * t531 + t416 * t822 + t603 * t705;
t538 = pkin(3) * t826 + pkin(12) * t793;
t877 = t538 * qJD(4) - t604 * t211 + t212 * t796 + t322 * t793;
t735 = t789 - t320;
t461 = t934 * t658;
t283 = -t415 * t796 + t416 * t604 + t461;
t738 = qJD(5) + t283;
t930 = pkin(13) * t607;
t981 = -pkin(14) * t882 - t958;
t980 = pkin(5) * t886 + pkin(14) * t887 + t952;
t512 = pkin(13) * t910 + t538;
t850 = t908 * pkin(3);
t513 = -pkin(4) * t793 - pkin(13) * t825 - t850;
t916 = -t512 * t868 + t513 * t867 + t982 * t603 - t607 * t983;
t172 = pkin(4) * t285 + pkin(13) * t283;
t197 = pkin(12) * t647 - t653;
t198 = -pkin(3) * t675 + t211;
t338 = -t413 * t909 + t911 * t484;
t228 = -t415 * pkin(3) - pkin(12) * t835 + t338;
t88 = -t604 * t197 + t198 * t796 + t228 * t793;
t71 = t603 * t172 + t607 * t88;
t979 = pkin(13) * t868 + t71;
t947 = pkin(4) * t835 + t877;
t598 = t601 ^ 2;
t977 = 0.2e1 * t598;
t239 = -t630 * t908 - t646 * t910 - qJDD(4);
t615 = pkin(12) * t618;
t811 = pkin(10) * t849;
t784 = qJD(2) * t817;
t761 = pkin(1) * t784;
t807 = pkin(1) * t812;
t860 = qJDD(1) * t608;
t840 = t601 * t860;
t857 = pkin(10) * t840 + t605 * t807 + t608 * t761;
t411 = -qJD(2) * t811 + t857;
t815 = qJDD(1) * t911;
t333 = (t695 + (t608 * t815 - t745) * t601) * pkin(11) + t411;
t577 = t608 * t807;
t700 = -t605 * t761 + t577;
t842 = t608 * t862;
t730 = -t842 - t861;
t701 = t730 * pkin(10);
t339 = t752 * pkin(2) + ((-t605 * t815 - t608 * t785) * pkin(11) + t701) * t601 + t700;
t694 = qJD(2) * t709;
t385 = (qJD(1) * t694 + qJDD(1) * t702) * t601;
t739 = -t333 * t932 + t339 * t798 + t385 * t795;
t747 = qJD(3) * t794;
t101 = pkin(3) * t646 - t410 * t847 - t413 * t751 - t484 * t747 + t615 * t910 + t739;
t254 = -t339 * t909 + t911 * t385;
t136 = pkin(3) * t630 + t615 * t908 + t254;
t846 = qJD(3) * t932;
t137 = t935 * t333 + t339 * t797 + t385 * t794 - t410 * t846 + t413 * t749 + t484 * t748;
t99 = -pkin(12) * t968 + t137;
t24 = t101 * t796 + t136 * t793 - t197 * t845 - t198 * t791 - t228 * t789 - t604 * t99;
t20 = pkin(4) * t239 - t24;
t130 = qJD(4) * t461 - t415 * t750 + t416 * t870 + t604 * t968 + t618 * t934;
t655 = -t675 * t910 - t408;
t652 = -qJD(4) - t655;
t336 = t607 * t652;
t76 = qJD(5) * t336 + t607 * t130 + t603 * t239 + t285 * t868;
t816 = -t603 * t130 + t607 * t239;
t192 = t607 * t285 - t603 * t652;
t869 = qJD(5) * t192;
t77 = t816 + t869;
t10 = pkin(5) * t77 + pkin(14) * t76 + t20;
t696 = -t101 * t826 - t136 * t825 + t197 * t870 - t198 * t750 - t228 * t746 - t934 * t99;
t19 = -pkin(13) * t239 - t696;
t68 = -t101 * t908 + t910 * t136;
t37 = t131 * pkin(4) + t130 * pkin(13) + t68;
t89 = t934 * t197 + (t198 * t910 + t228 * t908) * t604;
t84 = -pkin(13) * t652 + t89;
t134 = -t198 * t908 + t910 * t228;
t86 = t283 * pkin(4) - t285 * pkin(13) + t134;
t734 = -t607 * t19 - t603 * t37 + t84 * t868 - t86 * t867;
t5 = pkin(14) * t610 - t734;
t46 = t603 * t86 + t607 * t84;
t41 = pkin(14) * t738 + t46;
t190 = t285 * t603 + t336;
t83 = pkin(4) * t652 - t88;
t51 = t190 * pkin(5) - t192 * pkin(14) + t83;
t774 = t41 * t602 - t51 * t606;
t1 = -t774 * qJD(6) + t602 * t10 + t606 * t5;
t189 = qJD(6) + t190;
t975 = t774 * t189 + t1;
t949 = t603 * t258 + t607 * t267;
t957 = -qJD(5) * t949 + t603 * t985 + t984 * t607;
t974 = pkin(14) * t735 + t916;
t973 = pkin(5) * t879 + pkin(14) * t880 + t947;
t972 = t381 * t604;
t799 = t912 * t933;
t534 = -t605 * t799 + t608 * t936;
t685 = t605 * t936 + t608 * t799;
t664 = t685 * t911;
t627 = t534 * t932 + t664 * t935 - t760 * t933;
t967 = t685 * t909 + t933 * t828;
t969 = t627 * t908 + t910 * t967;
t45 = -t603 * t84 + t607 * t86;
t964 = -t45 * t738 - t734;
t499 = t684 * t601;
t414 = t499 * t908 + t729;
t963 = pkin(14) * t285 + t979;
t961 = -qJD(6) * t930 - t89 + t738 * (pkin(5) * t603 - pkin(14) * t607);
t16 = t41 * t606 + t51 * t602;
t2 = -qJD(6) * t16 + t606 * t10 - t602 * t5;
t960 = -t16 * t189 - t2;
t148 = pkin(14) * t452 + t949;
t297 = -t604 * t443 + t462 * t796 + t497 * t793;
t266 = t527 * pkin(4) - t297;
t375 = t455 * t603 + t607 * t527;
t173 = t375 * pkin(5) - t376 * pkin(14) + t266;
t90 = -t148 * t602 + t173 * t606;
t959 = qJD(6) * t90 + t602 * t980 - t606 * t981;
t91 = t148 * t606 + t173 * t602;
t928 = -qJD(6) * t91 + t602 * t981 + t606 * t980;
t927 = -pkin(5) * t882 - t957;
t948 = t607 * t512 + t603 * t513;
t915 = -qJD(5) * t948 + t603 * t983 + t982 * t607;
t779 = t912 * t909;
t444 = (t803 + t779) * pkin(11) + t875;
t463 = t593 + t672;
t876 = pkin(2) * t892 + pkin(11) * t806;
t498 = -pkin(1) * t601 - t876;
t298 = -t444 * t932 + t463 * t798 + t498 * t795;
t716 = t932 * t779;
t454 = t716 + t671;
t780 = t912 * t911;
t676 = t805 - t780;
t227 = -pkin(3) * t676 - t454 * t855 + t298;
t354 = -t463 * t909 + t911 * t498;
t683 = -t605 * t932 + t757;
t717 = t935 * t779;
t453 = -t601 * t683 - t717;
t834 = t454 * t908;
t742 = t453 * pkin(3) - pkin(12) * t834;
t257 = t742 + t354;
t146 = -t227 * t908 + t910 * t257;
t660 = t676 * t908;
t307 = t453 * t796 + t454 * t604 + t660 * t934;
t304 = t307 * pkin(4);
t649 = -t453 * t910 - t660;
t308 = t454 * t934 + t604 * t649;
t105 = -t308 * pkin(13) + t146 + t304;
t300 = t935 * t444 + t463 * t797 + t498 * t794;
t222 = pkin(12) * t649 + t300;
t112 = t934 * t222 + t227 * t826 + t257 * t825;
t370 = -t453 * t908 + t676 * t910;
t94 = -pkin(13) * t370 + t112;
t955 = t603 * t105 + t607 * t94;
t954 = t190 * t738;
t953 = t192 * t738;
t946 = -pkin(11) * t747 + t991;
t945 = -qJD(3) * t539 - t305;
t943 = t627 * t910 - t908 * t967;
t706 = t747 - t490;
t940 = (qJDD(2) + 0.2e1 * t812) * t601;
t141 = t606 * t192 + t602 * t738;
t866 = qJD(6) * t141;
t39 = -t602 * t76 - t606 * t610 + t866;
t522 = t875 * qJD(2);
t839 = t603 * t19 - t607 * t37;
t8 = -qJD(5) * t46 - t839;
t361 = qJD(3) * t716 + (qJD(2) * t684 + qJD(3) * t681) * t601;
t329 = -qJD(2) * t729 - t361 * t908;
t580 = qJD(2) * t593;
t468 = qJD(2) * t718 + t580;
t469 = t667 * qJD(2);
t506 = t601 * t694;
t183 = -t444 * t846 + t463 * t749 + t935 * t468 + t469 * t797 + t498 * t748 + t506 * t794;
t666 = qJD(2) * t728 - t361 * t910;
t152 = pkin(12) * t666 + t183;
t184 = -qJD(3) * t300 - t468 * t932 + t469 * t798 + t506 * t795;
t362 = qJD(3) * t717 + (qJD(2) * t682 + qJD(3) * t683) * t601;
t153 = pkin(3) * t753 - t362 * t855 + t184;
t367 = -t469 * t909 + t911 * t506;
t202 = t361 * pkin(3) - t362 * t854 + t367;
t47 = t934 * t152 + t153 * t826 + t202 * t825 - t222 * t870 + t227 * t750 + t257 * t746;
t43 = -pkin(13) * t329 + t47;
t110 = -t153 * t908 + t910 * t202;
t167 = -qJD(2) * t692 + qJD(4) * t308 + t361 * t796 + t362 * t604;
t168 = -qJD(4) * t307 + t362 * t934 + t604 * t666;
t59 = t167 * pkin(4) - t168 * pkin(13) + t110;
t14 = -qJD(5) * t955 - t43 * t603 + t59 * t607;
t937 = t870 * t908 - t320;
t609 = qJD(1) ^ 2;
t931 = pkin(10) * t601;
t708 = t606 * t738;
t865 = qJD(6) * t602;
t38 = -qJD(6) * t708 + t192 * t865 - t602 * t610 + t606 * t76;
t924 = t38 * t602;
t923 = t39 * t606;
t75 = qJDD(6) + t77;
t922 = t602 * t75;
t921 = t606 * t75;
t920 = t83 * t283;
t511 = -pkin(4) * t910 - t535;
t377 = t530 * pkin(5) - t531 * pkin(14) + t511;
t380 = -pkin(14) * t793 + t948;
t286 = t377 * t606 - t380 * t602;
t919 = qJD(6) * t286 + t602 * t973 + t606 * t974;
t287 = t377 * t602 + t380 * t606;
t918 = -qJD(6) * t287 - t602 * t974 + t606 * t973;
t917 = -pkin(5) * t735 - t915;
t571 = -pkin(5) * t607 - pkin(14) * t603 - pkin(4);
t864 = qJD(6) * t606;
t914 = t571 * t864 + t602 * t961 - t606 * t963;
t913 = -t571 * t865 + t602 * t963 + t606 * t961;
t139 = t192 * t602 - t708;
t907 = t139 * t189;
t906 = t139 * t602;
t905 = t141 * t139;
t904 = t189 * t141;
t903 = t192 * t190;
t233 = t934 * t941 + t972;
t902 = t233 * t603;
t901 = t233 * t607;
t384 = t534 * t935 - t664 * t932 + t759 * t933;
t237 = t384 * t604 + t934 * t943;
t900 = t237 * t603;
t899 = t237 * t607;
t898 = t283 * t285;
t897 = t307 * t603;
t896 = t307 * t607;
t895 = t415 * t416;
t894 = t598 * t609;
t891 = t602 * t607;
t890 = t606 * t607;
t889 = t376 * t865 - t452 * t864 - t602 * t882 + t606 * t887;
t312 = t376 * t606 + t452 * t602;
t888 = qJD(6) * t312 - t602 * t887 - t606 * t882;
t473 = t602 * t531 + t606 * t793;
t885 = -qJD(6) * t473 + t602 * t735 - t606 * t880;
t758 = t602 * t793;
t884 = -qJD(6) * t758 + t531 * t864 - t602 * t880 - t606 * t735;
t874 = t936 * pkin(1) + t933 * t931;
t599 = t605 ^ 2;
t600 = t608 ^ 2;
t873 = t599 - t600;
t858 = t608 * t894;
t848 = qJD(2) * t893;
t843 = pkin(1) * t977;
t838 = t381 * t908;
t836 = t384 * t908;
t830 = t533 * t909;
t829 = t534 * t909;
t814 = t189 * t606;
t813 = t738 * t603;
t810 = t605 * t858;
t808 = t605 * t842;
t801 = -pkin(1) * t933 + t936 * t931;
t790 = t601 * t609 * t912;
t238 = t384 * t934 - t604 * t943;
t164 = t238 * t603 - t607 * t969;
t787 = g(1) * t162 + g(2) * t164;
t235 = -t382 * t796 + t680 * t793 - t972;
t786 = g(1) * t235 + g(2) * t237;
t169 = -t283 * t891 - t606 * t285;
t783 = -t602 * t867 + t169;
t170 = -t283 * t890 + t285 * t602;
t782 = t606 * t867 - t170;
t775 = -t16 * t602 + t606 * t774;
t50 = pkin(14) * t307 + t955;
t223 = t308 * t603 + t370 * t607;
t224 = t308 * t607 - t370 * t603;
t111 = -t604 * t222 + t227 * t796 + t257 * t793;
t93 = t370 * pkin(4) - t111;
t69 = t223 * pkin(5) - t224 * pkin(14) + t93;
t22 = t50 * t606 + t602 * t69;
t21 = -t50 * t602 + t606 * t69;
t52 = t105 * t607 - t603 * t94;
t70 = t172 * t607 - t603 * t88;
t155 = t224 * t606 + t307 * t602;
t154 = t224 * t602 - t307 * t606;
t158 = t258 * t607 - t267 * t603;
t390 = -t603 * t512 + t607 * t513;
t762 = 0.2e1 * t817 + qJD(2);
t744 = -t382 * pkin(3) + pkin(12) * t838;
t743 = -t627 * pkin(3) + pkin(12) * t836;
t741 = -t532 * pkin(2) + pkin(11) * t830;
t740 = -t685 * pkin(2) + pkin(11) * t829;
t40 = -pkin(5) * t738 - t45;
t733 = -pkin(14) * t75 + t189 * t40;
t13 = t105 * t867 + t607 * t43 + t603 * t59 - t868 * t94;
t727 = g(1) * t936 + g(2) * t933;
t726 = g(1) * t164 - g(2) * t162 + g(3) * t223;
t165 = t238 * t607 + t603 * t969;
t725 = -g(1) * t165 + g(2) * t163 - g(3) * t224;
t404 = t532 * t932 - t533 * t798;
t405 = -t532 * t935 - t533 * t797;
t275 = t405 * t934 + (t404 * t910 + t533 * t776) * t604;
t347 = -t404 * t908 + t533 * t777;
t193 = t275 * t603 - t347 * t607;
t406 = -t534 * t798 + t685 * t932;
t407 = -t534 * t797 - t685 * t935;
t277 = t407 * t934 + (t406 * t910 + t534 * t776) * t604;
t348 = -t406 * t908 + t534 * t777;
t195 = t277 * t603 - t348 * t607;
t346 = t500 * t934 + (-t499 * t910 + t728) * t604;
t288 = t346 * t603 - t414 * t607;
t724 = -g(1) * t195 - g(2) * t193 - g(3) * t288;
t279 = -t381 * t826 - t382 * t934;
t203 = t279 * t603 - t381 * t822;
t281 = -t384 * t826 - t627 * t934;
t205 = t281 * t603 - t384 * t822;
t341 = -t453 * t934 - t454 * t826;
t291 = t341 * t603 - t454 * t822;
t723 = -g(1) * t205 - g(2) * t203 - g(3) * t291;
t722 = g(1) * t237 + g(2) * t233 + g(3) * t307;
t721 = g(1) * t238 - g(2) * t236 + g(3) * t308;
t274 = -t404 * t796 + t405 * t604 - t533 * t713;
t276 = -t406 * t796 + t407 * t604 - t534 * t713;
t345 = t499 * t796 + t500 * t604 - t692;
t720 = -g(1) * t276 - g(2) * t274 - g(3) * t345;
t278 = t381 * t796 - t382 * t604;
t280 = t384 * t796 - t604 * t627;
t340 = -t453 * t604 + t454 * t796;
t719 = -g(1) * t280 - g(2) * t278 - g(3) * t340;
t704 = t738 * qJD(5);
t6 = -pkin(5) * t610 - t8;
t703 = -t6 + t726;
t693 = t500 * pkin(3) + pkin(12) * t414 + t876;
t691 = t279 * pkin(4) + t278 * pkin(13) + t744;
t690 = t281 * pkin(4) + t280 * pkin(13) + t743;
t689 = t341 * pkin(4) + t340 * pkin(13) - t742;
t674 = pkin(14) * qJD(6) * t189 - t703;
t673 = -t533 * pkin(2) - pkin(11) * t680 + t801;
t668 = t346 * pkin(4) + t345 * pkin(13) + t693;
t659 = t675 * t909;
t656 = qJD(3) * t659;
t654 = t283 * t738 + t704;
t48 = -t604 * t152 + t153 * t796 + t202 * t793 - t222 * t845 - t227 * t791 - t257 * t789;
t650 = -t381 * pkin(3) - pkin(12) * t315 + t673;
t648 = t534 * pkin(2) + pkin(11) * t967 + t874;
t44 = t329 * pkin(4) - t48;
t645 = t405 * pkin(3) + pkin(12) * t347 + t741;
t644 = t407 * pkin(3) + pkin(12) * t348 + t740;
t642 = t652 * t908;
t641 = qJD(4) * t642;
t639 = t236 * pkin(4) + t235 * pkin(13) + t650;
t638 = t646 * t909;
t636 = t275 * pkin(4) + t274 * pkin(13) + t645;
t635 = t277 * pkin(4) + t276 * pkin(13) + t644;
t624 = t630 * t909;
t619 = t384 * pkin(3) + pkin(12) * t969 + t648;
t614 = t238 * pkin(4) + t237 * pkin(13) + t619;
t613 = t618 * t909;
t537 = -pkin(10) * t893 + t593;
t536 = -t766 + t591;
t521 = -pkin(10) * t848 + t580;
t514 = t579 - t811;
t508 = pkin(13) * t890 + t571 * t602;
t507 = -pkin(13) * t891 + t571 * t606;
t474 = t606 * t531 - t758;
t412 = t601 * t701 + t700;
t379 = pkin(5) * t793 - t390;
t311 = t376 * t602 - t606 * t452;
t292 = t341 * t607 + t454 * t827;
t289 = t346 * t607 + t414 * t603;
t231 = t237 * pkin(4);
t229 = t233 * pkin(4);
t206 = t281 * t607 + t384 * t827;
t204 = t279 * t607 + t381 * t827;
t196 = t277 * t607 + t348 * t603;
t194 = t275 * t607 + t347 * t603;
t147 = -pkin(5) * t452 - t158;
t138 = qJD(3) * t653 + t739;
t126 = pkin(5) * t192 + pkin(14) * t190;
t123 = t165 * t606 + t237 * t602;
t122 = -t165 * t602 + t237 * t606;
t104 = -qJD(5) * t223 + t168 * t607 - t329 * t603;
t103 = qJD(5) * t224 + t168 * t603 + t329 * t607;
t62 = -pkin(5) * t285 - t70;
t55 = -qJD(6) * t154 + t104 * t606 + t167 * t602;
t54 = qJD(6) * t155 + t104 * t602 - t167 * t606;
t49 = -pkin(5) * t307 - t52;
t34 = t126 * t602 + t45 * t606;
t33 = t126 * t606 - t45 * t602;
t17 = t103 * pkin(5) - t104 * pkin(14) + t44;
t12 = -pkin(5) * t167 - t14;
t11 = pkin(14) * t167 + t13;
t4 = -qJD(6) * t22 - t11 * t602 + t17 * t606;
t3 = qJD(6) * t21 + t11 * t606 + t17 * t602;
t7 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t933 - g(2) * t936, t727, 0, 0 (qJDD(1) * t599 + 0.2e1 * t808) * t598 (t605 * t860 - t862 * t873) * t977, t601 * t762 * t871 + t605 * t940 (qJDD(1) * t600 - 0.2e1 * t808) * t598, t608 * t940 - t762 * t848, t752 * t912, -t522 * t763 + t537 * t752 + t412 * t912 + g(1) * t533 - g(2) * t534 + (-t605 * t862 + t860) * t843, -g(1) * t532 + g(2) * t685 - t411 * t912 - t521 * t763 + t730 * t843 - t752 * t875 ((-t514 * qJD(2) + qJDD(1) * t875 + t411 + (-qJD(2) * t537 + t521) * qJD(1)) * t608 + (-t516 * qJD(2) - qJDD(1) * t537 - t412) * t605 - t727) * t601, t598 * qJDD(1) * pkin(1) ^ 2 - g(1) * t801 - g(2) * t874 + t411 * t875 + t412 * t537 - t514 * t522 + t516 * t521, t416 * t362 - t454 * t618, -t361 * t416 - t453 * t631 + t651 * t454 + t415 * t362 + (-t453 * t802 + (qJDD(1) * t454 * t798 - t453 * t938) * t608) * t601, -t362 * t675 + t416 * t753 + t454 * t646 + t618 * t676, -t415 * t361 + t453 * t630, t361 * t675 + t415 * t753 - t453 * t646 + t630 * t676, -t646 * t676 - t659 * t848, g(1) * t381 - g(2) * t384 - t138 * t676 - t184 * t675 + t254 * t453 + t262 * t753 + t298 * t646 + t338 * t361 + t354 * t630 - t367 * t415, -g(1) * t382 + g(2) * t627 + t137 * t676 + t183 * t675 + t254 * t454 - t300 * t646 + t338 * t362 - t354 * t618 + t367 * t416 + t653 * t753, g(1) * t680 - g(2) * t967 - t137 * t453 - t138 * t454 + t183 * t415 - t184 * t416 - t262 * t362 + t298 * t618 - t300 * t630 + t361 * t653, -g(1) * t673 - g(2) * t648 + t137 * t300 + t138 * t298 - t183 * t653 + t262 * t184 + t254 * t354 + t338 * t367, -t130 * t308 + t168 * t285, t130 * t307 - t131 * t308 - t167 * t285 - t168 * t283, t370 * t130 - t168 * t652 - t239 * t308 - t329 * t285, t131 * t307 + t167 * t283, t370 * t131 + t167 * t652 + t239 * t307 + t329 * t283, t239 * t370 + t329 * t652, -g(1) * t236 - g(2) * t238 + t110 * t283 - t111 * t239 + t146 * t131 + t134 * t167 - t24 * t370 + t68 * t307 - t88 * t329 - t48 * t652, t110 * t285 + t112 * t239 - t146 * t130 + t134 * t168 + t68 * t308 + t89 * t329 - t370 * t696 + t47 * t652 + t786, g(1) * t315 - g(2) * t969 + t111 * t130 - t112 * t131 - t89 * t167 - t88 * t168 - t24 * t308 - t47 * t283 - t48 * t285 + t307 * t696, -g(1) * t650 - g(2) * t619 + t134 * t110 + t24 * t111 - t112 * t696 + t68 * t146 + t89 * t47 + t88 * t48, t104 * t192 - t224 * t76, -t103 * t192 - t104 * t190 + t223 * t76 - t224 * t77, t104 * t738 + t167 * t192 + t224 * t610 - t307 * t76, t103 * t190 + t223 * t77, -t103 * t738 - t167 * t190 - t223 * t610 - t307 * t77, t167 * t738 + t307 * t610, -g(1) * t163 - g(2) * t165 + t83 * t103 + t14 * t738 + t45 * t167 + t44 * t190 + t20 * t223 + t8 * t307 + t52 * t610 + t93 * t77, t83 * t104 - t13 * t738 - t46 * t167 + t44 * t192 + t20 * t224 + t307 * t734 - t610 * t955 - t93 * t76 + t787, -t103 * t46 - t104 * t45 - t13 * t190 - t14 * t192 + t223 * t734 - t224 * t8 + t52 * t76 - t77 * t955 - t786, -g(1) * t639 - g(2) * t614 + t46 * t13 + t45 * t14 + t20 * t93 + t83 * t44 + t8 * t52 - t734 * t955, t141 * t55 - t155 * t38, -t139 * t55 - t141 * t54 + t154 * t38 - t155 * t39, t103 * t141 + t155 * t75 + t189 * t55 - t223 * t38, t139 * t54 + t154 * t39, -t103 * t139 - t154 * t75 - t189 * t54 - t223 * t39, t103 * t189 + t223 * t75, t4 * t189 + t21 * t75 + t2 * t223 - t774 * t103 + t12 * t139 + t49 * t39 + t6 * t154 + t40 * t54 - g(1) * (t235 * t602 + t986) - g(2) * t123, -t3 * t189 - t22 * t75 - t1 * t223 - t16 * t103 + t12 * t141 - t49 * t38 + t6 * t155 + t40 * t55 - g(1) * (t235 * t606 - t987) - g(2) * t122, -t1 * t154 - t139 * t3 - t141 * t4 - t155 * t2 - t16 * t54 + t21 * t38 - t22 * t39 + t55 * t774 - t787, t1 * t22 + t16 * t3 + t2 * t21 - t774 * t4 + t6 * t49 + t40 * t12 - g(1) * (t163 * pkin(5) + t162 * pkin(14) + t639) - g(2) * (t165 * pkin(5) + t164 * pkin(14) + t614); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t810, t873 * t894, -t608 * t790 + t841, t810, t605 * t790 + t840, t752, t577 + t516 * t763 + g(1) * t685 + g(2) * t532 + (-t784 + t894) * t605 * pkin(1) + (-g(3) * t608 + t701) * t601, t514 * t763 + pkin(1) * t858 + g(1) * t534 + g(2) * t533 + (pkin(10) * t862 + g(3)) * t893 - t857, 0, 0, t416 * t707 - t613 * t932, t415 * t707 - t416 * t706 - t613 * t935 - t624 * t932, -t416 * t754 + t491 * t675 - t618 * t911 + t638 * t932 - t656 * t935, -t415 * t706 - t624 * t935, -t415 * t754 - t490 * t675 - t630 * t911 + t638 * t935 + t656 * t932, t670 * t911 - (qJD(1) * t780 + qJD(3) - t560) * t754, -g(1) * t407 - g(2) * t405 - g(3) * t500 + t138 * t911 - t254 * t795 - t262 * t754 + t338 * t706 + t366 * t415 + t536 * t646 - t630 * t851 - t675 * t945, pkin(2) * t613 - g(1) * t406 - g(2) * t404 + g(3) * t499 - t137 * t911 + t254 * t794 + t338 * t707 - t366 * t416 - t539 * t646 - t653 * t754 + t675 * t946, -g(1) * t829 - g(2) * t830 - g(3) * t806 + t137 * t795 - t138 * t794 - t262 * t707 + t415 * t946 - t416 * t945 + t536 * t618 - t539 * t630 + t653 * t706, -g(1) * t740 - g(2) * t741 - g(3) * t876 + t137 * t539 + t138 * t536 - t254 * t851 + t262 * t945 - t338 * t366 - t653 * t946, -t130 * t455 - t285 * t883, t130 * t452 - t131 * t455 + t283 * t883 - t285 * t882, t527 * t130 - t239 * t455 + t285 * t687 + t652 * t883, t131 * t452 + t283 * t882, t527 * t131 + t239 * t452 - t283 * t687 + t652 * t882, t239 * t527 - t402 * t652 - t642 * t747, -g(1) * t277 - g(2) * t275 - g(3) * t346 + t353 * t131 + t134 * t882 - t297 * t239 - t24 * t527 + t283 * t881 + t68 * t452 - t652 * t950 + t687 * t88, -t353 * t130 - t134 * t883 + t299 * t239 + t285 * t881 + t68 * t455 - t527 * t696 - t652 * t951 - t687 * t89 - t720, -g(1) * t348 - g(2) * t347 - g(3) * t414 + t130 * t297 - t131 * t299 - t24 * t455 + t283 * t951 - t285 * t950 + t452 * t696 + t883 * t88 - t882 * t89, -g(1) * t644 - g(2) * t645 - g(3) * t693 + t134 * t881 + t24 * t297 - t299 * t696 + t68 * t353 + t88 * t950 - t89 * t951, -t192 * t887 - t376 * t76, t190 * t887 - t192 * t886 + t375 * t76 - t376 * t77, t192 * t882 + t376 * t610 - t452 * t76 - t738 * t887, t190 * t886 + t375 * t77, -t190 * t882 - t375 * t610 - t452 * t77 - t738 * t886, t452 * t610 + t738 * t882, -g(1) * t196 - g(2) * t194 - g(3) * t289 + t158 * t610 + t190 * t952 + t20 * t375 + t266 * t77 + t45 * t882 + t8 * t452 + t738 * t957 + t83 * t886, t192 * t952 + t20 * t376 - t266 * t76 + t452 * t734 - t46 * t882 - t610 * t949 - t738 * t958 - t83 * t887 - t724, t158 * t76 - t190 * t958 - t192 * t957 + t375 * t734 - t376 * t8 + t887 * t45 - t886 * t46 - t77 * t949 + t720, -g(1) * t635 - g(2) * t636 - g(3) * t668 + t8 * t158 + t20 * t266 + t45 * t957 + t46 * t958 - t734 * t949 + t83 * t952, -t141 * t889 - t312 * t38, t139 * t889 - t141 * t888 + t311 * t38 - t312 * t39, t141 * t886 - t189 * t889 + t312 * t75 - t375 * t38, t139 * t888 + t311 * t39, -t139 * t886 - t189 * t888 - t311 * t75 - t375 * t39, t189 * t886 + t375 * t75, t90 * t75 + t2 * t375 + t147 * t39 + t6 * t311 - g(1) * (t196 * t606 + t276 * t602) - g(2) * (t194 * t606 + t274 * t602) - g(3) * (t289 * t606 + t345 * t602) + t888 * t40 + t928 * t189 - t886 * t774 + t927 * t139, -t91 * t75 - t1 * t375 - t147 * t38 + t6 * t312 - g(1) * (-t196 * t602 + t276 * t606) - g(2) * (-t194 * t602 + t274 * t606) - g(3) * (-t289 * t602 + t345 * t606) - t889 * t40 - t959 * t189 - t886 * t16 + t927 * t141, -t1 * t311 - t139 * t959 - t141 * t928 - t16 * t888 - t2 * t312 + t38 * t90 - t39 * t91 - t774 * t889 + t724, t1 * t91 + t2 * t90 + t6 * t147 - g(1) * (t196 * pkin(5) + t195 * pkin(14) + t635) - g(2) * (t194 * pkin(5) + t193 * pkin(14) + t636) - g(3) * (t289 * pkin(5) + t288 * pkin(14) + t668) + t927 * t40 + t959 * t16 - t928 * t774; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t895, -t415 ^ 2 + t416 ^ 2, t415 * t675 - t618, t895, -t416 * t675 - t630, t646, g(1) * t627 + g(2) * t382 + g(3) * t453 - t338 * t416 - t653 * t686 + t739, g(1) * t384 + g(2) * t381 + g(3) * t454 - t262 * t675 - t338 * t415 - t137, 0, 0, -t130 * t825 + t285 * t705, -t130 * t793 - t131 * t825 + t283 * t321 + t320 * t285 + (-t283 * t793 - t285 * t825) * qJD(4), -t130 * t910 - t239 * t825 - t285 * t835 + t321 * t652 - t641 * t934, -t131 * t793 + t283 * t735, -t131 * t910 - t239 * t793 + t283 * t835 - t320 * t652 + t604 * t641, -t239 * t910 + t416 * t642, -g(1) * t281 - g(2) * t279 - g(3) * t341 - t131 * t850 + t134 * t735 - t157 * t283 - t535 * t239 + t24 * t910 + t652 * t877 - t68 * t793 - t835 * t88, t130 * t850 + t134 * t705 - t157 * t285 + t538 * t239 + t652 * t878 + t68 * t825 + t696 * t910 + t835 * t89 - t719, -t696 * t793 - g(1) * t836 - g(2) * t838 - g(3) * t834 - t24 * t825 + t535 * t130 - t538 * t131 + t89 * t320 + t88 * t321 + t877 * t285 - t878 * t283 + (-t793 * t88 - t825 * t89) * qJD(4), -g(1) * t743 - g(2) * t744 + g(3) * t742 - t134 * t157 + t24 * t535 - t538 * t696 - t68 * t850 - t877 * t88 + t878 * t89, -t192 * t880 - t531 * t76, t190 * t880 - t192 * t879 + t530 * t76 - t531 * t77, t937 * t192 + t531 * t610 - t880 * t738 + t76 * t793, t190 * t879 + t530 * t77, -t937 * t190 - t530 * t610 - t879 * t738 + t77 * t793, -t610 * t793 + t735 * t738, -g(1) * t206 - g(2) * t204 - g(3) * t292 + t190 * t947 + t20 * t530 + t390 * t610 + t45 * t735 + t511 * t77 + t738 * t915 - t793 * t8 + t83 * t879, t192 * t947 + t20 * t531 - t46 * t735 - t511 * t76 - t610 * t948 - t734 * t793 - t738 * t916 - t83 * t880 - t723, -t190 * t916 - t192 * t915 + t390 * t76 + t45 * t880 - t46 * t879 + t530 * t734 - t531 * t8 - t77 * t948 + t719, -g(1) * t690 - g(2) * t691 - g(3) * t689 + t20 * t511 + t8 * t390 + t915 * t45 + t916 * t46 - t734 * t948 + t83 * t947, t141 * t885 - t38 * t474, -t139 * t885 - t141 * t884 + t38 * t473 - t39 * t474, t141 * t879 + t189 * t885 - t38 * t530 + t474 * t75, t139 * t884 + t39 * t473, -t139 * t879 - t189 * t884 - t39 * t530 - t473 * t75, t189 * t879 + t530 * t75, t286 * t75 + t2 * t530 + t379 * t39 + t6 * t473 - g(1) * (t206 * t606 + t280 * t602) - g(2) * (t204 * t606 + t278 * t602) - g(3) * (t292 * t606 + t340 * t602) + t884 * t40 + t918 * t189 - t879 * t774 + t917 * t139, -t287 * t75 - t1 * t530 - t379 * t38 + t6 * t474 - g(1) * (-t206 * t602 + t280 * t606) - g(2) * (-t204 * t602 + t278 * t606) - g(3) * (-t292 * t602 + t340 * t606) + t885 * t40 - t919 * t189 - t879 * t16 + t917 * t141, -t1 * t473 - t139 * t919 - t141 * t918 - t16 * t884 - t2 * t474 + t286 * t38 - t287 * t39 + t774 * t885 + t723, t1 * t287 + t2 * t286 + t6 * t379 - g(1) * (t206 * pkin(5) + t205 * pkin(14) + t690) - g(2) * (t204 * pkin(5) + t203 * pkin(14) + t691) - g(3) * (t292 * pkin(5) + t291 * pkin(14) + t689) + t917 * t40 + t919 * t16 - t918 * t774; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t898, -t283 ^ 2 + t285 ^ 2, -t283 * t652 - t130, -t898, t655 * t285 - t612, -t239, -t134 * t285 - t652 * t89 + t24 + t722, t134 * t283 - t652 * t88 + t696 + t721, 0, 0, -t603 * t76 + t607 * t953 (-t76 - t954) * t607 + (-t77 - t953) * t603, -t285 * t192 + t603 * t610 + t607 * t654, t190 * t813 - t607 * t77, t285 * t190 - t603 * t654 + t607 * t610, -t738 * t285, -pkin(4) * t77 + g(1) * t899 + g(2) * t901 + g(3) * t896 - t89 * t190 - t20 * t607 - t45 * t285 - t70 * t738 - t704 * t930 + t83 * t868 + (-pkin(13) * t610 + t920) * t603, pkin(4) * t76 - g(1) * t900 - g(2) * t902 - g(3) * t897 - t89 * t192 + t20 * t603 + t46 * t285 + t607 * t920 - t610 * t930 + t738 * t979 + t83 * t867, t190 * t71 + t192 * t70 + ((-t77 + t869) * pkin(13) + t964) * t607 + (-t8 - t738 * t46 + (qJD(5) * t190 - t76) * pkin(13)) * t603 - t721, -t20 * pkin(4) + g(1) * t231 + g(2) * t229 + g(3) * t304 - t45 * t70 - t46 * t71 - t83 * t89 + (-t8 * t603 - t734 * t607 + (-t45 * t607 - t46 * t603) * qJD(5) - t721) * pkin(13), -t38 * t603 * t606 + (-t603 * t865 + t782) * t141, t139 * t170 + t141 * t169 + (-t139 * t606 - t141 * t602) * t867 + (t924 - t923 + (-t141 * t606 + t906) * qJD(6)) * t603, t38 * t607 + t782 * t189 + (t141 * t738 - t189 * t865 + t921) * t603, t39 * t602 * t603 + (t603 * t864 - t783) * t139, t39 * t607 + t783 * t189 + (-t139 * t738 - t189 * t864 - t922) * t603, t189 * t813 - t607 * t75, -t62 * t139 - t40 * t169 + t507 * t75 + t913 * t189 - t721 * t602 + (-t2 + (pkin(13) * t139 + t40 * t602) * qJD(5) + t722 * t606) * t607 + (pkin(13) * t39 + t40 * t864 + t6 * t602 - t738 * t774) * t603, -t62 * t141 - t40 * t170 - t508 * t75 - t914 * t189 - t721 * t606 + (t1 + (pkin(13) * t141 + t40 * t606) * qJD(5) - t722 * t602) * t607 + (-pkin(13) * t38 - t16 * t738 - t40 * t865 + t6 * t606) * t603, -t774 * t170 + t16 * t169 + t38 * t507 - t39 * t508 - t913 * t141 - t914 * t139 + t775 * t867 + (-t1 * t602 - t2 * t606 + (-t16 * t606 - t602 * t774) * qJD(6) + t722) * t603, t1 * t508 + t2 * t507 - t40 * t62 - g(1) * (-pkin(5) * t899 - pkin(14) * t900 - t231) - g(2) * (-pkin(5) * t901 - pkin(14) * t902 - t229) - g(3) * (-pkin(5) * t896 - pkin(14) * t897 - t304) + t914 * t16 - t913 * t774 + (t40 * t867 + t6 * t603 - t721) * pkin(13); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t903, -t190 ^ 2 + t192 ^ 2, -t76 + t954, -t903, t192 * t283 - t816, t610, -t83 * t192 + t283 * t46 + t726 - t839, t83 * t190 - t725 - t964, 0, 0, t141 * t814 - t924 (-t38 - t907) * t606 + (-t39 - t904) * t602, -t141 * t192 + t189 * t814 + t922, t189 * t906 - t923, -t189 ^ 2 * t602 + t139 * t192 + t921, -t189 * t192, -pkin(5) * t39 - t139 * t46 - t189 * t33 + t192 * t774 + t602 * t733 - t606 * t674, pkin(5) * t38 - t141 * t46 + t16 * t192 + t189 * t34 + t602 * t674 + t606 * t733, t139 * t34 + t141 * t33 + ((-t39 + t866) * pkin(14) + t975) * t606 + ((qJD(6) * t139 - t38) * pkin(14) + t960) * t602 + t725, t774 * t33 - t16 * t34 - t40 * t46 + t703 * pkin(5) + (qJD(6) * t775 + t1 * t606 - t2 * t602 + t725) * pkin(14); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t905, -t139 ^ 2 + t141 ^ 2, -t38 + t907, -t905, t904 - t39, t75, -t40 * t141 - g(1) * t122 - g(2) * (t233 * t606 + t987) + g(3) * t154 - t960, t40 * t139 + g(1) * t123 - g(2) * (-t233 * t602 + t986) + g(3) * t155 - t975, 0, 0;];
tau_reg  = t7;

% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:47
% EndTime: 2019-03-09 02:48:01
% DurationCPUTime: 66.48s
% Computational Cost: add. (33758->1213), mult. (51439->1548), div. (0->0), fcn. (52426->10), ass. (0->552)
t967 = Icges(5,1) + Icges(6,1);
t966 = Icges(5,4) - Icges(6,5);
t965 = Icges(6,4) + Icges(5,5);
t979 = Icges(5,2) + Icges(6,3);
t974 = Icges(6,2) + Icges(5,3);
t962 = Icges(5,6) - Icges(6,6);
t531 = pkin(9) + qJ(3);
t518 = cos(t531);
t532 = sin(pkin(10));
t540 = cos(qJ(1));
t796 = t540 * t532;
t534 = cos(pkin(10));
t538 = sin(qJ(1));
t797 = t538 * t534;
t440 = t518 * t796 - t797;
t798 = t534 * t540;
t799 = t532 * t538;
t441 = t518 * t798 + t799;
t438 = t518 * t799 + t798;
t439 = t518 * t797 - t796;
t517 = sin(t531);
t803 = t517 * t538;
t922 = t438 * t979 - t439 * t966 - t803 * t962;
t970 = -t966 * t438 + t439 * t967 + t965 * t803;
t982 = -t440 * t922 - t441 * t970;
t920 = -t438 * t962 + t439 * t965 + t803 * t974;
t978 = t532 * t979 - t534 * t966;
t977 = -t532 * t962 + t534 * t965;
t976 = -t966 * t532 + t534 * t967;
t957 = t438 * t922 + t439 * t970;
t935 = t532 * t922 + t534 * t970;
t802 = t517 * t540;
t973 = -t802 * t920 + t982;
t921 = t440 * t979 - t441 * t966 - t802 * t962;
t919 = -t440 * t962 + t441 * t965 + t802 * t974;
t918 = -t440 * t966 + t441 * t967 + t802 * t965;
t941 = t517 * t978 + t518 * t962;
t890 = t517 * t977 - t518 * t974;
t940 = -t517 * t976 + t518 * t965;
t537 = sin(qJ(6));
t539 = cos(qJ(6));
t277 = t438 * t539 - t439 * t537;
t278 = t438 * t537 + t439 * t539;
t124 = Icges(7,5) * t278 + Icges(7,6) * t277 - Icges(7,3) * t803;
t281 = t440 * t539 - t441 * t537;
t282 = t440 * t537 + t441 * t539;
t126 = Icges(7,5) * t282 + Icges(7,6) * t281 - Icges(7,3) * t802;
t829 = Icges(7,4) * t278;
t128 = -Icges(7,2) * t277 + Icges(7,6) * t803 - t829;
t266 = Icges(7,4) * t277;
t131 = -Icges(7,1) * t278 + Icges(7,5) * t803 - t266;
t794 = t128 * t281 + t131 * t282;
t828 = Icges(7,4) * t282;
t129 = Icges(7,2) * t281 - Icges(7,6) * t802 + t828;
t267 = Icges(7,4) * t281;
t132 = Icges(7,1) * t282 - Icges(7,5) * t802 + t267;
t795 = -t277 * t129 - t278 * t132;
t969 = t794 + t795 + (t124 * t540 + t126 * t538) * t517;
t793 = t281 * t129 + t282 * t132;
t908 = -t128 * t277 - t131 * t278;
t968 = -t908 + t517 * (t124 * t538 - t126 * t540) + t793;
t801 = t518 * t538;
t819 = Icges(4,3) * t540;
t378 = Icges(4,5) * t801 - Icges(4,6) * t803 - t819;
t489 = Icges(4,4) * t803;
t826 = Icges(4,5) * t540;
t382 = Icges(4,1) * t801 - t489 - t826;
t822 = Icges(4,6) * t540;
t380 = Icges(4,4) * t801 - Icges(4,2) * t803 - t822;
t808 = t380 * t517;
t619 = -t382 * t518 + t808;
t961 = -t378 * t540 - t538 * t619 + t803 * t920 + t957;
t507 = Icges(4,4) * t518;
t634 = -Icges(4,2) * t517 + t507;
t381 = Icges(4,6) * t538 + t540 * t634;
t831 = Icges(4,4) * t517;
t456 = Icges(4,1) * t518 - t831;
t383 = Icges(4,5) * t538 + t456 * t540;
t325 = t383 * t801;
t452 = Icges(4,5) * t518 - Icges(4,6) * t517;
t379 = Icges(4,3) * t538 + t452 * t540;
t678 = t379 * t540 - t325;
t156 = -t381 * t803 - t678;
t931 = t438 * t921 + t439 * t918 + t803 * t919;
t960 = t156 + t931;
t800 = t518 * t540;
t785 = t538 * t379 + t383 * t800;
t158 = -t381 * t802 + t785;
t930 = t440 * t921 + t441 * t918 + t802 * t919;
t959 = t158 + t930;
t891 = -t962 * t517 + t518 * t978;
t958 = t974 * t517 + t518 * t977;
t889 = t965 * t517 + t518 * t976;
t817 = t124 * t518;
t616 = t532 * t539 - t534 * t537;
t903 = t616 * t517;
t615 = t532 * t537 + t534 * t539;
t904 = t615 * t517;
t48 = -t128 * t903 - t131 * t904 + t817;
t954 = t532 * t941 - t534 * t940;
t936 = t532 * t921 + t534 * t918;
t522 = t540 * qJ(2);
t482 = pkin(1) * t538 - t522;
t535 = cos(pkin(9));
t510 = pkin(2) * t535 + pkin(1);
t536 = -pkin(7) - qJ(2);
t513 = t540 * t536;
t763 = -t538 * t510 - t513;
t373 = t482 + t763;
t363 = qJD(1) * t373;
t857 = pkin(3) * t518;
t459 = qJ(4) * t517 + t857;
t431 = t459 * t538;
t745 = qJD(3) * t540;
t700 = t518 * t745;
t469 = qJ(4) * t700;
t742 = qJD(4) * t540;
t480 = t517 * t742;
t463 = qJD(1) * t482;
t519 = qJD(2) * t538;
t758 = t519 - t463;
t764 = t480 + t519;
t953 = qJD(1) * t431 - t363 + t469 - t480 - t758 + t764;
t453 = Icges(4,2) * t518 + t831;
t952 = -t453 + t456 + t890;
t951 = t973 * t540;
t455 = Icges(4,1) * t517 + t507;
t950 = (t455 + t634 + t954 - t958) * qJD(1);
t617 = t453 * t517 - t455 * t518;
t451 = Icges(4,5) * t517 + Icges(4,6) * t518;
t804 = t451 * t540;
t949 = t438 * t941 - t439 * t940 - t538 * t617 + t803 * t890 - t804;
t805 = t451 * t538;
t948 = t440 * t941 - t441 * t940 - t540 * t617 + t802 * t890 + t805;
t133 = rSges(7,1) * t278 + rSges(7,2) * t277 - rSges(7,3) * t803;
t215 = rSges(7,1) * t904 + rSges(7,2) * t903 + rSges(7,3) * t518;
t739 = qJD(6) * t517;
t450 = t538 * t739 + t745;
t738 = qJD(6) * t518;
t497 = qJD(1) + t738;
t946 = -t133 * t497 - t215 * t450;
t701 = t517 * t745;
t317 = qJD(1) * t438 + t532 * t701;
t318 = -qJD(1) * t439 - t534 * t701;
t750 = qJD(1) * t538;
t705 = t517 * t750;
t580 = -t700 + t705;
t928 = -t317 * t979 - t318 * t966 + t580 * t962;
t746 = qJD(3) * t538;
t702 = t517 * t746;
t461 = t532 * t702;
t749 = qJD(1) * t540;
t704 = t518 * t749;
t319 = t532 * t704 - t534 * t750 - t461;
t320 = qJD(1) * t441 - t534 * t702;
t413 = t517 * t749 + t518 * t746;
t927 = -t319 * t979 + t320 * t966 + t413 * t962;
t926 = t317 * t962 + t318 * t965 - t580 * t974;
t925 = t319 * t962 - t320 * t965 - t413 * t974;
t924 = t966 * t317 + t318 * t967 - t965 * t580;
t923 = t966 * t319 - t320 * t967 - t965 * t413;
t944 = t891 * qJD(3);
t943 = t958 * qJD(3);
t942 = t889 * qJD(3);
t407 = qJD(5) * t440;
t169 = t318 * pkin(4) - qJ(5) * t317 + t407;
t287 = pkin(4) * t439 + t438 * qJ(5);
t939 = qJD(1) * t287 + t169 - t407 + t953;
t938 = -t540 * t890 - t936;
t937 = t538 * t890 + t935;
t846 = rSges(3,2) * sin(pkin(9));
t848 = rSges(3,1) * t535;
t613 = t538 * rSges(3,3) + (-t846 + t848) * t540;
t520 = qJD(2) * t540;
t521 = t538 * qJ(2);
t757 = t540 * pkin(1) + t521;
t909 = qJD(1) * t757 - t520;
t291 = qJD(1) * t613 + t909;
t208 = Icges(7,5) * t904 + Icges(7,6) * t903 + Icges(7,3) * t518;
t827 = Icges(7,4) * t904;
t210 = Icges(7,2) * t903 + Icges(7,6) * t518 + t827;
t375 = Icges(7,4) * t903;
t212 = Icges(7,1) * t904 + Icges(7,5) * t518 + t375;
t61 = -t208 * t803 + t210 * t277 + t212 * t278;
t929 = t497 * t61;
t854 = pkin(5) * t534;
t447 = pkin(8) * t518 + t517 * t854;
t858 = pkin(3) * t517;
t457 = -qJ(4) * t518 + t858;
t646 = pkin(4) * t534 + qJ(5) * t532;
t902 = t646 * t517;
t772 = -t902 - t457;
t709 = -t447 + t772;
t658 = t709 * t540;
t289 = t441 * pkin(4) + qJ(5) * t440;
t434 = pkin(3) * t800 + qJ(4) * t802;
t743 = qJD(4) * t538;
t478 = t517 * t743;
t491 = t540 * t510;
t676 = -t536 * t538 + t491;
t716 = qJD(1) * (t676 - t757) + t909;
t607 = qJD(1) * t434 - t457 * t746 + t478 + t716;
t741 = qJD(5) * t438;
t917 = qJD(1) * t289 - t746 * t902 + t607 + t741;
t786 = -t538 * t378 - t382 * t800;
t157 = -t380 * t802 - t786;
t916 = t157 - t973;
t915 = t948 * qJD(1);
t914 = (-t157 * t540 + t538 * t959 + t951) * qJD(3);
t913 = (t538 * t960 - t540 * t961) * qJD(3);
t912 = t949 * qJD(1);
t907 = 0.2e1 * qJD(3);
t747 = qJD(3) * t518;
t905 = t538 * t747;
t897 = t941 * t538;
t896 = t941 * t540;
t895 = t940 * t538;
t894 = t940 * t540;
t892 = -qJD(1) * t440 + t319 + t461;
t790 = t282 * rSges(7,1) + t281 * rSges(7,2);
t135 = -rSges(7,3) * t802 + t790;
t420 = t441 * pkin(5);
t341 = -pkin(8) * t802 + t420;
t449 = -t540 * t739 + t746;
t888 = -qJD(1) * t341 - t135 * t497 + t215 * t449 - t917;
t260 = t441 * rSges(6,1) + rSges(6,2) * t802 + t440 * rSges(6,3);
t887 = qJD(1) * t260 + t917;
t885 = t912 + t913;
t884 = t914 + t915;
t588 = qJD(3) * t453;
t253 = qJD(1) * t381 - t538 * t588;
t589 = qJD(3) * t455;
t255 = qJD(1) * t383 - t538 * t589;
t883 = (-t253 - t925) * t518 + (t532 * t927 + t534 * t923 - t255) * t517 + (-t517 * t920 - t518 * t935 + t619) * qJD(3);
t252 = -t540 * t588 + (-t538 * t634 + t822) * qJD(1);
t254 = -t540 * t589 + (-t456 * t538 + t826) * qJD(1);
t807 = t381 * t517;
t618 = -t383 * t518 + t807;
t882 = (t252 - t926) * t518 + (t532 * t928 + t534 * t924 + t254) * t517 + (t517 * t919 + t518 * t936 - t618) * qJD(3);
t443 = t634 * qJD(3);
t444 = t456 * qJD(3);
t547 = qJD(1) * t451 - t443 * t517 + t444 * t518 + (-t453 * t518 - t455 * t517) * qJD(3);
t873 = t617 * qJD(1) + t452 * qJD(3);
t881 = -t317 * t941 - t318 * t940 + t440 * t944 + t441 * t942 + t538 * t873 + t547 * t540 - t580 * t890 + t802 * t943;
t880 = t319 * t941 - t320 * t940 + t413 * t890 + t438 * t944 + t439 * t942 + t547 * t538 - t540 * t873 + t803 * t943;
t198 = t380 * t518 + t382 * t517;
t879 = t517 * t935 - t920 * t518 + t198;
t199 = t381 * t518 + t383 * t517;
t878 = t517 * t936 - t518 * t919 + t199;
t877 = ((t538 * t938 + t540 * t937) * qJD(3) - t950) * t517 + t952 * qJD(1) * t518;
t573 = t380 * t540 - t381 * t538;
t776 = -t453 * t540 + t383;
t777 = -Icges(4,2) * t801 + t382 - t489;
t876 = (-t538 * t776 + t540 * t777) * t517 + (t538 * t919 - t540 * t920 + t573) * t518;
t648 = rSges(6,1) * t534 + rSges(6,3) * t532;
t365 = -rSges(6,2) * t518 + t517 * t648;
t651 = rSges(5,1) * t534 - rSges(5,2) * t532;
t366 = -rSges(5,3) * t518 + t517 * t651;
t587 = qJD(3) * t451;
t875 = -t540 * t587 + (-t452 * t538 + t618 + t819) * qJD(1);
t754 = qJD(1) * t379;
t874 = qJD(1) * t619 - t538 * t587 + t754;
t872 = m(5) / 0.2e1;
t871 = m(6) / 0.2e1;
t870 = m(7) / 0.2e1;
t733 = qJD(3) * qJD(6);
t691 = t518 * t733;
t342 = qJD(1) * t449 - t538 * t691;
t869 = t342 / 0.2e1;
t343 = qJD(1) * t450 - t540 * t691;
t868 = t343 / 0.2e1;
t867 = -t449 / 0.2e1;
t866 = t449 / 0.2e1;
t865 = -t450 / 0.2e1;
t864 = t450 / 0.2e1;
t863 = -t497 / 0.2e1;
t862 = t497 / 0.2e1;
t861 = t538 / 0.2e1;
t860 = -t540 / 0.2e1;
t859 = -rSges(7,3) - pkin(8);
t855 = pkin(5) * t320;
t396 = t616 * t518;
t226 = qJD(3) * t396 - qJD(6) * t904;
t397 = t615 * t518;
t227 = qJD(3) * t397 + qJD(6) * t903;
t118 = -qJD(6) * t278 + t319 * t539 - t320 * t537;
t119 = qJD(6) * t277 + t319 * t537 + t320 * t539;
t66 = Icges(7,5) * t119 + Icges(7,6) * t118 - Icges(7,3) * t413;
t68 = Icges(7,4) * t119 + Icges(7,2) * t118 - Icges(7,6) * t413;
t70 = Icges(7,1) * t119 + Icges(7,4) * t118 - Icges(7,5) * t413;
t748 = qJD(3) * t517;
t8 = -t124 * t748 - t128 * t226 - t131 * t227 + t518 * t66 + t68 * t903 + t70 * t904;
t853 = t8 * t450;
t116 = -qJD(6) * t282 - t317 * t539 - t318 * t537;
t117 = qJD(6) * t281 - t317 * t537 + t318 * t539;
t65 = Icges(7,5) * t117 + Icges(7,6) * t116 + Icges(7,3) * t580;
t67 = Icges(7,4) * t117 + Icges(7,2) * t116 + Icges(7,6) * t580;
t69 = Icges(7,1) * t117 + Icges(7,4) * t116 + Icges(7,5) * t580;
t9 = -t126 * t748 + t129 * t226 + t132 * t227 + t518 * t65 + t67 * t903 + t69 * t904;
t852 = t9 * t449;
t849 = pkin(1) - t510;
t847 = rSges(4,1) * t518;
t845 = rSges(6,2) * t517;
t843 = rSges(5,3) * t517;
t123 = rSges(7,1) * t227 + rSges(7,2) * t226 - rSges(7,3) * t748;
t218 = -pkin(8) * t413 + t855;
t506 = qJD(4) * t517;
t735 = qJD(1) * qJD(3);
t695 = t538 * t735;
t736 = qJD(1) * qJD(2);
t512 = t540 * t736;
t734 = qJD(3) * qJD(4);
t693 = t518 * t734;
t710 = t457 * t695 + t540 * t693 + t512;
t638 = -qJD(5) * t317 + t695 * t902 + t710;
t170 = pkin(4) * t320 + t319 * qJ(5) + t741;
t765 = -pkin(3) * t702 + t478;
t223 = pkin(3) * t704 + qJ(4) * t413 + t765;
t503 = t536 * t750;
t783 = t503 - (-t540 * t849 - t521) * qJD(1) - t909;
t724 = -t223 + t783;
t672 = -t170 + t724;
t448 = -pkin(8) * t517 + t518 * t854;
t744 = qJD(4) * t518;
t391 = qJD(3) * t459 - t744;
t411 = t646 * t518;
t740 = qJD(5) * t532;
t476 = t517 * t740;
t787 = -qJD(3) * t411 - t391 - t476;
t718 = -t448 * qJD(3) + t787;
t647 = rSges(7,1) * t119 + rSges(7,2) * t118;
t72 = -rSges(7,3) * t413 + t647;
t12 = -t123 * t450 + t215 * t342 - t497 * t72 + (t133 * t739 + t540 * t718) * qJD(3) + (-t218 + (qJD(3) * t447 - t506) * t538 + t672) * qJD(1) + t638;
t841 = t12 * t538;
t788 = t318 * pkin(5) + pkin(8) * t705;
t217 = -pkin(8) * t700 + t788;
t581 = -t518 * t750 - t701;
t222 = pkin(3) * t581 - qJ(4) * t705 + t469 + t480;
t514 = qJ(2) * t749;
t759 = t514 + t519;
t770 = qJD(1) * (-pkin(1) * t750 + t759) + t538 * t736;
t717 = qJD(1) * (-t514 + (t538 * t849 - t513) * qJD(1)) + t770;
t609 = t538 * t693 + t717 + (t222 + t480) * qJD(1);
t575 = qJD(1) * t169 + qJD(5) * t319 + t609;
t726 = t117 * rSges(7,1) + t116 * rSges(7,2) + rSges(7,3) * t705;
t71 = -rSges(7,3) * t700 + t726;
t13 = qJD(1) * t217 - t123 * t449 - t215 * t343 + t497 * t71 + (qJD(1) * t658 - t135 * t739 + t538 * t718) * qJD(3) + t575;
t840 = t13 * t540;
t839 = t48 * t342;
t816 = t126 * t518;
t49 = t129 * t903 + t132 * t904 + t816;
t838 = t49 * t343;
t526 = t538 * rSges(4,3);
t649 = rSges(6,1) * t439 + rSges(6,3) * t438;
t256 = rSges(6,2) * t803 + t649;
t715 = -t365 + t772;
t659 = t540 * t715;
t779 = t373 - t482;
t713 = -t431 + t779;
t667 = -t287 + t713;
t711 = t407 + t764;
t74 = qJD(3) * t659 + (-t256 + t667) * qJD(1) + t711;
t835 = t74 * t365;
t834 = -rSges(6,2) - qJ(4);
t833 = -rSges(5,3) - qJ(4);
t392 = rSges(4,1) * t801 - rSges(4,2) * t803 - t540 * rSges(4,3);
t458 = rSges(4,1) * t517 + rSges(4,2) * t518;
t698 = t458 * t745;
t657 = t519 - t698;
t173 = (-t392 + t779) * qJD(1) + t657;
t813 = t173 * t538;
t393 = rSges(4,1) * t800 - rSges(4,2) * t802 + t526;
t703 = t458 * t746;
t174 = qJD(1) * t393 - t703 + t716;
t433 = t458 * t540;
t812 = t174 * t433;
t811 = t208 * t518;
t806 = t447 * t538;
t339 = pkin(5) * t439 - pkin(8) * t803;
t792 = t133 + t339;
t261 = t441 * rSges(5,1) - t440 * rSges(5,2) + rSges(5,3) * t802;
t791 = -t261 - t434;
t789 = -t289 - t434;
t368 = t518 * t651 + t843;
t784 = -t368 * qJD(3) - t391;
t782 = -t366 - t457;
t781 = -t368 - t459;
t371 = t538 * t902;
t429 = t457 * t538;
t780 = t371 + t429;
t437 = t457 * t750;
t778 = t750 * t902 + t437;
t774 = t538 * t431 + t540 * t434;
t432 = t457 * t540;
t773 = -qJD(1) * t432 + t518 * t743;
t771 = -t411 - t459;
t766 = rSges(4,2) * t705 + rSges(4,3) * t749;
t504 = t538 * t846;
t762 = rSges(3,3) * t749 + qJD(1) * t504;
t761 = t503 + t520;
t760 = t540 * rSges(3,3) + t504;
t752 = qJD(1) * t452;
t732 = t538 * t848;
t725 = t540 * t222 + t538 * t223 + t431 * t749;
t723 = -t260 + t789;
t722 = -t341 + t789;
t721 = t318 * rSges(5,1) + t317 * rSges(5,2) + rSges(5,3) * t700;
t720 = t318 * rSges(6,1) + rSges(6,2) * t700 - t317 * rSges(6,3);
t367 = t518 * t648 + t845;
t719 = -t367 * qJD(3) + t787;
t714 = -t367 + t771;
t712 = -t429 * t746 - t432 * t745 + t506;
t708 = -t448 + t771;
t706 = t491 + t434;
t477 = t518 * t740;
t697 = -t802 / 0.2e1;
t696 = -pkin(1) - t848;
t694 = t540 * t735;
t692 = t517 * t733;
t690 = t750 / 0.2e1;
t688 = -t748 / 0.2e1;
t687 = -t746 / 0.2e1;
t685 = -t745 / 0.2e1;
t684 = t745 / 0.2e1;
t682 = -t510 - t857;
t92 = qJD(1) * t261 - t366 * t746 + t607;
t681 = t92 * t782;
t679 = t540 * t782;
t677 = -t378 + t807;
t675 = pkin(3) * t701;
t42 = qJD(3) * t658 + (-t339 + t667) * qJD(1) + t711 + t946;
t674 = t42 * (-qJ(4) - t859);
t673 = -t123 + t718;
t671 = t222 * t745 + t223 * t746 + t431 * t694 + t517 * t734;
t670 = -t215 + t709;
t668 = t538 * t287 + t540 * t289 + t774;
t666 = t761 - t765;
t75 = -t365 * t746 + t887;
t664 = t75 * t715;
t660 = qJD(6) * t688;
t656 = t431 * t746 + t434 * t745 - t744;
t654 = -rSges(4,2) * t517 + t847;
t653 = rSges(5,1) * t320 - rSges(5,2) * t319;
t652 = rSges(5,1) * t439 - rSges(5,2) * t438;
t650 = rSges(6,1) * t320 + rSges(6,3) * t319;
t44 = -t124 * t803 + t908;
t45 = -t126 * t803 - t795;
t645 = -t44 * t540 + t45 * t538;
t644 = t44 * t538 + t45 * t540;
t46 = -t124 * t802 - t794;
t47 = -t126 * t802 + t793;
t643 = -t46 * t540 + t47 * t538;
t642 = t46 * t538 + t47 * t540;
t641 = -t48 * t540 + t49 * t538;
t640 = t48 * t538 + t49 * t540;
t629 = -t133 * t540 + t135 * t538;
t628 = -t173 * t540 - t174 * t538;
t614 = -t287 + t763;
t120 = Icges(7,5) * t227 + Icges(7,6) * t226 - Icges(7,3) * t748;
t121 = Icges(7,4) * t227 + Icges(7,2) * t226 - Icges(7,6) * t748;
t122 = Icges(7,1) * t227 + Icges(7,4) * t226 - Icges(7,5) * t748;
t19 = t120 * t518 + t121 * t903 + t122 * t904 - t208 * t748 + t210 * t226 + t212 * t227;
t78 = t210 * t903 + t212 * t904 + t811;
t612 = t19 * t497 - t692 * t78;
t481 = t518 * t742;
t611 = -t476 * t540 + t481;
t610 = t540 * t169 + t538 * t170 + t287 * t749 + t725;
t372 = t540 * t902;
t608 = -t371 * t746 - t372 * t745 + t477 + t712;
t430 = t458 * t538;
t594 = t289 + t706;
t214 = (t392 * t538 + t393 * t540) * qJD(3);
t583 = -t459 - t845;
t582 = -t459 - t843;
t579 = (Icges(7,5) * t277 - Icges(7,6) * t278) * t450 - (Icges(7,5) * t281 - Icges(7,6) * t282) * t449 - (Icges(7,5) * t903 - Icges(7,6) * t904) * t497;
t578 = -qJD(1) * t372 - t476 * t538 + t773;
t577 = t287 * t746 + t289 * t745 + t476 + t656;
t576 = qJD(3) * t477 + t169 * t745 + t170 * t746 + t287 * t694 + t671;
t257 = rSges(5,3) * t803 + t652;
t572 = -t170 + t666;
t571 = t517 * t579;
t559 = (Icges(7,1) * t281 - t129 - t828) * t449 - (Icges(7,1) * t277 + t128 - t829) * t450 + (Icges(7,1) * t903 - t210 - t827) * t497;
t558 = (-Icges(7,2) * t282 + t132 + t267) * t449 - (-Icges(7,2) * t278 - t131 + t266) * t450 + (-Icges(7,2) * t904 + t212 + t375) * t497;
t39 = t133 * t449 + t135 * t450 + (t339 * t538 + t341 * t540) * qJD(3) + t577;
t43 = -t447 * t746 - t888;
t550 = t39 * t629 + (-t42 * t538 + t43 * t540) * t215;
t549 = qJD(1) * t378 - qJD(3) * t198 - t253 * t517 + t255 * t518;
t548 = -qJD(3) * t199 - t252 * t517 + t254 * t518 + t754;
t344 = t538 * t903;
t345 = t538 * t904;
t175 = -Icges(7,5) * t345 - Icges(7,6) * t344 - Icges(7,3) * t801;
t346 = t540 * t903;
t347 = t540 * t904;
t176 = -Icges(7,5) * t347 - Icges(7,6) * t346 - Icges(7,3) * t800;
t209 = Icges(7,5) * t397 + Icges(7,6) * t396 - Icges(7,3) * t517;
t544 = (-t176 * t517 - t816) * t449 + (-t209 * t517 - t811) * t497 - (-t175 * t517 - t817) * t450;
t446 = t654 * qJD(3);
t419 = t732 - t760;
t412 = (t538 ^ 2 + t540 ^ 2) * t748;
t401 = t447 * t540;
t316 = t366 * t540;
t315 = t365 * t540;
t314 = t366 * t538;
t313 = t365 * t538;
t290 = t519 + (-t419 - t482) * qJD(1);
t265 = (t438 * t538 + t440 * t540) * qJD(3);
t263 = -qJD(3) * t430 + (t540 * t654 + t526) * qJD(1);
t262 = rSges(4,1) * t581 - rSges(4,2) * t700 + t766;
t249 = rSges(7,1) * t903 - rSges(7,2) * t904;
t225 = -qJD(1) * t291 + t512;
t224 = qJD(1) * (-qJD(1) * t732 + t762) + t770;
t216 = rSges(7,1) * t397 + rSges(7,2) * t396 - rSges(7,3) * t517;
t213 = Icges(7,1) * t397 + Icges(7,4) * t396 - Icges(7,5) * t517;
t211 = Icges(7,4) * t397 + Icges(7,2) * t396 - Icges(7,6) * t517;
t182 = -rSges(7,1) * t347 - rSges(7,2) * t346 - rSges(7,3) * t800;
t181 = -rSges(7,1) * t345 - rSges(7,2) * t344 - rSges(7,3) * t801;
t180 = -Icges(7,1) * t347 - Icges(7,4) * t346 - Icges(7,5) * t800;
t179 = -Icges(7,1) * t345 - Icges(7,4) * t344 - Icges(7,5) * t801;
t178 = -Icges(7,4) * t347 - Icges(7,2) * t346 - Icges(7,6) * t800;
t177 = -Icges(7,4) * t345 - Icges(7,2) * t344 - Icges(7,6) * t801;
t172 = rSges(7,1) * t281 - rSges(7,2) * t282;
t171 = rSges(7,1) * t277 - rSges(7,2) * t278;
t162 = rSges(5,3) * t413 + t653;
t161 = rSges(6,2) * t413 + t650;
t160 = -rSges(5,3) * t705 + t721;
t159 = -rSges(6,2) * t705 + t720;
t101 = -t446 * t745 + t512 + (-t263 + t703 + t783) * qJD(1);
t100 = -t446 * t746 + (t262 - t698) * qJD(1) + t717;
t99 = (t257 * t538 + t261 * t540) * qJD(3) + t656;
t91 = qJD(3) * t679 + (-t257 + t713) * qJD(1) + t764;
t73 = (t256 * t538 + t260 * t540) * qJD(3) + t577;
t62 = -t208 * t802 + t210 * t281 + t212 * t282;
t60 = t62 * t497;
t51 = t784 * t745 + (-t162 + (qJD(3) * t366 - t506) * t538 + t724) * qJD(1) + t710;
t50 = qJD(1) * t160 + (qJD(1) * t679 + t538 * t784) * qJD(3) + t609;
t34 = (t160 * t540 + t162 * t538 + (t257 * t540 + t538 * t791) * qJD(1)) * qJD(3) + t671;
t29 = qJD(1) * t159 + (qJD(1) * t659 + t538 * t719) * qJD(3) + t575;
t28 = t719 * t745 + (-t161 + (qJD(3) * t365 - t506) * t538 + t672) * qJD(1) + t638;
t17 = (t159 * t540 + t161 * t538 + (t256 * t540 + t538 * t723) * qJD(1)) * qJD(3) + t576;
t16 = t118 * t210 + t119 * t212 - t120 * t803 + t121 * t277 + t122 * t278 - t208 * t413;
t15 = t116 * t210 + t117 * t212 - t120 * t802 + t121 * t281 + t122 * t282 + t208 * t580;
t14 = t449 * t49 - t450 * t48 + t497 * t78;
t11 = t449 * t47 - t450 * t46 + t60;
t10 = -t44 * t450 + t449 * t45 + t929;
t7 = t133 * t343 - t135 * t342 + t449 * t72 + t450 * t71 + (t217 * t540 + t218 * t538 + (t339 * t540 + t538 * t722) * qJD(1)) * qJD(3) + t576;
t6 = t118 * t129 + t119 * t132 - t126 * t413 + t277 * t67 + t278 * t69 - t65 * t803;
t5 = -t118 * t128 - t119 * t131 - t124 * t413 + t277 * t68 + t278 * t70 - t66 * t803;
t4 = t116 * t129 + t117 * t132 + t126 * t580 + t281 * t67 + t282 * t69 - t65 * t802;
t3 = -t116 * t128 - t117 * t131 + t124 * t580 + t281 * t68 + t282 * t70 - t66 * t802;
t2 = t16 * t497 + t342 * t44 + t343 * t45 + t449 * t6 - t450 * t5 - t61 * t692;
t1 = t15 * t497 - t3 * t450 + t342 * t46 + t343 * t47 + t4 * t449 - t62 * t692;
t18 = [(t881 + t882) * t746 / 0.2e1 + t612 + (-(-qJD(1) * t392 - t173 + t363 - t463 + t657) * t174 + t101 * (-t392 + t763) + t173 * t761 + t100 * (t393 + t676) + t174 * (t519 + t766) + (t458 * t813 - t812) * qJD(3) + ((-t173 * rSges(4,3) + t174 * (-t510 - t847)) * t538 + (t173 * (-t510 - t654) - t174 * t536) * t540) * qJD(1)) * m(4) + (-t929 + (t47 - t968) * t450 + (t46 + t969) * t449 + t10) * t867 + (((t156 - t325 + (t379 + t808) * t540 + t786) * t540 + (t785 + t930) * t538 + t951) * qJD(3) + t915) * t684 + ((t443 - t943) * t518 + (t532 * t944 + t534 * t942 + t444) * t517 + (t890 * t517 + t518 * t954 - t617) * qJD(3)) * qJD(1) + (t16 + t11) * t865 + (t12 * (t614 - t792) + t13 * (t420 + t594 + t790) + (-t12 * t459 - t13 * t536 + t674 * t747) * t538 + (qJD(1) * t674 + t13 * t859) * t802 + (-qJD(3) * t806 + t682 * t749 + t572 - t647 - t855 - t888) * t42 + (t726 + t788 + ((t518 * t859 - t858) * t540 - t658) * qJD(3) + ((-t459 - t510) * t538 - t513 + t339) * qJD(1) + t939 - t946) * t43) * m(7) + (t880 - t883 + t884) * t685 + (t60 + (t45 + t969) * t450 + (t44 + t968) * t449) * t864 + t15 * t866 + t62 * t868 + t61 * t869 + (t225 * (t538 * t696 + t522 + t760) + t290 * t520 + t224 * (t613 + t757) + t291 * (t759 + t762) + (t290 * (t696 + t846) * t540 + (t290 * (-rSges(3,3) - qJ(2)) + t291 * t696) * t538) * qJD(1) - (-qJD(1) * t419 - t290 + t758) * t291) * m(3) + (t28 * (t614 - t649) + t29 * (t594 + t260) + (t28 * t583 - t29 * t536) * t538 - (t538 * t835 + t540 * t664) * qJD(3) + (t572 - t650 + t834 * t905 + (-t510 + t583) * t749 + t887) * t74 + (-t675 + t720 + ((t517 * t834 + t682) * t538 - t513 + t256) * qJD(1) + t939) * t75) * m(6) + (t51 * (-t652 + t763) + t50 * (t706 + t261) + (-t50 * t536 + t51 * t582) * t538 - t681 * t745 + (-t653 + t666 + t833 * t905 + (-t510 + t582) * t749) * t91 + (-t675 + t721 + t91 + ((t517 * t833 + t682) * t538 - t513 + t257) * qJD(1) + t953) * t92) * m(5) + t852 / 0.2e1 - t853 / 0.2e1 + t838 / 0.2e1 + t839 / 0.2e1 + (((t540 * t677 + t158 - t785 + t957) * t540 + (t538 * t677 + t678 + t916 - t931 + t982) * t538) * qJD(3) + t885 - t912) * t687 + ((t878 + t948) * t540 + (t879 + t949) * t538) * t735 / 0.2e1; 0.2e1 * (t841 / 0.2e1 - t840 / 0.2e1) * m(7) + 0.2e1 * (t28 * t861 + t29 * t860) * m(6) + 0.2e1 * (t50 * t860 + t51 * t861) * m(5) + 0.2e1 * (t100 * t860 + t101 * t861) * m(4) + 0.2e1 * (t224 * t860 + t225 * t861) * m(3); t641 * t660 + t14 * t739 / 0.2e1 + (qJD(1) * t640 + t538 * t9 - t540 * t8) * t862 + ((-t126 * t517 + t129 * t396 + t132 * t397 + t176 * t518 + t178 * t903 + t180 * t904) * t449 - (-t124 * t517 - t128 * t396 - t131 * t397 + t175 * t518 + t177 * t903 + t179 * t904) * t450 + (-t208 * t517 + t209 * t518 + t210 * t396 + t211 * t903 + t212 * t397 + t213 * t904) * t497 + (-t78 * t517 - t518 * t640) * qJD(6)) * t863 + ((-t129 * t344 - t132 * t345 + t178 * t277 + t180 * t278) * t449 - (t128 * t344 + t131 * t345 + t177 * t277 + t179 * t278) * t450 + (-t210 * t344 + t211 * t277 - t212 * t345 + t213 * t278) * t497 + (-t45 * t800 - t517 * t61) * qJD(6) + (-t44 * t738 + t544) * t538) * t864 + (qJD(1) * t644 - t5 * t540 + t538 * t6) * t865 + (qJD(1) * t642 - t3 * t540 + t4 * t538) * t866 + ((-t129 * t346 - t132 * t347 + t178 * t281 + t180 * t282) * t449 - (t128 * t346 + t131 * t347 + t177 * t281 + t179 * t282) * t450 + (-t210 * t346 + t211 * t281 - t212 * t347 + t213 * t282) * t497 + (-t46 * t801 - t517 * t62) * qJD(6) + (-t47 * t738 + t544) * t540) * t867 + t643 * t868 + t645 * t869 + (t10 * t538 + t11 * t540) * t738 / 0.2e1 + (t42 * t778 + t7 * t668 + t39 * t610 + (t12 * t670 + t42 * t673 + t7 * (t135 + t341) + t39 * (t217 + t71) + (t39 * t792 + t43 * t670) * qJD(1)) * t540 + (t13 * t670 + t43 * t673 + t7 * t792 + t39 * (t218 + t72) + (t42 * (t215 + t447) + t39 * (-t135 + t722)) * qJD(1)) * t538 - t42 * (-t181 * t497 - t216 * t450 + t611) - t43 * (t182 * t497 - t216 * t449 + t578) - t39 * (t181 * t449 + t182 * t450 + t608) - (t42 * (t806 + t780) - t43 * t401) * qJD(1) - ((t133 * t42 - t135 * t43) * t517 + t550 * t518) * qJD(6) - ((-t39 * t401 + t42 * t708) * t540 + (-t39 * t806 + t43 * t708) * t538) * qJD(3)) * m(7) + (-t74 * t611 - t75 * t578 - t73 * t608 - (t74 * (t313 + t780) - t75 * t315) * qJD(1) - ((-t73 * t315 + t714 * t74) * t540 + (-t73 * t313 + t714 * t75) * t538) * qJD(3) + t74 * t778 + t17 * t668 + t73 * t610 + (t28 * t715 + t74 * t719 + t17 * t260 + t73 * t159 + (t73 * t256 + t664) * qJD(1)) * t540 + (t29 * t715 + t75 * t719 + t17 * t256 + t73 * t161 + (t723 * t73 + t835) * qJD(1)) * t538) * m(6) + (t91 * t437 + t34 * t774 + t99 * t725 + (t51 * t782 + t91 * t784 + t34 * t261 + t99 * t160 + (t99 * t257 + t681) * qJD(1)) * t540 + (t50 * t782 + t92 * t784 + t34 * t257 + t99 * t162 + (t91 * t366 + t791 * t99) * qJD(1)) * t538 - t91 * (t481 + (t314 + t429) * qJD(1)) - t92 * (-qJD(1) * t316 + t773) - t99 * t712 - ((-t99 * t316 + t781 * t91) * t540 + (-t99 * t314 + t781 * t92) * t538) * qJD(3)) * m(5) + (0.2e1 * t214 * (t262 * t540 + t263 * t538 + (t392 * t540 - t393 * t538) * qJD(1)) + t628 * t446 + (-t100 * t538 - t101 * t540 + (-t174 * t540 + t813) * qJD(1)) * t458 - (t173 * t430 - t812) * qJD(1) - (t214 * (-t430 * t538 - t433 * t540) + t628 * t654) * qJD(3)) * m(4) - ((((-t777 - t937) * t540 + (t776 - t938) * t538) * qJD(3) + t950) * t518 + ((t573 + (t532 * t897 - t534 * t895 - t920) * t540 + (-t532 * t896 + t534 * t894 + t919) * t538) * qJD(3) + (t532 * t891 + t534 * t889 + t952) * qJD(1)) * t517) * qJD(1) / 0.2e1 + (t883 * t540 + t882 * t538 + (t879 * t538 + t878 * t540) * qJD(1)) * qJD(1) / 0.2e1 + (t752 * t538 + (-t440 * t896 + t441 * t894 - t538 * t804) * t746 + (t440 * t891 + t441 * t889) * qJD(1) + ((t440 * t897 - t441 * t895 + t538 * t805 + t876) * qJD(3) + t877) * t540) * t687 + (-t752 * t540 + (t438 * t897 - t439 * t895 - t540 * t805) * t745 + (t438 * t891 + t439 * t889) * qJD(1) + ((-t438 * t896 + t439 * t894 + t540 * t804 + t876) * qJD(3) + t877) * t538) * t684 + (t881 * qJD(1) + t1 + ((t922 * t317 - t318 * t970 + t927 * t440 + t923 * t441 - t549 * t540 + t920 * t580 + t925 * t802) * t540 + (t875 * t538 + t926 * t802 - t919 * t580 + (t548 - t874) * t540 + t924 * t441 + t928 * t440 + t918 * t318 - t921 * t317) * t538 + (t916 * t538 + t540 * t959) * qJD(1)) * t907) * t861 + (t880 * qJD(1) + t2 + ((-t922 * t319 - t320 * t970 - t920 * t413 + t927 * t438 + t923 * t439 + t874 * t540 + t925 * t803) * t540 + (t548 * t538 + t926 * t803 + (-t549 - t875) * t540 + t924 * t439 + t928 * t438 + t919 * t413 + t918 * t320 + t921 * t319) * t538 + (t538 * t961 + t540 * t960) * qJD(1)) * t907) * t860 + (t10 + t885 + t913) * t690 + (t11 + t884 + t914) * t749 / 0.2e1; -m(5) * (t412 * t99 + t413 * t92 - t580 * t91) - m(6) * (t412 * t73 + t413 * t75 - t580 * t74) - m(7) * (t39 * t412 + t413 * t43 - t42 * t580) + 0.2e1 * ((t745 * t91 + t746 * t92 - t34) * t872 + (t74 * t745 + t746 * t75 - t17) * t871 + (t42 * t745 + t43 * t746 - t7) * t870) * t518 + 0.2e1 * ((qJD(3) * t99 + t50 * t538 + t51 * t540 + t749 * t92 - t750 * t91) * t872 + (qJD(3) * t73 + t28 * t540 + t29 * t538 - t74 * t750 + t749 * t75) * t871 + (qJD(3) * t39 + t12 * t540 + t13 * t538 - t42 * t750 + t43 * t749) * t870) * t517; (t12 * t440 + t13 * t438 + (t39 * t747 + t517 * t7) * t532 - t265 * t39 + t892 * t43) * m(7) + (t28 * t440 + t29 * t438 + (t17 * t517 + t73 * t747) * t532 - t265 * t73 + t892 * t75) * m(6); t1 * t697 + (-t517 * t642 + t518 * t62) * t868 + ((-qJD(3) * t642 + t15) * t518 + (qJD(1) * t643 - qJD(3) * t62 - t3 * t538 - t4 * t540) * t517) * t866 - t2 * t803 / 0.2e1 + (-t517 * t644 + t518 * t61) * t869 + ((-qJD(3) * t644 + t16) * t518 + (qJD(1) * t645 - qJD(3) * t61 - t5 * t538 - t540 * t6) * t517) * t865 + t14 * t688 + t518 * (t612 + t838 + t839 + t852 - t853) / 0.2e1 + (-t517 * t640 + t518 * t78) * t660 + ((-qJD(3) * t640 + t19) * t518 + (qJD(1) * t641 - qJD(3) * t78 - t538 * t8 - t540 * t9) * t517) * t862 + (t281 * t558 + t282 * t559 + t540 * t571) * t867 + (t277 * t558 + t278 * t559 + t538 * t571) * t864 + (-t518 * t579 + t558 * t903 + t559 * t904) * t863 + (t517 * t690 + t518 * t685) * t11 + (qJD(1) * t697 + t518 * t687) * t10 + ((qJD(3) * t550 - t12 * t133 + t13 * t135 - t42 * t72 + t43 * t71) * t518 + (t42 * (qJD(3) * t133 - t123 * t538) + t43 * (-qJD(3) * t135 + t123 * t540) + t7 * t629 + t39 * (t133 * t750 + t135 * t749 + t538 * t71 - t540 * t72) + (-t841 + t840 + (-t42 * t540 - t43 * t538) * qJD(1)) * t215) * t517 - t42 * (-t171 * t497 - t249 * t450) - t43 * (t172 * t497 - t249 * t449) - t39 * (t171 * t449 + t172 * t450)) * m(7);];
tauc  = t18(:);

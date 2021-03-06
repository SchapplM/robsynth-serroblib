% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR10_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR10_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:03
% EndTime: 2019-12-31 19:44:23
% DurationCPUTime: 69.54s
% Computational Cost: add. (17784->1134), mult. (48297->1466), div. (0->0), fcn. (49974->8), ass. (0->505)
t903 = Icges(4,1) + Icges(5,1);
t902 = Icges(4,4) - Icges(5,5);
t901 = Icges(5,4) + Icges(4,5);
t914 = Icges(4,2) + Icges(5,3);
t907 = Icges(5,2) + Icges(4,3);
t898 = Icges(4,6) - Icges(5,6);
t494 = cos(pkin(8));
t497 = sin(qJ(1));
t499 = cos(qJ(2));
t493 = sin(pkin(8));
t500 = cos(qJ(1));
t736 = t500 * t493;
t397 = -t497 * t494 + t499 * t736;
t737 = t499 * t500;
t398 = t493 * t497 + t494 * t737;
t738 = t497 * t499;
t395 = t493 * t738 + t494 * t500;
t396 = t494 * t738 - t736;
t496 = sin(qJ(2));
t740 = t496 * t497;
t851 = t914 * t395 - t902 * t396 - t898 * t740;
t894 = -t902 * t395 + t903 * t396 + t901 * t740;
t917 = -t397 * t851 - t398 * t894;
t849 = -t898 * t395 + t901 * t396 + t907 * t740;
t913 = t914 * t493 - t902 * t494;
t912 = -t898 * t493 + t901 * t494;
t911 = -t902 * t493 + t903 * t494;
t892 = t395 * t851 + t396 * t894;
t879 = t493 * t851 + t494 * t894;
t739 = t496 * t500;
t906 = -t849 * t739 + t917;
t850 = t914 * t397 - t902 * t398 - t898 * t739;
t848 = -t898 * t397 + t901 * t398 + t907 * t739;
t846 = -t902 * t397 + t903 * t398 + t901 * t739;
t872 = t913 * t496 + t898 * t499;
t845 = t912 * t496 - t907 * t499;
t871 = t911 * t496 - t901 * t499;
t757 = Icges(3,3) * t500;
t358 = Icges(3,5) * t738 - Icges(3,6) * t740 - t757;
t469 = Icges(3,4) * t740;
t764 = Icges(3,5) * t500;
t362 = Icges(3,1) * t738 - t469 - t764;
t760 = Icges(3,6) * t500;
t360 = Icges(3,4) * t738 - Icges(3,2) * t740 - t760;
t745 = t360 * t496;
t575 = -t362 * t499 + t745;
t897 = -t358 * t500 - t497 * t575 + t740 * t849 + t892;
t487 = Icges(3,4) * t499;
t590 = -Icges(3,2) * t496 + t487;
t361 = Icges(3,6) * t497 + t500 * t590;
t769 = Icges(3,4) * t496;
t439 = Icges(3,1) * t499 - t769;
t363 = Icges(3,5) * t497 + t439 * t500;
t310 = t363 * t738;
t435 = Icges(3,5) * t499 - Icges(3,6) * t496;
t359 = Icges(3,3) * t497 + t435 * t500;
t632 = t359 * t500 - t310;
t162 = -t361 * t740 - t632;
t861 = t395 * t850 + t396 * t846 + t740 * t848;
t896 = t162 + t861;
t727 = t497 * t359 + t363 * t737;
t164 = -t361 * t739 + t727;
t860 = t397 * t850 + t398 * t846 + t739 * t848;
t895 = t164 + t860;
t819 = -t898 * t496 + t913 * t499;
t893 = t907 * t496 + t912 * t499;
t818 = t901 * t496 + t911 * t499;
t495 = sin(qJ(5));
t498 = cos(qJ(5));
t267 = t395 * t498 - t396 * t495;
t268 = t395 * t495 + t396 * t498;
t116 = Icges(6,5) * t268 + Icges(6,6) * t267 - Icges(6,3) * t740;
t271 = t397 * t498 - t398 * t495;
t272 = t397 * t495 + t398 * t498;
t118 = Icges(6,5) * t272 + Icges(6,6) * t271 - Icges(6,3) * t739;
t767 = Icges(6,4) * t268;
t120 = -Icges(6,2) * t267 + Icges(6,6) * t740 - t767;
t248 = Icges(6,4) * t267;
t123 = -Icges(6,1) * t268 + Icges(6,5) * t740 - t248;
t787 = t120 * t271 + t123 * t272;
t766 = Icges(6,4) * t272;
t121 = Icges(6,2) * t271 - Icges(6,6) * t739 + t766;
t249 = Icges(6,4) * t271;
t124 = Icges(6,1) * t272 - Icges(6,5) * t739 + t249;
t788 = -t267 * t121 - t268 * t124;
t891 = t787 + t788 + (t116 * t500 + t118 * t497) * t496;
t786 = t271 * t121 + t272 * t124;
t831 = -t120 * t267 - t123 * t268;
t890 = -t831 + t496 * (t116 * t497 - t118 * t500) + t786;
t889 = t906 * t500;
t436 = Icges(3,2) * t499 + t769;
t438 = Icges(3,1) * t496 + t487;
t573 = t436 * t496 - t438 * t499;
t434 = Icges(3,5) * t496 + Icges(3,6) * t499;
t741 = t434 * t500;
t888 = -t395 * t872 - t396 * t871 + t497 * t573 - t740 * t845 + t741;
t742 = t434 * t497;
t887 = t397 * t872 + t398 * t871 - t500 * t573 + t739 * t845 + t742;
t755 = t116 * t499;
t572 = t493 * t498 - t494 * t495;
t827 = t572 * t496;
t571 = t493 * t495 + t494 * t498;
t828 = t571 * t496;
t46 = -t120 * t827 - t123 * t828 + t755;
t693 = qJD(2) * t500;
t653 = t496 * t693;
t287 = qJD(1) * t395 + t493 * t653;
t288 = -qJD(1) * t396 - t494 * t653;
t652 = t499 * t693;
t699 = qJD(1) * t497;
t658 = t496 * t699;
t538 = -t652 + t658;
t857 = -t287 * t914 - t288 * t902 + t538 * t898;
t695 = qJD(2) * t497;
t654 = t496 * t695;
t433 = t493 * t654;
t697 = qJD(1) * t500;
t657 = t499 * t697;
t289 = t493 * t657 - t494 * t699 - t433;
t290 = qJD(1) * t398 - t494 * t654;
t694 = qJD(2) * t499;
t403 = t496 * t697 + t497 * t694;
t856 = -t289 * t914 + t290 * t902 + t403 * t898;
t855 = t287 * t898 + t288 * t901 - t538 * t907;
t854 = t289 * t898 - t290 * t901 - t403 * t907;
t853 = t287 * t902 + t288 * t903 - t538 * t901;
t852 = t289 * t902 - t290 * t903 - t403 * t901;
t884 = t819 * qJD(2);
t883 = t893 * qJD(2);
t882 = t818 * qJD(2);
t881 = t493 * t872 + t494 * t871;
t880 = t493 * t850 + t494 * t846;
t878 = -t436 + t439 + t845;
t877 = (t438 + t590 + t881 - t893) * qJD(1);
t125 = rSges(6,1) * t268 + rSges(6,2) * t267 - rSges(6,3) * t740;
t207 = rSges(6,1) * t828 + rSges(6,2) * t827 + rSges(6,3) * t499;
t688 = qJD(5) * t496;
t422 = t497 * t688 + t693;
t687 = qJD(5) * t499;
t475 = qJD(1) + t687;
t875 = -t125 * t475 - t207 * t422;
t544 = qJD(2) * t434;
t873 = t500 * (-t497 * t544 + (t359 + t575) * qJD(1));
t870 = -t500 * t845 - t880;
t869 = t497 * t845 + t879;
t441 = rSges(3,1) * t496 + rSges(3,2) * t499;
t705 = t500 * pkin(1) + t497 * pkin(6);
t843 = qJD(1) * t705;
t868 = t441 * t695 - t843;
t698 = qJD(1) * t499;
t867 = t887 * qJD(1);
t728 = -t497 * t358 - t362 * t737;
t163 = -t360 * t739 - t728;
t866 = (-t163 * t500 + t497 * t895 + t889) * qJD(2);
t865 = (t497 * t896 - t500 * t897) * qJD(2);
t863 = t888 * qJD(1);
t201 = Icges(6,5) * t828 + Icges(6,6) * t827 + Icges(6,3) * t499;
t765 = Icges(6,4) * t828;
t203 = Icges(6,2) * t827 + Icges(6,6) * t499 + t765;
t339 = Icges(6,4) * t827;
t205 = Icges(6,1) * t828 + Icges(6,5) * t499 + t339;
t68 = -t201 * t740 + t203 * t267 + t205 * t268;
t858 = t475 * t68;
t794 = pkin(4) * t494;
t428 = pkin(7) * t499 + t496 * t794;
t440 = pkin(2) * t496 - qJ(3) * t499;
t601 = pkin(3) * t494 + qJ(4) * t493;
t826 = t601 * t496;
t714 = -t826 - t440;
t664 = -t428 + t714;
t612 = t664 * t500;
t277 = t398 * pkin(3) + qJ(4) * t397;
t417 = pkin(2) * t737 + qJ(3) * t739;
t483 = qJD(3) * t496;
t462 = t497 * t483;
t622 = qJD(1) * t417 - t440 * t695 + t462 + t843;
t690 = qJD(4) * t395;
t844 = -qJD(1) * t277 + t695 * t826 - t622 - t690;
t842 = -t863 + t865;
t841 = t866 + t867;
t545 = qJD(2) * t436;
t546 = qJD(2) * t438;
t840 = (-qJD(1) * t361 + t497 * t545 - t854) * t499 + (-qJD(1) * t363 + t493 * t856 + t494 * t852 + t497 * t546) * t496 + (-t496 * t849 - t499 * t879 + t575) * qJD(2);
t744 = t361 * t496;
t574 = -t363 * t499 + t744;
t839 = (-t500 * t545 + (-t497 * t590 + t760) * qJD(1) - t855) * t499 + (-t500 * t546 + (-t439 * t497 + t764) * qJD(1) + t853 * t494 + t857 * t493) * t496 + (t496 * t848 + t499 * t880 - t574) * qJD(2);
t424 = t590 * qJD(2);
t425 = t439 * qJD(2);
t506 = qJD(1) * t434 - t424 * t496 + t425 * t499 + (-t436 * t499 - t438 * t496) * qJD(2);
t811 = t573 * qJD(1) + t435 * qJD(2);
t838 = -t287 * t872 + t288 * t871 + t397 * t884 + t398 * t882 + t497 * t811 + t506 * t500 - t538 * t845 + t739 * t883;
t837 = t289 * t872 + t290 * t871 + t395 * t884 + t396 * t882 + t403 * t845 + t506 * t497 - t500 * t811 + t740 * t883;
t836 = t163 - t906;
t835 = (t360 - t849) * t499 + (t362 + t879) * t496;
t834 = (t361 - t848) * t499 + (t363 + t880) * t496;
t830 = 0.2e1 * qJD(2);
t829 = pkin(6) * qJD(1);
t824 = -qJD(1) * t397 + t289 + t433;
t823 = t872 * t497;
t822 = t872 * t500;
t821 = t871 * t497;
t820 = t871 * t500;
t275 = pkin(3) * t396 + t395 * qJ(4);
t732 = t272 * rSges(6,1) + t271 * rSges(6,2);
t127 = -rSges(6,3) * t739 + t732;
t387 = t398 * pkin(4);
t317 = -pkin(7) * t739 + t387;
t421 = -t500 * t688 + t695;
t817 = -qJD(1) * t317 - t127 * t475 + t207 * t421 + t844;
t236 = t398 * rSges(5,1) + rSges(5,2) * t739 + t397 * rSges(5,3);
t816 = -qJD(1) * t236 + t844;
t815 = ((t497 * t870 + t500 * t869) * qJD(2) - t877) * t496 + t878 * t698;
t532 = t360 * t500 - t361 * t497;
t720 = -t436 * t500 + t363;
t721 = -Icges(3,2) * t738 + t362 - t469;
t814 = (-t497 * t720 + t500 * t721) * t496 + (t497 * t848 - t500 * t849 + t532) * t499;
t603 = rSges(5,1) * t494 + rSges(5,3) * t493;
t354 = -rSges(5,2) * t499 + t496 * t603;
t606 = rSges(4,1) * t494 - rSges(4,2) * t493;
t355 = -rSges(4,3) * t499 + t496 * t606;
t813 = -t500 * t544 + (-t435 * t497 + t574 + t757) * qJD(1);
t810 = m(4) / 0.2e1;
t809 = m(5) / 0.2e1;
t808 = m(6) / 0.2e1;
t683 = qJD(2) * qJD(5);
t644 = t499 * t683;
t324 = qJD(1) * t421 - t497 * t644;
t807 = t324 / 0.2e1;
t325 = qJD(1) * t422 - t500 * t644;
t806 = t325 / 0.2e1;
t805 = -t421 / 0.2e1;
t804 = t421 / 0.2e1;
t803 = -t422 / 0.2e1;
t802 = t422 / 0.2e1;
t801 = -t475 / 0.2e1;
t800 = t475 / 0.2e1;
t797 = -rSges(6,3) - pkin(7);
t795 = pkin(4) * t290;
t793 = t499 * pkin(2);
t375 = t572 * t499;
t239 = qJD(2) * t375 - qJD(5) * t828;
t376 = t571 * t499;
t240 = qJD(2) * t376 + qJD(5) * t827;
t104 = -qJD(5) * t268 + t289 * t498 - t290 * t495;
t105 = qJD(5) * t267 + t289 * t495 + t290 * t498;
t59 = Icges(6,5) * t105 + Icges(6,6) * t104 - Icges(6,3) * t403;
t61 = Icges(6,4) * t105 + Icges(6,2) * t104 - Icges(6,6) * t403;
t63 = Icges(6,1) * t105 + Icges(6,4) * t104 - Icges(6,5) * t403;
t696 = qJD(2) * t496;
t8 = -t116 * t696 - t120 * t239 - t123 * t240 + t499 * t59 + t61 * t827 + t63 * t828;
t792 = t8 * t422;
t102 = -qJD(5) * t272 - t287 * t498 - t288 * t495;
t103 = qJD(5) * t271 - t287 * t495 + t288 * t498;
t58 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t538;
t60 = Icges(6,4) * t103 + Icges(6,2) * t102 + Icges(6,6) * t538;
t62 = Icges(6,1) * t103 + Icges(6,4) * t102 + Icges(6,5) * t538;
t9 = -t118 * t696 + t121 * t239 + t124 * t240 + t499 * t58 + t60 * t827 + t62 * t828;
t791 = t9 * t421;
t785 = rSges(3,1) * t499;
t784 = rSges(5,2) * t496;
t782 = rSges(4,3) * t496;
t115 = rSges(6,1) * t240 + rSges(6,2) * t239 - rSges(6,3) * t696;
t193 = -pkin(7) * t403 + t795;
t685 = qJD(1) * qJD(2);
t648 = t497 * t685;
t684 = qJD(2) * qJD(3);
t646 = t499 * t684;
t715 = t440 * t648 + t500 * t646;
t623 = -qJD(4) * t287 + t648 * t826 + t715;
t602 = rSges(6,1) * t105 + rSges(6,2) * t104;
t65 = -rSges(6,3) * t403 + t602;
t429 = -pkin(7) * t496 + t499 * t794;
t443 = qJ(3) * t496 + t793;
t692 = qJD(3) * t499;
t391 = qJD(2) * t443 - t692;
t401 = t601 * t499;
t689 = qJD(4) * t493;
t459 = t496 * t689;
t729 = -qJD(2) * t401 - t391 - t459;
t668 = -t429 * qJD(2) + t729;
t144 = pkin(3) * t290 + t289 * qJ(4) + t690;
t707 = -pkin(2) * t654 + t462;
t243 = pkin(2) * t657 + qJ(3) * t403 + t707;
t733 = -t243 - t843;
t676 = -t144 + t733;
t12 = -t115 * t422 + t207 * t324 - t475 * t65 + (t125 * t688 + t500 * t668) * qJD(2) + (-t193 + (qJD(2) * t428 - t483) * t497 + t676) * qJD(1) + t623;
t780 = t12 * t497;
t315 = pkin(4) * t396 - pkin(7) * t740;
t414 = t443 * t497;
t491 = t500 * pkin(6);
t446 = pkin(1) * t497 - t491;
t712 = -t414 - t446;
t673 = -t275 + t712;
t381 = qJD(4) * t397;
t691 = qJD(3) * t500;
t464 = t496 * t691;
t717 = t381 + t464;
t44 = qJD(2) * t612 + (-t315 + t673) * qJD(1) + t717 + t875;
t779 = t44 * t497;
t778 = t44 * t500;
t777 = t46 * t324;
t753 = t118 * t499;
t47 = t121 * t827 + t124 * t828 + t753;
t776 = t47 * t325;
t488 = t497 * rSges(3,3);
t604 = rSges(5,1) * t396 + rSges(5,3) * t395;
t232 = rSges(5,2) * t740 + t604;
t667 = -t354 + t714;
t613 = t500 * t667;
t72 = qJD(2) * t613 + (-t232 + t673) * qJD(1) + t717;
t773 = t72 * t354;
t772 = -rSges(5,2) - qJ(3);
t771 = -rSges(4,3) - qJ(3);
t706 = rSges(3,2) * t740 + t500 * rSges(3,3);
t371 = rSges(3,1) * t738 - t706;
t655 = t441 * t693;
t195 = -t655 + (-t371 - t446) * qJD(1);
t751 = t195 * t497;
t750 = t195 * t500;
t372 = rSges(3,1) * t737 - rSges(3,2) * t739 + t488;
t196 = qJD(1) * t372 - t868;
t416 = t441 * t500;
t749 = t196 * t416;
t748 = t201 * t499;
t743 = t428 * t497;
t735 = t125 + t315;
t237 = t398 * rSges(4,1) - t397 * rSges(4,2) + rSges(4,3) * t739;
t734 = -t237 - t417;
t731 = -t277 - t417;
t730 = t288 * pkin(4) + pkin(7) * t658;
t357 = t499 * t606 + t782;
t726 = -t357 * qJD(2) - t391;
t337 = t497 * t826;
t412 = t440 * t497;
t725 = t337 + t412;
t419 = t440 * t699;
t724 = t699 * t826 + t419;
t723 = -t355 - t440;
t722 = -t357 - t443;
t718 = t497 * t414 + t500 * t417;
t415 = t440 * t500;
t716 = -qJD(1) * t415 + t497 * t692;
t713 = -t401 - t443;
t709 = qJ(3) * t652 + t464;
t708 = rSges(3,2) * t658 + rSges(3,3) * t697;
t700 = qJD(1) * t435;
t682 = -qJ(3) - t797;
t677 = t103 * rSges(6,1) + t102 * rSges(6,2) + rSges(6,3) * t658;
t539 = -t497 * t698 - t653;
t242 = pkin(2) * t539 - qJ(3) * t658 + t709;
t675 = t500 * t242 + t497 * t243 + t414 * t697;
t674 = -t236 + t731;
t672 = -t317 + t731;
t671 = t288 * rSges(4,1) + t287 * rSges(4,2) + rSges(4,3) * t652;
t670 = t288 * rSges(5,1) + rSges(5,2) * t652 - t287 * rSges(5,3);
t356 = t499 * t603 + t784;
t669 = -t356 * qJD(2) + t729;
t666 = -t356 + t713;
t665 = -t412 * t695 - t415 * t693 + t483;
t663 = -t429 + t713;
t482 = pkin(6) * t697;
t662 = t482 + t709;
t431 = qJD(1) * t446;
t661 = -qJD(1) * t414 - t431 + t464;
t660 = -pkin(1) - t793;
t659 = t797 * t500;
t460 = t499 * t689;
t650 = -t739 / 0.2e1;
t649 = -pkin(1) - t785;
t647 = t500 * t685;
t645 = t496 * t683;
t643 = t699 / 0.2e1;
t641 = -t696 / 0.2e1;
t640 = -t695 / 0.2e1;
t638 = -t693 / 0.2e1;
t637 = t693 / 0.2e1;
t99 = qJD(1) * t237 - t355 * t695 + t622;
t635 = t99 * t723;
t633 = t723 * t500;
t631 = -t358 + t744;
t630 = pkin(2) * t653;
t629 = -t115 + t668;
t628 = t242 * t693 + t243 * t695 + t414 * t647 + t496 * t684;
t627 = -t207 + t664;
t420 = qJD(1) * (-pkin(1) * t699 + t482);
t626 = t497 * t646 + t420 + (t242 + t464) * qJD(1);
t624 = t497 * t275 + t500 * t277 + t718;
t621 = t705 + t417;
t619 = -t275 + t491;
t73 = -t354 * t695 - t816;
t618 = t73 * t667;
t614 = qJD(5) * t641;
t611 = t414 * t695 + t417 * t693 - t692;
t609 = -rSges(3,2) * t496 + t785;
t608 = rSges(4,1) * t290 - rSges(4,2) * t289;
t607 = -rSges(4,1) * t396 + rSges(4,2) * t395;
t605 = rSges(5,1) * t290 + rSges(5,3) * t289;
t40 = -t116 * t740 + t831;
t41 = -t118 * t740 - t788;
t600 = -t40 * t500 + t41 * t497;
t599 = t40 * t497 + t41 * t500;
t42 = -t116 * t739 - t787;
t43 = -t118 * t739 + t786;
t598 = -t42 * t500 + t43 * t497;
t597 = t42 * t497 + t43 * t500;
t596 = -t46 * t500 + t47 * t497;
t595 = t46 * t497 + t47 * t500;
t594 = -qJD(1) * t275 + t381 + t661;
t585 = -t125 * t500 + t127 * t497;
t584 = -t196 * t497 - t750;
t112 = Icges(6,5) * t240 + Icges(6,6) * t239 - Icges(6,3) * t696;
t113 = Icges(6,4) * t240 + Icges(6,2) * t239 - Icges(6,6) * t696;
t114 = Icges(6,1) * t240 + Icges(6,4) * t239 - Icges(6,5) * t696;
t19 = t112 * t499 + t113 * t827 + t114 * t828 - t201 * t696 + t203 * t239 + t205 * t240;
t74 = t203 * t827 + t205 * t828 + t748;
t570 = t19 * t475 - t645 * t74;
t569 = -pkin(1) - t443;
t465 = t499 * t691;
t568 = -t459 * t500 + t465;
t143 = t288 * pkin(3) - qJ(4) * t287 + t381;
t567 = t500 * t143 + t497 * t144 + t275 * t697 + t675;
t566 = qJD(1) * t143 + qJD(4) * t289 + t626;
t338 = t500 * t826;
t564 = -t337 * t695 - t338 * t693 + t460 + t665;
t557 = -t144 - t707;
t413 = t441 * t497;
t191 = (t371 * t497 + t372 * t500) * qJD(2);
t541 = t496 * t772 + t660;
t540 = t496 * t771 + t660;
t537 = t277 + t621;
t536 = (Icges(6,5) * t267 - Icges(6,6) * t268) * t422 - (Icges(6,5) * t271 - Icges(6,6) * t272) * t421 - (Icges(6,5) * t827 - Icges(6,6) * t828) * t475;
t535 = -qJD(1) * t338 - t459 * t497 + t716;
t534 = t275 * t695 + t277 * t693 + t459 + t611;
t533 = qJD(2) * t460 + t143 * t693 + t144 * t695 + t275 * t647 + t628;
t233 = rSges(4,3) * t740 - t607;
t530 = t496 * t536;
t529 = t143 + t662;
t519 = (Icges(6,1) * t271 - t121 - t766) * t421 - (Icges(6,1) * t267 + t120 - t767) * t422 + (Icges(6,1) * t827 - t203 - t765) * t475;
t518 = (-Icges(6,2) * t272 + t124 + t249) * t421 - (-Icges(6,2) * t268 - t123 + t248) * t422 + (-Icges(6,2) * t828 + t205 + t339) * t475;
t39 = t125 * t421 + t127 * t422 + (t315 * t497 + t317 * t500) * qJD(2) + t534;
t45 = -t428 * t695 - t817;
t509 = t39 * t585 + (t45 * t500 - t779) * t207;
t328 = t497 * t827;
t329 = t497 * t828;
t177 = -Icges(6,5) * t329 - Icges(6,6) * t328 - Icges(6,3) * t738;
t330 = t500 * t827;
t331 = t500 * t828;
t178 = -Icges(6,5) * t331 - Icges(6,6) * t330 - Icges(6,3) * t737;
t202 = Icges(6,5) * t376 + Icges(6,6) * t375 - Icges(6,3) * t496;
t504 = (-t178 * t496 - t753) * t421 + (-t202 * t496 - t748) * t475 - (-t177 * t496 - t755) * t422;
t426 = t609 * qJD(2);
t402 = (t497 ^ 2 + t500 ^ 2) * t696;
t390 = t428 * t500;
t308 = t355 * t500;
t307 = t354 * t500;
t306 = t355 * t497;
t305 = t354 * t497;
t264 = -qJD(2) * t413 + (t500 * t609 + t488) * qJD(1);
t263 = rSges(3,1) * t539 - rSges(3,2) * t652 + t708;
t247 = rSges(6,1) * t827 - rSges(6,2) * t828;
t241 = (t497 * t395 + t397 * t500) * qJD(2);
t208 = rSges(6,1) * t376 + rSges(6,2) * t375 - rSges(6,3) * t496;
t206 = Icges(6,1) * t376 + Icges(6,4) * t375 - Icges(6,5) * t496;
t204 = Icges(6,4) * t376 + Icges(6,2) * t375 - Icges(6,6) * t496;
t192 = -pkin(7) * t652 + t730;
t184 = -rSges(6,1) * t331 - rSges(6,2) * t330 - rSges(6,3) * t737;
t183 = -rSges(6,1) * t329 - rSges(6,2) * t328 - rSges(6,3) * t738;
t182 = -Icges(6,1) * t331 - Icges(6,4) * t330 - Icges(6,5) * t737;
t181 = -Icges(6,1) * t329 - Icges(6,4) * t328 - Icges(6,5) * t738;
t180 = -Icges(6,4) * t331 - Icges(6,2) * t330 - Icges(6,6) * t737;
t179 = -Icges(6,4) * t329 - Icges(6,2) * t328 - Icges(6,6) * t738;
t160 = rSges(4,3) * t403 + t608;
t159 = rSges(5,2) * t403 + t605;
t158 = -rSges(4,3) * t658 + t671;
t157 = -rSges(5,2) * t658 + t670;
t141 = rSges(6,1) * t271 - rSges(6,2) * t272;
t140 = rSges(6,1) * t267 - rSges(6,2) * t268;
t129 = -t426 * t693 + (-t264 + t868) * qJD(1);
t128 = -t426 * t695 + t420 + (t263 - t655) * qJD(1);
t98 = t464 + qJD(2) * t633 + (-t233 + t712) * qJD(1);
t87 = (t233 * t497 + t237 * t500) * qJD(2) + t611;
t69 = -t201 * t739 + t203 * t271 + t205 * t272;
t67 = (t232 * t497 + t236 * t500) * qJD(2) + t534;
t66 = t69 * t475;
t64 = -rSges(6,3) * t652 + t677;
t57 = t726 * t693 + (-t160 + (qJD(2) * t355 - t483) * t497 + t733) * qJD(1) + t715;
t56 = qJD(1) * t158 + (qJD(1) * t633 + t497 * t726) * qJD(2) + t626;
t34 = (t158 * t500 + t160 * t497 + (t233 * t500 + t497 * t734) * qJD(1)) * qJD(2) + t628;
t33 = t669 * t693 + (-t159 + (qJD(2) * t354 - t483) * t497 + t676) * qJD(1) + t623;
t32 = qJD(1) * t157 + (qJD(1) * t613 + t497 * t669) * qJD(2) + t566;
t17 = (t157 * t500 + t159 * t497 + (t232 * t500 + t497 * t674) * qJD(1)) * qJD(2) + t533;
t16 = t104 * t203 + t105 * t205 - t112 * t740 + t113 * t267 + t114 * t268 - t201 * t403;
t15 = t102 * t203 + t103 * t205 - t112 * t739 + t113 * t271 + t114 * t272 + t201 * t538;
t14 = t421 * t47 - t422 * t46 + t475 * t74;
t13 = qJD(1) * t192 - t115 * t421 - t207 * t325 + t475 * t64 + (qJD(1) * t612 - t127 * t688 + t497 * t668) * qJD(2) + t566;
t11 = -t42 * t422 + t421 * t43 + t66;
t10 = -t40 * t422 + t41 * t421 + t858;
t7 = t125 * t325 - t127 * t324 + t421 * t65 + t422 * t64 + (t192 * t500 + t193 * t497 + (t315 * t500 + t497 * t672) * qJD(1)) * qJD(2) + t533;
t6 = t104 * t121 + t105 * t124 - t118 * t403 + t267 * t60 + t268 * t62 - t58 * t740;
t5 = -t104 * t120 - t105 * t123 - t116 * t403 + t267 * t61 + t268 * t63 - t59 * t740;
t4 = t102 * t121 + t103 * t124 + t118 * t538 + t271 * t60 + t272 * t62 - t58 * t739;
t3 = -t102 * t120 - t103 * t123 + t116 * t538 + t271 * t61 + t272 * t63 - t59 * t739;
t2 = t16 * t475 + t324 * t40 + t325 * t41 + t421 * t6 - t422 * t5 - t645 * t68;
t1 = t15 * t475 - t3 * t422 + t324 * t42 + t325 * t43 + t4 * t421 - t645 * t69;
t18 = [t570 + t776 / 0.2e1 + t777 / 0.2e1 - t792 / 0.2e1 + t791 / 0.2e1 + (t66 + (t41 + t891) * t422 + (t40 + t890) * t421) * t802 + t15 * t804 + t69 * t806 + t68 * t807 + (-t858 + (t43 - t890) * t422 + (t42 + t891) * t421 + t10) * t805 + (t16 + t11) * t803 + ((t424 - t883) * t499 + (t493 * t884 + t494 * t882 + t425) * t496 + (t845 * t496 + t499 * t881 - t573) * qJD(2)) * qJD(1) + (t12 * (t619 - t735) + t44 * (t557 - t602 - t795) + t13 * (t496 * t659 + t387 + t537 + t732) + t45 * (t529 + t677 + t730) + t569 * t780 + ((-t44 * pkin(6) + t45 * t569) * t497 + (t496 * t682 + t660) * t778) * qJD(1) - t44 * t817 - t45 * (-qJD(1) * t315 + t594 + t875) + (t682 * t779 * t499 - t44 * t743 + (-pkin(2) * t739 + t499 * t659 - t612) * t45) * qJD(2)) * m(6) + (t32 * (t537 + t236) - (t497 * t773 + t618 * t500) * qJD(2) + (t541 * t497 - t604 + t619) * t33 + (t529 - t630 + t670 - t594 + ((t569 - t784) * t497 + t232) * qJD(1)) * t73 + (t557 - t605 + t541 * t697 + (t772 * t694 - t829) * t497 - t816) * t72) * m(5) + (t56 * (t621 + t237) + (-t608 - t707 + t540 * t697 + (t771 * t694 - t829) * t497) * t98 - t635 * t693 + (t540 * t497 + t491 + t607) * t57 + (-t630 + t662 + t671 - t661 + t98 + ((t569 - t782) * t497 + t233) * qJD(1)) * t99) * m(4) + (-(-qJD(1) * t371 - t195 - t431 - t655) * t196 + t129 * (t497 * t649 + t491 + t706) + t128 * (t372 + t705) + t196 * (t482 + t708) + (t441 * t751 - t749) * qJD(2) + ((-pkin(1) - t609) * t750 + (t195 * (-rSges(3,3) - pkin(6)) + t196 * t649) * t497) * qJD(1)) * m(3) + (((t162 - t310 + (t359 + t745) * t500 + t728) * t500 + (t727 + t860) * t497 + t889) * qJD(2) + t867) * t637 + (t838 + t839) * t695 / 0.2e1 + (t837 - t840 + t841) * t638 + (((t500 * t631 + t164 - t727 + t892) * t500 + (t497 * t631 + t632 + t836 - t861 + t917) * t497) * qJD(2) + t842 + t863) * t640 + ((t835 - t888) * t497 + (t834 + t887) * t500) * t685 / 0.2e1; t596 * t614 + t14 * t688 / 0.2e1 + (qJD(1) * t595 + t497 * t9 - t500 * t8) * t800 + ((-t118 * t496 + t121 * t375 + t124 * t376 + t178 * t499 + t180 * t827 + t182 * t828) * t421 - (-t116 * t496 - t120 * t375 - t123 * t376 + t177 * t499 + t179 * t827 + t181 * t828) * t422 + (-t201 * t496 + t202 * t499 + t203 * t375 + t204 * t827 + t205 * t376 + t206 * t828) * t475 + (-t74 * t496 - t499 * t595) * qJD(5)) * t801 + ((-t121 * t328 - t124 * t329 + t180 * t267 + t182 * t268) * t421 - (t120 * t328 + t123 * t329 + t179 * t267 + t181 * t268) * t422 + (-t203 * t328 + t204 * t267 - t205 * t329 + t206 * t268) * t475 + (-t41 * t737 - t496 * t68) * qJD(5) + (-t40 * t687 + t504) * t497) * t802 + (qJD(1) * t599 + t497 * t6 - t5 * t500) * t803 + (qJD(1) * t597 - t3 * t500 + t4 * t497) * t804 + ((-t121 * t330 - t124 * t331 + t180 * t271 + t182 * t272) * t421 - (t120 * t330 + t123 * t331 + t179 * t271 + t181 * t272) * t422 + (-t203 * t330 + t204 * t271 - t205 * t331 + t206 * t272) * t475 + (-t42 * t738 - t496 * t69) * qJD(5) + (-t43 * t687 + t504) * t500) * t805 + t598 * t806 + t600 * t807 + (t10 * t497 + t11 * t500) * t687 / 0.2e1 + (t44 * t724 + t7 * t624 + t39 * t567 + (t12 * t627 + t44 * t629 + t7 * (t127 + t317) + t39 * (t192 + t64) + (t39 * t735 + t45 * t627) * qJD(1)) * t500 + (t13 * t627 + t45 * t629 + t7 * t735 + t39 * (t193 + t65) + (t44 * (t207 + t428) + t39 * (-t127 + t672)) * qJD(1)) * t497 - t44 * (-t183 * t475 - t208 * t422 + t568) - t45 * (t184 * t475 - t208 * t421 + t535) - t39 * (t183 * t421 + t184 * t422 + t564) - (t44 * (t743 + t725) - t45 * t390) * qJD(1) - ((t125 * t44 - t127 * t45) * t496 + t509 * t499) * qJD(5) - ((-t39 * t390 + t44 * t663) * t500 + (-t39 * t743 + t45 * t663) * t497) * qJD(2)) * m(6) + (-t72 * t568 - t73 * t535 - t67 * t564 - (t72 * (t305 + t725) - t73 * t307) * qJD(1) - ((-t67 * t307 + t666 * t72) * t500 + (-t67 * t305 + t666 * t73) * t497) * qJD(2) + t72 * t724 + t17 * t624 + t67 * t567 + (t33 * t667 + t72 * t669 + t17 * t236 + t67 * t157 + (t67 * t232 + t618) * qJD(1)) * t500 + (t32 * t667 + t73 * t669 + t17 * t232 + t67 * t159 + (t67 * t674 + t773) * qJD(1)) * t497) * m(5) + (t98 * t419 + t34 * t718 + t87 * t675 + (t57 * t723 + t98 * t726 + t34 * t237 + t87 * t158 + (t87 * t233 + t635) * qJD(1)) * t500 + (t56 * t723 + t99 * t726 + t34 * t233 + t87 * t160 + (t98 * t355 + t734 * t87) * qJD(1)) * t497 - t98 * (t465 + (t306 + t412) * qJD(1)) - t99 * (-qJD(1) * t308 + t716) - t87 * t665 - ((-t87 * t308 + t722 * t98) * t500 + (-t87 * t306 + t722 * t99) * t497) * qJD(2)) * m(4) + (-(t195 * t413 - t749) * qJD(1) - (t191 * (-t413 * t497 - t416 * t500) + t584 * t609) * qJD(2) + 0.2e1 * t191 * (t263 * t500 + t264 * t497 + (t371 * t500 - t372 * t497) * qJD(1)) + t584 * t426 + (-t128 * t497 - t129 * t500 + (-t196 * t500 + t751) * qJD(1)) * t441) * m(3) - ((((-t721 - t869) * t500 + (t720 - t870) * t497) * qJD(2) + t877) * t499 + ((t532 + (t493 * t823 + t494 * t821 - t849) * t500 + (-t493 * t822 - t494 * t820 + t848) * t497) * qJD(2) + (t493 * t819 + t494 * t818 + t878) * qJD(1)) * t496) * qJD(1) / 0.2e1 + (t840 * t500 + t839 * t497 + (t497 * t835 + t500 * t834) * qJD(1)) * qJD(1) / 0.2e1 + (t700 * t497 + (-t397 * t822 - t398 * t820 - t497 * t741) * t695 + (t397 * t819 + t398 * t818) * qJD(1) + ((t397 * t823 + t398 * t821 + t497 * t742 + t814) * qJD(2) + t815) * t500) * t640 + (-t700 * t500 + (t395 * t823 + t396 * t821 - t500 * t742) * t693 + (t395 * t819 + t396 * t818) * qJD(1) + ((-t395 * t822 - t396 * t820 + t500 * t741 + t814) * qJD(2) + t815) * t497) * t637 + (t1 + t838 * qJD(1) + ((t851 * t287 - t288 * t894 + t856 * t397 + t852 * t398 + t849 * t538 + t854 * t739) * t500 + (-t850 * t287 + t846 * t288 + t857 * t397 + t853 * t398 + t497 * t813 - t848 * t538 + t855 * t739 - t873) * t497 + (t836 * t497 + t500 * t895) * qJD(1)) * t830) * t497 / 0.2e1 - (t2 + t837 * qJD(1) + ((-t851 * t289 - t290 * t894 + t856 * t395 + t852 * t396 - t849 * t403 + t854 * t740 + t873) * t500 + (t850 * t289 + t846 * t290 + t857 * t395 + t853 * t396 + t848 * t403 - t500 * t813 + t855 * t740) * t497 + (t497 * t897 + t500 * t896) * qJD(1)) * t830) * t500 / 0.2e1 + (t10 + t842 + t865) * t643 + (t11 + t841 + t866) * t697 / 0.2e1; -m(4) * (t402 * t87 + t403 * t99 - t538 * t98) - m(5) * (t402 * t67 + t403 * t73 - t538 * t72) - m(6) * (t39 * t402 + t403 * t45 - t44 * t538) + 0.2e1 * ((t693 * t98 + t695 * t99 - t34) * t810 + (t693 * t72 + t695 * t73 - t17) * t809 + (t44 * t693 + t45 * t695 - t7) * t808) * t499 + 0.2e1 * ((qJD(2) * t87 + t497 * t56 + t500 * t57 + t697 * t99 - t699 * t98) * t810 + (qJD(2) * t67 + t32 * t497 + t33 * t500 + t697 * t73 - t699 * t72) * t809 + (qJD(2) * t39 + t12 * t500 + t13 * t497 - t44 * t699 + t45 * t697) * t808) * t496; (t12 * t397 + t13 * t395 + (t39 * t694 + t496 * t7) * t493 - t241 * t39 + t824 * t45) * m(6) + (t32 * t395 + t33 * t397 + (t17 * t496 + t67 * t694) * t493 - t241 * t67 + t824 * t73) * m(5); t1 * t650 + (-t496 * t597 + t499 * t69) * t806 + ((-qJD(2) * t597 + t15) * t499 + (qJD(1) * t598 - qJD(2) * t69 - t3 * t497 - t4 * t500) * t496) * t804 - t2 * t740 / 0.2e1 + (-t496 * t599 + t499 * t68) * t807 + ((-qJD(2) * t599 + t16) * t499 + (qJD(1) * t600 - qJD(2) * t68 - t497 * t5 - t500 * t6) * t496) * t803 + t14 * t641 + t499 * (t570 + t776 + t777 + t791 - t792) / 0.2e1 + (-t496 * t595 + t499 * t74) * t614 + ((-qJD(2) * t595 + t19) * t499 + (qJD(1) * t596 - qJD(2) * t74 - t497 * t8 - t500 * t9) * t496) * t800 + (t271 * t518 + t272 * t519 + t500 * t530) * t805 + (t267 * t518 + t268 * t519 + t497 * t530) * t802 + (-t499 * t536 + t518 * t827 + t519 * t828) * t801 + (t496 * t643 + t499 * t638) * t11 + (qJD(1) * t650 + t499 * t640) * t10 + ((qJD(2) * t509 - t12 * t125 + t13 * t127 - t44 * t65 + t45 * t64) * t499 + (t44 * (qJD(2) * t125 - t115 * t497) + t45 * (-qJD(2) * t127 + t115 * t500) + t7 * t585 + t39 * (t125 * t699 + t127 * t697 + t497 * t64 - t500 * t65) + (-t780 + t13 * t500 + (-t45 * t497 - t778) * qJD(1)) * t207) * t496 - t44 * (-t140 * t475 - t247 * t422) - t45 * (t141 * t475 - t247 * t421) - t39 * (t140 * t421 + t141 * t422)) * m(6);];
tauc = t18(:);

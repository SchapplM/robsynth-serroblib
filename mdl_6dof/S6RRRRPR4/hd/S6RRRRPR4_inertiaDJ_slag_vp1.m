% Calculate time derivative of joint inertia matrix for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:06:14
% EndTime: 2019-03-09 22:07:32
% DurationCPUTime: 51.10s
% Computational Cost: add. (100640->1358), mult. (91588->1790), div. (0->0), fcn. (86965->12), ass. (0->671)
t585 = qJ(4) + pkin(11);
t573 = sin(t585);
t574 = cos(t585);
t595 = cos(qJ(1));
t588 = qJ(2) + qJ(3);
t577 = cos(t588);
t592 = sin(qJ(1));
t838 = t577 * t592;
t486 = -t573 * t838 - t574 * t595;
t487 = -t573 * t595 + t574 * t838;
t576 = sin(t588);
t841 = t576 * t592;
t360 = Icges(6,5) * t487 + Icges(6,6) * t486 + Icges(6,3) * t841;
t590 = sin(qJ(4));
t833 = t590 * t592;
t774 = t577 * t833;
t593 = cos(qJ(4));
t830 = t593 * t595;
t512 = -t774 - t830;
t831 = t592 * t593;
t832 = t590 * t595;
t513 = t577 * t831 - t832;
t403 = Icges(5,5) * t513 + Icges(5,6) * t512 + Icges(5,3) * t841;
t955 = t360 + t403;
t837 = t577 * t595;
t488 = -t573 * t837 + t574 * t592;
t489 = t573 * t592 + t574 * t837;
t840 = t576 * t595;
t361 = Icges(6,5) * t489 + Icges(6,6) * t488 + Icges(6,3) * t840;
t514 = -t577 * t832 + t831;
t515 = t577 * t830 + t833;
t404 = Icges(5,5) * t515 + Icges(5,6) * t514 + Icges(5,3) * t840;
t954 = t361 + t404;
t362 = Icges(6,4) * t487 + Icges(6,2) * t486 + Icges(6,6) * t841;
t364 = Icges(6,1) * t487 + Icges(6,4) * t486 + Icges(6,5) * t841;
t405 = Icges(5,4) * t513 + Icges(5,2) * t512 + Icges(5,6) * t841;
t407 = Icges(5,1) * t513 + Icges(5,4) * t512 + Icges(5,5) * t841;
t953 = t362 * t488 + t364 * t489 + t405 * t514 + t407 * t515 + t840 * t955;
t363 = Icges(6,4) * t489 + Icges(6,2) * t488 + Icges(6,6) * t840;
t365 = Icges(6,1) * t489 + Icges(6,4) * t488 + Icges(6,5) * t840;
t406 = Icges(5,4) * t515 + Icges(5,2) * t514 + Icges(5,6) * t840;
t408 = Icges(5,1) * t515 + Icges(5,4) * t514 + Icges(5,5) * t840;
t952 = t363 * t488 + t365 * t489 + t406 * t514 + t408 * t515 + t840 * t954;
t733 = -qJD(4) * t577 + qJD(1);
t652 = t592 * t733;
t796 = qJD(1) * t577;
t732 = -qJD(4) + t796;
t649 = t732 * t595;
t584 = qJD(2) + qJD(3);
t836 = t584 * t592;
t778 = t576 * t836;
t926 = t649 - t778;
t339 = -t573 * t926 + t574 * t652;
t340 = t573 * t652 + t574 * t926;
t776 = t577 * t836;
t794 = qJD(1) * t595;
t620 = t576 * t794 + t776;
t226 = Icges(6,5) * t340 + Icges(6,6) * t339 + Icges(6,3) * t620;
t228 = Icges(6,4) * t340 + Icges(6,2) * t339 + Icges(6,6) * t620;
t230 = Icges(6,1) * t340 + Icges(6,4) * t339 + Icges(6,5) * t620;
t380 = -t590 * t926 + t593 * t652;
t650 = t733 * t590;
t835 = t584 * t593;
t381 = t593 * t649 + (-t576 * t835 + t650) * t592;
t257 = Icges(5,5) * t381 + Icges(5,6) * t380 + Icges(5,3) * t620;
t259 = Icges(5,4) * t381 + Icges(5,2) * t380 + Icges(5,6) * t620;
t261 = Icges(5,1) * t381 + Icges(5,4) * t380 + Icges(5,5) * t620;
t834 = t584 * t595;
t777 = t576 * t834;
t610 = t592 * t732 + t777;
t651 = t595 * t733;
t337 = t573 * t610 + t574 * t651;
t338 = t573 * t651 - t574 * t610;
t378 = t590 * t610 + t593 * t651;
t379 = -t593 * t610 + t595 * t650;
t795 = qJD(1) * t592;
t756 = t576 * t795;
t775 = t577 * t834;
t619 = -t756 + t775;
t951 = t228 * t488 + t230 * t489 + t259 * t514 + t261 * t515 + t337 * t362 + t338 * t364 + t378 * t405 + t379 * t407 + (t226 + t257) * t840 + t955 * t619;
t225 = Icges(6,5) * t338 + Icges(6,6) * t337 + Icges(6,3) * t619;
t227 = Icges(6,4) * t338 + Icges(6,2) * t337 + Icges(6,6) * t619;
t229 = Icges(6,1) * t338 + Icges(6,4) * t337 + Icges(6,5) * t619;
t256 = Icges(5,5) * t379 + Icges(5,6) * t378 + Icges(5,3) * t619;
t258 = Icges(5,4) * t379 + Icges(5,2) * t378 + Icges(5,6) * t619;
t260 = Icges(5,1) * t379 + Icges(5,4) * t378 + Icges(5,5) * t619;
t950 = t227 * t488 + t229 * t489 + t258 * t514 + t260 * t515 + t337 * t363 + t338 * t365 + t378 * t406 + t379 * t408 + (t225 + t256) * t840 + t954 * t619;
t686 = Icges(5,5) * t593 - Icges(5,6) * t590;
t839 = t577 * t584;
t366 = t686 * t839 + (Icges(5,3) * t584 + (-Icges(5,5) * t590 - Icges(5,6) * t593) * qJD(4)) * t576;
t865 = Icges(5,4) * t593;
t691 = -Icges(5,2) * t590 + t865;
t866 = Icges(5,4) * t590;
t367 = t691 * t839 + (Icges(5,6) * t584 + (-Icges(5,2) * t593 - t866) * qJD(4)) * t576;
t697 = Icges(5,1) * t593 - t866;
t368 = t697 * t839 + (Icges(5,5) * t584 + (-Icges(5,1) * t590 - t865) * qJD(4)) * t576;
t462 = -Icges(5,3) * t577 + t576 * t686;
t463 = -Icges(5,6) * t577 + t576 * t691;
t464 = -Icges(5,5) * t577 + t576 * t697;
t123 = t366 * t840 + t367 * t514 + t368 * t515 + t378 * t463 + t379 * t464 + t462 * t619;
t685 = Icges(6,5) * t574 - Icges(6,6) * t573;
t330 = t685 * t839 + (Icges(6,3) * t584 + (-Icges(6,5) * t573 - Icges(6,6) * t574) * qJD(4)) * t576;
t863 = Icges(6,4) * t574;
t690 = -Icges(6,2) * t573 + t863;
t864 = Icges(6,4) * t573;
t331 = t690 * t839 + (Icges(6,6) * t584 + (-Icges(6,2) * t574 - t864) * qJD(4)) * t576;
t696 = Icges(6,1) * t574 - t864;
t332 = t696 * t839 + (Icges(6,5) * t584 + (-Icges(6,1) * t573 - t863) * qJD(4)) * t576;
t447 = -Icges(6,3) * t577 + t576 * t685;
t448 = -Icges(6,6) * t577 + t576 * t690;
t449 = -Icges(6,5) * t577 + t576 * t696;
t96 = t330 * t840 + t331 * t488 + t332 * t489 + t337 * t448 + t338 * t449 + t447 * t619;
t949 = -t123 - t96;
t124 = t366 * t841 + t367 * t512 + t368 * t513 + t380 * t463 + t381 * t464 + t462 * t620;
t97 = t330 * t841 + t331 * t486 + t332 * t487 + t339 * t448 + t340 * t449 + t447 * t620;
t948 = -t124 - t97;
t246 = t447 * t841 + t448 * t486 + t449 * t487;
t272 = t462 * t841 + t463 * t512 + t464 * t513;
t947 = t246 + t272;
t247 = t447 * t840 + t448 * t488 + t449 * t489;
t273 = t462 * t840 + t463 * t514 + t464 * t515;
t946 = t247 + t273;
t945 = t447 + t462;
t944 = t592 * t953 + t595 * t952;
t940 = -t592 * t952 + t595 * t953;
t194 = t403 * t841 + t405 * t512 + t407 * t513;
t195 = t404 * t841 + t406 * t512 + t408 * t513;
t672 = t194 * t592 + t195 * t595;
t169 = t360 * t841 + t362 * t486 + t364 * t487;
t170 = t361 * t841 + t363 * t486 + t365 * t487;
t678 = t169 * t592 + t170 * t595;
t943 = t672 + t678;
t939 = (t169 + t194) * t595 + (-t170 - t195) * t592;
t575 = qJ(6) + t585;
t567 = sin(t575);
t568 = cos(t575);
t467 = -t567 * t838 - t568 * t595;
t468 = -t567 * t595 + t568 * t838;
t344 = Icges(7,5) * t468 + Icges(7,6) * t467 + Icges(7,3) * t841;
t346 = Icges(7,4) * t468 + Icges(7,2) * t467 + Icges(7,6) * t841;
t348 = Icges(7,1) * t468 + Icges(7,4) * t467 + Icges(7,5) * t841;
t669 = -t346 * t567 + t348 * t568;
t179 = -t344 * t577 + t576 * t669;
t684 = Icges(7,5) * t568 - Icges(7,6) * t567;
t442 = -Icges(7,3) * t577 + t576 * t684;
t861 = Icges(7,4) * t568;
t689 = -Icges(7,2) * t567 + t861;
t443 = -Icges(7,6) * t577 + t576 * t689;
t862 = Icges(7,4) * t567;
t695 = Icges(7,1) * t568 - t862;
t444 = -Icges(7,5) * t577 + t576 * t695;
t233 = t442 * t841 + t443 * t467 + t444 * t468;
t942 = -t233 - t179;
t470 = -t567 * t837 + t568 * t592;
t471 = t567 * t592 + t568 * t837;
t345 = Icges(7,5) * t471 + Icges(7,6) * t470 + Icges(7,3) * t840;
t347 = Icges(7,4) * t471 + Icges(7,2) * t470 + Icges(7,6) * t840;
t349 = Icges(7,1) * t471 + Icges(7,4) * t470 + Icges(7,5) * t840;
t668 = -t347 * t567 + t349 * t568;
t180 = -t345 * t577 + t576 * t668;
t234 = t442 * t840 + t443 * t470 + t444 * t471;
t941 = -t234 - t180;
t666 = -t363 * t573 + t365 * t574;
t189 = -t361 * t577 + t576 * t666;
t664 = -t406 * t590 + t408 * t593;
t220 = -t404 * t577 + t576 * t664;
t826 = t189 + t220;
t667 = -t362 * t573 + t364 * t574;
t188 = -t360 * t577 + t576 * t667;
t665 = -t405 * t590 + t407 * t593;
t219 = -t403 * t577 + t576 * t665;
t827 = t188 + t219;
t938 = t592 * t827 + t595 * t826;
t937 = (t584 * t944 + t949) * t577 + (qJD(1) * t940 + t584 * t946 + t592 * t951 + t595 * t950) * t576;
t54 = t226 * t841 + t228 * t486 + t230 * t487 + t339 * t362 + t340 * t364 + t360 * t620;
t55 = t225 * t841 + t227 * t486 + t229 * t487 + t339 * t363 + t340 * t365 + t361 * t620;
t64 = t257 * t841 + t259 * t512 + t261 * t513 + t380 * t405 + t381 * t407 + t403 * t620;
t65 = t256 * t841 + t258 * t512 + t260 * t513 + t380 * t406 + t381 * t408 + t404 * t620;
t936 = (t584 * t943 + t948) * t577 + ((t55 + t65) * t595 + (t54 + t64) * t592 + t947 * t584 + t939 * qJD(1)) * t576;
t935 = qJD(1) * t944 + t592 * t950 - t595 * t951;
t26 = qJD(1) * t678 - t54 * t595 + t55 * t592;
t33 = qJD(1) * t672 + t592 * t65 - t595 * t64;
t934 = t33 + t26;
t933 = t576 * t943 - t577 * t947;
t932 = t576 * t944 - t577 * t946;
t583 = qJD(4) + qJD(6);
t735 = -t577 * t583 + qJD(1);
t654 = t595 * t735;
t734 = -t583 + t796;
t908 = t592 * t734 + t777;
t311 = t567 * t908 + t568 * t654;
t312 = t567 * t654 - t568 * t908;
t769 = rSges(7,1) * t312 + rSges(7,2) * t311 + rSges(7,3) * t775;
t210 = -rSges(7,3) * t756 + t769;
t789 = qJD(4) * t590;
t783 = pkin(4) * t789;
t589 = -qJ(5) - pkin(9);
t582 = -pkin(10) + t589;
t800 = t582 - t589;
t616 = -t584 * t800 + t783;
t890 = pkin(4) * t590;
t542 = pkin(5) * t573 + t890;
t519 = t542 * qJD(4);
t580 = t593 * pkin(4);
t889 = pkin(5) * t574;
t520 = (t580 + t889) * qJD(4);
t724 = -t519 * t837 + t520 * t592 + t542 * t794 + t582 * t756;
t566 = pkin(4) * t832;
t788 = qJD(4) * t593;
t782 = pkin(4) * t788;
t757 = qJD(1) * t566 + t589 * t756 + t592 * t782;
t569 = t580 + pkin(3);
t536 = t569 + t889;
t805 = t536 - t569;
t929 = -t805 * t777 + (t595 * t616 - t795 * t805) * t577 + t724 - t757 + t210;
t851 = t463 * t590;
t852 = t449 * t574;
t928 = t945 * t577 + (t448 * t573 - t464 * t593 + t851 - t852) * t576;
t742 = t805 * t577;
t626 = -t576 * t582 + t742;
t804 = t589 * t841 + t566;
t301 = -t542 * t595 + t592 * t626 + t804;
t703 = -rSges(7,1) * t468 - rSges(7,2) * t467;
t350 = rSges(7,3) * t841 - t703;
t927 = t301 + t350;
t749 = t839 / 0.2e1;
t925 = -t595 * t749 + t756 / 0.2e1;
t746 = t794 / 0.2e1;
t924 = -t576 * t746 - t592 * t749;
t596 = -pkin(8) - pkin(7);
t591 = sin(qJ(2));
t792 = qJD(2) * t591;
t785 = pkin(2) * t792;
t923 = qJD(1) * t596 + t785;
t594 = cos(qJ(2));
t548 = rSges(3,1) * t591 + rSges(3,2) * t594;
t631 = qJD(2) * t548;
t922 = t592 * t631;
t869 = Icges(3,4) * t594;
t694 = -Icges(3,2) * t591 + t869;
t501 = Icges(3,6) * t592 + t595 * t694;
t870 = Icges(3,4) * t591;
t700 = Icges(3,1) * t594 - t870;
t503 = Icges(3,5) * t592 + t595 * t700;
t657 = t501 * t591 - t503 * t594;
t921 = t592 * t657;
t867 = Icges(4,4) * t577;
t692 = -Icges(4,2) * t576 + t867;
t477 = Icges(4,6) * t592 + t595 * t692;
t868 = Icges(4,4) * t576;
t698 = Icges(4,1) * t577 - t868;
t479 = Icges(4,5) * t592 + t595 * t698;
t659 = t477 * t576 - t479 * t577;
t920 = t592 * t659;
t570 = pkin(2) * t594 + pkin(1);
t886 = pkin(1) - t570;
t919 = t592 * t886;
t500 = -Icges(3,6) * t595 + t592 * t694;
t502 = -Icges(3,5) * t595 + t592 * t700;
t658 = t500 * t591 - t502 * t594;
t918 = t595 * t658;
t476 = -Icges(4,6) * t595 + t592 * t692;
t478 = -Icges(4,5) * t595 + t592 * t698;
t660 = t476 * t576 - t478 * t577;
t917 = t595 * t660;
t289 = t584 * t742 + (-t519 + t616) * t576;
t702 = rSges(7,1) * t568 - rSges(7,2) * t567;
t319 = t702 * t839 + (rSges(7,3) * t584 + (-rSges(7,1) * t567 - rSges(7,2) * t568) * t583) * t576;
t916 = -t289 - t319;
t351 = rSges(7,1) * t471 + rSges(7,2) * t470 + rSges(7,3) * t840;
t915 = -t350 * t592 - t351 * t595;
t914 = -t584 * t851 - t366;
t709 = -rSges(5,1) * t513 - rSges(5,2) * t512;
t409 = rSges(5,3) * t841 - t709;
t410 = rSges(5,1) * t515 + rSges(5,2) * t514 + rSges(5,3) * t840;
t913 = -t409 * t592 - t410 * t595;
t412 = t576 * t805 + t577 * t800;
t446 = -rSges(7,3) * t577 + t576 * t702;
t912 = -t412 - t446;
t790 = qJD(4) * t576;
t842 = t576 * t584;
t911 = t839 * t852 + (t332 * t576 - t448 * t790) * t574 + t576 * t593 * t368 + t945 * t842 + (t464 * t835 - t330) * t577;
t687 = Icges(4,5) * t577 - Icges(4,6) * t576;
t474 = -Icges(4,3) * t595 + t592 * t687;
t910 = qJD(1) * t474;
t909 = -t577 * t826 + t932;
t688 = Icges(3,5) * t594 - Icges(3,6) * t591;
t498 = -Icges(3,3) * t595 + t592 * t688;
t907 = t595 * t734 - t778;
t523 = Icges(4,2) * t577 + t868;
t524 = Icges(4,1) * t576 + t867;
t656 = t523 * t576 - t524 * t577;
t906 = qJD(1) * t656 + t584 * t687;
t905 = 2 * m(3);
t904 = 2 * m(4);
t903 = 2 * m(5);
t902 = 2 * m(6);
t901 = 2 * m(7);
t586 = t592 ^ 2;
t587 = t595 ^ 2;
t900 = m(6) / 0.2e1;
t899 = m(7) / 0.2e1;
t898 = -t577 / 0.2e1;
t897 = t592 / 0.2e1;
t896 = -t595 / 0.2e1;
t895 = -rSges(5,3) - pkin(9);
t894 = m(3) * t548;
t526 = rSges(4,1) * t576 + rSges(4,2) * t577;
t893 = m(4) * t526;
t892 = pkin(2) * t591;
t891 = pkin(3) * t577;
t888 = pkin(9) * t576;
t887 = t592 * pkin(7);
t581 = t595 * pkin(7);
t885 = pkin(3) - t569;
t884 = -pkin(7) - t596;
t883 = pkin(9) + t589;
t882 = rSges(3,1) * t594;
t881 = rSges(4,1) * t577;
t880 = rSges(3,2) * t591;
t879 = rSges(3,3) * t595;
t655 = t592 * t735;
t313 = -t567 * t907 + t568 * t655;
t314 = t567 * t655 + t568 * t907;
t205 = Icges(7,5) * t314 + Icges(7,6) * t313 + Icges(7,3) * t620;
t207 = Icges(7,4) * t314 + Icges(7,2) * t313 + Icges(7,6) * t620;
t209 = Icges(7,1) * t314 + Icges(7,4) * t313 + Icges(7,5) * t620;
t50 = (t584 * t669 - t205) * t577 + (t344 * t584 + (-t346 * t583 + t209) * t568 + (-t348 * t583 - t207) * t567) * t576;
t878 = t50 * t595;
t204 = Icges(7,5) * t312 + Icges(7,6) * t311 + Icges(7,3) * t619;
t206 = Icges(7,4) * t312 + Icges(7,2) * t311 + Icges(7,6) * t619;
t208 = Icges(7,1) * t312 + Icges(7,4) * t311 + Icges(7,5) * t619;
t51 = (t584 * t668 - t204) * t577 + (t345 * t584 + (-t347 * t583 + t208) * t568 + (-t349 * t583 - t206) * t567) * t576;
t877 = t51 * t592;
t579 = t592 * rSges(3,3);
t578 = t592 * rSges(4,3);
t60 = (t584 * t667 - t226) * t577 + (-t228 * t573 + t230 * t574 + t360 * t584 + (-t362 * t574 - t364 * t573) * qJD(4)) * t576;
t876 = t60 * t595;
t61 = (t584 * t666 - t225) * t577 + (-t227 * t573 + t229 * t574 + t361 * t584 + (-t363 * t574 - t365 * t573) * qJD(4)) * t576;
t875 = t61 * t592;
t70 = (t584 * t665 - t257) * t577 + (-t259 * t590 + t261 * t593 + t403 * t584 + (-t405 * t593 - t407 * t590) * qJD(4)) * t576;
t874 = t70 * t595;
t71 = (t584 * t664 - t256) * t577 + (-t258 * t590 + t260 * t593 + t404 * t584 + (-t406 * t593 - t408 * t590) * qJD(4)) * t576;
t873 = t71 * t592;
t872 = -rSges(7,3) + t582;
t854 = t367 * t590;
t853 = t444 * t568;
t711 = -rSges(4,2) * t576 + t881;
t495 = t711 * t584;
t850 = t495 * t592;
t849 = t500 * t594;
t848 = t501 * t594;
t847 = t502 * t591;
t846 = t503 * t591;
t845 = t523 * t584;
t844 = t524 * t584;
t843 = t576 * t583;
t829 = t595 * t596;
t780 = t448 * t839;
t828 = (-t780 + (-qJD(4) * t449 - t331) * t576) * t573 + t914 * t577 + (-t854 + (-t463 * t593 - t464 * t590) * qJD(4)) * t576 + t911;
t704 = rSges(7,1) * t314 + rSges(7,2) * t313;
t211 = rSges(7,3) * t620 + t704;
t825 = t211 * t840 + t350 * t775;
t768 = rSges(6,1) * t338 + rSges(6,2) * t337 + rSges(6,3) * t775;
t231 = -rSges(6,3) * t756 + t768;
t539 = pkin(9) * t775;
t640 = -t584 * t589 - t783;
t787 = qJD(5) * t576;
t549 = t595 * t787;
t723 = t549 + t757;
t248 = -t539 + (pkin(9) * t795 + t834 * t885) * t576 + (t595 * t640 + t795 * t885) * t577 + t723;
t824 = -t231 - t248;
t538 = pkin(3) * t778;
t565 = pkin(4) * t833;
t748 = t885 * t577;
t638 = -t748 - t888;
t701 = qJD(4) * pkin(4) * t774 + t569 * t778 + t589 * t620 + t595 * t782;
t249 = t538 + (-pkin(9) * t839 + t787) * t592 + (t595 * t638 + t565) * qJD(1) - t701;
t399 = t592 * t638 - t804;
t823 = t249 * t840 + t399 * t775;
t822 = t928 * t842;
t439 = t446 * t795;
t821 = t351 * t842 + t439 * t576;
t705 = rSges(6,1) * t574 - rSges(6,2) * t573;
t333 = t705 * t839 + (rSges(6,3) * t584 + (-rSges(6,1) * t573 - rSges(6,2) * t574) * qJD(4)) * t576;
t753 = t576 * t789;
t343 = -pkin(4) * t753 - qJD(5) * t577 + (-t576 * t883 - t748) * t584;
t820 = -t333 - t343;
t556 = pkin(3) * t837;
t509 = pkin(9) * t840 + t556;
t806 = -t569 * t837 - t565;
t645 = -t589 * t840 - t806;
t400 = t645 - t509;
t445 = -t576 * t885 + t577 * t883;
t438 = t445 * t795;
t819 = t400 * t842 + t438 * t576;
t706 = -rSges(6,1) * t487 - rSges(6,2) * t486;
t369 = rSges(6,3) * t841 - t706;
t818 = -t369 - t399;
t370 = rSges(6,1) * t489 + rSges(6,2) * t488 + rSges(6,3) * t840;
t817 = -t370 - t400;
t708 = rSges(5,1) * t593 - rSges(5,2) * t590;
t371 = t708 * t839 + (rSges(5,3) * t584 + (-rSges(5,1) * t590 - rSges(5,2) * t593) * qJD(4)) * t576;
t714 = t888 + t891;
t497 = t714 * t584;
t816 = -t371 - t497;
t815 = t399 * t577 + t445 * t841;
t814 = -t410 - t509;
t277 = t350 * t577 + t446 * t841;
t527 = pkin(3) * t576 - pkin(9) * t577;
t510 = t527 * t795;
t813 = t438 + t510;
t451 = -rSges(6,3) * t577 + t576 * t705;
t812 = -t445 - t451;
t465 = -rSges(5,3) * t577 + t576 * t708;
t450 = t465 * t795;
t811 = t450 + t510;
t472 = t581 + t829 - t919;
t550 = t595 * t570;
t473 = -pkin(1) * t595 + t592 * t884 + t550;
t810 = t472 * t592 + t473 * t595;
t482 = -rSges(4,3) * t595 + t592 * t711;
t483 = rSges(4,1) * t837 - rSges(4,2) * t840 + t578;
t382 = t482 * t592 + t483 * t595;
t809 = -t465 - t527;
t508 = t714 * t592;
t808 = t508 * t592 + t509 * t595;
t807 = t536 * t837 + t592 * t542;
t803 = rSges(4,2) * t756 + rSges(4,3) * t794;
t802 = t923 * t592;
t801 = t595 * t882 + t579;
t799 = t586 + t587;
t475 = Icges(4,3) * t592 + t595 * t687;
t798 = qJD(1) * t475;
t499 = Icges(3,3) * t592 + t595 * t688;
t797 = qJD(1) * t499;
t791 = qJD(2) * t594;
t786 = t595 * t880;
t784 = pkin(2) * t791;
t781 = t443 * t839;
t773 = -t248 - t929;
t772 = -t343 + t916;
t771 = -t399 - t927;
t302 = -t800 * t840 + t806 + t807;
t770 = -t302 - t351 - t400;
t767 = -t497 + t820;
t765 = -t509 + t817;
t621 = -t577 * t795 - t777;
t633 = t526 * t584;
t764 = t592 * (-t592 * t633 + (t595 * t711 + t578) * qJD(1)) + t595 * (rSges(4,1) * t621 - rSges(4,2) * t775 + t803) + t482 * t794;
t763 = rSges(5,1) * t379 + rSges(5,2) * t378 + rSges(5,3) * t775;
t762 = t592 * (pkin(9) * t620 + qJD(1) * t556 - t538) + t595 * (pkin(3) * t621 - pkin(9) * t756 + t539) + t508 * t794;
t761 = -t445 + t912;
t760 = t592 * ((-t595 * t886 - t887) * qJD(1) - t802) + t595 * (-t595 * t785 + (t595 * t884 + t919) * qJD(1)) + t472 * t794;
t441 = t451 * t795;
t759 = t441 + t813;
t758 = -t527 + t812;
t754 = t591 * t795;
t396 = t412 * t795;
t751 = t841 / 0.2e1;
t750 = t840 / 0.2e1;
t747 = t795 / 0.2e1;
t745 = -t526 - t892;
t744 = -t527 - t892;
t743 = t812 * t595;
t415 = t809 * t595;
t385 = -qJD(1) * t476 - t595 * t845;
t741 = t479 * t584 + t385;
t386 = qJD(1) * t477 - t592 * t845;
t740 = t478 * t584 + t386;
t387 = -qJD(1) * t478 - t595 * t844;
t739 = -t477 * t584 + t387;
t388 = qJD(1) * t479 - t592 * t844;
t738 = t476 * t584 - t388;
t737 = -t536 * t577 - t570;
t736 = -t592 * t596 + t550;
t731 = t577 * t211 + t319 * t841 + t446 * t620;
t730 = t577 * t249 + t343 * t841 + t445 * t620;
t729 = -t497 + t772;
t728 = -t509 + t770;
t727 = t399 * t592 + t400 * t595 + t808;
t726 = t396 + t439 + t813;
t725 = -t527 + t761;
t243 = t808 - t913;
t718 = -t445 + t744;
t717 = -t465 + t744;
t716 = -t497 - t784;
t715 = t595 * t761;
t294 = t758 * t595;
t268 = -t442 * t577 + (-t443 * t567 + t853) * t576;
t41 = t205 * t840 + t207 * t470 + t209 * t471 + t311 * t346 + t312 * t348 + t344 * t619;
t42 = t204 * t840 + t206 * t470 + t208 * t471 + t311 * t347 + t312 * t349 + t345 * t619;
t162 = t344 * t840 + t346 * t470 + t348 * t471;
t163 = t345 * t840 + t347 * t470 + t349 * t471;
t680 = t162 * t592 + t163 * t595;
t681 = t162 * t595 - t163 * t592;
t308 = t684 * t839 + (Icges(7,3) * t584 + (-Icges(7,5) * t567 - Icges(7,6) * t568) * t583) * t576;
t309 = t689 * t839 + (Icges(7,6) * t584 + (-Icges(7,2) * t568 - t862) * t583) * t576;
t310 = t695 * t839 + (Icges(7,5) * t584 + (-Icges(7,1) * t567 - t861) * t583) * t576;
t90 = t308 * t840 + t309 * t470 + t310 * t471 + t311 * t443 + t312 * t444 + t442 * t619;
t5 = (t584 * t680 - t90) * t577 + (qJD(1) * t681 + t234 * t584 + t41 * t592 + t42 * t595) * t576;
t43 = t205 * t841 + t207 * t467 + t209 * t468 + t313 * t346 + t314 * t348 + t344 * t620;
t44 = t204 * t841 + t206 * t467 + t208 * t468 + t313 * t347 + t314 * t349 + t345 * t620;
t160 = t344 * t841 + t346 * t467 + t348 * t468;
t161 = t345 * t841 + t347 * t467 + t349 * t468;
t682 = t160 * t592 + t161 * t595;
t683 = t160 * t595 - t161 * t592;
t91 = t308 * t841 + t309 * t467 + t310 * t468 + t313 * t443 + t314 * t444 + t442 * t620;
t6 = (t584 * t682 - t91) * t577 + (qJD(1) * t683 + t233 * t584 + t43 * t592 + t44 * t595) * t576;
t674 = t179 * t592 + t180 * t595;
t77 = -t233 * t577 + t576 * t682;
t78 = -t234 * t577 + t576 * t680;
t713 = t5 * t840 + t6 * t841 + t78 * t775 + (-t268 * t577 + t576 * t674) * t842 + t620 * t77;
t712 = -t880 + t882;
t710 = rSges(5,1) * t381 + rSges(5,2) * t380;
t707 = rSges(6,1) * t340 + rSges(6,2) * t339;
t699 = Icges(3,1) * t591 + t869;
t693 = Icges(3,2) * t594 + t870;
t522 = Icges(4,5) * t576 + Icges(4,6) * t577;
t675 = t179 * t595 - t180 * t592;
t663 = t409 * t595 - t410 * t592;
t653 = -t451 + t718;
t648 = -t343 + t716;
t647 = -t371 + t716;
t252 = t725 * t595;
t646 = -pkin(1) - t712;
t393 = t717 * t595;
t644 = -t577 * t827 + t933;
t643 = -t570 - t711;
t642 = t248 * t595 + t249 * t592 + t399 * t794 + t762;
t262 = -rSges(5,3) * t756 + t763;
t263 = rSges(5,3) * t620 + t710;
t641 = t262 * t595 + t263 * t592 + t409 * t794 + t762;
t155 = t369 * t592 + t370 * t595 + t727;
t639 = -rSges(6,3) * t576 - t569 * t577 - t570;
t636 = t188 / 0.2e1 + t219 / 0.2e1 + t246 / 0.2e1 + t272 / 0.2e1;
t635 = t189 / 0.2e1 + t220 / 0.2e1 + t247 / 0.2e1 + t273 / 0.2e1;
t634 = t718 + t912;
t632 = -t333 + t648;
t628 = t584 * t522;
t627 = t369 * t595 + t592 * t817;
t625 = qJD(2) * t699;
t624 = qJD(2) * t693;
t623 = qJD(2) * (-Icges(3,5) * t591 - Icges(3,6) * t594);
t284 = t653 * t595;
t622 = t576 * t895 - t570 - t891;
t618 = t648 + t916;
t617 = t639 * t592;
t615 = t576 * t872 + t737;
t614 = t301 * t595 + t592 * t770;
t245 = t634 * t595;
t125 = t301 * t592 + t302 * t595 + t727 - t915;
t232 = rSges(6,3) * t620 + t707;
t613 = t231 * t595 + t232 * t592 + t369 * t794 + t642;
t609 = rSges(3,2) * t754 + rSges(3,3) * t794 - t595 * t631;
t608 = -t577 * t308 + t442 * t842 + t839 * t853 + (t310 * t576 - t443 * t843) * t568;
t19 = qJD(1) * t680 - t41 * t595 + t42 * t592;
t20 = qJD(1) * t682 - t43 * t595 + t44 * t592;
t606 = t19 * t750 + t20 * t751 + t5 * t897 + t6 * t896 - t675 * t842 / 0.2e1 + (qJD(1) * t674 + t877 - t878) * t898 + t77 * t747 + t78 * t746 + t925 * t681 + t924 * t683;
t285 = -t474 * t595 - t592 * t660;
t286 = -t475 * t595 - t920;
t287 = t474 * t592 - t917;
t288 = t475 * t592 - t595 * t659;
t383 = -t595 * t628 - t910;
t384 = -t592 * t628 + t798;
t605 = (-t285 * t595 - t683 - t939) * t795 + (-t287 * t595 - t681 - t940) * t794 + (t19 + t286 * t795 + t288 * t794 + (t288 * qJD(1) + (t386 * t576 - t388 * t577 + t476 * t839 + t478 * t842 - t910) * t595) * t595 + ((t287 + t920) * qJD(1) + (-t384 + t739 * t577 - t741 * t576 + (t475 - t660) * qJD(1)) * t595 + t383 * t592) * t592 + t935) * t592;
t604 = t592 * t622 - t829;
t111 = (-t781 + (-t444 * t583 - t309) * t576) * t567 + t608;
t255 = t268 * t842;
t11 = t255 + (t584 * t674 - t111) * t577 + (qJD(1) * t675 + t50 * t592 + t51 * t595) * t576;
t603 = -t11 * t577 - t756 * t78 + t713;
t178 = -t520 * t595 + (-t536 * t842 + (-t582 * t584 - t519) * t577) * t592 + ((t542 - t890) * t592 + t626 * t595) * qJD(1) + t701;
t602 = t642 + t927 * t794 + t929 * t595 + (t178 + t211) * t592;
t601 = t255 + (t50 + t91) * t751 + (t51 + t90) * t750 + t941 * t925 + t942 * t924;
t493 = t692 * t584;
t494 = t698 * t584;
t600 = qJD(1) * t522 + (t494 - t845) * t577 + (-t493 - t844) * t576;
t57 = (t384 * t595 + (t286 + t917) * qJD(1)) * t595 + (t285 * qJD(1) + (-t385 * t576 + t387 * t577 - t477 * t839 - t479 * t842 + t798) * t592 + (-t383 + t738 * t577 + t740 * t576 + (-t474 - t659) * qJD(1)) * t595) * t592;
t599 = t605 + (-t20 - t57 - t934) * t595;
t598 = t877 / 0.2e1 - t874 / 0.2e1 + t873 / 0.2e1 - t876 / 0.2e1 - t878 / 0.2e1 + t875 / 0.2e1 + (t576 * t739 + t577 * t741 + t592 * t906 + t595 * t600 + t90 - t949) * t897 + (-t576 * t738 + t577 * t740 + t592 * t600 - t595 * t906 + t91 - t948) * t896 + (t476 * t577 + t478 * t576 - t522 * t595 - t592 * t656 + t827 - t942 + t947) * t747 + (t477 * t577 + t479 * t576 + t522 * t592 - t595 * t656 + t826 - t941 + t946) * t746;
t597 = t606 + (qJD(1) * t938 + t873 - t874 + t875 - t876) * t898 + t937 * t897 + t936 * t896 + (t592 * t826 - t595 * t827) * t842 / 0.2e1 + t934 * t751 + t935 * t750 + t933 * t747 + t932 * t746 + t940 * t925 + t939 * t924;
t561 = pkin(2) * t754;
t537 = t712 * qJD(2);
t505 = -t786 + t801;
t504 = t592 * t712 - t879;
t469 = t745 * t595;
t466 = t745 * t592;
t455 = t887 + (pkin(1) - t880) * t595 + t801;
t454 = t592 * t646 + t581 + t879;
t436 = t483 + t736;
t435 = (rSges(4,3) - t596) * t595 + t643 * t592;
t427 = t592 * t623 + t797;
t426 = -qJD(1) * t498 + t595 * t623;
t414 = t809 * t592;
t395 = t922 + ((-rSges(3,3) - pkin(7)) * t592 + t646 * t595) * qJD(1);
t394 = (t581 + (-pkin(1) - t882) * t592) * qJD(1) + t609;
t392 = t717 * t592;
t359 = t399 * t840;
t357 = -t526 * t794 - t850 + (-t591 * t794 - t592 * t791) * pkin(2);
t356 = t526 * t795 + t561 + (-t495 - t784) * t595;
t323 = t350 * t840;
t318 = t499 * t592 - t595 * t657;
t317 = t498 * t592 - t918;
t316 = -t499 * t595 - t921;
t315 = -t498 * t595 - t592 * t658;
t304 = t526 * t836 + (t595 * t643 - t578) * qJD(1) + t802;
t303 = (-t570 - t881) * t795 + (-t633 - t923) * t595 + t803;
t300 = t736 - t814;
t299 = t604 + t709;
t293 = t758 * t592;
t292 = -t410 * t577 - t465 * t840;
t291 = t409 * t577 + t465 * t841;
t283 = t653 * t592;
t281 = t645 + t736 + t370;
t280 = t617 + t706 + t804 - t829;
t278 = -t351 * t577 - t446 * t840;
t276 = t382 + t810;
t275 = t663 * t576;
t271 = -t582 * t840 + t351 + t736 + t807;
t270 = (t542 - t596) * t595 + t615 * t592 + t703;
t265 = qJD(1) * t415 + t592 * t816;
t264 = t595 * t816 + t811;
t251 = t725 * t592;
t250 = -t351 * t841 + t323;
t244 = t634 * t592;
t240 = qJD(1) * t393 + t592 * t647;
t239 = t595 * t647 + t561 + t811;
t216 = -t483 * t795 + t764;
t191 = t576 * t743 + t577 * t817;
t190 = t369 * t577 + t451 * t841 + t815;
t187 = t243 + t810;
t182 = t622 * t794 + t776 * t895 + t538 - t710 + t802;
t181 = t539 + (-pkin(3) * t842 - t785) * t595 + t604 * qJD(1) + t763;
t164 = t576 * t627 + t359;
t159 = qJD(1) * t294 + t592 * t767;
t158 = t595 * t767 + t759;
t157 = qJD(1) * t284 + t592 * t632;
t156 = t595 * t632 + t561 + t759;
t154 = (-rSges(6,3) * t839 - t787) * t592 + (t595 * t639 - t565) * qJD(1) + t701 - t707 + t802;
t153 = qJD(1) * t617 + (-t569 * t842 + t577 * t640 - t923) * t595 + t723 + t768;
t152 = (-t473 - t483) * t795 + t760 + t764;
t151 = t155 + t810;
t150 = (t465 * t836 + t263) * t577 + (t371 * t592 - t409 * t584 + t465 * t794) * t576;
t149 = (-t465 * t834 - t262) * t577 + (-t371 * t595 + t410 * t584 + t450) * t576;
t148 = t576 * t715 + t577 * t770;
t147 = t301 * t577 + t412 * t841 + t277 + t815;
t146 = (qJD(1) * t615 + t520) * t595 + (-qJD(1) * t542 + (t536 * t584 - qJD(5)) * t576 + (t584 * t872 + t519) * t577) * t592 - t704 + t802;
t145 = t549 + (-t785 + (-t536 * t576 - t577 * t582) * t584) * t595 + (-t829 + (-rSges(7,3) * t576 + t737) * t592) * qJD(1) + t724 + t769;
t134 = qJD(1) * t252 + t592 * t729;
t133 = t595 * t729 + t726;
t131 = qJD(1) * t245 + t592 * t618;
t130 = t595 * t618 + t561 + t726;
t129 = -t350 * t842 + t731;
t128 = -t319 * t840 + (-t446 * t834 - t210) * t577 + t821;
t127 = t576 * t614 + t323 + t359;
t108 = t125 + t810;
t101 = t663 * t839 + (qJD(1) * t913 - t262 * t592 + t263 * t595) * t576;
t98 = t795 * t814 + t641;
t80 = -t351 * t776 + (qJD(1) * t915 - t210 * t592) * t576 + t825;
t79 = (-t473 + t814) * t795 + t641 + t760;
t67 = (t451 * t836 + t232) * t577 + (t333 * t592 + t451 * t794 + t584 * t818) * t576 + t730;
t66 = (t584 * t743 + t824) * t577 + (t370 * t584 + t595 * t820 + t441) * t576 + t819;
t45 = t765 * t795 + t613;
t40 = (-t473 + t765) * t795 + t613 + t760;
t39 = t627 * t839 + (t232 * t595 + t824 * t592 + (t592 * t818 + t595 * t817) * qJD(1)) * t576 + t823;
t38 = (t412 * t836 + t178) * t577 + (t289 * t592 + t412 * t794 + t584 * t771) * t576 + t730 + t731;
t37 = (t302 * t584 + t595 * t772 + t396) * t576 + (t584 * t715 + t773) * t577 + t819 + t821;
t36 = t728 * t795 + t602;
t34 = t602 + (-t473 + t728) * t795 + t760;
t28 = t614 * t839 + (t178 * t595 + t773 * t592 + (t592 * t771 + t595 * t770) * qJD(1)) * t576 + t823 + t825;
t1 = [-t523 * t842 + (t145 * t271 + t146 * t270) * t901 + (t153 * t281 + t154 * t280) * t902 + (t181 * t300 + t182 * t299) * t903 + (t303 * t436 + t304 * t435) * t904 + (t394 * t455 + t395 * t454) * t905 + t608 + t524 * t839 - t464 * t753 + (-t693 + t700) * t792 + (t694 + t699) * t791 + (-t449 * t790 - t780) * t573 + (-t444 * t843 - t781) * t567 + (t493 + t914) * t577 + (-t309 * t567 - t331 * t573 - t463 * t788 + t494 - t854) * t576 + t911; (t586 / 0.2e1 + t587 / 0.2e1) * t688 * qJD(2) + m(3) * ((-t394 * t592 - t395 * t595) * t548 + (-t454 * t595 - t455 * t592) * t537) + m(7) * (t130 * t270 + t131 * t271 + t145 * t244 + t146 * t245) + m(6) * (t153 * t283 + t154 * t284 + t156 * t280 + t157 * t281) + m(5) * (t181 * t392 + t182 * t393 + t239 * t299 + t240 * t300) + m(4) * (t303 * t466 + t304 * t469 + t356 * t435 + t357 * t436) + (-qJD(2) * t658 + (qJD(1) * t501 - t592 * t624) * t594 + (qJD(1) * t503 - t592 * t625) * t591) * t896 + (-qJD(2) * t657 + (-qJD(1) * t500 - t595 * t624) * t594 + (-qJD(1) * t502 - t595 * t625) * t591) * t897 + ((t848 / 0.2e1 + t846 / 0.2e1 - t455 * t894) * t595 + (t454 * t894 + t849 / 0.2e1 + t847 / 0.2e1) * t592) * qJD(1) + t598; t605 + (-t315 * t595 + t316 * t592) * t795 + (-t317 * t595 + t318 * t592) * t794 + t592 * ((t426 * t592 + (t317 + t921) * qJD(1)) * t592 + (t318 * qJD(1) + (t500 * t791 + t502 * t792) * t595 + (-t427 + (-t846 - t848) * qJD(2) + (t499 - t658) * qJD(1)) * t592) * t595) - t595 * t20 - t595 * t26 - t595 * t33 - t595 * t57 + ((t504 * t592 + t505 * t595) * ((qJD(1) * t504 + t609) * t595 + (-t922 + (-t505 - t786 + t579) * qJD(1)) * t592) + t799 * t548 * t537) * t905 + (t108 * t34 + t130 * t245 + t131 * t244) * t901 + (t151 * t40 + t156 * t284 + t157 * t283) * t902 + (t187 * t79 + t239 * t393 + t240 * t392) * t903 + (t152 * t276 + t356 * t469 + t357 * t466) * t904 - t595 * ((t427 * t595 + (t316 + t918) * qJD(1)) * t595 + (t315 * qJD(1) + (-t501 * t791 - t503 * t792 + t797) * t592 + (-t426 + (t847 + t849) * qJD(2) - t657 * qJD(1)) * t595) * t592); m(7) * (t133 * t270 + t134 * t271 + t145 * t251 + t146 * t252) + m(6) * (t153 * t293 + t154 * t294 + t158 * t280 + t159 * t281) + m(5) * (t181 * t414 + t182 * t415 + t264 * t299 + t265 * t300) + (-t303 * t592 - t304 * t595 + (t435 * t592 - t436 * t595) * qJD(1)) * t893 + m(4) * (-t435 * t595 - t436 * t592) * t495 + t598; t599 + (-t356 * t595 - t357 * t592 + (-t466 * t595 + t469 * t592) * qJD(1)) * t893 + m(4) * (-t469 * t495 * t595 + t152 * t382 + t216 * t276 - t466 * t850) + m(7) * (t108 * t36 + t125 * t34 + t130 * t252 + t131 * t251 + t133 * t245 + t134 * t244) + m(6) * (t151 * t45 + t155 * t40 + t156 * t294 + t157 * t293 + t158 * t284 + t159 * t283) + m(5) * (t187 * t98 + t239 * t415 + t240 * t414 + t243 * t79 + t264 * t393 + t265 * t392); t599 + (t495 * t526 * t799 + t216 * t382) * t904 + (t125 * t36 + t133 * t252 + t134 * t251) * t901 + (t155 * t45 + t158 * t294 + t159 * t293) * t902 + (t243 * t98 + t264 * t415 + t265 * t414) * t903; t601 + (-t111 + (t592 * t636 + t595 * t635) * t584 - t828) * t577 + m(7) * (t145 * t148 + t146 * t147 + t270 * t38 + t271 * t37) + m(6) * (t153 * t191 + t154 * t190 + t280 * t67 + t281 * t66) + m(5) * (t149 * t300 + t150 * t299 + t181 * t292 + t182 * t291) + ((t61 / 0.2e1 + t71 / 0.2e1 + t96 / 0.2e1 + t123 / 0.2e1) * t595 + (t60 / 0.2e1 + t70 / 0.2e1 + t97 / 0.2e1 + t124 / 0.2e1) * t592 + (-t592 * t635 + t595 * t636) * qJD(1)) * t576 - t822; t597 + m(7) * (t108 * t28 + t127 * t34 + t130 * t147 + t131 * t148 + t244 * t37 + t245 * t38) + m(6) * (t151 * t39 + t156 * t190 + t157 * t191 + t164 * t40 + t283 * t66 + t284 * t67) + m(5) * (t101 * t187 + t149 * t392 + t150 * t393 + t239 * t291 + t240 * t292 + t275 * t79); t597 + m(7) * (t125 * t28 + t127 * t36 + t133 * t147 + t134 * t148 + t251 * t37 + t252 * t38) + m(6) * (t155 * t39 + t158 * t190 + t159 * t191 + t164 * t45 + t293 * t66 + t294 * t67) + m(5) * (t101 * t243 + t149 * t414 + t150 * t415 + t264 * t291 + t265 * t292 + t275 * t98); (t127 * t28 + t147 * t38 + t148 * t37) * t901 + (t164 * t39 + t190 * t67 + t191 * t66) * t902 + (t101 * t275 + t149 * t292 + t150 * t291) * t903 + (-t11 + t828 * t577 + (t644 * t592 + t595 * t909) * t584 + t822) * t577 + (t937 * t595 + t936 * t592 + t938 * t842 + ((-t61 - t71) * t595 + (-t60 - t70) * t592 + t928 * t584) * t577 + (t644 * t595 + (-t78 - t909) * t592) * qJD(1)) * t576 + t713; 0.2e1 * ((t270 * t595 + t271 * t592) * t899 + (t280 * t595 + t281 * t592) * t900) * t839 + 0.2e1 * ((t145 * t592 + t146 * t595 - t270 * t795 + t271 * t794) * t899 + (t153 * t592 + t154 * t595 - t280 * t795 + t281 * t794) * t900) * t576; 0.2e1 * ((t244 * t836 + t245 * t834 - t34) * t899 + (t283 * t836 + t284 * t834 - t40) * t900) * t577 + 0.2e1 * ((t108 * t584 + t130 * t595 + t131 * t592 + t244 * t794 - t245 * t795) * t899 + (t151 * t584 + t156 * t595 + t157 * t592 + t283 * t794 - t284 * t795) * t900) * t576; 0.2e1 * ((t251 * t836 + t252 * t834 - t36) * t899 + (t293 * t836 + t294 * t834 - t45) * t900) * t577 + 0.2e1 * ((t125 * t584 + t133 * t595 + t134 * t592 + t251 * t794 - t252 * t795) * t899 + (t155 * t584 + t158 * t595 + t159 * t592 + t293 * t794 - t294 * t795) * t900) * t576; 0.2e1 * ((t147 * t834 + t148 * t836 - t28) * t899 + (t190 * t834 + t191 * t836 - t39) * t900) * t577 + 0.2e1 * ((t127 * t584 - t147 * t795 + t148 * t794 + t37 * t592 + t38 * t595) * t899 + (t164 * t584 - t190 * t795 + t191 * t794 + t592 * t66 + t595 * t67) * t900) * t576; 0.4e1 * (t900 + t899) * (-0.1e1 + t799) * t576 * t839; t601 - t111 * t577 + m(7) * (t128 * t271 + t129 * t270 + t145 * t278 + t146 * t277); t606 + m(7) * (t108 * t80 + t128 * t244 + t129 * t245 + t130 * t277 + t131 * t278 + t250 * t34); t606 + m(7) * (t125 * t80 + t128 * t251 + t129 * t252 + t133 * t277 + t134 * t278 + t250 * t36); m(7) * (t127 * t80 + t128 * t148 + t129 * t147 + t250 * t28 + t277 * t38 + t278 * t37) + t603; m(7) * ((-t80 + (t277 * t595 + t278 * t592) * t584) * t577 + (t128 * t592 + t129 * t595 + t250 * t584 + (-t277 * t592 + t278 * t595) * qJD(1)) * t576); (t128 * t278 + t129 * t277 + t250 * t80) * t901 + t603;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
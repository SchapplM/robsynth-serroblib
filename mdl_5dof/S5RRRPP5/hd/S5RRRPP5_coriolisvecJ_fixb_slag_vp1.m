% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:36
% EndTime: 2019-12-31 20:58:22
% DurationCPUTime: 38.14s
% Computational Cost: add. (17008->862), mult. (23375->1050), div. (0->0), fcn. (18245->6), ass. (0->477)
t924 = -Icges(5,1) - Icges(6,1);
t476 = qJ(2) + qJ(3);
t453 = cos(t476);
t445 = Icges(4,4) * t453;
t452 = sin(t476);
t370 = Icges(4,1) * t452 + t445;
t764 = Icges(5,5) * t453;
t767 = Icges(6,4) * t453;
t917 = t924 * t452 + t764 + t767;
t923 = t370 - t917;
t478 = sin(qJ(1));
t719 = t453 * t478;
t722 = t452 * t478;
t480 = cos(qJ(1));
t757 = Icges(4,6) * t480;
t263 = Icges(4,4) * t719 - Icges(4,2) * t722 - t757;
t357 = Icges(5,3) * t452 + t764;
t255 = -Icges(5,6) * t480 + t357 * t478;
t361 = Icges(6,2) * t452 + t767;
t259 = Icges(6,6) * t480 + t361 * t478;
t913 = -t255 - t259;
t922 = t263 + t913;
t563 = -Icges(4,2) * t452 + t445;
t264 = Icges(4,6) * t478 + t480 * t563;
t718 = t453 * t480;
t420 = Icges(5,5) * t718;
t721 = t452 * t480;
t756 = Icges(5,6) * t478;
t256 = Icges(5,3) * t721 + t420 + t756;
t421 = Icges(6,4) * t718;
t754 = Icges(6,6) * t478;
t260 = Icges(6,2) * t721 + t421 - t754;
t902 = t256 + t260;
t921 = t264 - t902;
t443 = Icges(5,5) * t452;
t566 = Icges(5,1) * t453 + t443;
t267 = -Icges(5,4) * t480 + t478 * t566;
t444 = Icges(6,4) * t452;
t565 = Icges(6,1) * t453 + t444;
t265 = Icges(6,5) * t480 + t478 * t565;
t422 = Icges(4,4) * t722;
t765 = Icges(4,5) * t480;
t269 = Icges(4,1) * t719 - t422 - t765;
t905 = -t265 - t269;
t920 = t267 - t905;
t268 = Icges(5,4) * t478 + t480 * t566;
t266 = -Icges(6,5) * t478 + t480 * t565;
t769 = Icges(4,4) * t452;
t371 = Icges(4,1) * t453 - t769;
t270 = Icges(4,5) * t478 + t371 * t480;
t904 = t266 + t270;
t901 = t268 + t904;
t919 = t443 + t444 + (-Icges(6,2) - Icges(5,3)) * t453;
t364 = Icges(4,2) * t453 + t769;
t915 = t364 - t919;
t918 = -t357 - t361;
t916 = -t370 * t480 - t921;
t912 = t919 * t478 + t920;
t911 = (Icges(4,6) - Icges(5,6) + Icges(6,6)) * t453 + (Icges(5,4) + Icges(4,5) - Icges(6,5)) * t452;
t910 = t371 + t565 + t566 - t915;
t909 = -t563 - t923 - t918;
t908 = t915 * t480 - t901;
t907 = t923 * t478 + t922;
t359 = Icges(4,5) * t453 - Icges(4,6) * t452;
t258 = Icges(4,3) * t478 + t359 * t480;
t363 = Icges(5,4) * t453 + Icges(5,6) * t452;
t262 = Icges(5,2) * t478 + t363 * t480;
t906 = -t258 - t262;
t903 = t915 * t452 - t453 * t923;
t355 = Icges(6,5) * t453 + Icges(6,6) * t452;
t900 = t355 - t359 - t363;
t555 = t263 * t452 - t269 * t453;
t557 = t259 * t452 + t265 * t453;
t559 = t255 * t452 + t267 * t453;
t899 = -t555 + t557 + t559;
t473 = qJD(2) + qJD(3);
t898 = t910 * t473;
t897 = t909 * t473;
t896 = -t901 * qJD(1) + t907 * t473;
t398 = t473 * t480;
t895 = t917 * t398 + t916 * t473 + (-t371 * t478 - t265 - t267 + t765) * qJD(1);
t397 = t473 * t478;
t894 = -t364 * t397 + t912 * t473 + (t918 * t480 + t264 + t754 - t756) * qJD(1);
t893 = -t908 * t473 + (-t478 * t563 + t757 - t913) * qJD(1);
t892 = t911 * t480;
t891 = t911 * t478;
t832 = rSges(6,1) + pkin(4);
t708 = -t256 * t722 - t268 * t719;
t198 = t270 * t719;
t594 = t258 * t480 - t198;
t254 = -Icges(6,3) * t478 + t355 * t480;
t96 = t480 * t254 + t260 * t722 + t266 * t719;
t870 = -t264 * t722 - t594 + t96;
t840 = -t262 * t480 - t708 + t870;
t890 = t478 * t258 + t260 * t721 + t904 * t718;
t751 = Icges(4,3) * t480;
t257 = Icges(4,5) * t719 - Icges(4,6) * t722 - t751;
t889 = -t478 * t257 - t259 * t721 + t905 * t718;
t253 = Icges(6,3) * t480 + t355 * t478;
t261 = -Icges(5,2) * t480 + t363 * t478;
t888 = (-t253 + t261) * qJD(1);
t887 = t903 * t478 + t892;
t886 = -t903 * t480 + t891;
t885 = (t254 + t906) * qJD(1);
t532 = t555 * t478;
t738 = t257 * t480;
t739 = t253 * t480;
t871 = t478 * t557 - t532 - t738 + t739;
t740 = t253 * t478;
t884 = t263 * t721 + t740 + t889;
t869 = -t254 * t478 - t264 * t721 + t890;
t883 = t911 * qJD(1) + t897 * t452 + t898 * t453;
t882 = -t893 * t452 + t895 * t453 - t885;
t881 = -qJD(1) * t257 + t894 * t452 + t896 * t453 - t888;
t880 = t907 * t398 + (t924 * t721 + t420 + t421 + t916) * t397 + t910 * qJD(1);
t879 = (-Icges(4,2) * t719 - t422 + t912) * t398 + t908 * t397 + t909 * qJD(1);
t878 = -t903 * qJD(1) + t900 * t473;
t877 = t899 * qJD(1) + t891 * t473 + t885;
t736 = t264 * t452;
t876 = t888 + t892 * t473 + (t359 * t478 + t902 * t452 + t901 * t453 - t736 - t751) * qJD(1);
t874 = rSges(6,2) * t721 + t718 * t832;
t230 = t478 * t261;
t103 = t255 * t721 + t267 * t718 + t230;
t104 = t256 * t721 + t478 * t262 + t268 * t718;
t873 = qJD(1) * t886 - t103 * t398 + t104 * t397;
t97 = -t261 * t480 + t478 * t559;
t872 = qJD(1) * t887 + t398 * t97;
t783 = rSges(6,2) * t452;
t377 = rSges(6,1) * t453 + t783;
t272 = rSges(6,3) * t480 + t377 * t478;
t640 = t398 * t453;
t389 = rSges(6,2) * t640;
t641 = t452 * t398;
t659 = qJD(1) * t478;
t523 = -t453 * t659 - t641;
t653 = qJD(5) * t478;
t658 = qJD(1) * t480;
t712 = -rSges(6,1) * t641 + pkin(4) * t523 - qJ(5) * t658 - qJD(1) * t272 + t389 - t653;
t373 = rSges(6,1) * t452 - rSges(6,2) * t453;
t328 = t373 * t478;
t626 = t453 * t658;
t642 = t452 * t397;
t652 = qJD(5) * t480;
t782 = rSges(6,3) * t478;
t711 = -t473 * t328 + (t377 * t480 - t782) * qJD(1) - qJ(5) * t659 + t652 + (t626 - t642) * pkin(4);
t688 = pkin(4) * t719 + qJ(5) * t480 + t272;
t859 = -qJ(5) * t478 - t782 + t874;
t852 = -t478 * (t254 + t557) - t739 + t890;
t851 = t877 * t478 + t881 * t480;
t850 = -t876 * t478 + t882 * t480;
t849 = t881 * t478 - t877 * t480;
t848 = t882 * t478 + t876 * t480;
t847 = t397 * t840 - t871 * t398 - t872;
t846 = t869 * t397 + t884 * t398 + t873;
t845 = -t878 * t478 + t883 * t480;
t844 = t883 * t478 + t878 * t480;
t843 = t896 * t452 - t894 * t453;
t842 = t895 * t452 + t893 * t453;
t841 = t97 + t871;
t839 = t103 - t884;
t838 = t104 + t869;
t837 = t920 * t452 + t922 * t453;
t836 = t901 * t452 + t921 * t453;
t835 = t879 * t452 + t880 * t453;
t834 = -t900 * qJD(1) - t892 * t397 + t891 * t398;
t833 = 2 * qJD(2);
t628 = t452 * t659;
t674 = rSges(5,2) * t658 + rSges(5,3) * t640;
t159 = rSges(5,1) * t523 - rSges(5,3) * t628 + t674;
t580 = qJD(1) * t473;
t381 = t480 * t580;
t477 = sin(qJ(2));
t647 = qJD(1) * qJD(2);
t617 = t480 * t647;
t479 = cos(qJ(2));
t713 = t479 * (qJD(2) ^ 2);
t639 = t478 * t713;
t481 = -pkin(7) - pkin(6);
t450 = t480 * t481;
t451 = pkin(6) * t658;
t656 = qJD(2) * t480;
t621 = t477 * t656;
t581 = pkin(2) * t621;
t446 = pkin(2) * t479 + pkin(1);
t787 = pkin(1) - t446;
t185 = -t581 - t451 + (t478 * t787 - t450) * qJD(1);
t353 = qJD(1) * (-pkin(1) * t659 + t451);
t710 = qJD(1) * t185 + t353;
t509 = (-t477 * t617 - t639) * pkin(2) + t710;
t372 = pkin(3) * t452 - qJ(4) * t453;
t374 = rSges(5,1) * t452 - rSges(5,3) * t453;
t677 = -t372 - t374;
t723 = t452 * t473;
t226 = qJ(4) * t723 + (pkin(3) * t473 - qJD(4)) * t453;
t378 = rSges(5,1) * t453 + rSges(5,3) * t452;
t297 = t378 * t473;
t700 = -t226 - t297;
t654 = qJD(4) * t480;
t413 = t452 * t654;
t675 = qJ(4) * t640 + t413;
t132 = pkin(3) * t523 - qJ(4) * t628 + t675;
t655 = qJD(4) * t478;
t412 = t453 * t655;
t822 = t473 * t412 + (t132 + t413) * qJD(1);
t39 = qJD(1) * t159 + t381 * t677 + t397 * t700 + t509 + t822;
t831 = t39 * t478;
t329 = t374 * t478;
t466 = t478 * rSges(5,2);
t162 = -t473 * t329 + (t378 * t480 + t466) * qJD(1);
t380 = t478 * t580;
t414 = t453 * t654;
t780 = pkin(2) * qJD(2);
t645 = t477 * t780;
t434 = t478 * t645;
t547 = -pkin(2) * t480 * t713 + qJD(1) * t434;
t517 = -t398 * t226 + t380 * t372 + t473 * t414 + t547;
t627 = t452 * t658;
t243 = t397 * t453 + t627;
t392 = pkin(3) * t642;
t619 = t452 * t655;
t133 = pkin(3) * t626 + qJ(4) * t243 - t392 + t619;
t470 = t478 * pkin(6);
t666 = t481 * t659 + t434;
t186 = (-t480 * t787 - t470) * qJD(1) - t666;
t416 = t480 * pkin(1) + t470;
t393 = t416 * qJD(1);
t709 = -t186 - t393;
t531 = -t133 - t619 + t709;
t40 = -t297 * t398 + t374 * t380 + (-t162 + t531) * qJD(1) + t517;
t541 = -t397 * t372 - t434 + t619;
t424 = t480 * t446;
t584 = -t478 * t481 + t424;
t252 = t584 - t416;
t276 = rSges(5,1) * t718 + rSges(5,3) * t721 + t466;
t336 = pkin(3) * t718 + qJ(4) * t721;
t687 = t276 + t336;
t633 = -t252 - t687;
t87 = -t374 * t397 + (t416 - t633) * qJD(1) + t541;
t830 = qJD(1) * t87 + t40;
t471 = t480 * pkin(6);
t415 = pkin(1) * t478 - t471;
t669 = -t478 * t446 - t450;
t251 = t415 + t669;
t394 = qJD(1) * t415;
t826 = qJD(1) * t251 - t394;
t463 = Icges(3,4) * t479;
t564 = -Icges(3,2) * t477 + t463;
t405 = Icges(3,1) * t477 + t463;
t399 = pkin(4) * t628;
t821 = pkin(4) * t640 + t398 * t377 - t399;
t375 = rSges(4,1) * t452 + rSges(4,2) * t453;
t330 = t375 * t478;
t335 = t375 * t480;
t785 = rSges(4,1) * t453;
t379 = -rSges(4,2) * t452 + t785;
t274 = rSges(4,1) * t719 - rSges(4,2) * t722 - t480 * rSges(4,3);
t464 = t478 * rSges(4,3);
t277 = rSges(4,1) * t718 - rSges(4,2) * t721 + t464;
t657 = qJD(2) * t478;
t702 = -t251 * t657 + t252 * t656;
t92 = t274 * t397 + t277 * t398 + t702;
t525 = -t375 * t398 - t581;
t696 = t251 - t415;
t93 = (-t274 + t696) * qJD(1) + t525;
t695 = -t252 - t277;
t94 = -t375 * t397 - t434 + (t416 - t695) * qJD(1);
t820 = -(qJD(1) * t330 - t398 * t379) * t93 - t92 * (-t397 * t330 - t335 * t398) - t94 * (-qJD(1) * t335 - t379 * t397);
t715 = t478 * t479;
t717 = t477 * t478;
t752 = Icges(3,3) * t480;
t299 = Icges(3,5) * t715 - Icges(3,6) * t717 - t752;
t439 = Icges(3,4) * t717;
t766 = Icges(3,5) * t480;
t303 = Icges(3,1) * t715 - t439 - t766;
t758 = Icges(3,6) * t480;
t301 = Icges(3,4) * t715 - Icges(3,2) * t717 - t758;
t734 = t301 * t477;
t553 = -t303 * t479 + t734;
t109 = -t299 * t480 - t478 * t553;
t402 = Icges(3,5) * t479 - Icges(3,6) * t477;
t401 = Icges(3,5) * t477 + Icges(3,6) * t479;
t528 = qJD(2) * t401;
t770 = Icges(3,4) * t477;
t406 = Icges(3,1) * t479 - t770;
t304 = Icges(3,5) * t478 + t406 * t480;
t302 = Icges(3,6) * t478 + t480 * t564;
t733 = t302 * t477;
t552 = -t304 * t479 + t733;
t812 = -t480 * t528 + (-t402 * t478 + t552 + t752) * qJD(1);
t300 = Icges(3,3) * t478 + t402 * t480;
t661 = qJD(1) * t300;
t811 = qJD(1) * t553 - t478 * t528 + t661;
t403 = Icges(3,2) * t479 + t770;
t548 = t403 * t477 - t405 * t479;
t807 = t548 * qJD(1) + t402 * qJD(2);
t806 = t478 * (-t403 * t480 + t304) - t480 * (-Icges(3,2) * t715 + t303 - t439);
t646 = -pkin(3) - t832;
t748 = qJ(4) * t452;
t805 = -qJD(5) + (t646 * t453 - t446 - t748 - t783) * qJD(1);
t801 = m(5) / 0.2e1;
t800 = m(6) / 0.2e1;
t799 = t380 / 0.2e1;
t798 = t381 / 0.2e1;
t797 = -t397 / 0.2e1;
t796 = t397 / 0.2e1;
t795 = -t398 / 0.2e1;
t794 = t398 / 0.2e1;
t793 = t478 / 0.2e1;
t792 = -t480 / 0.2e1;
t791 = -rSges(5,1) - pkin(3);
t790 = pkin(2) * t477;
t789 = -qJD(1) / 0.2e1;
t788 = qJD(1) / 0.2e1;
t786 = rSges(3,1) * t479;
t784 = rSges(3,2) * t479;
t465 = t478 * rSges(3,3);
t778 = t93 * t375;
t777 = -rSges(6,2) - qJ(4);
t776 = -rSges(5,3) - qJ(4);
t775 = rSges(6,3) + qJ(5);
t376 = pkin(3) * t453 + t748;
t331 = t376 * t478;
t539 = -qJD(4) * t453 + t397 * t331 + t702;
t632 = t336 + t859;
t60 = t397 * t688 + t398 * t632 + t539;
t747 = qJD(1) * t60;
t469 = t480 * rSges(5,2);
t273 = t378 * t478 - t469;
t75 = t273 * t397 + t398 * t687 + t539;
t746 = qJD(1) * t75;
t665 = rSges(3,2) * t717 + t480 * rSges(3,3);
t305 = rSges(3,1) * t715 - t665;
t408 = rSges(3,1) * t477 + t784;
t622 = t408 * t656;
t166 = -t622 + (-t305 - t415) * qJD(1);
t743 = t166 * t478;
t742 = t166 * t480;
t623 = t408 * t657;
t714 = t479 * t480;
t716 = t477 * t480;
t306 = rSges(3,1) * t714 - rSges(3,2) * t716 + t465;
t682 = t306 + t416;
t167 = qJD(1) * t682 - t623;
t350 = t408 * t480;
t741 = t167 * t350;
t725 = t401 * t478;
t724 = t401 * t480;
t720 = t453 * t473;
t327 = t372 * t478;
t442 = qJD(4) * t452;
t704 = -t397 * t327 + t442;
t701 = -t478 * t251 + t480 * t252;
t699 = t478 * t274 + t480 * t277;
t698 = -t478 * t299 - t303 * t714;
t697 = t478 * t300 + t304 * t714;
t686 = t478 * t331 + t480 * t336;
t332 = t372 * t480;
t685 = -qJD(1) * t332 + t412;
t337 = t372 * t659;
t681 = t374 * t659 + t337;
t676 = -t376 - t378;
t673 = rSges(4,2) * t628 + rSges(4,3) * t658;
t672 = -t403 + t406;
t671 = t405 + t564;
t670 = -t398 * t376 + t414;
t625 = t477 * t659;
t667 = rSges(3,2) * t625 + rSges(3,3) * t658;
t660 = qJD(1) * t402;
t170 = -t478 * t548 - t724;
t648 = t170 * qJD(1);
t644 = t479 * t780;
t643 = -t481 - t775;
t638 = t480 * t132 + t478 * t133 + t331 * t658;
t160 = rSges(4,1) * t523 - rSges(4,2) * t640 + t673;
t163 = -t473 * t330 + (t379 * t480 + t464) * qJD(1);
t637 = t480 * t160 + t478 * t163 + t274 * t658;
t636 = t480 * t185 + t478 * t186 - t251 * t658;
t634 = -t331 + t696;
t631 = t373 * t659 + t337 + t399;
t630 = t392 + t666;
t629 = t424 + t336;
t624 = t477 * t658;
t618 = -pkin(1) - t786;
t616 = t659 / 0.2e1;
t615 = t658 / 0.2e1;
t614 = -t657 / 0.2e1;
t611 = t656 / 0.2e1;
t610 = -t375 - t790;
t609 = -pkin(4) * t452 - t373;
t607 = (-t478 ^ 2 - t480 ^ 2) * t477;
t240 = t304 * t715;
t593 = t300 * t480 - t240;
t592 = -t257 + t736;
t585 = -t299 + t733;
t579 = -t252 - t632;
t578 = t478 * t273 + t480 * t276 + t686;
t577 = t677 - t790;
t576 = -t372 + t609;
t298 = t379 * t473;
t573 = -t298 - t644;
t570 = -rSges(3,2) * t477 + t786;
t569 = -t478 * t94 - t480 * t93;
t560 = -t167 * t478 - t742;
t168 = t301 * t479 + t303 * t477;
t169 = t302 * t479 + t304 * t477;
t296 = t377 * t473;
t546 = -pkin(4) * t720 - t226 - t296;
t545 = t413 - t581;
t544 = -t644 + t700;
t543 = -t96 + t740;
t77 = t652 + t609 * t397 + (t416 - t579) * qJD(1) + t541;
t542 = t77 * t576;
t435 = pkin(2) * t625;
t540 = qJD(1) * t327 + t435 + t670;
t538 = t480 * t159 + t478 * t162 + t273 * t658 + t638;
t537 = t688 * t478 + t859 * t480 + t686;
t349 = t408 * t478;
t533 = t576 - t790;
t530 = qJD(2) * t405;
t529 = qJD(2) * t403;
t110 = -t302 * t717 - t593;
t527 = (-t109 * t480 + t110 * t478) * qJD(2);
t111 = -t301 * t716 - t698;
t112 = -t302 * t716 + t697;
t526 = (-t111 * t480 + t112 * t478) * qJD(2);
t164 = (t305 * t478 + t306 * t480) * qJD(2);
t524 = -t581 - t653;
t519 = -t252 * t478 * t647 + t185 * t656 + t186 * t657 - t251 * t617;
t518 = t301 * t480 - t302 * t478;
t516 = -qJD(1) * t331 + t545 + t826;
t515 = t546 - t644;
t514 = t711 * t478 + t712 * t480 + t688 * t658 + t638;
t513 = (-t477 * t671 + t479 * t672) * qJD(1);
t512 = -t446 + t676;
t511 = t397 * t133 + t381 * t331 + t473 * t442 + t519;
t179 = qJD(1) * t302 - t478 * t529;
t181 = qJD(1) * t304 - t478 * t530;
t493 = qJD(1) * t299 - qJD(2) * t168 - t179 * t477 + t181 * t479;
t178 = -t480 * t529 + (-t478 * t564 + t758) * qJD(1);
t180 = -t480 * t530 + (-t406 * t478 + t766) * qJD(1);
t492 = -qJD(2) * t169 - t178 * t477 + t180 * t479 + t661;
t385 = t564 * qJD(2);
t386 = t406 * qJD(2);
t491 = qJD(1) * t401 - t385 * t477 + t386 * t479 + (-t403 * t479 - t405 * t477) * qJD(2);
t334 = t374 * t480;
t490 = t75 * (-t397 * t329 + (-t332 - t334) * t398 + t704) + t87 * (-qJD(1) * t334 + t397 * t676 + t685);
t489 = -t477 * t806 + t518 * t479;
t485 = (t840 * t478 - t841 * t480) * t799 + (t838 * t478 - t839 * t480) * t798 + (t834 * t478 + t835 * t480) * t797 + (t851 * t480 + t850 * t478 + (t839 * t478 + t838 * t480) * qJD(1)) * t796 + (t849 * t480 + t848 * t478 + (t841 * t478 + t840 * t480) * qJD(1)) * t795 + (t835 * t478 - t834 * t480) * t794 + (t845 * qJD(1) + t839 * t380 + t838 * t381 + t850 * t397 + t851 * t398) * t793 + (t844 * qJD(1) + t841 * t380 + t840 * t381 + t848 * t397 + t849 * t398) * t792 + (t880 * t452 - t879 * t453) * t789 + (t843 * t480 + t842 * t478 + (t837 * t478 + t836 * t480) * qJD(1)) * t788 + t847 * t616 + t846 * t615;
t333 = t373 * t480;
t484 = t60 * (-t397 * t328 + t704 + (-pkin(4) * t721 - t332 - t333) * t398) + t77 * (-pkin(4) * t627 - qJD(1) * t333 + t685) + (t77 * (-pkin(4) * t453 - t376 - t377) - t60 * pkin(4) * t722) * t397;
t390 = t570 * qJD(2);
t249 = t398 * t378;
t246 = t398 * t372;
t244 = -t628 + t640;
t190 = (t397 * t478 + t398 * t480) * t452;
t183 = -qJD(2) * t349 + (t480 * t570 + t465) * qJD(1);
t182 = -t656 * t784 + (-t479 * t659 - t621) * rSges(3,1) + t667;
t171 = -t480 * t548 + t725;
t165 = t171 * qJD(1);
t108 = -t390 * t656 + (-t183 - t393 + t623) * qJD(1);
t107 = -t390 * t657 + t353 + (t182 - t622) * qJD(1);
t91 = t491 * t478 - t480 * t807;
t90 = t478 * t807 + t491 * t480;
t89 = -qJD(2) * t552 + t178 * t479 + t180 * t477;
t88 = -qJD(2) * t553 + t179 * t479 + t181 * t477;
t86 = -t374 * t398 - t246 + (-t273 + t634) * qJD(1) + t545;
t79 = -t298 * t398 + t375 * t380 + (-t163 + t709) * qJD(1) + t547;
t78 = qJD(1) * t160 - t298 * t397 - t375 * t381 + t509;
t76 = -t246 + t413 + t609 * t398 + (t634 - t688) * qJD(1) + t524;
t62 = t165 + t526;
t61 = t527 + t648;
t26 = t160 * t398 + t163 * t397 + t274 * t381 - t277 * t380 + t519;
t25 = -t296 * t398 + t373 * t380 + (t380 * t452 - t398 * t720) * pkin(4) + (t531 - t652 - t711) * qJD(1) + t517;
t24 = -pkin(2) * t639 + t546 * t397 + t576 * t381 + (t524 + t712) * qJD(1) + t710 + t822;
t14 = t162 * t397 + t273 * t381 + (t132 + t159) * t398 - t687 * t380 + t511;
t13 = t711 * t397 + t688 * t381 + (t132 + t712) * t398 - t632 * t380 + t511;
t1 = [(t165 + ((t110 - t240 + (t300 + t734) * t480 + t698) * t480 + t697 * t478) * qJD(2)) * t611 + (t61 - t648 + ((t480 * t585 + t112 - t697) * t480 + (t478 * t585 + t111 + t593) * t478) * qJD(2)) * t614 + (t89 + t90) * t657 / 0.2e1 + (-t548 * qJD(2) + t385 * t479 + t386 * t477 + t898 * t452 - t897 * t453) * qJD(1) + (t25 * t669 + t76 * t630 + t24 * (t629 + t874) + (-t25 * t775 + t805 * t76) * t480 + (t24 * t643 + (t473 * t76 * t777 + t25 * t646) * t453 + (t25 * t777 + t76 * (t473 * t832 - qJD(4))) * t452 + t76 * t775 * qJD(1)) * t478 + (t389 + t675 + (qJD(1) * t643 + t646 * t723 - t645) * t480 + t805 * t478 + qJD(1) * t688 - t516 + t653 + t76) * t77 - t542 * t398) * m(6) + (-(-qJD(1) * t273 + t398 * t677 + t516 - t86) * t87 + t40 * (t469 + t669) + t86 * t630 + t39 * (t276 + t629) + t87 * (t674 + t675) + (t87 * (t723 * t791 - t645) + (-t87 * t481 + t512 * t86) * qJD(1)) * t480 + (-t39 * t481 + (t473 * t776 * t86 + t40 * t791) * t453 + (t40 * t776 + t86 * (rSges(5,1) * t473 - qJD(4))) * t452 + (-t86 * rSges(5,2) + t512 * t87) * qJD(1)) * t478) * m(5) + (t79 * (-t274 + t669) + t93 * t666 + t78 * (t277 + t584) + t94 * (-t581 + t673) + (-t335 * t94 + t478 * t778) * t473 + ((-t93 * rSges(4,3) + t94 * (-t446 - t785)) * t478 + (t93 * (-t379 - t446) - t94 * t481) * t480) * qJD(1) - (-qJD(1) * t274 + t525 + t826 - t93) * t94) * m(4) + (t108 * (t478 * t618 + t471 + t665) + t107 * t682 + t167 * (t451 + t667) + (t408 * t743 - t741) * qJD(2) + ((-pkin(1) - t570) * t742 + (t166 * (-rSges(3,3) - pkin(6)) + t167 * t618) * t478) * qJD(1) - (-qJD(1) * t305 - t166 - t394 - t622) * t167) * m(3) + ((t543 + (t263 * t480 + t264 * t478) * t452 + t594 + t870 + t889) * t398 + (-t269 * t719 + t738 + (t263 * t478 - t264 * t480) * t452 + t852 + t871) * t397 + t873) * t794 - (t88 + t91 + t62) * t656 / 0.2e1 + (t837 - t887) * t799 + (t836 + t886) * t798 + ((t592 * t480 - t532 - t852 + t869) * t398 + (t592 * t478 - t198 - t230 + t543 + t708 + (-t899 - t906) * t480 + t839) * t397 + t847 + t872) * t797 + (t842 + t845) * t796 + (-t843 + t844 + t846) * t795 + ((t168 + t170) * t478 + (t169 + t171) * t480) * t647 / 0.2e1; ((-t657 * t724 + t660) * t478 + (t513 + (t478 * t725 + t489) * qJD(2)) * t480) * t614 + ((-t656 * t725 - t660) * t480 + (t513 + (t480 * t724 + t489) * qJD(2)) * t478) * t611 + ((t477 * t672 + t479 * t671) * qJD(1) + (t518 * t477 + t479 * t806) * qJD(2)) * t789 + (t478 * t89 - t480 * t88 + (t168 * t478 + t169 * t480) * qJD(1)) * t788 + t485 + (qJD(1) * t90 + (t478 * (t478 * t812 + t492 * t480) - t480 * (t478 * t811 + t493 * t480) + (t111 * t478 + t112 * t480) * qJD(1)) * t833) * t793 + (qJD(1) * t91 + (t478 * (t492 * t478 - t480 * t812) - t480 * (t493 * t478 - t480 * t811) + (t109 * t478 + t110 * t480) * qJD(1)) * t833) * t792 + (t61 + t527) * t616 + (t62 + t526) * t615 + (t76 * (t435 + t631) + t13 * (t537 + t701) + t60 * (t514 + t636) + (t515 * t76 + (qJD(1) * t77 + t25) * t533) * t480 + (t24 * t533 + t515 * t77 + t579 * t747) * t478 - t76 * (qJD(1) * t328 + t540 - t821) - (-t77 * t624 + ((-t478 * t77 - t480 * t76) * t479 + t60 * t607) * qJD(2)) * pkin(2) - t484) * m(6) + (t86 * (t435 + t681) + t14 * (t578 + t701) + t75 * (t538 + t636) + (t544 * t86 + t577 * t830) * t480 + (t39 * t577 + t544 * t87 + t633 * t746) * t478 - t86 * (qJD(1) * t329 - t249 + t540) - (-t87 * t624 + ((-t478 * t87 - t480 * t86) * t479 + t75 * t607) * qJD(2)) * pkin(2) - t490) * m(5) + (t26 * (t699 + t701) + t92 * (t636 + t637) + (t573 * t93 + (qJD(1) * t94 + t79) * t610) * t480 + (t78 * t610 + t94 * t573 + (t695 * t92 + t778) * qJD(1)) * t478 - (-t94 * t624 + (t479 * t569 + t607 * t92) * qJD(2)) * pkin(2) + t820) * m(4) + (0.2e1 * t164 * (t182 * t480 + t183 * t478 + (t305 * t480 - t306 * t478) * qJD(1)) + t560 * t390 + (-t107 * t478 - t108 * t480 + (-t167 * t480 + t743) * qJD(1)) * t408 - (t166 * t349 - t741) * qJD(1) - (t164 * (-t349 * t478 - t350 * t480) + t560 * t570) * qJD(2)) * m(3); t485 + (t13 * t537 + t60 * t514 + (qJD(1) * t542 + t25 * t576) * t480 + (t24 * t576 + t546 * t77 - t632 * t747) * t478 - t484 + (t631 + t546 * t480 - (t327 + t328) * qJD(1) - t670 + t821) * t76) * m(6) + (-t490 + t14 * t578 + t75 * t538 + (-t687 * t746 + t700 * t87) * t478 + (t480 * t830 + t831) * t677 + (t249 - (t327 + t329) * qJD(1) - t670 + t681 + t700 * t480) * t86) * m(5) + (t26 * t699 + t92 * (-t277 * t659 + t637) + t569 * t298 + (-t78 * t478 - t79 * t480 + (t478 * t93 - t480 * t94) * qJD(1)) * t375 + t820) * m(4); -m(5) * (t190 * t75 + t243 * t87 + t244 * t86) - m(6) * (t190 * t60 + t243 * t77 + t244 * t76) + 0.2e1 * ((t397 * t87 + t398 * t86 - t14) * t801 + (t397 * t77 + t398 * t76 - t13) * t800) * t453 + 0.2e1 * ((t40 * t480 + t473 * t75 + t658 * t87 - t659 * t86 + t831) * t801 + (t24 * t478 + t25 * t480 + t473 * t60 + t658 * t77 - t659 * t76) * t800) * t452; 0.2e1 * (t24 * t480 / 0.2e1 - t25 * t478 / 0.2e1 - t60 * (t397 * t480 - t398 * t478) / 0.2e1) * m(6);];
tauc = t1(:);
% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:48
% EndTime: 2019-12-05 18:38:43
% DurationCPUTime: 43.45s
% Computational Cost: add. (26976->1004), mult. (25729->1255), div. (0->0), fcn. (20188->10), ass. (0->557)
t943 = Icges(4,3) + Icges(5,3);
t512 = qJ(2) + qJ(3);
t475 = pkin(9) + t512;
t457 = sin(t475);
t458 = cos(t475);
t480 = sin(t512);
t481 = cos(t512);
t936 = Icges(4,5) * t481 + Icges(5,5) * t458 - Icges(4,6) * t480 - Icges(5,6) * t457;
t942 = Icges(4,5) * t480 + Icges(5,5) * t457 + Icges(4,6) * t481 + Icges(5,6) * t458;
t516 = cos(qJ(1));
t941 = t943 * t516;
t514 = sin(qJ(1));
t763 = t481 * t514;
t766 = t480 * t514;
t768 = t458 * t514;
t771 = t457 * t514;
t937 = -Icges(4,5) * t763 - Icges(5,5) * t768 + Icges(4,6) * t766 + Icges(5,6) * t771 + t941;
t940 = t943 * t514 + t936 * t516;
t818 = Icges(5,4) * t457;
t360 = Icges(5,2) * t458 + t818;
t445 = Icges(5,4) * t458;
t362 = Icges(5,1) * t457 + t445;
t819 = Icges(4,4) * t480;
t386 = Icges(4,2) * t481 + t819;
t456 = Icges(4,4) * t481;
t388 = Icges(4,1) * t480 + t456;
t939 = t360 * t457 - t362 * t458 + t386 * t480 - t388 * t481;
t806 = Icges(5,6) * t516;
t266 = Icges(5,4) * t768 - Icges(5,2) * t771 - t806;
t413 = Icges(5,4) * t771;
t814 = Icges(5,5) * t516;
t268 = Icges(5,1) * t768 - t413 - t814;
t807 = Icges(4,6) * t516;
t298 = Icges(4,4) * t763 - Icges(4,2) * t766 - t807;
t436 = Icges(4,4) * t766;
t815 = Icges(4,5) * t516;
t300 = Icges(4,1) * t763 - t436 - t815;
t935 = t266 * t457 - t268 * t458 + t298 * t480 - t300 * t481;
t363 = Icges(5,1) * t458 - t818;
t269 = Icges(5,5) * t514 + t363 * t516;
t389 = Icges(4,1) * t481 - t819;
t301 = Icges(4,5) * t514 + t389 * t516;
t938 = -t269 * t768 - t301 * t763;
t934 = t942 * t516;
t933 = t942 * t514;
t594 = -Icges(5,2) * t457 + t445;
t267 = Icges(5,6) * t514 + t516 * t594;
t595 = -Icges(4,2) * t480 + t456;
t299 = Icges(4,6) * t514 + t516 * t595;
t932 = t267 * t457 + t299 * t480;
t931 = t939 * t514 + t934;
t930 = -t939 * t516 + t933;
t929 = t935 * t514;
t928 = t940 * t516 + t938;
t927 = t940 * qJD(1);
t762 = t481 * t516;
t767 = t458 * t516;
t874 = -t269 * t767 - t301 * t762 - t940 * t514;
t926 = -t268 * t767 - t300 * t762 + t937 * t514;
t925 = t937 * t516;
t906 = t925 - t929;
t905 = -t267 * t771 - t299 * t766 - t928;
t765 = t480 * t516;
t770 = t457 * t516;
t904 = -t266 * t770 - t298 * t765 - t926;
t903 = -t267 * t770 - t299 * t765 - t874;
t509 = qJD(2) + qJD(3);
t888 = -t386 + t389;
t617 = t888 * t509;
t887 = t388 + t595;
t618 = t887 * t509;
t890 = -t360 + t363;
t621 = t890 * t509;
t889 = t362 + t594;
t622 = t889 * t509;
t924 = t942 * qJD(1) - t457 * t622 + t458 * t621 - t480 * t618 + t481 * t617;
t893 = -t388 * t516 - t299;
t635 = (-t389 * t514 + t815) * qJD(1) + t893 * t509;
t892 = -t386 * t516 + t301;
t637 = (-t514 * t595 + t807) * qJD(1) + t892 * t509;
t896 = -t362 * t516 - t267;
t639 = (-t363 * t514 + t814) * qJD(1) + t896 * t509;
t895 = -t360 * t516 + t269;
t641 = (-t514 * t594 + t806) * qJD(1) + t895 * t509;
t923 = -t457 * t641 + t458 * t639 - t480 * t637 + t481 * t635 + t927;
t894 = -t388 * t514 - t298;
t634 = qJD(1) * t301 + t509 * t894;
t760 = t509 * t514;
t636 = qJD(1) * t299 + t300 * t509 - t386 * t760;
t897 = -t362 * t514 - t266;
t638 = qJD(1) * t269 + t509 * t897;
t640 = qJD(1) * t267 + t268 * t509 - t360 * t760;
t922 = qJD(1) * t937 + t457 * t640 - t458 * t638 + t480 * t636 - t481 * t634;
t921 = t939 * qJD(1) + t936 * t509;
t920 = qJD(1) * t935 - t509 * t933 + t927;
t919 = -t934 * t509 + (-t269 * t458 - t301 * t481 - t514 * t936 + t932 + t941) * qJD(1);
t918 = t930 * qJD(1);
t917 = t931 * qJD(1);
t916 = -t920 * t514 + t922 * t516;
t915 = t919 * t514 + t923 * t516;
t914 = t922 * t514 + t920 * t516;
t913 = t923 * t514 - t919 * t516;
t479 = qJD(2) * t514;
t417 = qJD(3) * t514 + t479;
t418 = t509 * t516;
t912 = t905 * t417 - t906 * t418 - t917;
t911 = t903 * t417 - t904 * t418 + t918;
t910 = t921 * t514 + t924 * t516;
t909 = t924 * t514 - t921 * t516;
t908 = -t457 * t638 - t458 * t640 - t480 * t634 - t481 * t636;
t907 = t457 * t639 + t458 * t641 + t480 * t635 + t481 * t637;
t902 = t266 * t458 + t268 * t457 + t298 * t481 + t300 * t480;
t901 = t267 * t458 + t269 * t457 + t299 * t481 + t301 * t480;
t466 = qJ(5) + t475;
t448 = cos(t466);
t773 = t448 * t514;
t447 = sin(t466);
t775 = t447 * t514;
t805 = Icges(6,6) * t516;
t239 = Icges(6,4) * t773 - Icges(6,2) * t775 - t805;
t431 = Icges(6,4) * t448;
t337 = Icges(6,1) * t447 + t431;
t900 = -t337 * t514 - t239;
t593 = -Icges(6,2) * t447 + t431;
t240 = Icges(6,6) * t514 + t516 * t593;
t899 = -t337 * t516 - t240;
t817 = Icges(6,4) * t447;
t338 = Icges(6,1) * t448 - t817;
t242 = Icges(6,5) * t514 + t338 * t516;
t335 = Icges(6,2) * t448 + t817;
t898 = -t335 * t516 + t242;
t891 = t337 + t593;
t539 = qJD(1) * t888 + t417 * t893 - t418 * t894;
t540 = qJD(1) * t890 + t417 * t896 - t418 * t897;
t855 = qJD(1) * t889 + t417 * t895 - t418 * (-Icges(5,2) * t768 + t268 - t413);
t856 = qJD(1) * t887 + t417 * t892 - t418 * (-Icges(4,2) * t763 + t300 - t436);
t886 = -t457 * t855 + t540 * t458 - t480 * t856 + t539 * t481;
t885 = qJD(1) * t936 - t417 * t934 + t418 * t933;
t474 = qJD(5) + t509;
t379 = t474 * t516;
t694 = qJD(1) * t514;
t407 = rSges(6,2) * t775;
t693 = qJD(1) * t516;
t711 = rSges(6,3) * t693 + qJD(1) * t407;
t827 = rSges(6,2) * t448;
t145 = -t379 * t827 + (-t379 * t447 - t448 * t694) * rSges(6,1) + t711;
t351 = rSges(6,1) * t447 + t827;
t285 = t351 * t514;
t829 = rSges(6,1) * t448;
t352 = -rSges(6,2) * t447 + t829;
t496 = t514 * rSges(6,3);
t146 = -t474 * t285 + (t352 * t516 + t496) * qJD(1);
t252 = rSges(6,1) * t773 - t516 * rSges(6,3) - t407;
t286 = t351 * t516;
t378 = qJD(5) * t514 + t417;
t884 = t516 * t145 + t514 * t146 + t252 * t693 + t378 * t285 + t286 * t379;
t883 = 0.2e1 * qJD(2);
t844 = t417 / 0.2e1;
t843 = -t418 / 0.2e1;
t513 = sin(qJ(2));
t882 = rSges(3,2) * t513;
t837 = pkin(4) * t457;
t838 = pkin(3) * t480;
t382 = -t837 - t838;
t826 = pkin(2) * qJD(2);
t682 = t513 * t826;
t261 = t382 * t509 - t682;
t230 = t516 * t261;
t369 = -t509 * t838 - t682;
t339 = t516 * t369;
t517 = -pkin(7) - pkin(6);
t508 = -qJ(4) + t517;
t482 = -pkin(8) + t508;
t701 = -t482 + t508;
t465 = pkin(3) * t481;
t515 = cos(qJ(2));
t467 = t515 * pkin(2) + pkin(1);
t400 = t465 + t467;
t836 = pkin(4) * t458;
t326 = t400 + t836;
t721 = t326 - t400;
t103 = t230 - t339 + (-t514 * t721 + t516 * t701) * qJD(1);
t246 = t352 * t474;
t687 = qJD(1) * qJD(2);
t464 = t516 * t687;
t686 = qJD(1) * qJD(3);
t394 = t516 * t686 + t464;
t684 = qJD(1) * qJD(5);
t343 = t516 * t684 + t394;
t692 = qJD(2) * t516;
t661 = t513 * t692;
t613 = pkin(2) * t661;
t755 = t515 * qJD(2) ^ 2;
t680 = t514 * t755;
t472 = t516 * t517;
t473 = pkin(6) * t693;
t833 = pkin(1) - t467;
t209 = -t613 - t473 + (t514 * t833 - t472) * qJD(1);
t383 = qJD(1) * (-pkin(1) * t694 + t473);
t747 = qJD(1) * t209 + t383;
t769 = t458 * t509;
t700 = -t508 + t517;
t710 = t400 - t467;
t476 = qJD(4) * t514;
t717 = t339 + t476;
t124 = t613 + (-t514 * t710 + t516 * t700) * qJD(1) + t717;
t685 = qJD(1) * qJD(4);
t764 = t481 * t509;
t873 = qJD(1) * t124 + t514 * t685 + (-t394 * t480 - t417 * t764) * pkin(3);
t19 = -pkin(2) * t680 - t246 * t378 - t343 * t351 + (-t394 * t457 - t417 * t769) * pkin(4) + (t103 + t145 - t613) * qJD(1) + t747 + t873;
t881 = t19 * t514;
t364 = rSges(5,1) * t457 + rSges(5,2) * t458;
t317 = t364 * t514;
t830 = rSges(5,1) * t458;
t365 = -rSges(5,2) * t457 + t830;
t497 = t514 * rSges(5,3);
t160 = -t509 * t317 + (t365 * t516 + t497) * qJD(1);
t302 = t365 * t509;
t463 = t514 * t687;
t393 = t514 * t686 + t463;
t368 = t418 * t465;
t443 = t514 * t682;
t578 = -pkin(2) * t516 * t755 + qJD(1) * t443;
t551 = -t509 * t368 + t393 * t838 + t516 * t685 + t578;
t442 = t508 * t694;
t477 = qJD(4) * t516;
t605 = t369 * t514 - t442 - t477;
t703 = t517 * t694 + t443;
t125 = t693 * t710 + t605 + t703;
t504 = t514 * pkin(6);
t210 = (-t516 * t833 - t504) * qJD(1) - t703;
t433 = t516 * pkin(1) + t504;
t401 = t433 * qJD(1);
t744 = -t210 - t401;
t678 = -t125 + t744;
t47 = -t302 * t418 + t364 * t393 + (-t160 + t678) * qJD(1) + t551;
t366 = t417 * t838;
t668 = -t366 - t443 - t477;
t438 = t516 * t467;
t615 = -t514 * t517 + t438;
t289 = t615 - t433;
t371 = t516 * t400;
t212 = t514 * t700 + t371 - t438;
t274 = rSges(5,1) * t767 - rSges(5,2) * t770 + t497;
t743 = t212 + t274;
t672 = -t289 - t743;
t91 = -t364 * t417 + (t433 - t672) * qJD(1) + t668;
t880 = qJD(1) * t91 + t47;
t713 = t514 * t400 + t516 * t508;
t723 = -t514 * t326 - t516 * t482;
t167 = t713 + t723;
t772 = t448 * t516;
t774 = t447 * t516;
t253 = rSges(6,1) * t772 - rSges(6,2) * t774 + t496;
t705 = t514 * t467 + t472;
t211 = t705 - t713;
t506 = t516 * pkin(6);
t432 = pkin(1) * t514 - t506;
t288 = t432 - t705;
t734 = -t288 * t479 + t289 * t692;
t670 = -t417 * t211 + t734;
t308 = t516 * t326;
t168 = t514 * t701 + t308 - t371;
t749 = t168 + t212;
t36 = -t167 * t417 + t252 * t378 + t253 * t379 + t418 * t749 + t670;
t879 = t36 * t694;
t409 = qJD(1) * t432;
t878 = qJD(1) * t288 - t409;
t495 = Icges(3,4) * t515;
t596 = -Icges(3,2) * t513 + t495;
t424 = Icges(3,1) * t513 + t495;
t761 = t482 * t514;
t104 = t442 + (t261 - t369) * t514 + (t516 * t721 - t761) * qJD(1);
t679 = t516 * t124 + t514 * t125 - t211 * t693;
t877 = t516 * t103 + t514 * t104 - t167 * t693 + t679 + t884;
t876 = t932 + t937;
t218 = t379 * t352;
t666 = t480 * t694;
t419 = pkin(3) * t666;
t667 = t457 * t694;
t778 = t418 * t458;
t875 = t351 * t694 + t218 + t368 + t419 + (t667 + t778) * pkin(4);
t303 = rSges(4,1) * t763 - rSges(4,2) * t766 - t516 * rSges(4,3);
t828 = rSges(4,2) * t481;
t390 = rSges(4,1) * t480 + t828;
t561 = -t390 * t418 - t613;
t727 = t288 - t432;
t109 = (-t303 + t727) * qJD(1) + t561;
t498 = t514 * rSges(4,3);
t304 = rSges(4,1) * t762 - rSges(4,2) * t765 + t498;
t726 = -t289 - t304;
t110 = -t390 * t417 - t443 + (t433 - t726) * qJD(1);
t355 = t390 * t514;
t356 = t390 * t516;
t831 = rSges(4,1) * t481;
t391 = -rSges(4,2) * t480 + t831;
t96 = t303 * t417 + t304 * t418 + t734;
t872 = -t109 * (qJD(1) * t355 - t418 * t391) - t110 * (-qJD(1) * t356 - t391 * t417) - t96 * (-t417 * t355 - t356 * t418);
t757 = t514 * t515;
t759 = t513 * t514;
t804 = Icges(3,3) * t516;
t327 = Icges(3,5) * t757 - Icges(3,6) * t759 - t804;
t452 = Icges(3,4) * t759;
t816 = Icges(3,5) * t516;
t331 = Icges(3,1) * t757 - t452 - t816;
t808 = Icges(3,6) * t516;
t329 = Icges(3,4) * t757 - Icges(3,2) * t759 - t808;
t789 = t329 * t513;
t584 = -t331 * t515 + t789;
t129 = -t516 * t327 - t514 * t584;
t334 = Icges(6,5) * t448 - Icges(6,6) * t447;
t333 = Icges(6,5) * t447 + Icges(6,6) * t448;
t786 = t333 * t516;
t792 = t240 * t447;
t801 = Icges(6,3) * t516;
t871 = -t474 * t786 + (-t242 * t448 - t334 * t514 + t792 + t801) * qJD(1);
t405 = Icges(6,4) * t775;
t813 = Icges(6,5) * t516;
t241 = Icges(6,1) * t773 - t405 - t813;
t590 = t239 * t447 - t241 * t448;
t238 = Icges(6,3) * t514 + t334 * t516;
t699 = qJD(1) * t238;
t787 = t333 * t514;
t870 = qJD(1) * t590 - t474 * t787 + t699;
t421 = Icges(3,5) * t515 - Icges(3,6) * t513;
t420 = Icges(3,5) * t513 + Icges(3,6) * t515;
t564 = qJD(2) * t420;
t820 = Icges(3,4) * t513;
t425 = Icges(3,1) * t515 - t820;
t332 = Icges(3,5) * t514 + t425 * t516;
t330 = Icges(3,6) * t514 + t516 * t596;
t788 = t330 * t513;
t583 = -t332 * t515 + t788;
t865 = -t516 * t564 + (-t421 * t514 + t583 + t804) * qJD(1);
t328 = Icges(3,3) * t514 + t421 * t516;
t696 = qJD(1) * t328;
t864 = qJD(1) * t584 - t514 * t564 + t696;
t582 = t335 * t447 - t337 * t448;
t863 = qJD(1) * t582 + t334 * t474;
t422 = Icges(3,2) * t515 + t820;
t579 = t422 * t513 - t424 * t515;
t860 = t579 * qJD(1) + t421 * qJD(2);
t859 = t514 * (-t422 * t516 + t332) - t516 * (-Icges(3,2) * t757 + t331 - t452);
t342 = t514 * t684 + t393;
t20 = -t246 * t379 + t342 * t351 + (t393 * t457 - t418 * t769) * pkin(4) + (-t104 - t146 + t678) * qJD(1) + t551;
t675 = -t253 - t749;
t612 = -t289 + t675;
t68 = -t417 * t837 - t351 * t378 + (t433 - t612) * qJD(1) + t668;
t858 = (qJD(1) * t68 + t20) * t516 + t881;
t857 = qJD(1) * t891 + t378 * t898 - t379 * (-Icges(6,2) * t773 + t241 - t405);
t854 = t514 * t843 + t516 * t844;
t853 = t342 / 0.2e1;
t852 = t343 / 0.2e1;
t851 = -t378 / 0.2e1;
t850 = t378 / 0.2e1;
t849 = -t379 / 0.2e1;
t848 = t379 / 0.2e1;
t847 = t393 / 0.2e1;
t846 = t394 / 0.2e1;
t845 = -t417 / 0.2e1;
t842 = t418 / 0.2e1;
t841 = t514 / 0.2e1;
t840 = -t516 / 0.2e1;
t839 = pkin(2) * t513;
t835 = -qJD(1) / 0.2e1;
t834 = qJD(1) / 0.2e1;
t832 = rSges(3,1) * t515;
t499 = t514 * rSges(3,3);
t367 = t418 * t838;
t577 = t476 - t613;
t571 = -t367 + t577;
t673 = t211 + t727;
t750 = t167 - t252;
t782 = t379 * t351;
t67 = -t418 * t837 - t782 + (t673 + t750) * qJD(1) + t571;
t825 = t514 * t67;
t273 = rSges(5,1) * t768 - rSges(5,2) * t771 - t516 * rSges(5,3);
t73 = t273 * t417 + t418 * t743 + t670;
t798 = qJD(1) * t73;
t796 = t109 * t514;
t702 = rSges(3,2) * t759 + t516 * rSges(3,3);
t340 = rSges(3,1) * t757 - t702;
t427 = rSges(3,1) * t513 + rSges(3,2) * t515;
t662 = t427 * t692;
t182 = -t662 + (-t340 - t432) * qJD(1);
t795 = t182 * t514;
t794 = t182 * t516;
t663 = t427 * t479;
t756 = t515 * t516;
t758 = t513 * t516;
t341 = rSges(3,1) * t756 - rSges(3,2) * t758 + t499;
t716 = t341 + t433;
t183 = qJD(1) * t716 - t663;
t381 = t427 * t516;
t793 = t183 * t381;
t785 = t335 * t474;
t779 = t417 * t481;
t777 = t420 * t514;
t776 = t420 * t516;
t237 = Icges(6,5) * t773 - Icges(6,6) * t775 - t801;
t754 = t516 * t237;
t748 = -t514 * t211 + t516 * t212;
t746 = -t514 * t237 - t241 * t772;
t745 = t514 * t238 + t242 * t772;
t740 = t514 * t252 + t516 * t253;
t733 = -t514 * t288 + t516 * t289;
t730 = t514 * t303 + t516 * t304;
t729 = -t514 * t327 - t331 * t756;
t728 = t514 * t328 + t332 * t756;
t722 = t364 * t694 + t419;
t354 = t382 - t839;
t410 = -t838 - t839;
t715 = t354 - t410;
t709 = rSges(5,2) * t667 + rSges(5,3) * t693;
t708 = rSges(4,2) * t666 + rSges(4,3) * t693;
t707 = -t422 + t425;
t706 = t424 + t596;
t704 = rSges(3,3) * t693 + t694 * t882;
t695 = qJD(1) * t421;
t119 = -t514 * t582 - t786;
t691 = t119 * qJD(1);
t186 = -t514 * t579 - t776;
t688 = t186 * qJD(1);
t683 = pkin(3) * t764;
t681 = t515 * t826;
t175 = -t418 * t828 + (-t418 * t480 - t481 * t694) * rSges(4,1) + t708;
t176 = -t509 * t355 + (t391 * t516 + t498) * qJD(1);
t676 = t516 * t175 + t514 * t176 + t303 * t693;
t674 = t516 * t209 + t514 * t210 - t288 * t693;
t664 = t513 * t693;
t660 = t515 * t692;
t658 = -pkin(1) - t832;
t657 = t694 / 0.2e1;
t656 = t693 / 0.2e1;
t655 = -t479 / 0.2e1;
t652 = t692 / 0.2e1;
t651 = -t390 - t839;
t650 = t410 + t839;
t649 = -t364 - t838;
t648 = t382 + t838;
t646 = t513 * (-t514 ^ 2 - t516 ^ 2);
t645 = (-t514 * t593 + t805) * qJD(1) + t898 * t474;
t644 = qJD(1) * t240 + t241 * t474 - t514 * t785;
t643 = (-t338 * t514 + t813) * qJD(1) + t899 * t474;
t642 = qJD(1) * t242 + t474 * t900;
t205 = t242 * t773;
t633 = t516 * t238 - t205;
t318 = t364 * t516;
t630 = -t417 * t317 - t318 * t418;
t628 = -t237 + t792;
t627 = t891 * t474;
t626 = t338 * t474 - t785;
t624 = -qJD(1) * t286 - t352 * t378;
t276 = t332 * t757;
t623 = t516 * t328 - t276;
t620 = -qJD(1) * t318 - t365 * t417;
t616 = -t327 + t788;
t611 = t514 * t273 + t516 * t274 + t748;
t316 = t650 * t516;
t610 = -pkin(3) * t779 + qJD(1) * t316;
t609 = -t302 - t683;
t324 = t391 * t509;
t606 = -t324 - t681;
t604 = qJD(1) * t317 - t418 * t365 - t368;
t602 = t832 - t882;
t601 = t68 * t624;
t592 = -t109 * t516 - t110 * t514;
t591 = -t183 * t514 - t794;
t127 = t239 * t448 + t241 * t447;
t184 = t329 * t515 + t331 * t513;
t185 = t330 * t515 + t332 * t513;
t576 = -t364 + t410;
t159 = -rSges(5,2) * t778 + (-t418 * t457 - t458 * t694) * rSges(5,1) + t709;
t573 = t516 * t159 + t514 * t160 + t273 * t693 + t679;
t572 = -t514 * t167 + t516 * t168 + t740 + t748;
t380 = t427 * t514;
t569 = t590 * t514;
t566 = qJD(2) * t424;
t565 = qJD(2) * t422;
t130 = -t330 * t759 - t623;
t563 = (-t129 * t516 + t130 * t514) * qJD(2);
t131 = -t329 * t758 - t729;
t132 = -t330 * t758 + t728;
t562 = (-t131 * t516 + t132 * t514) * qJD(2);
t180 = (t340 * t514 + t341 * t516) * qJD(2);
t559 = -pkin(4) * t769 - t246 - t683;
t558 = t609 - t681;
t557 = qJD(1) * t334 - t378 * t786 + t379 * t787;
t553 = t209 * t692 + t210 * t479 - t288 * t464 - t289 * t463;
t552 = t329 * t516 - t330 * t514;
t550 = qJD(1) * t211 + t577 + t878;
t100 = -t240 * t774 + t745;
t128 = t240 * t448 + t242 * t447;
t536 = qJD(1) * t237 - t447 * t644 + t448 * t642;
t24 = t514 * t870 + t536 * t516;
t535 = -t447 * t645 + t448 * t643 + t699;
t25 = t514 * t871 + t535 * t516;
t26 = t536 * t514 - t516 * t870;
t27 = t535 * t514 - t516 * t871;
t97 = -t569 - t754;
t98 = -t240 * t775 - t633;
t40 = t378 * t98 - t379 * t97 + t691;
t120 = -t516 * t582 + t787;
t118 = t120 * qJD(1);
t99 = -t239 * t774 - t746;
t41 = t100 * t378 - t379 * t99 + t118;
t541 = t899 * t378 - t900 * t379 + (-t335 + t338) * qJD(1);
t523 = -t447 * t857 + t541 * t448;
t530 = qJD(1) * t333 - t447 * t627 + t448 * t626;
t59 = t514 * t863 + t530 * t516;
t60 = t530 * t514 - t516 * t863;
t69 = t447 * t642 + t448 * t644;
t70 = t447 * t643 + t448 * t645;
t548 = (qJD(1) * t59 + t100 * t343 - t24 * t379 + t25 * t378 + t342 * t99) * t841 + (t514 * t70 - t516 * t69 + (t127 * t514 + t128 * t516) * qJD(1)) * t834 + (t514 * t557 + t516 * t523) * t851 + (t514 * t523 - t516 * t557) * t848 + (qJD(1) * t60 - t26 * t379 + t27 * t378 + t342 * t97 + t343 * t98) * t840 + (t541 * t447 + t448 * t857) * t835 + t40 * t657 + t41 * t656 + (-t24 * t516 + t25 * t514 + (t100 * t516 + t514 * t99) * qJD(1)) * t850 + (-t26 * t516 + t27 * t514 + (t514 * t97 + t516 * t98) * qJD(1)) * t849 + (t514 * t98 - t516 * t97) * t853 + (t100 * t514 - t516 * t99) * t852;
t547 = t417 * t125 - t394 * t211 + t553;
t546 = (-t513 * t706 + t515 * t707) * qJD(1);
t545 = t559 - t681;
t543 = (-t464 * t513 - t680) * pkin(2) + t747;
t198 = qJD(1) * t330 - t514 * t565;
t200 = qJD(1) * t332 - t514 * t566;
t527 = qJD(1) * t327 - qJD(2) * t184 - t198 * t513 + t200 * t515;
t197 = -t516 * t565 + (-t514 * t596 + t808) * qJD(1);
t199 = -t516 * t566 + (-t425 * t514 + t816) * qJD(1);
t526 = -qJD(2) * t185 - t197 * t513 + t199 * t515 + t696;
t396 = t596 * qJD(2);
t397 = t425 * qJD(2);
t525 = qJD(1) * t420 - t396 * t513 + t397 * t515 + (-t422 * t515 - t424 * t513) * qJD(2);
t524 = -t513 * t859 + t552 * t515;
t520 = t548 + (t514 * t905 - t516 * t906) * t847 + (t514 * t903 - t516 * t904) * t846 + (t514 * t885 + t516 * t886) * t845 + (t916 * t516 + t915 * t514 + (t514 * t904 + t516 * t903) * qJD(1)) * t844 + (t914 * t516 + t913 * t514 + (t514 * t906 + t516 * t905) * qJD(1)) * t843 + (t514 * t886 - t516 * t885) * t842 + (t910 * qJD(1) + t904 * t393 + t903 * t394 + t915 * t417 + t418 * t916) * t841 + (qJD(1) * t909 + t393 * t906 + t394 * t905 + t417 * t913 + t418 * t914) * t840 + (t540 * t457 + t458 * t855 + t539 * t480 + t481 * t856) * t835 + (t908 * t516 + t907 * t514 + (t514 * t902 + t516 * t901) * qJD(1)) * t834 + t912 * t657 + t911 * t656;
t398 = t602 * qJD(2);
t353 = t516 * t367;
t315 = t650 * t514;
t257 = t648 * t516;
t256 = t648 * t514;
t220 = t417 * t315;
t214 = t715 * t516;
t213 = t715 * t514;
t202 = -qJD(2) * t380 + (t516 * t602 + t499) * qJD(1);
t201 = -rSges(3,2) * t660 + (-t515 * t694 - t661) * rSges(3,1) + t704;
t187 = -t516 * t579 + t777;
t181 = t187 * qJD(1);
t116 = -t398 * t692 + (-t202 - t401 + t663) * qJD(1);
t115 = -t398 * t479 + t383 + (t201 - t662) * qJD(1);
t95 = t525 * t514 - t516 * t860;
t94 = t514 * t860 + t525 * t516;
t93 = -qJD(2) * t583 + t197 * t515 + t199 * t513;
t92 = -qJD(2) * t584 + t198 * t515 + t200 * t513;
t90 = -t364 * t418 + (-t273 + t673) * qJD(1) + t571;
t85 = -t324 * t418 + t390 * t393 + (-t176 + t744) * qJD(1) + t578;
t84 = qJD(1) * t175 - t324 * t417 - t390 * t394 + t543;
t77 = t181 + t562;
t76 = t563 + t688;
t46 = qJD(1) * t159 - t302 * t417 - t364 * t394 + t543 + t873;
t39 = t175 * t418 + t176 * t417 + t303 * t394 - t304 * t393 + t553;
t14 = t160 * t417 + t273 * t394 + (t124 + t159) * t418 - t743 * t393 + t547;
t9 = t104 * t417 + t145 * t379 + t146 * t378 - t167 * t394 + t252 * t343 - t253 * t342 + (t103 + t124) * t418 - t749 * t393 + t547;
t1 = [(t181 + ((t130 - t276 + (t328 + t789) * t516 + t729) * t516 + t728 * t514) * qJD(2)) * t652 + (t118 + (t98 + (t239 * t516 + t240 * t514) * t447 + t633 + t746) * t379 + (-t241 * t773 + t754 + t97 + (t239 * t514 - t240 * t516) * t447 + t745) * t378) * t848 + (t119 + t127) * t853 + (t120 + t128) * t852 + (-t691 + (t100 - t569 - t745) * t379 + (t514 * t628 - t205 + t99) * t378 + ((t238 + t590) * t378 + t628 * t379) * t516 + t40) * t851 + (t70 + t59) * t850 + (((t298 * t516 + t299 * t514) * t480 + (t266 * t516 + t267 * t514) * t457 + t905 + t926 + t928) * t418 + (-t300 * t763 + (t298 * t514 - t299 * t516) * t480 - t268 * t768 + (t266 * t514 - t267 * t516) * t457 - t874 + t906 - t925) * t417 + t918) * t842 + (t76 - t688 + ((t516 * t616 + t132 - t728) * t516 + (t514 * t616 + t131 + t623) * t514) * qJD(2)) * t655 + (t93 + t94) * t479 / 0.2e1 + (t20 * (-t252 + t723) + t67 * (-t261 * t514 + t477) + t19 * (t253 + t308 - t761) + t351 * t825 * t474 + (-t286 * t474 - t382 * t418 + t230 + t476 - t550 + t67 + t711 + t782) * t68) * m(6) + (-qJD(2) * t579 + t396 * t515 + t397 * t513 + t447 * t626 + t448 * t627 + t457 * t621 + t458 * t622 + t480 * t617 + t481 * t618 + (-t750 * t68 + (t67 * (-t326 - t352) - t68 * t482) * t516 + (t67 * (-rSges(6,3) + t482) + t68 * (-t326 - t829)) * t514) * m(6)) * qJD(1) + (-(-qJD(1) * t273 + t418 * t649 + t550 - t90) * t91 + t47 * (-t273 - t713) - t90 * t605 + t46 * (-t508 * t514 + t274 + t371) + t91 * (t709 + t717) + (t317 * t90 - t318 * t91) * t509 + ((-t90 * rSges(5,3) + t91 * (-t400 - t830)) * t514 + (t90 * (-t365 - t400) - t91 * t508) * t516) * qJD(1)) * m(5) + (-(-qJD(1) * t303 - t109 + t561 + t878) * t110 + t85 * (-t303 - t705) + t109 * t703 + t84 * (t304 + t615) + t110 * (-t613 + t708) + (-t110 * t356 + t390 * t796) * t509 + ((-t109 * rSges(4,3) + t110 * (-t467 - t831)) * t514 + (t109 * (-t391 - t467) - t110 * t517) * t516) * qJD(1)) * m(4) + (t116 * (t514 * t658 + t506 + t702) + t115 * t716 + t183 * (t473 + t704) + (t427 * t795 - t793) * qJD(2) + ((-pkin(1) - t602) * t794 + (t182 * (-rSges(3,3) - pkin(6)) + t183 * t658) * t514) * qJD(1) - (-qJD(1) * t340 - t182 - t409 - t662) * t183) * m(3) + (t69 + t60 + t41) * t849 - (t92 + t95 + t77) * t692 / 0.2e1 + (t902 - t931) * t847 + (t901 + t930) * t846 + ((t516 * t876 + t874 + t903 - t929) * t418 + ((t935 + t940) * t516 + t876 * t514 + t904 + t938) * t417 + t912 + t917) * t845 + (t907 + t910) * t844 + (-t908 + t909 + t911) * t843 + ((t184 + t186) * t514 + (t185 + t187) * t516) * t687 / 0.2e1; t520 + ((t513 * t707 + t515 * t706) * qJD(1) + (t552 * t513 + t515 * t859) * qJD(2)) * t835 + ((-t479 * t776 + t695) * t514 + (t546 + (t514 * t777 + t524) * qJD(2)) * t516) * t655 + ((-t692 * t777 - t695) * t516 + (t546 + (t516 * t776 + t524) * qJD(2)) * t514) * t652 + (t514 * t93 - t516 * t92 + (t184 * t514 + t185 * t516) * qJD(1)) * t834 + (qJD(1) * t94 + (t514 * (t514 * t865 + t526 * t516) - t516 * (t514 * t864 + t527 * t516) + (t131 * t514 + t132 * t516) * qJD(1)) * t883) * t841 + (qJD(1) * t95 + (t514 * (t526 * t514 - t516 * t865) - t516 * (t527 * t514 - t516 * t864) + (t129 * t514 + t130 * t516) * qJD(1)) * t883) * t840 + (t563 + t76) * t657 + (t77 + t562) * t656 + (t9 * (t572 + t733) + t612 * t879 + (t417 * t836 - t610 - t624 - (-pkin(2) * t758 + t214) * qJD(1) + (t681 + t545) * t514) * t68 + t858 * (-t351 + t354) + (pkin(2) * t660 - (-t213 + t285 - t315) * qJD(1) + t545 * t516 + t875) * t67 + (-t213 * t417 - t220 - (t214 + t316) * t418 - t646 * t826 + t674 + t877) * t36) * m(6) + (-t90 * (-qJD(1) * t315 + t604) - t91 * (t610 + t620) - t73 * (t316 * t418 + t220 + t630) - (-t91 * t664 + ((-t514 * t91 - t516 * t90) * t515 + t73 * t646) * qJD(2)) * pkin(2) + t90 * t722 + t14 * (t611 + t733) + t73 * (t573 + t674) + (t558 * t90 + t576 * t880) * t516 + (t46 * t576 + t558 * t91 + t672 * t798) * t514) * m(5) + (-(-t110 * t664 + (t515 * t592 + t646 * t96) * qJD(2)) * pkin(2) + t39 * (t730 + t733) + t96 * (t674 + t676) + (t109 * t606 + (qJD(1) * t110 + t85) * t651) * t516 + (t84 * t651 + t110 * t606 + (t109 * t390 + t726 * t96) * qJD(1)) * t514 + t872) * m(4) + (-(t182 * t380 - t793) * qJD(1) - (t180 * (-t380 * t514 - t381 * t516) + t591 * t602) * qJD(2) + 0.2e1 * t180 * (t201 * t516 + t202 * t514 + (t340 * t516 - t341 * t514) * qJD(1)) + t591 * t398 + (-t115 * t514 - t116 * t516 + (-t183 * t516 + t795) * qJD(1)) * t427) * m(3); t520 + (t675 * t879 + t9 * t572 - t601 + (-(-t836 - t465) * t417 - (-pkin(3) * t765 + t257) * qJD(1) + t559 * t514) * t68 + t858 * (-t351 + t382) + (-t419 - (-t256 + t285) * qJD(1) + t559 * t516 + t875) * t67 + (-t257 * t418 + t353 - (-pkin(3) * t766 + t256) * t417 + t877) * t36) * m(6) + (-t91 * ((-t480 * t693 - t779) * pkin(3) + t620) + t14 * t611 + (t609 * t91 - t743 * t798) * t514 + (t46 * t514 + t516 * t880) * t649 + (t516 * t609 - t419 - t604 + t722) * t90 + (t366 * t514 + t353 + t573 - t630) * t73) * m(5) + (t39 * t730 + t96 * (-t304 * t694 + t676) + t592 * t324 + (-t84 * t514 - t85 * t516 + (-t110 * t516 + t796) * qJD(1)) * t390 + t872) * m(4); 0.2e1 * (t19 * t840 + t20 * t841 + t36 * t854) * m(6) + 0.2e1 * (t46 * t840 + t47 * t841 + t73 * t854) * m(5); t548 + (t9 * t740 + (-t514 * t68 - t516 * t67) * t246 + (-t881 - t20 * t516 + (-t516 * t68 + t825) * qJD(1)) * t351 - t67 * (qJD(1) * t285 - t218) - t601 + (-t253 * t694 + t884) * t36) * m(6);];
tauc = t1(:);
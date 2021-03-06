% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR7_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:56:14
% EndTime: 2019-03-09 02:57:14
% DurationCPUTime: 56.61s
% Computational Cost: add. (19169->1089), mult. (29411->1357), div. (0->0), fcn. (25448->8), ass. (0->545)
t947 = -Icges(4,3) - Icges(5,3);
t470 = qJ(3) + pkin(9);
t449 = sin(t470);
t450 = cos(t470);
t475 = sin(qJ(3));
t478 = cos(qJ(3));
t945 = Icges(4,5) * t475 + Icges(5,5) * t449 + Icges(4,6) * t478 + Icges(5,6) * t450;
t442 = Icges(6,6) * t449;
t348 = -Icges(6,2) * t450 + t442;
t777 = Icges(5,4) * t449;
t356 = Icges(5,1) * t450 - t777;
t946 = t348 - t356;
t776 = Icges(5,4) * t450;
t354 = -Icges(5,2) * t449 + t776;
t766 = Icges(6,6) * t450;
t585 = -Icges(6,3) * t449 + t766;
t924 = -t354 - t585;
t476 = sin(qJ(1));
t479 = cos(qJ(1));
t939 = -t945 * t476 + t479 * t947;
t351 = Icges(6,4) * t449 + Icges(6,5) * t450;
t250 = Icges(6,1) * t476 + t351 * t479;
t884 = t476 * t947 + t945 * t479;
t944 = t250 - t884;
t347 = Icges(6,2) * t449 + t766;
t595 = Icges(5,1) * t449 + t776;
t943 = -t347 - t595;
t942 = Icges(4,5) * t478 + Icges(5,5) * t450 - Icges(4,6) * t475 - Icges(5,6) * t449;
t584 = Icges(6,3) * t450 + t442;
t592 = Icges(5,2) * t450 + t777;
t941 = t584 + t592;
t243 = -Icges(5,6) * t476 + t479 * t592;
t245 = -Icges(5,5) * t476 + t479 * t595;
t779 = Icges(4,4) * t475;
t593 = Icges(4,2) * t478 + t779;
t276 = -Icges(4,6) * t476 + t479 * t593;
t778 = Icges(4,4) * t478;
t596 = Icges(4,1) * t475 + t778;
t278 = -Icges(4,5) * t476 + t479 * t596;
t246 = Icges(6,5) * t476 + t479 * t584;
t248 = Icges(6,4) * t476 + t347 * t479;
t575 = t246 * t450 + t248 * t449;
t913 = -t243 * t450 - t245 * t449 - t276 * t478 - t278 * t475 - t575;
t591 = Icges(6,4) * t450 - Icges(6,5) * t449;
t940 = t591 - t942;
t242 = Icges(5,6) * t479 + t476 * t592;
t247 = Icges(6,5) * t479 - t476 * t584;
t938 = t242 - t247;
t937 = t243 + t246;
t755 = t450 * t476;
t413 = Icges(5,4) * t755;
t757 = t449 * t476;
t770 = Icges(5,5) * t479;
t244 = Icges(5,1) * t757 + t413 + t770;
t249 = Icges(6,4) * t479 - t347 * t476;
t936 = -t244 + t249;
t935 = t245 + t248;
t934 = t924 * t479;
t933 = t946 * t476;
t932 = t946 * t479;
t397 = -Icges(4,2) * t475 + t778;
t399 = Icges(4,1) * t478 - t779;
t931 = t354 * t450 + t356 * t449 + t397 * t478 + t399 * t475;
t275 = Icges(4,6) * t479 + t476 * t593;
t748 = t476 * t478;
t435 = Icges(4,4) * t748;
t751 = t475 * t476;
t771 = Icges(4,5) * t479;
t277 = Icges(4,1) * t751 + t435 + t771;
t883 = t242 * t450 + t244 * t449 + t275 * t478 + t277 * t475;
t570 = t348 * t449 - t450 * t585;
t921 = -t570 + t931;
t885 = t939 * t476;
t902 = t479 * t883 + t885;
t251 = Icges(6,1) * t479 - t351 * t476;
t754 = t450 * t479;
t756 = t449 * t479;
t97 = t247 * t754 + t249 * t756 + t476 * t251;
t930 = t97 - t902;
t903 = -t242 * t755 - t244 * t757 - t275 * t748 - t277 * t751 + t479 * t939;
t224 = t479 * t251;
t574 = t247 * t450 + t249 * t449;
t99 = -t476 * t574 + t224;
t929 = t99 - t903;
t928 = t942 * t476;
t927 = t941 * qJD(3);
t926 = t943 * qJD(3);
t923 = t574 - t883;
t922 = t351 - t945;
t865 = t940 * t479;
t920 = t913 * t479;
t919 = qJD(1) * t938 + qJD(3) * t934;
t279 = t585 * t476;
t693 = qJD(3) * t476;
t918 = -qJD(1) * t937 - qJD(3) * t279 - t354 * t693;
t917 = t932 * qJD(3) + (t476 * t595 - t249 + t770) * qJD(1);
t916 = qJD(1) * t935 - qJD(3) * t933;
t915 = t479 * t931 - t928;
t223 = t479 * t250;
t870 = -t243 * t755 - t245 * t757 - t276 * t748 - t278 * t751 - t476 * t575 - t479 * t884 + t223;
t899 = t476 * t944 - t920;
t914 = t476 * t921 - t865;
t869 = t276 * t475 - t278 * t478 + t449 * t937 - t450 * t935;
t287 = t591 * t476;
t868 = -t287 + t928;
t867 = t275 * t475 - t277 * t478 + t449 * t938 + t450 * t936;
t810 = pkin(3) * t478;
t912 = rSges(4,2) * t478;
t911 = (Icges(5,2) * t757 - t279 - t413 + t936) * t479 + (-t934 + t935) * t476;
t910 = (-t933 - t938) * t479 + (t932 + t937) * t476;
t909 = t924 + t943;
t908 = -t946 - t941;
t681 = qJD(1) * qJD(3);
t907 = t921 * qJD(1) + qJD(3) * t922;
t906 = t944 * qJD(1);
t905 = (t251 - t939) * qJD(1);
t368 = t593 * qJD(3);
t369 = t596 * qJD(3);
t904 = -t368 * t478 - t369 * t475 - t927 * t450 + t926 * t449 + (-t397 * t475 + t399 * t478 + t449 * t924 - t450 * t946) * qJD(3) + t940 * qJD(1);
t474 = sin(qJ(6));
t477 = cos(qJ(6));
t587 = Icges(7,5) * t474 + Icges(7,6) * t477;
t228 = Icges(7,3) * t450 + t449 * t587;
t229 = -Icges(7,3) * t449 + t450 * t587;
t687 = qJD(6) * t476;
t692 = qJD(3) * t479;
t342 = t449 * t687 + t692;
t877 = qJD(6) * t450;
t418 = qJD(1) + t877;
t686 = qJD(6) * t479;
t543 = t449 * t686 - t693;
t773 = Icges(7,4) * t474;
t590 = Icges(7,2) * t477 + t773;
t230 = Icges(7,6) * t450 + t449 * t590;
t772 = Icges(7,4) * t477;
t412 = t449 * t772;
t758 = t449 * t474;
t768 = Icges(7,5) * t450;
t232 = Icges(7,1) * t758 + t412 + t768;
t579 = t230 * t477 + t232 * t474;
t749 = t476 * t477;
t752 = t474 * t479;
t323 = -t450 * t749 - t752;
t747 = t477 * t479;
t753 = t474 * t476;
t324 = -t450 * t753 + t747;
t774 = Icges(7,4) * t324;
t162 = Icges(7,2) * t323 + Icges(7,6) * t757 + t774;
t296 = Icges(7,4) * t323;
t165 = Icges(7,1) * t324 + Icges(7,5) * t757 + t296;
t581 = t162 * t477 + t165 * t474;
t322 = t450 * t752 + t749;
t295 = Icges(7,4) * t322;
t321 = -t450 * t747 + t753;
t160 = -Icges(7,2) * t321 - Icges(7,6) * t756 + t295;
t294 = Icges(7,4) * t321;
t164 = -Icges(7,1) * t322 + Icges(7,5) * t756 + t294;
t582 = t160 * t477 - t164 * t474;
t829 = t543 * (-t228 * t479 + t582) - t342 * (t228 * t476 + t581) - t418 * (t229 + t579);
t901 = t829 * t449;
t715 = t397 + t596;
t716 = -t593 + t399;
t900 = (t449 * t909 + t450 * t908 - t475 * t715 + t478 * t716) * qJD(1);
t898 = t914 * qJD(1);
t897 = -qJD(1) * t923 + t868 * qJD(3) - t906;
t896 = qJD(1) * t913 + qJD(3) * t865 + t905;
t895 = (t899 * t476 + t479 * t930) * qJD(3);
t894 = (t870 * t476 + t479 * t929) * qJD(3);
t114 = t479 * t570 - t287;
t893 = (-t114 + t915) * qJD(1);
t178 = qJD(1) * t276 + t397 * t693;
t339 = t399 * t476;
t180 = qJD(1) * t278 + qJD(3) * t339;
t892 = qJD(3) * t867 - t178 * t478 - t180 * t475 - t449 * t916 + t450 * t918 + t905;
t338 = t397 * t479;
t177 = qJD(1) * t275 - qJD(3) * t338;
t340 = t399 * t479;
t179 = -qJD(3) * t340 + (t476 * t596 + t771) * qJD(1);
t891 = qJD(3) * t869 + t177 * t478 + t179 * t475 + t449 * t917 + t450 * t919 - t906;
t158 = -Icges(7,5) * t322 + Icges(7,6) * t321 + Icges(7,3) * t756;
t890 = -t323 * t160 + t164 * t324;
t48 = -t158 * t757 - t890;
t159 = Icges(7,5) * t324 + Icges(7,6) * t323 + Icges(7,3) * t757;
t49 = t159 * t757 + t323 * t162 + t324 * t165;
t82 = t228 * t757 + t230 * t323 + t232 * t324;
t13 = t342 * t49 + t82 * t418 - t48 * t543;
t60 = -t450 * t158 + t449 * t582;
t46 = t158 * t756 - t160 * t321 - t164 * t322;
t889 = 0.2e1 * qJD(3);
t607 = rSges(5,1) * t449 + rSges(5,2) * t450;
t886 = t479 * t607;
t604 = rSges(6,2) * t450 - rSges(6,3) * t449;
t304 = t604 * t479;
t471 = t476 ^ 2;
t707 = t479 ^ 2 + t471;
t537 = t707 * t810;
t811 = pkin(3) * t475;
t532 = t607 + t811;
t643 = t449 * t692;
t698 = qJD(1) * t476;
t654 = t450 * t698;
t882 = t643 + t654;
t649 = t450 * t693;
t697 = qJD(1) * t479;
t655 = t449 * t697;
t531 = t649 + t655;
t362 = pkin(4) * t450 + qJ(5) * t449;
t689 = qJD(5) * t479;
t750 = t476 * t362;
t881 = -qJD(1) * t750 + t362 * t698 + t449 * t689;
t438 = qJD(5) * t449;
t762 = qJ(5) * t450;
t809 = pkin(4) * t449;
t358 = t762 - t809;
t878 = qJD(3) * t358;
t254 = t438 + t878;
t303 = t362 * t479;
t690 = qJD(5) * t476;
t880 = -qJD(1) * t303 + t476 * t254 - t358 * t693 + t362 * t697 - t449 * t690;
t645 = t450 * t689;
t656 = t449 * t698;
t648 = t450 * t692;
t663 = pkin(4) * t648 + qJ(5) * t882;
t131 = pkin(4) * t656 + t645 - t663;
t473 = -qJ(4) - pkin(7);
t434 = t473 * t697;
t440 = pkin(3) * t751;
t428 = t692 * t810;
t691 = qJD(4) * t476;
t624 = -t428 + t691;
t807 = pkin(7) * t479;
t197 = -t434 + (t440 - t807) * qJD(1) + t624;
t194 = t479 * t197;
t439 = qJD(5) * t450;
t879 = t479 * t131 + t303 * t692 + t194 - t439;
t688 = qJD(6) * t449;
t876 = t894 + t898;
t875 = -t893 + t895;
t874 = t476 * t907 - t479 * t904;
t873 = t476 * t904 + t479 * t907;
t872 = qJD(3) * t923 - t178 * t475 + t180 * t478 + t918 * t449 + t916 * t450;
t871 = -qJD(3) * t913 - t177 * t475 + t179 * t478 - t449 * t919 + t450 * t917;
t866 = t922 * qJD(1);
t166 = rSges(7,1) * t322 - rSges(7,2) * t321 - rSges(7,3) * t756;
t385 = pkin(8) * t649;
t863 = -t166 * t418 + t385;
t499 = t228 * t756 + t230 * t321 - t232 * t322;
t862 = t499 * t418 + t543 * t46;
t510 = t476 * (t278 + t338) - t479 * (-Icges(4,2) * t751 + t277 + t435);
t511 = t476 * (t276 - t340) - t479 * (t275 - t339);
t861 = -t449 * t910 + t450 * t911 - t511 * t475 + t510 * t478;
t852 = (-qJD(3) * t604 - t439) * t479;
t849 = t321 * t162 - t322 * t165;
t696 = qJD(3) * t449;
t674 = rSges(5,2) * t696;
t552 = -rSges(5,3) * qJD(1) - t674;
t653 = t450 * t697;
t662 = rSges(5,1) * t531 + rSges(5,2) * t653;
t154 = t476 * t552 + t662;
t652 = t475 * t697;
t451 = qJD(4) * t479;
t646 = t478 * t693;
t713 = pkin(3) * t646 + t451;
t618 = pkin(3) * t652 + t473 * t698 + t713;
t196 = pkin(7) * t698 + t618;
t848 = -t154 - t196;
t409 = t479 * pkin(1) + t476 * qJ(2);
t453 = qJD(2) * t479;
t332 = qJD(1) * t409 - t453;
t847 = -t197 - t332;
t799 = pkin(7) + t473;
t318 = t476 * t799 + t479 * t811;
t455 = t479 * qJ(2);
t405 = pkin(1) * t476 - t455;
t372 = qJD(1) * t405;
t846 = qJD(1) * t318 - t372;
t808 = pkin(7) * t476;
t629 = -t405 - t808;
t625 = -rSges(3,2) * t479 + t476 * rSges(3,3);
t844 = t409 + t625;
t796 = rSges(6,2) * t449;
t359 = rSges(6,3) * t450 + t796;
t467 = t479 * rSges(6,1);
t258 = -t359 * t476 + t467;
t365 = t479 * pkin(5) + pkin(8) * t757;
t257 = rSges(6,1) * t476 + t359 * t479;
t594 = Icges(7,1) * t474 + t772;
t842 = t449 * t594 + t768;
t605 = rSges(7,1) * t474 + rSges(7,2) * t477;
t234 = rSges(7,3) * t450 + t449 * t605;
t828 = t543 * (-Icges(7,2) * t322 - t164 - t294) - t342 * (-Icges(7,2) * t324 + t165 + t296) - t418 * (-Icges(7,2) * t758 + t232 + t412);
t827 = m(6) / 0.2e1;
t826 = m(7) / 0.2e1;
t825 = -pkin(1) - pkin(5);
t680 = qJD(3) * qJD(6);
t638 = t450 * t680;
t214 = qJD(1) * t543 + t476 * t638;
t824 = t214 / 0.2e1;
t215 = qJD(1) * t342 - t479 * t638;
t823 = t215 / 0.2e1;
t822 = t543 / 0.2e1;
t821 = -t543 / 0.2e1;
t820 = -t342 / 0.2e1;
t819 = t342 / 0.2e1;
t818 = -t418 / 0.2e1;
t817 = t418 / 0.2e1;
t816 = t450 / 0.2e1;
t815 = t476 / 0.2e1;
t813 = -rSges(6,1) - pkin(1);
t812 = rSges(3,2) - pkin(1);
t806 = pkin(7) * qJD(1) ^ 2;
t805 = pkin(8) * t450;
t623 = -qJD(1) * t450 - qJD(6);
t495 = t476 * t623 - t643;
t138 = -t418 * t752 + t477 * t495;
t139 = t418 * t747 + t474 * t495;
t530 = -t648 + t656;
t74 = Icges(7,5) * t139 + Icges(7,6) * t138 + Icges(7,3) * t530;
t76 = Icges(7,4) * t139 + Icges(7,2) * t138 + Icges(7,6) * t530;
t78 = Icges(7,1) * t139 + Icges(7,4) * t138 + Icges(7,5) * t530;
t8 = (qJD(3) * t582 + t74) * t450 + (qJD(3) * t158 + t474 * t78 + t477 * t76 + (-t160 * t474 - t164 * t477) * qJD(6)) * t449;
t804 = t8 * t543;
t560 = t623 * t479;
t136 = t477 * t560 + (t418 * t474 + t477 * t696) * t476;
t642 = t449 * t693;
t137 = -t418 * t749 + (t560 + t642) * t474;
t73 = Icges(7,5) * t137 + Icges(7,6) * t136 + Icges(7,3) * t531;
t75 = Icges(7,4) * t137 + Icges(7,2) * t136 + Icges(7,6) * t531;
t77 = Icges(7,1) * t137 + Icges(7,4) * t136 + Icges(7,5) * t531;
t9 = (qJD(3) * t581 + t73) * t450 + (-qJD(3) * t159 + t474 * t77 + t477 * t75 + (-t162 * t474 + t165 * t477) * qJD(6)) * t449;
t803 = t9 * t342;
t800 = -pkin(1) + t473;
t795 = rSges(3,3) * t479;
t235 = -rSges(7,3) * t449 + t450 * t605;
t297 = (rSges(7,1) * t477 - rSges(7,2) * t474) * t449;
t140 = qJD(3) * t235 + qJD(6) * t297;
t388 = pkin(8) * t648;
t226 = qJD(1) * t365 - t388;
t480 = qJD(3) ^ 2;
t682 = qJD(1) * qJD(2);
t445 = t479 * t682;
t617 = -t479 * t806 + t445;
t622 = t810 * t681;
t563 = t479 * t622 + t617;
t640 = t362 * t681;
t509 = qJD(5) * t642 + t254 * t693 + t479 * t640 + t563;
t556 = -t131 - t691 + t847;
t606 = rSges(7,1) * t139 + rSges(7,2) * t138;
t80 = rSges(7,3) * t530 + t606;
t10 = -t480 * t440 - t543 * t140 + t215 * t234 - t418 * t80 + (-pkin(8) * t476 * t480 + t166 * t680) * t449 + (-t226 + (pkin(8) * qJD(3) - qJD(5)) * t754 + t556) * qJD(1) + t509;
t793 = t10 * t476;
t168 = t324 * rSges(7,1) + t323 * rSges(7,2) + rSges(7,3) * t757;
t717 = pkin(8) * t655 + t385;
t225 = -pkin(5) * t698 + t717;
t421 = pkin(8) * t756;
t664 = pkin(4) * t531 + qJ(5) * t642;
t132 = (-qJ(5) * t697 - t690) * t450 + t664;
t452 = qJD(2) * t476;
t710 = qJ(2) * t697 + t452;
t723 = qJD(1) * (-pkin(1) * t698 + t710) + t476 * t682;
t564 = -t476 * t806 + t723;
t678 = t480 * t811;
t497 = t476 * t622 + t479 * t678 + t564 + (t196 + t451) * qJD(1);
t493 = qJD(1) * t132 + t476 * t640 + t497;
t644 = t450 * t690;
t79 = t137 * rSges(7,1) + t136 * rSges(7,2) + rSges(7,3) * t531;
t11 = (pkin(8) * t654 - t254 * t479 + (-qJD(6) * t168 - t689) * t449) * qJD(3) + (t225 - t644) * qJD(1) - t214 * t234 + t418 * t79 - t140 * t342 + t480 * t421 + t493;
t792 = t11 * t479;
t379 = rSges(6,3) * t642;
t155 = -rSges(6,2) * t649 - qJD(1) * t257 + t379;
t333 = t359 * qJD(3);
t732 = -t254 - t333;
t28 = (t155 - t644) * qJD(1) + (-t604 * t698 + (-t438 + t732) * t479) * qJD(3) + t493;
t791 = t28 * t479;
t156 = qJD(1) * t258 + qJD(3) * t304;
t29 = (qJD(3) * t333 - t678) * t476 + (-t156 + t556 + t852) * qJD(1) + t509;
t790 = t29 * t476;
t555 = -t362 * t692 - t453 + t624;
t420 = pkin(4) * t757;
t298 = -qJ(5) * t755 + t420;
t319 = -t479 * t799 + t440;
t628 = t409 + t807;
t615 = t319 + t628;
t562 = t298 + t615;
t43 = t645 + t168 * t418 - t234 * t342 - t388 + (t365 + t562) * qJD(1) + t555;
t787 = t43 * t234;
t786 = t43 * t479;
t466 = t479 * rSges(4,3);
t465 = t479 * rSges(5,3);
t364 = rSges(5,1) * t450 - rSges(5,2) * t449;
t314 = t364 * t693;
t334 = t607 * qJD(3);
t54 = t334 * t692 + (t154 + t314) * qJD(1) + t497;
t785 = t54 * t479;
t784 = t60 * t215;
t61 = t159 * t450 + t449 * t581;
t783 = t61 * t214;
t661 = t697 * t912 + (t646 + t652) * rSges(4,1);
t694 = qJD(3) * t475;
t184 = (-rSges(4,2) * t694 - rSges(4,3) * qJD(1)) * t476 + t661;
t408 = rSges(4,1) * t478 - rSges(4,2) * t475;
t357 = t408 * t693;
t608 = rSges(4,1) * t475 + t912;
t370 = t608 * qJD(3);
t83 = t370 * t692 + (t184 + t357) * qJD(1) + t564;
t782 = t83 * t479;
t344 = t408 * t479;
t183 = -qJD(3) * t344 + (t476 * t608 + t466) * qJD(1);
t650 = t408 * t692;
t84 = -t370 * t693 + (-t183 - t332 + t650) * qJD(1) + t617;
t781 = t84 * t476;
t463 = t476 * rSges(4,3);
t309 = t479 * t608 - t463;
t122 = t357 + t452 + (t309 + t629) * qJD(1);
t761 = t122 * t479;
t301 = t364 * t476;
t745 = -t132 - t196;
t255 = rSges(5,1) * t757 + rSges(5,2) * t755 + t465;
t731 = -t255 - t319;
t266 = t479 * t318;
t411 = qJ(5) * t754;
t302 = pkin(4) * t756 - t411;
t730 = -t479 * t302 - t266;
t725 = -t298 - t319;
t724 = t302 + t318;
t441 = pkin(3) * t748;
t722 = t750 + t441;
t714 = -t411 + t455;
t711 = t428 + t453;
t709 = rSges(3,2) * t698 + rSges(3,3) * t697;
t708 = t452 - t372;
t695 = qJD(3) * t450;
t679 = -rSges(4,3) - pkin(1) - pkin(7);
t676 = pkin(3) * t694;
t185 = t197 * t692;
t673 = qJD(3) * t439 + t131 * t692 + t185;
t672 = -t155 + t745;
t671 = -t225 + t745;
t669 = -t257 + t724;
t668 = -t258 + t725;
t667 = -t365 + t725;
t361 = pkin(5) * t476 - t421;
t666 = -t361 + t724;
t660 = t452 + t713;
t659 = t434 + t711;
t308 = rSges(4,1) * t751 + rSges(4,2) * t748 + t466;
t262 = t318 * t692;
t658 = -t302 * t692 - t262 + t438;
t657 = t440 + t409;
t647 = t475 * t693;
t641 = t757 / 0.2e1;
t639 = t449 * t680;
t636 = t697 / 0.2e1;
t635 = -t696 / 0.2e1;
t634 = -t693 / 0.2e1;
t633 = t693 / 0.2e1;
t632 = -t692 / 0.2e1;
t630 = -t362 - t810;
t627 = t234 + t805;
t621 = -t168 + t667;
t620 = t314 + t660;
t619 = t420 + t657;
t616 = t318 + t629;
t611 = qJD(6) * t635;
t610 = t809 + t811;
t603 = t10 * t479 + t11 * t476;
t47 = -t159 * t756 - t849;
t602 = t46 * t479 - t47 * t476;
t601 = t46 * t476 + t47 * t479;
t600 = t476 * t49 - t479 * t48;
t599 = t476 * t48 + t479 * t49;
t598 = t476 * t61 - t479 * t60;
t597 = t476 * t60 + t479 * t61;
t123 = -t650 - t453 + (t308 + t628) * qJD(1);
t583 = t122 * t476 - t123 * t479;
t580 = t166 * t476 + t168 * t479;
t561 = t302 + t616;
t283 = (Icges(7,5) * t477 - Icges(7,6) * t474) * t449;
t133 = qJD(3) * t229 + qJD(6) * t283;
t231 = -Icges(7,6) * t449 + t450 * t590;
t134 = (-Icges(7,2) * t474 + t772) * t688 + t231 * qJD(3);
t233 = -Icges(7,5) * t449 + t450 * t594;
t291 = (Icges(7,1) * t477 - t773) * t449;
t135 = qJD(3) * t233 + qJD(6) * t291;
t19 = (qJD(3) * t579 + t133) * t450 + (-qJD(3) * t228 + t134 * t477 + t135 * t474 + (-t230 * t474 + t232 * t477) * qJD(6)) * t449;
t89 = t228 * t450 + t449 * t579;
t558 = t19 * t418 - t639 * t89;
t554 = t659 + t663;
t553 = t618 + t710;
t546 = -t250 - t574;
t305 = t364 * t479;
t544 = -t364 * t692 + t691;
t462 = t476 * rSges(5,3);
t256 = -t462 + t886;
t536 = -t256 * t479 + t476 * t731;
t169 = (-t308 * t476 - t309 * t479) * qJD(3);
t529 = -t266 + t536;
t528 = -t158 * t543 - t159 * t342 - t228 * t418;
t527 = -(-Icges(7,5) * t321 - Icges(7,6) * t322) * t543 + (Icges(7,5) * t323 - Icges(7,6) * t324) * t342 + t283 * t418;
t526 = t362 * t693 - t644 + t660;
t498 = t553 + t664;
t496 = -t359 + t610;
t494 = qJD(1) * t302 + t526 + t846;
t492 = (Icges(7,1) * t323 - t162 - t774) * t342 - (-Icges(7,1) * t321 - t160 - t295) * t543 + (-t230 + t291) * t418;
t34 = t166 * t342 + t168 * t543 + (t361 * t479 + t476 * t667) * qJD(3) + t658;
t42 = -t234 * t543 + (-t361 + t561) * qJD(1) + t526 + t863;
t483 = t34 * t580 + (-t42 * t479 - t43 * t476) * t234;
t406 = rSges(3,2) * t476 + t795;
t343 = t408 * t476;
t313 = t604 * t693;
t270 = t707 * t695;
t269 = t642 - t653;
t217 = qJD(1) * t844 - t453;
t216 = t452 + (-t405 + t406) * qJD(1);
t205 = t234 * t479;
t204 = t234 * t476;
t203 = t842 * t479;
t202 = t842 * t476;
t201 = t230 * t479;
t200 = t230 * t476;
t193 = rSges(7,1) * t323 - rSges(7,2) * t324;
t192 = -rSges(7,1) * t321 - rSges(7,2) * t322;
t182 = t445 + (-qJD(1) * t625 - t332) * qJD(1);
t181 = qJD(1) * t709 + t723;
t153 = -qJD(3) * t305 + (t476 * t607 + t465) * qJD(1);
t91 = (t255 + t615) * qJD(1) + t544 - t711;
t90 = (t256 + t616) * qJD(1) + t620;
t88 = qJD(3) * t536 - t262;
t67 = -t852 + (t258 + t562) * qJD(1) + t555;
t66 = -t313 + (-t257 + t561) * qJD(1) + t526;
t62 = (t257 * t479 + t476 * t668) * qJD(3) + t658;
t55 = (-qJD(3) * t334 - t678) * t476 + (-t153 - t544 + t847) * qJD(1) + t563;
t17 = -t133 * t756 - t134 * t321 + t135 * t322 + t138 * t230 + t139 * t232 + t228 * t530;
t16 = t133 * t757 + t134 * t323 + t135 * t324 + t136 * t230 + t137 * t232 + t228 * t531;
t15 = (t156 * t479 + t672 * t476 + (t476 * t669 + t479 * t668) * qJD(1)) * qJD(3) + t673;
t14 = t342 * t61 + t418 * t89 - t543 * t60;
t12 = t342 * t47 - t862;
t7 = t138 * t162 + t139 * t165 + t159 * t530 - t321 * t75 + t322 * t77 - t73 * t756;
t6 = t138 * t160 - t139 * t164 - t158 * t530 - t321 * t76 + t322 * t78 - t74 * t756;
t5 = t136 * t162 + t137 * t165 + t159 * t531 + t323 * t75 + t324 * t77 + t73 * t757;
t4 = t136 * t160 - t137 * t164 - t158 * t531 + t323 * t76 + t324 * t78 + t74 * t757;
t3 = t166 * t214 - t168 * t215 + t543 * t79 + t342 * t80 + (t226 * t479 + t671 * t476 + (t476 * t666 + t479 * t667) * qJD(1)) * qJD(3) + t673;
t2 = t17 * t418 + t214 * t47 + t215 * t46 + t342 * t7 + t499 * t639 - t543 * t6;
t1 = t16 * t418 + t214 * t49 + t215 * t48 + t342 * t5 - t4 * t543 - t639 * t82;
t18 = [(-(-t313 - t66 + (-t257 - t808) * qJD(1) + t494) * t67 + t29 * t714 + t66 * t554 + t28 * (t467 + t619) + t67 * (t379 + t498) + (t29 * t496 + t66 * (-rSges(6,2) * t695 + rSges(6,3) * t696 - t439) - t28 * t473 + (t66 * t813 + t67 * (-t359 - t762)) * qJD(1)) * t479 + (t29 * (-rSges(6,1) + t800) - t66 * qJD(4) - t28 * t796 + (t28 * (-rSges(6,3) - qJ(5)) + t67 * (-rSges(6,2) * qJD(3) - qJD(5))) * t450 + (t66 * (-qJ(2) - t496) + t67 * t813) * qJD(1)) * t476) * m(6) - t499 * t823 - t804 / 0.2e1 + (t84 * (-t463 + t629) + t122 * t453 + t83 * (t308 + t409) + t123 * (-rSges(4,2) * t647 + t661 + t710) + (qJD(3) * t122 * t408 + t83 * pkin(7) + t608 * t84) * t479 + (t679 * t761 + (t122 * (-qJ(2) - t608) + t123 * t679) * t476) * qJD(1) - (-t122 + t357 + (t309 - t808) * qJD(1) + t708) * t123) * m(4) + t558 + t803 / 0.2e1 + (t17 + t13) * t821 + ((-t473 * t479 + t255 + t657) * t54 + (t476 * t800 + t479 * t532 + t455 - t462) * t55 + (t659 + (rSges(5,1) * t695 - pkin(1) * qJD(1) + t552) * t479 + (-qJD(4) + (-qJ(2) - t532) * qJD(1)) * t476) * t90 + (t553 + t662 + (-t674 + (-rSges(5,3) - pkin(1)) * qJD(1)) * t476 + t90 - (t256 - t808) * qJD(1) - t620 - t846) * t91) * m(5) + (-qJD(3) * t921 + t368 * t475 - t369 * t478 + t449 * t927 + t450 * t926) * qJD(1) + t783 / 0.2e1 + t784 / 0.2e1 + t13 * t822 + (t182 * (t476 * t812 + t455 + t795) + t216 * t453 + t181 * t844 + t217 * (t709 + t710) + (t216 * t812 * t479 + (t216 * (-rSges(3,3) - qJ(2)) - t217 * pkin(1)) * t476) * qJD(1) - (qJD(1) * t406 - t216 + t708) * t217) * m(3) + t16 * t819 + t82 * t824 + ((t48 + (t158 * t476 + t159 * t479) * t449 + t849 + t890) * t342 + t12 + t862) * t820 + (t10 * (-t166 + t421 + t714) + t42 * (t388 + t554 - t606) + t11 * (t619 + t168 + t365) + t43 * (t498 + t79 + t717) + (t10 * t610 - t11 * t473 + t42 * (rSges(7,3) * qJD(3) - qJD(5)) * t450 + (t42 * t825 - t43 * t762) * qJD(1)) * t479 + (t10 * (-pkin(5) + t800) - t42 * qJD(4) + (-qJ(5) * t11 - qJD(5) * t43) * t450 + (t43 * t825 + (-qJ(2) - t811 + (-rSges(7,3) - pkin(4) - pkin(8)) * t449) * t42) * qJD(1)) * t476 - (-t42 + (-t361 - t808) * qJD(1) + t494 + t863) * t43 + t787 * t543) * m(7) + (t871 + t874 + t876) * t633 + (((t224 + t899 - t903 + t920) * t479 + ((t546 - t883 + t884) * t479 + t870 - t885 + t902) * t476) * qJD(3) + t898) * t634 + (((t546 * t476 + t224 - t99) * t476 + t884 * t471 + (-t97 - t223 + (t883 + t884) * t479 + t870 + t885) * t479) * qJD(3) + t875 + t893) * t632 + ((t114 + t869) * qJD(1) + t872 + t873) * t692 / 0.2e1 - ((-t867 + t914) * t476 + t915 * t479) * t681 / 0.2e1; 0.2e1 * (t793 / 0.2e1 - t792 / 0.2e1) * m(7) + 0.2e1 * (-t791 / 0.2e1 + t790 / 0.2e1) * m(6) + 0.2e1 * (t55 * t815 - t785 / 0.2e1) * m(5) + 0.2e1 * (t781 / 0.2e1 - t782 / 0.2e1) * m(4) + 0.2e1 * (-t181 * t479 / 0.2e1 + t182 * t815) * m(3); t597 * t611 + ((-t200 * t321 + t202 * t322) * t342 - (t201 * t321 - t203 * t322) * t543 + (-t231 * t321 + t233 * t322) * t418 + (t449 * t499 + t47 * t755) * qJD(6) + ((-qJD(6) * t46 + t528) * t450 + t901) * t479) * t822 + (((t200 * t477 + t202 * t474 - t159) * t342 - (-t201 * t477 - t203 * t474 + t158) * t543 + (t231 * t477 + t233 * t474 - t228) * t418 - t89 * qJD(6)) * t449 + (qJD(6) * t598 - t829) * t450) * t818 + ((t200 * t323 + t202 * t324) * t342 - (-t201 * t323 - t203 * t324) * t543 + (t231 * t323 + t233 * t324) * t418 + (-t449 * t82 - t48 * t754) * qJD(6) + ((qJD(6) * t49 - t528) * t450 - t901) * t476) * t820 + t12 * t686 * t816 + t14 * t688 / 0.2e1 + (-qJD(1) * t598 + t476 * t8 + t479 * t9) * t817 + (-qJD(1) * t600 + t4 * t476 + t479 * t5) * t819 + (qJD(1) * t602 + t476 * t6 + t479 * t7) * t821 + t601 * t823 + t599 * t824 - t450 * t13 * t687 / 0.2e1 + (-t483 * t877 + t786 * t878 + t10 * t722 + t3 * t730 + (t11 * (-t627 + t630) + t3 * (t166 + t361)) * t479 + (t787 * qJD(1) + t10 * t627 + t3 * t621) * t476 + (t234 * t697 - t166 * t688 - t205 * t418 + t235 * t543 + (-(-pkin(8) * t449 - t811) * qJD(3) - pkin(8) * t696 + t140 - t676) * t476 + t880) * t42 + (-t204 * t418 + t235 * t342 + t168 * t688 + (-t140 - t254) * t479 + t881) * t43 + (-t204 * t543 + t205 * t342 - (-t476 * t750 - t707 * t805 - t537) * qJD(3) + (qJD(1) * t621 + t226 + t80) * t479 + (-t79 + t671 + (-t166 + t666) * qJD(1)) * t476 + t879) * t34) * m(7) + (t29 * t722 + t15 * t730 + (t28 * (t604 + t630) + t15 * t257) * t479 + (t15 * t668 - t29 * t604) * t476 + ((-(-t358 - t359) * qJD(3) + t732) * t479 + t881) * t67 + ((-(t359 - t811) * qJD(3) + t333 - t676) * t476 + t880) * t66 + (-(-t537 + t304 * t479 + (t604 * t476 - t750) * t476) * qJD(3) + (qJD(1) * t668 + t156) * t479 + (qJD(1) * t669 + t672) * t476 + t879) * t62) * m(6) + (t55 * (t441 + t301) - t90 * pkin(3) * t647 + (-t364 - t810) * t785 + (t153 * t692 + t693 * t848 + t185) * t529 + t88 * (t153 * t479 + t476 * t848 + t194) - (t91 * t886 + t88 * (-t305 * t479 - t537) + (-t88 * t301 - t532 * t90) * t476) * qJD(3) + (-t476 * t90 + t479 * t91) * t334 + ((t476 * t91 + t479 * t90) * t364 + (qJD(3) * t529 + t88) * (t731 * t479 + (t256 + t318) * t476) - t90 * t305 - t91 * t301) * qJD(1)) * m(5) + (-(t122 * t344 + t123 * t343) * qJD(1) - (t169 * (-t343 * t476 - t344 * t479) - t583 * t608) * qJD(3) + 0.2e1 * t169 * (t183 * t479 - t184 * t476 + (-t308 * t479 + t309 * t476) * qJD(1)) - t583 * t370 + (t781 - t782 + (t123 * t476 + t761) * qJD(1)) * t408) * m(4) - ((t449 * t911 + t450 * t910 + t475 * t510 + t478 * t511) * qJD(3) + (-t908 * t449 + t909 * t450 - t475 * t716 - t478 * t715) * qJD(1)) * qJD(1) / 0.2e1 + (t872 * t479 + t871 * t476 + (t867 * t476 + t869 * t479) * qJD(1)) * qJD(1) / 0.2e1 + ((t693 * t865 + t866) * t476 + ((t476 * t868 + t861) * qJD(3) - t900) * t479) * t634 + ((t692 * t868 + t866) * t479 + ((t479 * t865 - t861) * qJD(3) + t900) * t476) * t632 + (t874 * qJD(1) + t2 + ((t899 * qJD(1) + t892 * t479) * t479 + (t896 * t476 - t930 * qJD(1) + (-t891 + t897) * t479) * t476) * t889) * t815 + (t873 * qJD(1) + t1 + ((t870 * qJD(1) + t897 * t479) * t479 + (t891 * t476 - t929 * qJD(1) + (-t892 + t896) * t479) * t476) * t889) * t479 / 0.2e1 - (t13 + t876 + t894) * t698 / 0.2e1 + (t12 + t875 + t895) * t636; m(5) * (t476 * t54 + t479 * t55) + m(6) * (t28 * t476 + t29 * t479) + m(7) * t603; -m(6) * (t269 * t66 + t270 * t62 - t67 * t882) - m(7) * (t269 * t42 + t270 * t34 - t43 * t882) + 0.2e1 * ((t66 * t693 - t67 * t692 + t15) * t827 + (t42 * t693 - t43 * t692 + t3) * t826) * t449 + 0.2e1 * ((qJD(3) * t62 - t66 * t697 - t67 * t698 - t790 + t791) * t827 + (qJD(3) * t34 - t42 * t697 - t43 * t698 + t792 - t793) * t826) * t450; t1 * t641 + (t449 * t600 + t450 * t82) * t824 + ((qJD(3) * t600 + t16) * t450 + (qJD(1) * t599 - qJD(3) * t82 - t4 * t479 + t476 * t5) * t449) * t819 - t2 * t756 / 0.2e1 + (-t449 * t602 - t450 * t499) * t823 + ((-qJD(3) * t602 + t17) * t450 + (qJD(1) * t601 + qJD(3) * t499 + t476 * t7 - t479 * t6) * t449) * t821 + t14 * t635 + (t558 + t783 + t784 + t803 - t804) * t816 + (t449 * t598 + t450 * t89) * t611 + ((qJD(3) * t598 + t19) * t450 + (qJD(1) * t597 - qJD(3) * t89 + t476 * t9 - t479 * t8) * t449) * t817 + (-t323 * t828 + t324 * t492 + t527 * t757) * t820 + (t321 * t828 + t492 * t322 - t527 * t756) * t822 + (t527 * t450 + (t474 * t492 - t477 * t828) * t449) * t818 + (t449 * t636 + t450 * t633) * t13 + (qJD(1) * t641 + t450 * t632) * t12 + ((qJD(3) * t483 - t10 * t166 + t11 * t168 - t42 * t80 + t43 * t79) * t450 + (t42 * (qJD(3) * t166 - t140 * t479) + t43 * (-qJD(3) * t168 - t140 * t476) + t3 * t580 + t34 * (t166 * t697 - t168 * t698 + t476 * t80 + t479 * t79) + ((t42 * t476 - t786) * qJD(1) - t603) * t234) * t449 - t42 * (-t192 * t418 - t297 * t543) - t43 * (t193 * t418 - t297 * t342) - t34 * (t192 * t342 + t193 * t543)) * m(7);];
tauc  = t18(:);

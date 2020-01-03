% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR7
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR7_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:57
% EndTime: 2019-12-31 19:36:10
% DurationCPUTime: 57.70s
% Computational Cost: add. (18356->1014), mult. (27693->1310), div. (0->0), fcn. (24370->8), ass. (0->517)
t934 = Icges(5,4) - Icges(4,5);
t933 = Icges(5,5) - Icges(4,6);
t932 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t455 = qJ(2) + pkin(8);
t435 = sin(t455);
t436 = cos(t455);
t460 = sin(qJ(2));
t463 = cos(qJ(2));
t908 = Icges(3,5) * t463 - Icges(3,6) * t460 + t933 * t435 - t436 * t934;
t785 = Icges(4,4) * t435;
t346 = Icges(4,2) * t436 + t785;
t771 = Icges(5,6) * t435;
t565 = Icges(5,3) * t436 + t771;
t912 = -t346 - t565;
t425 = Icges(4,4) * t436;
t348 = Icges(4,1) * t435 + t425;
t770 = Icges(5,6) * t436;
t567 = Icges(5,2) * t435 + t770;
t931 = t348 + t567;
t464 = cos(qJ(1));
t930 = t932 * t464;
t929 = Icges(3,5) * t460 + Icges(4,5) * t435 + Icges(3,6) * t463 + Icges(4,6) * t436;
t349 = Icges(4,1) * t436 - t785;
t568 = Icges(5,2) * t436 - t771;
t928 = t349 + t568;
t566 = -Icges(5,3) * t435 + t770;
t573 = -Icges(4,2) * t435 + t425;
t927 = t566 + t573;
t461 = sin(qJ(1));
t738 = t461 * t463;
t741 = t460 * t461;
t745 = t436 * t461;
t747 = t435 * t461;
t896 = -Icges(3,5) * t738 + Icges(3,6) * t741 + t745 * t934 - t933 * t747 + t930;
t909 = t932 * t461 + t908 * t464;
t571 = Icges(5,4) * t435 + Icges(5,5) * t436;
t915 = -t571 + t929;
t772 = Icges(4,6) * t464;
t244 = Icges(4,4) * t745 - Icges(4,2) * t747 - t772;
t395 = Icges(5,6) * t747;
t783 = Icges(5,4) * t464;
t251 = Icges(5,2) * t745 - t395 + t783;
t773 = Icges(3,6) * t464;
t279 = Icges(3,4) * t738 - Icges(3,2) * t741 - t773;
t926 = t244 * t435 - t251 * t436 + t279 * t460;
t777 = Icges(5,5) * t464;
t249 = Icges(5,6) * t745 - Icges(5,3) * t747 + t777;
t925 = t244 + t249;
t245 = Icges(4,6) * t461 + t464 * t573;
t248 = Icges(5,5) * t461 - t464 * t566;
t924 = t245 - t248;
t400 = Icges(4,4) * t747;
t778 = Icges(4,5) * t464;
t246 = Icges(4,1) * t745 - t400 - t778;
t923 = t246 + t251;
t247 = Icges(4,5) * t461 + t349 * t464;
t746 = t435 * t464;
t396 = Icges(5,6) * t746;
t744 = t436 * t464;
t784 = Icges(5,4) * t461;
t250 = -Icges(5,2) * t744 + t396 + t784;
t922 = t247 - t250;
t920 = t912 * qJD(2);
t919 = t931 * qJD(2);
t786 = Icges(3,4) * t460;
t380 = Icges(3,2) * t463 + t786;
t445 = Icges(3,4) * t463;
t382 = Icges(3,1) * t460 + t445;
t918 = t346 * t435 - t348 * t436 + t380 * t460 - t382 * t463;
t383 = Icges(3,1) * t463 - t786;
t282 = Icges(3,5) * t461 + t383 * t464;
t916 = -t247 * t745 - t248 * t747 - t282 * t738;
t907 = t435 * t565 - t436 * t567 + t918;
t421 = Icges(3,4) * t741;
t779 = Icges(3,5) * t464;
t281 = Icges(3,1) * t738 - t421 - t779;
t892 = -t246 * t436 + t249 * t435 - t281 * t463 + t926;
t914 = t927 * qJD(2);
t913 = t928 * qJD(2);
t910 = t929 * t464;
t906 = t464 * t909 + t916;
t736 = t463 * t464;
t852 = t247 * t744 + t248 * t746 + t282 * t736 + t461 * t909;
t905 = -t246 * t744 + t249 * t746 - t281 * t736 + t461 * t896;
t850 = t915 * t461;
t574 = -Icges(3,2) * t460 + t445;
t280 = Icges(3,6) * t461 + t464 * t574;
t904 = t245 * t435 + t250 * t436 + t280 * t460;
t903 = t918 * t461 + t910;
t902 = t920 * t464 + (-t461 * t927 + t772 - t777) * qJD(1);
t901 = t924 * qJD(1) + t920 * t461;
t900 = -t919 * t464 + (-t461 * t928 + t778 - t783) * qJD(1);
t899 = t919 * t461 + (-t464 * t568 - t247 + t784) * qJD(1);
t898 = -t892 * t461 + t896 * t464;
t883 = -t245 * t747 - t250 * t745 - t280 * t741 - t906;
t740 = t460 * t464;
t882 = -t244 * t746 + t251 * t744 - t279 * t740 - t905;
t881 = -t245 * t746 - t250 * t744 - t280 * t740 + t852;
t897 = -t907 * t464 + t850;
t858 = t279 * t463 + t281 * t460 + t923 * t435 + t925 * t436;
t857 = t280 * t463 + t282 * t460 + t922 * t435 + t924 * t436;
t894 = t915 * qJD(2);
t893 = t247 * t436 + t248 * t435 + t282 * t463 - t904;
t891 = -t924 * t461 + t925 * t464;
t890 = t931 + t927;
t889 = t912 + t928;
t888 = (t395 + t400 + (Icges(4,2) + Icges(5,3)) * t745 - t923) * t464 + (-Icges(5,3) * t744 - t346 * t464 - t396 + t922) * t461;
t360 = t574 * qJD(2);
t361 = t383 * qJD(2);
t887 = -t360 * t460 + t361 * t463 + t913 * t436 - t914 * t435 + (-t380 * t463 - t382 * t460 - t435 * t931 + t912 * t436) * qJD(2) + t915 * qJD(1);
t886 = t909 * qJD(1);
t885 = t907 * qJD(1) + t908 * qJD(2);
t884 = t903 * qJD(1);
t880 = t897 * qJD(1);
t522 = qJD(2) * t380;
t183 = -t464 * t522 + (-t461 * t574 + t773) * qJD(1);
t524 = qJD(2) * t382;
t185 = -t464 * t524 + (-t383 * t461 + t779) * qJD(1);
t879 = -qJD(2) * t857 - t183 * t460 + t185 * t463 - t435 * t902 + t436 * t900 + t886;
t184 = qJD(1) * t280 - t461 * t522;
t186 = qJD(1) * t282 - t461 * t524;
t878 = qJD(1) * t896 + qJD(2) * t858 + t184 * t460 - t186 * t463 + t435 * t901 + t436 * t899;
t877 = (t461 * t881 - t464 * t882) * qJD(2);
t876 = (t461 * t883 - t464 * t898) * qJD(2);
t875 = qJD(1) * t892 - t461 * t894 + t886;
t874 = -t894 * t464 + (-t461 * t908 - t893 + t930) * qJD(1);
t671 = qJD(5) * t436;
t677 = qJD(2) * t461;
t334 = t464 * t671 + t677;
t676 = qJD(2) * t464;
t335 = -t461 * t671 + t676;
t672 = qJD(5) * t435;
t408 = qJD(1) + t672;
t462 = cos(qJ(5));
t737 = t462 * t464;
t459 = sin(qJ(5));
t743 = t459 * t461;
t316 = t435 * t737 - t743;
t739 = t461 * t462;
t742 = t459 * t464;
t317 = t435 * t742 + t739;
t154 = Icges(6,5) * t317 + Icges(6,6) * t316 + Icges(6,3) * t744;
t782 = Icges(6,4) * t317;
t157 = Icges(6,2) * t316 + Icges(6,6) * t744 + t782;
t298 = Icges(6,4) * t316;
t160 = Icges(6,1) * t317 + Icges(6,5) * t744 + t298;
t45 = t154 * t744 + t316 * t157 + t317 * t160;
t318 = t435 * t739 + t742;
t319 = -t435 * t743 + t737;
t156 = -Icges(6,5) * t319 + Icges(6,6) * t318 + Icges(6,3) * t745;
t300 = Icges(6,4) * t319;
t159 = Icges(6,2) * t318 + Icges(6,6) * t745 - t300;
t299 = Icges(6,4) * t318;
t161 = Icges(6,1) * t319 - Icges(6,5) * t745 - t299;
t46 = t156 * t744 + t316 * t159 - t161 * t317;
t569 = Icges(6,5) * t459 + Icges(6,6) * t462;
t494 = -Icges(6,3) * t435 + t436 * t569;
t781 = Icges(6,4) * t459;
t570 = Icges(6,2) * t462 + t781;
t495 = -Icges(6,6) * t435 + t436 * t570;
t780 = Icges(6,4) * t462;
t575 = Icges(6,1) * t459 + t780;
t496 = -Icges(6,5) * t435 + t436 * t575;
t79 = -t316 * t495 - t317 * t496 - t494 * t744;
t12 = t334 * t45 - t335 * t46 + t79 * t408;
t47 = t154 * t745 + t318 * t157 - t319 * t160;
t48 = t156 * t745 + t159 * t318 + t161 * t319;
t80 = -t318 * t495 + t319 * t496 - t494 * t745;
t13 = t334 * t47 - t335 * t48 + t408 * t80;
t563 = t159 * t462 - t161 * t459;
t61 = t156 * t435 - t436 * t563;
t765 = qJ(4) * t436;
t350 = pkin(3) * t435 - t765;
t301 = t350 * t461;
t680 = qJD(1) * t461;
t643 = t460 * t680;
t417 = pkin(2) * t643;
t673 = qJD(4) * t464;
t873 = -qJD(1) * t301 - t436 * t673 - t417;
t872 = 0.2e1 * qJD(2);
t262 = rSges(4,1) * t745 - rSges(4,2) * t747 - t464 * rSges(4,3);
t446 = t461 * rSges(4,3);
t263 = rSges(4,1) * t744 - rSges(4,2) * t746 + t446;
t554 = t262 * t461 + t263 * t464;
t453 = t464 * pkin(6);
t406 = pkin(1) * t461 - t453;
t458 = -qJ(3) - pkin(6);
t428 = t464 * t458;
t816 = pkin(2) * t463;
t429 = pkin(1) + t816;
t691 = -t461 * t429 - t428;
t240 = t406 + t691;
t451 = t461 * pkin(6);
t407 = t464 * pkin(1) + t451;
t413 = t464 * t429;
t606 = -t458 * t461 + t413;
t241 = t606 - t407;
t721 = -t461 * t240 + t464 * t241;
t511 = t554 + t721;
t722 = -t240 * t677 + t241 * t676;
t81 = qJD(2) * t554 + t722;
t871 = qJD(2) * t511 + t81;
t292 = t571 * t464;
t793 = t435 * rSges(5,2);
t585 = rSges(5,3) * t436 + t793;
t302 = t585 * t461;
t817 = pkin(2) * t460;
t617 = -t350 - t817;
t598 = t585 + t617;
t540 = t598 * t464;
t647 = t435 * t680;
t636 = t436 * t676;
t369 = qJ(4) * t636;
t387 = t435 * t673;
t698 = t369 + t387;
t638 = t435 * t676;
t645 = t436 * t680;
t866 = t638 + t645;
t128 = -pkin(3) * t866 - qJ(4) * t647 + t698;
t637 = t436 * t677;
t679 = qJD(1) * t464;
t646 = t435 * t679;
t275 = t637 + t646;
t639 = t435 * t677;
t376 = pkin(3) * t639;
t674 = qJD(4) * t461;
t631 = t435 * t674;
t644 = t436 * t679;
t129 = pkin(3) * t644 + qJ(4) * t275 - t376 + t631;
t354 = pkin(3) * t436 + qJ(4) * t435;
t304 = t354 * t461;
t306 = t350 * t464;
t424 = qJD(4) * t435;
t433 = pkin(6) * t679;
t437 = qJD(3) * t461;
t635 = t460 * t676;
t544 = -pkin(2) * t635 + t437;
t810 = pkin(1) - t429;
t179 = -t433 + (t461 * t810 - t428) * qJD(1) + t544;
t416 = t677 * t817;
t689 = qJD(3) * t464 + t416;
t648 = t458 * t680 + t689;
t180 = (-t464 * t810 - t451) * qJD(1) - t648;
t657 = t464 * t179 + t461 * t180 - t240 * t679;
t867 = t464 * t128 + t461 * t129 + t301 * t677 + t304 * t679 + t306 * t676 - t424 + t657;
t613 = rSges(5,1) * t464 - rSges(5,3) * t747;
t265 = rSges(5,2) * t745 + t613;
t865 = qJD(1) * t265;
t669 = (t565 * t747 - t567 * t745 - t292) * qJD(1);
t864 = -t669 + t876 - t884;
t863 = t877 + t880;
t862 = t461 * t885 + t464 * t887;
t861 = t461 * t887 - t464 * t885;
t860 = qJD(2) * t892 - t184 * t463 - t186 * t460 + t435 * t899 - t436 * t901;
t859 = qJD(2) * t893 + t183 * t463 + t185 * t460 + t435 * t900 + t436 * t902;
t502 = t279 * t464 - t280 * t461;
t832 = t461 * (-t380 * t464 + t282) - t464 * (-Icges(3,2) * t738 + t281 - t421);
t856 = -t888 * t435 + t436 * t891 - t460 * t832 + t502 * t463;
t695 = t382 + t574;
t696 = -t380 + t383;
t855 = (-t435 * t890 + t436 * t889 - t460 * t695 + t463 * t696) * qJD(1);
t854 = t896 + t904;
t853 = t908 * qJD(1);
t851 = -t292 + t910;
t587 = rSges(6,1) * t319 - rSges(6,2) * t318;
t165 = rSges(6,3) * t745 - t587;
t586 = rSges(6,1) * t459 + rSges(6,2) * t462;
t804 = rSges(6,3) * t435;
t498 = t436 * t586 - t804;
t815 = pkin(7) * t435;
t542 = t617 - t815;
t512 = t542 * t464;
t849 = qJD(2) * t512 - t165 * t408 + t335 * t498;
t847 = (-qJD(2) * t585 - t424) * t461;
t365 = qJD(1) * t406;
t846 = qJD(1) * t240 - t365;
t452 = t461 * pkin(4);
t353 = pkin(7) * t744 + t452;
t357 = pkin(4) * t464 - pkin(7) * t745;
t716 = t240 - t406;
t655 = -t304 + t716;
t693 = t387 + t437;
t41 = (t357 + t655) * qJD(1) + t693 + t849;
t163 = t317 * rSges(6,1) + t316 * rSges(6,2) + rSges(6,3) * t744;
t375 = pkin(7) * t639;
t314 = t350 * t677;
t650 = -t314 - t689;
t393 = qJ(4) * t746;
t309 = pkin(3) * t744 + t393;
t713 = t241 + t407;
t652 = t309 + t713;
t42 = t631 + t163 * t408 + t498 * t334 - t375 + (t353 + t652) * qJD(1) + t650;
t844 = t41 * t464 + t42 * t461;
t290 = (Icges(6,2) * t459 - t780) * t436;
t489 = t334 * (-Icges(6,2) * t317 + t160 + t298) - t335 * (Icges(6,2) * t319 - t161 + t299) + t408 * (-t496 + t290);
t295 = (-Icges(6,1) * t462 + t781) * t436;
t490 = t334 * (-Icges(6,1) * t316 + t157 + t782) - t335 * (-Icges(6,1) * t318 + t159 - t300) + t408 * (-t495 - t295);
t829 = m(5) / 0.2e1;
t828 = m(6) / 0.2e1;
t664 = qJD(2) * qJD(5);
t627 = t435 * t664;
t217 = qJD(1) * t334 - t461 * t627;
t827 = t217 / 0.2e1;
t218 = qJD(1) * t335 - t464 * t627;
t826 = t218 / 0.2e1;
t825 = -t334 / 0.2e1;
t824 = t334 / 0.2e1;
t823 = -t335 / 0.2e1;
t822 = t335 / 0.2e1;
t821 = -t408 / 0.2e1;
t820 = t408 / 0.2e1;
t819 = t461 / 0.2e1;
t818 = -t464 / 0.2e1;
t564 = t157 * t462 + t160 * t459;
t605 = qJD(1) * t435 + qJD(5);
t497 = -t461 * t605 + t636;
t547 = t408 * t459;
t135 = t462 * t497 - t464 * t547;
t548 = t462 * t408;
t136 = t459 * t497 + t464 * t548;
t72 = Icges(6,5) * t136 + Icges(6,6) * t135 - Icges(6,3) * t866;
t74 = Icges(6,4) * t136 + Icges(6,2) * t135 - Icges(6,6) * t866;
t76 = Icges(6,1) * t136 + Icges(6,4) * t135 - Icges(6,5) * t866;
t8 = (qJD(2) * t564 + t72) * t435 + (qJD(2) * t154 - t459 * t76 - t462 * t74 + (t157 * t459 - t160 * t462) * qJD(5)) * t436;
t814 = t8 * t334;
t546 = t605 * t464;
t678 = qJD(2) * t436;
t133 = t462 * t546 + (t462 * t678 - t547) * t461;
t134 = t461 * t548 + (t546 + t637) * t459;
t510 = -t639 + t644;
t71 = Icges(6,5) * t134 + Icges(6,6) * t133 + Icges(6,3) * t510;
t73 = Icges(6,4) * t134 + Icges(6,2) * t133 + Icges(6,6) * t510;
t75 = Icges(6,1) * t134 + Icges(6,4) * t133 + Icges(6,5) * t510;
t9 = (qJD(2) * t563 + t71) * t435 + (qJD(2) * t156 - t459 * t75 - t462 * t73 + (t159 * t459 + t161 * t462) * qJD(5)) * t436;
t813 = t9 * t335;
t227 = Icges(6,3) * t436 + t435 * t569;
t287 = (-Icges(6,5) * t462 + Icges(6,6) * t459) * t436;
t130 = qJD(2) * t227 + qJD(5) * t287;
t229 = Icges(6,6) * t436 + t435 * t570;
t131 = qJD(2) * t229 + qJD(5) * t290;
t231 = Icges(6,5) * t436 + t435 * t575;
t132 = qJD(2) * t231 + qJD(5) * t295;
t560 = -t459 * t496 - t462 * t495;
t19 = (qJD(2) * t560 + t130) * t435 + (-qJD(2) * t494 - t131 * t462 - t132 * t459 + (-t459 * t495 + t462 * t496) * qJD(5)) * t436;
t626 = t436 * t664;
t83 = -t435 * t494 - t436 * t560;
t809 = t19 * t408 + t83 * t626;
t808 = rSges(3,1) * t463;
t807 = rSges(4,1) * t436;
t806 = rSges(5,2) * t436;
t805 = rSges(5,3) * t435;
t803 = rSges(6,3) * t436;
t234 = t435 * t586 + t803;
t305 = (-rSges(6,1) * t462 + rSges(6,2) * t459) * t436;
t137 = qJD(2) * t234 + qJD(5) * t305;
t434 = pkin(4) * t679;
t226 = -pkin(7) * t866 + t434;
t675 = qJD(4) * t436;
t261 = qJD(2) * t354 - t675;
t465 = qJD(2) ^ 2;
t592 = -pkin(7) * t436 - t816;
t535 = t592 * t465;
t665 = qJD(2) * qJD(4);
t628 = t436 * t665;
t358 = qJD(1) * (-pkin(1) * t680 + t433);
t666 = qJD(1) * qJD(3);
t656 = qJD(1) * t179 + t461 * t666 + t358;
t537 = t461 * t628 + t656 + (t128 + t387) * qJD(1);
t735 = t136 * rSges(6,1) + t135 * rSges(6,2);
t78 = -rSges(6,3) * t866 + t735;
t10 = qJD(1) * t226 - t137 * t334 + t218 * t498 + t408 * t78 + t461 * t535 + (qJD(1) * t512 + t163 * t671 - t261 * t461) * qJD(2) + t537;
t802 = t10 * t464;
t225 = qJD(1) * t353 - t375;
t694 = qJD(1) * t416 + t464 * t666;
t602 = qJD(1) * t314 + t464 * t628 + t694;
t364 = t407 * qJD(1);
t728 = -t180 - t364;
t659 = -t129 + t728;
t588 = rSges(6,1) * t134 + rSges(6,2) * t133;
t77 = rSges(6,3) * t510 + t588;
t11 = -t165 * t626 - t137 * t335 - t217 * t498 - t408 * t77 + (-qJD(2) * t261 + t535) * t464 + (-t225 + (pkin(7) * qJD(2) - qJD(4)) * t747 + t659) * qJD(1) + t602;
t801 = t11 * t461;
t797 = t41 * t461;
t794 = t42 * t464;
t448 = t461 * rSges(5,1);
t447 = t461 * rSges(3,3);
t60 = t154 * t435 - t436 * t564;
t792 = t60 * t218;
t791 = t61 * t217;
t352 = rSges(4,1) * t435 + rSges(4,2) * t436;
t616 = -t352 - t817;
t591 = t616 * t464;
t539 = qJD(2) * t591;
t508 = t437 + t539;
t86 = (-t262 + t716) * qJD(1) + t508;
t790 = t86 * t352;
t688 = rSges(3,2) * t741 + t464 * rSges(3,3);
t312 = rSges(3,1) * t738 - t688;
t389 = rSges(3,1) * t460 + rSges(3,2) * t463;
t640 = t389 * t676;
t170 = -t640 + (-t312 - t406) * qJD(1);
t764 = t170 * t461;
t763 = t170 * t464;
t641 = t389 * t677;
t313 = rSges(3,1) * t736 - rSges(3,2) * t740 + t447;
t704 = t313 + t407;
t171 = qJD(1) * t704 - t641;
t337 = t389 * t464;
t762 = t171 * t337;
t730 = t163 + t353;
t729 = t165 - t357;
t714 = -t241 - t309;
t355 = t805 - t806;
t326 = t355 * qJD(2);
t708 = -t261 - t326;
t707 = -qJD(1) * t306 + t436 * t674;
t703 = t350 * t680 + t417;
t697 = rSges(4,2) * t647 + rSges(4,3) * t679;
t690 = rSges(3,2) * t643 + rSges(3,3) * t679;
t687 = t461 ^ 2 + t464 ^ 2;
t667 = qJD(1) * qJD(2);
t663 = -rSges(6,3) - pkin(3) - pkin(7);
t662 = pkin(2) * t740;
t661 = t465 * t816;
t660 = qJD(2) * t816;
t629 = t464 * t667;
t658 = t179 * t676 + t180 * t677 - t240 * t629;
t264 = -rSges(5,2) * t744 + rSges(5,3) * t746 + t448;
t654 = -t264 + t714;
t653 = -t353 + t714;
t649 = t413 + t309;
t642 = t352 * t677;
t634 = t463 * t676;
t630 = -pkin(1) - t808;
t624 = t679 / 0.2e1;
t623 = t678 / 0.2e1;
t622 = -t677 / 0.2e1;
t620 = -t676 / 0.2e1;
t619 = t676 / 0.2e1;
t356 = -rSges(4,2) * t435 + t807;
t615 = -t356 - t816;
t603 = t461 * t304 + t464 * t309 + t721;
t601 = rSges(5,1) * t679 + rSges(5,2) * t866 + rSges(5,3) * t636;
t600 = t376 + t648;
t599 = t687 * t817;
t594 = qJD(5) * t623;
t593 = -t261 - t660;
t589 = -rSges(3,2) * t460 + t808;
t584 = t45 * t464 + t46 * t461;
t583 = t45 * t461 - t46 * t464;
t582 = t461 * t48 + t464 * t47;
t581 = t461 * t47 - t464 * t48;
t580 = t461 * t61 + t464 * t60;
t579 = t461 * t60 - t464 * t61;
t578 = -qJD(1) * t304 + t693 + t846;
t562 = t163 * t461 - t165 * t464;
t561 = -t171 * t461 - t763;
t543 = -t326 + t593;
t327 = t356 * qJD(2);
t536 = -qJD(2) * t327 - t661;
t336 = t389 * t461;
t303 = t352 * t461;
t526 = t498 + t542;
t525 = t304 * t677 + t309 * t676 - t675 + t722;
t166 = (t312 * t461 + t313 * t464) * qJD(2);
t514 = qJD(2) * t540;
t513 = t128 * t676 + t129 * t677 + t304 * t629 + t435 * t665 + t658;
t507 = -t154 * t334 + t156 * t335 + t408 * t494;
t506 = (Icges(6,5) * t316 - Icges(6,6) * t317) * t334 - (Icges(6,5) * t318 + Icges(6,6) * t319) * t335 + t287 * t408;
t501 = -t354 - t429 - t805;
t500 = -pkin(7) * t678 - t137 + t593;
t499 = t436 * t506;
t30 = t163 * t335 + t165 * t334 + (t353 * t464 - t357 * t461) * qJD(2) + t525;
t481 = t30 * t562 - (t794 - t797) * t498;
t468 = (t494 * t464 + t564) * t334 - (t494 * t461 + t563) * t335 + (t227 + t560) * t408;
t467 = t468 * t436;
t362 = t589 * qJD(2);
t308 = t352 * t464;
t307 = t585 * t464;
t276 = t636 - t647;
t274 = t687 * t435 * qJD(2);
t204 = t498 * t464;
t203 = t498 * t461;
t202 = t496 * t464;
t201 = t496 * t461;
t200 = t495 * t464;
t199 = t495 * t461;
t196 = rSges(6,1) * t318 + rSges(6,2) * t319;
t195 = rSges(6,1) * t316 - rSges(6,2) * t317;
t188 = -qJD(2) * t336 + (t464 * t589 + t447) * qJD(1);
t187 = -rSges(3,2) * t634 + (-t463 * t680 - t635) * rSges(3,1) + t690;
t153 = -rSges(5,3) * t647 + t601;
t152 = qJD(2) * t302 + (t355 * t464 + t448) * qJD(1);
t151 = -qJD(2) * t303 + (t356 * t464 + t446) * qJD(1);
t150 = -rSges(4,1) * t866 - rSges(4,2) * t636 + t697;
t97 = -t362 * t676 + (-t188 - t364 + t641) * qJD(1);
t96 = -t362 * t677 + t358 + (t187 - t640) * qJD(1);
t87 = -t642 + (t263 + t713) * qJD(1) - t689;
t65 = -t847 + (t264 + t652) * qJD(1) + t650;
t64 = t514 + (t265 + t655) * qJD(1) + t693;
t59 = (t264 * t464 - t265 * t461) * qJD(2) + t525;
t54 = t536 * t464 + (-t151 + t642 + t728) * qJD(1) + t694;
t53 = t536 * t461 + (t150 + t539) * qJD(1) + t656;
t29 = (qJD(2) * t708 - t661) * t464 + (-t152 + t659 + t847) * qJD(1) + t602;
t28 = -t461 * t661 + qJD(1) * t153 + (qJD(1) * t540 + t461 * t708) * qJD(2) + t537;
t17 = t130 * t744 + t131 * t316 + t132 * t317 - t135 * t495 - t136 * t496 + t494 * t866;
t16 = t130 * t745 + t131 * t318 - t132 * t319 - t133 * t495 - t134 * t496 - t494 * t510;
t15 = t334 * t60 - t335 * t61 + t408 * t83;
t14 = (t152 * t461 + t153 * t464 + (-t265 * t464 + t461 * t654) * qJD(1)) * qJD(2) + t513;
t7 = t135 * t159 - t136 * t161 - t156 * t866 + t316 * t73 + t317 * t75 + t71 * t744;
t6 = t135 * t157 + t136 * t160 - t154 * t866 + t316 * t74 + t317 * t76 + t72 * t744;
t5 = t133 * t159 - t134 * t161 + t156 * t510 + t318 * t73 - t319 * t75 + t71 * t745;
t4 = t133 * t157 + t134 * t160 + t154 * t510 + t318 * t74 - t319 * t76 + t72 * t745;
t3 = -t163 * t217 + t165 * t218 + t334 * t77 + t335 * t78 + (t225 * t461 + t226 * t464 + (-t357 * t464 + t461 * t653) * qJD(1)) * qJD(2) + t513;
t2 = t17 * t408 + t217 * t46 + t218 * t45 + t334 * t6 - t335 * t7 + t626 * t79;
t1 = t16 * t408 + t217 * t48 + t218 * t47 + t334 * t4 - t335 * t5 + t626 * t80;
t18 = [t12 * t822 + t17 * t824 + t79 * t826 + t80 * t827 + t792 / 0.2e1 + t791 / 0.2e1 + t814 / 0.2e1 - t813 / 0.2e1 + t809 + (t16 + t12) * t823 + (-t907 * qJD(2) + t360 * t463 + t361 * t460 + t913 * t435 + t914 * t436) * qJD(1) + (t11 * (t357 + t587 + t691) + t41 * (t375 - t588 + t600) + t10 * (t649 + t730) + (t11 * (-t354 - t803) - t41 * t424 - t10 * t458) * t461 + ((-t765 + t804) * t797 + (t435 * t663 - t817) * t794) * qJD(2) + (t41 * (-t393 - t413 - t452) + t844 * t436 * t663) * qJD(1) + (t369 + t41 + t434 - t578 + t693 + t735 + (-qJ(4) * t747 - t357 + t691) * qJD(1) - t849) * t42) * m(6) + (-(t514 + t578 - t64 + t865) * t65 + t29 * (t613 + t691) + t64 * t600 + t28 * (t264 + t649) + t65 * (-pkin(3) * t638 + t544 + t601 + t698) + (t29 * (-t354 + t806) - t28 * t458 + (-t424 + (-t793 + (-rSges(5,3) - qJ(4)) * t436) * qJD(2)) * t64) * t461 + ((-t64 * rSges(5,1) + t501 * t65) * t461 + (t64 * (t501 + t806) - t65 * t458) * t464) * qJD(1)) * m(5) + (-(-qJD(1) * t262 + t508 + t846 - t86) * t87 + t54 * (-t262 + t691) + t86 * t648 + t53 * (t263 + t606) + t87 * (t437 + t697) + (t461 * t790 + t591 * t87) * qJD(2) + ((-t86 * rSges(4,3) + t87 * (-t429 - t807)) * t461 + (t86 * (-t356 - t429) - t87 * t458) * t464) * qJD(1)) * m(4) + (t97 * (t461 * t630 + t453 + t688) + t96 * t704 + t171 * (t433 + t690) + (t389 * t764 - t762) * qJD(2) + ((-pkin(1) - t589) * t763 + (t170 * (-rSges(3,3) - pkin(6)) + t171 * t630) * t461) * qJD(1) - (-qJD(1) * t312 - t170 - t365 - t640) * t171) * m(3) + ((t852 * t461 + ((t909 + t926) * t464 + t883 + t905 + t916) * t464) * qJD(2) + t880) * t619 + (t859 + t862) * t677 / 0.2e1 + (0.2e1 * t669 + ((t464 * t854 - t852 + t881) * t464 + (t461 * t854 + t882 + t906) * t461) * qJD(2) + t864 + t884) * t622 + (-t860 + t861 + t863) * t620 + ((t858 - t903) * t461 + (t857 + t897) * t464) * t667 / 0.2e1; (qJD(1) * t580 + t461 * t8 - t464 * t9) * t820 + (((-t200 * t462 - t202 * t459 + t154) * t334 - (-t199 * t462 - t201 * t459 + t156) * t335 + (-t229 * t462 - t231 * t459 - t494) * t408 + t83 * qJD(5)) * t436 + (-qJD(5) * t580 + t468) * t435) * t821 + ((t200 * t318 - t202 * t319) * t334 - (t199 * t318 - t201 * t319) * t335 + (t229 * t318 - t231 * t319) * t408 + (t436 * t80 - t47 * t746) * qJD(5) + ((-qJD(5) * t48 + t507) * t435 + t467) * t461) * t822 + (qJD(1) * t582 + t4 * t461 - t464 * t5) * t823 + (qJD(1) * t584 + t461 * t6 - t464 * t7) * t824 + ((t200 * t316 + t202 * t317) * t334 - (t199 * t316 + t201 * t317) * t335 + (t229 * t316 + t231 * t317) * t408 + (t436 * t79 - t46 * t747) * qJD(5) + ((-qJD(5) * t45 + t507) * t435 + t467) * t464) * t825 + t583 * t826 + t581 * t827 - t15 * t671 / 0.2e1 + t579 * t594 + (t12 * t464 + t13 * t461) * t672 / 0.2e1 + (-t41 * (-t203 * t408 - t234 * t335 - t873) - t42 * (-pkin(7) * t646 - qJD(1) * t662 + t204 * t408 - t234 * t334 + t707) - ((t163 * t42 - t165 * t41) * t436 + t481 * t435) * qJD(5) - t844 * (-t354 + t592) * qJD(2) + t41 * t703 + t3 * t603 + (-qJD(1) * t41 * t498 + t10 * t526 + t3 * t729 + t42 * t500) * t461 + (t3 * t730 + t41 * t500 + (qJD(1) * t42 + t11) * t526) * t464 + (-t203 * t334 - t204 * t335 - (-t687 * t815 - t599) * qJD(2) + (t225 + t77 + (-t163 + t653) * qJD(1)) * t461 + (qJD(1) * t729 + t226 + t78) * t464 + t867) * t30) * m(6) + (t14 * t603 + (t14 * t264 + t29 * t598) * t464 + (-t14 * t265 + t28 * t598) * t461 + (t543 * t461 - t707 + (-t307 + t662 + t540) * qJD(1)) * t65 + (t543 * t464 + t703 + t873) * t64 - (t65 * t461 + t64 * t464) * qJD(2) * (-t354 - t355 - t816) + ((t153 - t865) * t464 + (qJD(1) * t654 + t152) * t461 - (t302 * t461 + t307 * t464 - t599) * qJD(2) + t867) * t59) * m(5) + (-(t86 * t303 + t87 * (-t308 - t662)) * qJD(1) - (-t81 * t599 + (-t81 * t308 + t615 * t86) * t464 + (-t81 * t303 + t615 * t87) * t461) * qJD(2) + t54 * t591 - t86 * pkin(2) * t634 + (t150 * t676 + t151 * t677 + t658) * t511 + t81 * t657 + (-t86 * t327 + t81 * t150 + (t262 * t871 + t616 * t87) * qJD(1)) * t464 + (t53 * t616 + t87 * (-t327 - t660) + t81 * t151 + (t790 + t871 * (-t241 - t263)) * qJD(1)) * t461) * m(4) + (0.2e1 * t166 * (t187 * t464 + t188 * t461 + (t312 * t464 - t313 * t461) * qJD(1)) + t561 * t362 + (-t96 * t461 - t97 * t464 + (-t171 * t464 + t764) * qJD(1)) * t389 - (t170 * t336 - t762) * qJD(1) - (t166 * (-t336 * t461 - t337 * t464) + t561 * t589) * qJD(2)) * m(3) - ((t435 * t891 + t888 * t436 + t502 * t460 + t463 * t832) * qJD(2) + (t889 * t435 + t890 * t436 + t460 * t696 + t463 * t695) * qJD(1)) * qJD(1) / 0.2e1 + (t860 * t464 + t859 * t461 + (t858 * t461 + t857 * t464) * qJD(1)) * qJD(1) / 0.2e1 + ((-t677 * t851 + t853) * t461 + ((t461 * t850 + t856) * qJD(2) + t855) * t464) * t622 + ((-t676 * t850 - t853) * t464 + ((t464 * t851 + t856) * qJD(2) + t855) * t461) * t619 + (t862 * qJD(1) + t2 + ((t881 * qJD(1) + t464 * t878) * t464 + (t874 * t461 + t882 * qJD(1) + (-t875 + t879) * t464) * t461) * t872) * t819 + (t861 * qJD(1) + t1 + ((t883 * qJD(1) + t875 * t464) * t464 + (t879 * t461 + t898 * qJD(1) + (-t874 + t878) * t464) * t461) * t872) * t818 + (t13 + t864 + t876) * t680 / 0.2e1 + (t12 + t863 + t877) * t624; 0.2e1 * (-t802 / 0.2e1 + t801 / 0.2e1) * m(6) + 0.2e1 * (t28 * t818 + t29 * t819) * m(5) + 0.2e1 * (t53 * t818 + t54 * t819) * m(4); -m(5) * (t274 * t59 + t275 * t65 + t276 * t64) - m(6) * (t274 * t30 + t275 * t42 + t276 * t41) + 0.2e1 * ((t64 * t676 + t65 * t677 - t14) * t829 + (t41 * t676 + t42 * t677 - t3) * t828) * t436 + 0.2e1 * ((qJD(2) * t59 + t28 * t461 + t29 * t464 - t64 * t680 + t65 * t679) * t829 + (qJD(2) * t30 + t10 * t461 + t11 * t464 - t41 * t680 + t42 * t679) * t828) * t435; t2 * t744 / 0.2e1 + (t435 * t79 + t436 * t584) * t826 + ((-qJD(2) * t584 + t17) * t435 + (-qJD(1) * t583 + qJD(2) * t79 + t461 * t7 + t464 * t6) * t436) * t824 + t1 * t745 / 0.2e1 + (t435 * t80 + t436 * t582) * t827 + ((-qJD(2) * t582 + t16) * t435 + (-qJD(1) * t581 + qJD(2) * t80 + t4 * t464 + t461 * t5) * t436) * t823 + t15 * t623 + t435 * (t791 + t792 + t809 - t813 + t814) / 0.2e1 + (t435 * t83 + t436 * t580) * t594 + ((-qJD(2) * t580 + t19) * t435 + (-qJD(1) * t579 + qJD(2) * t83 + t461 * t9 + t464 * t8) * t436) * t820 + (t489 * t316 - t317 * t490 + t464 * t499) * t825 + (t318 * t489 + t319 * t490 + t461 * t499) * t822 + (t506 * t435 + (t490 * t459 - t462 * t489) * t436) * t821 + (t435 * t622 + t436 * t624) * t13 + (-t645 / 0.2e1 + t435 * t620) * t12 + ((qJD(2) * t481 + t10 * t163 - t11 * t165 - t41 * t77 + t42 * t78) * t435 + (t41 * (-qJD(2) * t165 + t137 * t461) + t42 * (qJD(2) * t163 - t137 * t464) - t3 * t562 + t30 * (-t163 * t679 - t165 * t680 - t461 * t78 + t464 * t77) - (qJD(1) * t844 + t801 - t802) * t498) * t436 - t41 * (-t196 * t408 - t305 * t335) - t42 * (t195 * t408 - t305 * t334) - t30 * (t195 * t335 + t196 * t334)) * m(6);];
tauc = t18(:);

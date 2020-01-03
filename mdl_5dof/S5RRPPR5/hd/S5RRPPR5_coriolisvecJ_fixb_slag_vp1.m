% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR5
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:54
% EndTime: 2019-12-31 19:29:46
% DurationCPUTime: 44.03s
% Computational Cost: add. (16947->872), mult. (25192->1085), div. (0->0), fcn. (22323->8), ass. (0->455)
t910 = Icges(3,3) + Icges(4,3);
t463 = qJ(2) + pkin(8);
t439 = sin(t463);
t440 = cos(t463);
t468 = sin(qJ(2));
t470 = cos(qJ(2));
t909 = Icges(3,5) * t470 + Icges(4,5) * t440 - Icges(3,6) * t468 - Icges(4,6) * t439;
t348 = Icges(5,4) * t440 + Icges(5,6) * t439;
t469 = sin(qJ(1));
t471 = cos(qJ(1));
t258 = Icges(5,2) * t469 + t348 * t471;
t891 = t910 * t469 + t909 * t471;
t908 = t258 + t891;
t749 = Icges(4,4) * t439;
t349 = Icges(4,2) * t440 + t749;
t429 = Icges(5,5) * t439;
t907 = Icges(5,3) * t440 + t349 - t429;
t743 = Icges(5,5) * t440;
t351 = Icges(5,1) * t439 - t743;
t430 = Icges(4,4) * t440;
t906 = Icges(4,1) * t439 + t351 + t430;
t905 = t910 * t471;
t705 = t469 * t470;
t707 = t468 * t469;
t709 = t440 * t469;
t711 = t439 * t469;
t876 = -Icges(3,5) * t705 - Icges(4,5) * t709 + Icges(3,6) * t707 + Icges(4,6) * t711 + t905;
t565 = Icges(5,1) * t440 + t429;
t261 = -Icges(5,4) * t471 + t469 * t565;
t404 = Icges(4,4) * t711;
t744 = Icges(4,5) * t471;
t263 = Icges(4,1) * t709 - t404 - t744;
t898 = t261 + t263;
t262 = Icges(5,4) * t469 + t471 * t565;
t354 = Icges(4,1) * t440 - t749;
t264 = Icges(4,5) * t469 + t354 * t471;
t897 = t262 + t264;
t344 = Icges(5,3) * t439 + t743;
t562 = -Icges(4,2) * t439 + t430;
t904 = -t344 + t562;
t903 = t354 + t565;
t735 = Icges(4,6) * t471;
t259 = Icges(4,4) * t709 - Icges(4,2) * t711 - t735;
t736 = Icges(3,6) * t471;
t295 = Icges(3,4) * t705 - Icges(3,2) * t707 - t736;
t902 = t259 * t439 + t295 * t468;
t892 = Icges(3,5) * t468 + Icges(3,6) * t470 + (Icges(4,6) - Icges(5,6)) * t440 + (Icges(5,4) + Icges(4,5)) * t439;
t750 = Icges(3,4) * t468;
t391 = Icges(3,1) * t470 - t750;
t298 = Icges(3,5) * t469 + t391 * t471;
t901 = -t264 * t709 - t298 * t705;
t253 = -Icges(5,6) * t471 + t344 * t469;
t900 = -t253 + t259;
t708 = t440 * t471;
t403 = Icges(5,5) * t708;
t710 = t439 * t471;
t734 = Icges(5,6) * t469;
t254 = Icges(5,3) * t710 + t403 + t734;
t260 = Icges(4,6) * t469 + t471 * t562;
t899 = -t254 + t260;
t896 = t348 + t909;
t895 = t907 * qJD(2);
t894 = t906 * qJD(2);
t425 = Icges(3,4) * t707;
t745 = Icges(3,5) * t471;
t297 = Icges(3,1) * t705 - t425 - t745;
t893 = t263 * t440 + t297 * t470 - t902;
t388 = Icges(3,2) * t470 + t750;
t453 = Icges(3,4) * t470;
t390 = Icges(3,1) * t468 + t453;
t881 = t388 * t468 - t390 * t470 + t907 * t439 - t906 * t440;
t890 = t904 * qJD(2);
t889 = t903 * qJD(2);
t575 = -t254 * t711 + t258 * t471 - t262 * t709;
t563 = -Icges(3,2) * t468 + t453;
t296 = Icges(3,6) * t469 + t471 * t563;
t885 = -t471 * t891 - t901;
t856 = -t260 * t711 - t296 * t707 + t885;
t886 = -t575 + t856;
t884 = t260 * t439 + t296 * t468;
t704 = t470 * t471;
t883 = t254 * t710 + t298 * t704 + t908 * t469 + t897 * t708;
t257 = -Icges(5,2) * t471 + t348 * t469;
t236 = t469 * t257;
t882 = -t253 * t710 - t297 * t704 + t469 * t876 - t898 * t708 - t236;
t819 = t892 * t471;
t818 = t892 * t469;
t880 = t895 * t471 + (t469 * t562 - t253 - t735) * qJD(1);
t879 = t895 * t469 + (t344 * t471 - t260 + t734) * qJD(1);
t878 = -t894 * t471 + (-t354 * t469 - t261 + t744) * qJD(1);
t877 = -qJD(1) * t897 + t469 * t894;
t724 = t257 * t471;
t557 = t253 * t439 + t261 * t440;
t814 = t469 * t557;
t88 = -t724 + t814;
t874 = t469 * t893 + t471 * t876 + t88;
t706 = t468 * t471;
t853 = -t259 * t710 - t295 * t706 - t882;
t852 = -t260 * t710 - t296 * t706 + t883;
t873 = t469 * t881 + t819;
t872 = -t471 * t881 + t818;
t871 = t892 * qJD(2);
t870 = t254 * t439 + t298 * t470 + t440 * t897 - t884;
t869 = -t557 - t893;
t822 = t296 * t470 + t298 * t468 + t439 * t897 + t440 * t899;
t821 = t295 * t470 + t297 * t468 + t439 * t898 + t440 * t900;
t467 = sin(qJ(5));
t775 = cos(qJ(5));
t536 = t439 * t467 + t440 * t775;
t637 = qJD(2) - qJD(5);
t204 = t637 * t536;
t620 = t439 * t775;
t332 = -t440 * t467 + t620;
t650 = qJD(1) * t471;
t116 = t204 * t469 + t332 * t650;
t286 = t536 * t471;
t837 = t637 * t332;
t117 = qJD(1) * t286 - t469 * t837;
t285 = t467 * t708 - t471 * t620;
t136 = Icges(6,5) * t286 - Icges(6,6) * t285 - Icges(6,3) * t469;
t252 = Icges(6,4) * t286;
t139 = -Icges(6,2) * t285 - Icges(6,6) * t469 + t252;
t251 = Icges(6,4) * t285;
t142 = Icges(6,1) * t286 - Icges(6,5) * t469 - t251;
t283 = t467 * t709 - t469 * t620;
t284 = t536 * t469;
t651 = qJD(1) * t469;
t114 = t204 * t471 - t332 * t651;
t115 = -t471 * t837 - t536 * t651;
t65 = Icges(6,5) * t115 + Icges(6,6) * t114 - Icges(6,3) * t650;
t67 = Icges(6,4) * t115 + Icges(6,2) * t114 - Icges(6,6) * t650;
t69 = Icges(6,1) * t115 + Icges(6,4) * t114 - Icges(6,5) * t650;
t10 = t116 * t139 + t117 * t142 - t136 * t651 - t283 * t67 + t284 * t69 + t471 * t65;
t384 = t637 * t469;
t647 = qJD(2) * t471;
t385 = -qJD(5) * t471 + t647;
t134 = Icges(6,5) * t284 - Icges(6,6) * t283 + Icges(6,3) * t471;
t746 = Icges(6,4) * t284;
t137 = -Icges(6,2) * t283 + Icges(6,6) * t471 + t746;
t747 = Icges(6,4) * t283;
t141 = -Icges(6,1) * t284 - Icges(6,5) * t471 + t747;
t47 = -t134 * t469 - t285 * t137 - t141 * t286;
t576 = t136 * t469 + t285 * t139 - t286 * t142;
t206 = Icges(6,5) * t332 - Icges(6,6) * t536;
t322 = Icges(6,4) * t332;
t209 = -Icges(6,2) * t536 + t322;
t321 = Icges(6,4) * t536;
t212 = Icges(6,1) * t332 - t321;
t64 = -t206 * t469 - t209 * t285 + t212 * t286;
t12 = t64 * qJD(1) - t384 * t576 - t385 * t47;
t100 = Icges(6,1) * t204 + Icges(6,4) * t837;
t98 = Icges(6,5) * t204 + Icges(6,6) * t837;
t99 = Icges(6,4) * t204 + Icges(6,2) * t837;
t13 = t100 * t286 + t114 * t209 + t115 * t212 - t206 * t650 - t285 * t99 - t469 * t98;
t14 = t100 * t284 + t116 * t209 + t117 * t212 - t206 * t651 - t283 * t99 + t471 * t98;
t68 = Icges(6,4) * t117 + Icges(6,2) * t116 - Icges(6,6) * t651;
t70 = Icges(6,1) * t117 + Icges(6,4) * t116 - Icges(6,5) * t651;
t17 = t137 * t837 - t141 * t204 + t332 * t70 - t536 * t68;
t18 = t139 * t837 + t142 * t204 + t332 * t69 - t536 * t67;
t589 = qJD(1) * t637;
t364 = t469 * t589;
t365 = t471 * t589;
t45 = t134 * t471 - t137 * t283 - t141 * t284;
t46 = t471 * t136 - t283 * t139 + t284 * t142;
t496 = qJD(1) * (Icges(6,1) * t536 + t209 + t322) + t384 * (Icges(6,1) * t285 + t139 + t252) - t385 * (Icges(6,1) * t283 + t137 + t746);
t510 = qJD(1) * (-Icges(6,5) * t536 - Icges(6,6) * t332) + (Icges(6,5) * t283 + Icges(6,6) * t284) * t385 - (Icges(6,5) * t285 + Icges(6,6) * t286) * t384;
t63 = t206 * t471 - t209 * t283 + t212 * t284;
t66 = Icges(6,5) * t117 + Icges(6,6) * t116 - Icges(6,3) * t651;
t7 = t114 * t137 - t115 * t141 - t134 * t650 - t285 * t68 + t286 * t70 - t469 * t66;
t769 = -qJD(1) / 0.2e1;
t77 = -t137 * t536 - t141 * t332;
t78 = -t139 * t536 + t142 * t332;
t780 = t385 / 0.2e1;
t789 = qJD(1) * (Icges(6,2) * t332 - t212 + t321) - t384 * (-Icges(6,2) * t286 + t142 - t251) + t385 * (-Icges(6,2) * t284 - t141 - t747);
t8 = t114 * t139 + t115 * t142 - t136 * t650 - t285 * t67 + t286 * t69 - t469 * t65;
t826 = t47 / 0.2e1 + t46 / 0.2e1;
t9 = t116 * t137 - t117 * t141 - t134 * t651 - t283 * t68 + t284 * t70 + t471 * t66;
t868 = (-t17 * t471 + t18 * t469 + (t469 * t77 + t471 * t78) * qJD(1)) * t769 - t469 * (qJD(1) * t13 + t384 * t8 - t385 * t7) / 0.2e1 + t471 * (qJD(1) * t14 + t10 * t384 - t385 * t9) / 0.2e1 - (qJD(1) * t63 + t384 * t46 - t385 * t45) * t651 / 0.2e1 - t12 * t650 / 0.2e1 + (t469 * t576 + t471 * t826) * t365 + (t45 * t471 - t469 * t826) * t364 - (t285 * t789 - t286 * t496 + (-t576 * qJD(1) - t7) * t471 + (t47 * qJD(1) - t510 + t8) * t469) * t384 / 0.2e1 + (t283 * t789 - t284 * t496 + (t45 * qJD(1) + t10) * t469 + (t46 * qJD(1) + t510 - t9) * t471) * t780;
t867 = rSges(3,2) * t468;
t866 = t900 * t471 + (-Icges(5,1) * t710 + t351 * t471 + t403 - t899) * t469;
t865 = t906 + t904;
t864 = -t907 + t903;
t863 = (Icges(4,2) * t709 + t404 - t898) * t471 + (-t349 * t471 + t897) * t469;
t367 = t563 * qJD(2);
t368 = t391 * qJD(2);
t862 = -t367 * t468 + t368 * t470 + t889 * t440 - t890 * t439 + (-t388 * t470 - t390 * t468 - t439 * t906 - t440 * t907) * qJD(2) + t892 * qJD(1);
t861 = t908 * qJD(1);
t860 = t881 * qJD(1) + qJD(2) * t896;
t859 = -t496 * t332 + t536 * t789;
t851 = t872 * qJD(1);
t525 = qJD(2) * t388;
t197 = -t471 * t525 + (-t469 * t563 + t736) * qJD(1);
t528 = qJD(2) * t390;
t199 = -t471 * t528 + (-t391 * t469 + t745) * qJD(1);
t850 = -t822 * qJD(2) - t197 * t468 + t199 * t470 + t880 * t439 + t878 * t440 + t861;
t198 = qJD(1) * t296 - t469 * t525;
t200 = qJD(1) * t298 - t469 * t528;
t808 = qJD(1) * t257;
t849 = t876 * qJD(1) + t821 * qJD(2) + t198 * t468 - t200 * t470 - t879 * t439 + t877 * t440 - t808;
t848 = (t852 * t469 - t853 * t471) * qJD(2);
t847 = (t886 * t469 - t874 * t471) * qJD(2);
t846 = t873 * qJD(1);
t845 = t869 * qJD(1) - t871 * t469 + t861;
t844 = -t808 - t871 * t471 + (-t469 * t909 - t870 + t905) * qJD(1);
t843 = 0.2e1 * qJD(2);
t276 = rSges(4,1) * t709 - rSges(4,2) * t711 - t471 * rSges(4,3);
t454 = t469 * rSges(4,3);
t278 = rSges(4,1) * t708 - rSges(4,2) * t710 + t454;
t552 = t276 * t469 + t278 * t471;
t461 = t471 * pkin(6);
t409 = pkin(1) * t469 - t461;
t466 = -qJ(3) - pkin(6);
t433 = t471 * t466;
t773 = pkin(2) * t470;
t434 = pkin(1) + t773;
t662 = -t469 * t434 - t433;
t245 = t409 + t662;
t460 = t469 * pkin(6);
t410 = t471 * pkin(1) + t460;
t414 = t471 * t434;
t590 = -t466 * t469 + t414;
t246 = t590 - t410;
t690 = -t469 * t245 + t471 * t246;
t514 = t552 + t690;
t648 = qJD(2) * t469;
t691 = -t245 * t648 + t246 * t647;
t83 = qJD(2) * t552 + t691;
t842 = qJD(2) * t514 + t83;
t214 = -rSges(6,1) * t536 - rSges(6,2) * t332;
t839 = t214 * t384;
t838 = t214 * t385;
t764 = rSges(5,1) * t439;
t358 = -rSges(5,3) * t440 + t764;
t312 = t358 * t469;
t357 = pkin(3) * t439 - qJ(4) * t440;
t774 = pkin(2) * t468;
t598 = -t357 - t774;
t582 = -t358 + t598;
t541 = t471 * t582;
t570 = rSges(6,1) * t284 - rSges(6,2) * t283;
t143 = rSges(6,3) * t471 + t570;
t355 = pkin(4) * t709 + pkin(7) * t471;
t698 = t143 + t355;
t613 = t440 * t647;
t376 = qJ(4) * t613;
t644 = qJD(4) * t471;
t394 = t439 * t644;
t513 = -t439 * t647 - t440 * t651;
t619 = t439 * t651;
t146 = pkin(3) * t513 - qJ(4) * t619 + t376 + t394;
t291 = t439 * t650 + t440 * t648;
t614 = t439 * t648;
t382 = pkin(3) * t614;
t412 = pkin(3) * t708;
t645 = qJD(4) * t469;
t608 = t439 * t645;
t147 = qJ(4) * t291 + qJD(1) * t412 - t382 + t608;
t311 = t357 * t469;
t360 = pkin(3) * t440 + qJ(4) * t439;
t314 = t360 * t469;
t315 = t357 * t471;
t428 = qJD(4) * t439;
t438 = pkin(6) * t650;
t442 = qJD(3) * t469;
t612 = t468 * t647;
t767 = pkin(1) - t434;
t193 = -pkin(2) * t612 - t438 + t442 + (t469 * t767 - t433) * qJD(1);
t417 = t648 * t774;
t660 = qJD(3) * t471 + t417;
t621 = t466 * t651 + t660;
t194 = (-t471 * t767 - t460) * qJD(1) - t621;
t630 = t471 * t193 + t469 * t194 - t245 * t650;
t836 = t471 * t146 + t469 * t147 + t311 * t648 + t314 * t650 + t315 * t647 - t428 + t630;
t835 = t876 + t884;
t834 = t357 * t651 - t440 * t644;
t361 = rSges(5,1) * t440 + rSges(5,3) * t439;
t459 = t471 * rSges(5,2);
t275 = t361 * t469 - t459;
t833 = qJD(1) * t275;
t832 = -t846 + t847;
t831 = t848 + t851;
t830 = t469 * t860 + t471 * t862;
t829 = t469 * t862 - t471 * t860;
t828 = t869 * qJD(2) - t198 * t470 - t200 * t468 + t877 * t439 + t879 * t440;
t827 = t870 * qJD(2) + t197 * t470 + t199 * t468 + t878 * t439 - t880 * t440;
t505 = t295 * t471 - t296 * t469;
t792 = t469 * (-t388 * t471 + t298) - t471 * (-Icges(3,2) * t705 + t297 - t425);
t825 = -t863 * t439 + t440 * t866 - t468 * t792 + t505 * t470;
t666 = t390 + t563;
t667 = -t388 + t391;
t824 = (-t439 * t865 + t440 * t864 - t468 * t666 + t470 * t667) * qJD(1);
t823 = t724 + t883;
t820 = t896 * qJD(1);
t176 = rSges(6,1) * t283 + rSges(6,2) * t284;
t178 = rSges(6,1) * t285 + rSges(6,2) * t286;
t816 = t176 * t384 + t178 * t385;
t813 = (qJD(2) * t358 - t428) * t469;
t372 = qJD(1) * t409;
t812 = qJD(1) * t245 - t372;
t811 = -t360 - t361;
t411 = pkin(4) * t708;
t356 = -pkin(7) * t469 + t411;
t215 = rSges(6,1) * t332 - rSges(6,2) * t536;
t772 = pkin(4) * t439;
t544 = t598 - t772;
t508 = t544 * t647;
t810 = -t385 * t215 + t508;
t687 = t245 - t409;
t628 = -t314 + t687;
t664 = t394 + t442;
t43 = (t628 - t698) * qJD(1) + t664 + t810;
t381 = pkin(4) * t614;
t680 = t286 * rSges(6,1) - t285 * rSges(6,2);
t145 = -rSges(6,3) * t469 + t680;
t318 = qJ(4) * t710 + t412;
t685 = -t246 - t318;
t626 = -t356 + t685;
t588 = -t145 + t626;
t323 = t357 * t648;
t624 = -t323 - t660;
t44 = t608 - t215 * t384 - t381 + (t410 - t588) * qJD(1) + t624;
t809 = t43 * t471 + t44 * t469;
t787 = m(5) / 0.2e1;
t786 = m(6) / 0.2e1;
t779 = t469 / 0.2e1;
t778 = -t471 / 0.2e1;
t777 = -rSges(5,1) - pkin(3);
t776 = rSges(6,3) + pkin(7);
t771 = pkin(4) * t440;
t766 = rSges(3,1) * t470;
t765 = rSges(4,1) * t440;
t456 = t469 * rSges(5,2);
t455 = t469 * rSges(3,3);
t359 = rSges(4,1) * t439 + rSges(4,2) * t440;
t597 = -t359 - t774;
t574 = t471 * t597;
t540 = qJD(2) * t574;
t512 = t442 + t540;
t86 = (-t276 + t687) * qJD(1) + t512;
t760 = t86 * t359;
t759 = -rSges(5,3) - qJ(4);
t238 = pkin(4) * t513 - pkin(7) * t650;
t703 = t115 * rSges(6,1) + t114 * rSges(6,2);
t71 = -rSges(6,3) * t650 + t703;
t758 = t238 + t71;
t239 = qJD(1) * t356 - t381;
t571 = rSges(6,1) * t117 + rSges(6,2) * t116;
t72 = -rSges(6,3) * t651 + t571;
t757 = t239 + t72;
t659 = rSges(3,2) * t707 + t471 * rSges(3,3);
t319 = rSges(3,1) * t705 - t659;
t396 = rSges(3,1) * t468 + rSges(3,2) * t470;
t615 = t396 * t647;
t184 = -t615 + (-t319 - t409) * qJD(1);
t728 = t184 * t469;
t727 = t184 * t471;
t616 = t396 * t648;
t320 = rSges(3,1) * t704 - rSges(3,2) * t706 + t455;
t675 = t320 + t410;
t185 = qJD(1) * t675 - t616;
t342 = t396 * t471;
t726 = t185 * t342;
t371 = t410 * qJD(1);
t697 = -t194 - t371;
t686 = -t246 - t278;
t646 = qJD(4) * t440;
t274 = qJD(2) * t360 - t646;
t333 = t361 * qJD(2);
t679 = -t274 - t333;
t678 = -qJD(1) * t315 + t440 * t645;
t669 = rSges(5,2) * t650 + rSges(5,3) * t613;
t668 = rSges(4,2) * t619 + rSges(4,3) * t650;
t639 = qJD(1) * qJD(3);
t665 = qJD(1) * t417 + t471 * t639;
t661 = rSges(3,3) * t650 + t651 * t867;
t658 = t469 ^ 2 + t471 ^ 2;
t649 = qJD(2) * t440;
t640 = qJD(1) * qJD(2);
t638 = qJD(2) * qJD(4);
t636 = -t466 - t776;
t635 = pkin(2) * t706;
t472 = qJD(2) ^ 2;
t634 = t472 * t773;
t633 = qJD(2) * t773;
t632 = -t147 + t697;
t606 = t471 * t640;
t631 = t193 * t647 + t194 * t648 - t245 * t606;
t363 = qJD(1) * (-pkin(1) * t651 + t438);
t629 = qJD(1) * t193 + t469 * t639 + t363;
t277 = rSges(5,1) * t708 + rSges(5,3) * t710 + t456;
t627 = -t277 + t685;
t623 = t376 + t664;
t622 = t414 + t318;
t617 = t359 * t648;
t611 = t470 * t647;
t607 = -pkin(1) - t766;
t605 = t440 * t638;
t602 = -t648 / 0.2e1;
t599 = t647 / 0.2e1;
t362 = -rSges(4,2) * t439 + t765;
t596 = -t362 - t773;
t586 = t469 * t314 + t471 * t318 + t690;
t585 = qJD(1) * t323 + t471 * t605 + t665;
t584 = t382 + t621;
t583 = t658 * t774;
t578 = -t274 - t633;
t577 = -t771 - t773;
t572 = t766 - t867;
t101 = rSges(6,1) * t204 + rSges(6,2) * t837;
t503 = -qJD(2) * t274 + t472 * t577;
t538 = t469 * t605 + t629 + (t146 + t394) * qJD(1);
t15 = -t101 * t384 - t215 * t365 + t503 * t469 + (t508 + t758) * qJD(1) + t538;
t16 = -t101 * t385 + t215 * t364 + t503 * t471 + ((pkin(4) * qJD(2) - qJD(4)) * t711 + t632 - t757) * qJD(1) + t585;
t569 = t15 * t469 + t16 * t471;
t568 = -qJD(1) * t314 + t664 + t812;
t559 = -t185 * t469 - t727;
t545 = -t333 + t578;
t334 = t362 * qJD(2);
t537 = -qJD(2) * t334 - t634;
t341 = t396 * t469;
t313 = t359 * t469;
t530 = -t215 + t544;
t529 = t314 * t648 + t318 * t647 - t646 + t691;
t180 = (t319 * t469 + t320 * t471) * qJD(2);
t517 = qJD(2) * t541;
t516 = t146 * t647 + t147 * t648 + t314 * t606 + t439 * t638 + t631;
t515 = -t360 - t771;
t509 = -t434 + t515;
t504 = -pkin(4) * t649 - t101 + t578;
t499 = -t434 + t811;
t369 = t572 * qJD(2);
t317 = t359 * t471;
t316 = t358 * t471;
t292 = t613 - t619;
t290 = t658 * t439 * qJD(2);
t202 = -qJD(2) * t341 + (t471 * t572 + t455) * qJD(1);
t201 = -rSges(3,2) * t611 + (-t470 * t651 - t612) * rSges(3,1) + t661;
t175 = -qJD(2) * t313 + (t362 * t471 + t454) * qJD(1);
t174 = -qJD(2) * t312 + (t361 * t471 + t456) * qJD(1);
t173 = rSges(4,1) * t513 - rSges(4,2) * t613 + t668;
t172 = rSges(5,1) * t513 - rSges(5,3) * t619 + t669;
t97 = -t369 * t647 + (-t202 - t371 + t616) * qJD(1);
t96 = -t369 * t648 + t363 + (t201 - t615) * qJD(1);
t87 = -t617 + (t410 - t686) * qJD(1) - t660;
t76 = -t813 + (t410 - t627) * qJD(1) + t624;
t75 = t517 + (-t275 + t628) * qJD(1) + t664;
t61 = (t275 * t469 + t277 * t471) * qJD(2) + t529;
t56 = t537 * t471 + (-t175 + t617 + t697) * qJD(1) + t665;
t55 = t537 * t469 + (t173 + t540) * qJD(1) + t629;
t32 = t143 * t384 + t145 * t385 + (t355 * t469 + t356 * t471) * qJD(2) + t529;
t31 = (qJD(2) * t679 - t634) * t471 + (-t174 + t632 + t813) * qJD(1) + t585;
t30 = -t469 * t634 + qJD(1) * t172 + (qJD(1) * t541 + t469 * t679) * qJD(2) + t538;
t19 = (t172 * t471 + t174 * t469 + (t275 * t471 + t469 * t627) * qJD(1)) * qJD(2) + t516;
t6 = t143 * t365 - t145 * t364 + t384 * t72 + t385 * t71 + (t238 * t471 + t239 * t469 + (t355 * t471 + t469 * t626) * qJD(1)) * qJD(2) + t516;
t1 = [t12 * t780 + (t77 + t63) * t364 / 0.2e1 + (t78 + t64) * t365 / 0.2e1 + (t13 + t18) * t384 / 0.2e1 + (t43 * (t381 - t571 + t584) + t15 * (t411 + t622 + t680) + (t43 * (-qJ(4) * t649 - t428) + t15 * t636) * t469 + (t515 * t469 - t776 * t471 - t570 + t662) * t16 + (t43 - t568 - t810 + t623 + t703 + (-t774 + (-pkin(3) - pkin(4)) * t439) * t647) * t44) * m(6) + (t100 * t332 + t837 * t209 + t204 * t212 - t536 * t99 + t367 * t470 + t368 * t468 + t890 * t440 + t889 * t439 + (t698 * t44 + (t43 * t776 + t44 * t509) * t469 + (t43 * t509 + t44 * t636) * t471) * m(6) - t881 * qJD(2)) * qJD(1) + (-(t517 + t568 - t75 - t833) * t76 + t31 * (t459 + t662) + t75 * t584 + t30 * (t277 + t622) + t76 * (t623 + t669) + (t76 * (t439 * t777 - t774) * qJD(2) + (-t76 * t466 + t499 * t75) * qJD(1)) * t471 + (-t30 * t466 + t31 * t777 * t440 + (-t75 * qJD(4) + t31 * t759) * t439 + t75 * (t440 * t759 + t764) * qJD(2) + (-t75 * rSges(5,2) + t499 * t76) * qJD(1)) * t469) * m(5) + (-(-qJD(1) * t276 + t512 + t812 - t86) * t87 + t56 * (-t276 + t662) + t86 * t621 + t55 * (t278 + t590) + t87 * (t442 + t668) + (t469 * t760 + t574 * t87) * qJD(2) + ((-t86 * rSges(4,3) + t87 * (-t434 - t765)) * t469 + (t86 * (-t362 - t434) - t87 * t466) * t471) * qJD(1)) * m(4) + (t97 * (t469 * t607 + t461 + t659) + t96 * t675 + t185 * (t438 + t661) + (t396 * t728 - t726) * qJD(2) + ((-pkin(1) - t572) * t727 + (t184 * (-rSges(3,3) - pkin(6)) + t185 * t607) * t469) * qJD(1) - (-qJD(1) * t319 - t184 - t372 - t615) * t185) * m(3) - (t14 + t12 + t17) * t385 / 0.2e1 + (((t88 - t814 + t823) * t469 + ((t891 + t902) * t471 + t856 + t882 + t901) * t471) * qJD(2) + t851) * t599 + (t827 + t830) * t648 / 0.2e1 + (((t471 * t835 - t823 + t852) * t471 + (t469 * t835 - t236 + t575 + t853 - t885) * t469) * qJD(2) + t832 + t846) * t602 - (-t828 + t829 + t831) * t647 / 0.2e1 + ((t821 - t873) * t469 + (t822 + t872) * t471) * t640 / 0.2e1; (-t809 * (-t360 + t577) * qJD(2) + t6 * t586 + (t15 * t530 + t6 * t698) * t469 + (t16 * t530 + t6 * (t145 + t356)) * t471 + (-t839 + t504 * t469 - t678 + (pkin(4) * t710 + t471 * t530 - t178 + t635) * qJD(1)) * t44 + (-t838 + t504 * t471 + (t215 * t469 + t176 - t311) * qJD(1) + t834) * t43 + (-(-t658 * t772 - t583) * qJD(2) + (qJD(1) * t588 + t757) * t469 + (qJD(1) * t698 + t758) * t471 - t816 + t836) * t32) * m(6) + (t19 * t586 + (t19 * t277 + t31 * t582) * t471 + (t19 * t275 + t30 * t582) * t469 + (t545 * t469 - t678 + (t316 + t635 + t541) * qJD(1)) * t76 + (-t311 * qJD(1) + t545 * t471 + t834) * t75 - (t469 * t76 + t471 * t75) * qJD(2) * (-t773 + t811) + (-(-t312 * t469 - t316 * t471 - t583) * qJD(2) + (t172 + t833) * t471 + (qJD(1) * t627 + t174) * t469 + t836) * t61) * m(5) + (t56 * t574 - t86 * pkin(2) * t611 + (t173 * t647 + t175 * t648 + t631) * t514 + t83 * t630 + (-t86 * t334 + t83 * t173 + (t276 * t842 + t597 * t87) * qJD(1)) * t471 + (t55 * t597 + t87 * (-t334 - t633) + t83 * t175 + (t686 * t842 + t760) * qJD(1)) * t469 - (t86 * t313 + t87 * (-t317 - t635)) * qJD(1) - (-t83 * t583 + (-t83 * t317 + t596 * t86) * t471 + (-t83 * t313 + t596 * t87) * t469) * qJD(2)) * m(4) + (0.2e1 * t180 * (t201 * t471 + t202 * t469 + (t319 * t471 - t320 * t469) * qJD(1)) + t559 * t369 + (-t96 * t469 - t97 * t471 + (-t185 * t471 + t728) * qJD(1)) * t396 - (t184 * t341 - t726) * qJD(1) - (t180 * (-t341 * t469 - t342 * t471) + t559 * t572) * qJD(2)) * m(3) + (t828 * t471 + t827 * t469 + (t821 * t469 + t822 * t471) * qJD(1)) * qJD(1) / 0.2e1 + ((-t648 * t819 + t820) * t469 + ((t469 * t818 + t825) * qJD(2) + t824) * t471) * t602 + ((-t647 * t818 - t820) * t471 + ((t471 * t819 + t825) * qJD(2) + t824) * t469) * t599 + ((t866 * t439 + t863 * t440 + t505 * t468 + t470 * t792) * qJD(2) + (t439 * t864 + t440 * t865 + t468 * t667 + t470 * t666) * qJD(1) - t859) * t769 + (t830 * qJD(1) + ((t852 * qJD(1) + t849 * t471) * t471 + (t844 * t469 + t853 * qJD(1) + (-t845 + t850) * t471) * t469) * t843) * t779 + (t829 * qJD(1) + ((t886 * qJD(1) + t845 * t471) * t471 + (t850 * t469 + t874 * qJD(1) + (-t844 + t849) * t471) * t469) * t843) * t778 + (t832 + t847) * t651 / 0.2e1 + (t831 + t848) * t650 / 0.2e1 - t868; 0.2e1 * (t15 * t778 + t16 * t779) * m(6) + 0.2e1 * (t30 * t778 + t31 * t779) * m(5) + 0.2e1 * (t55 * t778 + t56 * t779) * m(4); -m(5) * (t290 * t61 + t291 * t76 + t292 * t75) - m(6) * (t290 * t32 + t291 * t44 + t292 * t43) + 0.2e1 * ((t647 * t75 + t648 * t76 - t19) * t787 + (t43 * t647 + t44 * t648 - t6) * t786) * t440 + 0.2e1 * ((qJD(2) * t61 + t30 * t469 + t31 * t471 + t650 * t76 - t651 * t75) * t787 + (qJD(2) * t32 - t43 * t651 + t44 * t650 + t569) * t786) * t439; t859 * t769 + (t6 * (-t143 * t469 - t145 * t471) + t809 * t101 + ((-t43 * t469 + t44 * t471) * qJD(1) + t569) * t215 - t43 * (qJD(1) * t176 - t838) - t44 * (-qJD(1) * t178 - t839) + (-t469 * t72 - t471 * t71 + (-t143 * t471 + t145 * t469) * qJD(1) + t816) * t32) * m(6) + t868;];
tauc = t1(:);

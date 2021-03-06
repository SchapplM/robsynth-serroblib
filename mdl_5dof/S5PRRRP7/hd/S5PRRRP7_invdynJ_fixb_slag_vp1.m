% Calculate vector of inverse dynamics joint torques for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP7_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:08
% EndTime: 2019-12-05 16:56:24
% DurationCPUTime: 117.20s
% Computational Cost: add. (51891->1478), mult. (141804->2016), div. (0->0), fcn. (164620->10), ass. (0->501)
t895 = Icges(5,4) + Icges(6,4);
t888 = Icges(5,1) + Icges(6,1);
t887 = Icges(5,5) + Icges(6,5);
t886 = Icges(5,2) + Icges(6,2);
t885 = Icges(5,6) + Icges(6,6);
t884 = Icges(5,3) + Icges(6,3);
t558 = sin(pkin(9));
t560 = cos(pkin(9));
t566 = cos(qJ(2));
t564 = sin(qJ(2));
t775 = cos(pkin(5));
t658 = t564 * t775;
t521 = t558 * t566 + t560 * t658;
t559 = sin(pkin(5));
t563 = sin(qJ(3));
t744 = t559 * t563;
t782 = cos(qJ(3));
t449 = t521 * t782 - t560 * t744;
t562 = sin(qJ(4));
t565 = cos(qJ(4));
t657 = t566 * t775;
t600 = -t558 * t564 + t560 * t657;
t344 = -t449 * t562 - t600 * t565;
t910 = t895 * t344;
t643 = t558 * t658;
t523 = t560 * t566 - t643;
t451 = t523 * t782 + t558 * t744;
t604 = -t558 * t657 - t560 * t564;
t346 = -t451 * t562 - t604 * t565;
t909 = t895 * t346;
t677 = t559 * t782;
t525 = t563 * t775 + t564 * t677;
t742 = t559 * t566;
t454 = -t525 * t562 - t565 * t742;
t908 = t895 * t454;
t690 = t562 * t742;
t619 = -t525 * t565 + t690;
t907 = t895 * t619;
t752 = t604 * t562;
t347 = t451 * t565 - t752;
t906 = t895 * t347;
t756 = t600 * t562;
t345 = t449 * t565 - t756;
t905 = t895 * t345;
t496 = t600 * qJD(2);
t611 = -t521 * t563 - t560 * t677;
t340 = qJD(3) * t611 + t496 * t782;
t497 = t521 * qJD(2);
t154 = -qJD(4) * t345 - t340 * t562 + t497 * t565;
t586 = qJD(4) * t344 + t497 * t562;
t155 = t340 * t565 + t586;
t339 = qJD(3) * t449 + t496 * t563;
t904 = t885 * t154 + t887 * t155 + t884 * t339;
t498 = t604 * qJD(2);
t612 = -t523 * t563 + t558 * t677;
t342 = qJD(3) * t612 + t498 * t782;
t704 = qJD(2) * t560;
t499 = -qJD(2) * t643 + t566 * t704;
t156 = -qJD(4) * t347 - t342 * t562 + t499 * t565;
t585 = qJD(4) * t346 + t499 * t562;
t157 = t342 * t565 + t585;
t341 = qJD(3) * t451 + t498 * t563;
t903 = t885 * t156 + t887 * t157 + t884 * t341;
t902 = t886 * t154 + t895 * t155 + t885 * t339;
t901 = t886 * t156 + t895 * t157 + t885 * t341;
t900 = t895 * t154 + t888 * t155 + t887 * t339;
t899 = t895 * t156 + t888 * t157 + t887 * t341;
t743 = t559 * t564;
t524 = t563 * t743 - t775 * t782;
t674 = t566 * t782;
t645 = t559 * t674;
t453 = qJD(2) * t645 - qJD(3) * t524;
t703 = qJD(2) * t564;
t672 = t559 * t703;
t263 = qJD(4) * t619 - t453 * t562 + t565 * t672;
t575 = qJD(4) * t454 + t562 * t672;
t264 = t453 * t565 + t575;
t691 = t563 * t742;
t452 = qJD(2) * t691 + qJD(3) * t525;
t898 = t885 * t263 + t887 * t264 + t884 * t452;
t897 = t886 * t263 + t895 * t264 + t885 * t452;
t896 = t895 * t263 + t888 * t264 + t887 * t452;
t835 = t885 * t344 + t887 * t345 - t884 * t611;
t834 = t885 * t346 + t887 * t347 - t884 * t612;
t878 = t886 * t344 - t885 * t611 + t905;
t877 = t886 * t346 - t885 * t612 + t906;
t876 = t888 * t345 - t887 * t611 + t910;
t875 = t888 * t347 - t887 * t612 + t909;
t825 = t885 * t454 + t884 * t524 - t887 * t619;
t868 = t886 * t454 + t885 * t524 - t907;
t867 = t887 * t524 - t888 * t619 + t908;
t858 = t878 * t154 + t876 * t155 + t835 * t339 + t902 * t344 + t900 * t345 - t904 * t611;
t857 = t877 * t154 + t875 * t155 + t834 * t339 + t901 * t344 + t899 * t345 - t903 * t611;
t856 = t878 * t156 + t876 * t157 + t835 * t341 + t902 * t346 + t900 * t347 - t904 * t612;
t855 = t877 * t156 + t875 * t157 + t834 * t341 + t901 * t346 + t899 * t347 - t903 * t612;
t854 = t878 * t263 + t876 * t264 + t835 * t452 + t902 * t454 + t904 * t524 - t900 * t619;
t853 = t877 * t263 + t875 * t264 + t834 * t452 + t901 * t454 + t903 * t524 - t899 * t619;
t849 = t868 * t154 + t867 * t155 + t825 * t339 + t897 * t344 + t896 * t345 - t898 * t611;
t848 = t868 * t156 + t867 * t157 + t825 * t341 + t897 * t346 + t896 * t347 - t898 * t612;
t847 = t868 * t263 + t867 * t264 + t825 * t452 + t897 * t454 + t898 * t524 - t896 * t619;
t846 = t878 * t344 + t876 * t345 - t835 * t611;
t845 = t877 * t344 + t875 * t345 - t834 * t611;
t844 = t878 * t346 + t876 * t347 - t835 * t612;
t843 = t877 * t346 + t875 * t347 - t834 * t612;
t842 = t878 * t454 + t835 * t524 - t876 * t619;
t841 = t877 * t454 + t834 * t524 - t875 * t619;
t840 = t868 * t344 + t867 * t345 - t825 * t611;
t839 = t868 * t346 + t867 * t347 - t825 * t612;
t561 = -qJ(5) - pkin(8);
t890 = rSges(6,3) - t561;
t836 = t868 * t454 + t825 * t524 - t867 * t619;
t893 = rSges(6,1) + pkin(4);
t380 = Icges(4,5) * t525 - Icges(4,6) * t524 - Icges(4,3) * t742;
t769 = Icges(4,4) * t525;
t381 = -Icges(4,2) * t524 - Icges(4,6) * t742 + t769;
t516 = Icges(4,4) * t524;
t382 = Icges(4,1) * t525 - Icges(4,5) * t742 - t516;
t112 = -t380 * t600 + t381 * t611 + t382 * t449;
t693 = qJDD(2) * t559;
t548 = t558 * t693;
t351 = qJD(3) * t499 - qJDD(3) * t604 + t548;
t666 = t560 * t693;
t352 = qJD(3) * t497 - qJDD(3) * t600 - t666;
t158 = Icges(4,5) * t340 - Icges(4,6) * t339 + Icges(4,3) * t497;
t160 = Icges(4,4) * t340 - Icges(4,2) * t339 + Icges(4,6) * t497;
t162 = Icges(4,1) * t340 - Icges(4,4) * t339 + Icges(4,5) * t497;
t238 = Icges(4,5) * t449 + Icges(4,6) * t611 - Icges(4,3) * t600;
t771 = Icges(4,4) * t449;
t240 = Icges(4,2) * t611 - Icges(4,6) * t600 + t771;
t435 = Icges(4,4) * t611;
t242 = Icges(4,1) * t449 - Icges(4,5) * t600 + t435;
t45 = -t158 * t600 + t160 * t611 + t162 * t449 + t238 * t497 - t240 * t339 + t242 * t340;
t159 = Icges(4,5) * t342 - Icges(4,6) * t341 + Icges(4,3) * t499;
t161 = Icges(4,4) * t342 - Icges(4,2) * t341 + Icges(4,6) * t499;
t163 = Icges(4,1) * t342 - Icges(4,4) * t341 + Icges(4,5) * t499;
t239 = Icges(4,5) * t451 + Icges(4,6) * t612 - Icges(4,3) * t604;
t770 = Icges(4,4) * t451;
t241 = Icges(4,2) * t612 - Icges(4,6) * t604 + t770;
t436 = Icges(4,4) * t612;
t243 = Icges(4,1) * t451 - Icges(4,5) * t604 + t436;
t46 = -t159 * t600 + t161 * t611 + t163 * t449 + t239 * t497 - t241 * t339 + t243 * t340;
t705 = qJD(2) * t559;
t549 = t558 * t705;
t460 = -qJD(3) * t604 + t549;
t673 = t559 * t704;
t461 = -qJD(3) * t600 - t673;
t556 = qJDD(2) * t775;
t700 = qJD(3) * t564;
t671 = t559 * t700;
t464 = qJD(2) * t671 - qJDD(3) * t742 + t556;
t557 = qJD(2) * t775;
t533 = -qJD(3) * t742 + t557;
t266 = Icges(4,5) * t453 - Icges(4,6) * t452 + Icges(4,3) * t672;
t267 = Icges(4,4) * t453 - Icges(4,2) * t452 + Icges(4,6) * t672;
t268 = Icges(4,1) * t453 - Icges(4,4) * t452 + Icges(4,5) * t672;
t64 = -t266 * t600 + t267 * t611 + t268 * t449 - t339 * t381 + t340 * t382 + t380 * t497;
t152 = qJD(4) * t341 - qJDD(4) * t612 + t351;
t153 = qJD(4) * t339 - qJDD(4) * t611 + t352;
t265 = qJD(4) * t452 + qJDD(4) * t524 + t464;
t316 = -qJD(4) * t612 + t460;
t317 = -qJD(4) * t611 + t461;
t434 = qJD(4) * t524 + t533;
t861 = t845 * t152 + t846 * t153 + t840 * t265 + t857 * t316 + t317 * t858 + t849 * t434;
t95 = -t238 * t600 + t240 * t611 + t242 * t449;
t96 = -t239 * t600 + t241 * t611 + t243 * t449;
t892 = t112 * t464 + t351 * t96 + t352 * t95 + t461 * t45 + t46 * t460 + t533 * t64 + t861;
t113 = -t380 * t604 + t381 * t612 + t382 * t451;
t47 = -t158 * t604 + t160 * t612 + t162 * t451 + t238 * t499 - t240 * t341 + t242 * t342;
t48 = -t159 * t604 + t161 * t612 + t163 * t451 + t239 * t499 - t241 * t341 + t243 * t342;
t65 = -t266 * t604 + t267 * t612 + t268 * t451 - t341 * t381 + t342 * t382 + t380 * t499;
t860 = t152 * t843 + t153 * t844 + t265 * t839 + t316 * t855 + t317 * t856 + t434 * t848;
t97 = -t238 * t604 + t240 * t612 + t242 * t451;
t98 = -t239 * t604 + t241 * t612 + t243 * t451;
t891 = t113 * t464 + t351 * t98 + t352 * t97 + t460 * t48 + t461 * t47 + t533 * t65 + t860;
t880 = -pkin(8) + t890;
t106 = -t238 * t742 - t524 * t240 + t525 * t242;
t107 = -t239 * t742 - t524 * t241 + t525 * t243;
t130 = -t380 * t742 - t524 * t381 + t525 * t382;
t51 = -t524 * t160 + t525 * t162 - t452 * t240 + t453 * t242 + (-t158 * t566 + t238 * t703) * t559;
t52 = -t524 * t161 + t525 * t163 - t452 * t241 + t453 * t243 + (-t159 * t566 + t239 * t703) * t559;
t72 = -t524 * t267 + t525 * t268 - t452 * t381 + t453 * t382 + (-t266 * t566 + t380 * t703) * t559;
t859 = t152 * t841 + t153 * t842 + t265 * t836 + t316 * t853 + t317 * t854 + t434 * t847;
t889 = t106 * t352 + t107 * t351 + t130 * t464 + t460 * t52 + t461 * t51 + t533 * t72 + t859;
t883 = t885 * t562 - t887 * t565;
t882 = t886 * t562 - t895 * t565;
t881 = t895 * t562 - t888 * t565;
t852 = t316 * t845 + t317 * t846 + t434 * t840;
t851 = t316 * t843 + t317 * t844 + t434 * t839;
t850 = t316 * t841 + t317 * t842 + t434 * t836;
t513 = (pkin(2) * t566 + pkin(7) * t564) * t705;
t297 = t453 * pkin(3) + t452 * pkin(8);
t696 = qJD(5) * t524;
t555 = pkin(4) * t565 + pkin(3);
t780 = -pkin(3) + t555;
t739 = rSges(6,1) * t264 + rSges(6,2) * t263 + pkin(4) * t575 + t452 * t880 + t453 * t780 + t696;
t687 = -t297 - t739;
t879 = t559 * (-t513 + t687);
t676 = t562 * t782;
t392 = t521 * t565 - t600 * t676;
t675 = t565 * t782;
t753 = t521 * t562;
t393 = t600 * t675 + t753;
t755 = t600 * t563;
t874 = t392 * t885 + t393 * t887 + t755 * t884;
t394 = t523 * t565 - t604 * t676;
t749 = t523 * t562;
t395 = t604 * t675 + t749;
t751 = t604 * t563;
t873 = t394 * t885 + t395 * t887 + t751 * t884;
t872 = t392 * t886 + t393 * t895 + t755 * t885;
t871 = t394 * t886 + t395 * t895 + t751 * t885;
t870 = t392 * t895 + t393 * t888 + t755 * t887;
t869 = t394 * t895 + t395 * t888 + t751 * t887;
t470 = (-t562 * t674 + t564 * t565) * t559;
t740 = t562 * t564;
t471 = (t565 * t674 + t740) * t559;
t866 = t470 * t885 + t471 * t887 + t691 * t884;
t865 = t470 * t886 + t471 * t895 + t691 * t885;
t864 = t470 * t895 + t471 * t888 + t691 * t887;
t863 = t555 * t782 - t561 * t563;
t176 = pkin(3) * t342 + pkin(8) * t341;
t290 = pkin(3) * t451 - pkin(8) * t612;
t427 = t525 * pkin(3) + t524 * pkin(8);
t399 = pkin(2) * t498 + pkin(7) * t499;
t419 = pkin(2) * t523 - pkin(7) * t604;
t530 = (pkin(2) * t564 - pkin(7) * t566) * t559;
t746 = t558 * t559;
t583 = t399 * t557 + t419 * t556 + (-qJD(2) * t513 - qJDD(2) * t530) * t746;
t572 = t533 * t176 + t464 * t290 - t297 * t460 - t351 * t427 + t583;
t726 = -rSges(6,1) * t619 + rSges(6,2) * t454 - pkin(4) * t690 + t524 * t880 + t525 * t780;
t736 = rSges(6,1) * t347 + rSges(6,2) * t346 - pkin(4) * t752 + t451 * t780 - t612 * t880;
t697 = qJD(5) * t612;
t777 = rSges(6,1) * t157 + rSges(6,2) * t156 + pkin(4) * t585 + t341 * t880 + t342 * t780 - t697;
t11 = qJD(5) * t339 - qJDD(5) * t611 - t152 * t726 + t265 * t736 - t316 * t739 + t434 * t777 + t572;
t175 = pkin(3) * t340 + pkin(8) * t339;
t288 = pkin(3) * t449 - pkin(8) * t611;
t417 = pkin(2) * t521 - pkin(7) * t600;
t655 = t775 * t417;
t745 = t559 * t560;
t601 = -t530 * t745 - t655;
t398 = pkin(2) * t496 + pkin(7) * t497;
t656 = t775 * t398;
t571 = (-t513 * t745 - t656) * qJD(2) + t601 * qJDD(2);
t568 = -t533 * t175 - t464 * t288 + t461 * t297 + t352 * t427 + t571;
t737 = rSges(6,1) * t345 + rSges(6,2) * t344 - pkin(4) * t756 + t449 * t780 - t611 * t880;
t698 = qJD(5) * t611;
t778 = rSges(6,1) * t155 + rSges(6,2) * t154 + pkin(4) * t586 + t339 * t880 + t340 * t780 - t698;
t12 = qJD(5) * t341 - qJDD(5) * t612 + t153 * t726 - t265 * t737 + t317 * t739 - t434 * t778 + t568;
t837 = t12 * t737;
t862 = t11 * t736 - t837;
t150 = rSges(5,1) * t347 + rSges(5,2) * t346 - rSges(5,3) * t612;
t735 = -t150 - t290;
t833 = t393 * rSges(6,1) + t392 * rSges(6,2) + rSges(6,3) * t755 + pkin(4) * t753 + t600 * t863;
t832 = t395 * rSges(6,1) + t394 * rSges(6,2) + rSges(6,3) * t751 + pkin(4) * t749 + t604 * t863;
t831 = t449 * t885 - t611 * t882;
t830 = t451 * t885 - t612 * t882;
t829 = t449 * t887 - t611 * t881;
t828 = t451 * t887 - t612 * t881;
t761 = t611 * t565;
t762 = t611 * t562;
t827 = rSges(6,1) * t761 - rSges(6,2) * t762 + t449 * t890 + t555 * t611;
t758 = t612 * t565;
t759 = t612 * t562;
t826 = rSges(6,1) * t758 - rSges(6,2) * t759 + t451 * t890 + t555 * t612;
t255 = -rSges(5,1) * t619 + rSges(5,2) * t454 + rSges(5,3) * t524;
t725 = -t255 - t427;
t824 = t524 * t882 + t525 * t885;
t823 = t524 * t881 + t525 * t887;
t747 = t524 * t565;
t748 = t524 * t562;
t822 = -rSges(6,1) * t747 + rSges(6,2) * t748 - t524 * t555 + t525 * t890;
t821 = t559 * pkin(4) * t740 + t471 * rSges(6,1) + t470 * rSges(6,2) + t555 * t645 + t691 * t890;
t820 = (t524 * t883 + t525 * t884 + t562 * t868 - t565 * t867) * t434 + (t449 * t884 + t562 * t878 - t565 * t876 - t611 * t883) * t317 + (t451 * t884 + t562 * t877 - t565 * t875 - t612 * t883) * t316;
t819 = (t619 * t886 + t867 + t908) * t434 + (-t345 * t886 + t876 + t910) * t317 + (-t347 * t886 + t875 + t909) * t316;
t818 = (t454 * t888 - t868 + t907) * t434 + (t344 * t888 - t878 - t905) * t317 + (t346 * t888 - t877 - t906) * t316;
t817 = (t454 * t887 + t619 * t885) * t434 + (t344 * t887 - t345 * t885) * t317 + (t346 * t887 - t347 * t885) * t316;
t287 = pkin(3) * t611 + pkin(8) * t449;
t289 = pkin(3) * t612 + pkin(8) * t451;
t816 = -t175 * t604 - t460 * t287 + t499 * t288 + t289 * t461;
t667 = t417 * t549 + t419 * t673 + qJD(1);
t613 = t460 * t288 - t290 * t461 + t667;
t38 = t316 * t737 - t317 * t736 + t613 + t696;
t624 = t398 * t549 + t399 * t673 + t417 * t548 + t419 * t666 + qJDD(1);
t576 = t460 * t175 - t176 * t461 + t351 * t288 - t290 * t352 + t624;
t7 = qJD(5) * t452 + qJDD(5) * t524 + t152 * t737 - t153 * t736 + t316 * t778 - t317 * t777 + t576;
t815 = -t38 * t777 - t7 * t736;
t616 = Icges(4,5) * t782 - Icges(4,6) * t563;
t814 = t460 * (-Icges(4,3) * t523 - t241 * t563 + t243 * t782 - t604 * t616) + t461 * (-Icges(4,3) * t521 - t240 * t563 + t242 * t782 - t600 * t616) + t533 * (-t381 * t563 + t382 * t782 - (Icges(4,3) * t564 + t566 * t616) * t559);
t813 = t152 / 0.2e1;
t812 = t153 / 0.2e1;
t811 = t265 / 0.2e1;
t810 = -t316 / 0.2e1;
t809 = t316 / 0.2e1;
t808 = -t317 / 0.2e1;
t807 = t317 / 0.2e1;
t804 = t351 / 0.2e1;
t803 = t352 / 0.2e1;
t800 = -t434 / 0.2e1;
t799 = t434 / 0.2e1;
t795 = -t460 / 0.2e1;
t794 = t460 / 0.2e1;
t793 = -t461 / 0.2e1;
t792 = t461 / 0.2e1;
t791 = t464 / 0.2e1;
t784 = -t533 / 0.2e1;
t783 = t533 / 0.2e1;
t94 = rSges(5,1) * t157 + rSges(5,2) * t156 + rSges(5,3) * t341;
t776 = -t176 - t94;
t774 = Icges(3,4) * t521;
t773 = Icges(3,4) * t523;
t772 = Icges(3,4) * t564;
t123 = rSges(5,1) * t264 + rSges(5,2) * t263 + rSges(5,3) * t452;
t738 = -t123 - t297;
t734 = -t287 + t827;
t733 = -t289 + t826;
t357 = t775 * t399;
t731 = t775 * t176 + t357;
t730 = -t345 * rSges(6,2) + t344 * t893;
t729 = -t347 * rSges(6,2) + t346 * t893;
t681 = t600 * t782;
t707 = -pkin(3) * t681 - pkin(8) * t755;
t728 = t707 + t833;
t679 = t604 * t782;
t706 = -pkin(3) * t679 - pkin(8) * t751;
t727 = t706 + t832;
t426 = -t524 * pkin(3) + pkin(8) * t525;
t724 = -t426 + t822;
t723 = t288 * t742 - t427 * t600;
t407 = t775 * t419;
t722 = t775 * t290 + t407;
t721 = t619 * rSges(6,2) + t454 * t893;
t491 = pkin(3) * t645 + pkin(8) * t691;
t720 = -t491 + t821;
t719 = t398 * t746 + t399 * t745;
t360 = Icges(3,2) * t600 - Icges(3,6) * t745 + t774;
t412 = Icges(3,1) * t600 - t774;
t718 = -t360 + t412;
t361 = Icges(3,2) * t604 + Icges(3,6) * t746 + t773;
t413 = Icges(3,1) * t604 - t773;
t717 = -t361 + t413;
t509 = Icges(3,4) * t600;
t362 = Icges(3,1) * t521 - Icges(3,5) * t745 + t509;
t410 = -Icges(3,2) * t521 + t509;
t716 = -t362 - t410;
t510 = Icges(3,4) * t604;
t363 = Icges(3,1) * t523 + Icges(3,5) * t746 + t510;
t411 = -Icges(3,2) * t523 + t510;
t715 = -t363 - t411;
t416 = pkin(2) * t600 + pkin(7) * t521;
t418 = pkin(2) * t604 + pkin(7) * t523;
t714 = t416 * t549 + t418 * t673;
t713 = t417 * t746 + t419 * t745;
t475 = Icges(3,6) * t775 + (Icges(3,2) * t566 + t772) * t559;
t528 = (Icges(3,1) * t566 - t772) * t559;
t709 = -t475 + t528;
t550 = Icges(3,4) * t742;
t476 = Icges(3,1) * t743 + Icges(3,5) * t775 + t550;
t527 = -Icges(3,2) * t743 + t550;
t708 = -t476 - t527;
t531 = pkin(2) * t742 + pkin(7) * t743;
t702 = qJD(3) * t521;
t701 = qJD(3) * t523;
t699 = qJD(4) * t563;
t695 = qJD(5) * t563;
t692 = -t176 - t777;
t686 = -t513 + t738;
t685 = -t290 - t736;
t684 = t175 * t742 - t297 * t600 + t497 * t427;
t683 = -t427 - t726;
t682 = -t530 + t725;
t209 = t393 * rSges(5,1) + t392 * rSges(5,2) + rSges(5,3) * t755;
t211 = t395 * rSges(5,1) + t394 * rSges(5,2) + rSges(5,3) * t751;
t326 = t471 * rSges(5,1) + t470 * rSges(5,2) + rSges(5,3) * t691;
t665 = t705 / 0.2e1;
t661 = t38 * t737;
t620 = t419 * t557 - t530 * t549;
t589 = t533 * t290 - t427 * t460 + t620;
t43 = -t316 * t726 + t434 * t736 + t589 - t698;
t660 = t43 * t736;
t590 = t601 * qJD(2);
t574 = -t533 * t288 + t461 * t427 + t590;
t44 = t317 * t726 - t434 * t737 + t574 - t697;
t659 = t44 * t726;
t269 = rSges(4,1) * t453 - rSges(4,2) * t452 + rSges(4,3) * t672;
t654 = t559 * (-t269 - t513);
t383 = t525 * rSges(4,1) - t524 * rSges(4,2) - rSges(4,3) * t742;
t653 = t559 * (-t383 - t530);
t651 = t533 * t289 - t426 * t460;
t650 = -t287 * t533 + t461 * t426;
t648 = t175 * t746 + t176 * t745 + t719;
t647 = -t530 + t683;
t646 = t288 * t746 + t290 * t745 + t713;
t148 = rSges(5,1) * t345 + rSges(5,2) * t344 - rSges(5,3) * t611;
t92 = rSges(5,1) * t155 + rSges(5,2) * t154 + rSges(5,3) * t339;
t13 = t148 * t152 - t150 * t153 + t316 * t92 - t317 * t94 + t576;
t53 = t148 * t316 - t150 * t317 + t613;
t642 = t13 * t148 + t53 * t92;
t386 = Icges(3,5) * t496 - Icges(3,6) * t497;
t387 = Icges(3,5) * t498 - Icges(3,6) * t499;
t388 = Icges(3,4) * t496 - Icges(3,2) * t497;
t389 = Icges(3,4) * t498 - Icges(3,2) * t499;
t390 = Icges(3,1) * t496 - Icges(3,4) * t497;
t391 = Icges(3,1) * t498 - Icges(3,4) * t499;
t634 = -(-t360 * t499 + t362 * t498 + t386 * t746 + t388 * t604 + t390 * t523) * t560 + (-t361 * t499 + t363 * t498 + t387 * t746 + t389 * t604 + t391 * t523) * t558;
t628 = -t361 * t564 + t363 * t566;
t629 = -t360 * t564 + t362 * t566;
t633 = -(t775 * t386 + (qJD(2) * t629 + t388 * t566 + t390 * t564) * t559) * t560 + (t775 * t387 + (qJD(2) * t628 + t389 * t566 + t391 * t564) * t559) * t558;
t358 = Icges(3,5) * t521 + Icges(3,6) * t600 - Icges(3,3) * t745;
t359 = Icges(3,5) * t523 + Icges(3,6) * t604 + Icges(3,3) * t746;
t632 = -(-t358 * t745 + t360 * t600 + t362 * t521) * t560 + (-t359 * t745 + t361 * t600 + t363 * t521) * t558;
t631 = -(t358 * t746 + t360 * t604 + t362 * t523) * t560 + (t359 * t746 + t361 * t604 + t363 * t523) * t558;
t630 = -(t775 * t358 + (t360 * t566 + t362 * t564) * t559) * t560 + (t775 * t359 + (t361 * t566 + t363 * t564) * t559) * t558;
t366 = rSges(3,1) * t521 + rSges(3,2) * t600 - rSges(3,3) * t745;
t367 = rSges(3,1) * t523 + rSges(3,2) * t604 + rSges(3,3) * t746;
t627 = t366 * t558 + t367 * t560;
t396 = rSges(3,1) * t496 - rSges(3,2) * t497;
t397 = rSges(3,1) * t498 - rSges(3,2) * t499;
t626 = t396 * t558 + t397 * t560;
t625 = -t475 * t564 + t476 * t566;
t306 = rSges(4,1) * t681 - rSges(4,2) * t755 + t521 * rSges(4,3);
t307 = rSges(4,1) * t679 - rSges(4,2) * t751 + t523 * rSges(4,3);
t228 = rSges(5,1) * t761 - rSges(5,2) * t762 + t449 * rSges(5,3);
t230 = rSges(5,1) * t758 - rSges(5,2) * t759 + t451 * rSges(5,3);
t315 = -rSges(5,1) * t747 + rSges(5,2) * t748 + t525 * rSges(5,3);
t621 = t418 * t557 - t531 * t549;
t618 = Icges(4,1) * t782 - Icges(4,4) * t563;
t617 = Icges(4,4) * t782 - Icges(4,2) * t563;
t529 = (rSges(3,1) * t566 - rSges(3,2) * t564) * t559;
t459 = rSges(4,1) * t645 - rSges(4,2) * t691 + rSges(4,3) * t743;
t614 = t559 * t627;
t526 = (Icges(3,5) * t566 - Icges(3,6) * t564) * t559;
t607 = t38 * t778 + t7 * t737;
t494 = t775 * rSges(3,3) + (rSges(3,1) * t564 + rSges(3,2) * t566) * t559;
t606 = t367 * t775 - t494 * t746;
t508 = qJD(2) * t529;
t605 = t397 * t775 - t508 * t746;
t603 = -t366 * t775 - t494 * t745;
t602 = -t396 * t775 - t508 * t745;
t597 = (Icges(4,5) * t611 - Icges(4,6) * t449) * t461 + (Icges(4,5) * t612 - Icges(4,6) * t451) * t460 + (-Icges(4,5) * t524 - Icges(4,6) * t525) * t533;
t596 = t12 * t726 + t44 * t739;
t595 = -t43 * t726 + t661;
t594 = -t44 * t737 + t660;
t593 = -t38 * t736 + t659;
t592 = -t175 * t775 - t656;
t591 = -t288 * t775 - t655;
t588 = -t288 * t671 + t427 * t702 + t461 * t491 + t533 * t707;
t587 = -t416 * t557 - t531 * t673;
t584 = t288 * t701 - t290 * t702 - t460 * t707 + t461 * t706 + t714;
t578 = (Icges(4,1) * t612 - t241 - t770) * t460 + (Icges(4,1) * t611 - t240 - t771) * t461 + (-Icges(4,1) * t524 - t381 - t769) * t533;
t577 = (Icges(4,2) * t451 - t243 - t436) * t460 + (Icges(4,2) * t449 - t242 - t435) * t461 + (Icges(4,2) * t525 - t382 + t516) * t533;
t573 = t290 * t671 - t427 * t701 - t460 * t491 - t533 * t706 + t621;
t503 = qJD(2) * t528;
t502 = (Icges(3,4) * t566 - Icges(3,2) * t564) * t705;
t501 = qJD(2) * t526;
t500 = (t566 * t699 + t700) * t559;
t474 = Icges(3,3) * t775 + (Icges(3,5) * t564 + Icges(3,6) * t566) * t559;
t458 = (Icges(4,5) * t564 + t566 * t618) * t559;
t457 = (Icges(4,6) * t564 + t566 * t617) * t559;
t423 = -rSges(4,1) * t524 - rSges(4,2) * t525;
t415 = rSges(3,1) * t604 - rSges(3,2) * t523;
t414 = rSges(3,1) * t600 - rSges(3,2) * t521;
t409 = Icges(3,5) * t604 - Icges(3,6) * t523;
t408 = Icges(3,5) * t600 - Icges(3,6) * t521;
t406 = t604 * t699 + t701;
t405 = t600 * t699 + t702;
t305 = Icges(4,5) * t523 + t604 * t618;
t304 = Icges(4,5) * t521 + t600 * t618;
t303 = Icges(4,6) * t523 + t604 * t617;
t302 = Icges(4,6) * t521 + t600 * t617;
t299 = rSges(5,1) * t454 + rSges(5,2) * t619;
t286 = rSges(4,1) * t612 - rSges(4,2) * t451;
t285 = rSges(4,1) * t611 - rSges(4,2) * t449;
t270 = t290 * t672;
t259 = t603 * qJD(2);
t258 = t606 * qJD(2);
t253 = t604 * t288;
t245 = rSges(4,1) * t451 + rSges(4,2) * t612 - rSges(4,3) * t604;
t244 = rSges(4,1) * t449 + rSges(4,2) * t611 - rSges(4,3) * t600;
t214 = t474 * t746 + t475 * t604 + t476 * t523;
t213 = -t474 * t745 + t475 * t600 + t476 * t521;
t193 = rSges(5,1) * t346 - rSges(5,2) * t347;
t191 = rSges(5,1) * t344 - rSges(5,2) * t345;
t165 = rSges(4,1) * t342 - rSges(4,2) * t341 + rSges(4,3) * t499;
t164 = rSges(4,1) * t340 - rSges(4,2) * t339 + rSges(4,3) * t497;
t125 = -t475 * t499 + t476 * t498 + t501 * t746 + t502 * t604 + t503 * t523;
t124 = -t475 * t497 + t476 * t496 - t501 * t745 + t502 * t600 + t503 * t521;
t114 = qJDD(1) + (qJD(2) * t626 + qJDD(2) * t627) * t559;
t111 = -t533 * t244 + t461 * t383 + t590;
t110 = t245 * t533 - t383 * t460 + t620;
t101 = t244 * t460 - t245 * t461 + t667;
t67 = -t434 * t148 + t317 * t255 + t574;
t66 = t150 * t434 - t255 * t316 + t589;
t55 = -t533 * t164 - t464 * t244 + t461 * t269 + t352 * t383 + t571;
t54 = t165 * t533 + t245 * t464 - t269 * t460 - t351 * t383 + t583;
t50 = t164 * t460 - t165 * t461 + t244 * t351 - t245 * t352 + t624;
t49 = t106 * t461 + t107 * t460 + t533 * t130;
t42 = t113 * t533 + t460 * t98 + t461 * t97;
t41 = t112 * t533 + t460 * t96 + t461 * t95;
t31 = t317 * t123 - t265 * t148 + t153 * t255 - t434 * t92 + t568;
t30 = -t123 * t316 + t150 * t265 - t152 * t255 + t434 * t94 + t572;
t1 = [m(2) * qJDD(1) + (-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3) + m(3) * t114 + m(4) * t50 + m(5) * t13 + m(6) * t7; (-g(1) * (t418 + t211 - t706) - g(2) * (t416 + t209 - t707) - g(3) * (t326 + t531 + t491) + t31 * (-t148 * t775 + t591) + t67 * (-t775 * t92 + t592) + t30 * (t150 * t775 + t722) + t66 * (t775 * t94 + t731) + t13 * t646 + t53 * t648 + ((t13 * t150 + t31 * t682 + t53 * t94 + t67 * t686) * t560 + (t30 * t682 + t66 * t686 + t642) * t558) * t559 - t67 * (-t500 * t148 - t434 * t209 + t405 * t255 + t317 * t326 + (-t416 * t775 - t531 * t745) * qJD(2) + t588) - t66 * (t150 * t500 + t211 * t434 - t255 * t406 - t316 * t326 + t573) - t53 * (t148 * t406 - t150 * t405 + t209 * t316 - t211 * t317 + t584)) * m(5) + (t125 * t557 + t214 * t556 + (qJD(2) * t634 + qJDD(2) * t631) * t559 + t891) * t746 / 0.2e1 - (t632 * t693 + t124 * t557 + t213 * t556 + (t124 * t775 + 0.2e1 * t559 * (-(-t360 * t497 + t362 * t496 - t386 * t745 + t388 * t600 + t390 * t521) * t560 + (-t361 * t497 + t363 * t496 - t387 * t745 + t389 * t600 + t391 * t521) * t558)) * qJD(2) + t892) * t745 / 0.2e1 + (-t558 * ((t409 * t746 + t523 * t717 - t604 * t715) * t746 - (t408 * t746 + t523 * t718 - t604 * t716) * t745 + (t523 * t709 + t526 * t746 - t604 * t708) * t775) / 0.2e1 + t560 * ((-t409 * t745 + t521 * t717 - t600 * t715) * t746 - (-t408 * t745 + t521 * t718 - t600 * t716) * t745 + (t521 * t709 - t526 * t745 - t600 * t708) * t775) / 0.2e1) * t559 * qJD(2) ^ 2 + (t836 * t500 + (t454 * t865 + t470 * t868 + t471 * t867 + t524 * t866 - t619 * t864 + t691 * t825) * t434 + t841 * t406 + t842 * t405 + (t454 * t872 + t470 * t878 + t471 * t876 + t524 * t874 - t619 * t870 + t691 * t835) * t317 + (t454 * t871 + t470 * t877 + t471 * t875 + t524 * t873 - t619 * t869 + t691 * t834) * t316) * t800 + (t836 * t775 + (t558 * t841 - t560 * t842) * t559) * t811 + (t839 * t500 + (t346 * t865 + t347 * t864 + t394 * t868 + t395 * t867 - t612 * t866 + t751 * t825) * t434 + t843 * t406 + t844 * t405 + (t346 * t872 + t347 * t870 + t394 * t878 + t395 * t876 - t612 * t874 + t751 * t835) * t317 + (t346 * t871 + t347 * t869 + t394 * t877 + t395 * t875 - t612 * t873 + t751 * t834) * t316) * t810 + (t839 * t775 + (t558 * t843 - t560 * t844) * t559) * t813 + (t840 * t500 + (t344 * t865 + t345 * t864 + t392 * t868 + t393 * t867 - t611 * t866 + t755 * t825) * t434 + t845 * t406 + t846 * t405 + (t344 * t872 + t345 * t870 + t392 * t878 + t393 * t876 - t611 * t874 + t755 * t835) * t317 + (t344 * t871 + t345 * t869 + t392 * t877 + t393 * t875 - t611 * t873 + t755 * t834) * t316) * t808 + (t840 * t775 + (t558 * t845 - t560 * t846) * t559) * t812 + (t630 * t556 + t633 * t557) * t559 / 0.2e1 + (-(t259 * (-t414 * t775 - t529 * t745) + t258 * (t415 * t775 - t529 * t746)) * qJD(2) + (qJD(2) * t602 + qJDD(2) * t603) * t603 + t259 * t602 + (qJD(2) * t605 + qJDD(2) * t606) * t606 + t258 * t605 + t114 * t614 - g(1) * t415 - g(2) * t414 - g(3) * t529 + (-(t414 * t558 + t415 * t560) * t705 + t626 * t559) * (qJD(2) * t614 + qJD(1))) * m(3) + (t55 * (-t244 * t775 + t560 * t653 - t655) + t54 * (t245 * t775 + t558 * t653 + t407) + t50 * ((t244 * t558 + t245 * t560) * t559 + t713) - g(1) * (t307 + t418) - g(2) * (t306 + t416) - g(3) * (t459 + t531) + (-t164 * t775 + t560 * t654 - t656 + t533 * t306 - t461 * t459 - t587 - (-t244 * t743 + t383 * t521) * qJD(3)) * t111 + (t165 * t775 + t558 * t654 + t357 - t307 * t533 + t459 * t460 - t621 - (t245 * t743 - t383 * t523) * qJD(3)) * t110 + ((t164 * t558 + t165 * t560) * t559 + t719 - t306 * t460 + t307 * t461 - t714 - (t244 * t523 - t245 * t521) * qJD(3)) * t101) * m(4) + (t12 * t591 + t11 * t722 + t7 * t646 + ((t12 * t647 - t815) * t560 + (t11 * t647 + t607) * t558) * t559 + t862 * t775 - g(1) * (t418 + t832) - g(2) * (t416 + t833) - g(3) * (t531 + t821) - t594 * t500 - t595 * t406 - t593 * t405 + (-t720 * t317 + t728 * t434 + t560 * t879 - t604 * t695 - t778 * t775 - t587 - t588 + t592) * t44 + (t720 * t316 - t727 * t434 + t558 * t879 - t600 * t695 + t777 * t775 - t573 + t731) * t43 + (-qJD(5) * t691 - t728 * t316 + t727 * t317 - t584 + t648) * t38) * m(6) - t850 * t500 / 0.2e1 - t851 * t406 / 0.2e1 - t852 * t405 / 0.2e1 + (t847 * t775 + (t558 * t853 - t560 * t854) * t559) * t799 + (t848 * t775 + (t558 * t855 - t560 * t856) * t559) * t809 + (t849 * t775 + (t558 * t857 - t560 * t858) * t559) * t807 - (t213 * t775 + t559 * t632) * t666 / 0.2e1 + (t214 * t775 + t559 * t631) * t548 / 0.2e1 - t42 * t701 / 0.2e1 - t41 * t702 / 0.2e1 + t556 * t775 * (t775 * t474 + (t475 * t566 + t476 * t564) * t559) - t49 * t671 / 0.2e1 + t557 * t775 * (t775 * t501 + (qJD(2) * t625 + t502 * t566 + t503 * t564) * t559) + ((qJD(2) * t633 + qJDD(2) * t630) * t559 + t889) * t775 / 0.2e1 + ((-t524 * t303 + t525 * t305) * t460 + (-t524 * t302 + t525 * t304) * t461 + (-t524 * t457 + t525 * t458) * t533 + (t106 * t521 + t107 * t523) * qJD(3) + ((t130 * qJD(3) + t238 * t461 + t239 * t460 + t380 * t533) * t564 + t814 * t566) * t559) * t784 + ((t523 * t239 + t303 * t612 + t451 * t305) * t460 + (t523 * t238 + t302 * t612 + t451 * t304) * t461 + (t523 * t380 + t451 * t458 + t457 * t612) * t533 + (t113 * t743 + t521 * t97 + t523 * t98) * qJD(3) + t814 * t604) * t795 + ((t521 * t239 + t303 * t611 + t449 * t305) * t460 + (t521 * t238 + t302 * t611 + t449 * t304) * t461 + (t521 * t380 + t449 * t458 + t457 * t611) * t533 + (t112 * t743 + t521 * t95 + t523 * t96) * qJD(3) + t814 * t600) * t793 + t558 * (t125 * t775 + t559 * t634) * t665 + (t112 * t775 + (t558 * t96 - t560 * t95) * t559) * t803 + (t113 * t775 + (t558 * t98 - t560 * t97) * t559) * t804 + (t130 * t775 + (-t106 * t560 + t107 * t558) * t559) * t791 + (t64 * t775 + (-t45 * t560 + t46 * t558) * t559) * t792 + (t65 * t775 + (-t47 * t560 + t48 * t558) * t559) * t794 - (((t411 * t566 + t413 * t564 + t628) * t558 - (t410 * t566 + t412 * t564 + t629) * t560) * t559 * t705 + ((-t408 * t560 + t409 * t558 + t527 * t566 + t528 * t564 + t625) * t705 + t526 * t557) * t775) * t557 / 0.2e1 + (t72 * t775 + (-t51 * t560 + t52 * t558) * t559) * t783; (t42 + t851) * t499 / 0.2e1 + (t41 + t852) * t497 / 0.2e1 - t891 * t604 / 0.2e1 - t892 * t600 / 0.2e1 + (-t600 * t846 - t604 * t845 - t742 * t840) * t812 + ((-t566 * t847 + t703 * t836) * t559 - t853 * t604 - t854 * t600 + t841 * t499 + t842 * t497) * t799 + (-t600 * t842 - t604 * t841 - t742 * t836) * t811 + (-t600 * t844 - t604 * t843 - t742 * t839) * t813 + ((-t566 * t848 + t703 * t839) * t559 - t855 * t604 - t856 * t600 + t843 * t499 + t844 * t497) * t809 + ((-t566 * t849 + t703 * t840) * t559 - t857 * t604 - t858 * t600 + t845 * t499 + t846 * t497) * t807 + (t106 * t497 + t107 * t499 - t51 * t600 - t52 * t604 + (t130 * t703 - t566 * t72) * t559) * t783 + (t50 * (-t244 * t604 + t245 * t600) + (t110 * t604 - t111 * t600) * t269 + (-t110 * t499 + t111 * t497 + t54 * t604 - t55 * t600) * t383 + ((t110 * t245 - t111 * t244) * t703 + (-t110 * t165 + t111 * t164 + t244 * t55 - t245 * t54) * t566) * t559 - t111 * (-t285 * t533 + t423 * t461) - t110 * (t286 * t533 - t423 * t460) - g(1) * t286 - g(2) * t285 - g(3) * t423 + (-t164 * t604 + t165 * t600 + t244 * t499 - t245 * t497 - t285 * t460 + t286 * t461) * t101) * m(4) + (-t112 * t742 - t600 * t95 - t604 * t96) * t803 + (-t113 * t742 - t600 * t97 - t604 * t98) * t804 + (-t106 * t600 - t107 * t604 - t130 * t742) * t791 + (-t45 * t600 - t46 * t604 + t95 * t497 + t96 * t499 + (t112 * t703 - t566 * t64) * t559) * t792 + (-t47 * t600 - t48 * t604 + t97 * t497 + t98 * t499 + (t113 * t703 - t566 * t65) * t559) * t794 + (t451 * t578 - t577 * t612 - t597 * t604) * t795 + (-(t449 * t593 + t451 * t595 + t525 * t594) * qJD(4) - g(1) * t826 - g(2) * t827 - g(3) * t822 + t12 * t723 - t7 * t253 + t661 * t499 + t659 * t497 - (t685 * t7 + t596) * t600 - (t11 * t683 + t607) * t604 + (t660 * t703 + (t11 * t685 + t837) * t566) * t559 + (-qJD(5) * t451 - t650 + t734 * t434 - t724 * t317 + t684 + ((-t288 - t737) * t703 + t778 * t566) * t559) * t44 + (-qJD(5) * t449 + t316 * t724 - t434 * t733 + t499 * t683 - t604 * t687 + t692 * t742 + t270 - t651) * t43 + (-qJD(5) * t525 - t316 * t734 + t317 * t733 + t497 * t685 - t600 * t692 + t816) * t38) * m(6) + (t31 * t723 + t67 * (t497 * t255 + t684) + t66 * (t499 * t725 + t270) - t13 * t253 - (t67 * t123 + t13 * t735 + t31 * t255) * t600 - (t30 * t725 + t66 * t738 + t642) * t604 + ((t67 * (-t148 - t288) + t66 * t150) * t703 + (t31 * t148 + t30 * t735 + t66 * t776 + t67 * t92) * t566) * t559 - t67 * (-t228 * t434 + t315 * t317 + t650) - t66 * (t230 * t434 - t315 * t316 + t651) - (t67 * (-t148 * t525 + t255 * t449) + t66 * (t150 * t525 - t255 * t451)) * qJD(4) - g(1) * (t230 + t289) - g(2) * (t228 + t287) - g(3) * (t315 + t426) + (t148 * t499 - t776 * t600 - t228 * t316 + t230 * t317 - (t148 * t451 - t150 * t449) * qJD(4) + t735 * t497 + t816) * t53) * m(5) + (t820 * t524 + (t454 * t824 + t525 * t825 - t619 * t823) * t434 + (t454 * t831 + t525 * t835 - t619 * t829) * t317 + (t454 * t830 + t525 * t834 - t619 * t828) * t316 + (t449 * t842 + t451 * t841 + t525 * t836) * qJD(4)) * t800 + (t449 * t578 - t577 * t611 - t597 * t600) * t793 + (-t820 * t611 + (t344 * t824 + t345 * t823 + t449 * t825) * t434 + (t344 * t831 + t345 * t829 + t449 * t835) * t317 + (t344 * t830 + t345 * t828 + t449 * t834) * t316 + (t449 * t846 + t451 * t845 + t525 * t840) * qJD(4)) * t808 + (-t820 * t612 + (t346 * t824 + t347 * t823 + t451 * t825) * t434 + (t346 * t831 + t347 * t829 + t451 * t835) * t317 + (t346 * t830 + t347 * t828 + t451 * t834) * t316 + (t449 * t844 + t451 * t843 + t525 * t839) * qJD(4)) * t810 - (t449 * t852 + t451 * t851 + t525 * t850) * qJD(4) / 0.2e1 + (t49 + t850) * t564 * t665 - t889 * t742 / 0.2e1 + (t524 * t577 + t525 * t578 - t597 * t742) * t784; (t524 * t839 - t611 * t844 - t612 * t843) * t813 + (t524 * t840 - t611 * t846 - t612 * t845) * t812 + (t524 * t836 - t611 * t842 - t612 * t841) * t811 + (t346 * t819 + t347 * t818 - t612 * t817) * t810 + (t339 * t844 + t341 * t843 + t452 * t839 + t524 * t848 - t611 * t856 - t612 * t855) * t809 + (t344 * t819 + t345 * t818 - t611 * t817) * t808 + (t339 * t846 + t341 * t845 + t452 * t840 + t524 * t849 - t611 * t858 - t612 * t857) * t807 + t852 * t339 / 0.2e1 + t851 * t341 / 0.2e1 + (t454 * t819 + t524 * t817 - t619 * t818) * t800 + (t339 * t842 + t341 * t841 + t452 * t836 + t524 * t847 - t611 * t854 - t612 * t853) * t799 - t861 * t611 / 0.2e1 - t860 * t612 / 0.2e1 + t850 * t452 / 0.2e1 + t859 * t524 / 0.2e1 + (-g(1) * t729 - g(2) * t730 - g(3) * t721 - (t43 * t729 - t44 * t730) * t434 - (-t38 * t729 + t44 * t721) * t317 - (t38 * t730 - t43 * t721) * t316 + t594 * t452 + t595 * t341 + t593 * t339 + (t43 * t777 - t44 * t778 + t862) * t524 - (-t11 * t726 - t43 * t739 + t607) * t612 - (t596 + t815) * t611) * m(6) + (-g(1) * t193 - g(2) * t191 - g(3) * t299 + t31 * (-t148 * t524 - t255 * t611) + t30 * (t150 * t524 + t255 * t612) + t13 * (-t148 * t612 + t150 * t611) + (-t123 * t611 - t148 * t452 + t191 * t434 + t255 * t339 - t299 * t317 - t524 * t92) * t67 + (t123 * t612 + t150 * t452 - t193 * t434 - t255 * t341 + t299 * t316 + t524 * t94) * t66 + (t148 * t341 - t150 * t339 - t191 * t316 + t193 * t317 + t611 * t94 - t612 * t92) * t53) * m(5); ((-g(3) + t7) * t524 - (t12 - g(1)) * t612 - (-g(2) + t11) * t611 + (-t317 * t524 - t434 * t611 + t341) * t44 + (t316 * t524 + t434 * t612 + t339) * t43 + (t316 * t611 - t317 * t612 + t452) * t38) * m(6);];
tau = t1;

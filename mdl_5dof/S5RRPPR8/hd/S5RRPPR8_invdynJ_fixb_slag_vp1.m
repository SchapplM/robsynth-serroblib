% Calculate vector of inverse dynamics joint torques for
% S5RRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR8_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR8_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:38:08
% EndTime: 2019-12-31 19:39:13
% DurationCPUTime: 56.22s
% Computational Cost: add. (16195->1089), mult. (31008->1288), div. (0->0), fcn. (28876->8), ass. (0->504)
t532 = cos(qJ(2));
t530 = sin(qJ(2));
t807 = Icges(3,4) * t530;
t433 = Icges(3,2) * t532 + t807;
t512 = Icges(4,5) * t530;
t968 = Icges(4,3) * t532 + t433 - t512;
t800 = Icges(4,5) * t532;
t435 = Icges(4,1) * t530 - t800;
t515 = Icges(3,4) * t532;
t966 = Icges(3,1) * t530 + t435 + t515;
t430 = Icges(3,5) * t532 - Icges(3,6) * t530;
t531 = sin(qJ(1));
t533 = cos(qJ(1));
t323 = Icges(3,3) * t531 + t430 * t533;
t432 = Icges(4,4) * t532 + Icges(4,6) * t530;
t325 = Icges(4,2) * t531 + t432 * t533;
t967 = t323 + t325;
t616 = Icges(4,1) * t532 + t512;
t328 = -Icges(4,4) * t533 + t531 * t616;
t772 = t530 * t531;
t482 = Icges(3,4) * t772;
t770 = t531 * t532;
t801 = Icges(3,5) * t533;
t330 = Icges(3,1) * t770 - t482 - t801;
t963 = t328 + t330;
t329 = Icges(4,4) * t531 + t533 * t616;
t438 = Icges(3,1) * t532 - t807;
t331 = Icges(3,5) * t531 + t438 * t533;
t962 = t329 + t331;
t428 = Icges(4,3) * t530 + t800;
t320 = -Icges(4,6) * t533 + t428 * t531;
t794 = Icges(3,6) * t533;
t326 = Icges(3,4) * t770 - Icges(3,2) * t772 - t794;
t965 = -t320 + t326;
t767 = t532 * t533;
t481 = Icges(4,5) * t767;
t771 = t530 * t533;
t793 = Icges(4,6) * t531;
t321 = Icges(4,3) * t771 + t481 + t793;
t613 = -Icges(3,2) * t530 + t515;
t327 = Icges(3,6) * t531 + t533 * t613;
t964 = t321 - t327;
t961 = t428 - t613;
t960 = (Icges(3,6) - Icges(4,6)) * t532 + (Icges(4,4) + Icges(3,5)) * t530;
t957 = t438 + t616;
t956 = t968 * qJD(2);
t955 = t966 * qJD(2);
t528 = sin(pkin(8));
t773 = t530 * t528;
t813 = cos(pkin(8));
t582 = t532 * t813 + t773;
t349 = t582 * t531;
t310 = Icges(5,4) * t349;
t658 = t530 * t813;
t768 = t532 * t528;
t689 = t531 * t768;
t348 = -t531 * t658 + t689;
t182 = Icges(5,2) * t348 - Icges(5,6) * t533 - t310;
t309 = Icges(5,4) * t348;
t184 = Icges(5,1) * t349 + Icges(5,5) * t533 - t309;
t324 = -Icges(4,2) * t533 + t432 * t531;
t688 = t528 * t767;
t350 = -t533 * t658 + t688;
t351 = t582 * t533;
t954 = -t182 * t350 - t351 * t184 - t531 * t324;
t953 = -t321 * t771 - t967 * t531 - t962 * t767;
t940 = -t530 * t968 + t966 * t532;
t311 = Icges(5,4) * t351;
t183 = -Icges(5,2) * t350 - Icges(5,6) * t531 + t311;
t805 = Icges(5,4) * t350;
t186 = Icges(5,1) * t351 - Icges(5,5) * t531 - t805;
t952 = t348 * t183 - t349 * t186 - t321 * t772 + t325 * t533 - t329 * t770;
t791 = Icges(3,3) * t533;
t322 = Icges(3,5) * t770 - Icges(3,6) * t772 - t791;
t951 = -t320 * t771 - t531 * t322 - t963 * t767 + t954;
t180 = Icges(5,5) * t351 - Icges(5,6) * t350 - Icges(5,3) * t531;
t764 = -t350 * t183 + t351 * t186;
t60 = -t180 * t531 + t764;
t906 = -t327 * t771 - t953;
t888 = t60 + t906;
t950 = t963 * t530 + t965 * t532;
t949 = -t962 * t530 + t964 * t532;
t948 = t956 * t533 + (t531 * t613 - t320 - t794) * qJD(1);
t947 = t956 * t531 + (t428 * t533 - t327 + t793) * qJD(1);
t946 = -t955 * t533 + (-t438 * t531 - t328 + t801) * qJD(1);
t945 = -qJD(1) * t962 + t531 * t955;
t944 = t961 * qJD(2);
t943 = t957 * qJD(2);
t942 = -t430 - t432;
t941 = t960 * qJD(2);
t939 = -t530 * t966 - t532 * t968;
t780 = t327 * t530;
t938 = t321 * t530 + t532 * t962 - t780;
t781 = t326 * t530;
t598 = -t330 * t532 + t781;
t601 = t320 * t530 + t328 * t532;
t937 = t598 - t601;
t899 = t960 * t533;
t898 = t960 * t531;
t782 = t324 * t533;
t874 = t531 * t601;
t119 = -t782 + t874;
t178 = Icges(5,5) * t349 - Icges(5,6) * t348 + Icges(5,3) * t533;
t910 = t182 * t348 + t184 * t349;
t57 = t178 * t533 + t910;
t890 = -t322 * t533 - t531 * t598 + t119 + t57;
t889 = -t178 * t531 - t326 * t771 - t951;
t274 = t331 * t770;
t656 = t323 * t533 - t274;
t122 = -t327 * t772 - t656;
t883 = t533 * t180 + t122 - t952;
t936 = t967 * qJD(1);
t694 = pkin(8) + qJ(5);
t502 = sin(t694);
t657 = cos(t694);
t619 = t530 * t657;
t355 = -t532 * t502 + t619;
t302 = t355 * t531;
t562 = t530 * t502 + t532 * t657;
t303 = t562 * t531;
t630 = t303 * rSges(6,1) + t302 * rSges(6,2);
t147 = rSges(6,3) * t533 + t630;
t529 = -pkin(7) - qJ(4);
t692 = pkin(3) * t770;
t475 = pkin(4) * t773;
t493 = pkin(4) * t813 + pkin(3);
t725 = t532 * t493 + t475;
t733 = t725 * t531;
t231 = t692 + (qJ(4) + t529) * t533 - t733;
t759 = t147 - t231;
t935 = qJD(1) * t759;
t402 = t658 - t768;
t253 = Icges(5,5) * t402 - Icges(5,6) * t582;
t373 = Icges(5,4) * t402;
t255 = -Icges(5,2) * t582 + t373;
t804 = Icges(5,4) * t582;
t257 = Icges(5,1) * t402 - t804;
t887 = t253 * t533 - t255 * t348 + t257 * t349 + t531 * t940 - t899;
t886 = -t253 * t531 - t255 * t350 + t257 * t351 + t533 * t940 + t898;
t867 = qJD(2) - qJD(5);
t234 = t867 * t562;
t711 = qJD(1) * t533;
t129 = t234 * t531 + t355 * t711;
t305 = t562 * t533;
t907 = t867 * t355;
t130 = qJD(1) * t305 - t531 * t907;
t304 = t502 * t767 - t533 * t619;
t140 = Icges(6,5) * t305 - Icges(6,6) * t304 - Icges(6,3) * t531;
t287 = Icges(6,4) * t305;
t143 = -Icges(6,2) * t304 - Icges(6,6) * t531 + t287;
t286 = Icges(6,4) * t304;
t146 = Icges(6,1) * t305 - Icges(6,5) * t531 - t286;
t713 = qJD(1) * t531;
t127 = t234 * t533 - t355 * t713;
t128 = -t533 * t907 - t562 * t713;
t67 = Icges(6,5) * t128 + Icges(6,6) * t127 - Icges(6,3) * t711;
t69 = Icges(6,4) * t128 + Icges(6,2) * t127 - Icges(6,6) * t711;
t71 = Icges(6,1) * t128 + Icges(6,4) * t127 - Icges(6,5) * t711;
t10 = t129 * t143 + t130 * t146 - t140 * t713 + t302 * t69 + t303 * t71 + t533 * t67;
t425 = t867 * t531;
t707 = qJD(2) * t533;
t426 = -qJD(5) * t533 + t707;
t138 = Icges(6,5) * t303 + Icges(6,6) * t302 + Icges(6,3) * t533;
t802 = Icges(6,4) * t303;
t141 = Icges(6,2) * t302 + Icges(6,6) * t533 + t802;
t803 = Icges(6,4) * t302;
t145 = -Icges(6,1) * t303 - Icges(6,5) * t533 - t803;
t51 = -t138 * t531 - t304 * t141 - t145 * t305;
t639 = t140 * t531 + t304 * t143 - t305 * t146;
t236 = Icges(6,5) * t355 - Icges(6,6) * t562;
t344 = Icges(6,4) * t355;
t239 = -Icges(6,2) * t562 + t344;
t343 = Icges(6,4) * t562;
t242 = Icges(6,1) * t355 - t343;
t64 = -t236 * t531 - t239 * t304 + t242 * t305;
t14 = t64 * qJD(1) - t425 * t639 - t426 * t51;
t104 = Icges(6,5) * t234 + Icges(6,6) * t907;
t105 = Icges(6,4) * t234 + Icges(6,2) * t907;
t106 = Icges(6,1) * t234 + Icges(6,4) * t907;
t16 = -t104 * t531 - t105 * t304 + t106 * t305 + t127 * t239 + t128 * t242 - t236 * t711;
t17 = t104 * t533 + t105 * t302 + t106 * t303 + t129 * t239 + t130 * t242 - t236 * t713;
t70 = Icges(6,4) * t130 + Icges(6,2) * t129 - Icges(6,6) * t713;
t72 = Icges(6,1) * t130 + Icges(6,4) * t129 - Icges(6,5) * t713;
t22 = t141 * t907 - t145 * t234 + t355 * t72 - t562 * t70;
t23 = t143 * t907 + t146 * t234 + t355 * t71 - t562 * t69;
t699 = qJD(1) * qJD(2);
t418 = qJDD(2) * t531 + t533 * t699;
t698 = qJD(1) * qJD(5);
t300 = -qJDD(5) * t531 - t533 * t698 + t418;
t495 = t531 * t699;
t301 = -t531 * t698 + t495 + (-qJDD(2) + qJDD(5)) * t533;
t49 = t138 * t533 + t141 * t302 - t145 * t303;
t50 = t533 * t140 + t302 * t143 + t303 * t146;
t548 = qJD(1) * (Icges(6,1) * t562 + t239 + t344) + t425 * (Icges(6,1) * t304 + t143 + t287) - t426 * (-Icges(6,1) * t302 + t141 + t802);
t564 = qJD(1) * (-Icges(6,5) * t562 - Icges(6,6) * t355) + (-Icges(6,5) * t302 + Icges(6,6) * t303) * t426 - (Icges(6,5) * t304 + Icges(6,6) * t305) * t425;
t63 = t236 * t533 + t239 * t302 + t242 * t303;
t68 = Icges(6,5) * t130 + Icges(6,6) * t129 - Icges(6,3) * t713;
t7 = t127 * t141 - t128 * t145 - t138 * t711 - t304 * t70 + t305 * t72 - t531 * t68;
t75 = -t141 * t562 - t145 * t355;
t76 = -t143 * t562 + t146 * t355;
t8 = t127 * t143 + t128 * t146 - t140 * t711 - t304 * t69 + t305 * t71 - t531 * t67;
t828 = -qJD(1) / 0.2e1;
t837 = t426 / 0.2e1;
t851 = qJD(1) * (Icges(6,2) * t355 - t242 + t343) - t425 * (-Icges(6,2) * t305 + t146 - t286) + t426 * (-Icges(6,2) * t303 - t145 + t803);
t891 = t51 / 0.2e1 + t50 / 0.2e1;
t9 = t129 * t141 - t130 * t145 - t138 * t713 + t302 * t70 + t303 * t72 + t533 * t68;
t934 = (-t22 * t533 + t23 * t531 + (t531 * t75 + t533 * t76) * qJD(1)) * t828 - qJDD(1) * (t531 * t76 - t533 * t75) / 0.2e1 - t531 * (qJD(1) * t16 + qJDD(1) * t64 + t425 * t8 - t426 * t7) / 0.2e1 + t533 * (qJD(1) * t17 + qJDD(1) * t63 + t10 * t425 - t426 * t9) / 0.2e1 - (qJD(1) * t63 + t425 * t50 - t426 * t49) * t713 / 0.2e1 - t14 * t711 / 0.2e1 + (t49 * t533 - t531 * t891) * t301 + (t639 * t531 + t533 * t891) * t300 - (t304 * t851 - t305 * t548 + (-t639 * qJD(1) - t7) * t533 + (t51 * qJD(1) - t564 + t8) * t531) * t425 / 0.2e1 + (-t302 * t851 - t303 * t548 + (t49 * qJD(1) + t10) * t531 + (t50 * qJD(1) + t564 - t9) * t533) * t837;
t933 = t949 * qJD(2) + t948 * t530 + t946 * t532 + t936;
t869 = qJD(1) * t324;
t932 = -qJD(1) * t322 + t950 * qJD(2) - t947 * t530 + t945 * t532 - t869;
t931 = t965 * t533 + (-Icges(4,1) * t771 + t435 * t533 + t481 + t964) * t531;
t930 = t966 - t961;
t929 = t957 - t968;
t928 = (Icges(3,2) * t770 + t482 - t963) * t533 + (-t433 * t533 + t962) * t531;
t632 = -t349 * rSges(5,1) + t348 * rSges(5,2);
t188 = -rSges(5,3) * t533 + t632;
t927 = qJD(1) * t188;
t366 = t582 * qJD(2);
t201 = t366 * t531 + t402 * t711;
t868 = t402 * qJD(2);
t202 = qJD(1) * t351 - t531 * t868;
t926 = -Icges(5,5) * t202 - Icges(5,6) * t201 + Icges(5,3) * t713 + t937 * qJD(1) - t941 * t531 + t936;
t199 = t366 * t533 - t402 * t713;
t200 = -t533 * t868 - t582 * t713;
t925 = -Icges(5,5) * t200 - Icges(5,6) * t199 + Icges(5,3) * t711 - t869 - t941 * t533 + (-t430 * t531 + t791 - t938) * qJD(1);
t924 = Icges(5,5) * t366 + Icges(5,6) * t868 + t940 * qJD(1) + t942 * qJD(2);
t923 = t943 * t532 + t944 * t530 + t939 * qJD(2) + (-t253 + t960) * qJD(1);
t922 = t888 * t531 - t889 * t533;
t921 = t883 * t531 - t890 * t533;
t920 = -t548 * t355 + t562 * t851;
t246 = rSges(6,1) * t562 + t355 * rSges(6,2);
t917 = t246 * t425;
t916 = t246 * t426;
t913 = t886 * qJD(1);
t912 = t887 * qJD(1);
t164 = t302 * rSges(6,1) - t303 * rSges(6,2);
t166 = -t304 * rSges(6,1) - t305 * rSges(6,2);
t911 = -t164 * t425 - t166 * t426;
t905 = -t530 * t928 + t532 * t931;
t904 = (-t530 * t930 + t532 * t929) * qJD(1);
t712 = qJD(1) * t532;
t568 = -t530 * t707 - t531 * t712;
t704 = qJD(4) * t531;
t266 = pkin(3) * t568 - qJ(4) * t711 - t704;
t709 = qJD(2) * t531;
t674 = t530 * t709;
t463 = pkin(3) * t674;
t721 = qJ(4) * t713 + t463;
t267 = (pkin(3) * t712 + qJD(4)) * t533 - t721;
t789 = qJ(4) * t533;
t413 = t692 + t789;
t476 = qJ(3) * t770;
t388 = -pkin(2) * t772 + t476;
t478 = qJ(3) * t767;
t392 = -pkin(2) * t771 + t478;
t506 = qJD(3) * t530;
t683 = t388 * t709 + t392 * t707 + t506;
t676 = t530 * t713;
t705 = qJD(3) * t533;
t473 = t530 * t705;
t673 = t532 * t707;
t724 = qJ(3) * t673 + t473;
t203 = pkin(2) * t568 - qJ(3) * t676 + t724;
t708 = qJD(2) * t532;
t369 = t530 * t711 + t531 * t708;
t490 = pkin(2) * t767;
t464 = pkin(2) * t674;
t670 = t531 * t506;
t637 = t464 - t670;
t204 = qJ(3) * t369 + qJD(1) * t490 - t637;
t507 = t530 * qJ(3);
t871 = t532 * pkin(2) + t507;
t391 = t871 * t531;
t686 = t533 * t203 + t531 * t204 + t391 * t711;
t903 = t533 * t266 + t531 * t267 + t413 * t711 - t683 + t686;
t902 = -t782 + t953;
t697 = qJD(2) * qJD(3);
t901 = qJDD(3) * t530 + t532 * t697;
t443 = pkin(2) * t530 - qJ(3) * t532;
t398 = t443 * t713;
t474 = t532 * t705;
t900 = t398 - t474;
t96 = Icges(5,4) * t202 + Icges(5,2) * t201 - Icges(5,6) * t713;
t98 = Icges(5,1) * t202 + Icges(5,4) * t201 - Icges(5,5) * t713;
t897 = t937 * qJD(2) + t182 * t868 - t184 * t366 - t402 * t98 + t945 * t530 + t947 * t532 + t582 * t96;
t95 = Icges(5,4) * t200 + Icges(5,2) * t199 - Icges(5,6) * t711;
t97 = Icges(5,1) * t200 + Icges(5,4) * t199 - Icges(5,5) * t711;
t896 = t938 * qJD(2) + t183 * t868 + t186 * t366 + t402 * t97 + t946 * t530 - t948 * t532 - t582 * t95;
t248 = Icges(5,4) * t366 + Icges(5,2) * t868;
t249 = Icges(5,1) * t366 + Icges(5,4) * t868;
t895 = t199 * t255 + t200 * t257 - t248 * t350 + t249 * t351 - t531 * t924 + t533 * t923;
t894 = t201 * t255 + t202 * t257 - t248 * t348 + t249 * t349 + t531 * t923 + t533 * t924;
t893 = qJD(2) * t921 + t912;
t892 = qJD(2) * t922 + t913;
t885 = -t183 * t582 + t186 * t402 - t949;
t884 = -t182 * t582 - t184 * t402 - t950;
t882 = (-t180 * t713 + t183 * t201 + t186 * t202 - t348 * t95 + t349 * t97 + t933 * t531) * t531 + (t178 * t713 + t182 * t201 - t184 * t202 + t348 * t96 - t349 * t98 + t926 * t533 + (-t925 + t932) * t531) * t533;
t881 = (t178 * t711 + t182 * t199 - t184 * t200 + t350 * t96 - t351 * t98 + t932 * t533) * t533 + (-t180 * t711 + t183 * t199 + t186 * t200 - t350 * t95 + t351 * t97 + t925 * t531 + (-t926 + t933) * t533) * t531;
t107 = rSges(6,1) * t234 + rSges(6,2) * t907;
t825 = pkin(3) - t493;
t860 = t532 * t825 - t475;
t271 = t860 * qJD(2);
t706 = qJD(3) * t532;
t356 = qJD(2) * t871 - t706;
t640 = -pkin(3) * t708 - t356;
t832 = pkin(3) * t532;
t659 = -t871 - t832;
t879 = -qJD(2) * (-t725 + t832 + t659) - t107 + t271 + t640;
t589 = -pkin(1) - t871;
t878 = (t589 - t725) * qJD(1) - qJD(4);
t877 = ((Icges(5,5) * t348 + Icges(5,6) * t349) * t533 - (Icges(5,5) * t350 + Icges(5,6) * t351) * t531) * qJD(2) + (-Icges(5,5) * t582 - Icges(5,6) * t402 - t942) * qJD(1);
t876 = t57 + t60;
t668 = t530 * t825;
t524 = t533 * pkin(6);
t452 = pkin(1) * t531 - t524;
t422 = qJD(1) * t452;
t872 = -qJD(1) * t391 - t422;
t634 = t532 * rSges(4,1) + t530 * rSges(4,3);
t453 = t533 * pkin(1) + t531 * pkin(6);
t245 = rSges(6,1) * t355 - rSges(6,2) * t562;
t691 = pkin(4) * t768;
t307 = -t691 - t668;
t833 = pkin(3) * t530;
t660 = -t443 - t833;
t643 = -t307 + t660;
t870 = -t426 * t245 + t643 * t707;
t518 = t531 * rSges(4,2);
t341 = rSges(4,1) * t767 + rSges(4,3) * t771 + t518;
t395 = qJ(3) * t771 + t490;
t648 = t395 + t453;
t212 = t648 + t341;
t829 = g(2) * t531;
t861 = (g(1) * t533 + t829) * t530;
t559 = t531 * (-Icges(5,2) * t351 + t186 - t805) - t533 * (-Icges(5,2) * t349 + t184 - t309);
t849 = m(4) / 0.2e1;
t848 = m(5) / 0.2e1;
t847 = m(6) / 0.2e1;
t846 = -m(5) - m(6);
t845 = -pkin(2) - pkin(3);
t842 = t418 / 0.2e1;
t419 = -qJDD(2) * t533 + t495;
t841 = t419 / 0.2e1;
t834 = -rSges(4,1) - pkin(2);
t831 = g(1) * t531;
t826 = -pkin(2) - t493;
t824 = rSges(3,1) * t532;
t823 = rSges(4,1) * t530;
t517 = t531 * rSges(3,3);
t818 = -rSges(4,3) - qJ(3);
t817 = -rSges(5,3) - qJ(4);
t816 = -rSges(6,3) + t529;
t646 = t533 * t668;
t442 = pkin(4) * t688;
t732 = qJD(2) * t442 + t529 * t711;
t108 = qJD(2) * t646 + (t531 * t860 + t789) * qJD(1) + t732;
t766 = t128 * rSges(6,1) + t127 * rSges(6,2);
t73 = -rSges(6,3) * t711 + t766;
t815 = t108 + t73;
t497 = t531 * t529;
t774 = t530 * t493;
t109 = (t691 - t774) * t709 + (-t533 * t860 + t497) * qJD(1) + t721;
t631 = -rSges(6,1) * t130 - rSges(6,2) * t129;
t74 = -rSges(6,3) * t713 - t631;
t814 = t109 + t74;
t445 = rSges(3,1) * t530 + rSges(3,2) * t532;
t675 = t445 * t707;
t719 = rSges(3,2) * t772 + t533 * rSges(3,3);
t340 = rSges(3,1) * t770 - t719;
t742 = -t340 - t452;
t170 = qJD(1) * t742 - t675;
t786 = t170 * t531;
t785 = t170 * t533;
t342 = rSges(3,1) * t767 - rSges(3,2) * t771 + t517;
t269 = t342 + t453;
t171 = qJD(1) * t269 - t445 * t709;
t394 = t445 * t533;
t784 = t171 * t394;
t149 = t305 * rSges(6,1) - t304 * rSges(6,2) - rSges(6,3) * t531;
t489 = pkin(3) * t767;
t414 = -qJ(4) * t531 + t489;
t680 = t533 * t475 + t493 * t767 + t497;
t232 = -t414 + t680;
t758 = t149 + t232;
t755 = t200 * rSges(5,1) + t199 * rSges(5,2);
t752 = Icges(5,1) * t582 + t255 + t373;
t751 = Icges(5,2) * t402 - t257 + t804;
t213 = t348 * rSges(5,1) + t349 * rSges(5,2);
t214 = t350 * rSges(5,1) + t351 * rSges(5,2);
t748 = t351 * rSges(5,1) - t350 * rSges(5,2);
t741 = -t341 - t395;
t740 = t531 * t391 + t533 * t395;
t739 = -t634 * qJD(2) - t356;
t738 = qJD(1) * t392 + t531 * t706;
t737 = -t391 - t452;
t736 = -t395 - t414;
t501 = pkin(6) * t711;
t734 = qJD(1) * (-pkin(1) * t713 + t501) + qJDD(1) * t453;
t444 = -rSges(4,3) * t532 + t823;
t727 = -t443 - t444;
t726 = -t871 - t634;
t723 = rSges(4,2) * t711 + rSges(4,3) * t673;
t722 = rSges(3,2) * t676 + rSges(3,3) * t711;
t718 = t531 ^ 2 + t533 ^ 2;
t710 = qJD(2) * t530;
t703 = qJD(4) * t533;
t696 = qJD(4) * qJD(1);
t693 = pkin(3) * t771;
t690 = qJD(2) ^ 2 * t832;
t189 = -rSges(5,3) * t531 + t748;
t687 = -t189 + t736;
t685 = -t232 + t736;
t684 = t419 * t443 + t533 * t901;
t521 = t533 * rSges(4,2);
t339 = t531 * t634 - t521;
t682 = -t339 + t737;
t681 = -t413 + t737;
t679 = t501 + t724;
t678 = t524 - t733;
t677 = t533 * t834;
t669 = -pkin(1) - t824;
t664 = -t709 / 0.2e1;
t663 = t709 / 0.2e1;
t662 = -t707 / 0.2e1;
t661 = t707 / 0.2e1;
t655 = -t322 + t780;
t654 = qJD(2) * t739;
t653 = t473 - t704;
t652 = -t149 + t685;
t651 = t188 + t681;
t649 = t531 * t413 + t533 * t414 + t740;
t647 = t718 * t833;
t261 = rSges(5,1) * t402 - rSges(5,2) * t582;
t645 = -t261 + t660;
t636 = t391 * t709 + t395 * t707 - t706;
t397 = t443 * t709;
t635 = -t397 - t463 + t703;
t450 = rSges(2,1) * t533 - rSges(2,2) * t531;
t446 = rSges(2,1) * t531 + rSges(2,2) * t533;
t449 = -rSges(3,2) * t530 + t824;
t633 = rSges(5,1) * t202 + rSges(5,2) * t201;
t417 = t453 * qJD(1);
t586 = -t204 - t417 - t670;
t555 = -t267 + t586 - t703;
t566 = -t690 + (t271 - t356) * qJD(2);
t572 = -qJDD(4) * t531 + t419 * t833 + t684;
t618 = t681 - t759;
t11 = -t426 * t107 + t301 * t245 + t419 * t307 + t566 * t533 + t618 * qJDD(1) + (t555 - t814) * qJD(1) + t572;
t570 = qJDD(1) * t395 + t734 + t901 * t531 + (t203 + t473) * qJD(1);
t551 = qJD(1) * t266 + qJDD(1) * t414 + qJDD(4) * t533 + t570;
t12 = t551 + (t566 - t696) * t531 + t815 * qJD(1) + t643 * t418 + t758 * qJDD(1) - t300 * t245 - t425 * t107;
t628 = t11 * t533 + t12 * t531;
t605 = -t171 * t531 - t785;
t228 = rSges(3,1) * t568 - rSges(3,2) * t673 + t722;
t390 = t445 * t531;
t230 = -qJD(2) * t390 + (t449 * t533 + t517) * qJD(1);
t603 = t228 * t533 + t230 * t531;
t596 = t340 * t531 + t342 * t533;
t591 = -t245 + t643;
t250 = rSges(5,1) * t366 + rSges(5,2) * t868;
t590 = -t250 + t640;
t66 = (-qJD(2) * t261 + t506) * t531 + (t453 - t687) * qJD(1) + t635;
t587 = t66 * t645;
t584 = t707 * t727 + t473;
t581 = t413 * t709 + t414 * t707 + t636;
t580 = -qJD(1) * t413 + t653 + t872;
t573 = -qJDD(3) * t532 + t203 * t707 + t204 * t709 + t418 * t391 + t530 * t697;
t567 = -t690 + (-t250 - t356) * qJD(2);
t563 = t589 - t832;
t560 = (Icges(5,1) * t350 + t183 + t311) * t531 - (Icges(5,1) * t348 - t182 + t310) * t533;
t556 = t530 * t818 + t532 * t834 - pkin(1);
t554 = t266 * t707 + t267 * t709 + t418 * t413 + t573;
t486 = rSges(4,3) * t767;
t484 = rSges(4,3) * t770;
t440 = pkin(4) * t689;
t410 = t449 * qJD(2);
t393 = -rSges(4,1) * t771 + t486;
t389 = -rSges(4,1) * t772 + t484;
t370 = t673 - t676;
t368 = t718 * t710;
t265 = t442 + t646;
t264 = t531 * t668 + t440;
t229 = -t444 * t709 + (t533 * t634 + t518) * qJD(1);
t227 = rSges(4,1) * t568 - rSges(4,3) * t676 + t723;
t167 = t596 * qJD(2);
t103 = -t397 + (-qJD(2) * t444 + t506) * t531 + t212 * qJD(1);
t102 = qJD(1) * t682 + t584;
t101 = (t339 * t531 + t341 * t533) * qJD(2) + t636;
t100 = -rSges(5,3) * t713 + t633;
t99 = -rSges(5,3) * t711 + t755;
t89 = qJD(1) * t228 + qJDD(1) * t342 - t410 * t709 - t418 * t445 + t734;
t88 = -t410 * t707 + t419 * t445 + t742 * qJDD(1) + (-t230 - t417) * qJD(1);
t65 = qJD(1) * t651 + t645 * t707 + t653;
t61 = (-t188 * t531 + t189 * t533) * qJD(2) + t581;
t48 = -t245 * t425 + (-qJD(2) * t307 + t506) * t531 + (t453 - t652) * qJD(1) + t635;
t47 = qJD(1) * t618 + t653 + t870;
t46 = qJD(1) * t227 + qJDD(1) * t341 + t418 * t727 + t531 * t654 + t570;
t45 = t419 * t444 + t533 * t654 + t682 * qJDD(1) + (-t229 + t586) * qJD(1) + t684;
t35 = t147 * t425 + t149 * t426 + (-t231 * t531 + t232 * t533) * qJD(2) + t581;
t34 = t339 * t418 + t741 * t419 + (t227 * t533 + t229 * t531) * qJD(2) + t573;
t27 = t551 + (t567 - t696) * t531 + t645 * t418 + qJD(1) * t99 + qJDD(1) * t189;
t26 = t261 * t419 + t567 * t533 + t651 * qJDD(1) + (-t100 + t555) * qJD(1) + t572;
t15 = -t188 * t418 + (t100 * t531 + t533 * t99) * qJD(2) + t687 * t419 + t554;
t6 = t147 * t300 - t149 * t301 - t231 * t418 + t425 * t74 + t426 * t73 + (t108 * t533 + t109 * t531) * qJD(2) + t685 * t419 + t554;
t1 = [t14 * t837 - m(2) * (-g(1) * t446 + g(2) * t450) + (t76 + t64) * t300 / 0.2e1 + (t75 + t63) * t301 / 0.2e1 + (t16 + t23) * t425 / 0.2e1 - (t17 + t14 + t22) * t426 / 0.2e1 + (((t122 - t274 + (t323 + t781) * t533 + t951) * t533 + (t119 - t874 + t876 - t902 - t910) * t531) * qJD(2) + t913) * t661 + (t940 * qJD(2) - t105 * t562 + t106 * t355 + t234 * t242 + t907 * t239 - t248 * t582 + t249 * t402 + t255 * t868 + t257 * t366 + t943 * t530 - t944 * t532) * qJD(1) + (-g(1) * (t529 * t533 - t147 + t678) - t589 * t831 + (t12 - g(2)) * (t149 + t648 + t680) + (t531 * t589 + t533 * t816 - t630 + t678) * t11 + (t464 + t631 + t878 * t533 + (-t506 + (t774 + (-pkin(4) * t528 - qJ(3)) * t532) * qJD(2) + (-pkin(6) - t816) * qJD(1)) * t531) * t47 + (t679 + t732 + t766 + (-rSges(6,3) * qJD(1) + t710 * t826) * t533 + t878 * t531 + t47 + t935 - t580 - t870) * t48) * m(6) + (-(t580 - t65 + t927) * t66 - t587 * t707 + t26 * (t524 + t632) + t65 * (t464 - t633 + t721) + t66 * (t679 + t755) + (t66 * t710 * t845 - t65 * qJD(4) + t26 * t817) * t533 + (t26 * t563 + t65 * (-qJ(3) * t708 - t506) - t66 * qJD(4)) * t531 + ((t563 * t65 + t66 * t817) * t533 + (t65 * (rSges(5,3) - pkin(6)) + t66 * t563) * t531) * qJD(1) - g(1) * (t188 + t524 - t789) - (t532 * t845 - pkin(1) - t507) * t831 + (t27 - g(2)) * (t531 * t817 + t489 + t648 + t748)) * m(5) + (-(-qJD(1) * t339 - t102 + t584 + t872) * t103 + t102 * t637 + t103 * (t679 + t723) + (t103 * t530 * t677 + t102 * (t532 * t818 + t823) * t531) * qJD(2) + (t102 * t556 * t533 + (t102 * (-rSges(4,2) - pkin(6)) + t103 * (t589 - t634)) * t531) * qJD(1) + (-g(2) + t46) * t212 + (-g(1) + t45) * (t531 * t556 + t521 + t524)) * m(4) + (t171 * (t501 + t722) + (t445 * t786 - t784) * qJD(2) + ((-pkin(1) - t449) * t785 + (t170 * (-rSges(3,3) - pkin(6)) + t171 * t669) * t531) * qJD(1) - (-qJD(1) * t340 - t170 - t422 - t675) * t171 + (t89 - g(2)) * t269 + (t88 - g(1)) * (t669 * t531 + t524 + t719)) * m(3) + (t885 + t886) * t842 + (-t884 + t887) * t841 + (t895 + t896) * t663 + (-t255 * t582 + t257 * t402 - t239 * t562 + t242 * t355 + m(2) * (t446 ^ 2 + t450 ^ 2) + Icges(2,3) - t939) * qJDD(1) + (t892 + t894 - t897) * t662 + (((t533 * t655 - t764 + t876 + t902 + t906) * t533 + (t656 + (t655 + t178) * t531 + t889 + t952 + t954) * t531) * qJD(2) + t893 - t912) * t664; t922 * t842 + t921 * t841 + (t895 * qJD(1) + t881 * qJD(2) + t886 * qJDD(1) + t888 * t418 + t889 * t419) * t531 / 0.2e1 - (t894 * qJD(1) + t882 * qJD(2) + t887 * qJDD(1) + t883 * t418 + t890 * t419) * t533 / 0.2e1 + (t897 * t533 + t896 * t531 + (-t884 * t531 + t885 * t533) * qJD(1)) * qJD(1) / 0.2e1 + (t885 * t531 + t884 * t533) * qJDD(1) / 0.2e1 + t893 * t713 / 0.2e1 + t892 * t711 / 0.2e1 + ((-t350 * t751 + t351 * t752) * qJD(1) + (t350 * t559 + t351 * t560) * qJD(2) + (-t709 * t899 + t877) * t531 + ((t531 * t898 + t905) * qJD(2) + t904) * t533) * t664 + ((t889 * t531 + t533 * t888) * qJD(1) + t881) * t663 + ((t531 * t890 + t883 * t533) * qJD(1) + t882) * t662 + ((-t348 * t751 + t349 * t752) * qJD(1) + (t348 * t559 + t349 * t560) * qJD(2) + (-t707 * t898 - t877) * t533 + ((t533 * t899 + t905) * qJD(2) + t904) * t531) * t661 + (t6 * t649 + (t11 * t591 + t6 * t758) * t533 + (t12 * t591 + t6 * t759) * t531 - g(1) * (t442 + t478 - t166) - g(2) * (t440 + t476 - t164) - g(3) * (t871 + t725 + t246) - t826 * t861 + (t917 - t738 + t879 * t531 + (t533 * t591 + t166 - t265 + t693) * qJD(1)) * t48 + (t916 + t879 * t533 + (-t164 + t264 + t388 + (t245 + t307) * t531) * qJD(1) + t900) * t47 + (-(t264 * t531 + t265 * t533 - t647) * qJD(2) + (t815 + t935) * t533 + (qJD(1) * t652 + t814) * t531 + t903 - t911) * t35) * m(6) + (-g(1) * (t478 + t214) - g(2) * (t476 + t213) - t845 * t861 + t15 * t649 + (qJD(1) * t587 + t15 * t189 + t26 * t645) * t533 + (-t15 * t188 + t27 * t645) * t531 + (-t738 - (t214 - t693) * qJD(1) + t590 * t531) * t66 + (t590 * t533 + (t261 * t531 + t213 + t388) * qJD(1) + t900) * t65 + (-(t66 * t531 + t65 * t533) * qJD(2) + g(3)) * (-rSges(5,1) * t582 - t402 * rSges(5,2) + t659) + (-(t213 * t531 + t214 * t533 - t647) * qJD(2) + (t99 - t927) * t533 + (qJD(1) * t687 + t100) * t531 + t903) * t61) * m(5) + (-t102 * (t474 + (-t388 - t389) * qJD(1)) - t103 * (qJD(1) * t393 + t738) - t101 * t683 - ((t101 * t393 + t102 * t726) * t533 + (t101 * t389 + t103 * t726) * t531) * qJD(2) + t102 * t398 + t34 * t740 + t101 * t686 + (t45 * t727 + t102 * t739 + t34 * t341 + t101 * t227 + (t101 * t339 + t103 * t727) * qJD(1)) * t533 + (t46 * t727 + t103 * t739 + t34 * t339 + t101 * t229 + (t101 * t741 + t102 * t444) * qJD(1)) * t531 - g(1) * (t478 + t486) - g(2) * (t476 + t484) + g(3) * t726 - (g(1) * t677 + t829 * t834) * t530) * m(4) + (g(1) * t394 + g(2) * t390 - g(3) * t449 - (t170 * t390 - t784) * qJD(1) - (t167 * (-t390 * t531 - t394 * t533) + t605 * t449) * qJD(2) + (qJD(2) * t603 + t340 * t418 - t342 * t419) * t596 + t167 * ((t340 * t533 - t342 * t531) * qJD(1) + t603) + t605 * t410 + (-t89 * t531 - t88 * t533 + (-t171 * t533 + t786) * qJD(1)) * t445) * m(3) + ((t560 * t402 + t530 * t931 + t532 * t928 + t582 * t559) * qJD(2) + (t402 * t752 + t530 * t929 + t532 * t930 - t582 * t751) * qJD(1) - t920) * t828 - t934; (-m(4) + t846) * (-g(3) * t532 + t861) - m(4) * (t101 * t368 + t102 * t370 + t103 * t369) - m(5) * (t368 * t61 + t369 * t66 + t370 * t65) - m(6) * (t35 * t368 + t369 * t48 + t370 * t47) + 0.2e1 * ((t102 * t707 + t103 * t709 - t34) * t849 + (t65 * t707 + t66 * t709 - t15) * t848 + (t47 * t707 + t48 * t709 - t6) * t847) * t532 + 0.2e1 * ((qJD(2) * t101 - t102 * t713 + t103 * t711 + t45 * t533 + t46 * t531) * t849 + (qJD(2) * t61 + t26 * t533 + t27 * t531 - t65 * t713 + t66 * t711) * t848 + (qJD(2) * t35 - t47 * t713 + t48 * t711 + t628) * t847) * t530; t846 * (g(2) * t533 - t831) + m(5) * (-t26 * t531 + t27 * t533) + m(6) * (-t11 * t531 + t12 * t533); t920 * t828 + (t6 * (-t147 * t531 - t149 * t533) + (t47 * t533 + t48 * t531) * t107 + ((-t47 * t531 + t48 * t533) * qJD(1) + t628) * t245 - t47 * (-qJD(1) * t164 + t916) - t48 * (qJD(1) * t166 + t917) - g(1) * t166 - g(2) * t164 + g(3) * t246 + (-t531 * t74 - t533 * t73 + (-t147 * t533 + t149 * t531) * qJD(1) + t911) * t35) * m(6) + t934;];
tau = t1;

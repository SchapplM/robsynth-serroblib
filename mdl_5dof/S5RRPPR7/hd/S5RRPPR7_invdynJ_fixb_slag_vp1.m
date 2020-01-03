% Calculate vector of inverse dynamics joint torques for
% S5RRPPR7
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR7_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:57
% EndTime: 2019-12-31 19:36:21
% DurationCPUTime: 68.91s
% Computational Cost: add. (19479->1109), mult. (28938->1390), div. (0->0), fcn. (25462->8), ass. (0->543)
t1000 = Icges(5,4) - Icges(4,5);
t999 = Icges(5,5) - Icges(4,6);
t998 = Icges(5,1) + Icges(3,3) + Icges(4,3);
t497 = qJ(2) + pkin(8);
t474 = sin(t497);
t475 = cos(t497);
t502 = sin(qJ(2));
t505 = cos(qJ(2));
t967 = Icges(3,5) * t505 - Icges(3,6) * t502 - t1000 * t475 + t999 * t474;
t833 = Icges(4,4) * t474;
t360 = Icges(4,2) * t475 + t833;
t819 = Icges(5,6) * t474;
t602 = Icges(5,3) * t475 + t819;
t989 = -t360 - t602;
t458 = Icges(4,4) * t475;
t362 = Icges(4,1) * t474 + t458;
t818 = Icges(5,6) * t475;
t604 = Icges(5,2) * t474 + t818;
t997 = t362 + t604;
t506 = cos(qJ(1));
t996 = t998 * t506;
t503 = sin(qJ(1));
t786 = t503 * t505;
t789 = t502 * t503;
t793 = t475 * t503;
t795 = t474 * t503;
t968 = -Icges(3,5) * t786 + Icges(3,6) * t789 + t1000 * t793 - t999 * t795 + t996;
t980 = t998 * t503 + t967 * t506;
t820 = Icges(4,6) * t506;
t254 = Icges(4,4) * t793 - Icges(4,2) * t795 - t820;
t420 = Icges(5,6) * t795;
t831 = Icges(5,4) * t506;
t261 = Icges(5,2) * t793 - t420 + t831;
t821 = Icges(3,6) * t506;
t289 = Icges(3,4) * t786 - Icges(3,2) * t789 - t821;
t995 = t254 * t474 - t261 * t475 + t289 * t502;
t825 = Icges(5,5) * t506;
t259 = Icges(5,6) * t793 - Icges(5,3) * t795 + t825;
t994 = t254 + t259;
t610 = -Icges(4,2) * t474 + t458;
t255 = Icges(4,6) * t503 + t506 * t610;
t603 = -Icges(5,3) * t474 + t818;
t258 = Icges(5,5) * t503 - t506 * t603;
t993 = t255 - t258;
t425 = Icges(4,4) * t795;
t826 = Icges(4,5) * t506;
t256 = Icges(4,1) * t793 - t425 - t826;
t992 = t256 + t261;
t363 = Icges(4,1) * t475 - t833;
t257 = Icges(4,5) * t503 + t363 * t506;
t794 = t474 * t506;
t421 = Icges(5,6) * t794;
t792 = t475 * t506;
t832 = Icges(5,4) * t503;
t260 = -Icges(5,2) * t792 + t421 + t832;
t991 = t257 - t260;
t990 = Icges(3,5) * t502 + Icges(4,5) * t474 + Icges(3,6) * t505 + Icges(4,6) * t475;
t605 = Icges(5,2) * t475 - t819;
t987 = t363 + t605;
t986 = t989 * qJD(2);
t985 = t997 * qJD(2);
t984 = t603 + t610;
t834 = Icges(3,4) * t502;
t405 = Icges(3,1) * t505 - t834;
t292 = Icges(3,5) * t503 + t405 * t506;
t983 = -t257 * t793 - t258 * t795 - t292 * t786;
t608 = Icges(5,4) * t474 + Icges(5,5) * t475;
t982 = -t608 + t990;
t453 = Icges(3,4) * t789;
t827 = Icges(3,5) * t506;
t291 = Icges(3,1) * t786 - t453 - t827;
t962 = -t256 * t475 + t259 * t474 - t291 * t505 + t995;
t402 = Icges(3,2) * t505 + t834;
t486 = Icges(3,4) * t505;
t404 = Icges(3,1) * t502 + t486;
t981 = -t360 * t474 + t362 * t475 - t402 * t502 + t404 * t505;
t965 = t474 * t602 - t475 * t604 - t981;
t979 = t506 * t980 + t983;
t784 = t505 * t506;
t978 = -t256 * t792 + t259 * t794 - t291 * t784 + t503 * t968;
t913 = t257 * t792 + t258 * t794 + t292 * t784 + t503 * t980;
t611 = -Icges(3,2) * t502 + t486;
t290 = Icges(3,6) * t503 + t506 * t611;
t977 = t255 * t474 + t260 * t475 + t290 * t502;
t302 = t608 * t506;
t114 = t602 * t795 - t604 * t793 - t302;
t969 = t990 * t506;
t938 = t981 * t503 - t969;
t976 = t114 - t938;
t975 = t986 * t506 + (-t503 * t984 + t820 - t825) * qJD(1);
t974 = qJD(1) * t993 + t503 * t986;
t973 = -t985 * t506 + (-t503 * t987 + t826 - t831) * qJD(1);
t972 = t985 * t503 + (-t506 * t605 - t257 + t832) * qJD(1);
t971 = t984 * qJD(2);
t970 = t987 * qJD(2);
t926 = -t503 * t962 + t506 * t968;
t925 = -t255 * t795 - t260 * t793 - t290 * t789 - t979;
t788 = t502 * t506;
t924 = -t254 * t794 + t261 * t792 - t289 * t788 - t978;
t923 = -t255 * t794 - t260 * t792 - t290 * t788 + t913;
t921 = t289 * t505 + t291 * t502 + t992 * t474 + t475 * t994;
t920 = t290 * t505 + t292 * t502 + t474 * t991 + t475 * t993;
t966 = t982 * qJD(2);
t964 = -t402 * t505 - t404 * t502 - t474 * t997 + t475 * t989;
t963 = t257 * t475 + t258 * t474 + t292 * t505 - t977;
t911 = t982 * t503;
t922 = -t506 * t965 + t911;
t961 = t980 * qJD(1);
t504 = cos(qJ(5));
t787 = t503 * t504;
t501 = sin(qJ(5));
t790 = t501 * t506;
t331 = t474 * t787 + t790;
t785 = t504 * t506;
t791 = t501 * t503;
t332 = -t474 * t791 + t785;
t752 = t332 * rSges(6,1) - t331 * rSges(6,2);
t168 = rSges(6,3) * t793 - t752;
t852 = rSges(6,3) * t474;
t636 = rSges(6,1) * t501 + rSges(6,2) * t504;
t934 = t475 * t636;
t245 = t852 - t934;
t717 = qJD(5) * t475;
t722 = qJD(2) * t506;
t349 = -t503 * t717 + t722;
t718 = qJD(5) * t474;
t440 = qJD(1) + t718;
t813 = qJ(4) * t475;
t868 = pkin(3) * t474;
t364 = -t813 + t868;
t869 = pkin(2) * t502;
t664 = -t364 - t869;
t867 = pkin(7) * t474;
t573 = t664 - t867;
t960 = t168 * t440 + t245 * t349 - t573 * t722;
t959 = t506 ^ 2;
t958 = -t993 * t503 + t506 * t994;
t957 = t984 + t997;
t956 = t987 + t989;
t955 = (t420 + t425 + (Icges(4,2) + Icges(5,3)) * t793 - t992) * t506 + (-Icges(5,3) * t792 - t360 * t506 - t421 + t991) * t503;
t374 = t611 * qJD(2);
t375 = t405 * qJD(2);
t954 = qJD(1) * t982 + t964 * qJD(2) - t374 * t502 + t375 * t505 - t971 * t474 + t970 * t475;
t559 = qJD(2) * t402;
t188 = -t506 * t559 + (-t503 * t611 + t821) * qJD(1);
t561 = qJD(2) * t404;
t190 = -t506 * t561 + (-t405 * t503 + t827) * qJD(1);
t953 = -t920 * qJD(2) - t188 * t502 + t190 * t505 - t975 * t474 + t973 * t475 + t961;
t189 = qJD(1) * t290 - t503 * t559;
t191 = qJD(1) * t292 - t503 * t561;
t952 = t968 * qJD(1) + t921 * qJD(2) + t189 * t502 - t191 * t505 + t974 * t474 + t972 * t475;
t951 = t503 * t923 - t506 * t924;
t950 = t925 * t503 - t506 * t926;
t949 = qJD(1) * t965 + qJD(2) * t967;
t948 = qJD(1) * t962 - t503 * t966 + t961;
t947 = -t966 * t506 + (-t503 * t967 - t963 + t996) * qJD(1);
t366 = rSges(4,1) * t474 + rSges(4,2) * t475;
t313 = t366 * t503;
t487 = t503 * rSges(4,3);
t461 = t475 * rSges(4,1);
t903 = -rSges(4,2) * t474 + t461;
t154 = -qJD(2) * t313 + (t506 * t903 + t487) * qJD(1);
t340 = t903 * qJD(2);
t713 = qJD(1) * qJD(2);
t385 = -qJDD(2) * t506 + t503 * t713;
t712 = qJD(1) * qJD(3);
t692 = qJDD(3) * t503 + t385 * t869 + t506 * t712;
t272 = rSges(4,1) * t793 - rSges(4,2) * t795 - t506 * rSges(4,3);
t495 = t506 * pkin(6);
t437 = pkin(1) * t503 - t495;
t500 = -qJ(3) - pkin(6);
t467 = t506 * t500;
t494 = t505 * pkin(2);
t468 = t494 + pkin(1);
t737 = -t503 * t468 - t467;
t250 = t437 + t737;
t763 = t250 - t437;
t698 = -t272 + t763;
t507 = qJD(2) ^ 2;
t783 = t505 * t507;
t705 = pkin(2) * t783;
t492 = t503 * pkin(6);
t726 = qJD(1) * t503;
t450 = t500 * t726;
t707 = pkin(2) * t789;
t735 = qJD(2) * t707 + qJD(3) * t506;
t690 = t450 + t735;
t859 = pkin(1) - t468;
t185 = (-t506 * t859 - t492) * qJD(1) - t690;
t438 = t506 * pkin(1) + t492;
t379 = t438 * qJD(1);
t775 = -t185 - t379;
t39 = t366 * t385 + (-qJD(2) * t340 - t705) * t506 + t698 * qJDD(1) + (-t154 + t775) * qJD(1) + t692;
t946 = t39 - g(1);
t460 = t475 * rSges(6,3);
t244 = t474 * t636 + t460;
t371 = pkin(4) * t506 - pkin(7) * t793;
t945 = t371 + t737 + t752;
t944 = t922 * qJD(1);
t943 = t976 * qJD(1);
t723 = qJD(2) * t503;
t348 = t506 * t717 + t723;
t329 = t474 * t785 - t791;
t330 = t474 * t790 + t787;
t157 = Icges(6,5) * t330 + Icges(6,6) * t329 + Icges(6,3) * t792;
t830 = Icges(6,4) * t330;
t160 = Icges(6,2) * t329 + Icges(6,6) * t792 + t830;
t308 = Icges(6,4) * t329;
t163 = Icges(6,1) * t330 + Icges(6,5) * t792 + t308;
t47 = t157 * t792 + t329 * t160 + t330 * t163;
t159 = -Icges(6,5) * t332 + Icges(6,6) * t331 + Icges(6,3) * t793;
t310 = Icges(6,4) * t332;
t162 = Icges(6,2) * t331 + Icges(6,6) * t793 - t310;
t309 = Icges(6,4) * t331;
t164 = Icges(6,1) * t332 - Icges(6,5) * t793 - t309;
t48 = t159 * t792 + t329 * t162 - t164 * t330;
t606 = Icges(6,5) * t501 + Icges(6,6) * t504;
t531 = -Icges(6,3) * t474 + t475 * t606;
t829 = Icges(6,4) * t501;
t607 = Icges(6,2) * t504 + t829;
t532 = -Icges(6,6) * t474 + t475 * t607;
t828 = Icges(6,4) * t504;
t612 = Icges(6,1) * t501 + t828;
t533 = -Icges(6,5) * t474 + t475 * t612;
t81 = -t329 * t532 - t330 * t533 - t531 * t792;
t12 = t348 * t47 - t349 * t48 + t81 * t440;
t49 = t157 * t793 + t331 * t160 - t332 * t163;
t50 = t159 * t793 + t162 * t331 + t164 * t332;
t82 = -t331 * t532 + t332 * t533 - t531 * t793;
t13 = t348 * t49 - t349 * t50 + t440 * t82;
t598 = t162 * t504 - t164 * t501;
t61 = t159 * t474 - t475 * t598;
t674 = t475 * t722;
t391 = qJ(4) * t674;
t719 = qJD(4) * t506;
t410 = t474 * t719;
t688 = t474 * t726;
t686 = t475 * t726;
t935 = t474 * t722 + t686;
t131 = -pkin(3) * t935 - qJ(4) * t688 + t391 + t410;
t675 = t475 * t723;
t725 = qJD(1) * t506;
t687 = t474 * t725;
t285 = t675 + t687;
t682 = t474 * t723;
t398 = pkin(3) * t682;
t720 = qJD(4) * t503;
t676 = t474 * t720;
t685 = t475 * t725;
t132 = pkin(3) * t685 + qJ(4) * t285 - t398 + t676;
t416 = qJ(4) * t793;
t311 = -pkin(3) * t795 + t416;
t457 = t474 * qJ(4);
t902 = t475 * pkin(3) + t457;
t314 = t902 * t503;
t418 = qJ(4) * t792;
t316 = -pkin(3) * t794 + t418;
t456 = qJD(4) * t474;
t472 = pkin(6) * t725;
t478 = qJD(3) * t503;
t680 = t502 * t722;
t184 = -pkin(2) * t680 - t472 + t478 + (t503 * t859 - t467) * qJD(1);
t699 = t506 * t184 + t503 * t185 - t250 * t725;
t937 = t506 * t131 + t503 * t132 - t311 * t723 + t314 * t725 - t316 * t722 - t456 + t699;
t711 = qJD(2) * qJD(4);
t936 = qJDD(4) * t474 + t475 * t711;
t662 = rSges(5,1) * t506 - rSges(5,3) * t795;
t275 = rSges(5,2) * t793 + t662;
t933 = qJD(1) * t275;
t932 = qJD(2) * t950 - t943;
t931 = qJD(2) * t951 + t944;
t930 = t503 * t949 + t506 * t954;
t929 = t503 * t954 - t506 * t949;
t928 = t962 * qJD(2) - t189 * t505 - t191 * t502 + t972 * t474 - t974 * t475;
t927 = t963 * qJD(2) + t188 * t505 + t190 * t502 + t973 * t474 + t975 * t475;
t540 = t289 * t506 - t290 * t503;
t889 = t503 * (-t402 * t506 + t292) - t506 * (-Icges(3,2) * t786 + t291 - t453);
t919 = -t474 * t955 + t475 * t958 - t502 * t889 + t540 * t505;
t740 = t404 + t611;
t741 = -t402 + t405;
t918 = (-t474 * t957 + t475 * t956 - t502 * t740 + t505 * t741) * qJD(1);
t917 = t948 * t959 + (t953 * t503 + (-t947 + t952) * t506) * t503;
t916 = t952 * t959 + (t947 * t503 + (-t948 + t953) * t506) * t503;
t915 = t968 + t977;
t914 = t967 * qJD(1);
t912 = -t302 + t969;
t459 = t474 * rSges(5,3);
t854 = rSges(5,2) * t475;
t904 = t459 - t854;
t339 = t904 * qJD(2);
t721 = qJD(4) * t475;
t271 = qJD(2) * t902 - t721;
t701 = qJD(2) * t494;
t643 = -t271 - t701;
t689 = t494 + t902;
t906 = t904 + t689;
t910 = qJD(2) * t906 - t339 + t643;
t551 = -t366 - t869;
t908 = t506 * t551;
t387 = qJD(1) * t437;
t907 = qJD(1) * t250 - t387;
t905 = t903 + t494;
t493 = t503 * pkin(4);
t367 = pkin(7) * t792 + t493;
t697 = -t314 + t763;
t647 = t371 + t697;
t739 = t410 + t478;
t43 = qJD(1) * t647 + t739 - t960;
t166 = t330 * rSges(6,1) + t329 * rSges(6,2) + rSges(6,3) * t792;
t326 = t364 * t723;
t397 = pkin(7) * t682;
t570 = t676 - t735;
t417 = qJ(4) * t794;
t319 = pkin(3) * t792 + t417;
t445 = t506 * t468;
t655 = -t500 * t503 + t445;
t251 = t655 - t438;
t760 = t251 + t438;
t694 = t319 + t760;
t44 = t166 * t440 - t245 * t348 - t326 - t397 + (t367 + t694) * qJD(1) + t570;
t899 = t43 * t506 + t44 * t503;
t300 = (Icges(6,2) * t501 - t828) * t475;
t525 = t348 * (-Icges(6,2) * t330 + t163 + t308) - t349 * (Icges(6,2) * t332 - t164 + t309) + t440 * (-t533 + t300);
t305 = (-Icges(6,1) * t504 + t829) * t475;
t526 = t348 * (-Icges(6,1) * t329 + t160 + t830) - t349 * (-Icges(6,1) * t331 + t162 - t310) + t440 * (-t532 - t305);
t886 = m(5) / 0.2e1;
t885 = m(6) / 0.2e1;
t884 = -m(5) - m(6);
t384 = qJDD(2) * t503 + t506 * t713;
t709 = qJDD(5) * t475;
t182 = -qJD(5) * t935 + t506 * t709 + t384;
t883 = t182 / 0.2e1;
t550 = -t682 + t685;
t183 = qJD(5) * t550 + t503 * t709 + t385;
t882 = t183 / 0.2e1;
t328 = qJD(2) * t717 + qJDD(5) * t474 + qJDD(1);
t881 = t328 / 0.2e1;
t880 = -t348 / 0.2e1;
t879 = t348 / 0.2e1;
t878 = -t349 / 0.2e1;
t877 = t349 / 0.2e1;
t876 = t384 / 0.2e1;
t875 = t385 / 0.2e1;
t874 = -t440 / 0.2e1;
t873 = t440 / 0.2e1;
t872 = t503 / 0.2e1;
t871 = -t506 / 0.2e1;
t870 = -rSges(6,3) - pkin(3);
t866 = pkin(7) * t475;
t865 = g(1) * t503;
t864 = g(2) * t503;
t315 = (-rSges(6,1) * t504 + rSges(6,2) * t501) * t475;
t140 = qJD(2) * t244 + qJD(5) * t315;
t473 = pkin(4) * t725;
t234 = -pkin(7) * t935 + t473;
t745 = qJD(1) * (-pkin(1) * t726 + t472) + qJDD(1) * t438;
t544 = qJD(1) * t184 + qJDD(1) * t251 - qJDD(3) * t506 + t503 * t712 + t745;
t524 = qJDD(1) * t319 + t544 + t936 * t503 + (t131 + t410) * qJD(1);
t641 = -t494 - t866;
t535 = -qJD(2) * t271 + t507 * t641;
t652 = qJD(1) * t474 + qJD(5);
t534 = -t503 * t652 + t674;
t576 = t440 * t501;
t138 = t504 * t534 - t506 * t576;
t577 = t504 * t440;
t139 = t501 * t534 + t506 * t577;
t782 = t139 * rSges(6,1) + t138 * rSges(6,2);
t78 = -rSges(6,3) * t935 + t782;
t8 = qJD(1) * t234 + qJDD(1) * t367 - t348 * t140 + t328 * t166 - t182 * t245 + t384 * t573 + t440 * t78 + t503 * t535 + t524;
t863 = t8 * t506;
t233 = qJD(1) * t367 - t397;
t563 = -t132 - t676 + t775;
t568 = t385 * t364 + t506 * t936 + t692;
t575 = t652 * t506;
t724 = qJD(2) * t475;
t136 = t504 * t575 + (t504 * t724 - t576) * t503;
t137 = t503 * t577 + (t575 + t675) * t501;
t638 = rSges(6,1) * t137 + rSges(6,2) * t136;
t77 = rSges(6,3) * t550 + t638;
t9 = t385 * t867 - t349 * t140 - t328 * t168 + t183 * t245 - t440 * t77 + t535 * t506 + t647 * qJDD(1) + (-t233 + t563) * qJD(1) + t568;
t862 = t9 * t503;
t236 = Icges(6,3) * t475 + t474 * t606;
t297 = (-Icges(6,5) * t504 + Icges(6,6) * t501) * t475;
t133 = qJD(2) * t236 + qJD(5) * t297;
t238 = Icges(6,6) * t475 + t474 * t607;
t134 = qJD(2) * t238 + qJD(5) * t300;
t240 = Icges(6,5) * t475 + t474 * t612;
t135 = qJD(2) * t240 + qJD(5) * t305;
t594 = -t501 * t533 - t504 * t532;
t21 = (qJD(2) * t594 + t133) * t474 + (-qJD(2) * t531 - t134 * t504 - t135 * t501 + (-t501 * t532 + t504 * t533) * qJD(5)) * t475;
t85 = -t474 * t531 - t475 * t594;
t858 = t21 * t440 + t85 * t328;
t857 = rSges(3,1) * t505;
t599 = t160 * t504 + t163 * t501;
t72 = Icges(6,5) * t139 + Icges(6,6) * t138 - Icges(6,3) * t935;
t74 = Icges(6,4) * t139 + Icges(6,2) * t138 - Icges(6,6) * t935;
t76 = Icges(6,1) * t139 + Icges(6,4) * t138 - Icges(6,5) * t935;
t10 = (qJD(2) * t599 + t72) * t474 + (qJD(2) * t157 - t501 * t76 - t504 * t74 + (t160 * t501 - t163 * t504) * qJD(5)) * t475;
t851 = t10 * t348;
t71 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t550;
t73 = Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t550;
t75 = Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t550;
t11 = (qJD(2) * t598 + t71) * t474 + (qJD(2) * t159 - t501 * t75 - t504 * t73 + (t162 * t501 + t164 * t504) * qJD(5)) * t475;
t850 = t11 * t349;
t847 = t43 * t503;
t844 = t44 * t506;
t842 = t474 * rSges(5,2);
t489 = t503 * rSges(5,1);
t488 = t503 * rSges(3,3);
t60 = t157 * t474 - t475 * t599;
t841 = t60 * t182;
t840 = t61 * t183;
t548 = qJD(2) * t908 + t478;
t88 = qJD(1) * t698 + t548;
t839 = t88 * t366;
t412 = rSges(3,1) * t502 + rSges(3,2) * t505;
t683 = t412 * t722;
t734 = rSges(3,2) * t789 + t506 * rSges(3,3);
t324 = rSges(3,1) * t786 - t734;
t751 = -t324 - t437;
t173 = qJD(1) * t751 - t683;
t812 = t173 * t503;
t811 = t173 * t506;
t325 = rSges(3,1) * t784 - rSges(3,2) * t788 + t488;
t228 = t325 + t438;
t174 = qJD(1) * t228 - t412 * t723;
t351 = t412 * t506;
t810 = t174 * t351;
t777 = t166 + t367;
t776 = t168 - t371;
t769 = -t250 * t723 + t251 * t722;
t768 = -t503 * t250 + t506 * t251;
t273 = rSges(4,1) * t792 - rSges(4,2) * t794 + t487;
t762 = -t251 - t273;
t761 = -t251 - t319;
t755 = qJD(1) * t316 + t475 * t720;
t684 = t502 * t726;
t449 = pkin(2) * t684;
t750 = t364 * t726 + t449;
t744 = t934 * t503;
t743 = t934 * t506;
t742 = rSges(4,2) * t688 + rSges(4,3) * t725;
t738 = t475 * t719 + t449;
t312 = rSges(5,2) * t795 + rSges(5,3) * t793;
t317 = rSges(5,2) * t794 + rSges(5,3) * t792;
t736 = rSges(3,2) * t684 + rSges(3,3) * t725;
t733 = t503 ^ 2 + t959;
t708 = -pkin(7) + t870;
t706 = pkin(2) * t788;
t700 = t184 * t722 + t185 * t723 - t384 * t250;
t274 = -rSges(5,2) * t792 + rSges(5,3) * t794 + t489;
t696 = -t274 + t761;
t695 = -t367 + t761;
t691 = t391 + t739;
t679 = t505 * t722;
t673 = -pkin(1) - t857;
t670 = t725 / 0.2e1;
t669 = -t723 / 0.2e1;
t668 = t723 / 0.2e1;
t667 = -t722 / 0.2e1;
t666 = t722 / 0.2e1;
t654 = t503 * t708;
t653 = t506 * t708;
t650 = t503 * t314 + t506 * t319 + t768;
t649 = t275 + t697;
t648 = rSges(5,1) * t725 + rSges(5,2) * t935 + rSges(5,3) * t674;
t646 = t733 * t869;
t634 = rSges(5,3) * t475 + t842;
t645 = t634 + t664;
t642 = -t868 - t869;
t415 = rSges(2,1) * t506 - rSges(2,2) * t503;
t413 = rSges(2,1) * t503 + rSges(2,2) * t506;
t414 = -rSges(3,2) * t502 + t857;
t626 = t47 * t506 + t48 * t503;
t625 = t47 * t503 - t48 * t506;
t624 = t49 * t506 + t50 * t503;
t623 = t49 * t503 - t50 * t506;
t622 = t503 * t61 + t506 * t60;
t621 = t503 * t60 - t506 * t61;
t616 = -qJD(1) * t314 + t739 + t907;
t153 = -rSges(4,1) * t935 - rSges(4,2) * t674 + t742;
t600 = t153 * t506 + t154 * t503;
t597 = t166 * t503 - t168 * t506;
t596 = -t174 * t503 - t811;
t192 = -rSges(3,2) * t679 + (-t505 * t726 - t680) * rSges(3,1) + t736;
t350 = t412 * t503;
t193 = -qJD(2) * t350 + (t414 * t506 + t488) * qJD(1);
t595 = t192 * t506 + t193 * t503;
t588 = t272 * t503 + t273 * t506;
t585 = t324 * t503 + t325 * t506;
t571 = t655 + t319;
t567 = t642 * t506;
t564 = -t245 + t573;
t562 = t314 * t723 + t319 * t722 - t721 + t769;
t552 = t645 * t722;
t547 = -t705 + (-t271 - t339) * qJD(2);
t546 = -t157 * t348 + t159 * t349 + t440 * t531;
t545 = (Icges(6,5) * t329 - Icges(6,6) * t330) * t348 - (Icges(6,5) * t331 + Icges(6,6) * t332) * t349 + t297 * t440;
t539 = -t468 - t902 - t459;
t537 = -pkin(7) * t724 - t140 + t643;
t536 = t475 * t545;
t530 = -qJDD(4) * t475 + t131 * t722 + t132 * t723 + t384 * t314 + t474 * t711 + t700;
t30 = t166 * t349 + t168 * t348 + (t367 * t506 - t371 * t503) * qJD(2) + t562;
t523 = t30 * t597 + (t844 - t847) * t245;
t510 = (t531 * t506 + t599) * t348 - (t531 * t503 + t598) * t349 + (t236 + t594) * t440;
t509 = t510 * t475;
t376 = t414 * qJD(2);
t318 = t366 * t506;
t286 = t674 - t688;
t284 = t733 * t474 * qJD(2);
t212 = -rSges(6,3) * t794 + t743;
t211 = -rSges(6,3) * t795 + t744;
t210 = t533 * t506;
t209 = t533 * t503;
t208 = t532 * t506;
t207 = t532 * t503;
t201 = rSges(6,1) * t331 + rSges(6,2) * t332;
t200 = rSges(6,1) * t329 - rSges(6,2) * t330;
t169 = t585 * qJD(2);
t156 = -rSges(5,3) * t688 + t648;
t155 = t634 * t723 + (t506 * t904 + t489) * qJD(1);
t89 = -t366 * t723 + (t273 + t760) * qJD(1) - t735;
t84 = qJD(2) * t588 + t769;
t80 = qJD(1) * t192 + qJDD(1) * t325 - t376 * t723 - t384 * t412 + t745;
t79 = -t376 * t722 + t385 * t412 + t751 * qJDD(1) + (-t193 - t379) * qJD(1);
t65 = -t326 + (qJD(2) * t634 + t456) * t503 + (t274 + t694) * qJD(1) - t735;
t64 = qJD(1) * t649 + t552 + t739;
t59 = (t274 * t506 - t275 * t503) * qJD(2) + t562;
t40 = -t340 * t723 + qJD(1) * t153 + qJDD(1) * t273 - t366 * t384 + (-t384 * t502 - t503 * t783) * pkin(2) + t544;
t19 = qJD(1) * t156 + qJDD(1) * t274 + t384 * t645 + t503 * t547 + t524;
t18 = -t634 * t385 + t547 * t506 + t649 * qJDD(1) + (-t155 + t563) * qJD(1) + t568;
t17 = t133 * t792 + t134 * t329 + t135 * t330 - t138 * t532 - t139 * t533 + t531 * t935;
t16 = t133 * t793 + t134 * t331 - t135 * t332 - t136 * t532 - t137 * t533 - t531 * t550;
t15 = t348 * t60 - t349 * t61 + t440 * t85;
t14 = -t275 * t384 + (t155 * t503 + t156 * t506) * qJD(2) + t696 * t385 + t530;
t7 = t138 * t162 - t139 * t164 - t159 * t935 + t329 * t73 + t330 * t75 + t71 * t792;
t6 = t138 * t160 + t139 * t163 - t157 * t935 + t329 * t74 + t330 * t76 + t72 * t792;
t5 = t136 * t162 - t137 * t164 + t159 * t550 + t331 * t73 - t332 * t75 + t71 * t793;
t4 = t136 * t160 + t137 * t163 + t157 * t550 + t331 * t74 - t332 * t76 + t72 * t793;
t3 = -t166 * t183 + t168 * t182 + t348 * t77 + t349 * t78 - t371 * t384 + (t233 * t503 + t234 * t506) * qJD(2) + t695 * t385 + t530;
t2 = t17 * t440 + t182 * t47 + t183 * t48 + t328 * t81 + t348 * t6 - t349 * t7;
t1 = t16 * t440 + t182 * t49 + t183 * t50 + t328 * t82 + t348 * t4 - t349 * t5;
t20 = [(t16 + t12) * t878 + (((-t813 + t852) * t847 + (t474 * t708 - t869) * t844) * qJD(2) - g(1) * t945 - (t475 * t870 - t457) * t865 + (t8 - g(2)) * (t571 + t777) + ((-t902 - t460) * t503 + t945) * t9 + (t397 + t398 - t638 + t690 - t456 * t503 + (t653 * t475 - t417 - t445 - t493) * qJD(1)) * t43 + (t43 + t473 - t616 + t691 + t782 + (-qJ(4) * t795 + t654 * t475 - t371 + t737) * qJD(1) + t960) * t44) * m(6) + (m(2) * (t413 ^ 2 + t415 ^ 2) + Icges(2,3) - t964) * qJDD(1) + (-t965 * qJD(2) + t374 * t505 + t375 * t502 + t970 * t474 + t971 * t475) * qJD(1) + t82 * t882 + t81 * t883 + t17 * t879 + t858 + t851 / 0.2e1 - t850 / 0.2e1 + (t921 + t938) * t875 + (t64 * (t398 + t450 - t570) + t65 * (t648 + t691) + (t65 * t567 + t64 * (-t842 + (-rSges(5,3) - qJ(4)) * t475) * t503) * qJD(2) + ((-t64 * rSges(5,1) + t539 * t65) * t503 + (t64 * (t539 + t854) - t65 * t500) * t506) * qJD(1) - (t552 + t616 - t64 + t933) * t65 + (t19 - g(2)) * (t274 + t571) + (t18 - g(1)) * ((-t457 + (rSges(5,2) - pkin(3)) * t475) * t503 + t662 + t737)) * m(5) + t12 * t877 + t840 / 0.2e1 + t841 / 0.2e1 + (t920 + t922) * t876 + ((t913 * t503 + ((t980 + t995) * t506 + t925 + t978 + t983) * t506) * qJD(2) + t944) * t666 + (t927 + t930) * t668 + (-t928 + t929 + t931) * t667 + (((t506 * t915 - t913 + t923) * t506 + (t503 * t915 + t924 + t979) * t503) * qJD(2) + t932 + t943) * t669 + (-(-qJD(1) * t324 - t173 - t387 - t683) * t174 + t174 * (t472 + t736) + (t412 * t812 - t810) * qJD(2) + ((-pkin(1) - t414) * t811 + (t173 * (-rSges(3,3) - pkin(6)) + t174 * t673) * t503) * qJD(1) + (t80 - g(2)) * t228 + (t79 - g(1)) * (t673 * t503 + t495 + t734)) * m(3) - t385 * t114 / 0.2e1 + (t88 * t690 + t89 * (t478 + t742) + (t503 * t839 + t89 * t908) * qJD(2) + ((-t88 * rSges(4,3) + t89 * (-t468 - t461)) * t503 + (t88 * (-t468 - t903) - t89 * t500) * t506) * qJD(1) - (-qJD(1) * t272 + t548 - t88 + t907) * t89 + (t40 - g(2)) * (t273 + t655) + t946 * (-t272 + t737)) * m(4) - m(2) * (-g(1) * t413 + g(2) * t415); (-t43 * (-qJD(1) * t311 - t211 * t440 - t244 * t349 + t738) - t44 * (-pkin(7) * t687 - qJD(1) * t706 + t212 * t440 - t244 * t348 + t755) - ((t166 * t44 - t168 * t43) * t475 + t523 * t474) * qJD(5) - t899 * (-t902 + t641) * qJD(2) + t43 * t750 + t3 * t650 + (qJD(1) * t245 * t43 + t3 * t776 + t44 * t537 + t8 * t564) * t503 + (t3 * t777 + t43 * t537 + (qJD(1) * t44 + t9) * t564) * t506 - g(1) * (t418 - t706 + t743) - g(2) * (t416 - t707 + t744) - g(3) * (t244 + t689 + t866) - (g(1) * t653 + g(2) * t654) * t474 + (-t211 * t348 - t212 * t349 - (-t733 * t867 - t646) * qJD(2) + (t233 + t77 + (-t166 + t695) * qJD(1)) * t503 + (qJD(1) * t776 + t234 + t78) * t506 + t937) * t30) * m(6) - ((t474 * t958 + t475 * t955 + t540 * t502 + t505 * t889) * qJD(2) + (t956 * t474 + t957 * t475 + t502 * t741 + t505 * t740) * qJD(1)) * qJD(1) / 0.2e1 + (((-t208 * t504 - t210 * t501 + t157) * t348 - (-t207 * t504 - t209 * t501 + t159) * t349 + (-t238 * t504 - t240 * t501 - t531) * t440 + t85 * qJD(5)) * t475 + (-qJD(5) * t622 + t510) * t474) * t874 + t950 * t875 + t951 * t876 + (g(1) * t351 + g(2) * t350 - g(3) * t414 - (t173 * t350 - t810) * qJD(1) - (t169 * (-t350 * t503 - t351 * t506) + t596 * t414) * qJD(2) + (qJD(2) * t595 + t324 * t384 - t325 * t385) * t585 + t169 * ((t324 * t506 - t325 * t503) * qJD(1) + t595) + t596 * t376 + (-t80 * t503 - t79 * t506 + (-t174 * t506 + t812) * qJD(1)) * t412) * m(3) + t621 * t881 + t623 * t882 + t625 * t883 + (qJD(1) * t622 + t10 * t503 - t11 * t506) * t873 + ((t208 * t331 - t210 * t332) * t348 - (t207 * t331 - t209 * t332) * t349 + (t238 * t331 - t240 * t332) * t440 + (t475 * t82 - t49 * t794) * qJD(5) + ((-qJD(5) * t50 + t546) * t474 + t509) * t503) * t877 + (qJD(1) * t624 + t4 * t503 - t5 * t506) * t878 + (qJD(1) * t626 + t503 * t6 - t506 * t7) * t879 + ((t208 * t329 + t210 * t330) * t348 - (t207 * t329 + t209 * t330) * t349 + (t238 * t329 + t240 * t330) * t440 + (t475 * t81 - t48 * t795) * qJD(5) + ((-qJD(5) * t47 + t546) * t474 + t509) * t506) * t880 + ((-t722 * t911 - t914) * t506 + ((t506 * t912 + t919) * qJD(2) + t918) * t503) * t666 + ((-t723 * t912 + t914) * t503 + ((t503 * t911 + t919) * qJD(2) + t918) * t506) * t669 - t15 * t717 / 0.2e1 + (t929 * qJD(1) + t917 * qJD(2) - t976 * qJDD(1) + t925 * t384 + t926 * t385 + t1) * t871 + (-g(1) * (t418 + t567 + t317) - g(2) * (t503 * t642 + t312 + t416) - g(3) * t906 + t14 * t650 + (t14 * t274 + t18 * t645) * t506 + (-t14 * t275 + t19 * t645) * t503 + (-t755 + t910 * t503 + (t506 * t645 - t317 + t706) * qJD(1)) * t65 + (-t738 + t750 + t910 * t506 + (-t503 * t634 + t311 + t312) * qJD(1)) * t64 + ((t156 - t933) * t506 + (qJD(1) * t696 + t155) * t503 - (t312 * t503 + t317 * t506 - t646) * qJD(2) + t937) * t59) * m(5) + (t13 + t932) * t726 / 0.2e1 + (-g(3) * t905 - t551 * t864 + t88 * (-pkin(2) * t679 - t340 * t506) + (qJD(2) * t600 + t272 * t384 + t385 * t762 + t700) * (t588 + t768) + t84 * (t600 + t699) + (t84 * t272 + t551 * t89) * t725 + (t40 * t551 + t89 * (-t340 - t701) + (t762 * t84 + t839) * qJD(1)) * t503 - (t88 * t313 + t89 * (-t318 - t706)) * qJD(1) - (-t84 * t646 + (-t84 * t318 - t88 * t905) * t506 + (-t84 * t313 - t89 * t905) * t503) * qJD(2) + t946 * t908) * m(4) + ((t503 * t924 + t506 * t923) * qJD(1) + t916) * t668 + ((t503 * t926 + t506 * t925) * qJD(1) + t917) * t667 + (t928 * t506 + t927 * t503 + (t503 * t921 + t506 * t920) * qJD(1)) * qJD(1) / 0.2e1 + (qJD(1) * t930 + qJD(2) * t916 + qJDD(1) * t922 + t384 * t923 + t385 * t924 + t2) * t872 + (t12 + t931) * t670 + (t503 * t920 - t506 * t921) * qJDD(1) / 0.2e1 + (t12 * t506 + t13 * t503) * t718 / 0.2e1; (-m(4) + t884) * (-g(2) * t506 + t865) + 0.2e1 * (t862 / 0.2e1 - t863 / 0.2e1) * m(6) + 0.2e1 * (t18 * t872 + t19 * t871) * m(5) + 0.2e1 * (t39 * t872 + t40 * t871) * m(4); t884 * (-g(3) * t475 + (g(1) * t506 + t864) * t474) - m(5) * (t284 * t59 + t285 * t65 + t286 * t64) - m(6) * (t284 * t30 + t285 * t44 + t286 * t43) + 0.2e1 * ((t64 * t722 + t65 * t723 - t14) * t886 + (t43 * t722 + t44 * t723 - t3) * t885) * t475 + 0.2e1 * ((qJD(2) * t59 + t18 * t506 + t19 * t503 - t64 * t726 + t65 * t725) * t886 + (qJD(2) * t30 - t43 * t726 + t44 * t725 + t503 * t8 + t506 * t9) * t885) * t474; t2 * t792 / 0.2e1 + (t474 * t81 + t475 * t626) * t883 + ((-qJD(2) * t626 + t17) * t474 + (-qJD(1) * t625 + qJD(2) * t81 + t503 * t7 + t506 * t6) * t475) * t879 + t1 * t793 / 0.2e1 + (t474 * t82 + t475 * t624) * t882 + ((-qJD(2) * t624 + t16) * t474 + (-qJD(1) * t623 + qJD(2) * t82 + t4 * t506 + t5 * t503) * t475) * t878 + t15 * t724 / 0.2e1 + t474 * (t840 + t841 - t850 + t851 + t858) / 0.2e1 + (t474 * t85 + t475 * t622) * t881 + ((-qJD(2) * t622 + t21) * t474 + (-qJD(1) * t621 + qJD(2) * t85 + t10 * t506 + t11 * t503) * t475) * t873 + (t525 * t329 - t330 * t526 + t506 * t536) * t880 + (t331 * t525 + t332 * t526 + t503 * t536) * t877 + (t545 * t474 + (t526 * t501 - t504 * t525) * t475) * t874 + (t474 * t669 + t475 * t670) * t13 + (-t686 / 0.2e1 + t474 * t667) * t12 + ((qJD(2) * t523 + t8 * t166 - t9 * t168 - t43 * t77 + t44 * t78) * t474 + (t43 * (-qJD(2) * t168 + t140 * t503) + t44 * (qJD(2) * t166 - t140 * t506) - t3 * t597 + t30 * (-t166 * t725 - t168 * t726 - t503 * t78 + t506 * t77) + (qJD(1) * t899 + t862 - t863) * t245) * t475 - t43 * (-t201 * t440 - t315 * t349) - t44 * (t200 * t440 - t315 * t348) - t30 * (t200 * t349 + t201 * t348) - g(1) * t200 - g(2) * t201 - g(3) * t315) * m(6);];
tau = t20;

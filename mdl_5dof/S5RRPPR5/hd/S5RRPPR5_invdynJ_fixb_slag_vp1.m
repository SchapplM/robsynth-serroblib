% Calculate vector of inverse dynamics joint torques for
% S5RRPPR5
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:54
% EndTime: 2019-12-31 19:29:56
% DurationCPUTime: 53.26s
% Computational Cost: add. (18322->972), mult. (26773->1175), div. (0->0), fcn. (23822->8), ass. (0->480)
t972 = Icges(3,3) + Icges(4,3);
t503 = qJ(2) + pkin(8);
t476 = sin(t503);
t477 = cos(t503);
t508 = sin(qJ(2));
t510 = cos(qJ(2));
t971 = Icges(3,5) * t510 + Icges(4,5) * t477 - Icges(3,6) * t508 - Icges(4,6) * t476;
t370 = Icges(5,4) * t477 + Icges(5,6) * t476;
t509 = sin(qJ(1));
t511 = cos(qJ(1));
t269 = Icges(5,2) * t509 + t370 * t511;
t952 = t972 * t509 + t971 * t511;
t970 = t269 + t952;
t792 = Icges(4,4) * t476;
t371 = Icges(4,2) * t477 + t792;
t460 = Icges(5,5) * t476;
t967 = Icges(5,3) * t477 + t371 - t460;
t786 = Icges(5,5) * t477;
t373 = Icges(5,1) * t476 - t786;
t461 = Icges(4,4) * t477;
t969 = Icges(4,1) * t476 + t373 + t461;
t968 = t972 * t511;
t750 = t509 * t510;
t752 = t508 * t509;
t754 = t477 * t509;
t756 = t476 * t509;
t943 = -Icges(3,5) * t750 - Icges(4,5) * t754 + Icges(3,6) * t752 + Icges(4,6) * t756 + t968;
t601 = Icges(5,1) * t477 + t460;
t272 = -Icges(5,4) * t511 + t509 * t601;
t431 = Icges(4,4) * t756;
t787 = Icges(4,5) * t511;
t274 = Icges(4,1) * t754 - t431 - t787;
t962 = -t272 - t274;
t273 = Icges(5,4) * t509 + t511 * t601;
t376 = Icges(4,1) * t477 - t792;
t275 = Icges(4,5) * t509 + t376 * t511;
t961 = t273 + t275;
t778 = Icges(4,6) * t511;
t270 = Icges(4,4) * t754 - Icges(4,2) * t756 - t778;
t779 = Icges(3,6) * t511;
t314 = Icges(3,4) * t750 - Icges(3,2) * t752 - t779;
t966 = t270 * t476 + t314 * t508;
t793 = Icges(3,4) * t508;
t415 = Icges(3,1) * t510 - t793;
t317 = Icges(3,5) * t509 + t415 * t511;
t965 = -t275 * t754 - t317 * t750;
t366 = Icges(5,3) * t476 + t786;
t264 = -Icges(5,6) * t511 + t366 * t509;
t964 = t264 - t270;
t753 = t477 * t511;
t430 = Icges(5,5) * t753;
t755 = t476 * t511;
t777 = Icges(5,6) * t509;
t265 = Icges(5,3) * t755 + t430 + t777;
t598 = -Icges(4,2) * t476 + t461;
t271 = Icges(4,6) * t509 + t511 * t598;
t963 = -t265 + t271;
t960 = t366 - t598;
t957 = t376 + t601;
t956 = t967 * qJD(2);
t955 = t969 * qJD(2);
t455 = Icges(3,4) * t752;
t788 = Icges(3,5) * t511;
t316 = Icges(3,1) * t750 - t455 - t788;
t954 = t274 * t477 + t316 * t510 - t966;
t953 = Icges(3,5) * t508 + Icges(3,6) * t510 + (Icges(4,6) - Icges(5,6)) * t477 + (Icges(5,4) + Icges(4,5)) * t476;
t623 = -t265 * t756 + t269 * t511 - t273 * t754;
t492 = Icges(3,4) * t510;
t599 = -Icges(3,2) * t508 + t492;
t315 = Icges(3,6) * t509 + t511 * t599;
t951 = -t952 * t511 - t965;
t905 = -t271 * t756 - t315 * t752 + t951;
t885 = -t623 + t905;
t950 = t271 * t476 + t315 * t508;
t749 = t510 * t511;
t949 = t265 * t755 + t317 * t749 + t970 * t509 + t961 * t753;
t268 = -Icges(5,2) * t511 + t370 * t509;
t245 = t509 * t268;
t948 = -t264 * t755 - t316 * t749 + t943 * t509 + t962 * t753 - t245;
t412 = Icges(3,2) * t510 + t793;
t414 = Icges(3,1) * t508 + t492;
t937 = t412 * t508 - t414 * t510 + t967 * t476 - t477 * t969;
t947 = t956 * t511 + (t509 * t598 - t264 - t778) * qJD(1);
t946 = t956 * t509 + (t366 * t511 - t271 + t777) * qJD(1);
t945 = -t955 * t511 + (-t376 * t509 - t272 + t787) * qJD(1);
t944 = -t961 * qJD(1) + t955 * t509;
t942 = t960 * qJD(2);
t941 = t957 * qJD(2);
t746 = t511 * t268;
t589 = t264 * t476 + t272 * t477;
t870 = t509 * t589;
t92 = -t746 + t870;
t889 = t954 * t509 + t943 * t511 + t92;
t751 = t508 * t511;
t888 = -t270 * t755 - t314 * t751 - t948;
t887 = -t271 * t755 - t315 * t751 + t949;
t939 = t370 + t971;
t938 = t953 * qJD(2);
t936 = -t412 * t510 - t414 * t508 - t476 * t969 - t477 * t967;
t935 = t265 * t476 + t317 * t510 + t961 * t477 - t950;
t934 = -t589 - t954;
t879 = t315 * t510 + t317 * t508 + t961 * t476 + t963 * t477;
t878 = -t314 * t510 - t316 * t508 + t962 * t476 + t964 * t477;
t874 = t953 * t511;
t873 = t953 * t509;
t507 = sin(qJ(5));
t823 = cos(qJ(5));
t661 = t476 * t823;
t302 = t507 * t754 - t509 * t661;
t567 = t476 * t507 + t477 * t823;
t303 = t567 * t509;
t618 = t303 * rSges(6,1) - t302 * rSges(6,2);
t147 = rSges(6,3) * t511 + t618;
t818 = pkin(7) * t511;
t377 = pkin(4) * t754 + t818;
t740 = t147 + t377;
t933 = qJD(1) * t740;
t884 = -t509 * t937 - t874;
t883 = -t511 * t937 + t873;
t932 = t970 * qJD(1);
t860 = qJD(2) - qJD(5);
t208 = t860 * t567;
t353 = -t477 * t507 + t661;
t692 = qJD(1) * t511;
t118 = t208 * t509 + t353 * t692;
t305 = t567 * t511;
t902 = t860 * t353;
t119 = qJD(1) * t305 - t509 * t902;
t304 = t507 * t753 - t511 * t661;
t140 = Icges(6,5) * t305 - Icges(6,6) * t304 - Icges(6,3) * t509;
t263 = Icges(6,4) * t305;
t143 = -Icges(6,2) * t304 - Icges(6,6) * t509 + t263;
t262 = Icges(6,4) * t304;
t146 = Icges(6,1) * t305 - Icges(6,5) * t509 - t262;
t693 = qJD(1) * t509;
t116 = t208 * t511 - t353 * t693;
t117 = -t511 * t902 - t567 * t693;
t66 = Icges(6,5) * t117 + Icges(6,6) * t116 - Icges(6,3) * t692;
t68 = Icges(6,4) * t117 + Icges(6,2) * t116 - Icges(6,6) * t692;
t70 = Icges(6,1) * t117 + Icges(6,4) * t116 - Icges(6,5) * t692;
t10 = t118 * t143 + t119 * t146 - t140 * t693 - t302 * t68 + t303 * t70 + t511 * t66;
t408 = t860 * t509;
t689 = qJD(2) * t511;
t409 = -qJD(5) * t511 + t689;
t138 = Icges(6,5) * t303 - Icges(6,6) * t302 + Icges(6,3) * t511;
t789 = Icges(6,4) * t303;
t141 = -Icges(6,2) * t302 + Icges(6,6) * t511 + t789;
t790 = Icges(6,4) * t302;
t145 = -Icges(6,1) * t303 - Icges(6,5) * t511 + t790;
t50 = -t138 * t509 - t304 * t141 - t145 * t305;
t624 = t140 * t509 + t304 * t143 - t305 * t146;
t212 = Icges(6,5) * t353 - Icges(6,6) * t567;
t341 = Icges(6,4) * t353;
t215 = -Icges(6,2) * t567 + t341;
t340 = Icges(6,4) * t567;
t218 = Icges(6,1) * t353 - t340;
t65 = -t212 * t509 - t215 * t304 + t218 * t305;
t14 = t65 * qJD(1) - t408 * t624 - t409 * t50;
t100 = Icges(6,5) * t208 + Icges(6,6) * t902;
t101 = Icges(6,4) * t208 + Icges(6,2) * t902;
t102 = Icges(6,1) * t208 + Icges(6,4) * t902;
t15 = -t100 * t509 - t101 * t304 + t102 * t305 + t116 * t215 + t117 * t218 - t212 * t692;
t16 = t100 * t511 - t101 * t302 + t102 * t303 + t118 * t215 + t119 * t218 - t212 * t693;
t69 = Icges(6,4) * t119 + Icges(6,2) * t118 - Icges(6,6) * t693;
t71 = Icges(6,1) * t119 + Icges(6,4) * t118 - Icges(6,5) * t693;
t18 = t141 * t902 - t145 * t208 + t353 * t71 - t567 * t69;
t19 = t143 * t902 + t146 * t208 + t353 * t70 - t567 * t68;
t682 = qJD(1) * qJD(2);
t393 = qJDD(2) * t509 + t511 * t682;
t680 = qJD(1) * qJD(5);
t298 = -qJDD(5) * t509 - t511 * t680 + t393;
t468 = t509 * t682;
t299 = -t509 * t680 + t468 + (-qJDD(2) + qJDD(5)) * t511;
t48 = t138 * t511 - t141 * t302 - t145 * t303;
t49 = t511 * t140 - t302 * t143 + t303 * t146;
t531 = qJD(1) * (Icges(6,1) * t567 + t215 + t341) + t408 * (Icges(6,1) * t304 + t143 + t263) - t409 * (Icges(6,1) * t302 + t141 + t789);
t547 = qJD(1) * (-Icges(6,5) * t567 - Icges(6,6) * t353) + (Icges(6,5) * t302 + Icges(6,6) * t303) * t409 - (Icges(6,5) * t304 + Icges(6,6) * t305) * t408;
t64 = t212 * t511 - t215 * t302 + t218 * t303;
t67 = Icges(6,5) * t119 + Icges(6,6) * t118 - Icges(6,3) * t693;
t7 = t116 * t141 - t117 * t145 - t138 * t692 - t304 * t69 + t305 * t71 - t509 * t67;
t78 = -t141 * t567 - t145 * t353;
t79 = -t143 * t567 + t146 * t353;
t8 = t116 * t143 + t117 * t146 - t140 * t692 - t304 * t68 + t305 * t70 - t509 * t66;
t814 = -qJD(1) / 0.2e1;
t828 = t409 / 0.2e1;
t841 = qJD(1) * (Icges(6,2) * t353 - t218 + t340) - t408 * (-Icges(6,2) * t305 + t146 - t262) + t409 * (-Icges(6,2) * t303 - t145 - t790);
t886 = t50 / 0.2e1 + t49 / 0.2e1;
t9 = t118 * t141 - t119 * t145 - t138 * t693 - t302 * t69 + t303 * t71 + t511 * t67;
t931 = (-t18 * t511 + t19 * t509 + (t509 * t78 + t511 * t79) * qJD(1)) * t814 - qJDD(1) * (t509 * t79 - t511 * t78) / 0.2e1 - t509 * (qJD(1) * t15 + qJDD(1) * t65 + t408 * t8 - t409 * t7) / 0.2e1 + t511 * (qJD(1) * t16 + qJDD(1) * t64 + t10 * t408 - t409 * t9) / 0.2e1 - (qJD(1) * t64 + t408 * t49 - t409 * t48) * t693 / 0.2e1 - t14 * t692 / 0.2e1 + (t48 * t511 - t509 * t886) * t299 + (t509 * t624 + t511 * t886) * t298 - (t304 * t841 - t305 * t531 + (-qJD(1) * t624 - t7) * t511 + (qJD(1) * t50 - t547 + t8) * t509) * t408 / 0.2e1 + (t302 * t841 - t303 * t531 + (qJD(1) * t48 + t10) * t509 + (qJD(1) * t49 + t547 - t9) * t511) * t828;
t930 = t511 ^ 2;
t929 = rSges(3,2) * t508;
t928 = -t964 * t511 + (-Icges(5,1) * t755 + t373 * t511 + t430 - t963) * t509;
t927 = t969 - t960;
t926 = t957 - t967;
t925 = (Icges(4,2) * t754 + t431 + t962) * t511 + (-t371 * t511 + t961) * t509;
t387 = t599 * qJD(2);
t388 = t415 * qJD(2);
t924 = t953 * qJD(1) + t936 * qJD(2) - t387 * t508 + t388 * t510 + t942 * t476 + t941 * t477;
t560 = qJD(2) * t412;
t201 = -t511 * t560 + (-t509 * t599 + t779) * qJD(1);
t563 = qJD(2) * t414;
t203 = -t511 * t563 + (-t415 * t509 + t788) * qJD(1);
t923 = -qJD(2) * t879 - t201 * t508 + t203 * t510 + t476 * t947 + t477 * t945 + t932;
t202 = qJD(1) * t315 - t509 * t560;
t204 = qJD(1) * t317 - t509 * t563;
t861 = qJD(1) * t268;
t922 = qJD(1) * t943 - qJD(2) * t878 + t202 * t508 - t204 * t510 - t476 * t946 + t477 * t944 - t861;
t921 = t509 * t887 - t511 * t888;
t920 = t509 * t885 - t889 * t511;
t919 = t937 * qJD(1) + qJD(2) * t939;
t918 = qJD(1) * t934 - t509 * t938 + t932;
t917 = -t861 - t938 * t511 + (-t509 * t971 - t935 + t968) * qJD(1);
t916 = -t531 * t353 + t567 * t841;
t381 = rSges(4,1) * t476 + rSges(4,2) * t477;
t332 = t381 * t509;
t493 = t509 * rSges(4,3);
t464 = t477 * rSges(4,1);
t865 = -rSges(4,2) * t476 + t464;
t179 = -qJD(2) * t332 + (t511 * t865 + t493) * qJD(1);
t355 = t865 * qJD(2);
t394 = -qJDD(2) * t511 + t468;
t681 = qJD(1) * qJD(3);
t822 = pkin(2) * t508;
t665 = qJDD(3) * t509 + t394 * t822 + t511 * t681;
t293 = rSges(4,1) * t754 - rSges(4,2) * t756 - t511 * rSges(4,3);
t501 = t511 * pkin(6);
t438 = pkin(1) * t509 - t501;
t506 = -qJ(3) - pkin(6);
t470 = t511 * t506;
t500 = t510 * pkin(2);
t471 = t500 + pkin(1);
t704 = -t509 * t471 - t470;
t256 = t438 + t704;
t728 = t256 - t438;
t670 = -t293 + t728;
t512 = qJD(2) ^ 2;
t748 = t510 * t512;
t674 = pkin(2) * t748;
t499 = t509 * pkin(6);
t452 = t506 * t693;
t676 = pkin(2) * t752;
t702 = qJD(2) * t676 + qJD(3) * t511;
t663 = t452 + t702;
t812 = pkin(1) - t471;
t198 = (-t511 * t812 - t499) * qJD(1) - t663;
t439 = t511 * pkin(1) + t499;
t392 = t439 * qJD(1);
t738 = -t198 - t392;
t42 = t381 * t394 + (-qJD(2) * t355 - t674) * t511 + t670 * qJDD(1) + (-t179 + t738) * qJD(1) + t665;
t913 = t42 - g(1);
t222 = rSges(6,1) * t567 + t353 * rSges(6,2);
t912 = t222 * t408;
t911 = t222 * t409;
t908 = t883 * qJD(1);
t907 = t884 * qJD(1);
t181 = -t302 * rSges(6,1) - t303 * rSges(6,2);
t183 = -t304 * rSges(6,1) - t305 * rSges(6,2);
t906 = -t181 * t408 - t183 * t409;
t656 = t477 * t689;
t400 = qJ(4) * t656;
t686 = qJD(4) * t511;
t419 = t476 * t686;
t551 = -t476 * t689 - t477 * t693;
t660 = t476 * t693;
t150 = pkin(3) * t551 - qJ(4) * t660 + t400 + t419;
t690 = qJD(2) * t509;
t310 = t476 * t692 + t477 * t690;
t657 = t476 * t690;
t406 = pkin(3) * t657;
t442 = pkin(3) * t753;
t687 = qJD(4) * t509;
t651 = t476 * t687;
t151 = qJ(4) * t310 + qJD(1) * t442 - t406 + t651;
t425 = qJ(4) * t754;
t330 = -pkin(3) * t756 + t425;
t459 = t476 * qJ(4);
t864 = t477 * pkin(3) + t459;
t333 = t864 * t509;
t427 = qJ(4) * t753;
t334 = -pkin(3) * t755 + t427;
t458 = qJD(4) * t476;
t475 = pkin(6) * t692;
t481 = qJD(3) * t509;
t655 = t508 * t689;
t197 = -pkin(2) * t655 - t475 + t481 + (t509 * t812 - t470) * qJD(1);
t671 = t511 * t197 + t509 * t198 - t256 * t692;
t901 = t511 * t150 + t509 * t151 - t330 * t690 + t333 * t692 - t334 * t689 - t458 + t671;
t900 = t943 + t950;
t679 = qJD(2) * qJD(4);
t899 = qJDD(4) * t476 + t477 * t679;
t379 = pkin(3) * t476 - qJ(4) * t477;
t897 = t379 * t693 - t477 * t686;
t498 = t511 * rSges(5,2);
t866 = t477 * rSges(5,1) + t476 * rSges(5,3);
t292 = t509 * t866 - t498;
t896 = qJD(1) * t292;
t895 = qJD(2) * t920 + t907;
t894 = qJD(2) * t921 + t908;
t893 = t509 * t919 + t511 * t924;
t892 = t509 * t924 - t511 * t919;
t891 = qJD(2) * t934 - t202 * t510 - t204 * t508 + t476 * t944 + t477 * t946;
t890 = qJD(2) * t935 + t201 * t510 + t203 * t508 + t476 * t945 - t477 * t947;
t541 = t314 * t511 - t315 * t509;
t844 = t509 * (-t412 * t511 + t317) - t511 * (-Icges(3,2) * t750 + t316 - t455);
t882 = -t476 * t925 + t477 * t928 - t508 * t844 + t541 * t510;
t707 = t414 + t599;
t708 = -t412 + t415;
t881 = (-t476 * t927 + t477 * t926 - t508 * t707 + t510 * t708) * qJD(1);
t880 = t746 + t949;
t877 = t918 * t930 + (t923 * t509 + (-t917 + t922) * t511) * t509;
t876 = t922 * t930 + (t917 * t509 + (-t918 + t923) * t511) * t509;
t875 = t939 * qJD(1);
t552 = -t381 - t822;
t869 = t511 * t552;
t396 = qJD(1) * t438;
t868 = qJD(1) * t256 - t396;
t867 = t865 + t500;
t441 = pkin(4) * t753;
t378 = -pkin(7) * t509 + t441;
t662 = t500 + t864;
t221 = rSges(6,1) * t353 - rSges(6,2) * t567;
t642 = -t379 - t822;
t821 = pkin(4) * t476;
t573 = t642 - t821;
t863 = -t409 * t221 + t573 * t689;
t669 = -t333 + t728;
t604 = t669 - t740;
t706 = t419 + t481;
t46 = qJD(1) * t604 + t706 + t863;
t344 = t379 * t690;
t405 = pkin(4) * t657;
t570 = t651 - t702;
t721 = t305 * rSges(6,1) - t304 * rSges(6,2);
t149 = -rSges(6,3) * t509 + t721;
t444 = t511 * t471;
t636 = -t506 * t509 + t444;
t257 = t636 - t439;
t337 = qJ(4) * t755 + t442;
t726 = -t257 - t337;
t667 = -t378 + t726;
t635 = -t149 + t667;
t47 = -t221 * t408 - t344 - t405 + (t439 - t635) * qJD(1) + t570;
t862 = t46 * t511 + t47 * t509;
t815 = g(2) * t509;
t854 = (g(1) * t511 + t815) * t476;
t839 = m(5) / 0.2e1;
t838 = m(6) / 0.2e1;
t837 = -m(5) - m(6);
t836 = -pkin(3) - pkin(4);
t833 = t393 / 0.2e1;
t832 = t394 / 0.2e1;
t827 = t509 / 0.2e1;
t826 = -t511 / 0.2e1;
t825 = -rSges(5,1) - pkin(3);
t824 = rSges(6,3) + pkin(7);
t820 = pkin(4) * t477;
t817 = g(1) * t509;
t811 = rSges(3,1) * t510;
t805 = t476 * rSges(5,1);
t495 = t509 * rSges(5,2);
t494 = t509 * rSges(3,3);
t550 = qJD(2) * t869 + t481;
t90 = qJD(1) * t670 + t550;
t804 = t90 * t381;
t802 = -rSges(5,3) - qJ(4);
t247 = pkin(4) * t551 - pkin(7) * t692;
t745 = t117 * rSges(6,1) + t116 * rSges(6,2);
t72 = -rSges(6,3) * t692 + t745;
t801 = t247 + t72;
t248 = qJD(1) * t378 - t405;
t619 = rSges(6,1) * t119 + rSges(6,2) * t118;
t73 = -rSges(6,3) * t693 + t619;
t800 = t248 + t73;
t421 = rSges(3,1) * t508 + rSges(3,2) * t510;
t658 = t421 * t689;
t701 = rSges(3,2) * t752 + t511 * rSges(3,3);
t338 = rSges(3,1) * t750 - t701;
t717 = -t338 - t438;
t188 = qJD(1) * t717 - t658;
t771 = t188 * t509;
t770 = t188 * t511;
t339 = rSges(3,1) * t749 - rSges(3,2) * t751 + t494;
t242 = t339 + t439;
t189 = qJD(1) * t242 - t421 * t690;
t364 = t421 * t511;
t769 = t189 * t364;
t739 = t149 + t378;
t732 = -t256 * t690 + t257 * t689;
t731 = -t509 * t256 + t511 * t257;
t295 = rSges(4,1) * t753 - rSges(4,2) * t755 + t493;
t727 = -t257 - t295;
t720 = qJD(1) * t334 + t477 * t687;
t711 = qJD(1) * (-pkin(1) * t693 + t475) + qJDD(1) * t439;
t710 = rSges(5,2) * t692 + rSges(5,3) * t656;
t709 = rSges(4,2) * t660 + rSges(4,3) * t692;
t703 = rSges(3,3) * t692 + t693 * t929;
t700 = t509 ^ 2 + t930;
t691 = qJD(2) * t477;
t688 = qJD(4) * t477;
t677 = -t506 - t824;
t675 = pkin(2) * t751;
t673 = qJD(2) * t500;
t672 = t197 * t689 + t198 * t690 - t393 * t256;
t294 = rSges(5,1) * t753 + rSges(5,3) * t755 + t495;
t668 = -t294 + t726;
t664 = t400 + t706;
t654 = t510 * t689;
t650 = -pkin(1) - t811;
t646 = -t690 / 0.2e1;
t645 = t690 / 0.2e1;
t644 = -t689 / 0.2e1;
t643 = t689 / 0.2e1;
t633 = t509 * t333 + t511 * t337 + t731;
t632 = -t292 + t669;
t631 = t700 * t822;
t630 = t425 - t676;
t629 = t427 - t675;
t380 = -rSges(5,3) * t477 + t805;
t628 = -t380 + t642;
t291 = qJD(2) * t864 - t688;
t626 = -t291 - t673;
t625 = -t500 - t820;
t424 = rSges(2,1) * t511 - rSges(2,2) * t509;
t422 = rSges(2,1) * t509 + rSges(2,2) * t511;
t423 = t811 - t929;
t103 = rSges(6,1) * t208 + rSges(6,2) * t902;
t546 = qJD(1) * t197 + qJDD(1) * t257 - qJDD(3) * t511 + t509 * t681 + t711;
t528 = qJDD(1) * t337 + t546 + t899 * t509 + (t150 + t419) * qJD(1);
t538 = -qJD(2) * t291 + t512 * t625;
t11 = qJD(1) * t801 + qJDD(1) * t739 - t408 * t103 - t298 * t221 + t393 * t573 + t509 * t538 + t528;
t565 = -t151 - t651 + t738;
t568 = t394 * t379 + t511 * t899 + t665;
t12 = t394 * t821 - t409 * t103 + t299 * t221 + t538 * t511 + t604 * qJDD(1) + (t565 - t800) * qJD(1) + t568;
t616 = t11 * t509 + t12 * t511;
t605 = -qJD(1) * t333 + t706 + t868;
t177 = rSges(4,1) * t551 - rSges(4,2) * t656 + t709;
t593 = t177 * t511 + t179 * t509;
t592 = -t189 * t509 - t770;
t205 = -rSges(3,2) * t654 + (-t510 * t693 - t655) * rSges(3,1) + t703;
t363 = t421 * t509;
t206 = -qJD(2) * t363 + (t423 * t511 + t494) * qJD(1);
t591 = t205 * t511 + t206 * t509;
t584 = t293 * t509 + t295 * t511;
t581 = t338 * t509 + t339 * t511;
t354 = t866 * qJD(2);
t574 = -t354 + t626;
t566 = -t221 + t573;
t564 = t333 * t690 + t337 * t689 - t688 + t732;
t554 = t628 * t689;
t553 = -t864 - t820;
t549 = -t674 + (-t291 - t354) * qJD(2);
t545 = -t471 + t553;
t539 = -pkin(4) * t691 - t103 + t626;
t537 = -qJDD(4) * t477 + t150 * t689 + t151 * t690 + t393 * t333 + t476 * t679 + t672;
t533 = -t471 - t864 - t866;
t435 = rSges(5,3) * t753;
t433 = rSges(5,3) * t754;
t389 = t423 * qJD(2);
t336 = t381 * t511;
t335 = -rSges(5,1) * t755 + t435;
t331 = -rSges(5,1) * t756 + t433;
t311 = t656 - t660;
t309 = t700 * t476 * qJD(2);
t184 = t581 * qJD(2);
t178 = -t380 * t690 + (t511 * t866 + t495) * qJD(1);
t176 = rSges(5,1) * t551 - rSges(5,3) * t660 + t710;
t91 = -t381 * t690 + (t439 - t727) * qJD(1) - t702;
t86 = qJD(2) * t584 + t732;
t85 = qJD(1) * t205 + qJDD(1) * t339 - t389 * t690 - t393 * t421 + t711;
t84 = -t389 * t689 + t394 * t421 + t717 * qJDD(1) + (-t206 - t392) * qJD(1);
t77 = -t344 + (-qJD(2) * t380 + t458) * t509 + (t439 - t668) * qJD(1) - t702;
t76 = qJD(1) * t632 + t554 + t706;
t62 = (t292 * t509 + t294 * t511) * qJD(2) + t564;
t43 = -t355 * t690 + qJD(1) * t177 + qJDD(1) * t295 - t381 * t393 + (-t393 * t508 - t509 * t748) * pkin(2) + t546;
t33 = t147 * t408 + t149 * t409 + (t377 * t509 + t378 * t511) * qJD(2) + t564;
t23 = qJD(1) * t176 + qJDD(1) * t294 + t393 * t628 + t509 * t549 + t528;
t22 = t380 * t394 + t549 * t511 + t632 * qJDD(1) + (-t178 + t565) * qJD(1) + t568;
t17 = t292 * t393 + (t176 * t511 + t178 * t509) * qJD(2) + t668 * t394 + t537;
t6 = t147 * t298 - t149 * t299 + t377 * t393 + t408 * t73 + t409 * t72 + (t247 * t511 + t248 * t509) * qJD(2) + t667 * t394 + t537;
t1 = [t14 * t828 - m(2) * (-g(1) * t422 + g(2) * t424) + (t79 + t65) * t298 / 0.2e1 + (t78 + t64) * t299 / 0.2e1 + (t15 + t19) * t408 / 0.2e1 - (t16 + t14 + t18) * t409 / 0.2e1 + (((t92 - t870 + t880) * t509 + ((t952 + t966) * t511 + t905 + t948 + t965) * t511) * qJD(2) + t908) * t643 + (-t937 * qJD(2) - t101 * t567 + t102 * t353 + t208 * t218 + t902 * t215 + t387 * t510 + t388 * t508 + t941 * t476 - t942 * t477) * qJD(1) + (-(-t46 + t605 + t863 - t933) * t47 + t12 * (-t618 + t704) + t46 * (t405 + t406 - t619 + t663) + t47 * (t664 + t745) + (-t12 * t824 + t47 * (t476 * t836 - t822) * qJD(2)) * t511 + (t12 * t553 + t46 * (-qJ(4) * t691 - t458)) * t509 + ((t46 * t824 + t47 * t545) * t509 + (t46 * t545 + t47 * t677) * t511) * qJD(1) - g(1) * (-t147 + t704 - t818) - (t477 * t836 - t459) * t817 + (t11 - g(2)) * (t509 * t677 + t337 + t441 + t444 + t721)) * m(6) + (-(t554 + t605 - t76 - t896) * t77 + t76 * (t406 + t452 - t570) + t77 * (t664 + t710) + (t77 * (t476 * t825 - t822) * t511 + t76 * (t477 * t802 + t805) * t509) * qJD(2) + ((-t77 * t506 + t533 * t76) * t511 + (-t76 * rSges(5,2) + t533 * t77) * t509) * qJD(1) + (t23 - g(2)) * (t636 + t294 + t337) + (t22 - g(1)) * (t498 + (t476 * t802 + t477 * t825) * t509 + t704)) * m(5) + (-(-qJD(1) * t293 + t550 + t868 - t90) * t91 + t90 * t663 + t91 * (t481 + t709) + (t509 * t804 + t869 * t91) * qJD(2) + ((-t90 * rSges(4,3) + t91 * (-t471 - t464)) * t509 + (t90 * (-t471 - t865) - t91 * t506) * t511) * qJD(1) + (t43 - g(2)) * (t295 + t636) + t913 * (-t293 + t704)) * m(4) + (t189 * (t475 + t703) + (t421 * t771 - t769) * qJD(2) + ((-pkin(1) - t423) * t770 + (t188 * (-rSges(3,3) - pkin(6)) + t189 * t650) * t509) * qJD(1) - (-qJD(1) * t338 - t188 - t396 - t658) * t189 + (-g(2) + t85) * t242 + (-g(1) + t84) * (t650 * t509 + t501 + t701)) * m(3) + (t879 + t883) * t833 + (-t878 + t884) * t832 + (t890 + t893) * t645 + (-t215 * t567 + t218 * t353 + m(2) * (t422 ^ 2 + t424 ^ 2) + Icges(2,3) - t936) * qJDD(1) + (((t511 * t900 - t880 + t887) * t511 + (t509 * t900 - t245 + t623 + t888 - t951) * t509) * qJD(2) + t895 - t907) * t646 + (-t891 + t892 + t894) * t644; t921 * t833 + t920 * t832 + (qJD(1) * t893 + qJD(2) * t876 + qJDD(1) * t883 + t393 * t887 + t394 * t888) * t827 + (qJD(1) * t892 + qJD(2) * t877 + qJDD(1) * t884 + t393 * t885 + t394 * t889) * t826 + (t891 * t511 + t890 * t509 + (-t878 * t509 + t879 * t511) * qJD(1)) * qJD(1) / 0.2e1 + (t879 * t509 + t878 * t511) * qJDD(1) / 0.2e1 + t895 * t693 / 0.2e1 + t894 * t692 / 0.2e1 + ((-t690 * t874 + t875) * t509 + ((t509 * t873 + t882) * qJD(2) + t881) * t511) * t646 + ((t509 * t888 + t511 * t887) * qJD(1) + t876) * t645 + ((t509 * t889 + t511 * t885) * qJD(1) + t877) * t644 + ((-t689 * t873 - t875) * t511 + ((t511 * t874 + t882) * qJD(2) + t881) * t509) * t643 + (-t862 * (-t864 + t625) * qJD(2) - g(1) * (t629 - t183) - g(2) * (t630 - t181) - g(3) * (t222 + t662 + t820) - t836 * t854 + t6 * t633 + (t11 * t566 + t6 * t740) * t509 + (t12 * t566 + t6 * t739) * t511 + (t912 + t539 * t509 - t720 + (pkin(4) * t755 + t511 * t566 + t183 + t675) * qJD(1)) * t47 + (t911 + t539 * t511 + (t221 * t509 - t181 + t330) * qJD(1) + t897) * t46 + (-(-t700 * t821 - t631) * qJD(2) + (qJD(1) * t635 + t800) * t509 + (t801 + t933) * t511 + t901 - t906) * t33) * m(6) + (-g(1) * (t435 + t629) - g(2) * (t433 + t630) - t825 * t854 + t17 * t633 + (t17 * t294 + t22 * t628) * t511 + (t17 * t292 + t23 * t628) * t509 + (t574 * t509 - t720 + (t511 * t628 - t335 + t675) * qJD(1)) * t77 + (t574 * t511 + (t380 * t509 + t330 + t331) * qJD(1) + t897) * t76 + (-(t331 * t509 + t335 * t511 - t631) * qJD(2) + (t176 + t896) * t511 + (qJD(1) * t668 + t178) * t509 + t901) * t62 + (-g(3) + (t77 * t509 + t76 * t511) * qJD(2)) * (t662 + t866)) * m(5) + (-g(3) * t867 - t552 * t815 + t90 * (-pkin(2) * t654 - t355 * t511) + (qJD(2) * t593 + t293 * t393 + t394 * t727 + t672) * (t584 + t731) + t86 * (t593 + t671) + (t86 * t293 + t552 * t91) * t692 + (t43 * t552 + t91 * (-t355 - t673) + (t727 * t86 + t804) * qJD(1)) * t509 - (t90 * t332 + t91 * (-t336 - t675)) * qJD(1) - (-t86 * t631 + (-t86 * t336 - t867 * t90) * t511 + (-t86 * t332 - t867 * t91) * t509) * qJD(2) + t913 * t869) * m(4) + (g(1) * t364 + g(2) * t363 - g(3) * t423 - (t188 * t363 - t769) * qJD(1) - (t184 * (-t363 * t509 - t364 * t511) + t592 * t423) * qJD(2) + (qJD(2) * t591 + t338 * t393 - t339 * t394) * t581 + t184 * ((t338 * t511 - t339 * t509) * qJD(1) + t591) + t592 * t389 + (-t85 * t509 - t84 * t511 + (-t189 * t511 + t771) * qJD(1)) * t421) * m(3) + ((t476 * t928 + t477 * t925 + t541 * t508 + t510 * t844) * qJD(2) + (t476 * t926 + t477 * t927 + t508 * t708 + t510 * t707) * qJD(1) - t916) * t814 - t931; (-m(4) + t837) * (-g(2) * t511 + t817) + 0.2e1 * (t11 * t826 + t12 * t827) * m(6) + 0.2e1 * (t22 * t827 + t23 * t826) * m(5) + 0.2e1 * (t42 * t827 + t43 * t826) * m(4); t837 * (-g(3) * t477 + t854) - m(5) * (t309 * t62 + t310 * t77 + t311 * t76) - m(6) * (t309 * t33 + t310 * t47 + t311 * t46) + 0.2e1 * ((t689 * t76 + t690 * t77 - t17) * t839 + (t46 * t689 + t47 * t690 - t6) * t838) * t477 + 0.2e1 * ((qJD(2) * t62 + t22 * t511 + t23 * t509 + t692 * t77 - t693 * t76) * t839 + (qJD(2) * t33 - t46 * t693 + t47 * t692 + t616) * t838) * t476; t916 * t814 + (-g(1) * t183 - g(2) * t181 + g(3) * t222 + t6 * (-t147 * t509 - t149 * t511) + t862 * t103 + ((-t46 * t509 + t47 * t511) * qJD(1) + t616) * t221 - t46 * (-qJD(1) * t181 + t911) - t47 * (qJD(1) * t183 + t912) + (-t509 * t73 - t511 * t72 + (-t147 * t511 + t149 * t509) * qJD(1) + t906) * t33) * m(6) + t931;];
tau = t1;

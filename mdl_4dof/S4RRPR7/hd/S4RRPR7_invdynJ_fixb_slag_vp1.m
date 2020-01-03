% Calculate vector of inverse dynamics joint torques for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR7_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:57
% EndTime: 2019-12-31 17:06:43
% DurationCPUTime: 40.28s
% Computational Cost: add. (14471->875), mult. (22084->1169), div. (0->0), fcn. (20186->8), ass. (0->419)
t724 = Icges(3,3) + Icges(4,3);
t364 = qJ(2) + pkin(7);
t343 = sin(t364);
t344 = cos(t364);
t369 = sin(qJ(2));
t372 = cos(qJ(2));
t712 = Icges(3,5) * t372 + Icges(4,5) * t344 - Icges(3,6) * t369 - Icges(4,6) * t343;
t373 = cos(qJ(1));
t723 = t724 * t373;
t370 = sin(qJ(1));
t565 = t370 * t372;
t568 = t369 * t370;
t571 = t344 * t370;
t573 = t343 * t370;
t713 = -Icges(3,5) * t565 - Icges(4,5) * t571 + Icges(3,6) * t568 + Icges(4,6) * t573 + t723;
t718 = t724 * t370 + t712 * t373;
t589 = Icges(4,6) * t373;
t202 = Icges(4,4) * t571 - Icges(4,2) * t573 - t589;
t590 = Icges(3,6) * t373;
t212 = Icges(3,4) * t565 - Icges(3,2) * t568 - t590;
t722 = t202 * t343 + t212 * t369;
t302 = Icges(4,4) * t573;
t594 = Icges(4,5) * t373;
t204 = Icges(4,1) * t571 - t302 - t594;
t325 = Icges(3,4) * t568;
t595 = Icges(3,5) * t373;
t214 = Icges(3,1) * t565 - t325 - t595;
t701 = -t204 * t344 - t214 * t372 + t722;
t683 = -t370 * t701 + t713 * t373;
t599 = Icges(4,4) * t343;
t267 = Icges(4,1) * t344 - t599;
t205 = Icges(4,5) * t370 + t267 * t373;
t600 = Icges(3,4) * t369;
t294 = Icges(3,1) * t372 - t600;
t215 = Icges(3,5) * t370 + t294 * t373;
t720 = -t205 * t571 - t215 * t565;
t719 = Icges(3,5) * t369 + Icges(4,5) * t343 + Icges(3,6) * t372 + Icges(4,6) * t344;
t264 = Icges(4,2) * t344 + t599;
t328 = Icges(4,4) * t344;
t266 = Icges(4,1) * t343 + t328;
t291 = Icges(3,2) * t372 + t600;
t355 = Icges(3,4) * t372;
t293 = Icges(3,1) * t369 + t355;
t710 = t264 * t343 - t266 * t344 + t291 * t369 - t293 * t372;
t717 = t718 * t373 + t720;
t563 = t372 * t373;
t570 = t344 * t373;
t655 = -t205 * t570 - t215 * t563 - t718 * t370;
t716 = -t204 * t570 - t214 * t563 + t713 * t370;
t440 = -Icges(4,2) * t343 + t328;
t203 = Icges(4,6) * t370 + t373 * t440;
t441 = -Icges(3,2) * t369 + t355;
t213 = Icges(3,6) * t370 + t373 * t441;
t714 = t203 * t343 + t213 * t369;
t682 = -t203 * t573 - t213 * t568 - t717;
t567 = t369 * t373;
t572 = t343 * t373;
t681 = -t202 * t572 - t212 * t567 - t716;
t680 = -t203 * t572 - t213 * t567 - t655;
t677 = t202 * t344 + t204 * t343 + t212 * t372 + t214 * t369;
t676 = t203 * t344 + t205 * t343 + t213 * t372 + t215 * t369;
t711 = t719 * qJD(2);
t709 = -t264 * t344 - t266 * t343 - t291 * t372 - t293 * t369;
t708 = t205 * t344 + t215 * t372 - t714;
t664 = t719 * t373;
t663 = t719 * t370;
t675 = -t370 * t710 - t664;
t674 = -t373 * t710 + t663;
t707 = t718 * qJD(1);
t371 = cos(qJ(4));
t564 = t371 * t373;
t368 = sin(qJ(4));
t569 = t368 * t370;
t242 = t344 * t569 + t564;
t559 = t373 * t368;
t566 = t370 * t371;
t243 = t344 * t566 - t559;
t460 = t243 * rSges(5,1) - t242 * rSges(5,2);
t123 = -rSges(5,3) * t573 - t460;
t613 = rSges(5,2) * t368;
t615 = rSges(5,1) * t371;
t459 = -t613 + t615;
t193 = -rSges(5,3) * t344 + t343 * t459;
t513 = qJD(4) * t343;
t514 = qJD(2) * t373;
t259 = -t370 * t513 + t514;
t512 = qJD(4) * t344;
t312 = qJD(1) - t512;
t347 = qJD(3) * t370;
t623 = t343 * pkin(3);
t271 = -pkin(6) * t344 + t623;
t626 = pkin(2) * t369;
t474 = -t271 - t626;
t706 = t312 * t123 - t259 * t193 + t474 * t514 + t347;
t113 = Icges(5,5) * t243 - Icges(5,6) * t242 + Icges(5,3) * t573;
t226 = Icges(5,4) * t243;
t116 = -Icges(5,2) * t242 + Icges(5,6) * t573 + t226;
t225 = Icges(5,4) * t242;
t120 = -Icges(5,1) * t243 - Icges(5,5) * t573 + t225;
t690 = t116 * t368 + t120 * t371;
t46 = -t113 * t344 - t343 * t690;
t705 = t373 ^ 2;
t247 = t440 * qJD(2);
t248 = t267 * qJD(2);
t274 = t441 * qJD(2);
t275 = t294 * qJD(2);
t704 = t719 * qJD(1) + t709 * qJD(2) - t247 * t343 + t248 * t344 - t274 * t369 + t275 * t372;
t408 = qJD(2) * t264;
t107 = -t373 * t408 + (-t370 * t440 + t589) * qJD(1);
t410 = qJD(2) * t266;
t109 = -t373 * t410 + (-t267 * t370 + t594) * qJD(1);
t409 = qJD(2) * t291;
t146 = -t373 * t409 + (-t370 * t441 + t590) * qJD(1);
t411 = qJD(2) * t293;
t148 = -t373 * t411 + (-t294 * t370 + t595) * qJD(1);
t703 = -t676 * qJD(2) - t107 * t343 + t109 * t344 - t146 * t369 + t148 * t372 + t707;
t108 = qJD(1) * t203 - t370 * t408;
t110 = qJD(1) * t205 - t370 * t410;
t147 = qJD(1) * t213 - t370 * t409;
t149 = qJD(1) * t215 - t370 * t411;
t702 = t713 * qJD(1) + t677 * qJD(2) + t108 * t343 - t110 * t344 + t147 * t369 - t149 * t372;
t700 = t680 * t370 - t681 * t373;
t699 = t682 * t370 - t683 * t373;
t698 = t710 * qJD(1) + t712 * qJD(2);
t697 = t701 * qJD(1) - t711 * t370 + t707;
t696 = -t711 * t373 + (-t712 * t370 - t708 + t723) * qJD(1);
t268 = rSges(4,1) * t343 + rSges(4,2) * t344;
t229 = t268 * t370;
t356 = t370 * rSges(4,3);
t330 = t344 * rSges(4,1);
t658 = -rSges(4,2) * t343 + t330;
t112 = -qJD(2) * t229 + (t373 * t658 + t356) * qJD(1);
t249 = t658 * qJD(2);
t510 = qJD(1) * qJD(2);
t283 = -qJDD(2) * t373 + t370 * t510;
t509 = qJD(1) * qJD(3);
t495 = qJDD(3) * t370 + t283 * t626 + t373 * t509;
t206 = rSges(4,1) * t571 - rSges(4,2) * t573 - t373 * rSges(4,3);
t362 = t373 * pkin(5);
t309 = pkin(1) * t370 - t362;
t367 = -qJ(3) - pkin(5);
t337 = t373 * t367;
t361 = t372 * pkin(2);
t339 = t361 + pkin(1);
t526 = -t370 * t339 - t337;
t198 = t309 + t526;
t544 = t198 - t309;
t497 = -t206 + t544;
t562 = t372 * qJD(2) ^ 2;
t505 = pkin(2) * t562;
t360 = t370 * pkin(5);
t518 = qJD(1) * t370;
t507 = pkin(2) * t568;
t524 = qJD(2) * t507 + qJD(3) * t373;
t494 = t367 * t518 + t524;
t618 = pkin(1) - t339;
t141 = (-t373 * t618 - t360) * qJD(1) - t494;
t310 = t373 * pkin(1) + t360;
t279 = t310 * qJD(1);
t553 = -t141 - t279;
t29 = t268 * t283 + (-qJD(2) * t249 - t505) * t373 + t497 * qJDD(1) + (-t112 + t553) * qJD(1) + t495;
t695 = t29 - g(1);
t694 = t674 * qJD(1);
t329 = t343 * rSges(5,3);
t194 = t344 * t459 + t329;
t692 = t675 * qJD(1);
t332 = t344 * pkin(3);
t657 = t343 * pkin(6) + t332;
t238 = t657 * t370;
t496 = -t238 + t544;
t40 = qJD(1) * t496 + t706;
t244 = -t344 * t559 + t566;
t245 = t344 * t564 + t569;
t124 = t245 * rSges(5,1) + t244 * rSges(5,2) + rSges(5,3) * t572;
t240 = pkin(3) * t570 + pkin(6) * t572;
t515 = qJD(2) * t370;
t258 = t373 * t513 + t515;
t468 = t373 * t339 - t367 * t370;
t199 = t468 - t310;
t541 = t199 + t310;
t41 = -t271 * t515 + t124 * t312 - t193 * t258 + (t240 + t541) * qJD(1) - t524;
t607 = t370 * t41;
t691 = t373 * t40 + t607;
t36 = t113 * t573 - t116 * t242 - t120 * t243;
t115 = Icges(5,5) * t245 + Icges(5,6) * t244 + Icges(5,3) * t572;
t598 = Icges(5,4) * t245;
t118 = Icges(5,2) * t244 + Icges(5,6) * t572 + t598;
t227 = Icges(5,4) * t244;
t121 = Icges(5,1) * t245 + Icges(5,5) * t572 + t227;
t37 = t115 * t573 - t242 * t118 + t243 * t121;
t438 = Icges(5,5) * t371 - Icges(5,6) * t368;
t186 = -Icges(5,3) * t344 + t343 * t438;
t596 = Icges(5,4) * t371;
t439 = -Icges(5,2) * t368 + t596;
t188 = -Icges(5,6) * t344 + t343 * t439;
t597 = Icges(5,4) * t368;
t442 = Icges(5,1) * t371 - t597;
t190 = -Icges(5,5) * t344 + t343 * t442;
t63 = t186 * t573 - t188 * t242 + t190 * t243;
t12 = t258 * t37 - t259 * t36 + t312 * t63;
t38 = t113 * t572 + t244 * t116 - t120 * t245;
t39 = t115 * t572 + t244 * t118 + t245 * t121;
t64 = t186 * t572 + t188 * t244 + t190 * t245;
t13 = t258 * t39 - t259 * t38 + t64 * t312;
t689 = qJD(2) * t699 + t692;
t688 = qJD(2) * t700 + t694;
t687 = t370 * t698 + t373 * t704;
t686 = t370 * t704 - t373 * t698;
t685 = qJD(2) * t701 - t108 * t344 - t110 * t343 - t147 * t372 - t149 * t369;
t684 = t708 * qJD(2) + t107 * t344 + t109 * t343 + t146 * t372 + t148 * t369;
t679 = rSges(3,2) * t369;
t397 = t212 * t373 - t213 * t370;
t398 = t202 * t373 - t203 * t370;
t643 = t370 * (-t264 * t373 + t205) - t373 * (-Icges(4,2) * t571 + t204 - t302);
t644 = t370 * (-t291 * t373 + t215) - t373 * (-Icges(3,2) * t565 + t214 - t325);
t669 = -t343 * t643 + t398 * t344 - t369 * t644 + t397 * t372;
t527 = t293 + t441;
t528 = -t291 + t294;
t533 = t266 + t440;
t534 = -t264 + t267;
t668 = (-t343 * t533 + t344 * t534 - t369 * t527 + t372 * t528) * qJD(1);
t667 = t697 * t705 + (t703 * t370 + (-t696 + t702) * t373) * t370;
t666 = t702 * t705 + (t696 * t370 + (-t697 + t703) * t373) * t370;
t665 = t712 * qJD(1);
t405 = -t268 - t626;
t661 = t373 * t405;
t285 = qJD(1) * t309;
t660 = qJD(1) * t198 - t285;
t659 = t658 + t361;
t656 = t713 + t714;
t228 = (-rSges(5,1) * t368 - rSges(5,2) * t371) * t343;
t104 = qJD(2) * t194 + qJD(4) * t228;
t517 = qJD(1) * t373;
t403 = t343 * t517 + t344 * t515;
t508 = qJDD(4) * t343;
t139 = qJD(4) * t403 + t370 * t508 + t283;
t486 = t343 * t515;
t143 = t403 * pkin(6) + (t344 * t517 - t486) * pkin(3);
t241 = qJD(2) * t513 - qJDD(4) * t344 + qJDD(1);
t250 = t657 * qJD(2);
t467 = qJD(1) * t344 - qJD(4);
t102 = t312 * t566 + (-t373 * t467 + t486) * t368;
t516 = qJD(2) * t343;
t103 = t467 * t564 + (t312 * t368 - t371 * t516) * t370;
t461 = rSges(5,1) * t103 + rSges(5,2) * t102;
t60 = rSges(5,3) * t403 + t461;
t10 = -t104 * t259 + t123 * t241 + t139 * t193 + t271 * t283 - t312 * t60 + (-qJD(2) * t250 - t505) * t373 + t496 * qJDD(1) + (-t143 + t553) * qJD(1) + t495;
t282 = qJDD(2) * t370 + t373 * t510;
t490 = t344 * t514;
t493 = t343 * t518;
t402 = t490 - t493;
t138 = qJD(4) * t402 + t373 * t508 + t282;
t288 = pkin(6) * t490;
t485 = t343 * t514;
t404 = -t344 * t518 - t485;
t142 = pkin(3) * t404 - pkin(6) * t493 + t288;
t342 = pkin(5) * t517;
t489 = t369 * t514;
t419 = -pkin(2) * t489 + t347;
t140 = -t342 + (t370 * t618 - t337) * qJD(1) + t419;
t532 = qJD(1) * (-pkin(1) * t518 + t342) + qJDD(1) * t310;
t379 = qJD(1) * t140 + qJDD(1) * t199 + t370 * t509 + (-t282 * t369 - t370 * t562) * pkin(2) - qJDD(3) * t373 + t532;
t420 = t312 * t373;
t651 = t370 * t467 + t485;
t100 = t368 * t651 + t371 * t420;
t101 = t368 * t420 - t371 * t651;
t504 = t101 * rSges(5,1) + t100 * rSges(5,2) + rSges(5,3) * t490;
t59 = -rSges(5,3) * t493 + t504;
t11 = qJD(1) * t142 + qJDD(1) * t240 - t104 * t258 + t124 * t241 - t138 * t193 - t250 * t515 - t271 * t282 + t312 * t59 + t379;
t653 = t10 * t370 - t11 * t373;
t625 = g(1) * t370;
t652 = g(2) * t373 - t625;
t187 = Icges(5,3) * t343 + t344 * t438;
t431 = -t188 * t368 + t190 * t371;
t435 = -t118 * t368 + t121 * t371;
t642 = t258 * (-t186 * t373 - t435) - t259 * (-t186 * t370 + t690) + t312 * (t187 - t431);
t219 = (-Icges(5,2) * t371 - t597) * t343;
t641 = t258 * (-Icges(5,2) * t245 + t121 + t227) - t259 * (-Icges(5,2) * t243 - t120 - t225) + t312 * (t190 + t219);
t640 = t138 / 0.2e1;
t639 = t139 / 0.2e1;
t638 = t241 / 0.2e1;
t637 = -t258 / 0.2e1;
t636 = t258 / 0.2e1;
t635 = -t259 / 0.2e1;
t634 = t259 / 0.2e1;
t633 = t282 / 0.2e1;
t632 = t283 / 0.2e1;
t631 = -t312 / 0.2e1;
t630 = t312 / 0.2e1;
t629 = t370 / 0.2e1;
t628 = -t373 / 0.2e1;
t627 = -rSges(5,3) - pkin(6);
t624 = g(2) * t370;
t54 = Icges(5,5) * t103 + Icges(5,6) * t102 + Icges(5,3) * t403;
t56 = Icges(5,4) * t103 + Icges(5,2) * t102 + Icges(5,6) * t403;
t58 = Icges(5,1) * t103 + Icges(5,4) * t102 + Icges(5,5) * t403;
t8 = (-qJD(2) * t690 - t54) * t344 + (qJD(2) * t113 - t368 * t56 + t371 * t58 + (-t116 * t371 + t120 * t368) * qJD(4)) * t343;
t622 = t8 * t259;
t53 = Icges(5,5) * t101 + Icges(5,6) * t100 + Icges(5,3) * t402;
t55 = Icges(5,4) * t101 + Icges(5,2) * t100 + Icges(5,6) * t402;
t57 = Icges(5,1) * t101 + Icges(5,4) * t100 + Icges(5,5) * t402;
t9 = (qJD(2) * t435 - t53) * t344 + (qJD(2) * t115 - t368 * t55 + t371 * t57 + (-t118 * t371 - t121 * t368) * qJD(4)) * t343;
t621 = t9 * t258;
t216 = (-Icges(5,5) * t368 - Icges(5,6) * t371) * t343;
t97 = qJD(2) * t187 + qJD(4) * t216;
t189 = Icges(5,6) * t343 + t344 * t439;
t98 = qJD(2) * t189 + qJD(4) * t219;
t191 = Icges(5,5) * t343 + t344 * t442;
t222 = (-Icges(5,1) * t368 - t596) * t343;
t99 = qJD(2) * t191 + qJD(4) * t222;
t18 = (qJD(2) * t431 - t97) * t344 + (qJD(2) * t186 - t368 * t98 + t371 * t99 + (-t188 * t371 - t190 * t368) * qJD(4)) * t343;
t67 = -t186 * t344 + t343 * t431;
t617 = t18 * t312 + t67 * t241;
t616 = rSges(3,1) * t372;
t357 = t370 * rSges(3,3);
t606 = t46 * t139;
t47 = -t115 * t344 + t343 * t435;
t605 = t47 * t138;
t401 = qJD(2) * t661 + t347;
t70 = qJD(1) * t497 + t401;
t604 = t70 * t268;
t296 = rSges(3,1) * t369 + rSges(3,2) * t372;
t491 = t296 * t514;
t523 = rSges(3,2) * t568 + t373 * rSges(3,3);
t235 = rSges(3,1) * t565 - t523;
t535 = -t235 - t309;
t129 = qJD(1) * t535 - t491;
t584 = t129 * t370;
t583 = t129 * t373;
t236 = rSges(3,1) * t563 - rSges(3,2) * t567 + t357;
t182 = t236 + t310;
t130 = qJD(1) * t182 - t296 * t515;
t261 = t296 * t373;
t582 = t130 * t261;
t555 = -t123 + t238;
t554 = t124 + t240;
t549 = -t198 * t515 + t199 * t514;
t548 = -t370 * t198 + t373 * t199;
t207 = rSges(4,1) * t570 - rSges(4,2) * t572 + t356;
t543 = -t199 - t207;
t542 = -t199 - t240;
t502 = t343 * t613;
t531 = rSges(5,3) * t571 + t370 * t502;
t530 = rSges(5,3) * t570 + t373 * t502;
t529 = rSges(4,2) * t493 + rSges(4,3) * t517;
t525 = rSges(3,3) * t517 + t518 * t679;
t506 = pkin(2) * t567;
t503 = t343 * t615;
t500 = qJD(2) * t361;
t499 = t140 * t514 + t141 * t515 - t282 * t198;
t498 = t373 * t140 + t370 * t141 - t198 * t517;
t488 = t372 * t514;
t484 = -pkin(1) - t616;
t481 = t517 / 0.2e1;
t480 = -t515 / 0.2e1;
t479 = t515 / 0.2e1;
t478 = -t514 / 0.2e1;
t477 = t514 / 0.2e1;
t466 = (-t370 ^ 2 - t705) * t626;
t465 = -t193 + t474;
t299 = rSges(2,1) * t373 - rSges(2,2) * t370;
t297 = rSges(2,1) * t370 + rSges(2,2) * t373;
t298 = t616 - t679;
t454 = t36 * t373 - t37 * t370;
t453 = t36 * t370 + t37 * t373;
t452 = t370 * t39 - t373 * t38;
t451 = t370 * t38 + t373 * t39;
t450 = t370 * t47 - t373 * t46;
t449 = t370 * t46 + t373 * t47;
t111 = rSges(4,1) * t404 - rSges(4,2) * t490 + t529;
t437 = t111 * t373 + t112 * t370;
t434 = -t123 * t373 - t124 * t370;
t433 = -t130 * t370 - t583;
t150 = -rSges(3,2) * t488 + (-t372 * t518 - t489) * rSges(3,1) + t525;
t260 = t296 * t370;
t151 = -qJD(2) * t260 + (t298 * t373 + t357) * qJD(1);
t432 = t150 * t373 + t151 * t370;
t428 = t206 * t370 + t207 * t373;
t425 = t235 * t370 + t236 * t373;
t418 = -t104 - t250 - t500;
t417 = t343 * t627 - t332;
t400 = -t113 * t259 + t115 * t258 + t186 * t312;
t399 = (-Icges(5,5) * t242 - Icges(5,6) * t243) * t259 - (Icges(5,5) * t244 - Icges(5,6) * t245) * t258 - t216 * t312;
t395 = t343 * t399;
t388 = (Icges(5,1) * t244 - t118 - t598) * t258 - (-Icges(5,1) * t242 - t116 - t226) * t259 + (-t188 + t222) * t312;
t31 = -t123 * t258 + t124 * t259 + (t238 * t370 + t240 * t373) * qJD(2) + t549;
t386 = t31 * t434 + (t370 * t40 - t373 * t41) * t193;
t376 = t642 * t343;
t315 = pkin(6) * t570;
t313 = pkin(6) * t571;
t276 = t298 * qJD(2);
t239 = -pkin(3) * t572 + t315;
t237 = -pkin(3) * t573 + t313;
t230 = t268 * t373;
t170 = -t373 * t503 + t530;
t169 = -t370 * t503 + t531;
t168 = t190 * t373;
t167 = t190 * t370;
t166 = t188 * t373;
t165 = t188 * t370;
t159 = rSges(5,1) * t244 - rSges(5,2) * t245;
t158 = -rSges(5,1) * t242 - rSges(5,2) * t243;
t125 = t425 * qJD(2);
t71 = -t268 * t515 + (t207 + t541) * qJD(1) - t524;
t66 = qJD(2) * t428 + t549;
t62 = qJD(1) * t150 + qJDD(1) * t236 - t276 * t515 - t282 * t296 + t532;
t61 = -t276 * t514 + t283 * t296 + t535 * qJDD(1) + (-t151 - t279) * qJD(1);
t30 = qJD(1) * t111 + qJDD(1) * t207 - t249 * t515 - t268 * t282 + t379;
t16 = t102 * t188 + t103 * t190 + t186 * t403 - t242 * t98 + t243 * t99 + t573 * t97;
t15 = t100 * t188 + t101 * t190 + t186 * t402 + t244 * t98 + t245 * t99 + t572 * t97;
t14 = t258 * t47 - t259 * t46 + t312 * t67;
t7 = t102 * t118 + t103 * t121 + t115 * t403 - t242 * t55 + t243 * t57 + t53 * t573;
t6 = t102 * t116 - t103 * t120 + t113 * t403 - t242 * t56 + t243 * t58 + t54 * t573;
t5 = t100 * t118 + t101 * t121 + t115 * t402 + t244 * t55 + t245 * t57 + t53 * t572;
t4 = t100 * t116 - t101 * t120 + t113 * t402 + t244 * t56 + t245 * t58 + t54 * t572;
t3 = -t123 * t138 - t124 * t139 + t238 * t282 + t258 * t60 + t259 * t59 + t542 * t283 + (t142 * t373 + t143 * t370) * qJD(2) + t499;
t2 = t138 * t37 + t139 * t36 + t16 * t312 + t241 * t63 + t258 * t7 - t259 * t6;
t1 = t138 * t39 + t139 * t38 + t15 * t312 + t241 * t64 + t258 * t5 - t259 * t4;
t17 = [t13 * t634 + t15 * t636 + t63 * t639 + t64 * t640 + t605 / 0.2e1 + t606 / 0.2e1 - m(2) * (-g(1) * t297 + g(2) * t299) - t622 / 0.2e1 + t621 / 0.2e1 + t617 + (t16 + t13) * t635 + ((((t718 + t722) * t373 + t682 + t716 + t720) * t373 - t655 * t370) * qJD(2) + t694) * t477 + (-t710 * qJD(2) + t247 * t344 + t248 * t343 + t274 * t372 + t275 * t369) * qJD(1) + (t40 * (-t461 + t494) + t41 * (-pkin(3) * t485 + t288 + t419 + t504) + (t10 * t417 + t40 * (t344 * t627 + t623) * qJD(2)) * t370 + ((-t339 + t417) * t607 + (t40 * (-t339 - t657 - t329) - t41 * t367) * t373) * qJD(1) - t417 * t625 - (-qJD(1) * t238 - t40 + t660 + t706) * t41 + (t11 - g(2)) * (t468 + t554) + (t10 - g(1)) * (-t460 + t526)) * m(5) + (t70 * t494 + t71 * (t347 + t529) + (t604 * t370 + t661 * t71) * qJD(2) + ((-t70 * rSges(4,3) + t71 * (-t339 - t330)) * t370 + (t70 * (-t339 - t658) - t71 * t367) * t373) * qJD(1) - (-qJD(1) * t206 + t401 + t660 - t70) * t71 + (t30 - g(2)) * (t207 + t468) + t695 * (-t206 + t526)) * m(4) + (t130 * (t342 + t525) + (t296 * t584 - t582) * qJD(2) + ((-pkin(1) - t298) * t583 + (t129 * (-rSges(3,3) - pkin(5)) + t130 * t484) * t370) * qJD(1) - (-qJD(1) * t235 - t129 - t285 - t491) * t130 + (t62 - g(2)) * t182 + (t61 - g(1)) * (t484 * t370 + t362 + t523)) * m(3) + (t674 + t676) * t633 + (t675 + t677) * t632 + (((t373 * t656 + t655 + t680) * t373 + (t370 * t656 + t681 + t717) * t370) * qJD(2) + t689 - t692) * t480 + (t684 + t687) * t479 + (m(2) * (t297 ^ 2 + t299 ^ 2) + Icges(2,3) - t709) * qJDD(1) + (-t685 + t686 + t688) * t478; (g(1) * t261 + g(2) * t260 - g(3) * t298 - (t129 * t260 - t582) * qJD(1) - (t125 * (-t260 * t370 - t261 * t373) + t433 * t298) * qJD(2) + (qJD(2) * t432 + t235 * t282 - t236 * t283) * t425 + t125 * ((t235 * t373 - t236 * t370) * qJD(1) + t432) + t433 * t276 + (-t62 * t370 - t61 * t373 + (-t130 * t373 + t584) * qJD(1)) * t296) * m(3) - (t12 * t370 + t13 * t373) * t512 / 0.2e1 - ((t398 * t343 + t344 * t643 + t397 * t369 + t372 * t644) * qJD(2) + (t343 * t534 + t344 * t533 + t369 * t528 + t372 * t527) * qJD(1)) * qJD(1) / 0.2e1 + (qJD(1) * t449 + t370 * t9 - t373 * t8) * t630 + (qJD(1) * t453 + t370 * t7 - t373 * t6) * t635 + (qJD(1) * t451 + t370 * t5 - t373 * t4) * t636 + t450 * t638 + t452 * t640 + ((t166 * t242 - t168 * t243) * t258 - (t165 * t242 - t167 * t243) * t259 + (-t189 * t242 + t191 * t243) * t312 + (t343 * t63 + t37 * t570) * qJD(4) + ((qJD(4) * t36 + t400) * t344 + t376) * t370) * t634 + ((-t166 * t244 - t168 * t245) * t258 - (-t165 * t244 - t167 * t245) * t259 + (t189 * t244 + t191 * t245) * t312 + (t343 * t64 + t38 * t571) * qJD(4) + ((qJD(4) * t39 + t400) * t344 + t376) * t373) * t637 + (((t166 * t368 - t168 * t371 + t115) * t258 - (t165 * t368 - t167 * t371 + t113) * t259 + (-t189 * t368 + t191 * t371 + t186) * t312 + t67 * qJD(4)) * t343 + (qJD(4) * t449 - t642) * t344) * t631 + (-g(3) * t659 - t405 * t624 - (t70 * t229 + t71 * (-t230 - t506)) * qJD(1) - (t66 * t466 + (-t66 * t230 - t659 * t70) * t373 + (-t66 * t229 - t659 * t71) * t370) * qJD(2) + t70 * (-pkin(2) * t488 - t249 * t373) + (qJD(2) * t437 + t206 * t282 + t283 * t543 + t499) * (t428 + t548) + t66 * (t437 + t498) + (t66 * t206 + t405 * t71) * t517 + (t30 * t405 + t71 * (-t249 - t500) + (t543 * t66 + t604) * qJD(1)) * t370 + t695 * t661) * m(4) + (-g(1) * (t315 - t506 + t530) - g(2) * (t313 - t507 + t531) - g(3) * t194 - (g(1) * t373 + t624) * t343 * (-pkin(3) - t615) + t3 * t548 + (t3 * t554 + t40 * t418 + (qJD(1) * t41 + t10) * t465) * t373 + (t11 * t465 + t41 * t418 + t3 * t555 + t40 * (t193 + t271) * qJD(1)) * t370 - t40 * (-qJD(1) * t237 - t169 * t312 - t194 * t259) - t41 * (t170 * t312 - t194 * t258 + (t239 - t506) * qJD(1)) - ((t123 * t40 + t124 * t41) * t343 + t386 * t344) * qJD(4) + (-qJD(2) * t691 + g(3)) * (-t657 - t361) + (-t169 * t258 - t170 * t259 - (t237 * t370 + t239 * t373 + t466) * qJD(2) + t498 + (t555 * qJD(1) + t142 + t59) * t373 + (t143 + t60 + (-t124 + t542) * qJD(1)) * t370) * t31) * m(5) - t139 * t454 / 0.2e1 - t14 * t513 / 0.2e1 + t699 * t632 + t700 * t633 + (t12 + t689) * t518 / 0.2e1 + ((t370 * t681 + t373 * t680) * qJD(1) + t666) * t479 + ((t370 * t683 + t373 * t682) * qJD(1) + t667) * t478 + (t685 * t373 + t684 * t370 + (t677 * t370 + t676 * t373) * qJD(1)) * qJD(1) / 0.2e1 + (qJD(1) * t686 + qJD(2) * t667 + qJDD(1) * t675 + t282 * t682 + t283 * t683 + t2) * t628 + (qJD(1) * t687 + qJD(2) * t666 + qJDD(1) * t674 + t282 * t680 + t283 * t681 + t1) * t629 + (t13 + t688) * t481 + (t676 * t370 - t677 * t373) * qJDD(1) / 0.2e1 + ((-t514 * t663 - t665) * t373 + ((t373 * t664 + t669) * qJD(2) + t668) * t370) * t477 + ((-t515 * t664 + t665) * t370 + ((t663 * t370 + t669) * qJD(2) + t668) * t373) * t480; (0.2e1 * t29 * t629 + 0.2e1 * t30 * t628 + t652) * m(4) + (t652 + t653) * m(5); t1 * t572 / 0.2e1 + (t343 * t451 - t344 * t64) * t640 + ((qJD(2) * t451 - t15) * t344 + (-qJD(1) * t452 + qJD(2) * t64 + t370 * t4 + t373 * t5) * t343) * t636 + t2 * t573 / 0.2e1 + (t343 * t453 - t344 * t63) * t639 + ((qJD(2) * t453 - t16) * t344 + (qJD(1) * t454 + qJD(2) * t63 + t370 * t6 + t373 * t7) * t343) * t635 + t14 * t516 / 0.2e1 - t344 * (t605 + t606 + t617 + t621 - t622) / 0.2e1 + (t343 * t449 - t344 * t67) * t638 + ((qJD(2) * t449 - t18) * t344 + (-qJD(1) * t450 + qJD(2) * t67 + t370 * t8 + t373 * t9) * t343) * t630 + (t244 * t641 + t388 * t245 - t373 * t395) * t637 + (-t242 * t641 + t243 * t388 - t370 * t395) * t634 + (t399 * t344 + (-t368 * t641 + t371 * t388) * t343) * t631 + (-t493 / 0.2e1 + t344 * t477) * t13 + (t343 * t481 + t344 * t479) * t12 + ((qJD(2) * t386 - t10 * t123 - t11 * t124 + t40 * t60 - t41 * t59) * t344 + (t40 * (qJD(2) * t123 + t104 * t370) + t41 * (qJD(2) * t124 - t104 * t373) + t3 * t434 + t31 * (t123 * t518 - t124 * t517 - t370 * t59 + t373 * t60) + (qJD(1) * t691 + t653) * t193) * t343 - t40 * (-t158 * t312 - t228 * t259) - t41 * (t159 * t312 - t228 * t258) - t31 * (t158 * t258 + t159 * t259) - g(1) * t159 - g(2) * t158 - g(3) * t228) * m(5);];
tau = t17;

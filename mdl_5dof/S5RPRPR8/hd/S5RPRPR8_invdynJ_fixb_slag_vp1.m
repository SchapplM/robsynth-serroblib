% Calculate vector of inverse dynamics joint torques for
% S5RPRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR8_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR8_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:18
% EndTime: 2019-12-31 18:22:10
% DurationCPUTime: 45.38s
% Computational Cost: add. (25440->1074), mult. (27399->1405), div. (0->0), fcn. (25552->10), ass. (0->493)
t408 = qJ(1) + pkin(8);
t397 = sin(t408);
t399 = cos(t408);
t410 = cos(pkin(9));
t409 = sin(pkin(9));
t414 = cos(qJ(3));
t622 = t409 * t414;
t268 = t397 * t622 + t399 * t410;
t621 = t410 * t414;
t625 = t399 * t409;
t269 = t397 * t621 - t625;
t412 = sin(qJ(3));
t629 = t397 * t412;
t141 = Icges(5,5) * t269 - Icges(5,6) * t268 + Icges(5,3) * t629;
t628 = t397 * t414;
t643 = Icges(4,3) * t399;
t215 = Icges(4,5) * t628 - Icges(4,6) * t629 - t643;
t364 = Icges(4,4) * t629;
t650 = Icges(4,5) * t399;
t219 = Icges(4,1) * t628 - t364 - t650;
t646 = Icges(4,6) * t399;
t217 = Icges(4,4) * t628 - Icges(4,2) * t629 - t646;
t635 = t217 * t412;
t473 = -t219 * t414 + t635;
t144 = Icges(5,4) * t269 - Icges(5,2) * t268 + Icges(5,6) * t629;
t147 = Icges(5,1) * t269 - Icges(5,4) * t268 + Icges(5,5) * t629;
t477 = -t144 * t268 + t147 * t269;
t738 = t141 * t629 - t215 * t399 - t397 * t473 + t477;
t270 = t397 * t410 - t399 * t622;
t630 = t397 * t409;
t271 = t399 * t621 + t630;
t624 = t399 * t412;
t143 = Icges(5,5) * t271 + Icges(5,6) * t270 + Icges(5,3) * t624;
t146 = Icges(5,4) * t271 + Icges(5,2) * t270 + Icges(5,6) * t624;
t149 = Icges(5,1) * t271 + Icges(5,4) * t270 + Icges(5,5) * t624;
t46 = t143 * t629 - t268 * t146 + t269 * t149;
t402 = Icges(4,4) * t414;
t486 = -Icges(4,2) * t412 + t402;
t218 = Icges(4,6) * t397 + t399 * t486;
t654 = Icges(4,4) * t412;
t345 = Icges(4,1) * t414 - t654;
t220 = Icges(4,5) * t397 + t345 * t399;
t207 = t220 * t628;
t341 = Icges(4,5) * t414 - Icges(4,6) * t412;
t216 = Icges(4,3) * t397 + t341 * t399;
t524 = t216 * t399 - t207;
t82 = -t218 * t629 - t524;
t737 = t46 + t82;
t483 = Icges(5,5) * t410 - Icges(5,6) * t409;
t272 = -Icges(5,3) * t414 + t412 * t483;
t485 = Icges(5,4) * t410 - Icges(5,2) * t409;
t274 = -Icges(5,6) * t414 + t412 * t485;
t488 = Icges(5,1) * t410 - Icges(5,4) * t409;
t276 = -Icges(5,5) * t414 + t412 * t488;
t342 = Icges(4,2) * t414 + t654;
t344 = Icges(4,1) * t412 + t402;
t467 = t342 * t412 - t344 * t414;
t340 = Icges(4,5) * t412 + Icges(4,6) * t414;
t632 = t340 * t399;
t743 = -t268 * t274 + t269 * t276 + t272 * t629 - t397 * t467 - t632;
t633 = t340 * t397;
t742 = t270 * t274 + t271 * t276 + t272 * t624 - t399 * t467 + t633;
t753 = t737 * t397 - t738 * t399;
t415 = cos(qJ(1));
t406 = t415 * pkin(1);
t721 = t270 * t144 + t271 * t147;
t47 = t141 * t624 + t721;
t48 = t143 * t624 + t270 * t146 + t271 * t149;
t494 = t397 * t48 - t399 * t47;
t752 = t742 * qJD(1) + qJD(3) * t494;
t316 = rSges(3,1) * t397 + rSges(3,2) * t399;
t413 = sin(qJ(1));
t682 = pkin(1) * t413;
t300 = -t316 - t682;
t751 = t743 * qJD(1);
t407 = pkin(9) + qJ(5);
t396 = sin(t407);
t398 = cos(t407);
t627 = t398 * t399;
t245 = t396 * t628 + t627;
t626 = t398 * t414;
t246 = -t399 * t396 + t397 * t626;
t506 = t246 * rSges(6,1) - t245 * rSges(6,2);
t132 = -rSges(6,3) * t629 - t506;
t667 = rSges(6,2) * t396;
t670 = rSges(6,1) * t398;
t505 = -t667 + t670;
t256 = -rSges(6,3) * t414 + t412 * t505;
t571 = qJD(5) * t412;
t575 = qJD(3) * t399;
t303 = -t397 * t571 + t575;
t400 = qJD(4) * t412;
t350 = t399 * t400;
t570 = qJD(5) * t414;
t383 = qJD(1) - t570;
t411 = -pkin(7) - qJ(4);
t616 = qJ(4) + t411;
t715 = t414 * t616;
t393 = pkin(4) * t410 + pkin(3);
t674 = pkin(3) - t393;
t716 = t412 * t674;
t238 = t715 - t716;
t352 = pkin(3) * t412 - qJ(4) * t414;
t601 = -t238 - t352;
t750 = t383 * t132 - t303 * t256 + t575 * t601 + t350;
t122 = Icges(6,5) * t246 - Icges(6,6) * t245 + Icges(6,3) * t629;
t230 = Icges(6,4) * t246;
t125 = -Icges(6,2) * t245 + Icges(6,6) * t629 + t230;
t229 = Icges(6,4) * t245;
t129 = -Icges(6,1) * t246 - Icges(6,5) * t629 + t229;
t739 = t125 * t396 + t129 * t398;
t49 = -t122 * t414 - t739 * t412;
t416 = qJD(1) ^ 2;
t562 = t416 * t406;
t749 = t753 * qJD(3) + t751;
t623 = t399 * t414;
t610 = -t397 * t215 - t219 * t623;
t83 = -t217 * t624 - t610;
t609 = t397 * t216 + t220 * t623;
t84 = -t218 * t624 + t609;
t490 = t397 * t84 - t399 * t83;
t748 = qJD(3) * t490 + t752;
t447 = qJD(3) * t342;
t169 = qJD(1) * t218 - t397 * t447;
t448 = qJD(3) * t344;
t171 = qJD(1) * t220 - t397 * t448;
t720 = t144 * t409 - t147 * t410;
t574 = qJD(3) * t412;
t542 = t409 * t574;
t199 = qJD(1) * t270 + t397 * t542;
t541 = t410 * t574;
t200 = qJD(1) * t271 - t397 * t541;
t568 = t412 * qJD(1);
t573 = qJD(3) * t414;
t279 = t397 * t573 + t399 * t568;
t86 = Icges(5,5) * t200 + Icges(5,6) * t199 + Icges(5,3) * t279;
t88 = Icges(5,4) * t200 + Icges(5,2) * t199 + Icges(5,6) * t279;
t90 = Icges(5,1) * t200 + Icges(5,4) * t199 + Icges(5,5) * t279;
t747 = (t86 - t169) * t414 + (t409 * t88 - t410 * t90 - t171) * t412 + (-t141 * t412 + t414 * t720 + t473) * qJD(3);
t168 = -t399 * t447 + (-t397 * t486 + t646) * qJD(1);
t170 = -t399 * t448 + (-t345 * t397 + t650) * qJD(1);
t634 = t218 * t412;
t472 = -t220 * t414 + t634;
t475 = -t146 * t409 + t149 * t410;
t197 = qJD(1) * t268 + t399 * t542;
t198 = -qJD(1) * t269 - t399 * t541;
t543 = t399 * t573;
t548 = t397 * t568;
t280 = t543 - t548;
t85 = Icges(5,5) * t198 + Icges(5,6) * t197 + Icges(5,3) * t280;
t87 = Icges(5,4) * t198 + Icges(5,2) * t197 + Icges(5,6) * t280;
t89 = Icges(5,1) * t198 + Icges(5,4) * t197 + Icges(5,5) * t280;
t746 = (-t85 + t168) * t414 + (-t409 * t87 + t410 * t89 + t170) * t412 + (t143 * t412 + t414 * t475 - t472) * qJD(3);
t273 = Icges(5,3) * t412 + t414 * t483;
t242 = t273 * qJD(3);
t275 = Icges(5,6) * t412 + t414 * t485;
t243 = t275 * qJD(3);
t277 = Icges(5,5) * t412 + t414 * t488;
t244 = t277 * qJD(3);
t322 = t486 * qJD(3);
t323 = t345 * qJD(3);
t468 = t342 * t414 + t344 * t412;
t419 = qJD(1) * t340 - qJD(3) * t468 - t322 * t412 + t323 * t414;
t703 = t467 * qJD(1) + t341 * qJD(3);
t745 = t197 * t274 + t198 * t276 + t242 * t624 + t243 * t270 + t244 * t271 + t272 * t280 + t397 * t703 + t419 * t399;
t744 = t199 * t274 + t200 * t276 + t242 * t629 - t243 * t268 + t244 * t269 + t272 * t279 + t419 * t397 - t399 * t703;
t736 = t47 + t83;
t134 = t217 * t414 + t219 * t412;
t734 = -t141 * t414 - t412 * t720 + t134;
t135 = t218 * t414 + t220 * t412;
t733 = -t143 * t414 + t412 * t475 + t135;
t509 = t269 * rSges(5,1) - t268 * rSges(5,2);
t151 = -rSges(5,3) * t629 - t509;
t576 = qJD(3) * t397;
t302 = t399 * t571 + t576;
t41 = t122 * t629 - t125 * t245 - t129 * t246;
t631 = t397 * t398;
t247 = -t396 * t623 + t631;
t248 = t396 * t397 + t398 * t623;
t124 = Icges(6,5) * t248 + Icges(6,6) * t247 + Icges(6,3) * t624;
t653 = Icges(6,4) * t248;
t127 = Icges(6,2) * t247 + Icges(6,6) * t624 + t653;
t231 = Icges(6,4) * t247;
t130 = Icges(6,1) * t248 + Icges(6,5) * t624 + t231;
t42 = t124 * t629 - t245 * t127 + t246 * t130;
t482 = Icges(6,5) * t398 - Icges(6,6) * t396;
t249 = -Icges(6,3) * t414 + t412 * t482;
t651 = Icges(6,4) * t398;
t484 = -Icges(6,2) * t396 + t651;
t251 = -Icges(6,6) * t414 + t412 * t484;
t652 = Icges(6,4) * t396;
t487 = Icges(6,1) * t398 - t652;
t253 = -Icges(6,5) * t414 + t412 * t487;
t72 = -t245 * t251 + t246 * t253 + t249 * t629;
t12 = t302 * t42 - t303 * t41 + t383 * t72;
t43 = t122 * t624 + t247 * t125 - t129 * t248;
t44 = t124 * t624 + t247 * t127 + t248 * t130;
t73 = t247 * t251 + t248 * t253 + t249 * t624;
t13 = t302 * t44 - t303 * t43 + t73 * t383;
t735 = t48 + t84;
t387 = t397 * rSges(4,3);
t237 = rSges(4,1) * t623 - rSges(4,2) * t624 + t387;
t319 = t399 * pkin(2) + t397 * pkin(6);
t712 = t406 + t319;
t185 = t237 + t712;
t365 = t414 * t393;
t522 = -t411 * t412 + t365;
t404 = t412 * rSges(5,3);
t401 = t412 * qJ(4);
t405 = t414 * pkin(3);
t504 = t405 + t401;
t728 = -t404 - t504;
t581 = qJD(1) * t216;
t420 = -qJD(3) * t135 - t168 * t412 + t170 * t414 + t581;
t421 = qJD(1) * t215 - qJD(3) * t134 - t169 * t412 + t171 * t414;
t446 = qJD(3) * t340;
t704 = qJD(1) * t473 - t397 * t446 + t581;
t705 = -t399 * t446 + (-t341 * t397 + t472 + t643) * qJD(1);
t727 = (-t141 * t279 - t144 * t199 - t147 * t200 + t268 * t88 - t269 * t90 + t399 * t704 - t629 * t86) * t399 + (t143 * t279 + t146 * t199 + t149 * t200 - t268 * t87 + t269 * t89 + t420 * t397 + t629 * t85 + (-t421 - t705) * t399) * t397;
t726 = (-t141 * t280 - t144 * t197 - t147 * t198 - t270 * t88 - t271 * t90 - t421 * t399 - t624 * t86) * t399 + (t143 * t280 + t146 * t197 + t149 * t198 + t270 * t87 + t271 * t89 + t397 * t705 + t624 * t85 + (t420 - t704) * t399) * t397;
t566 = qJD(3) * qJD(4);
t725 = qJDD(4) * t412 + t414 * t566;
t724 = t12 * t397 + t13 * t399;
t469 = -t274 * t409 + t276 * t410;
t722 = qJD(3) * (-t397 * t475 - t399 * t720) + (t273 - t469) * qJD(1);
t719 = -rSges(5,3) - qJ(4);
t717 = t412 * t616;
t714 = t414 * t674;
t294 = t504 * t397;
t391 = t399 * pkin(6);
t318 = pkin(2) * t397 - t391;
t311 = qJD(1) * t318;
t713 = -qJD(1) * t294 - t311;
t317 = t399 * rSges(3,1) - rSges(3,2) * t397;
t301 = t317 + t406;
t710 = g(1) * t399 + g(2) * t397;
t707 = t710 * t412;
t577 = qJD(1) * t414;
t518 = -qJD(5) + t577;
t544 = t399 * t574;
t706 = t397 * t518 + t544;
t605 = -Icges(4,2) * t628 + t219 - t364;
t607 = t344 * t397 + t217;
t701 = -t412 * t605 - t414 * t607;
t250 = Icges(6,3) * t412 + t414 * t482;
t470 = -t251 * t396 + t253 * t398;
t479 = -t127 * t396 + t130 * t398;
t700 = t302 * (-t249 * t399 - t479) - t303 * (-t249 * t397 + t739) + t383 * (t250 - t470);
t286 = (-Icges(6,2) * t398 - t652) * t412;
t699 = t302 * (-Icges(6,2) * t248 + t130 + t231) - t303 * (-Icges(6,2) * t246 - t129 - t229) + t383 * (t253 + t286);
t698 = m(5) / 0.2e1;
t697 = m(6) / 0.2e1;
t696 = -m(5) - m(6);
t567 = qJD(1) * qJD(3);
t308 = qJDD(3) * t397 + t399 * t567;
t564 = qJDD(5) * t412;
t180 = qJD(5) * t280 + t399 * t564 + t308;
t695 = t180 / 0.2e1;
t309 = -qJDD(3) * t399 + t397 * t567;
t181 = qJD(5) * t279 + t397 * t564 + t309;
t694 = t181 / 0.2e1;
t693 = -t302 / 0.2e1;
t692 = t302 / 0.2e1;
t691 = -t303 / 0.2e1;
t690 = t303 / 0.2e1;
t689 = t308 / 0.2e1;
t688 = t309 / 0.2e1;
t320 = qJD(3) * t571 - qJDD(5) * t414 + qJDD(1);
t687 = t320 / 0.2e1;
t686 = -t383 / 0.2e1;
t685 = t383 / 0.2e1;
t681 = g(1) * t397;
t547 = t397 * t574;
t114 = t383 * t631 + (-t399 * t518 + t547) * t396;
t115 = t518 * t627 + (t383 * t396 - t398 * t574) * t397;
t60 = Icges(6,5) * t115 + Icges(6,6) * t114 + Icges(6,3) * t279;
t62 = Icges(6,4) * t115 + Icges(6,2) * t114 + Icges(6,6) * t279;
t64 = Icges(6,1) * t115 + Icges(6,4) * t114 + Icges(6,5) * t279;
t8 = (-qJD(3) * t739 - t60) * t414 + (qJD(3) * t122 - t396 * t62 + t398 * t64 + (-t125 * t398 + t129 * t396) * qJD(5)) * t412;
t678 = t8 * t303;
t464 = t383 * t399;
t112 = t396 * t706 + t398 * t464;
t113 = t396 * t464 - t398 * t706;
t59 = Icges(6,5) * t113 + Icges(6,6) * t112 + Icges(6,3) * t280;
t61 = Icges(6,4) * t113 + Icges(6,2) * t112 + Icges(6,6) * t280;
t63 = Icges(6,1) * t113 + Icges(6,4) * t112 + Icges(6,5) * t280;
t9 = (qJD(3) * t479 - t59) * t414 + (qJD(3) * t124 - t396 * t61 + t398 * t63 + (-t127 * t398 - t130 * t396) * qJD(5)) * t412;
t677 = t9 * t302;
t283 = (-Icges(6,5) * t396 - Icges(6,6) * t398) * t412;
t176 = qJD(3) * t250 + qJD(5) * t283;
t252 = Icges(6,6) * t412 + t414 * t484;
t177 = qJD(3) * t252 + qJD(5) * t286;
t254 = Icges(6,5) * t412 + t414 * t487;
t289 = (-Icges(6,1) * t396 - t651) * t412;
t178 = qJD(3) * t254 + qJD(5) * t289;
t34 = (qJD(3) * t470 - t176) * t414 + (qJD(3) * t249 - t177 * t396 + t178 * t398 + (-t251 * t398 - t253 * t396) * qJD(5)) * t412;
t93 = -t249 * t414 + t412 * t470;
t673 = t93 * t320 + t34 * t383;
t672 = rSges(4,1) * t414;
t671 = rSges(5,1) * t410;
t668 = rSges(5,2) * t409;
t563 = pkin(4) * t625;
t516 = -t397 * t365 + t563;
t154 = (t405 + t717) * t397 + t516;
t528 = -t318 - t682;
t514 = -t294 + t528;
t465 = t154 + t514;
t39 = qJD(1) * t465 + t750;
t661 = t39 * t399;
t403 = t412 * rSges(6,3);
t659 = t49 * t181;
t50 = -t124 * t414 + t412 * t479;
t658 = t50 * t180;
t656 = -rSges(6,3) + t411;
t582 = rSges(4,2) * t629 + t399 * rSges(4,3);
t236 = rSges(4,1) * t628 - t582;
t515 = -t236 + t528;
t353 = rSges(4,1) * t412 + rSges(4,2) * t414;
t545 = t353 * t575;
t117 = qJD(1) * t515 - t545;
t639 = t117 * t397;
t638 = t117 * t399;
t118 = qJD(1) * t185 - t353 * t576;
t297 = t353 * t399;
t637 = t118 * t297;
t619 = t412 * t393;
t618 = t414 * t411;
t613 = -t132 - t154;
t152 = t271 * rSges(5,1) + t270 * rSges(5,2) + rSges(5,3) * t624;
t375 = pkin(3) * t623;
t298 = t399 * t401 + t375;
t612 = -t152 - t298;
t358 = pkin(4) * t630;
t463 = t399 * t522 + t358;
t155 = t463 - t298;
t611 = -t155 - t298;
t572 = qJD(4) * t414;
t304 = qJD(3) * t504 - t572;
t439 = -t714 - t717;
t608 = -t439 * qJD(3) - t304;
t606 = -t344 * t399 - t218;
t604 = -t342 * t399 + t220;
t600 = t238 + t256;
t598 = t397 * t294 + t399 * t298;
t508 = -t668 + t671;
t596 = -(t414 * t508 + t404) * qJD(3) - t304;
t361 = qJ(4) * t623;
t296 = -pkin(3) * t624 + t361;
t594 = qJD(1) * t296 + t397 * t572;
t281 = -rSges(5,3) * t414 + t412 * t508;
t593 = -t281 - t352;
t592 = -rSges(5,1) * t621 + rSges(5,2) * t622 + t728;
t558 = t412 * t667;
t591 = rSges(6,3) * t628 + t397 * t558;
t590 = rSges(6,3) * t623 + t399 * t558;
t589 = qJD(1) * t563 + t411 * t548;
t559 = t412 * t668;
t588 = rSges(5,3) * t628 + t397 * t559;
t587 = rSges(5,3) * t623 + t399 * t559;
t578 = qJD(1) * t399;
t586 = rSges(4,2) * t548 + rSges(4,3) * t578;
t585 = -t342 + t345;
t584 = t344 + t486;
t380 = pkin(6) * t578;
t583 = t350 + t380;
t580 = qJD(1) * t341;
t579 = qJD(1) * t397;
t561 = t412 * t671;
t560 = t412 * t670;
t555 = t113 * rSges(6,1) + t112 * rSges(6,2) + rSges(6,3) * t543;
t333 = qJ(4) * t543;
t444 = -t397 * t577 - t544;
t162 = pkin(3) * t444 - qJ(4) * t548 + t333 + t350;
t338 = pkin(3) * t547;
t540 = t397 * t400;
t163 = qJ(4) * t279 + qJD(1) * t375 - t338 + t540;
t554 = t399 * t162 + t397 * t163 + t294 * t578;
t295 = (-rSges(6,1) * t396 - rSges(6,2) * t398) * t412;
t179 = qJD(5) * t295 + (t414 * t505 + t403) * qJD(3);
t553 = -t179 + t608;
t552 = t198 * rSges(5,1) + t197 * rSges(5,2) + rSges(5,3) * t543;
t359 = qJ(4) * t628;
t292 = -pkin(3) * t629 + t359;
t551 = t292 * t576 + t296 * t575 + t400;
t133 = t248 * rSges(6,1) + t247 * rSges(6,2) + rSges(6,3) * t624;
t550 = -t352 - t600;
t539 = -pkin(2) - t672;
t535 = t578 / 0.2e1;
t534 = -t576 / 0.2e1;
t533 = t576 / 0.2e1;
t532 = -t575 / 0.2e1;
t531 = t575 / 0.2e1;
t526 = t391 - t682;
t306 = t352 * t576;
t513 = t298 + t712;
t70 = -t306 + (-qJD(3) * t281 + t400) * t397 + (t152 + t513) * qJD(1);
t525 = t70 * t593;
t523 = -t215 + t634;
t521 = qJD(3) * t608;
t520 = qJD(3) * t596;
t519 = -qJD(1) * t292 + t399 * t572;
t517 = qJDD(1) * t406 - t416 * t682;
t512 = t656 * t412 - pkin(2);
t357 = rSges(2,1) * t415 - rSges(2,2) * t413;
t354 = rSges(2,1) * t413 + rSges(2,2) * t415;
t356 = -rSges(4,2) * t412 + t672;
t510 = rSges(5,1) * t200 + rSges(5,2) * t199;
t507 = -rSges(6,1) * t115 - rSges(6,2) * t114;
t499 = t397 * t42 - t399 * t41;
t498 = t397 * t41 + t399 * t42;
t497 = t397 * t44 - t399 * t43;
t496 = t397 * t43 + t399 * t44;
t493 = t397 * t50 - t399 * t49;
t492 = t397 * t49 + t399 * t50;
t481 = -t118 * t397 - t638;
t478 = -t132 * t399 - t133 * t397;
t174 = rSges(4,1) * t444 - rSges(4,2) * t543 + t586;
t293 = t353 * t397;
t175 = -qJD(3) * t293 + (t356 * t399 + t387) * qJD(1);
t474 = t174 * t399 + t175 * t397;
t471 = t236 * t397 + t237 * t399;
t466 = t151 + t514;
t257 = rSges(6,1) * t626 - t414 * t667 + t403;
t305 = t319 * qJD(1);
t462 = -t163 - t305 - t540;
t461 = t294 * t576 + t298 * t575 + qJD(2) - t572;
t452 = t309 * t352 + t399 * t725 - t562;
t451 = qJD(1) * (-pkin(2) * t579 + t380) + qJDD(1) * t319 + t517;
t450 = -t618 + t716;
t449 = (-t141 * t399 + t143 * t397) * t414;
t445 = t412 * t719 - pkin(2) - t405;
t442 = -t122 * t303 + t124 * t302 + t249 * t383;
t441 = (-Icges(6,5) * t245 - Icges(6,6) * t246) * t303 - (Icges(6,5) * t247 - Icges(6,6) * t248) * t302 - t283 * t383;
t440 = t450 * t399;
t438 = -pkin(2) + t728;
t436 = -qJDD(4) * t414 + t162 * t575 + t163 * t576 + t308 * t294 + t412 * t566 + qJDD(2);
t435 = -t412 * t604 + t414 * t606;
t433 = t412 * t441;
t426 = (-t412 * t584 + t414 * t585) * qJD(1);
t425 = qJDD(1) * t298 + t451 + t725 * t397 + (t162 + t350) * qJD(1);
t424 = (Icges(6,1) * t247 - t127 - t653) * t302 - (-Icges(6,1) * t245 - t125 - t230) * t303 + (-t251 + t289) * t383;
t33 = -t132 * t302 + t133 * t303 + (-t154 * t397 + t155 * t399) * qJD(3) + t461;
t40 = t133 * t383 - t256 * t302 - t306 + (-qJD(3) * t238 + t400) * t397 + (t155 + t513) * qJD(1);
t422 = t33 * t478 + (t39 * t397 - t399 * t40) * t256;
t418 = t700 * t412;
t417 = t272 * t577 + t412 * t722;
t325 = t356 * qJD(3);
t307 = t352 * t579;
t278 = (t397 ^ 2 + t399 ^ 2) * t574;
t211 = -t399 * t561 + t587;
t210 = -t397 * t561 + t588;
t206 = t276 * t399;
t205 = t276 * t397;
t204 = t274 * t399;
t203 = t274 * t397;
t193 = -t399 * t560 + t590;
t192 = -t397 * t560 + t591;
t191 = t253 * t399;
t190 = t253 * t397;
t189 = t251 * t399;
t188 = t251 * t397;
t183 = -t361 + t440;
t182 = t397 * t450 - t359;
t165 = rSges(6,1) * t247 - rSges(6,2) * t248;
t164 = -rSges(6,1) * t245 - rSges(6,2) * t246;
t116 = qJD(3) * t471 + qJD(2);
t92 = rSges(5,3) * t279 + t510;
t91 = -rSges(5,3) * t548 + t552;
t79 = t338 + (-t619 - t715) * t576 + (t399 * t439 + t358) * qJD(1);
t78 = -t333 + qJD(3) * t440 + (t401 + t714) * t579 + t589;
t69 = qJD(1) * t466 + t575 * t593 + t350;
t66 = rSges(6,3) * t279 - t507;
t65 = -rSges(6,3) * t548 + t555;
t58 = qJD(1) * t174 + qJDD(1) * t237 - t308 * t353 - t325 * t576 + t451;
t57 = -t562 - t325 * t575 + t309 * t353 + (-t175 - t305) * qJD(1) + t515 * qJDD(1);
t52 = (-t151 * t397 + t152 * t399) * qJD(3) + t461;
t51 = qJD(3) * t474 + t236 * t308 - t237 * t309 + qJDD(2);
t27 = qJD(1) * t91 + qJDD(1) * t152 + t308 * t593 + t397 * t520 + t425;
t26 = t281 * t309 + t399 * t520 + t466 * qJDD(1) + (t462 - t92) * qJD(1) + t452;
t25 = t114 * t251 + t115 * t253 + t176 * t629 - t177 * t245 + t178 * t246 + t249 * t279;
t24 = t112 * t251 + t113 * t253 + t176 * t624 + t177 * t247 + t178 * t248 + t249 * t280;
t21 = t302 * t50 - t303 * t49 + t383 * t93;
t20 = -t151 * t308 + t612 * t309 + (t397 * t92 + t399 * t91) * qJD(3) + t436;
t11 = qJD(1) * t78 + qJDD(1) * t155 + t133 * t320 - t179 * t302 - t180 * t256 + t308 * t601 + t383 * t65 + t397 * t521 + t425;
t10 = t320 * t132 - t303 * t179 + t181 * t256 + t309 * t238 - t383 * t66 + t399 * t521 + t465 * qJDD(1) + (t462 - t79) * qJD(1) + t452;
t7 = t114 * t127 + t115 * t130 + t124 * t279 - t245 * t61 + t246 * t63 + t59 * t629;
t6 = t114 * t125 - t115 * t129 + t122 * t279 - t245 * t62 + t246 * t64 + t60 * t629;
t5 = t112 * t127 + t113 * t130 + t124 * t280 + t247 * t61 + t248 * t63 + t59 * t624;
t4 = t112 * t125 - t113 * t129 + t122 * t280 + t247 * t62 + t248 * t64 + t60 * t624;
t3 = -t132 * t180 - t133 * t181 - t154 * t308 + t302 * t66 + t303 * t65 + t611 * t309 + (t397 * t79 + t399 * t78) * qJD(3) + t436;
t2 = t180 * t42 + t181 * t41 + t25 * t383 + t302 * t7 - t303 * t6 + t320 * t72;
t1 = t180 * t44 + t181 * t43 + t24 * t383 + t302 * t5 - t303 * t4 + t320 * t73;
t14 = [-t678 / 0.2e1 + t677 / 0.2e1 + t673 + t658 / 0.2e1 + t659 / 0.2e1 + t25 * t691 + t24 * t692 + t72 * t694 + t73 * t695 - m(2) * (-g(1) * t354 + g(2) * t357) + (((t82 - t207 + (t216 + t635) * t399 + t610) * t399 + t609 * t397) * qJD(3) + t752) * t531 + (t690 + t691) * t13 + ((t322 - t242) * t414 + (-t243 * t409 + t244 * t410 + t323) * t412 + (t272 * t412 + t414 * t469 - t467) * qJD(3)) * qJD(1) + ((-t316 * t416 - g(2) + t517) * t301 + (-t562 + (-0.2e1 * t317 - t406 + t301) * t416 - g(1)) * t300) * m(3) + (-(-t39 + (t154 - t682) * qJD(1) + t713 + t750) * t40 - t512 * t681 + t39 * t507 + t40 * (-t393 * t544 - t411 * t543 + t555 + t583 + t589) + (-t10 * pkin(2) + (-t39 * qJD(4) + t10 * t656) * t412 + t39 * (t414 * t656 + t619) * qJD(3)) * t397 + ((-t39 * t415 - t40 * t413) * pkin(1) + (t512 - t365) * t661 + (t39 * (-pkin(4) * t409 - pkin(6)) + t40 * (-pkin(2) - t365 - t403)) * t397) * qJD(1) + (-g(2) + t11) * (t463 + t712 + t133) + (t10 - g(1)) * (t516 + t526 - t506)) * m(6) + (-(t350 - t69 + (t151 - t682) * qJD(1) + t713) * t70 - t525 * t575 - t445 * t681 + t69 * (t338 - t510) + t70 * (-pkin(3) * t544 + t333 + t552 + t583) + (t26 * t438 + t69 * (t573 * t719 - t400)) * t397 + ((-t413 * t70 - t415 * t69) * pkin(1) + t69 * t445 * t399 + (-t69 * pkin(6) + t438 * t70) * t397) * qJD(1) + (t27 - g(2)) * (t712 - t612) + (t26 - g(1)) * (-t509 + t526)) * m(5) + (t118 * (t380 + t586) + (t353 * t639 - t637) * qJD(3) + ((-t117 * t415 - t118 * t413) * pkin(1) + (-pkin(2) - t356) * t638 + (t117 * (-rSges(4,3) - pkin(6)) + t118 * t539) * t397) * qJD(1) - (-t545 - t117 - t311 + (-t236 - t682) * qJD(1)) * t118 + (t58 - g(2)) * t185 + (t57 - g(1)) * (t539 * t397 + t526 + t582)) * m(4) + (t742 + t733) * t689 + (t743 + t734) * t688 + (t745 + t746) * t533 + (m(3) * (t300 ^ 2 + t317 * t301) + t468 - t272 * t414 + t412 * t469 + m(2) * (t354 ^ 2 + t357 ^ 2) + Icges(2,3) + Icges(3,3)) * qJDD(1) + (((t399 * t523 + t477 - t609 + t84) * t399 + (t397 * t523 - t46 + t524 - t721 + t736) * t397) * qJD(3) + t749 - t751) * t534 + (t744 - t747 + t748) * t532; m(3) * qJDD(2) + m(4) * t51 + m(5) * t20 + m(6) * t3 + (-m(3) - m(4) + t696) * g(3); ((t189 * t245 - t191 * t246) * t302 - (t188 * t245 - t190 * t246) * t303 + (-t245 * t252 + t246 * t254) * t383 + (t412 * t72 + t42 * t623) * qJD(5) + ((qJD(5) * t41 + t442) * t414 + t418) * t397) * t690 + ((-t189 * t247 - t191 * t248) * t302 - (-t188 * t247 - t190 * t248) * t303 + (t247 * t252 + t248 * t254) * t383 + (t412 * t73 + t43 * t628) * qJD(5) + ((qJD(5) * t44 + t442) * t414 + t418) * t399) * t693 - ((t412 * t585 + t414 * t584) * qJD(1) + ((t397 * t604 - t399 * t605) * t414 + (t397 * t606 + t399 * t607) * t412) * qJD(3) - t722 * t414 + ((-t275 * t409 + t277 * t410 + t272) * qJD(1) + ((t204 * t409 - t206 * t410 + t143) * t397 - (t203 * t409 - t205 * t410 + t141) * t399) * qJD(3)) * t412) * qJD(1) / 0.2e1 + ((-t576 * t632 + t580) * t397 + (t426 + (-t701 * t399 + (t633 + t435) * t397) * qJD(3)) * t399 + (-t204 * t270 - t206 * t271) * t576 + (t270 * t275 + t271 * t277) * qJD(1) + ((t270 * t203 + t271 * t205 + t449) * qJD(3) + t417) * t399) * t534 + ((-t575 * t633 - t580) * t399 + (t426 + (t435 * t397 + (t632 - t701) * t399) * qJD(3)) * t397 - (t203 * t268 - t205 * t269) * t575 + (-t268 * t275 + t269 * t277) * qJD(1) + ((t268 * t204 - t269 * t206 + t449) * qJD(3) + t417) * t397) * t531 + (((t189 * t396 - t191 * t398 + t124) * t302 - (t188 * t396 - t190 * t398 + t122) * t303 + (-t252 * t396 + t254 * t398 + t249) * t383 + t93 * qJD(5)) * t412 + (qJD(5) * t492 - t700) * t414) * t686 + (t397 * t733 - t399 * t734) * qJDD(1) / 0.2e1 - t724 * t570 / 0.2e1 + t753 * t688 + ((t736 * t397 + t399 * t735) * qJD(1) + t726) * t533 + ((t397 * t738 + t737 * t399) * qJD(1) + t727) * t532 - t21 * t571 / 0.2e1 + (-t39 * (-qJD(1) * t182 - t192 * t383 - t257 * t303 + t519) - t40 * (qJD(1) * t183 + t193 * t383 - t257 * t302 + t594) - t33 * (t192 * t302 + t193 * t303 + t551) - ((t132 * t39 + t133 * t40) * t412 + t422 * t414) * qJD(5) - ((t33 * t183 - t39 * t522) * t399 + (t33 * t182 - t40 * t522) * t397) * qJD(3) - g(1) * t590 - g(2) * t591 - g(3) * (t257 + t522) - t710 * (-t618 + (-t393 - t670) * t412) + t39 * t307 + t3 * t598 + t33 * t554 + (t10 * t550 + t39 * t553 + t3 * (t133 + t155) + t33 * (t65 + t78) + (t33 * t613 + t40 * t550) * qJD(1)) * t399 + (t11 * t550 + t40 * t553 + t3 * t613 + t33 * (t66 + t79) + (t39 * t600 + t33 * (-t133 + t611)) * qJD(1)) * t397) * m(6) + (-g(1) * (t361 + t587) - g(2) * (t359 + t588) + g(3) * t592 - (-pkin(3) - t671) * t707 - t69 * (-qJD(1) * t210 + t519) - t70 * (qJD(1) * t211 + t594) - t52 * t551 - ((t52 * t211 + t592 * t69) * t399 + (t52 * t210 + t592 * t70) * t397) * qJD(3) + t69 * t307 + t20 * t598 + t52 * t554 + (t26 * t593 + t69 * t596 + t20 * t152 + t52 * t91 + (-t151 * t52 + t525) * qJD(1)) * t399 + (t27 * t593 + t70 * t596 - t20 * t151 + t52 * t92 + (t69 * t281 + t52 * t612) * qJD(1)) * t397) * m(5) + (g(1) * t297 + g(2) * t293 - g(3) * t356 - (t117 * t293 - t637) * qJD(1) - (t116 * (-t293 * t397 - t297 * t399) + t481 * t356) * qJD(3) + t51 * t471 + t116 * ((t236 * t399 - t237 * t397) * qJD(1) + t474) + t481 * t325 + (-t58 * t397 - t57 * t399 + (-t118 * t399 + t639) * qJD(1)) * t353) * m(4) + (t494 + t490) * t689 + (t747 * t399 + t746 * t397 + (t734 * t397 + t733 * t399) * qJD(1)) * qJD(1) / 0.2e1 + (t13 + t748) * t535 + (t12 + t749) * t579 / 0.2e1 - (t744 * qJD(1) + t727 * qJD(3) + t743 * qJDD(1) + t737 * t308 + t738 * t309 + t2) * t399 / 0.2e1 + (t745 * qJD(1) + t726 * qJD(3) + t742 * qJDD(1) + t735 * t308 + t736 * t309 + t1) * t397 / 0.2e1 + (qJD(1) * t492 + t397 * t9 - t399 * t8) * t685 + t493 * t687 + (qJD(1) * t498 + t397 * t7 - t399 * t6) * t691 + (qJD(1) * t496 + t397 * t5 - t399 * t4) * t692 + t499 * t694 + t497 * t695; t696 * (-g(3) * t414 + t707) - m(5) * (t278 * t52 + t279 * t70 + t280 * t69) - m(6) * (t278 * t33 + t279 * t40 + t280 * t39) + 0.2e1 * ((t575 * t69 + t576 * t70 - t20) * t698 + (t39 * t575 + t40 * t576 - t3) * t697) * t414 + 0.2e1 * ((qJD(3) * t52 + t26 * t399 + t27 * t397 + t578 * t70 - t579 * t69) * t698 + (qJD(3) * t33 + t10 * t399 + t11 * t397 - t39 * t579 + t40 * t578) * t697) * t412; -t13 * t548 / 0.2e1 + t1 * t624 / 0.2e1 + (t412 * t496 - t414 * t73) * t695 + ((qJD(3) * t496 - t24) * t414 + (-qJD(1) * t497 + qJD(3) * t73 + t397 * t4 + t399 * t5) * t412) * t692 + t412 * t12 * t535 + t2 * t629 / 0.2e1 + (t412 * t498 - t414 * t72) * t694 + ((qJD(3) * t498 - t25) * t414 + (-qJD(1) * t499 + qJD(3) * t72 + t397 * t6 + t399 * t7) * t412) * t691 + t21 * t574 / 0.2e1 - t414 * (t658 + t659 + t673 + t677 - t678) / 0.2e1 + (t412 * t492 - t414 * t93) * t687 + ((qJD(3) * t492 - t34) * t414 + (-qJD(1) * t493 + qJD(3) * t93 + t397 * t8 + t399 * t9) * t412) * t685 + (t247 * t699 + t424 * t248 - t399 * t433) * t693 + (-t245 * t699 + t246 * t424 - t397 * t433) * t690 + (t441 * t414 + (-t396 * t699 + t398 * t424) * t412) * t686 + t724 * t573 / 0.2e1 + ((qJD(3) * t422 - t10 * t132 - t11 * t133 + t39 * t66 - t40 * t65) * t414 + (t39 * (qJD(3) * t132 + t179 * t397) + t40 * (qJD(3) * t133 - t179 * t399) + t3 * t478 + t33 * (t132 * t579 - t133 * t578 - t397 * t65 + t399 * t66) + (t10 * t397 - t11 * t399 + (t397 * t40 + t661) * qJD(1)) * t256) * t412 - t39 * (-t164 * t383 - t295 * t303) - t40 * (t165 * t383 - t295 * t302) - t33 * (t164 * t302 + t165 * t303) - g(1) * t165 - g(2) * t164 - g(3) * t295) * m(6);];
tau = t14;

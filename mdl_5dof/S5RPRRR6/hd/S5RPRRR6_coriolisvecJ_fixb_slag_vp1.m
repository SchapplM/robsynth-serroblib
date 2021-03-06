% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR6_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:57
% EndTime: 2019-12-31 19:01:35
% DurationCPUTime: 30.24s
% Computational Cost: add. (32017->939), mult. (29291->1277), div. (0->0), fcn. (26816->10), ass. (0->495)
t405 = qJ(1) + pkin(9);
t399 = cos(t405);
t413 = -pkin(7) - pkin(6);
t376 = t399 * t413;
t411 = cos(qJ(3));
t395 = pkin(3) * t411 + pkin(2);
t398 = sin(t405);
t597 = -t398 * t395 - t376;
t409 = sin(qJ(1));
t693 = pkin(1) * t409;
t747 = t597 - t693;
t406 = qJ(3) + qJ(4);
t401 = cos(t406);
t636 = t398 * t401;
t571 = rSges(5,1) * t636;
t746 = -t571 + t747;
t407 = sin(qJ(5));
t623 = t401 * t407;
t410 = cos(qJ(5));
t627 = t399 * t410;
t287 = t398 * t623 + t627;
t622 = t401 * t410;
t629 = t399 * t407;
t288 = t398 * t622 - t629;
t400 = sin(t406);
t637 = t398 * t400;
t151 = Icges(6,5) * t288 - Icges(6,6) * t287 + Icges(6,3) * t637;
t269 = Icges(6,4) * t288;
t154 = -Icges(6,2) * t287 + Icges(6,6) * t637 + t269;
t268 = Icges(6,4) * t287;
t158 = -Icges(6,1) * t288 - Icges(6,5) * t637 + t268;
t745 = t154 * t407 + t158 * t410;
t73 = -t151 * t401 - t400 * t745;
t379 = qJD(3) * t398;
t314 = qJD(4) * t398 + t379;
t579 = qJD(5) * t400;
t254 = t399 * t579 + t314;
t404 = qJD(3) + qJD(4);
t315 = t399 * t404;
t255 = -t398 * t579 + t315;
t578 = qJD(5) * t401;
t368 = qJD(1) - t578;
t63 = t151 * t637 - t154 * t287 - t158 * t288;
t633 = t398 * t410;
t289 = -t399 * t623 + t633;
t635 = t398 * t407;
t290 = t399 * t622 + t635;
t631 = t399 * t400;
t153 = Icges(6,5) * t290 + Icges(6,6) * t289 + Icges(6,3) * t631;
t666 = Icges(6,4) * t290;
t156 = Icges(6,2) * t289 + Icges(6,6) * t631 + t666;
t270 = Icges(6,4) * t289;
t159 = Icges(6,1) * t290 + Icges(6,5) * t631 + t270;
t64 = t153 * t637 - t287 * t156 + t288 * t159;
t488 = Icges(6,5) * t410 - Icges(6,6) * t407;
t258 = -Icges(6,3) * t401 + t400 * t488;
t664 = Icges(6,4) * t410;
t490 = -Icges(6,2) * t407 + t664;
t260 = -Icges(6,6) * t401 + t400 * t490;
t665 = Icges(6,4) * t407;
t494 = Icges(6,1) * t410 - t665;
t262 = -Icges(6,5) * t401 + t400 * t494;
t92 = t258 * t637 - t260 * t287 + t262 * t288;
t26 = t254 * t64 - t255 * t63 + t368 * t92;
t65 = t151 * t631 + t289 * t154 - t158 * t290;
t66 = t153 * t631 + t289 * t156 + t290 * t159;
t93 = t258 * t631 + t260 * t289 + t262 * t290;
t27 = t254 * t66 - t255 * t65 + t93 * t368;
t415 = qJD(1) ^ 2;
t507 = rSges(6,1) * t288 - rSges(6,2) * t287;
t160 = rSges(6,3) * t637 + t507;
t690 = pkin(4) * t401;
t330 = pkin(8) * t400 + t690;
t282 = t330 * t398;
t741 = t160 + t282;
t162 = t290 * rSges(6,1) + t289 * rSges(6,2) + rSges(6,3) * t631;
t630 = t399 * t401;
t284 = pkin(4) * t630 + pkin(8) * t631;
t740 = t162 + t284;
t657 = Icges(5,6) * t399;
t237 = Icges(5,4) * t636 - Icges(5,2) * t637 - t657;
t394 = Icges(5,4) * t401;
t325 = Icges(5,1) * t400 + t394;
t739 = -t325 * t398 - t237;
t491 = -Icges(5,2) * t400 + t394;
t238 = Icges(5,6) * t398 + t399 * t491;
t738 = -t325 * t399 - t238;
t667 = Icges(5,4) * t400;
t326 = Icges(5,1) * t401 - t667;
t240 = Icges(5,5) * t398 + t326 * t399;
t323 = Icges(5,2) * t401 + t667;
t737 = -t323 * t399 + t240;
t506 = rSges(6,1) * t410 - rSges(6,2) * t407;
t265 = -rSges(6,3) * t401 + t400 * t506;
t691 = pkin(4) * t400;
t329 = -pkin(8) * t401 + t691;
t736 = t265 + t329;
t735 = t325 + t491;
t341 = Icges(5,4) * t637;
t662 = Icges(5,5) * t399;
t239 = Icges(5,1) * t636 - t341 - t662;
t647 = t237 * t400;
t482 = -t239 * t401 + t647;
t461 = t482 * t398;
t322 = Icges(5,5) * t401 - Icges(5,6) * t400;
t236 = Icges(5,3) * t398 + t322 * t399;
t614 = t398 * t236 + t240 * t630;
t734 = -t461 - t614;
t733 = t26 * t398 + t27 * t399;
t408 = sin(qJ(3));
t581 = qJD(3) * t408;
t555 = t399 * t581;
t523 = pkin(3) * t555;
t732 = -t160 * t368 - t255 * t265 - t315 * t329 - t523;
t731 = 0.2e1 * qJD(3);
t392 = t399 * pkin(6);
t319 = pkin(2) * t398 - t392;
t230 = t319 + t597;
t313 = qJD(1) * t319;
t729 = qJD(1) * t230 - t313;
t387 = t398 * rSges(4,3);
t626 = t399 * t411;
t628 = t399 * t408;
t257 = rSges(4,1) * t626 - rSges(4,2) * t628 + t387;
t391 = t398 * pkin(6);
t320 = t399 * pkin(2) + t391;
t412 = cos(qJ(1));
t403 = t412 * pkin(1);
t539 = t320 + t403;
t728 = t257 + t539;
t727 = -rSges(5,2) * t637 - t399 * rSges(5,3);
t536 = t399 * rSges(3,1) - rSges(3,2) * t398;
t402 = Icges(4,4) * t411;
t492 = -Icges(4,2) * t408 + t402;
t356 = Icges(4,1) * t408 + t402;
t726 = t403 + t536;
t214 = t265 * t398;
t466 = t506 * t401;
t683 = rSges(6,3) * t400;
t266 = t466 + t683;
t281 = t329 * t398;
t554 = t398 * t578;
t587 = qJD(1) * t398;
t725 = -qJD(1) * t281 + t160 * t579 - t214 * t368 + t255 * t266 - t265 * t554 + t315 * t330 + t587 * t736;
t327 = rSges(5,1) * t400 + rSges(5,2) * t401;
t279 = t327 * t398;
t280 = t327 * t399;
t684 = rSges(5,2) * t400;
t328 = rSges(5,1) * t401 - t684;
t241 = t571 + t727;
t386 = t398 * rSges(5,3);
t596 = rSges(5,1) * t630 + t386;
t242 = -rSges(5,2) * t631 + t596;
t343 = t399 * t395;
t525 = -t398 * t413 + t343;
t231 = t525 - t320;
t582 = qJD(3) * t399;
t549 = -t230 * t379 + t231 * t582 + qJD(2);
t78 = t241 * t314 + t242 * t315 + t549;
t455 = -t315 * t327 - t523;
t540 = -t319 - t693;
t518 = t230 + t540;
t94 = (-t241 + t518) * qJD(1) + t455;
t570 = pkin(3) * t581;
t350 = t398 * t570;
t517 = t231 + t539;
t95 = -t314 * t327 - t350 + (t242 + t517) * qJD(1);
t722 = -t94 * (qJD(1) * t279 - t315 * t328) - t78 * (-t314 * t279 - t280 * t315) - t95 * (-qJD(1) * t280 - t314 * t328);
t632 = t398 * t411;
t634 = t398 * t408;
t655 = Icges(4,3) * t399;
t248 = Icges(4,5) * t632 - Icges(4,6) * t634 - t655;
t365 = Icges(4,4) * t634;
t663 = Icges(4,5) * t399;
t252 = Icges(4,1) * t632 - t365 - t663;
t658 = Icges(4,6) * t399;
t250 = Icges(4,4) * t632 - Icges(4,2) * t634 - t658;
t644 = t250 * t408;
t480 = -t252 * t411 + t644;
t103 = -t248 * t399 - t398 * t480;
t321 = Icges(5,5) * t400 + Icges(5,6) * t401;
t641 = t321 * t399;
t646 = t238 * t400;
t654 = Icges(5,3) * t399;
t721 = -t404 * t641 + (-t240 * t401 - t322 * t398 + t646 + t654) * qJD(1);
t590 = qJD(1) * t236;
t642 = t321 * t398;
t720 = qJD(1) * t482 - t404 * t642 + t590;
t353 = Icges(4,5) * t411 - Icges(4,6) * t408;
t352 = Icges(4,5) * t408 + Icges(4,6) * t411;
t458 = qJD(3) * t352;
t668 = Icges(4,4) * t408;
t357 = Icges(4,1) * t411 - t668;
t253 = Icges(4,5) * t398 + t357 * t399;
t251 = Icges(4,6) * t398 + t399 * t492;
t643 = t251 * t408;
t479 = -t253 * t411 + t643;
t719 = -t399 * t458 + (-t353 * t398 + t479 + t655) * qJD(1);
t249 = Icges(4,3) * t398 + t353 * t399;
t589 = qJD(1) * t249;
t718 = qJD(1) * t480 - t398 * t458 + t589;
t476 = t323 * t400 - t325 * t401;
t717 = qJD(1) * t476 + t322 * t404;
t354 = Icges(4,2) * t411 + t668;
t475 = t354 * t408 - t356 * t411;
t716 = t475 * qJD(1) + t353 * qJD(3);
t603 = -Icges(4,2) * t632 + t252 - t365;
t605 = t356 * t398 + t250;
t715 = -t408 * t603 - t411 * t605;
t505 = -rSges(6,1) * t407 - rSges(6,2) * t410;
t168 = t404 * t466 + (rSges(6,3) * t404 + qJD(5) * t505) * t400;
t624 = t401 * t404;
t565 = t399 * t624;
t317 = pkin(8) * t565;
t625 = t400 * t404;
t566 = t399 * t625;
t584 = qJD(1) * t401;
t451 = -t398 * t584 - t566;
t585 = qJD(1) * t400;
t561 = t398 * t585;
t171 = pkin(4) * t451 - pkin(8) * t561 + t317;
t575 = qJD(1) * qJD(3);
t372 = t399 * t575;
t574 = qJD(1) * qJD(4);
t311 = t399 * t574 + t372;
t450 = -t561 + t565;
t198 = qJD(5) * t450 + t311;
t295 = t330 * t404;
t586 = qJD(1) * t399;
t377 = pkin(6) * t586;
t687 = pkin(2) - t395;
t186 = -t523 - t377 + (t398 * t687 - t376) * qJD(1);
t573 = t415 * t693;
t519 = qJD(1) * (-pkin(2) * t587 + t377) - t573;
t621 = t411 * qJD(3) ^ 2;
t420 = qJD(1) * t186 + (-t372 * t408 - t398 * t621) * pkin(3) + t519;
t552 = t404 * t579;
t441 = t368 * t410 + t407 * t625;
t522 = -qJD(5) + t584;
t133 = t399 * t441 + t522 * t635;
t440 = t368 * t407 - t410 * t625;
t134 = t399 * t440 - t522 * t633;
t563 = t134 * rSges(6,1) + t133 * rSges(6,2) + rSges(6,3) * t565;
t87 = -rSges(6,3) * t561 + t563;
t28 = qJD(1) * t171 + t162 * t552 - t168 * t254 - t198 * t265 - t295 * t314 - t311 * t329 + t368 * t87 + t420;
t567 = t398 * t624;
t452 = t399 * t585 + t567;
t568 = t398 * t625;
t172 = t452 * pkin(8) + (t399 * t584 - t568) * pkin(4);
t371 = t398 * t575;
t310 = t398 * t574 + t371;
t197 = qJD(5) * t452 + t310;
t572 = t415 * t403;
t444 = -pkin(3) * t399 * t621 + qJD(1) * t350 - t572;
t594 = t413 * t587 + t350;
t187 = (-t399 * t687 - t391) * qJD(1) - t594;
t312 = t320 * qJD(1);
t616 = -t187 - t312;
t135 = t398 * t441 - t522 * t629;
t136 = t398 * t440 + t522 * t627;
t508 = rSges(6,1) * t136 + rSges(6,2) * t135;
t88 = rSges(6,3) * t452 + t508;
t29 = -t160 * t552 - t168 * t255 + t197 * t265 - t295 * t315 + t310 * t329 - t368 * t88 + (-t172 + t616) * qJD(1) + t444;
t62 = t162 * t368 - t254 * t265 - t314 * t329 - t350 + (t284 + t517) * qJD(1);
t714 = (qJD(1) * t62 + t29) * t399 + t28 * t398;
t462 = t488 * t401;
t477 = -t260 * t407 + t262 * t410;
t484 = -t156 * t407 + t159 * t410;
t713 = t254 * (-t258 * t399 - t484) - t255 * (-t258 * t398 + t745) + t368 * (Icges(6,3) * t400 + t462 - t477);
t489 = -Icges(6,2) * t410 - t665;
t712 = t254 * (-Icges(6,2) * t290 + t159 + t270) - t255 * (-Icges(6,2) * t288 - t158 - t268) + t368 * (t489 * t400 + t262);
t711 = qJD(1) * t735 + t314 * t737 - t315 * (-Icges(5,2) * t636 + t239 - t341);
t710 = t197 / 0.2e1;
t709 = t198 / 0.2e1;
t708 = -t254 / 0.2e1;
t707 = t254 / 0.2e1;
t706 = -t255 / 0.2e1;
t705 = t255 / 0.2e1;
t704 = t310 / 0.2e1;
t703 = t311 / 0.2e1;
t702 = -t314 / 0.2e1;
t701 = t314 / 0.2e1;
t700 = -t315 / 0.2e1;
t699 = t315 / 0.2e1;
t698 = -t368 / 0.2e1;
t697 = t368 / 0.2e1;
t696 = t398 / 0.2e1;
t695 = -t399 / 0.2e1;
t694 = -rSges(6,3) - pkin(8);
t692 = pkin(3) * t408;
t689 = -qJD(1) / 0.2e1;
t688 = qJD(1) / 0.2e1;
t686 = rSges(4,1) * t411;
t82 = Icges(6,5) * t136 + Icges(6,6) * t135 + Icges(6,3) * t452;
t84 = Icges(6,4) * t136 + Icges(6,2) * t135 + Icges(6,6) * t452;
t86 = Icges(6,1) * t136 + Icges(6,4) * t135 + Icges(6,5) * t452;
t20 = (-t404 * t745 - t82) * t401 + (t151 * t404 - t407 * t84 + t410 * t86 + (-t154 * t410 + t158 * t407) * qJD(5)) * t400;
t681 = t20 * t255;
t81 = Icges(6,5) * t134 + Icges(6,6) * t133 + Icges(6,3) * t450;
t83 = Icges(6,4) * t134 + Icges(6,2) * t133 + Icges(6,6) * t450;
t85 = Icges(6,1) * t134 + Icges(6,4) * t133 + Icges(6,5) * t450;
t21 = (t404 * t484 - t81) * t401 + (t153 * t404 - t407 * t83 + t410 * t85 + (-t156 * t410 - t159 * t407) * qJD(5)) * t400;
t680 = t21 * t254;
t674 = t398 * t62;
t673 = t73 * t197;
t74 = -t153 * t401 + t400 * t484;
t672 = t74 * t198;
t107 = -t258 * t401 + t400 * t477;
t487 = -Icges(6,5) * t407 - Icges(6,6) * t410;
t163 = t404 * t462 + (Icges(6,3) * t404 + qJD(5) * t487) * t400;
t463 = t490 * t401;
t164 = t404 * t463 + (Icges(6,6) * t404 + qJD(5) * t489) * t400;
t464 = t494 * t401;
t493 = -Icges(6,1) * t407 - t664;
t165 = t404 * t464 + (Icges(6,5) * t404 + qJD(5) * t493) * t400;
t47 = (t404 * t477 - t163) * t401 + (-t164 * t407 + t165 * t410 + t258 * t404 + (-t260 * t410 - t262 * t407) * qJD(5)) * t400;
t671 = t107 * t552 + t47 * t368;
t56 = t160 * t254 + t162 * t255 + t282 * t314 + t284 * t315 + t549;
t652 = qJD(1) * t56;
t591 = rSges(4,2) * t634 + t399 * rSges(4,3);
t256 = rSges(4,1) * t632 - t591;
t358 = rSges(4,1) * t408 + rSges(4,2) * t411;
t556 = t358 * t582;
t126 = -t556 + (-t256 + t540) * qJD(1);
t650 = t126 * t398;
t649 = t126 * t399;
t558 = t358 * t379;
t127 = qJD(1) * t728 - t558;
t304 = t358 * t399;
t648 = t127 * t304;
t640 = t323 * t404;
t639 = t352 * t398;
t638 = t352 * t399;
t617 = -t168 - t295;
t235 = Icges(5,5) * t636 - Icges(5,6) * t637 - t654;
t615 = -t398 * t235 - t239 * t630;
t612 = -t398 * t230 + t399 * t231;
t611 = t398 * t241 + t399 * t242;
t610 = -t398 * t248 - t252 * t626;
t609 = t398 * t249 + t253 * t626;
t604 = -t356 * t399 - t251;
t602 = -t354 * t399 + t253;
t598 = rSges(5,2) * t561 + rSges(5,3) * t586;
t583 = qJD(1) * t408;
t595 = rSges(4,2) * t398 * t583 + rSges(4,3) * t586;
t593 = -t354 + t357;
t592 = t356 + t492;
t588 = qJD(1) * t353;
t580 = qJD(3) * t411;
t124 = -t398 * t476 - t641;
t577 = t124 * qJD(1);
t180 = -t398 * t475 - t638;
t576 = t180 * qJD(1);
t569 = pkin(3) * t580;
t145 = rSges(5,1) * t451 - rSges(5,2) * t565 + t598;
t146 = -t404 * t279 + (t328 * t399 + t386) * qJD(1);
t564 = t399 * t145 + t398 * t146 + t241 * t586;
t562 = t399 * t186 + t398 * t187 - t230 * t586;
t559 = t399 * t583;
t553 = t399 * t578;
t551 = t625 / 0.2e1;
t548 = -pkin(2) - t686;
t547 = t587 / 0.2e1;
t546 = t586 / 0.2e1;
t545 = -t379 / 0.2e1;
t542 = t582 / 0.2e1;
t538 = -t327 - t692;
t535 = (-t398 ^ 2 - t399 ^ 2) * t408;
t534 = (-t398 * t491 + t657) * qJD(1) + t737 * t404;
t533 = qJD(1) * t238 + t239 * t404 - t398 * t640;
t532 = (-t326 * t398 + t662) * qJD(1) + t738 * t404;
t531 = qJD(1) * t240 + t404 * t739;
t227 = t253 * t632;
t530 = t249 * t399 - t227;
t529 = -t235 + t646;
t528 = -t248 + t643;
t527 = t735 * t404;
t526 = t326 * t404 - t640;
t521 = t398 * t741 + t399 * t740;
t520 = qJD(5) * t551;
t294 = t328 * t404;
t513 = -t294 - t569;
t200 = t240 * t636;
t512 = t238 * t637 - t200;
t318 = rSges(3,1) * t398 + rSges(3,2) * t399;
t509 = -rSges(4,2) * t408 + t686;
t61 = (-t282 + t518) * qJD(1) + t732;
t504 = t399 * t61 + t674;
t503 = t398 * t64 - t399 * t63;
t502 = t398 * t63 + t399 * t64;
t501 = t398 * t66 - t399 * t65;
t500 = t398 * t65 + t399 * t66;
t499 = t398 * t74 - t399 * t73;
t498 = t398 * t73 + t399 * t74;
t497 = -t398 * t95 - t399 * t94;
t486 = -t127 * t398 - t649;
t483 = t160 * t399 - t162 * t398;
t114 = t237 * t401 + t239 * t400;
t147 = t250 * t411 + t252 * t408;
t148 = t251 * t411 + t253 * t408;
t478 = t256 * t398 + t257 * t399;
t474 = -t569 + t617;
t473 = t741 * t586 + (t171 + t87) * t399 + (t172 + t88) * t398;
t303 = t358 * t398;
t460 = qJD(3) * t356;
t459 = qJD(3) * t354;
t104 = -t251 * t634 - t530;
t457 = (-t103 * t399 + t104 * t398) * qJD(3);
t105 = -t250 * t628 - t610;
t106 = -t251 * t628 + t609;
t456 = (-t105 * t399 + t106 * t398) * qJD(3);
t454 = -t330 - t683;
t449 = -t151 * t255 + t153 * t254 + t258 * t368;
t448 = (-Icges(6,5) * t287 - Icges(6,6) * t288) * t255 - (Icges(6,5) * t289 - Icges(6,6) * t290) * t254 - t487 * t400 * t368;
t447 = qJD(1) * t322 - t314 * t641 + t315 * t642;
t446 = t186 * t582 + t187 * t379 - t230 * t372 - t231 * t371;
t445 = -t408 * t602 + t411 * t604;
t443 = t400 * t448;
t435 = (-t408 * t592 + t411 * t593) * qJD(1);
t432 = (Icges(6,1) * t289 - t156 - t666) * t254 - (-Icges(6,1) * t287 - t154 - t269) * t255 + (t493 * t400 - t260) * t368;
t183 = -rSges(4,2) * t399 * t580 + (-t411 * t587 - t555) * rSges(4,1) + t595;
t184 = -qJD(3) * t303 + (t399 * t509 + t387) * qJD(1);
t430 = t183 * t399 + t184 * t398 + (t256 * t399 - t257 * t398) * qJD(1);
t429 = t738 * t314 - t739 * t315 + (-t323 + t326) * qJD(1);
t426 = qJD(1) * t235 - t400 * t533 + t401 * t531;
t425 = -t400 * t534 + t401 * t532 + t590;
t424 = qJD(1) * t321 - t400 * t527 + t401 * t526;
t177 = qJD(1) * t251 - t398 * t459;
t179 = qJD(1) * t253 - t398 * t460;
t423 = qJD(1) * t248 - qJD(3) * t147 - t177 * t408 + t179 * t411;
t176 = -t399 * t459 + (-t398 * t492 + t658) * qJD(1);
t178 = -t399 * t460 + (-t357 * t398 + t663) * qJD(1);
t422 = -qJD(3) * t148 - t176 * t408 + t178 * t411 + t589;
t332 = t492 * qJD(3);
t333 = t357 * qJD(3);
t421 = qJD(1) * t352 - t332 * t408 + t333 * t411 + (-t354 * t411 - t356 * t408) * qJD(3);
t115 = t238 * t401 + t240 * t400;
t16 = t133 * t154 - t134 * t158 + t151 * t450 + t289 * t84 + t290 * t86 + t631 * t82;
t17 = t133 * t156 + t134 * t159 + t153 * t450 + t289 * t83 + t290 * t85 + t631 * t81;
t18 = t135 * t154 - t136 * t158 + t151 * t452 - t287 * t84 + t288 * t86 + t637 * t82;
t19 = t135 * t156 + t136 * t159 + t153 * t452 - t287 * t83 + t288 * t85 + t637 * t81;
t205 = t260 * t398;
t206 = t260 * t399;
t207 = t262 * t398;
t208 = t262 * t399;
t261 = Icges(6,6) * t400 + t463;
t263 = Icges(6,5) * t400 + t464;
t37 = t133 * t260 + t134 * t262 + t163 * t631 + t164 * t289 + t165 * t290 + t258 * t450;
t3 = -t16 * t255 + t17 * t254 + t197 * t65 + t198 * t66 + t368 * t37 + t552 * t93;
t32 = t107 * t368 + t254 * t74 - t255 * t73;
t38 = t135 * t260 + t136 * t262 + t163 * t637 - t164 * t287 + t165 * t288 + t258 * t452;
t4 = -t18 * t255 + t19 * t254 + t197 * t63 + t198 * t64 + t368 * t38 + t552 * t92;
t40 = t398 * t720 + t426 * t399;
t41 = t398 * t721 + t425 * t399;
t417 = -t400 * t711 + t429 * t401;
t418 = t713 * t400;
t42 = t426 * t398 - t399 * t720;
t43 = t425 * t398 - t399 * t721;
t96 = -t235 * t399 - t461;
t97 = -t236 * t399 - t512;
t54 = t314 * t97 - t315 * t96 + t577;
t125 = -t399 * t476 + t642;
t117 = t125 * qJD(1);
t98 = -t237 * t631 - t615;
t99 = -t238 * t631 + t614;
t55 = t314 * t99 - t315 * t98 + t117;
t69 = t400 * t531 + t401 * t533;
t70 = t400 * t532 + t401 * t534;
t71 = t398 * t717 + t424 * t399;
t72 = t424 * t398 - t399 * t717;
t419 = (t398 * t70 - t399 * t69 + (t114 * t398 + t115 * t399) * qJD(1)) * t688 + (t429 * t400 + t401 * t711) * t689 + (((t206 * t407 - t208 * t410 + t153) * t254 - (t205 * t407 - t207 * t410 + t151) * t255 + (-t261 * t407 + t263 * t410 + t258) * t368 + t107 * qJD(5)) * t400 + (qJD(5) * t498 - t713) * t401) * t698 + t499 * t520 - t32 * t579 / 0.2e1 + (qJD(1) * t498 - t20 * t399 + t21 * t398) * t697 + (t398 * t417 - t399 * t447) * t699 + (t398 * t43 - t399 * t42 + (t398 * t96 + t399 * t97) * qJD(1)) * t700 + (t398 * t41 - t399 * t40 + (t398 * t98 + t399 * t99) * qJD(1)) * t701 + (t398 * t447 + t399 * t417) * t702 + (t398 * t99 - t399 * t98) * t703 + (t398 * t97 - t399 * t96) * t704 + ((t206 * t287 - t208 * t288) * t254 - (t205 * t287 - t207 * t288) * t255 + (-t261 * t287 + t263 * t288) * t368 + (t400 * t92 + t630 * t64) * qJD(5) + ((qJD(5) * t63 + t449) * t401 + t418) * t398) * t705 + (qJD(1) * t502 - t18 * t399 + t19 * t398) * t706 + (qJD(1) * t500 - t16 * t399 + t17 * t398) * t707 + ((-t206 * t289 - t208 * t290) * t254 - (-t205 * t289 - t207 * t290) * t255 + (t261 * t289 + t263 * t290) * t368 + (t400 * t93 + t636 * t65) * qJD(5) + ((qJD(5) * t66 + t449) * t401 + t418) * t399) * t708 + t501 * t709 + t503 * t710 + (qJD(1) * t71 + t310 * t98 + t311 * t99 + t314 * t41 - t315 * t40 + t3) * t696 + (qJD(1) * t72 + t310 * t96 + t311 * t97 + t314 * t43 - t315 * t42 + t4) * t695 + (t54 + t26) * t547 + (t55 + t27) * t546 - t733 * t578 / 0.2e1;
t215 = t265 * t399;
t283 = t329 * t399;
t416 = t56 * (t160 * t553 - t162 * t554 - t254 * t214 - t215 * t255 - t314 * t281 - t283 * t315) + t62 * (-qJD(1) * t283 + t162 * t579 - t368 * t215 - t254 * t266 - t265 * t553 - t314 * t330);
t335 = t509 * qJD(3);
t308 = t505 * t400;
t196 = rSges(6,1) * t289 - rSges(6,2) * t290;
t195 = -rSges(6,1) * t287 - rSges(6,2) * t288;
t181 = -t399 * t475 + t639;
t173 = t181 * qJD(1);
t123 = qJD(3) * t478 + qJD(2);
t101 = -t572 - t335 * t582 + (-t184 - t312 + t558) * qJD(1);
t100 = -t335 * t379 + (t183 - t556) * qJD(1) + t519;
t90 = t421 * t398 - t399 * t716;
t89 = t398 * t716 + t421 * t399;
t77 = -qJD(3) * t479 + t176 * t411 + t178 * t408;
t76 = -t480 * qJD(3) + t177 * t411 + t179 * t408;
t75 = t430 * qJD(3);
t68 = -t294 * t315 + t310 * t327 + (-t146 + t616) * qJD(1) + t444;
t67 = qJD(1) * t145 - t294 * t314 - t311 * t327 + t420;
t58 = t173 + t456;
t57 = t457 + t576;
t39 = t145 * t315 + t146 * t314 + t241 * t311 - t242 * t310 + t446;
t15 = t160 * t198 - t162 * t197 + t171 * t315 + t172 * t314 + t254 * t88 + t255 * t87 + t282 * t311 - t284 * t310 + t446;
t1 = [(t705 + t706) * t27 + m(3) * ((-t318 * t415 - t573) * t726 + (-t572 + (-0.2e1 * t536 - t403 + t726) * t415) * (-t318 - t693)) + (-qJD(3) * t475 + t332 * t411 + t333 * t408 + t400 * t526 + t401 * t527) * qJD(1) + (t124 + t114) * t704 + (t70 + t71) * t701 + (t117 + (t97 + (t236 + t647) * t399 + t512 + t615) * t315 + (-t399 * t529 - t734 + t96) * t314) * t699 + (-t577 + (t99 + t734) * t315 + (t398 * t529 - t200 + t98) * t314 + ((t236 + t482) * t314 + t529 * t315) * t399 + t54) * t702 - (t76 + t90 + t58) * t582 / 0.2e1 + (t77 + t89) * t379 / 0.2e1 + (t173 + ((t104 - t227 + (t249 + t644) * t399 + t610) * t399 + t609 * t398) * qJD(3)) * t542 + (t125 + t115) * t703 + ((t148 + t181) * t399 + (t147 + t180) * t398) * t575 / 0.2e1 + (t29 * (-t507 + t747) + t61 * (-t508 + t594) + t28 * (t343 + t403 + t740) + t62 * (-pkin(4) * t566 + t317 - t523 + t563) + (t29 * t454 - t28 * t413 + t61 * (t401 * t694 + t691) * t404) * t398 + ((-t409 * t62 - t412 * t61) * pkin(1) + (t61 * (-t395 + t454) - t62 * t413) * t399 + (t400 * t694 - t395 - t690) * t674) * qJD(1) - (-t61 + (-t282 - t693) * qJD(1) + t729 + t732) * t62) * m(6) + (t57 - t576 + ((t399 * t528 + t106 - t609) * t399 + (t398 * t528 + t105 + t530) * t398) * qJD(3)) * t545 + (t69 + t72 + t55) * t700 + (t101 * (t398 * t548 + t392 + t591 - t693) + t100 * t728 + t127 * (t377 + t595) + (t358 * t650 - t648) * qJD(3) + ((-t126 * t412 - t127 * t409) * pkin(1) + (-pkin(2) - t509) * t649 + (t126 * (-rSges(4,3) - pkin(6)) + t127 * t548) * t398) * qJD(1) - (-t556 - t126 - t313 + (-t256 - t693) * qJD(1)) * t127) * m(4) + (t68 * (-t727 + t746) + (-t684 * t399 + t403 + t525 + t596) * t67 + (rSges(5,1) * t568 + rSges(5,2) * t567 + t594 + (-t386 - t403 + (-t328 - t395) * t399) * qJD(1)) * t94 + (t94 - t455 - t729 + t598 + (-rSges(5,1) * t625 - rSges(5,2) * t624 - t570) * t399 + (t241 + t693 + t746) * qJD(1)) * t95) * m(5) - t681 / 0.2e1 + t680 / 0.2e1 + t672 / 0.2e1 + t673 / 0.2e1 + t671 + t38 * t706 + t37 * t707 + t93 * t709 + t92 * t710; m(4) * t75 + m(5) * t39 + m(6) * t15; ((t408 * t593 + t411 * t592) * qJD(1) + ((t398 * t602 - t399 * t603) * t411 + (t398 * t604 + t399 * t605) * t408) * qJD(3)) * t689 + t419 + (t398 * t77 - t399 * t76 + (t147 * t398 + t148 * t399) * qJD(1)) * t688 + ((-t379 * t638 + t588) * t398 + (t435 + (-t715 * t399 + (t639 + t445) * t398) * qJD(3)) * t399) * t545 + ((-t582 * t639 - t588) * t399 + (t435 + (t445 * t398 + (t638 - t715) * t399) * qJD(3)) * t398) * t542 + (qJD(1) * t89 + (t398 * (t398 * t719 + t422 * t399) - t399 * (t398 * t718 + t423 * t399) + (t105 * t398 + t106 * t399) * qJD(1)) * t731) * t696 + (qJD(1) * t90 + (t398 * (t422 * t398 - t399 * t719) - t399 * (t423 * t398 - t399 * t718) + (t103 * t398 + t104 * t399) * qJD(1)) * t731) * t695 + (t57 + t457) * t547 + (t456 + t58) * t546 + (t15 * (t521 + t612) + t56 * (t473 + t562) + (t62 * t474 + (-t231 - t740) * t652) * t398 - (-t62 * t559 + (-t411 * t504 + t535 * t56) * qJD(3)) * pkin(3) - t416 + t714 * (-t736 - t692) + (t399 * t474 + t725) * t61) * m(6) + (t39 * (t611 + t612) + t78 * (t562 + t564) + (t513 * t94 + (qJD(1) * t95 + t68) * t538) * t399 + (t67 * t538 + t95 * t513 + (t94 * t327 + t78 * (-t231 - t242)) * qJD(1)) * t398 - (-t95 * t559 + (t411 * t497 + t535 * t78) * qJD(3)) * pkin(3) + t722) * m(5) + (-(t126 * t303 - t648) * qJD(1) - (t123 * (-t303 * t398 - t304 * t399) + t486 * t509) * qJD(3) + t75 * t478 + t123 * t430 + t486 * t335 + (-t100 * t398 - t101 * t399 + (-t127 * t399 + t650) * qJD(1)) * t358) * m(4); t419 + (t15 * t521 + t56 * t473 + (t617 * t62 - t652 * t740) * t398 - t416 - t714 * t736 + (t399 * t617 + t725) * t61) * m(6) + (t39 * t611 + t78 * (-t242 * t587 + t564) + t497 * t294 + (-t67 * t398 - t68 * t399 + (t398 * t94 - t399 * t95) * qJD(1)) * t327 + t722) * m(5); -t27 * t561 / 0.2e1 + t3 * t631 / 0.2e1 + (t400 * t500 - t401 * t93) * t709 + ((t404 * t500 - t37) * t401 + (-qJD(1) * t501 + t16 * t398 + t17 * t399 + t404 * t93) * t400) * t707 + t400 * t26 * t546 + t4 * t637 / 0.2e1 + (t400 * t502 - t401 * t92) * t710 + ((t404 * t502 - t38) * t401 + (-qJD(1) * t503 + t18 * t398 + t19 * t399 + t404 * t92) * t400) * t706 + t32 * t551 - t401 * (t671 + t672 + t673 + t680 - t681) / 0.2e1 + (-t107 * t401 + t400 * t498) * t520 + ((t404 * t498 - t47) * t401 + (-qJD(1) * t499 + t107 * t404 + t20 * t398 + t21 * t399) * t400) * t697 + (t289 * t712 + t432 * t290 - t399 * t443) * t708 + (-t287 * t712 + t288 * t432 - t398 * t443) * t705 + (t448 * t401 + (-t407 * t712 + t410 * t432) * t400) * t698 + t733 * t624 / 0.2e1 + ((t29 * t160 - t28 * t162 + t61 * t88 - t62 * t87 + (t56 * t483 + (t398 * t61 - t399 * t62) * t265) * t404) * t401 + (t61 * (-t160 * t404 + t168 * t398) + t62 * (t162 * t404 - t168 * t399) + t15 * t483 + t56 * (-t160 * t587 - t162 * t586 - t398 * t87 + t399 * t88) + (qJD(1) * t504 - t28 * t399 + t29 * t398) * t265) * t400 - t61 * (-t195 * t368 - t255 * t308) - t62 * (t196 * t368 - t254 * t308) - t56 * (t195 * t254 + t196 * t255)) * m(6);];
tauc = t1(:);

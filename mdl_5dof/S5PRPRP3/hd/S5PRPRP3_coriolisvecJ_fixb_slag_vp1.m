% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:40
% EndTime: 2019-12-05 15:33:27
% DurationCPUTime: 33.61s
% Computational Cost: add. (14962->644), mult. (22009->949), div. (0->0), fcn. (20917->8), ass. (0->358)
t637 = Icges(5,5) + Icges(6,5);
t610 = Icges(5,6) + Icges(6,6);
t281 = qJ(2) + pkin(8);
t278 = cos(t281);
t283 = cos(pkin(7));
t288 = cos(qJ(4));
t472 = t283 * t288;
t282 = sin(pkin(7));
t286 = sin(qJ(4));
t475 = t282 * t286;
t246 = -t278 * t475 - t472;
t473 = t283 * t286;
t474 = t282 * t288;
t247 = t278 * t474 - t473;
t277 = sin(t281);
t479 = t277 * t282;
t101 = Icges(5,5) * t247 + Icges(5,6) * t246 + Icges(5,3) * t479;
t99 = Icges(6,5) * t247 + Icges(6,6) * t246 + Icges(6,3) * t479;
t601 = t99 + t101;
t248 = -t278 * t473 + t474;
t249 = t278 * t472 + t475;
t478 = t277 * t283;
t100 = Icges(6,5) * t249 + Icges(6,6) * t248 + Icges(6,3) * t478;
t102 = Icges(5,5) * t249 + Icges(5,6) * t248 + Icges(5,3) * t478;
t598 = t100 + t102;
t491 = Icges(6,4) * t247;
t103 = Icges(6,2) * t246 + Icges(6,6) * t479 + t491;
t495 = Icges(5,4) * t247;
t105 = Icges(5,2) * t246 + Icges(5,6) * t479 + t495;
t623 = t103 + t105;
t490 = Icges(6,4) * t249;
t104 = Icges(6,2) * t248 + Icges(6,6) * t478 + t490;
t494 = Icges(5,4) * t249;
t106 = Icges(5,2) * t248 + Icges(5,6) * t478 + t494;
t622 = t104 + t106;
t220 = Icges(6,4) * t246;
t107 = Icges(6,1) * t247 + Icges(6,5) * t479 + t220;
t222 = Icges(5,4) * t246;
t109 = Icges(5,1) * t247 + Icges(5,5) * t479 + t222;
t621 = t107 + t109;
t221 = Icges(6,4) * t248;
t108 = Icges(6,1) * t249 + Icges(6,5) * t478 + t221;
t223 = Icges(5,4) * t248;
t110 = Icges(5,1) * t249 + Icges(5,5) * t478 + t223;
t620 = t108 + t110;
t585 = Icges(5,1) + Icges(6,1);
t636 = Icges(5,2) + Icges(6,2);
t645 = Icges(5,3) + Icges(6,3);
t368 = Icges(6,5) * t288 - Icges(6,6) * t286;
t166 = -Icges(6,3) * t278 + t277 * t368;
t369 = Icges(5,5) * t288 - Icges(5,6) * t286;
t168 = -Icges(5,3) * t278 + t277 * t369;
t635 = -t168 - t166;
t488 = Icges(6,4) * t288;
t372 = -Icges(6,2) * t286 + t488;
t170 = -Icges(6,6) * t278 + t277 * t372;
t492 = Icges(5,4) * t288;
t373 = -Icges(5,2) * t286 + t492;
t172 = -Icges(5,6) * t278 + t277 * t373;
t584 = t170 + t172;
t489 = Icges(6,4) * t286;
t378 = Icges(6,1) * t288 - t489;
t174 = -Icges(6,5) * t278 + t277 * t378;
t493 = Icges(5,4) * t286;
t379 = Icges(5,1) * t288 - t493;
t176 = -Icges(5,5) * t278 + t277 * t379;
t593 = t174 + t176;
t644 = (-t637 * t286 - t610 * t288) * t277;
t607 = t622 * t246 + t620 * t247 + t598 * t479;
t606 = t623 * t248 + t621 * t249 + t601 * t478;
t608 = t623 * t246 + t621 * t247 + t479 * t601;
t567 = t248 * t622 + t249 * t620 + t478 * t598;
t643 = t246 * t584 + t247 * t593 - t479 * t635;
t642 = t248 * t584 + t249 * t593 - t478 * t635;
t452 = qJD(2) * t277;
t434 = t286 * t452;
t150 = -qJD(4) * t247 + t282 * t434;
t433 = t288 * t452;
t151 = qJD(4) * t246 - t282 * t433;
t450 = qJD(2) * t282;
t432 = t278 * t450;
t641 = t150 * t610 + t151 * t637 + t432 * t645;
t152 = -qJD(4) * t249 + t283 * t434;
t153 = qJD(4) * t248 - t283 * t433;
t448 = qJD(2) * t283;
t431 = t278 * t448;
t640 = t152 * t610 + t153 * t637 + t431 * t645;
t167 = Icges(6,3) * t277 + t278 * t368;
t169 = Icges(5,3) * t277 + t278 * t369;
t639 = t644 * qJD(4) + (t167 + t169) * qJD(2);
t638 = Icges(5,4) + Icges(6,4);
t171 = Icges(6,6) * t277 + t278 * t372;
t173 = Icges(5,6) * t277 + t278 * t373;
t634 = t171 + t173;
t175 = Icges(6,5) * t277 + t278 * t378;
t177 = Icges(5,5) * t277 + t278 * t379;
t561 = -t175 - t177;
t633 = (-t288 * t636 - t489 - t493) * t277;
t632 = (t286 * t585 + t488 + t492) * t277;
t631 = t607 * t283;
t630 = t606 * t282;
t629 = t636 * t150 + t638 * t151 + t610 * t432;
t628 = t636 * t152 + t638 * t153 + t610 * t431;
t627 = t638 * t150 + t585 * t151 + t637 * t432;
t626 = t638 * t152 + t585 * t153 + t637 * t431;
t625 = t634 * qJD(2) + t633 * qJD(4);
t624 = t561 * qJD(2) + t632 * qJD(4);
t451 = qJD(2) * t278;
t619 = -t639 * t277 + t635 * t451;
t618 = t640 * t277 + t598 * t451;
t617 = t641 * t277 + t601 * t451;
t287 = sin(qJ(2));
t289 = cos(qJ(2));
t589 = -Icges(3,5) * t287 - Icges(4,5) * t277 - Icges(3,6) * t289 - Icges(4,6) * t278;
t616 = t567 * t283 + t630;
t615 = t608 * t282 + t631;
t614 = t642 * t277;
t613 = t643 * t277;
t445 = qJD(4) * t277;
t429 = t283 * t445;
t260 = t429 + t450;
t430 = t282 * t445;
t261 = t430 - t448;
t609 = -t260 * t282 + t261 * t283;
t528 = t260 / 0.2e1;
t526 = t261 / 0.2e1;
t573 = t623 * t150 + t621 * t151 + t629 * t246 + t627 * t247 + t617 * t282;
t572 = t622 * t150 + t620 * t151 + t628 * t246 + t626 * t247 + t618 * t282;
t571 = t623 * t152 + t621 * t153 + t629 * t248 + t627 * t249 + t617 * t283;
t570 = t622 * t152 + t620 * t153 + t628 * t248 + t626 * t249 + t618 * t283;
t366 = -t103 * t286 + t107 * t288;
t47 = t277 * t366 - t278 * t99;
t364 = -t105 * t286 + t109 * t288;
t49 = -t101 * t278 + t277 * t364;
t605 = t47 + t49;
t365 = -t104 * t286 + t108 * t288;
t48 = -t100 * t278 + t277 * t365;
t363 = -t106 * t286 + t110 * t288;
t50 = -t102 * t278 + t277 * t363;
t604 = t48 + t50;
t361 = -t170 * t286 + t174 * t288;
t481 = t166 * t278;
t66 = t277 * t361 - t481;
t360 = -t172 * t286 + t176 * t288;
t480 = t168 * t278;
t67 = t277 * t360 - t480;
t586 = t66 + t67;
t592 = (-t584 * t152 - t593 * t153 - t625 * t248 + t624 * t249 + t619 * t283) * t278 + (t616 * t278 + t614) * qJD(2);
t591 = (-t584 * t150 - t593 * t151 - t625 * t246 + t624 * t247 + t619 * t282) * t278 + (t615 * t278 + t613) * qJD(2);
t590 = t589 * qJD(2);
t499 = Icges(3,4) * t287;
t376 = -Icges(3,2) * t289 - t499;
t498 = Icges(3,4) * t289;
t382 = -Icges(3,1) * t287 - t498;
t496 = Icges(4,4) * t278;
t497 = Icges(4,4) * t277;
t588 = -(-t287 * t376 + t289 * t382) * qJD(2) + (-Icges(4,2) * t278 - t497) * t452 - (-Icges(4,1) * t277 - t496) * t451;
t444 = qJD(4) * t278;
t539 = t260 * (-t249 * t636 + t221 + t223 + t620) + t261 * (-t247 * t636 + t220 + t222 + t621) - t444 * (t593 + t633);
t587 = t571 * t526 + t570 * t528 + t592 * qJD(4) / 0.2e1;
t279 = t282 ^ 2;
t280 = t283 ^ 2;
t453 = t279 + t280;
t583 = t644 * t444 + (-t246 * t637 + t610 * t247) * t261 + (-t248 * t637 + t610 * t249) * t260;
t389 = t49 * t282 + t50 * t283;
t390 = t47 * t282 + t48 * t283;
t582 = t389 + t390;
t578 = qJD(4) * t591 + t260 * t572 + t261 * t573;
t546 = ((t364 + t366) * qJD(2) - t641) * t278 + (t627 * t288 - t629 * t286 + (-t286 * t621 - t288 * t623) * qJD(4) + t601 * qJD(2)) * t277;
t545 = ((t363 + t365) * qJD(2) - t640) * t278 + (t626 * t288 - t628 * t286 + (-t286 * t620 - t288 * t622) * qJD(4) + t598 * qJD(2)) * t277;
t576 = rSges(6,1) + pkin(4);
t575 = t260 * t607 + t261 * t608 - t444 * t643;
t574 = t260 * t567 + t261 * t606 - t444 * t642;
t569 = t260 * t604 + t261 * t605 - t444 * t586;
t138 = t170 * t282;
t140 = t172 * t282;
t566 = -t138 - t140;
t139 = t170 * t283;
t141 = t172 * t283;
t565 = -t139 - t141;
t142 = t174 * t282;
t144 = t176 * t282;
t564 = -t142 - t144;
t143 = t174 * t283;
t145 = t176 * t283;
t563 = -t143 - t145;
t519 = pkin(4) * t288;
t155 = qJ(5) * t277 + t278 * t519;
t397 = rSges(6,1) * t288 - rSges(6,2) * t286;
t462 = rSges(6,3) * t277 + t278 * t397 + t155;
t560 = t589 * t282;
t559 = t589 * t283;
t558 = t590 * t282;
t557 = t590 * t283;
t341 = t169 - t360;
t342 = t167 - t361;
t533 = -(-t166 * t283 - t365) * t260 - (-t166 * t282 - t366) * t261;
t534 = -(-t168 * t283 - t363) * t260 - (-t168 * t282 - t364) * t261;
t556 = (-t533 - t534 + (-t341 - t342) * t444) * t277;
t377 = -Icges(3,2) * t287 + t498;
t209 = Icges(3,6) * t282 + t283 * t377;
t383 = Icges(3,1) * t289 - t499;
t211 = Icges(3,5) * t282 + t283 * t383;
t375 = -Icges(4,2) * t277 + t496;
t381 = Icges(4,1) * t278 - t497;
t547 = -(Icges(4,6) * t282 + t283 * t375) * t278 - (Icges(4,5) * t282 + t283 * t381) * t277;
t555 = t588 * t283 + (t209 * t289 + t211 * t287 - t547) * qJD(2);
t208 = -Icges(3,6) * t283 + t282 * t377;
t210 = -Icges(3,5) * t283 + t282 * t383;
t548 = (-Icges(4,6) * t283 + t282 * t375) * t278 + (-Icges(4,5) * t283 + t282 * t381) * t277;
t554 = t588 * t282 + (t208 * t289 + t210 * t287 + t548) * qJD(2);
t325 = qJ(5) * t278 - t277 * t519;
t553 = rSges(6,3) * t278 - t277 * t397 + t325;
t552 = t260 * t598 + t261 * t601;
t551 = -t480 - t481;
t540 = t586 * t452 + (t639 * t278 + (t624 * t288 + t625 * t286 + (t286 * t593 + t584 * t288) * qJD(4)) * t277 + ((-t360 - t361) * t278 + t635 * t277 + t582) * qJD(2)) * t278;
t538 = (t584 + t632) * t444 + (t246 * t585 - t491 - t495 - t623) * t261 + (t248 * t585 - t490 - t494 - t622) * t260;
t537 = t583 * t277;
t271 = rSges(3,1) * t287 + rSges(3,2) * t289;
t337 = qJD(2) * t271;
t398 = rSges(5,1) * t288 - rSges(5,2) * t286;
t181 = -rSges(5,3) * t278 + t277 * t398;
t267 = pkin(3) * t278 + pkin(6) * t277;
t226 = t267 * t282;
t227 = t267 * t283;
t521 = pkin(2) * t289;
t184 = -qJ(3) * t283 + t282 * t521;
t185 = qJ(3) * t282 + t283 * t521;
t426 = t184 * t450 + t185 * t448 + qJD(1);
t384 = t226 * t450 + t227 * t448 + t426;
t442 = qJD(5) * t278;
t500 = rSges(6,1) * t249 + rSges(6,2) * t248 + rSges(6,3) * t478 + pkin(4) * t475 + t155 * t283;
t501 = rSges(6,1) * t247 + rSges(6,2) * t246 + rSges(6,3) * t479 - pkin(4) * t473 + t155 * t282;
t27 = t260 * t501 - t261 * t500 + t384 - t442;
t266 = pkin(3) * t277 - pkin(6) * t278;
t339 = qJD(2) * t266;
t206 = t282 * t339;
t207 = t283 * t339;
t511 = pkin(2) * qJD(2);
t439 = t287 * t511;
t446 = qJD(3) * t283;
t262 = -t282 * t439 - t446;
t276 = qJD(3) * t282;
t263 = -t283 * t439 + t276;
t455 = t262 * t450 + t263 * t448;
t410 = -t206 * t450 - t207 * t448 + t455;
t443 = qJD(5) * t277;
t270 = t283 * t443;
t510 = pkin(4) * qJD(4);
t437 = t286 * t510;
t295 = qJD(2) * t325 - t278 * t437;
t436 = t288 * t510;
t515 = rSges(6,1) * t153 + rSges(6,2) * t152 + rSges(6,3) * t431 + t282 * t436 + t283 * t295 + t270;
t269 = t282 * t443;
t516 = rSges(6,1) * t151 + rSges(6,2) * t150 + rSges(6,3) * t432 + t282 * t295 - t283 * t436 + t269;
t5 = -t515 * t261 + t516 * t260 + (t443 + (-t282 * t500 + t283 * t501) * t444) * qJD(2) + t410;
t535 = t27 * t515 + t5 * t500;
t532 = -t287 * (t376 * t282 + t210) - t289 * (-t382 * t282 + t208);
t290 = qJD(2) ^ 2;
t531 = -m(6) / 0.2e1;
t530 = m(6) / 0.2e1;
t529 = -t260 / 0.2e1;
t527 = -t261 / 0.2e1;
t522 = pkin(2) * t287;
t236 = (-rSges(6,1) * t286 - rSges(6,2) * t288) * t277;
t514 = qJD(2) * t462 + qJD(4) * t236 - t277 * t437 - t442;
t477 = t278 * t282;
t476 = t278 * t283;
t467 = t553 * t282;
t466 = t553 * t283;
t465 = -rSges(6,2) * t247 + t246 * t576;
t464 = rSges(6,2) * t249 - t248 * t576;
t461 = t282 * t184 + t283 * t185;
t458 = t453 * t339;
t454 = t282 * t262 + t283 * t263;
t441 = qJD(2) * qJD(4);
t440 = t290 * t521;
t438 = t289 * t511;
t423 = t448 / 0.2e1;
t421 = -t444 / 0.2e1;
t420 = t444 / 0.2e1;
t264 = rSges(4,1) * t277 + rSges(4,2) * t278;
t419 = -t264 - t522;
t265 = rSges(4,1) * t278 - rSges(4,2) * t277;
t418 = -t265 - t521;
t417 = -t266 - t522;
t416 = t441 / 0.2e1;
t415 = t453 * t287;
t414 = t282 * t440;
t413 = t283 * t440;
t253 = t267 * qJD(2);
t412 = -t253 + t442;
t411 = t282 * t226 + t283 * t227 + t461;
t409 = -t282 * t206 - t283 * t207 + t454;
t404 = pkin(4) * t277 * t286 - t236;
t403 = -t181 + t417;
t402 = t277 * t416;
t401 = t278 * t416;
t252 = t265 * qJD(2);
t400 = -t252 - t438;
t399 = -t253 - t438;
t272 = rSges(3,1) * t289 - rSges(3,2) * t287;
t114 = rSges(5,1) * t249 + rSges(5,2) * t248 + rSges(5,3) * t478;
t86 = rSges(5,1) * t153 + rSges(5,2) * t152 + rSges(5,3) * t431;
t237 = (-rSges(5,1) * t286 - rSges(5,2) * t288) * t277;
t320 = rSges(5,3) * t277 + t278 * t398;
t98 = qJD(2) * t320 + qJD(4) * t237;
t34 = -t414 - t86 * t444 - t260 * t98 + (-t253 * t282 + (t114 * t277 - t181 * t476) * qJD(4)) * qJD(2);
t112 = rSges(5,1) * t247 + rSges(5,2) * t246 + rSges(5,3) * t479;
t84 = rSges(5,1) * t151 + rSges(5,2) * t150 + rSges(5,3) * t432;
t35 = -t413 + t84 * t444 + t261 * t98 + (-t253 * t283 + (-t112 * t277 + t181 * t477) * qJD(4)) * qJD(2);
t396 = t282 * t35 - t283 * t34;
t386 = qJD(2) * t417;
t321 = t282 * t386 - t446;
t36 = t260 * t553 - t444 * t500 + t269 + t321;
t327 = t283 * t386 + t276;
t37 = -t261 * t553 + t444 * t501 + t270 + t327;
t395 = t282 * t36 + t283 * t37;
t55 = -t114 * t444 - t181 * t260 + t321;
t56 = t112 * t444 + t181 * t261 + t327;
t388 = -t282 * t55 - t283 * t56;
t387 = qJD(2) * t419;
t362 = t112 * t283 - t114 * t282;
t359 = t453 * t272;
t358 = t453 * t337;
t357 = t399 - t98;
t356 = t417 + t553;
t355 = t282 * t401;
t354 = t283 * t401;
t353 = -qJD(2) * t252 - t440;
t340 = t399 - t514;
t190 = -rSges(4,3) * t283 + t265 * t282;
t191 = rSges(4,3) * t282 + t265 * t283;
t65 = (t190 * t282 + t191 * t283) * qJD(2) + t426;
t338 = t65 * t264;
t336 = qJD(2) * t264;
t38 = t112 * t260 - t114 * t261 + t384;
t335 = t38 * t362;
t328 = t27 * t516 + t5 * t501;
t326 = t36 * t500 - t37 * t501;
t324 = (t382 * t283 - t209) * t289 + (-t376 * t283 - t211) * t287;
t294 = (t27 * t501 + t36 * t553) * t283 + (-t27 * t500 - t37 * t553) * t282;
t293 = t282 * t547 + t283 * t548;
t205 = t283 * t336;
t204 = t282 * t336;
t159 = t353 * t283;
t158 = t353 * t282;
t157 = t283 * t387 + t276;
t156 = t282 * t387 - t446;
t133 = rSges(5,1) * t248 - rSges(5,2) * t249;
t131 = rSges(5,1) * t246 - rSges(5,2) * t247;
t115 = t358 * qJD(2);
t90 = qJD(2) * t359 + qJD(1);
t68 = (-t204 * t282 - t205 * t283) * qJD(2) + t455;
t26 = t278 * t362 * t441 + t260 * t84 - t261 * t86 + t410;
t15 = -t413 + t514 * t261 + t516 * t444 + (t412 * t283 + (-t277 * t501 - t477 * t553) * qJD(4)) * qJD(2);
t14 = -t414 - t514 * t260 - t515 * t444 + (t412 * t282 + (t277 * t500 + t476 * t553) * qJD(4)) * qJD(2);
t1 = [-m(3) * t115 + m(4) * t68 + m(5) * t26 + m(6) * t5; (((-t248 * t634 + t249 * t561 + t630) * t278 + t614) * qJD(4) + (((t551 + t567) * qJD(4) + t552) * t278 + t556) * t283 + (t248 * t566 + t249 * t564) * t261 + (t248 * t565 + t249 * t563) * t260) * t529 + (t282 * t570 - t283 * t571) * t528 + (((-t246 * t634 + t247 * t561 + t631) * t278 + t613) * qJD(4) + (((t551 + t608) * qJD(4) + t552) * t278 + t556) * t282 + (t246 * t566 + t247 * t564) * t261 + (t246 * t565 + t247 * t563) * t260) * t527 + (t282 * t572 - t283 * t573) * t526 + t282 * t587 - t578 * t283 / 0.2e1 + (t554 * t280 + (t557 * t282 + (-t555 - t558) * t283) * t282) * t450 + (-t558 * t280 + (t555 * t282 + (-t554 + t557) * t283) * t282) * t448 - (t559 * qJD(2) * t279 + (-t532 * t283 + t293 + (t324 - t560) * t282) * t448) * t450 / 0.2e1 + ((t324 * t282 + t293 + (-t532 - t559) * t283) * t450 + t560 * qJD(2) * t280) * t423 - t569 * t445 / 0.2e1 + (((t139 * t286 - t143 * t288 + t100) * t260 + (t138 * t286 - t142 * t288 + t99) * t261 + t66 * qJD(4)) * t277 + ((t342 * t278 + (t171 * t286 - t175 * t288 - t166) * t277 + t390) * qJD(4) + t533) * t278 + ((t141 * t286 - t145 * t288 + t102) * t260 + (t140 * t286 - t144 * t288 + t101) * t261 + t67 * qJD(4)) * t277 + ((t341 * t278 + (t173 * t286 - t177 * t288 - t168) * t277 + t389) * qJD(4) + t534) * t278) * t420 + (t282 * t604 - t283 * t605) * t402 + (t282 * t607 - t283 * t608) * t355 + (t567 * t282 - t283 * t606) * t354 + (t5 * t411 + t27 * t409 + (t15 * t356 + t340 * t37 + t535) * t283 + (t14 * t356 + t340 * t36 + t328) * t282 + t27 * t458 - (t27 * t277 + t278 * t395) * qJD(5) - (-t27 * t466 + t37 * t462) * t261 - (t27 * t467 - t36 * t462) * t260 - (-t395 * t267 + (-t27 * t415 - t289 * t395) * pkin(2)) * qJD(2) - (t326 * t277 + (-t36 * t466 + t37 * t467 + t294) * t278) * qJD(4)) * m(6) + (t26 * t411 + t38 * t409 + (t26 * t114 + t35 * t403 + t357 * t56 + t38 * t86) * t283 + (t26 * t112 + t34 * t403 + t357 * t55 + t38 * t84) * t282 - t38 * (t181 * t609 - t458) - (-t260 * t55 + t261 * t56) * t320 - (t388 * t267 + (t289 * t388 - t38 * t415) * pkin(2)) * qJD(2) - ((-t112 * t56 + t114 * t55) * t277 + t335 * t278) * qJD(4)) * m(5) + (t68 * t461 + t65 * t454 + (t157 * t400 + t159 * t419 + t68 * t191 - t65 * t205) * t283 + (t156 * t400 + t158 * t419 + t68 * t190 - t65 * t204) * t282 - (-t65 * pkin(2) * t415 + (t157 * t418 - t283 * t338) * t283 + (t156 * t418 - t282 * t338) * t282) * qJD(2)) * m(4) + (-t115 * t359 - t358 * t90 + (t271 * t272 * t290 + t337 * t90) * t453) * m(3) + ((-t546 + t574) * t283 + (t545 + t575) * t282) * t421; m(4) * (-t158 * t283 + t159 * t282) + m(5) * t396 + m(6) * (-t14 * t283 + t15 * t282); (t248 * t539 + t249 * t538 - t283 * t537) * t529 + ((t282 * t571 + t283 * t570) * t277 + t592) * t528 + (t246 * t539 + t247 * t538 - t282 * t537) * t527 + ((t282 * t573 + t283 * t572) * t277 + t591) * t526 - (qJD(4) * t540 + t260 * t545 + t261 * t546) * t278 / 0.2e1 + t578 * t479 / 0.2e1 + t478 * t587 + t569 * t452 / 0.2e1 + ((t282 * t546 + t283 * t545) * t277 + t540) * t421 + (t583 * t278 + (-t286 * t539 + t538 * t288) * t277) * t420 + t575 * t432 / 0.2e1 + t574 * t278 * t423 + (t277 * t582 - t278 * t586) * t402 + (t615 * t277 - t278 * t643) * t355 + (t616 * t277 - t278 * t642) * t354 + (-(t27 * t464 - t37 * t404) * t261 - (t27 * t465 + t36 * t404) * t260 - (t36 * t464 + t37 * t465) * t444 + (qJD(2) * t294 - t14 * t500 + t15 * t501 - t36 * t515 + t37 * t516) * t278 + (t326 * qJD(2) + (t14 * t553 - t36 * t514 + t328) * t283 + (-t15 * t553 + t37 * t514 - t535) * t282) * t277) * m(6) + (-t56 * (t131 * t444 + t237 * t261) - t55 * (-t133 * t444 - t237 * t260) - t38 * (t131 * t260 - t133 * t261) + (t35 * t112 - t34 * t114 - t55 * t86 + t56 * t84 + (t335 + (t282 * t56 - t283 * t55) * t181) * qJD(2)) * t278 + (t56 * (-qJD(2) * t112 + t282 * t98) + t55 * (qJD(2) * t114 - t283 * t98) + t26 * t362 + t38 * (-t282 * t86 + t283 * t84) + t396 * t181) * t277) * m(5); 0.2e1 * ((qJD(2) * t27 + t14 * t282 + t15 * t283) * t530 - t27 * t609 * t531) * t277 + 0.2e1 * ((t36 * t450 + t37 * t448 - t5) * t530 + (t37 * (-t261 + t430) + t36 * (t260 - t429)) * t531) * t278;];
tauc = t1(:);

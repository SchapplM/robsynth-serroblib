% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR10_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:10
% EndTime: 2019-12-31 17:11:51
% DurationCPUTime: 36.90s
% Computational Cost: add. (8273->804), mult. (22065->1069), div. (0->0), fcn. (20136->6), ass. (0->398)
t701 = Icges(4,4) - Icges(3,5);
t700 = Icges(4,5) - Icges(3,6);
t699 = Icges(4,1) + Icges(3,3);
t342 = sin(qJ(2));
t345 = cos(qJ(2));
t669 = t700 * t342 - t701 * t345;
t346 = cos(qJ(1));
t698 = t699 * t346;
t343 = sin(qJ(1));
t546 = t343 * t345;
t549 = t342 * t343;
t693 = t701 * t546 - t700 * t549 + t698;
t679 = t699 * t343 + t669 * t346;
t575 = Icges(3,4) * t342;
t278 = Icges(3,1) * t345 - t575;
t564 = Icges(4,6) * t342;
t416 = Icges(4,2) * t345 - t564;
t697 = t278 + t416;
t563 = Icges(4,6) * t345;
t414 = -Icges(4,3) * t342 + t563;
t331 = Icges(3,4) * t345;
t421 = -Icges(3,2) * t342 + t331;
t696 = t414 + t421;
t271 = Icges(3,5) * t342 + Icges(3,6) * t345;
t419 = Icges(4,4) * t342 + Icges(4,5) * t345;
t685 = -t419 + t271;
t565 = Icges(3,6) * t346;
t181 = Icges(3,4) * t546 - Icges(3,2) * t549 - t565;
t307 = Icges(4,6) * t549;
t573 = Icges(4,4) * t346;
t190 = Icges(4,2) * t546 - t307 + t573;
t695 = t181 * t342 - t190 * t345;
t186 = Icges(3,5) * t343 + t278 * t346;
t187 = Icges(4,5) * t343 - t346 * t414;
t694 = -t186 * t546 - t187 * t549;
t275 = Icges(3,2) * t345 + t575;
t413 = Icges(4,3) * t345 + t564;
t688 = -t275 - t413;
t277 = Icges(3,1) * t342 + t331;
t415 = Icges(4,2) * t342 + t563;
t692 = t277 + t415;
t402 = t275 * t342 - t277 * t345;
t686 = -t342 * t413 + t345 * t415 - t402;
t312 = Icges(3,4) * t549;
t569 = Icges(3,5) * t346;
t185 = Icges(3,1) * t546 - t312 - t569;
t568 = Icges(4,5) * t346;
t188 = Icges(4,6) * t546 - Icges(4,3) * t549 + t568;
t666 = -t185 * t345 + t188 * t342 + t695;
t501 = qJD(2) * t342;
t499 = qJD(2) * t345;
t182 = Icges(3,6) * t343 + t346 * t421;
t677 = t182 - t187;
t691 = t677 * t345;
t690 = t696 * qJD(2);
t689 = t697 * qJD(2);
t684 = t679 * t346 + t694;
t544 = t345 * t346;
t548 = t342 * t346;
t641 = t186 * t544 + t187 * t548 + t679 * t343;
t683 = -t185 * t544 + t188 * t548 + t693 * t343;
t629 = t685 * t343;
t682 = -t688 * t499 + t692 * t501;
t681 = -t666 * t343 + t693 * t346;
t308 = Icges(4,6) * t548;
t574 = Icges(4,4) * t343;
t189 = -Icges(4,2) * t544 + t308 + t574;
t649 = -t182 * t549 - t189 * t546 - t684;
t648 = -t181 * t548 + t190 * t544 - t683;
t647 = -t182 * t548 - t189 * t544 + t641;
t680 = t686 * t346 + t629;
t678 = t181 + t188;
t676 = t185 + t190;
t675 = t186 - t189;
t670 = t182 * t342 + t189 * t345;
t668 = t689 * t345 - t690 * t342 + (-t342 * t692 + t688 * t345) * qJD(2) + t685 * qJD(1);
t667 = t685 * qJD(2);
t665 = t186 * t345 + t187 * t342 - t670;
t664 = t686 * qJD(1) - t669 * qJD(2);
t663 = t680 * qJD(1);
t662 = (-t667 * t343 + (t666 + t679) * qJD(1)) * t346;
t661 = -t343 * t677 + t346 * t678;
t660 = (t649 * t343 - t346 * t681) * qJD(2);
t659 = (t343 * t647 - t346 * t648) * qJD(2);
t658 = t696 + t692;
t657 = t697 + t688;
t656 = (t307 + t312 + (Icges(3,2) + Icges(4,3)) * t546 - t676) * t346 + (-Icges(4,3) * t544 - t275 * t346 - t308 + t675) * t343;
t502 = qJD(1) * t346;
t234 = t419 * t346;
t492 = (t413 * t549 - t415 * t546 - t234) * qJD(1);
t550 = t271 * t346;
t103 = -t343 * t402 - t550;
t493 = t103 * qJD(1);
t655 = t493 - t492 + t660;
t654 = t659 + t663;
t653 = t666 * qJD(2) + t682 * t343 + (-t691 + (-t346 * t416 - t186 + t574) * t342) * qJD(1);
t652 = t665 * qJD(2) - t682 * t346 + ((t565 - t568) * t345 + (t569 - t573) * t342 + (-t342 * t697 - t345 * t696) * t343) * qJD(1);
t651 = -t664 * t343 + t668 * t346;
t650 = t668 * t343 + t664 * t346;
t646 = t342 * t676 + t345 * t678;
t645 = t342 * t675 + t691;
t642 = t670 + t693;
t640 = -t667 * t346 + (-t669 * t343 - t665 + t698) * qJD(1);
t494 = qJD(4) * t345;
t500 = qJD(2) * t343;
t252 = t346 * t494 + t500;
t498 = qJD(2) * t346;
t253 = -t343 * t494 + t498;
t495 = qJD(4) * t342;
t318 = qJD(1) + t495;
t341 = sin(qJ(4));
t344 = cos(qJ(4));
t545 = t344 * t346;
t221 = -t341 * t343 + t342 * t545;
t547 = t343 * t344;
t222 = t341 * t548 + t547;
t110 = Icges(5,5) * t222 + Icges(5,6) * t221 + Icges(5,3) * t544;
t572 = Icges(5,4) * t222;
t113 = Icges(5,2) * t221 + Icges(5,6) * t544 + t572;
t208 = Icges(5,4) * t221;
t116 = Icges(5,1) * t222 + Icges(5,5) * t544 + t208;
t35 = t110 * t544 + t221 * t113 + t222 * t116;
t223 = t341 * t346 + t342 * t547;
t224 = -t341 * t549 + t545;
t112 = -Icges(5,5) * t224 + Icges(5,6) * t223 + Icges(5,3) * t546;
t210 = Icges(5,4) * t224;
t115 = Icges(5,2) * t223 + Icges(5,6) * t546 - t210;
t209 = Icges(5,4) * t223;
t117 = Icges(5,1) * t224 - Icges(5,5) * t546 - t209;
t36 = t112 * t544 + t221 * t115 - t117 * t222;
t417 = Icges(5,5) * t341 + Icges(5,6) * t344;
t370 = -Icges(5,3) * t342 + t345 * t417;
t571 = Icges(5,4) * t341;
t418 = Icges(5,2) * t344 + t571;
t371 = -Icges(5,6) * t342 + t345 * t418;
t570 = Icges(5,4) * t344;
t422 = Icges(5,1) * t341 + t570;
t372 = -Icges(5,5) * t342 + t345 * t422;
t60 = -t221 * t371 - t222 * t372 - t370 * t544;
t10 = t252 * t35 - t253 * t36 + t60 * t318;
t37 = t110 * t546 + t223 * t113 - t224 * t116;
t38 = t112 * t546 + t115 * t223 + t117 * t224;
t61 = -t223 * t371 + t224 * t372 - t370 * t546;
t11 = t252 * t37 - t253 * t38 + t318 * t61;
t411 = t115 * t344 - t117 * t341;
t42 = t112 * t342 - t345 * t411;
t638 = 0.2e1 * qJD(2);
t634 = -t656 * t342 + t661 * t345;
t633 = (-t658 * t342 + t657 * t345) * qJD(1);
t478 = t342 * t498;
t503 = qJD(1) * t343;
t483 = t345 * t503;
t632 = t478 + t483;
t631 = t669 * qJD(1);
t630 = t550 - t234;
t433 = rSges(5,1) * t224 - rSges(5,2) * t223;
t123 = rSges(5,3) * t546 - t433;
t496 = qJD(3) * t346;
t303 = t342 * t496;
t432 = rSges(5,1) * t341 + rSges(5,2) * t344;
t374 = -rSges(5,3) * t342 + t345 * t432;
t279 = pkin(2) * t342 - qJ(3) * t345;
t458 = -pkin(6) * t342 - t279;
t436 = t458 * t346;
t628 = qJD(2) * t436 - t123 * t318 + t253 * t374 + t303;
t327 = qJD(3) * t342;
t593 = rSges(4,2) * t342;
t431 = rSges(4,3) * t345 + t593;
t626 = (-qJD(2) * t431 - t327) * t343;
t333 = t343 * rSges(4,1);
t204 = -rSges(4,2) * t544 + rSges(4,3) * t548 + t333;
t248 = pkin(2) * t544 + qJ(3) * t548;
t288 = t346 * pkin(1) + t343 * pkin(5);
t443 = t248 + t288;
t625 = t204 + t443;
t559 = qJ(3) * t342;
t283 = pkin(2) * t345 + t559;
t243 = t283 * t343;
t214 = qJD(1) * t243;
t337 = t346 * pkin(5);
t287 = pkin(1) * t343 - t337;
t266 = qJD(1) * t287;
t624 = -t214 - t266;
t320 = pkin(6) * t544;
t264 = t343 * pkin(3) + t320;
t489 = -rSges(5,3) - pkin(2) - pkin(6);
t623 = t489 * t345 - pkin(1) - t559;
t319 = pkin(6) * t546;
t265 = pkin(3) * t346 - t319;
t519 = -t243 - t287;
t39 = (t265 + t519) * qJD(1) + t628;
t121 = t222 * rSges(5,1) + t221 * rSges(5,2) + rSges(5,3) * t544;
t249 = t279 * t500;
t479 = t342 * t500;
t298 = pkin(6) * t479;
t474 = t343 * t327;
t40 = t474 + t121 * t318 + t374 * t252 - t249 - t298 + (t264 + t443) * qJD(1);
t622 = t343 * t40 + t346 * t39;
t232 = (Icges(5,2) * t341 - t570) * t345;
t365 = t252 * (-Icges(5,2) * t222 + t116 + t208) - t253 * (Icges(5,2) * t224 - t117 + t209) + t318 * (-t372 + t232);
t237 = (-Icges(5,1) * t344 + t571) * t345;
t366 = t252 * (-Icges(5,1) * t221 + t113 + t572) - t253 * (-Icges(5,1) * t223 + t115 - t210) + t318 * (-t371 - t237);
t612 = m(4) / 0.2e1;
t611 = m(5) / 0.2e1;
t490 = qJD(2) * qJD(4);
t469 = t342 * t490;
t166 = qJD(1) * t252 - t343 * t469;
t610 = t166 / 0.2e1;
t167 = qJD(1) * t253 - t346 * t469;
t609 = t167 / 0.2e1;
t608 = -t252 / 0.2e1;
t607 = t252 / 0.2e1;
t606 = -t253 / 0.2e1;
t605 = t253 / 0.2e1;
t604 = -t318 / 0.2e1;
t603 = t318 / 0.2e1;
t412 = t113 * t344 + t116 * t341;
t446 = qJD(1) * t342 + qJD(4);
t476 = t345 * t498;
t373 = -t343 * t446 + t476;
t401 = t318 * t341;
t101 = t344 * t373 - t346 * t401;
t400 = t318 * t344;
t102 = t341 * t373 + t346 * t400;
t52 = Icges(5,5) * t102 + Icges(5,6) * t101 - Icges(5,3) * t632;
t54 = Icges(5,4) * t102 + Icges(5,2) * t101 - Icges(5,6) * t632;
t56 = Icges(5,1) * t102 + Icges(5,4) * t101 - Icges(5,5) * t632;
t8 = (qJD(2) * t412 + t52) * t342 + (qJD(2) * t110 - t341 * t56 - t344 * t54 + (t113 * t341 - t116 * t344) * qJD(4)) * t345;
t600 = t8 * t252;
t399 = t446 * t346;
t477 = t343 * t499;
t100 = t343 * t400 + (t399 + t477) * t341;
t482 = t345 * t502;
t381 = -t479 + t482;
t99 = t344 * t399 + (t344 * t499 - t401) * t343;
t51 = Icges(5,5) * t100 + Icges(5,6) * t99 + Icges(5,3) * t381;
t53 = Icges(5,4) * t100 + Icges(5,2) * t99 + Icges(5,6) * t381;
t55 = Icges(5,1) * t100 + Icges(5,4) * t99 + Icges(5,5) * t381;
t9 = (qJD(2) * t411 + t51) * t342 + (qJD(2) * t112 - t341 * t55 - t344 * t53 + (t115 * t341 + t117 * t344) * qJD(4)) * t345;
t599 = t9 * t253;
t597 = qJD(1) / 0.2e1;
t175 = Icges(5,3) * t345 + t342 * t417;
t229 = (-Icges(5,5) * t344 + Icges(5,6) * t341) * t345;
t124 = qJD(2) * t175 + qJD(4) * t229;
t179 = Icges(5,6) * t345 + t342 * t418;
t127 = qJD(2) * t179 + qJD(4) * t232;
t183 = Icges(5,5) * t345 + t342 * t422;
t130 = qJD(2) * t183 + qJD(4) * t237;
t409 = -t341 * t372 - t344 * t371;
t19 = (qJD(2) * t409 + t124) * t342 + (-qJD(2) * t370 - t127 * t344 - t130 * t341 + (-t341 * t371 + t344 * t372) * qJD(4)) * t345;
t468 = t345 * t490;
t70 = -t342 * t370 - t345 * t409;
t596 = t19 * t318 + t70 * t468;
t595 = t102 * rSges(5,1) + t101 * rSges(5,2);
t594 = rSges(3,1) * t345;
t592 = rSges(4,2) * t345;
t591 = rSges(4,3) * t342;
t589 = rSges(5,3) * t345;
t281 = rSges(3,1) * t342 + rSges(3,2) * t345;
t247 = t281 * t346;
t481 = t281 * t500;
t332 = t343 * rSges(3,3);
t203 = rSges(3,1) * t544 - rSges(3,2) * t548 + t332;
t525 = t203 + t288;
t90 = qJD(1) * t525 - t481;
t588 = t247 * t90;
t509 = rSges(3,2) * t549 + t346 * rSges(3,3);
t201 = rSges(3,1) * t546 - t509;
t480 = t281 * t498;
t89 = -t480 + (-t201 - t287) * qJD(1);
t583 = t343 * t89;
t581 = t346 * t89;
t41 = t110 * t342 - t345 * t412;
t580 = t41 * t167;
t579 = t42 * t166;
t578 = -rSges(4,3) - qJ(3);
t497 = qJD(3) * t345;
t437 = t243 * t500 + t248 * t498 - t497;
t34 = t121 * t253 + t123 * t252 + (t264 * t346 - t265 * t343) * qJD(2) + t437;
t558 = qJD(2) * t34;
t219 = t342 * t502 + t477;
t299 = pkin(2) * t479;
t120 = pkin(2) * t482 + qJ(3) * t219 - t299 + t474;
t263 = t288 * qJD(1);
t539 = -t120 - t263;
t538 = t121 + t264;
t537 = t123 - t265;
t524 = -t204 - t248;
t523 = t343 * t243 + t346 * t248;
t211 = qJD(2) * t283 - t497;
t284 = t591 - t592;
t522 = -t284 * qJD(2) - t211;
t245 = t279 * t346;
t521 = -qJD(1) * t245 + t343 * t497;
t491 = qJD(2) * qJD(3);
t470 = t345 * t491;
t520 = qJD(1) * t249 + t346 * t470;
t518 = -t248 - t264;
t513 = -t279 + t431;
t512 = -t283 - t284;
t511 = qJ(3) * t476 + t303;
t484 = t342 * t503;
t510 = rSges(3,2) * t484 + rSges(3,3) * t502;
t508 = t343 ^ 2 + t346 ^ 2;
t488 = t40 * t502;
t119 = -pkin(2) * t632 - qJ(3) * t484 + t511;
t487 = t346 * t119 + t343 * t120 + t243 * t502;
t240 = t279 * t343;
t486 = -t240 * t500 - t245 * t498 + t327;
t325 = pkin(5) * t502;
t485 = t325 + t511;
t473 = t121 * t494;
t472 = t123 * t494;
t471 = -pkin(1) - t594;
t466 = t502 / 0.2e1;
t465 = -t500 / 0.2e1;
t463 = t499 / 0.2e1;
t462 = -t498 / 0.2e1;
t461 = t498 / 0.2e1;
t455 = rSges(4,1) * t346 - rSges(4,3) * t549;
t454 = t513 * t346;
t449 = t120 * t500 + t342 * t491 + (t119 + t214) * t498;
t448 = qJD(1) * t240 + t345 * t496;
t251 = qJD(1) * (-pkin(1) * t503 + t325);
t445 = t343 * t470 + t251 + (t119 + t303) * qJD(1);
t444 = rSges(4,1) * t502 + rSges(4,2) * t632 + rSges(4,3) * t476;
t442 = t374 + t458;
t439 = qJD(4) * t463;
t438 = rSges(5,1) * t100 + rSges(5,2) * t99;
t434 = -rSges(3,2) * t342 + t594;
t430 = t343 * t36 + t346 * t35;
t429 = t343 * t35 - t346 * t36;
t428 = t343 * t38 + t346 * t37;
t427 = t343 * t37 - t346 * t38;
t426 = t343 * t42 + t346 * t41;
t425 = t343 * t41 - t346 * t42;
t424 = -t343 * t90 - t581;
t410 = t121 * t343 - t123 * t346;
t200 = t342 * t432 + t589;
t244 = (-rSges(5,1) * t344 + rSges(5,2) * t341) * t345;
t139 = qJD(2) * t200 + qJD(4) * t244;
t397 = -pkin(6) * t499 - t139 - t211;
t396 = -pkin(1) - t283;
t393 = qJD(2) * t454 + t303;
t242 = t281 * t343;
t241 = t431 * t343;
t83 = (t201 * t343 + t203 * t346) * qJD(2);
t379 = -t110 * t252 + t112 * t253 + t318 * t370;
t378 = (Icges(5,5) * t221 - Icges(5,6) * t222) * t252 - (Icges(5,5) * t223 + Icges(5,6) * t224) * t253 + t229 * t318;
t375 = t345 * t378;
t359 = t34 * t410 - (-t343 * t39 + t346 * t40) * t374;
t350 = (t370 * t346 + t412) * t252 - (t370 * t343 + t411) * t253 + (t175 + t409) * t318;
t349 = t350 * t345;
t347 = qJD(2) ^ 2;
t326 = pkin(3) * t502;
t261 = t434 * qJD(2);
t250 = t279 * t503;
t246 = t431 * t346;
t220 = t476 - t484;
t218 = t508 * t501;
t205 = rSges(4,2) * t546 + t455;
t173 = -pkin(6) * t632 + t326;
t172 = qJD(1) * t264 - t298;
t159 = t374 * t346;
t158 = t374 * t343;
t157 = t372 * t346;
t156 = t372 * t343;
t155 = t371 * t346;
t154 = t371 * t343;
t151 = rSges(5,1) * t223 + rSges(5,2) * t224;
t150 = rSges(5,1) * t221 - rSges(5,2) * t222;
t143 = -rSges(4,3) * t484 + t444;
t142 = qJD(2) * t241 + (t284 * t346 + t333) * qJD(1);
t141 = -qJD(2) * t242 + (t346 * t434 + t332) * qJD(1);
t140 = -rSges(3,1) * t632 - rSges(3,2) * t476 + t510;
t66 = qJD(1) * t625 - t249 - t626;
t65 = (t205 + t519) * qJD(1) + t393;
t64 = -t261 * t498 + (-t141 - t263 + t481) * qJD(1);
t63 = -t261 * t500 + t251 + (t140 - t480) * qJD(1);
t62 = (t204 * t346 - t205 * t343) * qJD(2) + t437;
t58 = -rSges(5,3) * t632 + t595;
t57 = rSges(5,3) * t381 + t438;
t33 = t522 * t498 + (-t142 + t539 + t626) * qJD(1) + t520;
t32 = qJD(1) * t143 + (qJD(1) * t454 + t343 * t522) * qJD(2) + t445;
t17 = (t142 * t343 + t143 * t346 + (-t205 * t346 + t343 * t524) * qJD(1)) * qJD(2) + t449;
t16 = -t101 * t371 - t102 * t372 + t124 * t544 + t127 * t221 + t130 * t222 + t370 * t632;
t15 = -t100 * t372 + t124 * t546 + t127 * t223 - t130 * t224 - t370 * t381 - t371 * t99;
t14 = t252 * t41 - t253 * t42 + t318 * t70;
t13 = -t347 * t320 - t139 * t253 - t166 * t374 - t318 * t57 + (-t211 * t346 - t472) * qJD(2) + (-t172 + (pkin(6) * qJD(2) - qJD(3)) * t549 + t539) * qJD(1) + t520;
t12 = -t347 * t319 + qJD(1) * t173 - t139 * t252 + t167 * t374 + t318 * t58 + (qJD(1) * t436 - t211 * t343 + t473) * qJD(2) + t445;
t7 = -t121 * t166 + t123 * t167 + t252 * t57 + t253 * t58 + (t172 * t343 + t173 * t346 + (-t265 * t346 + t343 * t518) * qJD(1)) * qJD(2) + t449;
t6 = t101 * t115 - t102 * t117 - t112 * t632 + t221 * t53 + t222 * t55 + t51 * t544;
t5 = t101 * t113 + t102 * t116 - t110 * t632 + t221 * t54 + t222 * t56 + t52 * t544;
t4 = -t100 * t117 + t112 * t381 + t115 * t99 + t223 * t53 - t224 * t55 + t51 * t546;
t3 = t100 * t116 + t110 * t381 + t113 * t99 + t223 * t54 - t224 * t56 + t52 * t546;
t2 = t16 * t318 + t166 * t36 + t167 * t35 + t252 * t5 - t253 * t6 + t468 * t60;
t1 = t15 * t318 + t166 * t38 + t167 * t37 + t252 * t3 - t253 * t4 + t468 * t61;
t18 = [t600 / 0.2e1 - t599 / 0.2e1 + t596 + t10 * t605 + t580 / 0.2e1 + t579 / 0.2e1 + t61 * t610 + t16 * t607 + t60 * t609 + (t15 + t10) * t606 + ((t641 * t343 + ((t679 + t695) * t346 + t649 + t683 + t694) * t346) * qJD(2) + t663) * t461 + (t686 * qJD(2) + t689 * t342 + t690 * t345) * qJD(1) + (-(qJD(1) * t265 - t39 + t624 + t628) * t40 + t13 * (-t319 + t337 + t433) + t39 * (t298 + t299 - t438) + t12 * (t443 + t538) + t40 * (t326 + t485 + t595) + (qJD(1) * t39 * t623 + t40 * t489 * t501 + t13 * pkin(3)) * t346 + (t13 * (t396 - t589) + t39 * (rSges(5,3) * t501 - qJ(3) * t499 - t327) + (t39 * (-pkin(3) - pkin(5)) + t623 * t40) * qJD(1)) * t343) * m(5) + (t32 * t625 + (t299 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t345 + t578 * t342) * t502 + (-t327 + (t345 * t578 - t593) * qJD(2) + (-rSges(4,1) - pkin(5)) * qJD(1)) * t343) * t65 + (t337 + t455 + (t396 + t592) * t343) * t33 + (-pkin(2) * t478 - t393 + t444 + t485 - t624 + t65 + (-t205 + (t396 - t591) * t343) * qJD(1)) * t66) * m(4) + (-(-qJD(1) * t201 - t266 - t480 - t89) * t90 + t64 * (t343 * t471 + t337 + t509) + t63 * t525 + t90 * (t325 + t510) + (t281 * t583 - t588) * qJD(2) + ((-pkin(1) - t434) * t581 + (t89 * (-rSges(3,3) - pkin(5)) + t90 * t471) * t343) * qJD(1)) * m(3) + (t651 + t652) * t500 / 0.2e1 + (-t493 + 0.2e1 * t492 + ((t346 * t642 - t641 + t647) * t346 + (t343 * t642 + t648 + t684) * t343) * qJD(2) + t655) * t465 + (t650 - t653 + t654) * t462 + ((t103 + t646) * t343 + (t645 + t680) * t346) * qJD(2) * t597; t427 * t610 + (qJD(1) * t426 + t343 * t8 - t346 * t9) * t603 + (((-t155 * t344 - t157 * t341 + t110) * t252 - (-t154 * t344 - t156 * t341 + t112) * t253 + (-t179 * t344 - t183 * t341 - t370) * t318 + t70 * qJD(4)) * t345 + (-qJD(4) * t426 + t350) * t342) * t604 + ((t155 * t223 - t157 * t224) * t252 - (t154 * t223 - t156 * t224) * t253 + (t179 * t223 - t183 * t224) * t318 + (t345 * t61 - t37 * t548) * qJD(4) + ((-qJD(4) * t38 + t379) * t342 + t349) * t343) * t605 + (qJD(1) * t428 + t3 * t343 - t346 * t4) * t606 + (qJD(1) * t430 + t343 * t5 - t346 * t6) * t607 + ((t155 * t221 + t157 * t222) * t252 - (t154 * t221 + t156 * t222) * t253 + (t179 * t221 + t183 * t222) * t318 + (t345 * t60 - t36 * t549) * qJD(4) + ((-qJD(4) * t35 + t379) * t342 + t349) * t346) * t608 + t429 * t609 - t14 * t494 / 0.2e1 + t425 * t439 - ((t661 * t342 + t656 * t345) * qJD(2) + (t657 * t342 + t658 * t345) * qJD(1)) * qJD(1) / 0.2e1 + (t653 * t346 + t652 * t343 + (t343 * t646 + t346 * t645) * qJD(1)) * t597 + ((-t500 * t630 + t631) * t343 + ((t343 * t629 + t634) * qJD(2) + t633) * t346) * t465 + ((-t498 * t629 - t631) * t346 + ((t346 * t630 + t634) * qJD(2) + t633) * t343) * t461 + (t346 * t10 + t343 * t11) * t495 / 0.2e1 + (-t40 * (t159 * t318 - t200 * t252 + t473 + t521) - ((-t508 * t558 - t488) * pkin(6) + t359 * qJD(4)) * t342 - t622 * qJD(2) * (-pkin(6) * t345 - t283) + t7 * t523 + (t12 * t442 + t40 * t397 + t7 * t537) * t343 + (t7 * t538 + (t40 * qJD(1) + t13) * t442) * t346 + (t158 * t318 + t200 * t253 + t397 * t346 - t374 * t503 + t250 - t448 + t472) * t39 + (-t158 * t252 - t159 * t253 - t486 + t487 + (t172 + t57 + (-t121 + t518) * qJD(1)) * t343 + (t537 * qJD(1) + t173 + t58) * t346) * t34) * m(5) + (-t65 * (-qJD(1) * t241 + t448) - t66 * (qJD(1) * t246 + t521) - t62 * t486 - ((t62 * t246 + t512 * t65) * t346 + (t62 * t241 + t512 * t66) * t343) * qJD(2) + t65 * t250 + t17 * t523 + t62 * t487 + (t33 * t513 + t65 * t522 + t17 * t204 + t62 * t143 + (-t62 * t205 + t513 * t66) * qJD(1)) * t346 + (t32 * t513 + t66 * t522 - t17 * t205 + t62 * t142 + (-t431 * t65 + t524 * t62) * qJD(1)) * t343) * m(4) + (0.2e1 * t83 * (t140 * t346 + t141 * t343 + (t201 * t346 - t203 * t343) * qJD(1)) + t424 * t261 + (-t63 * t343 - t64 * t346 + (-t346 * t90 + t583) * qJD(1)) * t281 - (t242 * t89 - t588) * qJD(1) - (t83 * (-t242 * t343 - t247 * t346) + t424 * t434) * qJD(2)) * m(3) + (t2 + t651 * qJD(1) + (t647 * t502 + (t648 * qJD(1) + t640 * t343 - t662) * t343) * t638) * t343 / 0.2e1 - (t1 + t650 * qJD(1) + ((t649 * qJD(1) + t662) * t346 + (qJD(1) * t681 - t640 * t346) * t343) * t638) * t346 / 0.2e1 + (t11 + t655 + t660) * t503 / 0.2e1 + (t10 + t654 + t659) * t466; -m(4) * (t218 * t62 + t219 * t66 + t220 * t65) - m(5) * (t218 * t34 + t219 * t40 + t220 * t39) + 0.2e1 * ((t498 * t65 + t500 * t66 - t17) * t612 + (t39 * t498 + t40 * t500 - t7) * t611) * t345 + 0.2e1 * ((qJD(2) * t62 + t32 * t343 + t33 * t346 + t502 * t66 - t503 * t65) * t612 + (t12 * t343 + t13 * t346 - t39 * t503 + t488 + t558) * t611) * t342; t2 * t544 / 0.2e1 + (t342 * t60 + t345 * t430) * t609 + ((-qJD(2) * t430 + t16) * t342 + (-qJD(1) * t429 + qJD(2) * t60 + t343 * t6 + t346 * t5) * t345) * t607 + t1 * t546 / 0.2e1 + (t342 * t61 + t345 * t428) * t610 + ((-qJD(2) * t428 + t15) * t342 + (-qJD(1) * t427 + qJD(2) * t61 + t3 * t346 + t343 * t4) * t345) * t606 + t14 * t463 + t342 * (t579 + t580 + t596 - t599 + t600) / 0.2e1 + (t342 * t70 + t345 * t426) * t439 + ((-qJD(2) * t426 + t19) * t342 + (-qJD(1) * t425 + qJD(2) * t70 + t343 * t9 + t346 * t8) * t345) * t603 + (t365 * t221 - t222 * t366 + t346 * t375) * t608 + (t223 * t365 + t224 * t366 + t343 * t375) * t605 + (t378 * t342 + (t366 * t341 - t344 * t365) * t345) * t604 + (t342 * t465 + t345 * t466) * t11 + (-t483 / 0.2e1 + t342 * t462) * t10 + ((qJD(2) * t359 + t12 * t121 - t13 * t123 - t39 * t57 + t40 * t58) * t342 + (t39 * (-qJD(2) * t123 + t139 * t343) + t40 * (qJD(2) * t121 - t139 * t346) - t7 * t410 + t34 * (-t121 * t502 - t123 * t503 - t343 * t58 + t346 * t57) - (qJD(1) * t622 - t12 * t346 + t13 * t343) * t374) * t345 - t39 * (-t151 * t318 - t244 * t253) - t40 * (t150 * t318 - t244 * t252) - t34 * (t150 * t253 + t151 * t252)) * m(5);];
tauc = t18(:);

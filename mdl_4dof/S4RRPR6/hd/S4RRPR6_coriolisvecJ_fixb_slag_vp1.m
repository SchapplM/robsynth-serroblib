% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:23
% EndTime: 2019-12-31 17:04:49
% DurationCPUTime: 21.73s
% Computational Cost: add. (10452->649), mult. (12802->844), div. (0->0), fcn. (10026->8), ass. (0->376)
t642 = Icges(3,3) + Icges(4,3);
t327 = qJ(2) + pkin(7);
t300 = sin(t327);
t301 = cos(t327);
t331 = sin(qJ(2));
t333 = cos(qJ(2));
t635 = Icges(3,5) * t333 + Icges(4,5) * t301 - Icges(3,6) * t331 - Icges(4,6) * t300;
t334 = cos(qJ(1));
t641 = t642 * t334;
t332 = sin(qJ(1));
t501 = t332 * t333;
t503 = t331 * t332;
t506 = t301 * t332;
t508 = t300 * t332;
t627 = -Icges(3,5) * t501 - Icges(4,5) * t506 + Icges(3,6) * t503 + Icges(4,6) * t508 + t641;
t636 = t642 * t332 + t635 * t334;
t637 = Icges(3,5) * t331 + Icges(4,5) * t300 + Icges(3,6) * t333 + Icges(4,6) * t301;
t535 = Icges(4,6) * t334;
t166 = Icges(4,4) * t506 - Icges(4,2) * t508 - t535;
t536 = Icges(3,6) * t334;
t191 = Icges(3,4) * t501 - Icges(3,2) * t503 - t536;
t640 = t166 * t300 + t191 * t331;
t271 = Icges(4,4) * t508;
t541 = Icges(4,5) * t334;
t168 = Icges(4,1) * t506 - t271 - t541;
t285 = Icges(3,4) * t503;
t542 = Icges(3,5) * t334;
t193 = Icges(3,1) * t501 - t285 - t542;
t622 = -t168 * t301 - t193 * t333 + t640;
t619 = -t622 * t332 + t627 * t334;
t544 = Icges(4,4) * t300;
t233 = Icges(4,1) * t301 - t544;
t169 = Icges(4,5) * t332 + t233 * t334;
t545 = Icges(3,4) * t331;
t263 = Icges(3,1) * t333 - t545;
t194 = Icges(3,5) * t332 + t263 * t334;
t638 = -t169 * t506 - t194 * t501;
t230 = Icges(4,2) * t301 + t544;
t289 = Icges(4,4) * t301;
t232 = Icges(4,1) * t300 + t289;
t260 = Icges(3,2) * t333 + t545;
t314 = Icges(3,4) * t333;
t262 = Icges(3,1) * t331 + t314;
t634 = t230 * t300 - t232 * t301 + t260 * t331 - t262 * t333;
t633 = t636 * t334 + t638;
t500 = t333 * t334;
t505 = t301 * t334;
t583 = -t169 * t505 - t194 * t500 - t636 * t332;
t632 = -t168 * t505 - t193 * t500 + t627 * t332;
t588 = t637 * t334;
t587 = t637 * t332;
t395 = -Icges(4,2) * t300 + t289;
t167 = Icges(4,6) * t332 + t334 * t395;
t396 = -Icges(3,2) * t331 + t314;
t192 = Icges(3,6) * t332 + t334 * t396;
t631 = t167 * t300 + t192 * t331;
t618 = -t167 * t508 - t192 * t503 - t633;
t502 = t331 * t334;
t507 = t300 * t334;
t617 = -t166 * t507 - t191 * t502 - t632;
t616 = -t167 * t507 - t192 * t502 - t583;
t599 = t166 * t301 + t168 * t300 + t191 * t333 + t193 * t331;
t598 = t167 * t301 + t169 * t300 + t192 * t333 + t194 * t331;
t629 = -t634 * t332 - t588;
t628 = -t634 * t334 + t587;
t625 = t637 * qJD(2);
t624 = t169 * t301 + t194 * t333 - t631;
t207 = t395 * qJD(2);
t208 = t233 * qJD(2);
t241 = t396 * qJD(2);
t242 = t263 * qJD(2);
t623 = -t207 * t300 + t208 * t301 - t241 * t331 + t242 * t333 + (-t230 * t301 - t232 * t300 - t260 * t333 - t262 * t331) * qJD(2) + t637 * qJD(1);
t621 = t636 * qJD(1);
t620 = t634 * qJD(1) + t635 * qJD(2);
t615 = t628 * qJD(1);
t370 = qJD(2) * t260;
t119 = -t334 * t370 + (-t332 * t396 + t536) * qJD(1);
t372 = qJD(2) * t262;
t121 = -t334 * t372 + (-t263 * t332 + t542) * qJD(1);
t369 = qJD(2) * t230;
t96 = -t334 * t369 + (-t332 * t395 + t535) * qJD(1);
t371 = qJD(2) * t232;
t98 = -t334 * t371 + (-t233 * t332 + t541) * qJD(1);
t614 = -qJD(2) * t598 - t119 * t331 + t121 * t333 - t300 * t96 + t301 * t98 + t621;
t120 = qJD(1) * t192 - t332 * t370;
t122 = qJD(1) * t194 - t332 * t372;
t97 = qJD(1) * t167 - t332 * t369;
t99 = qJD(1) * t169 - t332 * t371;
t613 = qJD(1) * t627 + qJD(2) * t599 + t120 * t331 - t122 * t333 + t300 * t97 - t301 * t99;
t612 = (t332 * t616 - t334 * t617) * qJD(2);
t611 = (t618 * t332 - t334 * t619) * qJD(2);
t610 = t629 * qJD(1);
t609 = qJD(1) * t622 - t332 * t625 + t621;
t608 = -t625 * t334 + (-t332 * t635 - t624 + t641) * qJD(1);
t607 = 0.2e1 * qJD(2);
t606 = t610 + t611;
t605 = t612 + t615;
t604 = t620 * t332 + t623 * t334;
t603 = t623 * t332 - t620 * t334;
t602 = t622 * qJD(2) - t120 * t333 - t122 * t331 - t300 * t99 - t301 * t97;
t601 = qJD(2) * t624 + t119 * t333 + t121 * t331 + t300 * t98 + t301 * t96;
t600 = rSges(3,2) * t331;
t272 = rSges(4,2) * t508;
t174 = rSges(4,1) * t506 - t334 * rSges(4,3) - t272;
t316 = t332 * rSges(4,3);
t175 = rSges(4,1) * t505 - rSges(4,2) * t507 + t316;
t388 = t174 * t332 + t175 * t334;
t323 = t334 * pkin(5);
t274 = pkin(1) * t332 - t323;
t330 = -qJ(3) - pkin(5);
t294 = t334 * t330;
t322 = t333 * pkin(2);
t295 = t322 + pkin(1);
t464 = t332 * t295 + t294;
t162 = t274 - t464;
t321 = t332 * pkin(5);
t275 = t334 * pkin(1) + t321;
t277 = t334 * t295;
t410 = -t330 * t332 + t277;
t163 = t410 - t275;
t486 = -t332 * t162 + t334 * t163;
t366 = t388 + t486;
t451 = qJD(2) * t334;
t452 = qJD(2) * t332;
t487 = -t162 * t452 + t163 * t451;
t53 = qJD(2) * t388 + t487;
t597 = qJD(2) * t366 + t53;
t302 = qJ(4) + t327;
t291 = cos(t302);
t510 = t291 * t332;
t290 = sin(t302);
t512 = t290 * t332;
t534 = Icges(5,6) * t334;
t152 = Icges(5,4) * t510 - Icges(5,2) * t512 - t534;
t281 = Icges(5,4) * t291;
t214 = Icges(5,1) * t290 + t281;
t596 = -t214 * t332 - t152;
t394 = -Icges(5,2) * t290 + t281;
t153 = Icges(5,6) * t332 + t334 * t394;
t595 = -t214 * t334 - t153;
t543 = Icges(5,4) * t290;
t215 = Icges(5,1) * t291 - t543;
t155 = Icges(5,5) * t332 + t215 * t334;
t212 = Icges(5,2) * t291 + t543;
t594 = -t212 * t334 + t155;
t593 = -t212 + t215;
t592 = t214 + t394;
t360 = t191 * t334 - t192 * t332;
t361 = t166 * t334 - t167 * t332;
t571 = t332 * (-t230 * t334 + t169) - t334 * (-Icges(4,2) * t506 + t168 - t271);
t572 = t332 * (-t260 * t334 + t194) - t334 * (-Icges(3,2) * t501 + t193 - t285);
t591 = -t300 * t571 + t361 * t301 - t331 * t572 + t360 * t333;
t466 = t262 + t396;
t467 = -t260 + t263;
t471 = t232 + t395;
t472 = -t230 + t233;
t590 = (-t300 * t471 + t301 * t472 - t331 * t466 + t333 * t467) * qJD(1);
t589 = t635 * qJD(1);
t156 = qJD(1) * t162;
t249 = qJD(1) * t274;
t586 = t156 - t249;
t303 = qJD(3) * t332;
t560 = pkin(3) * t300;
t561 = pkin(2) * t331;
t248 = -t560 - t561;
t234 = t248 * qJD(2);
t499 = t334 * t234;
t585 = t303 + t499;
t584 = t627 + t631;
t211 = Icges(5,5) * t291 - Icges(5,6) * t290;
t326 = qJD(2) + qJD(4);
t210 = Icges(5,5) * t290 + Icges(5,6) * t291;
t518 = t210 * t334;
t526 = t153 * t290;
t531 = Icges(5,3) * t334;
t579 = -t326 * t518 + (-t155 * t291 - t211 * t332 + t526 + t531) * qJD(1);
t253 = Icges(5,4) * t512;
t540 = Icges(5,5) * t334;
t154 = Icges(5,1) * t510 - t253 - t540;
t392 = t152 * t290 - t154 * t291;
t151 = Icges(5,3) * t332 + t211 * t334;
t459 = qJD(1) * t151;
t519 = t210 * t332;
t578 = qJD(1) * t392 - t326 * t519 + t459;
t385 = t212 * t290 - t214 * t291;
t575 = qJD(1) * t385 + t211 * t326;
t256 = t326 * t332;
t257 = t326 * t334;
t570 = qJD(1) * t592 + t256 * t594 - t257 * (-Icges(5,2) * t510 + t154 - t253);
t407 = qJD(1) * t326;
t238 = t332 * t407;
t569 = t238 / 0.2e1;
t239 = t334 * t407;
t568 = t239 / 0.2e1;
t567 = -t256 / 0.2e1;
t566 = t256 / 0.2e1;
t565 = -t257 / 0.2e1;
t564 = t257 / 0.2e1;
t563 = t332 / 0.2e1;
t562 = -t334 / 0.2e1;
t559 = pkin(3) * t301;
t558 = -qJD(1) / 0.2e1;
t557 = qJD(1) / 0.2e1;
t556 = pkin(1) - t295;
t555 = rSges(3,1) * t333;
t554 = rSges(4,1) * t301;
t553 = rSges(5,1) * t291;
t552 = rSges(4,2) * t301;
t551 = rSges(5,2) * t291;
t317 = t332 * rSges(3,3);
t315 = t332 * rSges(5,3);
t216 = rSges(5,1) * t290 + t551;
t355 = -t257 * t216 + t585;
t481 = t162 - t274;
t243 = t295 + t559;
t325 = -pkin(6) + t330;
t473 = -t332 * t243 - t334 * t325;
t125 = t464 + t473;
t254 = rSges(5,2) * t512;
t157 = rSges(5,1) * t510 - t334 * rSges(5,3) - t254;
t495 = t125 - t157;
t47 = (t481 + t495) * qJD(1) + t355;
t550 = t332 * t47;
t235 = rSges(4,1) * t300 + t552;
t426 = -t235 - t561;
t402 = t334 * t426;
t381 = qJD(2) * t402;
t364 = t303 + t381;
t60 = (-t174 + t481) * qJD(1) + t364;
t549 = t60 * t235;
t461 = rSges(3,2) * t503 + t334 * rSges(3,3);
t203 = rSges(3,1) * t501 - t461;
t265 = rSges(3,1) * t331 + rSges(3,2) * t333;
t437 = t265 * t451;
t106 = -t437 + (-t203 - t274) * qJD(1);
t530 = t106 * t332;
t529 = t106 * t334;
t438 = t265 * t452;
t204 = rSges(3,1) * t500 - rSges(3,2) * t502 + t317;
t475 = t204 + t275;
t107 = qJD(1) * t475 - t438;
t227 = t265 * t334;
t528 = t107 * t227;
t150 = Icges(5,5) * t510 - Icges(5,6) * t512 - t531;
t527 = t150 * t334;
t515 = t234 * t332;
t511 = t290 * t334;
t509 = t291 * t334;
t504 = t325 * t332;
t73 = -t332 * t385 - t518;
t498 = t73 * qJD(1);
t454 = qJD(1) * t332;
t282 = t330 * t454;
t279 = t452 * t561;
t304 = qJD(3) * t334;
t462 = t279 + t304;
t441 = t282 + t462;
t116 = (-t334 * t556 - t321) * qJD(1) - t441;
t246 = t275 * qJD(1);
t496 = -t116 - t246;
t219 = t334 * t243;
t460 = -t325 + t330;
t126 = t332 * t460 + t219 - t277;
t494 = -t126 - t163;
t493 = -t332 * t150 - t154 * t509;
t492 = t332 * t151 + t155 * t509;
t158 = rSges(5,1) * t509 - rSges(5,2) * t511 + t315;
t489 = t332 * t157 + t334 * t158;
t480 = -t163 - t175;
t470 = t243 - t295;
t453 = qJD(1) * t334;
t469 = rSges(5,3) * t453 + qJD(1) * t254;
t468 = rSges(4,3) * t453 + qJD(1) * t272;
t449 = qJD(1) * qJD(3);
t465 = qJD(1) * t279 + t334 * t449;
t463 = rSges(3,3) * t453 + t454 * t600;
t448 = pkin(2) * t502;
t92 = -t257 * t551 + (-t257 * t290 - t291 * t454) * rSges(5,1) + t469;
t187 = t216 * t332;
t217 = -rSges(5,2) * t290 + t553;
t93 = -t326 * t187 + (t217 * t334 + t315) * qJD(1);
t447 = t157 * t453 + t332 * t93 + t334 * t92;
t446 = qJD(2) * t322;
t299 = pkin(5) * t453;
t436 = t331 * t451;
t409 = pkin(2) * t436;
t115 = -t409 - t299 + t303 + (t332 * t556 - t294) * qJD(1);
t445 = t116 * t452 + (t115 - t156) * t451;
t444 = t334 * t115 + t332 * t116 - t162 * t453;
t237 = qJD(1) * (-pkin(1) * t454 + t299);
t443 = qJD(1) * t115 + t332 * t449 + t237;
t442 = -t158 + t494;
t439 = t235 * t452;
t435 = t333 * t451;
t433 = -pkin(1) - t555;
t432 = t454 / 0.2e1;
t431 = t453 / 0.2e1;
t430 = -t452 / 0.2e1;
t427 = t451 / 0.2e1;
t236 = -rSges(4,2) * t300 + t554;
t425 = -t236 - t322;
t424 = t248 + t561;
t422 = qJD(1) * t155 + t326 * t596;
t421 = (-t215 * t332 + t540) * qJD(1) + t595 * t326;
t420 = qJD(1) * t153 + t154 * t326 - t212 * t256;
t419 = (-t332 * t394 + t534) * qJD(1) + t594 * t326;
t127 = t155 * t510;
t418 = t151 * t334 - t127;
t416 = -t150 + t526;
t413 = t592 * t326;
t412 = t593 * t326;
t408 = t452 * t560;
t406 = (-t332 ^ 2 - t334 ^ 2) * t561;
t403 = -t322 - t559;
t400 = t555 - t600;
t393 = -t107 * t332 - t529;
t77 = t152 * t291 + t154 * t290;
t382 = -t216 + t248;
t209 = t236 * qJD(2);
t335 = qJD(2) ^ 2;
t380 = -qJD(2) * t209 - t322 * t335;
t379 = t403 * t335;
t226 = t265 * t332;
t201 = t235 * t332;
t188 = t216 * t334;
t373 = t392 * t332;
t102 = (t203 * t332 + t204 * t334) * qJD(2);
t176 = t217 * t326;
t363 = -qJD(2) * t559 - t176 - t446;
t362 = qJD(1) * t211 - t256 * t518 + t257 * t519;
t347 = -t290 * t419 + t291 * t421 + t459;
t10 = t332 * t579 + t347 * t334;
t348 = qJD(1) * t150 - t290 * t420 + t291 * t422;
t11 = t348 * t332 - t334 * t578;
t12 = t347 * t332 - t334 * t579;
t54 = -t373 - t527;
t55 = -t153 * t512 - t418;
t20 = t256 * t55 - t257 * t54 + t498;
t56 = -t152 * t511 - t493;
t57 = -t153 * t511 + t492;
t74 = -t334 * t385 + t519;
t68 = t74 * qJD(1);
t21 = t256 * t57 - t257 * t56 + t68;
t346 = qJD(1) * t210 - t290 * t413 + t291 * t412;
t30 = t332 * t575 + t346 * t334;
t31 = t346 * t332 - t334 * t575;
t353 = qJD(1) * t593 + t256 * t595 - t257 * t596;
t337 = -t290 * t570 + t353 * t291;
t37 = t290 * t422 + t291 * t420;
t38 = t290 * t421 + t291 * t419;
t78 = t153 * t291 + t155 * t290;
t9 = t332 * t578 + t348 * t334;
t359 = (qJD(1) * t30 + t10 * t256 + t238 * t56 + t239 * t57 - t257 * t9) * t563 + (t353 * t290 + t291 * t570) * t558 + t20 * t432 + t21 * t431 + (qJD(1) * t31 - t11 * t257 + t12 * t256 + t238 * t54 + t239 * t55) * t562 + (t10 * t332 - t334 * t9 + (t332 * t56 + t334 * t57) * qJD(1)) * t566 + (t332 * t55 - t334 * t54) * t569 + (t332 * t57 - t334 * t56) * t568 + (-t11 * t334 + t12 * t332 + (t332 * t54 + t334 * t55) * qJD(1)) * t565 + (t332 * t38 - t334 * t37 + (t332 * t77 + t334 * t78) * qJD(1)) * t557 + (t332 * t362 + t334 * t337) * t567 + (t332 * t337 - t334 * t362) * t564;
t34 = t157 * t256 + t158 * t257 + (-t125 * t332 + t126 * t334) * qJD(2) + t487;
t48 = -t408 - t216 * t256 + (t275 - t442) * qJD(1) - t462;
t356 = t34 * (-t256 * t187 - t188 * t257) + t48 * (-qJD(1) * t188 - t217 * t256);
t244 = t400 * qJD(2);
t202 = t235 * t334;
t180 = t424 * t334;
t179 = t424 * t332;
t149 = t257 * t217;
t124 = -qJD(2) * t226 + (t334 * t400 + t317) * qJD(1);
t123 = -rSges(3,2) * t435 + (-t333 * t454 - t436) * rSges(3,1) + t463;
t101 = -qJD(2) * t201 + (t236 * t334 + t316) * qJD(1);
t100 = -t451 * t552 + (-t300 * t451 - t301 * t454) * rSges(4,1) + t468;
t76 = t515 + t279 + t282 + (t334 * t470 - t504) * qJD(1);
t75 = t409 + t499 + (-t332 * t470 + t334 * t460) * qJD(1);
t67 = -t244 * t451 + (-t124 - t246 + t438) * qJD(1);
t66 = -t244 * t452 + t237 + (t123 - t437) * qJD(1);
t61 = -t439 + (t275 - t480) * qJD(1) - t462;
t42 = t380 * t334 + (-t101 + t439 + t496) * qJD(1) + t465;
t41 = t380 * t332 + (t100 + t381) * qJD(1) + t443;
t23 = -t176 * t257 + t216 * t238 + t334 * t379 + (-t76 - t93 + t408 + t496) * qJD(1) + t465;
t22 = -t176 * t256 - t216 * t239 + t332 * t379 + (t75 + t92 + t499) * qJD(1) + t443;
t5 = t157 * t239 - t158 * t238 + t256 * t93 + t257 * t92 + (t332 * t76 + t334 * t75 + (-t125 * t334 + t332 * t494) * qJD(1)) * qJD(2) + t445;
t1 = [(t68 + (t55 + (t152 * t334 + t153 * t332) * t290 + t418 + t493) * t257 + (-t154 * t510 + t527 + t54 + (t152 * t332 - t153 * t334) * t290 + t492) * t256) * t564 + (t77 + t73) * t569 + (t78 + t74) * t568 + (-t498 + (t57 - t373 - t492) * t257 + (t332 * t416 - t127 + t56) * t256 + ((t151 + t392) * t256 + t416 * t257) * t334 + t20) * t567 + (t38 + t30) * t566 + ((((t636 + t640) * t334 + t618 + t632 + t638) * t334 - t583 * t332) * qJD(2) + t615) * t427 + (-t634 * qJD(2) + t207 * t301 + t208 * t300 + t241 * t333 + t242 * t331 + t290 * t412 + t291 * t413) * qJD(1) + (t23 * (-t157 + t473) + t47 * (t304 - t515) + t22 * (t158 + t219 - t504) + t48 * (t469 + t585) + (-t188 * t48 + t216 * t550) * t326 + ((t47 * (-t217 - t243) - t48 * t325) * t334 + (t47 * (-rSges(5,3) + t325) + t48 * (-t243 - t553)) * t332) * qJD(1) - (qJD(1) * t495 + t355 - t47 + t586) * t48) * m(5) + (t42 * (-t174 - t464) + t60 * t441 + t41 * (t175 + t410) + t61 * (t303 + t468) + (t332 * t549 + t402 * t61) * qJD(2) + ((-t60 * rSges(4,3) + t61 * (-t295 - t554)) * t332 + (t60 * (-t236 - t295) - t61 * t330) * t334) * qJD(1) - (-qJD(1) * t174 + t364 + t586 - t60) * t61) * m(4) + (t67 * (t433 * t332 + t323 + t461) + t66 * t475 + t107 * (t299 + t463) + (t265 * t530 - t528) * qJD(2) + ((-pkin(1) - t400) * t529 + (t106 * (-rSges(3,3) - pkin(5)) + t107 * t433) * t332) * qJD(1) - (-qJD(1) * t203 - t106 - t249 - t437) * t107) * m(3) + (t37 + t31 + t21) * t565 + (((t334 * t584 + t583 + t616) * t334 + (t332 * t584 + t617 + t633) * t332) * qJD(2) + t606 - t610) * t430 + (t601 + t604) * t452 / 0.2e1 - (-t602 + t603 + t605) * t451 / 0.2e1 + ((t599 + t629) * t332 + (t598 + t628) * t334) * qJD(2) * t557; t359 + ((t361 * t300 + t301 * t571 + t360 * t331 + t333 * t572) * qJD(2) + (t300 * t472 + t301 * t471 + t331 * t467 + t333 * t466) * qJD(1)) * t558 + (t602 * t334 + t601 * t332 + (t599 * t332 + t598 * t334) * qJD(1)) * t557 + ((-t452 * t588 + t589) * t332 + ((t587 * t332 + t591) * qJD(2) + t590) * t334) * t430 + ((-t451 * t587 - t589) * t334 + ((t334 * t588 + t591) * qJD(2) + t590) * t332) * t427 + (t47 * t149 - (t47 * (-t179 + t187) + t48 * (t180 - t448)) * qJD(1) - (t34 * t406 + (t34 * t180 + t403 * t47) * t334 + (t34 * t179 + t403 * t48) * t332) * qJD(2) - t356 + t5 * (t486 + t489) + t34 * (t444 + t447) + (t23 * t382 + t47 * t363 + t5 * t126 + t34 * t75 + (-t34 * t125 + t382 * t48) * qJD(1)) * t334 + (t22 * t382 + t48 * t363 - t5 * t125 + t34 * t76 + (t47 * (t216 + t560) + t34 * t442) * qJD(1)) * t332) * m(5) + (-(t60 * t201 + t61 * (-t202 - t448)) * qJD(1) - (t53 * t406 + (-t53 * t202 + t425 * t60) * t334 + (-t53 * t201 + t425 * t61) * t332) * qJD(2) + t42 * t402 - t60 * pkin(2) * t435 + (t100 * t451 + t101 * t452 + t445) * t366 + t53 * t444 + (-t60 * t209 + t53 * t100 + (t174 * t597 + t426 * t61) * qJD(1)) * t334 + (t41 * t426 + t61 * (-t209 - t446) + t53 * t101 + (t480 * t597 + t549) * qJD(1)) * t332) * m(4) + (0.2e1 * t102 * (t123 * t334 + t124 * t332 + (t203 * t334 - t204 * t332) * qJD(1)) + t393 * t244 + (-t66 * t332 - t67 * t334 + (-t107 * t334 + t530) * qJD(1)) * t265 - (t106 * t226 - t528) * qJD(1) - (t102 * (-t226 * t332 - t227 * t334) + t393 * t400) * qJD(2)) * m(3) + (t604 * qJD(1) + ((t616 * qJD(1) + t613 * t334) * t334 + (t608 * t332 + t617 * qJD(1) + (-t609 + t614) * t334) * t332) * t607) * t563 + (t603 * qJD(1) + ((t618 * qJD(1) + t609 * t334) * t334 + (t614 * t332 + t619 * qJD(1) + (-t608 + t613) * t334) * t332) * t607) * t562 + (t606 + t611) * t432 + (t605 + t612) * t431; 0.2e1 * (t22 * t562 + t23 * t563) * m(5) + 0.2e1 * (t41 * t562 + t42 * t563) * m(4); t359 + (-t47 * (qJD(1) * t187 - t149) - t356 + t5 * t489 + t34 * (-t158 * t454 + t447) + (-t332 * t48 - t334 * t47) * t176 + (-t22 * t332 - t23 * t334 + (-t334 * t48 + t550) * qJD(1)) * t216) * m(5);];
tauc = t1(:);

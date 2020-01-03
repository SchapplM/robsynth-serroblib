% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:15
% EndTime: 2020-01-03 11:39:15
% DurationCPUTime: 34.22s
% Computational Cost: add. (16646->727), mult. (13380->913), div. (0->0), fcn. (10340->10), ass. (0->402)
t645 = Icges(4,3) + Icges(5,3);
t321 = qJ(3) + pkin(9);
t309 = sin(t321);
t311 = cos(t321);
t233 = Icges(5,5) * t311 - Icges(5,6) * t309;
t324 = sin(qJ(3));
t326 = cos(qJ(3));
t271 = Icges(4,5) * t326 - Icges(4,6) * t324;
t638 = t233 + t271;
t322 = qJ(1) + pkin(8);
t312 = cos(t322);
t499 = t311 * t312;
t503 = t309 * t312;
t310 = sin(t322);
t525 = Icges(5,6) * t310;
t160 = Icges(5,4) * t499 - Icges(5,2) * t503 + t525;
t496 = t312 * t326;
t497 = t312 * t324;
t526 = Icges(4,6) * t310;
t178 = Icges(4,4) * t496 - Icges(4,2) * t497 + t526;
t644 = -t160 * t309 - t178 * t324;
t643 = t645 * t310;
t614 = t638 * t310 - t645 * t312;
t613 = Icges(4,5) * t496 + Icges(5,5) * t499 - Icges(4,6) * t497 - Icges(5,6) * t503 + t643;
t535 = Icges(4,4) * t324;
t272 = Icges(4,2) * t326 + t535;
t495 = t324 * t272;
t314 = Icges(4,4) * t326;
t586 = Icges(4,1) * t324 + t314;
t370 = t326 * t586 - t495;
t534 = Icges(5,4) * t309;
t234 = Icges(5,2) * t311 + t534;
t298 = Icges(5,4) * t311;
t587 = Icges(5,1) * t309 + t298;
t372 = -t234 * t309 + t311 * t587;
t642 = t370 + t372;
t261 = Icges(5,4) * t503;
t531 = Icges(5,5) * t310;
t162 = Icges(5,1) * t499 - t261 + t531;
t282 = Icges(4,4) * t497;
t532 = Icges(4,5) * t310;
t180 = Icges(4,1) * t496 - t282 + t532;
t641 = -t162 * t311 - t180 * t326 - t644;
t237 = Icges(5,1) * t311 - t534;
t161 = -Icges(5,5) * t312 + t237 * t310;
t275 = Icges(4,1) * t326 - t535;
t179 = -Icges(4,5) * t312 + t275 * t310;
t500 = t310 * t326;
t502 = t310 * t311;
t616 = -t161 * t502 - t179 * t500;
t640 = t162 * t502 + t180 * t500;
t384 = -Icges(5,2) * t309 + t298;
t159 = -Icges(5,6) * t312 + t310 * t384;
t385 = -Icges(4,2) * t324 + t314;
t177 = -Icges(4,6) * t312 + t310 * t385;
t615 = t159 * t503 + t177 * t497;
t232 = Icges(5,5) * t309 + Icges(5,6) * t311;
t187 = t232 * t312;
t270 = Icges(4,5) * t324 + Icges(4,6) * t326;
t202 = t270 * t312;
t639 = t187 + t202;
t509 = t270 * t310;
t510 = t232 * t310;
t637 = -t509 - t510;
t501 = t310 * t324;
t504 = t309 * t310;
t621 = -t159 * t504 - t177 * t501 - t614 * t312 - t616;
t620 = -t160 * t504 - t178 * t501 - t613 * t312 + t640;
t619 = -t161 * t499 - t179 * t496 - t614 * t310 + t615;
t618 = -t613 * t310 + t641 * t312;
t636 = t642 * t310 - t639;
t635 = -t234 * t503 - t312 * t495 + t496 * t586 + t499 * t587 - t637;
t610 = -t161 * t311 - t179 * t326;
t609 = t159 * t309 + t177 * t324;
t634 = t642 * qJD(1) - t638 * qJD(3);
t633 = t609 + t610;
t217 = t384 * qJD(3);
t218 = t237 * qJD(3);
t254 = t385 * qJD(3);
t255 = t275 * qJD(3);
t632 = t217 * t309 - t218 * t311 + t254 * t324 - t255 * t326 + (t234 * t311 + t272 * t326 + t309 * t587 + t324 * t586) * qJD(3) + (-t232 - t270) * qJD(1);
t631 = t636 * qJD(1);
t630 = (t618 * t310 - t619 * t312) * qJD(3);
t629 = (t620 * t310 - t621 * t312) * qJD(3);
t628 = t635 * qJD(1);
t627 = t629 + t631;
t626 = -t628 + t630;
t188 = t234 * t310;
t101 = -qJD(3) * t188 + (t312 * t384 + t525) * qJD(1);
t190 = t587 * t310;
t103 = -qJD(3) * t190 + (t237 * t312 + t531) * qJD(1);
t203 = t272 * t310;
t115 = -qJD(3) * t203 + (t312 * t385 + t526) * qJD(1);
t205 = t586 * t310;
t117 = -qJD(3) * t205 + (t275 * t312 + t532) * qJD(1);
t625 = t633 * qJD(3) - t101 * t311 - t103 * t309 - t115 * t326 - t117 * t324;
t447 = qJD(3) * t312;
t100 = qJD(1) * t159 + t234 * t447;
t191 = t587 * t312;
t102 = qJD(1) * t161 + qJD(3) * t191;
t114 = qJD(1) * t177 + t272 * t447;
t206 = t586 * t312;
t116 = qJD(1) * t179 + qJD(3) * t206;
t624 = -qJD(3) * t641 - t100 * t311 - t102 * t309 - t114 * t326 - t116 * t324;
t623 = t634 * t310 + t632 * t312;
t622 = -t632 * t310 + t634 * t312;
t617 = t159 * t311 + t161 * t309 + t177 * t326 + t179 * t324;
t611 = t160 * t311 + t162 * t309 + t178 * t326 + t180 * t324;
t608 = t614 * qJD(1);
t313 = qJ(5) + t321;
t302 = sin(t313);
t303 = cos(t313);
t507 = t302 * t312;
t248 = Icges(6,4) * t507;
t505 = t303 * t312;
t530 = Icges(6,5) * t310;
t153 = Icges(6,1) * t505 - t248 + t530;
t320 = qJD(3) + qJD(5);
t230 = t310 * t320;
t498 = t312 * t320;
t293 = Icges(6,4) * t303;
t383 = -Icges(6,2) * t302 + t293;
t588 = Icges(6,1) * t302 + t293;
t596 = t588 + t383;
t533 = Icges(6,4) * t302;
t225 = Icges(6,1) * t303 - t533;
t152 = -Icges(6,5) * t312 + t225 * t310;
t222 = Icges(6,2) * t303 + t533;
t603 = -t222 * t310 + t152;
t335 = qJD(1) * t596 + t230 * (-Icges(6,2) * t505 + t153 - t248) - t498 * t603;
t602 = t222 - t225;
t524 = Icges(6,6) * t310;
t151 = Icges(6,4) * t505 - Icges(6,2) * t507 + t524;
t604 = t588 * t312 + t151;
t150 = -Icges(6,6) * t312 + t310 * t383;
t605 = t588 * t310 + t150;
t570 = qJD(1) * t602 + t230 * t604 - t498 * t605;
t607 = t335 * t302 + t303 * t570;
t349 = t310 * (-Icges(5,2) * t499 + t162 - t261) - t312 * (t161 - t188);
t575 = t310 * (t160 + t191) - t312 * (t159 + t190);
t606 = -t349 * t309 - t311 * t575;
t601 = t637 * qJD(3) + (t638 * t312 + t633 + t643) * qJD(1);
t600 = -qJD(1) * t641 + t639 * qJD(3) + t608;
t598 = qJD(1) * t613 - qJD(3) * t611 + t100 * t309 - t102 * t311 + t114 * t324 - t116 * t326;
t597 = -qJD(3) * t617 - t101 * t309 + t103 * t311 - t115 * t324 + t117 * t326 + t608;
t463 = t586 + t385;
t464 = t272 - t275;
t595 = (t324 * t463 + t326 * t464) * qJD(1);
t470 = t587 + t384;
t471 = t234 - t237;
t594 = (t309 * t470 + t311 * t471) * qJD(1);
t593 = 0.2e1 * qJD(3);
t221 = Icges(6,5) * t303 - Icges(6,6) * t302;
t148 = -Icges(6,3) * t312 + t221 * t310;
t519 = t151 * t302;
t382 = -t153 * t303 + t519;
t368 = -t148 + t382;
t592 = t498 * t368;
t317 = t326 * pkin(3);
t304 = t317 + pkin(2);
t266 = t312 * t304;
t323 = -qJ(4) - pkin(6);
t560 = pkin(2) * t312;
t147 = t560 - t266 + (pkin(6) + t323) * t310;
t142 = qJD(1) * t147;
t244 = pkin(6) * t310 + t560;
t229 = qJD(1) * t244;
t589 = t229 - t142;
t393 = -t310 * rSges(3,1) - t312 * rSges(3,2);
t553 = rSges(5,1) * t309;
t240 = rSges(5,2) * t311 + t553;
t559 = pkin(3) * t324;
t416 = t240 + t559;
t390 = qJD(3) * t416;
t547 = rSges(6,2) * t303;
t552 = rSges(6,1) * t302;
t227 = t547 + t552;
t181 = t227 * t310;
t182 = t227 * t312;
t548 = rSges(6,2) * t302;
t551 = rSges(6,1) * t303;
t228 = -t548 + t551;
t288 = t312 * t323;
t466 = t310 * t304 + t288;
t557 = pkin(4) * t311;
t256 = t304 + t557;
t319 = -pkin(7) + t323;
t474 = t310 * t256 + t312 * t319;
t122 = -t466 + t474;
t213 = t312 * t256;
t460 = t319 - t323;
t123 = t310 * t460 - t213 + t266;
t506 = t303 * t310;
t250 = rSges(6,1) * t506;
t508 = t302 * t310;
t544 = rSges(6,3) * t312;
t154 = -rSges(6,2) * t508 + t250 - t544;
t251 = rSges(6,2) * t507;
t434 = rSges(6,1) * t505;
t155 = rSges(6,3) * t310 - t251 + t434;
t301 = t310 * pkin(2);
t243 = -pkin(6) * t312 + t301;
t146 = -t243 + t466;
t449 = qJD(3) * t310;
t442 = t146 * t449 + qJD(2);
t33 = t154 * t230 + t155 * t498 + (t122 * t310 + (-t123 - t147) * t312) * qJD(3) + t442;
t297 = qJD(4) * t310;
t558 = pkin(4) * t309;
t396 = t558 + t559;
t316 = sin(qJ(1)) * pkin(1);
t417 = t243 + t316;
t402 = t146 + t417;
t47 = t227 * t498 - t297 + t396 * t447 + (t122 + t154 + t402) * qJD(1);
t584 = -t47 * (-qJD(1) * t181 + t228 * t498) - t33 * (-t181 * t230 - t182 * t498);
t583 = (-t310 ^ 2 - t312 ^ 2) * t559;
t220 = Icges(6,5) * t302 + Icges(6,6) * t303;
t406 = t602 * t320;
t407 = t596 * t320;
t573 = -qJD(1) * t220 + t302 * t407 + t303 * t406;
t521 = Icges(6,3) * t310;
t149 = Icges(6,5) * t505 - Icges(6,6) * t507 + t521;
t411 = -qJD(1) * t150 + t153 * t320 - t222 * t498;
t413 = qJD(1) * t152 + t320 * t604;
t572 = -qJD(1) * t149 + t302 * t411 + t303 * t413;
t412 = (t312 * t383 + t524) * qJD(1) + t603 * t320;
t414 = -(t225 * t312 + t530) * qJD(1) + t605 * t320;
t459 = qJD(1) * t148;
t571 = t302 * t412 + t303 * t414 - t459;
t329 = qJD(1) ^ 2;
t214 = qJD(1) * t230;
t569 = t214 / 0.2e1;
t215 = qJD(1) * t498;
t568 = -t215 / 0.2e1;
t567 = -t230 / 0.2e1;
t566 = t230 / 0.2e1;
t565 = t498 / 0.2e1;
t564 = -t498 / 0.2e1;
t563 = -t310 / 0.2e1;
t562 = -t312 / 0.2e1;
t561 = rSges(4,3) + pkin(6);
t318 = cos(qJ(1)) * pkin(1);
t556 = -qJD(1) / 0.2e1;
t555 = qJD(1) / 0.2e1;
t554 = rSges(4,1) * t326;
t550 = rSges(4,2) * t324;
t549 = rSges(5,2) * t309;
t546 = rSges(4,3) * t312;
t545 = rSges(5,3) * t312;
t284 = rSges(4,2) * t497;
t437 = rSges(4,1) * t496;
t184 = rSges(4,3) * t310 - t284 + t437;
t306 = qJD(1) * t318;
t276 = rSges(4,1) * t324 + rSges(4,2) * t326;
t426 = t276 * t449;
t394 = t306 - t426;
t95 = (t184 + t244) * qJD(1) + t394;
t543 = t310 * t95;
t541 = rSges(5,3) - t323;
t540 = rSges(6,3) - t319;
t450 = qJD(1) * t312;
t295 = pkin(6) * t450;
t446 = qJD(3) * t324;
t433 = pkin(3) * t446;
t109 = t312 * t433 + t295 - t297 + (t288 + (-pkin(2) + t304) * t310) * qJD(1);
t238 = t396 * qJD(3);
t539 = -t109 - (t238 - t433) * t312 - (t460 * t312 + (t256 - t304) * t310) * qJD(1);
t520 = t150 * t302;
t518 = t152 * t303;
t511 = t220 * t310;
t170 = t220 * t312;
t200 = t310 * t238;
t78 = -t222 * t507 + t505 * t588 + t511;
t494 = t78 * qJD(1);
t193 = t240 * t312;
t241 = rSges(5,1) * t311 - t549;
t492 = -qJD(3) * t193 - (t241 * t310 - t545) * qJD(1) - t109;
t491 = -t123 + t155;
t436 = rSges(5,1) * t499;
t398 = -rSges(5,2) * t503 + t436;
t164 = rSges(5,3) * t310 + t398;
t489 = -t147 + t164;
t480 = t177 + t205;
t479 = t178 + t206;
t478 = t179 - t203;
t477 = -Icges(4,2) * t496 + t180 - t282;
t476 = t256 * t450 - t200;
t305 = t329 * t318;
t451 = qJD(1) * t310;
t462 = pkin(2) * t450 + pkin(6) * t451;
t475 = qJD(1) * t462 + t305;
t469 = rSges(6,3) * t451 + qJD(1) * t434;
t468 = rSges(5,3) * t451 + qJD(1) * t436;
t465 = rSges(4,3) * t451 + qJD(1) * t437;
t456 = qJD(1) * t164;
t453 = qJD(1) * t233;
t452 = qJD(1) * t271;
t448 = qJD(3) * t311;
t445 = qJD(3) * t326;
t444 = qJD(4) * t312;
t441 = pkin(3) * t501;
t328 = qJD(3) ^ 2;
t440 = t328 * t317;
t286 = pkin(3) * t497;
t439 = t329 * t316;
t91 = t320 * t182 + (t228 * t310 - t544) * qJD(1);
t438 = -t91 + t539;
t435 = t320 * t552;
t432 = pkin(3) * t445;
t252 = t304 * t450;
t405 = t252 - t444;
t110 = (-qJD(1) * t323 - t433) * t310 + t405 - t462;
t137 = t146 * t450;
t430 = qJD(3) * t137 + (t110 + t142) * t449;
t429 = t310 * t110 + t147 * t451 + t137;
t92 = -t310 * t435 + (-t230 * t303 - t302 * t450) * rSges(6,2) + t469;
t428 = t154 * t450 - t155 * t451 + t310 * t92;
t427 = -t147 + t491;
t425 = t276 * t447;
t423 = t451 / 0.2e1;
t422 = -t450 / 0.2e1;
t420 = t449 / 0.2e1;
t418 = t447 / 0.2e1;
t410 = -t149 - t520;
t409 = -t149 + t518;
t404 = qJD(1) * t110 + t312 * t440 + t475;
t403 = qJD(1) * t297 - t439;
t395 = -t317 - t557;
t242 = rSges(3,1) * t312 - rSges(3,2) * t310;
t392 = -t550 + t554;
t283 = rSges(4,1) * t500;
t183 = -rSges(4,2) * t501 + t283 - t546;
t94 = t425 + (t183 + t417) * qJD(1);
t391 = t312 * t94 - t543;
t76 = -t151 * t303 - t153 * t302;
t375 = t183 * t310 + t184 * t312;
t374 = -t222 * t302 + t303 * t588;
t369 = -t227 - t396;
t208 = t276 * t312;
t367 = qJD(3) * t276;
t163 = rSges(5,1) * t502 - rSges(5,2) * t504 - t545;
t357 = -qJD(1) * t221 + t170 * t230 - t498 * t511;
t354 = qJD(1) * t382 - t170 * t320 - t459;
t353 = t320 * t511 + (-t221 * t312 + t518 - t520 - t521) * qJD(1);
t348 = t324 * t478 + t326 * t480;
t347 = t324 * t477 + t326 * t479;
t344 = qJD(1) * t374 - t221 * t320;
t341 = -t310 * t390 - t444;
t10 = t354 * t310 - t312 * t572;
t11 = -t310 * t571 + t353 * t312;
t12 = t310 * t572 + t354 * t312;
t124 = t152 * t506;
t56 = -t148 * t312 - t150 * t508 + t124;
t125 = t153 * t506;
t57 = t149 * t312 + t151 * t508 - t125;
t77 = t310 * t374 - t170;
t74 = t77 * qJD(1);
t21 = -t230 * t57 - t498 * t56 + t74;
t126 = t150 * t507;
t58 = -t148 * t310 - t152 * t505 + t126;
t59 = t149 * t310 - t312 * t382;
t22 = -t230 * t59 - t498 * t58 - t494;
t36 = t344 * t310 + t312 * t573;
t37 = -t310 * t573 + t344 * t312;
t38 = -t302 * t414 + t303 * t412;
t39 = t302 * t413 - t303 * t411;
t75 = t150 * t303 + t152 * t302;
t9 = t353 * t310 + t312 * t571;
t340 = (qJD(1) * t36 - t10 * t230 + t214 * t58 - t215 * t59 - t498 * t9) * t563 + (-t302 * t570 + t303 * t335) * t556 + t21 * t423 + (qJD(1) * t37 - t11 * t498 - t12 * t230 + t214 * t56 - t215 * t57) * t562 + t22 * t422 + (-t10 * t310 - t312 * t9 + (t310 * t58 - t312 * t59) * qJD(1)) * t567 + (-t310 * t57 - t312 * t56) * t569 + (-t310 * t59 - t312 * t58) * t568 + (-t11 * t312 - t12 * t310 + (t310 * t56 - t312 * t57) * qJD(1)) * t564 + (-t310 * t39 - t312 * t38 + (t310 * t75 - t312 * t76) * qJD(1)) * t555 + (t357 * t310 + t312 * t607) * t566 + (-t310 * t607 + t357 * t312) * t565;
t339 = t306 + t341;
t338 = -t444 - t200;
t120 = qJD(3) * t208 + (t310 * t392 - t546) * qJD(1);
t121 = -rSges(4,1) * t310 * t446 + (-t310 * t445 - t324 * t450) * rSges(4,2) + t465;
t337 = -t120 * t312 + t121 * t310 + (t183 * t312 - t184 * t310) * qJD(1);
t330 = -t230 * t227 + t306 + t338;
t257 = t392 * qJD(3);
t226 = pkin(2) * t451 - t295;
t219 = t241 * qJD(3);
t207 = t276 * t310;
t197 = t228 * t320;
t192 = t240 * t310;
t168 = t312 * t396 - t286;
t167 = (-t396 + t559) * t310;
t156 = t230 * t228;
t141 = t310 * t154;
t140 = t310 * t146;
t105 = -t449 * t553 + (-t309 * t450 - t310 * t448) * rSges(5,2) + t468;
t93 = qJD(3) * t375 + qJD(2);
t73 = -t252 + (-qJD(1) * t460 + t433) * t310 + t476;
t67 = t257 * t447 + (t121 - t426) * qJD(1) + t475;
t66 = -t439 - t257 * t449 + (-t120 - t226 - t425) * qJD(1);
t61 = (t244 + t489) * qJD(1) + t339;
t60 = -t297 + t416 * t447 + (t163 + t402) * qJD(1);
t52 = (t163 * t310 + t312 * t489) * qJD(3) + t442;
t49 = t337 * qJD(3);
t48 = (t244 + t427) * qJD(1) + t330;
t41 = t219 * t447 + (t105 + t341) * qJD(1) + t404;
t40 = (-qJD(3) * t219 - t440) * t310 + (-t312 * t390 - t226 + t492) * qJD(1) + t403;
t24 = pkin(4) * t328 * t499 + t197 * t498 - t214 * t227 + (t73 + t92 + t338) * qJD(1) + t404;
t23 = -t197 * t230 - t215 * t227 + t395 * t328 * t310 + (-t238 * t312 - t226 + t438) * qJD(1) + t403;
t14 = ((t105 - t456) * t310 + (qJD(1) * t163 + t492) * t312) * qJD(3) + t430;
t5 = t154 * t215 - t155 * t214 + t230 * t92 - t498 * t91 + ((qJD(1) * t123 + t73) * t310 + (qJD(1) * t122 + t539) * t312) * qJD(3) + t430;
t1 = [m(3) * ((t329 * t393 - t439) * (t242 + t318) + (-t305 + (-0.2e1 * rSges(3,1) * t450 + 0.2e1 * rSges(3,2) * t451 + qJD(1) * t242) * qJD(1)) * (t393 - t316)) + (t74 - (t310 * t410 + t124 + t59) * t498 + (t125 - t126 + t58 + (t148 - t519) * t310) * t230 + (t230 * t409 - t592) * t312) * t566 + (qJD(3) * t370 + t254 * t326 + t255 * t324) * qJD(1) + (qJD(3) * t372 + t217 * t311 + t218 * t309) * qJD(1) + t76 * t568 + (-t302 * t406 + t303 * t407) * qJD(1) + t215 * t78 / 0.2e1 + (t75 + t77) * t569 + (t494 - (-t126 + t57) * t498 + (-t124 + t56) * t230 + (-t230 * t368 - t409 * t498) * t312 + (-t230 * t410 + t592) * t310 + t22) * t565 + (t38 + t37) * t564 + ((((-t641 + t614) * t312 + t616 + t618) * t312 + ((t614 + t644) * t310 + (t609 - t610) * t312 - t615 + t619 + t640) * t310) * qJD(3) + t631) * t420 + (-(qJD(1) * t491 + t330 - t48 + t589) * t47 + t23 * (t213 - t251 + t318) + t48 * t297 + t24 * (t250 + t316 + t474) + t47 * (t306 + t469 + t476) + (t23 * t551 + t48 * (-t320 * t547 - t238 - t435) - t24 * rSges(6,3) - t47 * qJD(4)) * t312 + (-t227 * t320 * t47 + t23 * t540 - t24 * t548) * t310 + (-t48 * t316 + (-t47 * t548 + t48 * t540) * t312 + (t48 * (-t228 - t256) - t47 * t319) * t310) * qJD(1)) * m(6) + (-(t339 - t61 + t456 + t589) * t60 + t40 * (t310 * t541 + t266 + t318 + t398) + t61 * t297 + t41 * (t163 + t316 + t466) + t60 * (t306 + t405 + t468) + (-t61 * t316 + (t541 * t61 - t549 * t60) * t312 + (t61 * (-t241 - t304) - t60 * t323) * t310) * qJD(1) - (t310 * t60 + t312 * t61) * t390) * m(5) + (t66 * (-t284 + t318) + t95 * t295 + t67 * (t283 + t301 + t316) + (-t550 * t67 + t561 * t66) * t310 + (t66 * (pkin(2) + t554) - t67 * t561 - t95 * t367) * t312 + (t95 * (t546 - t316) + (-pkin(2) - t392) * t543) * qJD(1) + (-t367 * t310 - t229 + t306 - t394 + t462 + t465 + t95 + (-t184 - t284) * qJD(1)) * t94) * m(4) + (t39 + t36 + t21) * t567 + ((((t610 + t613) * t312 + t615 + t620) * t312 + ((t609 + t613) * t310 + t616 + t621) * t310) * qJD(3) + t626 + t628) * t418 - (t623 - t624 + t627) * t449 / 0.2e1 - (-qJD(1) * t611 + t622 - t625) * t447 / 0.2e1 + (t635 * t312 + (t617 + t636) * t310) * qJD(3) * t555; m(4) * t49 + m(5) * t14 + m(6) * t5; t340 + (((t310 * t477 - t312 * t478) * t326 + (-t310 * t479 + t312 * t480) * t324 - t309 * t575 + t311 * t349) * qJD(3) + (-t309 * t471 + t311 * t470 - t324 * t464 + t326 * t463) * qJD(1)) * t556 + (t625 * t312 + t624 * t310 + (t310 * t617 + t312 * t611) * qJD(1)) * t555 + ((t187 * t449 - t453) * t310 + (t594 + (-t310 * t510 - t606) * qJD(3)) * t312 + (t202 * t449 - t452) * t310 + (t595 + (-t348 * t312 + (-t509 + t347) * t310) * qJD(3)) * t312) * t420 + ((-t447 * t510 - t453) * t312 + (-t594 + (t312 * t187 + t606) * qJD(3)) * t310 + (-t447 * t509 - t452) * t312 + (-t595 + (-t347 * t310 + (t202 + t348) * t312) * qJD(3)) * t310) * t418 + (t5 * (t140 + t141) + t33 * (t428 + t429) + t24 * t286 + (t5 * t427 + t33 * t438 + t24 * (t227 + t558) + t47 * t197 + (t33 * t122 + t369 * t48) * qJD(1)) * t312 + (t5 * t122 + t33 * t73 + t23 * t369 + t48 * (-pkin(4) * t448 - t197 - t432) + (t33 * t123 + t369 * t47) * qJD(1)) * t310 + t48 * t156 - (t48 * (-t168 - t182 - t286) + t47 * (t167 - t441)) * qJD(1) - (t33 * (-t168 * t312 + t583) + (t33 * t167 + t395 * t48) * t310) * qJD(3) + t584) * m(6) + (-(-t207 * t94 - t208 * t95) * qJD(1) - (t93 * (-t207 * t310 - t208 * t312) + t391 * t392) * qJD(3) + t49 * t375 + t93 * t337 + t391 * t257 + (-t66 * t310 + t67 * t312 + (-t310 * t94 - t312 * t95) * qJD(1)) * t276) * m(4) + (t623 * qJD(1) + ((t618 * qJD(1) + t597 * t312) * t312 + (t600 * t310 + t619 * qJD(1) + (-t598 + t601) * t312) * t310) * t593) * t563 + (t622 * qJD(1) + ((t620 * qJD(1) + t601 * t312) * t312 + (t598 * t310 + t621 * qJD(1) + (-t597 + t600) * t312) * t310) * t593) * t562 + (-(t61 * (-t193 - t286) + t60 * (-t192 - t441)) * qJD(1) - (t52 * (-t193 * t312 + t583) + t60 * t241 * t312 + (-t52 * t192 + t61 * (-t241 - t317)) * t310) * qJD(3) + t14 * t140 + t52 * t429 + t41 * t286 + (t14 * t489 + t52 * t492 + t41 * t240 + t60 * t219 + (t52 * t163 - t416 * t61) * qJD(1)) * t312 + (t14 * t163 + t52 * t105 - t40 * t416 + t61 * (-t219 - t432) + (-t52 * t164 - t416 * t60) * qJD(1)) * t310) * m(5) + (t627 + t629) * t423 + (t626 + t630) * t422; m(5) * (-t310 * t41 - t312 * t40) + m(6) * (-t23 * t312 - t24 * t310); t340 + (t5 * (t155 * t312 + t141) + t33 * (-t312 * t91 + t428) + (-t310 * t48 + t312 * t47) * t197 + (-t23 * t310 + t24 * t312 + (-t310 * t47 - t312 * t48) * qJD(1)) * t227 - t48 * (-qJD(1) * t182 - t156) + t584) * m(6);];
tauc = t1(:);

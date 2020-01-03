% Calculate vector of inverse dynamics joint torques for
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR8_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR8_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:19
% EndTime: 2019-12-31 17:42:41
% DurationCPUTime: 19.09s
% Computational Cost: add. (19446->706), mult. (19856->1025), div. (0->0), fcn. (18404->10), ass. (0->380)
t344 = qJ(2) + qJ(3);
t336 = pkin(9) + t344;
t326 = sin(t336);
t327 = cos(t336);
t337 = sin(t344);
t338 = cos(t344);
t633 = -Icges(4,5) * t337 - Icges(5,5) * t326 - Icges(4,6) * t338 - Icges(5,6) * t327;
t632 = Icges(4,3) + Icges(5,3);
t345 = sin(pkin(8));
t631 = t633 * t345;
t346 = cos(pkin(8));
t630 = t633 * t346;
t629 = Icges(4,5) * t338 + Icges(5,5) * t327 - Icges(4,6) * t337 - Icges(5,6) * t326;
t343 = qJD(2) + qJD(3);
t628 = t631 * t343;
t627 = t630 * t343;
t626 = t629 * t345 - t632 * t346;
t625 = t632 * t345 + t629 * t346;
t547 = Icges(4,4) * t337;
t453 = Icges(4,1) * t338 - t547;
t227 = Icges(4,5) * t345 + t346 * t453;
t444 = -Icges(4,2) * t338 - t547;
t603 = -t444 * t346 - t227;
t546 = Icges(4,4) * t338;
t445 = -Icges(4,2) * t337 + t546;
t225 = Icges(4,6) * t345 + t346 * t445;
t452 = -Icges(4,1) * t337 - t546;
t605 = t452 * t346 - t225;
t545 = Icges(5,4) * t326;
t451 = Icges(5,1) * t327 - t545;
t204 = Icges(5,5) * t345 + t346 * t451;
t442 = -Icges(5,2) * t327 - t545;
t607 = -t442 * t346 - t204;
t544 = Icges(5,4) * t327;
t443 = -Icges(5,2) * t326 + t544;
t202 = Icges(5,6) * t345 + t346 * t443;
t450 = -Icges(5,1) * t326 - t544;
t609 = t450 * t346 - t202;
t624 = (t326 * t607 + t327 * t609 + t337 * t603 + t338 * t605) * t343;
t226 = -Icges(4,5) * t346 + t345 * t453;
t604 = -t444 * t345 - t226;
t224 = -Icges(4,6) * t346 + t345 * t445;
t606 = t452 * t345 - t224;
t203 = -Icges(5,5) * t346 + t345 * t451;
t608 = -t442 * t345 - t203;
t201 = -Icges(5,6) * t346 + t345 * t443;
t610 = t450 * t345 - t201;
t623 = (-t326 * t608 - t327 * t610 - t337 * t604 - t338 * t606) * t343;
t622 = -t202 * t326 + t204 * t327 - t225 * t337 + t227 * t338;
t621 = -t201 * t326 + t203 * t327 - t224 * t337 + t226 * t338;
t348 = sin(qJ(2));
t350 = cos(qJ(2));
t318 = rSges(3,1) * t348 + rSges(3,2) * t350;
t341 = t345 ^ 2;
t342 = t346 ^ 2;
t583 = t341 + t342;
t620 = t318 * t583;
t619 = t623 * t345 + t628 * t346;
t618 = t624 * t345 - t627 * t346;
t617 = -t628 * t345 + t623 * t346;
t616 = t627 * t345 + t624 * t346;
t615 = t621 * t345 - t626 * t346;
t614 = t622 * t345 - t625 * t346;
t613 = t626 * t345 + t621 * t346;
t612 = t625 * t345 + t622 * t346;
t286 = rSges(5,1) * t326 + rSges(5,2) * t327;
t404 = t286 * t345;
t192 = t343 * t404;
t239 = t286 * t346;
t193 = t343 * t239;
t335 = qJD(2) * t345;
t310 = qJD(3) * t345 + t335;
t611 = -t345 * t192 - t346 * t193 + t310 * t404;
t529 = t343 * t346;
t602 = (t310 * t605 - t529 * t606) * t338 + (t310 * t603 - t529 * t604) * t337 + (t310 * t609 - t529 * t610) * t327 + (t310 * t607 - t529 * t608) * t326;
t561 = pkin(3) * t337;
t562 = pkin(2) * t348;
t309 = -t561 - t562;
t299 = t345 * t309;
t466 = t343 * t561;
t505 = qJD(4) * t346;
t160 = -t345 * t466 - t505;
t333 = qJD(4) * t345;
t161 = -t346 * t466 + t333;
t519 = t345 * t160 + t346 * t161;
t601 = t519 - t310 * (t345 * t562 + t299);
t600 = t630 * t310 - t631 * t529;
t347 = sin(qJ(5));
t349 = cos(qJ(5));
t556 = rSges(6,1) * t349;
t461 = -rSges(6,2) * t347 + t556;
t527 = t345 * t349;
t495 = t326 * t527;
t528 = t345 * t347;
t496 = t326 * t528;
t532 = t327 * t345;
t511 = rSges(6,2) * t496 + rSges(6,3) * t532;
t147 = -rSges(6,1) * t495 + t511;
t525 = t346 * t349;
t493 = t326 * t525;
t526 = t346 * t347;
t494 = t326 * t526;
t531 = t327 * t346;
t510 = rSges(6,2) * t494 + rSges(6,3) * t531;
t148 = -rSges(6,1) * t493 + t510;
t288 = pkin(4) * t326 - pkin(7) * t327;
t409 = t343 * t288;
t194 = t345 * t409;
t195 = t346 * t409;
t504 = qJD(5) * t326;
t240 = t346 * t504 + t310;
t501 = qJD(5) * t345;
t241 = t326 * t501 - t529;
t314 = pkin(7) * t532;
t269 = -t327 * t528 - t525;
t270 = t327 * t527 - t526;
t535 = t326 * t345;
t120 = rSges(6,1) * t270 + rSges(6,2) * t269 + rSges(6,3) * t535;
t503 = qJD(5) * t327;
t488 = t120 * t503;
t149 = -qJD(5) * t270 + t343 * t496;
t150 = qJD(5) * t269 - t343 * t495;
t492 = t343 * t532;
t91 = rSges(6,1) * t150 + rSges(6,2) * t149 + rSges(6,3) * t492;
t272 = t327 * t525 + t528;
t151 = -qJD(5) * t272 + t343 * t494;
t271 = -t327 * t526 + t527;
t152 = qJD(5) * t271 - t343 * t493;
t491 = t327 * t529;
t92 = rSges(6,1) * t152 + rSges(6,2) * t151 + rSges(6,3) * t491;
t599 = (-t194 + t91) * t345 - t240 * t147 + t148 * t241 - t310 * (-pkin(4) * t535 + t314) + (-t195 + t92 - t488) * t346;
t548 = Icges(3,4) * t350;
t549 = Icges(3,4) * t348;
t598 = (-t348 * (-Icges(3,2) * t350 - t549) + t350 * (-Icges(3,1) * t348 - t548)) * qJD(2);
t114 = Icges(6,5) * t270 + Icges(6,6) * t269 + Icges(6,3) * t535;
t543 = Icges(6,4) * t270;
t116 = Icges(6,2) * t269 + Icges(6,6) * t535 + t543;
t249 = Icges(6,4) * t269;
t118 = Icges(6,1) * t270 + Icges(6,5) * t535 + t249;
t50 = t114 * t535 + t116 * t269 + t118 * t270;
t534 = t326 * t346;
t115 = Icges(6,5) * t272 + Icges(6,6) * t271 + Icges(6,3) * t534;
t542 = Icges(6,4) * t272;
t117 = Icges(6,2) * t271 + Icges(6,6) * t534 + t542;
t250 = Icges(6,4) * t271;
t119 = Icges(6,1) * t272 + Icges(6,5) * t534 + t250;
t51 = t115 * t535 + t117 * t269 + t119 * t270;
t52 = t114 * t534 + t116 * t271 + t118 * t272;
t53 = t115 * t534 + t117 * t271 + t119 * t272;
t433 = Icges(6,5) * t349 - Icges(6,6) * t347;
t175 = -Icges(6,3) * t327 + t326 * t433;
t540 = Icges(6,4) * t349;
t441 = -Icges(6,2) * t347 + t540;
t177 = -Icges(6,6) * t327 + t326 * t441;
t541 = Icges(6,4) * t347;
t449 = Icges(6,1) * t349 - t541;
t179 = -Icges(6,5) * t327 + t326 * t449;
t68 = t175 * t535 + t177 * t269 + t179 * t270;
t69 = t175 * t534 + t177 * t271 + t179 * t272;
t597 = (t240 * t53 + t241 * t52 - t503 * t69) * t346 + (t240 * t51 + t241 * t50 - t503 * t68) * t345;
t303 = rSges(4,1) * t337 + rSges(4,2) * t338;
t273 = t303 * t345;
t216 = t343 * t273;
t274 = t303 * t346;
t217 = t343 * t274;
t584 = t338 * rSges(4,1) - rSges(4,2) * t337;
t228 = -rSges(4,3) * t346 + t345 * t584;
t229 = rSges(4,3) * t345 + t346 * t584;
t339 = t350 * pkin(2);
t220 = -pkin(6) * t346 + t339 * t345;
t221 = pkin(6) * t345 + t339 * t346;
t506 = qJD(2) * t346;
t484 = t220 * t335 + t221 * t506 + qJD(1);
t596 = (-t345 * t216 - t346 * t217 + t310 * t273 + t274 * t529) * (t228 * t310 + t229 * t529 + t484);
t328 = pkin(3) * t338;
t586 = t327 * rSges(5,1) - rSges(5,2) * t326;
t465 = t586 + t328;
t593 = t310 * t465;
t460 = -rSges(6,1) * t347 - rSges(6,2) * t349;
t533 = t327 * t343;
t107 = t461 * t533 + (rSges(6,3) * t343 + qJD(5) * t460) * t326;
t585 = t327 * pkin(4) + t326 * pkin(7);
t231 = t585 * t343;
t589 = -t107 - t231;
t587 = -t585 - t328;
t284 = t529 * t561;
t582 = t346 * t284 + t519;
t581 = g(1) * t346 + g(2) * t345;
t429 = -t117 * t347 + t119 * t349;
t430 = -t116 * t347 + t118 * t349;
t579 = -(-t175 * t346 - t429) * t240 - (-t175 * t345 - t430) * t241;
t440 = -Icges(6,2) * t349 - t541;
t359 = t240 * (-Icges(6,2) * t272 + t119 + t250) + t241 * (-Icges(6,2) * t270 + t118 + t249) - t503 * (t440 * t326 + t179);
t352 = qJD(2) ^ 2;
t332 = qJDD(2) * t345;
t307 = qJDD(3) * t345 + t332;
t502 = qJD(5) * t343;
t399 = qJDD(5) * t326 + t327 * t502;
t155 = t346 * t399 + t307;
t577 = t155 / 0.2e1;
t308 = (-qJDD(2) - qJDD(3)) * t346;
t156 = t345 * t399 + t308;
t576 = t156 / 0.2e1;
t575 = -t240 / 0.2e1;
t574 = t240 / 0.2e1;
t573 = -t241 / 0.2e1;
t572 = t241 / 0.2e1;
t251 = -qJDD(5) * t327 + t326 * t502;
t571 = t251 / 0.2e1;
t564 = t345 / 0.2e1;
t563 = -t346 / 0.2e1;
t552 = pkin(2) * qJD(2);
t551 = t345 * t52;
t550 = t346 * t51;
t536 = t175 * t327;
t530 = t338 * t343;
t197 = t346 * t221;
t139 = -qJ(4) * t346 + t328 * t345;
t140 = qJ(4) * t345 + t328 * t346;
t522 = t345 * t139 + t346 * t140;
t209 = rSges(5,3) * t345 + t346 * t586;
t521 = -t140 - t209;
t244 = t585 * t346;
t520 = -t140 - t244;
t516 = t345 * t220 + t197;
t515 = t345 * t228 + t346 * t229;
t285 = t529 * t328;
t512 = -t529 * t586 - t285;
t500 = t345 * t561;
t499 = pkin(3) * t530;
t498 = t348 * t552;
t497 = t350 * t552;
t490 = t314 + t511;
t315 = pkin(7) * t531;
t489 = t315 + t510;
t121 = rSges(6,1) * t272 + rSges(6,2) * t271 + rSges(6,3) * t534;
t487 = t121 * t503;
t486 = t327 * t501;
t481 = -t503 / 0.2e1;
t480 = t503 / 0.2e1;
t385 = -t303 - t562;
t384 = -t286 - t561;
t478 = -t288 - t561;
t476 = t345 * t498;
t475 = t345 * t497;
t474 = t346 * t498;
t473 = t346 * t497;
t208 = -rSges(5,3) * t346 + t345 * t586;
t472 = t345 * t208 + t346 * t209 + t522;
t471 = t310 * t139 + t484;
t230 = t586 * t343;
t469 = -t230 - t499;
t190 = -rSges(6,3) * t327 + t326 * t461;
t468 = -t190 + t478;
t259 = t584 * t343;
t467 = -t259 - t497;
t319 = rSges(3,1) * t350 - rSges(3,2) * t348;
t397 = pkin(2) * (-qJDD(2) * t348 - t350 * t352);
t383 = t345 * t397;
t357 = (-t307 * t337 - t310 * t530) * pkin(3) + t383 - qJDD(4) * t346;
t30 = -t107 * t240 + t121 * t251 - t155 * t190 - t231 * t310 - t288 * t307 - t503 * t92 + t357;
t382 = t346 * t397;
t368 = qJDD(4) * t345 - t343 * t285 + t308 * t561 + t382;
t31 = t107 * t241 - t120 * t251 + t156 * t190 - t231 * t529 + t288 * t308 + t503 * t91 + t368;
t459 = -t30 * t346 + t31 * t345;
t458 = t345 * t50 + t550;
t457 = t346 * t53 + t551;
t60 = -t114 * t327 + t326 * t430;
t61 = -t115 * t327 + t326 * t429;
t456 = t60 * t345 + t61 * t346;
t455 = Icges(3,1) * t350 - t549;
t448 = -Icges(6,1) * t347 - t540;
t447 = -Icges(3,2) * t348 + t548;
t439 = Icges(3,5) * t350 - Icges(3,6) * t348;
t438 = -Icges(3,5) * t348 - Icges(3,6) * t350;
t432 = -Icges(6,5) * t347 - Icges(6,6) * t349;
t431 = t114 * t241 + t115 * t240;
t428 = t120 * t346 - t121 * t345;
t427 = -t177 * t347 + t179 * t349;
t254 = -Icges(3,6) * t346 + t345 * t447;
t256 = -Icges(3,5) * t346 + t345 * t455;
t421 = -t254 * t348 + t256 * t350;
t255 = Icges(3,6) * t345 + t346 * t447;
t257 = Icges(3,5) * t345 + t346 * t455;
t420 = -t255 * t348 + t257 * t350;
t419 = t583 * t319;
t417 = qJD(2) * t620;
t416 = -t499 + t589;
t191 = t326 * rSges(6,3) + t327 * t461;
t414 = -t286 + t309;
t242 = t585 * t345;
t412 = t522 + (t121 + t244) * t346 + (t120 + t242) * t345;
t85 = Icges(6,5) * t150 + Icges(6,6) * t149 + Icges(6,3) * t492;
t411 = t114 * t533 + t326 * t85;
t86 = Icges(6,5) * t152 + Icges(6,6) * t151 + Icges(6,3) * t491;
t410 = t115 * t533 + t326 * t86;
t394 = t433 * t327;
t104 = t343 * t394 + (Icges(6,3) * t343 + qJD(5) * t432) * t326;
t408 = t104 * t326 + t175 * t533;
t405 = Icges(6,3) * t326 + t394 - t427;
t403 = -t284 + t333 - t474;
t398 = -t190 - t288 + t309;
t396 = t449 * t327;
t395 = t441 * t327;
t392 = -t497 - t499;
t389 = qJD(2) * t438;
t386 = -t476 - t505;
t381 = -t230 + t392;
t378 = t191 - t587;
t377 = t392 + t589;
t376 = -t120 * t504 + t147 * t503 + t190 * t486 + t241 * t191 - t529 * t585 - t285;
t375 = -(Icges(6,5) * t269 - Icges(6,6) * t270) * t241 - (Icges(6,5) * t271 - Icges(6,6) * t272) * t240 + t432 * t326 * t503;
t374 = -t352 * t562 * t583 + qJDD(2) * t197 + t220 * t332 + qJDD(1);
t370 = t326 * t375;
t369 = t307 * t139 + t310 * t160 + t374;
t363 = (-t254 * t350 - t256 * t348) * qJD(2) + t598 * t345;
t362 = (-t255 * t350 - t257 * t348) * qJD(2) + t598 * t346;
t361 = t581 * (-pkin(4) - t556) * t326;
t360 = (Icges(6,1) * t271 - t117 - t542) * t240 + (Icges(6,1) * t269 - t116 - t543) * t241 - (t448 * t326 - t177) * t503;
t358 = -t191 * t240 + t121 * t504 + t587 * t310 + (-t190 * t346 - t148) * t503;
t87 = Icges(6,4) * t150 + Icges(6,2) * t149 + Icges(6,6) * t492;
t89 = Icges(6,1) * t150 + Icges(6,4) * t149 + Icges(6,5) * t492;
t14 = (t343 * t430 - t85) * t327 + (t114 * t343 - t347 * t87 + t349 * t89 + (-t116 * t349 - t118 * t347) * qJD(5)) * t326;
t143 = t177 * t345;
t144 = t177 * t346;
t145 = t179 * t345;
t146 = t179 * t346;
t88 = Icges(6,4) * t152 + Icges(6,2) * t151 + Icges(6,6) * t491;
t90 = Icges(6,1) * t152 + Icges(6,4) * t151 + Icges(6,5) * t491;
t15 = (t343 * t429 - t86) * t327 + (t115 * t343 - t347 * t88 + t349 * t90 + (-t117 * t349 - t119 * t347) * qJD(5)) * t326;
t178 = Icges(6,6) * t326 + t395;
t180 = Icges(6,5) * t326 + t396;
t20 = t116 * t149 + t118 * t150 + t269 * t87 + t270 * t89 + t345 * t411;
t21 = t117 * t149 + t119 * t150 + t269 * t88 + t270 * t90 + t345 * t410;
t22 = t116 * t151 + t118 * t152 + t271 * t87 + t272 * t89 + t346 * t411;
t23 = t117 * t151 + t119 * t152 + t271 * t88 + t272 * t90 + t346 * t410;
t77 = t326 * t427 - t536;
t27 = t240 * t61 + t241 * t60 - t503 * t77;
t105 = t343 * t395 + (Icges(6,6) * t343 + qJD(5) * t440) * t326;
t106 = t343 * t396 + (Icges(6,5) * t343 + qJD(5) * t448) * t326;
t34 = t105 * t269 + t106 * t270 + t149 * t177 + t150 * t179 + t345 * t408;
t3 = t155 * t51 + t156 * t50 + t20 * t241 + t21 * t240 + t251 * t68 - t34 * t503;
t353 = (-t405 * t503 - t579) * t326;
t35 = t105 * t271 + t106 * t272 + t151 * t177 + t152 * t179 + t346 * t408;
t4 = t155 * t53 + t156 * t52 + t22 * t241 + t23 * t240 + t251 * t69 - t35 * t503;
t354 = (((t144 * t347 - t146 * t349 + t115) * t240 + (t143 * t347 - t145 * t349 + t114) * t241 + t77 * qJD(5)) * t326 + ((t405 * t327 + (t178 * t347 - t180 * t349 - t175) * t326 + t456) * qJD(5) + t579) * t327) * t480 - t27 * t504 / 0.2e1 + (-t20 * t346 + t21 * t345) * t572 + ((-t144 * t269 - t146 * t270) * t240 + (-t143 * t269 - t145 * t270) * t241 + (t68 * t326 + (-t178 * t269 - t180 * t270 + t550) * t327) * qJD(5) + (((t50 - t536) * qJD(5) + t431) * t327 + t353) * t345) * t573 + (-t22 * t346 + t23 * t345) * t574 + ((-t144 * t271 - t146 * t272) * t240 + (-t143 * t271 - t145 * t272) * t241 + (t69 * t326 + (-t178 * t271 - t180 * t272 + t551) * t327) * qJD(5) + (((t53 - t536) * qJD(5) + t431) * t327 + t353) * t346) * t575 + (t345 * t51 - t346 * t50) * t576 + (t345 * t53 - t346 * t52) * t577 + (t345 * t61 - t346 * t60) * t571 + (t345 * t612 - t346 * t613) * t307 / 0.2e1 + (t345 * t614 - t346 * t615) * t308 / 0.2e1 - (t345 * t600 + t346 * t602) * t310 / 0.2e1 + (t345 * t616 + t346 * t617) * t310 / 0.2e1 + (t345 * t602 - t346 * t600) * t529 / 0.2e1 - (t345 * t618 + t346 * t619) * t529 / 0.2e1 + (t307 * t612 + t308 * t613 + t310 * t616 + t529 * t617 + t4) * t564 + (t307 * t614 + t308 * t615 + t310 * t618 + t529 * t619 + t3) * t563 + (-t14 * t346 + t15 * t345 + t597) * t481;
t317 = g(1) * t345 - g(2) * t346;
t300 = t346 * t309;
t298 = t318 * t346;
t297 = t318 * t345;
t292 = t438 * t346;
t291 = t438 * t345;
t276 = t346 * t389;
t275 = t345 * t389;
t253 = Icges(3,3) * t345 + t346 * t439;
t252 = -Icges(3,3) * t346 + t345 * t439;
t248 = t460 * t326;
t243 = -pkin(4) * t534 + t315;
t238 = t346 * t562 + t300;
t189 = -t310 * t346 + t345 * t529;
t158 = -t303 * t529 - t474;
t157 = -t303 * t310 - t476;
t135 = rSges(6,1) * t271 - rSges(6,2) * t272;
t134 = rSges(6,1) * t269 - rSges(6,2) * t270;
t111 = -t286 * t529 + t403;
t110 = t310 * t384 + t386;
t102 = -t259 * t529 + t303 * t308 + t382;
t101 = -t259 * t310 - t303 * t307 + t383;
t82 = -qJD(2) * t417 + qJDD(2) * t419 + qJDD(1);
t76 = -t230 * t529 + t286 * t308 + t368;
t75 = -t230 * t310 - t286 * t307 + t357;
t63 = t190 * t241 - t288 * t529 + t403 + t488;
t62 = -t190 * t240 + t310 * t478 + t386 - t487;
t49 = -t216 * t310 - t217 * t529 + t228 * t307 - t229 * t308 + t374;
t46 = t208 * t310 - t521 * t529 + t471;
t39 = t120 * t240 - t121 * t241 + t242 * t310 - t520 * t529 + t471;
t33 = (t343 * t427 - t104) * t327 + (-t105 * t347 + t106 * t349 + t175 * t343 + (-t177 * t349 - t179 * t347) * qJD(5)) * t326;
t32 = -t192 * t310 + t208 * t307 - (-t161 + t193) * t529 + t521 * t308 + t369;
t11 = t120 * t155 - t121 * t156 - t194 * t310 + t240 * t91 - t241 * t92 + t242 * t307 - (-t161 + t195) * t529 + t520 * t308 + t369;
t1 = [m(2) * qJDD(1) + (-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3) + m(3) * t82 + m(4) * t49 + m(5) * t32 + m(6) * t11; -(t292 * qJD(2) * t341 - t345 * t291 * t506) * t335 / 0.2e1 + (t291 * qJD(2) * t342 - t346 * t292 * t335) * t506 / 0.2e1 + t354 + (-t63 * (t376 - t473) - t62 * (t358 - t475) + t11 * (t412 + t516) + (t31 * t398 + t377 * t63) * t346 + (t30 * t398 + t377 * t62) * t345 - g(1) * (t300 + t489) - g(2) * (t299 + t490) - g(3) * (t339 + t378) - t361 + (t121 * t486 + (-t238 - t243) * t529 + t599 + t601) * t39) * m(6) + (-g(1) * (t300 - t239) - g(2) * (t299 - t404) - g(3) * (t339 + t465) - t111 * (-t473 + t512) - t110 * (-t475 - t593) + t32 * (t472 + t516) + (t111 * t381 + t414 * t76) * t346 + (t110 * t381 + t414 * t75) * t345 + ((-t238 + t239) * t529 + t601 + t611) * t46) * m(5) + (t49 * (t515 + t516) + (t102 * t385 + t158 * t467) * t346 + (t101 * t385 + t157 * t467) * t345 - t158 * (-t529 * t584 - t473) - t157 * (-t310 * t584 - t475) - g(3) * (t584 + t339) - t581 * t385 + t596) * m(4) + (g(1) * t298 + g(2) * t297 - g(3) * t319 + t82 * t419 + (qJDD(2) * t318 + t319 * t352) * t620 + (-t417 - (-t297 * t345 - t298 * t346) * qJD(2)) * (qJD(2) * t419 + qJD(1))) * m(3) + 0.2e1 * (((t253 * t345 + t346 * t420) * t345 - t346 * (t252 * t345 + t346 * t421)) * qJDD(2) + (t345 * (t276 * t345 + t346 * t362) - t346 * (t275 * t345 + t346 * t363)) * qJD(2)) * t564 + 0.2e1 * ((t345 * (-t253 * t346 + t345 * t420) - t346 * (-t252 * t346 + t345 * t421)) * qJDD(2) + (t345 * (-t276 * t346 + t345 * t362) - t346 * (-t275 * t346 + t345 * t363)) * qJD(2)) * t563; t354 + (-g(1) * (-t346 * t561 + t489) - g(2) * (t490 - t500) - g(3) * t378 - t361 + t11 * t412 + (t31 * t468 + t416 * t63) * t346 + (t30 * t468 + t416 * t62) * t345 - t63 * t376 - t62 * t358 + (-t243 * t529 - (-t310 * t561 - t487) * t345 + t582 + t599) * t39) * m(6) + (-g(3) * t465 + t32 * t472 + (t345 * t469 + t593) * t110 + (t239 * t529 + t310 * t500 + t582 + t611) * t46 + (t346 * t469 - t512) * t111 + (t75 * t345 + t76 * t346 - t581) * t384) * m(5) + (t49 * t515 + (-t101 * t345 - t102 * t346) * t303 + (-t157 * t345 - t158 * t346) * t259 + g(1) * t274 + g(2) * t273 + t596 + (t157 * t310 + t158 * t529 - g(3)) * t584) * m(4); (-t39 * t189 - t317 + t459) * m(6) + (-t46 * t189 + t345 * t76 - t346 * t75 - t317) * m(5); t4 * t534 / 0.2e1 + (t326 * t457 - t327 * t69) * t577 + ((t343 * t457 - t35) * t327 + (t22 * t345 + t23 * t346 + t343 * t69) * t326) * t574 + t3 * t535 / 0.2e1 + (t326 * t458 - t327 * t68) * t576 + ((t343 * t458 - t34) * t327 + (t20 * t345 + t21 * t346 + t343 * t68) * t326) * t572 + t343 * t326 * t27 / 0.2e1 - t327 * (t14 * t241 + t15 * t240 + t155 * t61 + t156 * t60 + t251 * t77 - t33 * t503) / 0.2e1 + (t326 * t456 - t327 * t77) * t571 + ((t343 * t456 - t33) * t327 + (t14 * t345 + t15 * t346 + t343 * t77) * t326) * t481 + (t271 * t359 + t272 * t360 - t346 * t370) * t575 + (t269 * t359 + t270 * t360 - t345 * t370) * t573 + (t375 * t327 + (-t347 * t359 + t360 * t349) * t326) * t480 + t597 * t533 / 0.2e1 + ((t31 * t120 - t30 * t121 - t62 * t92 + t63 * t91 + (t39 * t428 + (t345 * t63 - t346 * t62) * t190) * t343) * t327 + (t63 * (t107 * t345 - t120 * t343) + t62 * (-t107 * t346 + t121 * t343) + t11 * t428 + t39 * (-t345 * t92 + t346 * t91) + t459 * t190) * t326 - t63 * (t134 * t503 + t241 * t248) - t62 * (-t135 * t503 - t240 * t248) - t39 * (t134 * t240 - t135 * t241) - g(1) * t135 - g(2) * t134 - g(3) * t248) * m(6);];
tau = t1;

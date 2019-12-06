% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:21:59
% EndTime: 2019-12-05 18:22:18
% DurationCPUTime: 12.43s
% Computational Cost: add. (12262->426), mult. (9083->506), div. (0->0), fcn. (6938->8), ass. (0->278)
t562 = Icges(5,3) + Icges(6,3);
t254 = qJ(1) + qJ(2);
t246 = pkin(8) + t254;
t242 = sin(t246);
t243 = cos(t246);
t258 = cos(qJ(4));
t408 = t243 * t258;
t256 = sin(qJ(4));
t409 = t243 * t256;
t130 = Icges(6,4) * t408 - Icges(6,2) * t409 + Icges(6,6) * t242;
t132 = Icges(5,4) * t408 - Icges(5,2) * t409 + Icges(5,6) * t242;
t552 = t130 + t132;
t209 = Icges(6,4) * t409;
t134 = Icges(6,1) * t408 + Icges(6,5) * t242 - t209;
t211 = Icges(5,4) * t409;
t136 = Icges(5,1) * t408 + Icges(5,5) * t242 - t211;
t550 = t134 + t136;
t564 = Icges(5,5) + Icges(6,5);
t563 = -Icges(5,6) - Icges(6,6);
t432 = Icges(6,4) * t256;
t220 = Icges(6,2) * t258 + t432;
t433 = Icges(5,4) * t256;
t222 = Icges(5,2) * t258 + t433;
t559 = t220 + t222;
t250 = Icges(5,4) * t258;
t487 = Icges(5,1) * t256 + t250;
t249 = Icges(6,4) * t258;
t488 = Icges(6,1) * t256 + t249;
t556 = t487 + t488;
t219 = Icges(5,5) * t258 - Icges(5,6) * t256;
t217 = Icges(6,5) * t258 - Icges(6,6) * t256;
t298 = t217 * t242;
t561 = -t219 * t242 + t562 * t243 - t298;
t326 = -Icges(6,2) * t256 + t249;
t300 = t326 * t242;
t129 = Icges(6,6) * t243 - t300;
t327 = -Icges(5,2) * t256 + t250;
t131 = Icges(5,6) * t243 - t242 * t327;
t560 = t129 + t131;
t412 = t242 * t256;
t208 = Icges(6,4) * t412;
t411 = t242 * t258;
t133 = -Icges(6,1) * t411 + Icges(6,5) * t243 + t208;
t210 = Icges(5,4) * t412;
t135 = -Icges(5,1) * t411 + Icges(5,5) * t243 + t210;
t551 = -t133 - t135;
t225 = Icges(6,1) * t258 - t432;
t227 = Icges(5,1) * t258 - t433;
t558 = t225 + t227;
t557 = t552 * t256 - t550 * t258;
t253 = qJD(1) + qJD(2);
t555 = t253 * t258;
t554 = t562 * t242 + t564 * t408 + t563 * t409;
t216 = Icges(6,5) * t256 + Icges(6,6) * t258;
t218 = Icges(5,5) * t256 + Icges(5,6) * t258;
t549 = t216 + t218;
t545 = t559 * t256 - t556 * t258;
t548 = (t556 * qJD(4) - t564 * t253) * t256 + (t559 * qJD(4) + t563 * t253) * t258;
t547 = (t326 + t327) * qJD(4);
t546 = t558 * qJD(4);
t497 = -t561 * t242 + t551 * t408;
t496 = t561 * t243 + t560 * t412;
t544 = t551 * t258;
t495 = t560 * t256;
t541 = t557 * t243;
t540 = -t558 * t256 * t253 - t327 * t555;
t539 = t551 * t411 + t496;
t507 = t554 * t243 - t550 * t411 + t552 * t412;
t538 = -t409 * t560 - t497;
t526 = t554 * t242 - t541;
t416 = t218 * t243;
t418 = t216 * t243;
t537 = t545 * t242 + t416 + t418;
t143 = t216 * t242;
t145 = t218 * t242;
t536 = -t545 * t243 + t143 + t145;
t534 = t495 + t544;
t533 = t549 * qJD(4) - t562 * t253;
t410 = t243 * t253;
t415 = t219 * t253;
t532 = (-t217 * t410 + t533 * t242 - t415 * t243 + t534 * t253) * t243;
t531 = t545 * t253 + (t217 + t219) * qJD(4);
t215 = rSges(5,2) * t409;
t373 = rSges(5,1) * t408;
t440 = rSges(5,3) * t242;
t315 = -t373 - t440;
t138 = -t215 - t315;
t122 = t253 * t138;
t231 = rSges(5,1) * t256 + rSges(5,2) * t258;
t381 = qJD(4) * t242;
t166 = t231 * t381;
t449 = pkin(7) * t242;
t451 = pkin(3) * t243;
t169 = t449 + t451;
t248 = cos(t254);
t406 = t248 * t253;
t376 = pkin(2) * t406;
t337 = t253 * t169 + t376;
t379 = qJD(4) * t256;
t358 = t242 * t379;
t378 = qJD(4) * t258;
t498 = t242 * t378 + t253 * t409;
t362 = rSges(5,1) * t358 + t498 * rSges(5,2);
t530 = -t337 + t166 - t122 - t362;
t529 = -t546 * t258 + t547 * t256 - t549 * t253 + (t556 * t256 + t559 * t258) * qJD(4);
t525 = t537 * t253;
t522 = -t415 * t242 + (-t298 + t557) * t253 - t533 * t243;
t521 = (t526 * t242 + t538 * t243) * qJD(4);
t520 = (t507 * t242 + t539 * t243) * qJD(4);
t255 = -qJ(5) - pkin(7);
t446 = pkin(7) + t255;
t489 = rSges(6,2) * t412 + t243 * rSges(6,3);
t519 = t446 * t243 - t489;
t238 = t242 * rSges(4,2);
t443 = rSges(4,1) * t243;
t168 = -t238 + t443;
t413 = t242 * t253;
t200 = rSges(4,2) * t413;
t453 = pkin(2) * t248;
t336 = -t443 - t453;
t518 = -t376 - t200 + (-t168 - t336) * t253;
t517 = t536 * t253;
t514 = 0.2e1 * qJD(4);
t513 = t520 + t525;
t512 = t517 + t521;
t511 = -t326 * t410 * t258 - t534 * qJD(4) + t548 * t242 + t540 * t243;
t510 = -qJD(4) * t557 + t540 * t242 - t548 * t243 - t300 * t555;
t509 = t531 * t242 - t529 * t243;
t508 = t529 * t242 + t531 * t243;
t506 = t551 * t256 - t258 * t560;
t505 = t550 * t256 + t552 * t258;
t450 = pkin(4) * t258;
t244 = pkin(3) + t450;
t186 = t243 * t244;
t214 = rSges(6,2) * t409;
t372 = rSges(6,1) * t408;
t314 = rSges(6,3) * t242 + t372;
t400 = t242 * t446 - t186 + t214 - t314 + t451;
t504 = t253 * t400;
t385 = t488 + t326;
t386 = t220 - t225;
t503 = (t256 * t385 + t258 * t386) * t253;
t383 = t487 + t327;
t384 = t222 - t227;
t502 = (t256 * t383 + t258 * t384) * t253;
t247 = sin(t254);
t333 = rSges(3,1) * t247 + rSges(3,2) * t248;
t161 = t333 * t253;
t257 = sin(qJ(1));
t439 = pkin(1) * qJD(1);
t371 = t257 * t439;
t140 = t161 + t371;
t501 = -t554 + t544;
t441 = rSges(6,1) * t258;
t344 = -t244 - t441;
t500 = (-pkin(3) - t344) * t242 + t519;
t356 = t243 * t379;
t499 = -t253 * t411 - t356;
t454 = pkin(2) * t247;
t494 = -t489 + t454;
t490 = -t186 - t453;
t213 = rSges(5,2) * t412;
t387 = t243 * rSges(5,3) + t213;
t237 = qJD(5) * t243;
t390 = pkin(4) * t358 + t237;
t405 = t253 * t255;
t311 = rSges(6,1) * t358 + t498 * rSges(6,2) + t242 * t405 + t390;
t230 = rSges(6,1) * t256 + rSges(6,2) * t258;
t366 = t230 * t381 + t390;
t486 = -t311 - t337 + t366 + t504;
t236 = qJD(5) * t242;
t345 = pkin(4) * t256 + t230;
t330 = qJD(4) * t345;
t485 = t243 * t330 - t236;
t355 = t243 * t378;
t484 = t499 * rSges(6,1) - rSges(6,2) * t355 - t243 * t405 - t244 * t413 + t236;
t442 = rSges(5,1) * t258;
t331 = -rSges(5,2) * t256 + t442;
t188 = t331 * qJD(4);
t483 = -t188 * t243 + t231 * t413;
t392 = -Icges(5,2) * t408 + t136 - t211;
t396 = t243 * t487 + t132;
t466 = t256 * t392 + t258 * t396;
t393 = Icges(5,2) * t411 + t135 + t210;
t397 = -t242 * t487 + t131;
t465 = -t256 * t393 - t258 * t397;
t394 = -Icges(6,2) * t408 + t134 - t209;
t398 = t243 * t488 + t130;
t464 = t256 * t394 + t258 * t398;
t395 = Icges(6,2) * t411 + t133 + t208;
t399 = -t242 * t488 + t129;
t463 = -t256 * t395 - t258 * t399;
t252 = t253 ^ 2;
t458 = -rSges(5,3) - pkin(7);
t457 = pkin(1) * t257;
t259 = cos(qJ(1));
t456 = pkin(1) * t259;
t455 = pkin(1) * qJD(1) ^ 2;
t452 = pkin(2) * t252;
t448 = t243 * pkin(7);
t447 = pkin(3) - t244;
t203 = pkin(3) * t413;
t445 = t203 + (-pkin(4) * t379 - pkin(7) * t253) * t243 + t253 * t489 + t484;
t444 = -t311 + (-t243 * t447 + t314 - t449) * t253;
t346 = -t169 - t453;
t370 = t259 * t439;
t60 = -t370 + t166 + (-t138 + t346) * t253;
t438 = t242 * t60;
t417 = t217 * t253;
t156 = t231 * t242;
t414 = t231 * t243;
t407 = t247 * t253;
t389 = -rSges(4,1) * t413 - rSges(4,2) * t410;
t245 = t257 * t455;
t382 = t247 * t452 + t245;
t380 = qJD(4) * t243;
t377 = pkin(2) * t407;
t374 = t259 * t455;
t364 = t499 * rSges(5,1) - rSges(5,2) * t355;
t359 = t231 * t380;
t352 = -pkin(3) - t442;
t350 = -t381 / 0.2e1;
t348 = -t380 / 0.2e1;
t347 = t380 / 0.2e1;
t177 = rSges(3,1) * t248 - t247 * rSges(3,2);
t341 = t203 - t364;
t338 = -t253 * (-pkin(3) * t242 + t448) + t377;
t162 = -rSges(3,1) * t406 + rSges(3,2) * t407;
t332 = -rSges(4,1) * t242 - rSges(4,2) * t243;
t233 = -rSges(6,2) * t256 + t441;
t317 = rSges(5,1) * t411 - t387;
t316 = t238 + t336;
t141 = -t177 * t253 - t370;
t187 = t233 * qJD(4);
t310 = (t450 * qJD(4) + t187) * qJD(4);
t305 = -t248 * t452 - t374;
t304 = t387 + t448 - t454;
t297 = t371 + t377;
t295 = -t243 * t255 - t494;
t120 = t253 * t317;
t294 = t332 - t454;
t289 = -t252 * t169 + t305;
t286 = t214 - t372 + t490;
t277 = t215 - t373 - t451 - t453;
t276 = t120 + t338 + t359;
t274 = t243 * t138 + t242 * t317;
t273 = pkin(4) * t356 - t484;
t268 = t338 + t485 + (rSges(6,1) * t411 - t242 * t447 + t519) * t253;
t91 = t253 * t387 + t364;
t93 = t253 * t315 + t362;
t267 = (-t93 - t122) * t242 + (t91 + t120) * t243;
t266 = t500 * t242 - t400 * t243;
t265 = (((t496 + t526 + t541) * t243 + ((-t495 + t501) * t243 - t497 + t507 - t538) * t242) * qJD(4) + t525) * t350 + ((((t495 - t554) * t243 + t497 + t507) * t243 + (t501 * t242 + t496 - t539) * t242) * qJD(4) + t512 - t517) * t348 + (t508 + t511) * t347 + (t509 + t510 + t513) * t381 / 0.2e1 + (t547 * t258 + t546 * t256 - t545 * qJD(4) + (-t506 + t537) * t350 + (t505 + t536) * t347) * t253;
t264 = (t444 + t504) * t242 + (t500 * t253 + t445) * t243;
t142 = pkin(7) * t410 - t203;
t23 = t310 * t242 + (-t142 - t445 + t485) * t253 + t382;
t24 = -t310 * t243 + (t242 * t330 + t237 - t444) * t253 + t289;
t49 = t268 + t371;
t50 = -t370 + (t346 + t400) * t253 + t366;
t263 = (t23 * (-rSges(6,3) + t255) + t24 * t344) * t242 + (t50 * t494 - t49 * (-t314 + t490)) * t253;
t47 = t188 * t381 + (-t142 - t91 + t359) * t253 + t382;
t48 = qJD(4) * t483 + t253 * t93 + t289;
t59 = t276 + t371;
t262 = (t352 * t48 + t458 * t47) * t242 + (t60 * (-t213 + t454) - t59 * (-t440 - t449 - t453) + (-t352 * t59 + t458 * t60) * t243) * t253;
t159 = t253 * t332;
t157 = t230 * t243;
t155 = t230 * t242;
t124 = t162 * t253 - t374;
t123 = t161 * t253 + t245;
t118 = -t370 + (-t168 - t453) * t253;
t117 = -t159 + t297;
t101 = t253 * (-rSges(4,1) * t410 + t200) + t305;
t100 = -t253 * t389 + t382;
t65 = qJD(4) * t274 + qJD(3);
t42 = qJD(4) * t266 + qJD(3);
t25 = t267 * qJD(4);
t5 = t264 * qJD(4);
t1 = [t265 + m(3) * (t123 * (-t177 - t456) + t124 * (-t333 - t457) + (t141 - t162 + t370) * t140) + (t23 * (t286 - t456) + t50 * (t273 + t371) + t24 * (t295 - t457) + t263 + (-t50 + t486) * t49) * m(6) + (t47 * (t277 - t456) + t60 * (t341 + t371) + t48 * (t304 - t457) + t262 + (-t60 + t530) * t59) * m(5) + (t100 * (t316 - t456) + t118 * (t297 - t389) + t101 * (t294 - t457) + (-t118 + t518) * t117) * m(4); t265 + (t23 * t286 + t24 * t295 + t263 + (-t268 + t273) * t50 + t486 * t49) * m(6) + (t47 * t277 + t48 * t304 + t262 + (t341 - t276) * t60 + t530 * t59) * m(5) + (t100 * t316 + t101 * t294 + t518 * t117 + (-t389 + t159) * t118) * m(4) + (-t123 * t177 - t124 * t333 - t140 * t162 + t141 * t161 - (t140 * t177 + t141 * t333) * t253) * m(3); m(5) * t25 + m(6) * t5; -(((t383 + t385) * t258 + (-t384 - t386) * t256) * t253 + (((t393 + t395) * t243 + (t392 + t394) * t242) * t258 + ((-t397 - t399) * t243 + (-t396 - t398) * t242) * t256) * qJD(4)) * t253 / 0.2e1 + ((t505 * t253 + t511) * t243 + (t506 * t253 + t510) * t242) * t253 / 0.2e1 + ((-t381 * t418 + t417) * t242 + (-t503 + (t463 * t243 + (t143 - t464) * t242) * qJD(4)) * t243 + (-t381 * t416 + t415) * t242 + (-t502 + (t465 * t243 + (t145 - t466) * t242) * qJD(4)) * t243) * t350 + ((t143 * t380 + t417) * t243 + (t503 + (t464 * t242 + (-t418 - t463) * t243) * qJD(4)) * t242 + (t145 * t380 + t415) * t243 + (t502 + (t466 * t242 + (-t416 - t465) * t243) * qJD(4)) * t242) * t348 + (-((t50 * t233 + (-pkin(4) * t412 - t155) * t42) * t242 + (-t49 * (-t233 - t450) + (-pkin(4) * t409 - t157) * t42) * t243) * qJD(4) + t5 * t266 + t42 * t264 + t50 * (t187 * t242 + t230 * t410) - t49 * (t230 * t413 + (-pkin(4) * t378 - t187) * t243) + (t23 * t242 - t24 * t243) * t345 + (t49 * t155 - t50 * t157) * t253) * m(6) + (t60 * t231 * t410 + t47 * t156 + t188 * t438 + t25 * t274 + t65 * t267 - t48 * t414 - t483 * t59 - (-t156 * t59 + t414 * t60) * t253 - (t65 * (-t156 * t242 - t243 * t414) + (t243 * t59 + t438) * t331) * qJD(4)) * m(5) + (t509 * t253 + (t526 * t410 + (t522 * t242 - t253 * t538 + t532) * t242) * t514) * t242 / 0.2e1 + (t508 * t253 + ((t507 * t253 + t532) * t243 + (t522 * t243 - t253 * t539) * t242) * t514) * t243 / 0.2e1 - (t513 + t520) * t413 / 0.2e1 + (t512 + t521) * t410 / 0.2e1; m(6) * (t23 * t243 + t24 * t242);];
tauc = t1(:);

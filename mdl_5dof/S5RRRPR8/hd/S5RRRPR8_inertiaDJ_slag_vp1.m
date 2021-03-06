% Calculate time derivative of joint inertia matrix for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR8_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR8_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:08
% EndTime: 2019-12-31 21:19:35
% DurationCPUTime: 16.65s
% Computational Cost: add. (20653->789), mult. (28379->1078), div. (0->0), fcn. (26047->8), ass. (0->411)
t325 = qJ(2) + qJ(3);
t316 = cos(t325);
t315 = sin(t325);
t508 = Icges(4,4) * t315;
t277 = Icges(4,2) * t316 + t508;
t498 = Icges(5,6) * t315;
t382 = Icges(5,3) * t316 + t498;
t572 = -t277 - t382;
t507 = Icges(4,4) * t316;
t278 = Icges(4,1) * t315 + t507;
t497 = Icges(5,6) * t316;
t384 = Icges(5,2) * t315 + t497;
t571 = t384 + t278;
t328 = sin(qJ(1));
t331 = cos(qJ(1));
t387 = Icges(4,5) * t316 - Icges(4,6) * t315;
t224 = Icges(4,3) * t328 + t331 * t387;
t391 = Icges(5,4) * t316 - Icges(5,5) * t315;
t233 = Icges(5,1) * t328 - t331 * t391;
t570 = t224 + t233;
t392 = -Icges(4,2) * t315 + t507;
t226 = Icges(4,6) * t328 + t331 * t392;
t383 = -Icges(5,3) * t315 + t497;
t229 = Icges(5,5) * t328 - t331 * t383;
t569 = t226 - t229;
t396 = Icges(4,1) * t316 - t508;
t228 = Icges(4,5) * t328 + t331 * t396;
t385 = Icges(5,2) * t316 - t498;
t231 = Icges(5,4) * t328 - t331 * t385;
t568 = t228 - t231;
t322 = qJD(2) + qJD(3);
t567 = (-t383 - t392) * t322;
t566 = (t385 + t396) * t322;
t275 = Icges(4,5) * t315 + Icges(4,6) * t316;
t390 = Icges(5,4) * t315 + Icges(5,5) * t316;
t560 = t275 - t390;
t558 = t315 * t572 + t316 * t571;
t223 = -Icges(4,3) * t331 + t328 * t387;
t539 = Icges(5,1) * t331 + t328 * t391;
t225 = -Icges(4,6) * t331 + t328 * t392;
t227 = -Icges(4,5) * t331 + t328 * t396;
t376 = t225 * t315 - t227 * t316;
t546 = t376 * t331;
t540 = Icges(5,4) * t331 + t328 * t385;
t541 = Icges(5,5) * t331 + t328 * t383;
t373 = -t315 * t541 + t316 * t540;
t549 = t373 * t331;
t565 = t546 - t549 + (-t223 + t539) * t328;
t564 = t225 + t541;
t563 = t227 + t540;
t562 = t539 * t331;
t351 = t322 * t382;
t352 = t322 * t384;
t374 = t229 * t315 - t231 * t316;
t375 = t226 * t315 - t228 * t316;
t561 = (t374 - t375) * t331 + t570 * t328;
t485 = t278 * t322;
t486 = t277 * t322;
t559 = (-t351 - t486 + t566) * t316 + (-t352 - t485 + t567) * t315 + t560 * qJD(1);
t160 = -qJD(1) * t227 - t331 * t485;
t165 = qJD(1) * t540 + t331 * t352;
t557 = -t322 * t569 + t160 - t165;
t158 = -qJD(1) * t225 - t331 * t486;
t163 = qJD(1) * t541 + t331 * t351;
t556 = -t322 * t568 - t158 + t163;
t456 = qJD(1) * t328;
t432 = t316 * t456;
t477 = t322 * t331;
t555 = t315 * t477 + t432;
t332 = -pkin(7) - pkin(6);
t327 = sin(qJ(2));
t453 = qJD(2) * t327;
t447 = pkin(2) * t453;
t554 = qJD(1) * t332 + t447;
t553 = (-t387 + t391) * t322 + t558 * qJD(1);
t330 = cos(qJ(2));
t296 = rSges(3,1) * t327 + rSges(3,2) * t330;
t358 = qJD(2) * t296;
t552 = t328 * t358;
t509 = Icges(3,4) * t330;
t394 = -Icges(3,2) * t327 + t509;
t253 = Icges(3,6) * t328 + t331 * t394;
t510 = Icges(3,4) * t327;
t398 = Icges(3,1) * t330 - t510;
t255 = Icges(3,5) * t328 + t331 * t398;
t371 = t253 * t327 - t255 * t330;
t551 = t371 * t328;
t252 = -Icges(3,6) * t331 + t328 * t394;
t254 = -Icges(3,5) * t331 + t328 * t398;
t372 = t252 * t327 - t254 * t330;
t550 = t372 * t331;
t548 = t374 * t328;
t547 = t375 * t328;
t310 = pkin(2) * t330 + pkin(1);
t522 = pkin(1) - t310;
t545 = t522 * t328;
t329 = cos(qJ(5));
t473 = t329 * t331;
t326 = sin(qJ(5));
t476 = t326 * t328;
t265 = t315 * t473 - t476;
t474 = t328 * t329;
t475 = t326 * t331;
t266 = t315 * t475 + t474;
t480 = t316 * t331;
t179 = rSges(6,1) * t266 + rSges(6,2) * t265 + rSges(6,3) * t480;
t267 = t315 * t474 + t475;
t268 = t315 * t476 - t473;
t408 = -t268 * rSges(6,1) - t267 * rSges(6,2);
t481 = t316 * t328;
t180 = rSges(6,3) * t481 - t408;
t544 = -t179 * t331 - t180 * t328;
t320 = t328 * pkin(4);
t271 = pkin(8) * t480 + t320;
t543 = -t179 - t271;
t542 = qJD(1) * t223;
t388 = Icges(3,5) * t330 - Icges(3,6) * t327;
t250 = -Icges(3,3) * t331 + t328 * t388;
t536 = 2 * m(3);
t535 = 2 * m(4);
t534 = 2 * m(5);
t533 = 2 * m(6);
t323 = t328 ^ 2;
t324 = t331 ^ 2;
t532 = m(5) / 0.2e1;
t531 = m(6) / 0.2e1;
t530 = t328 / 0.2e1;
t529 = -t331 / 0.2e1;
t528 = rSges(6,3) + pkin(8);
t527 = m(3) * t296;
t281 = rSges(4,1) * t315 + rSges(4,2) * t316;
t526 = m(4) * t281;
t525 = pkin(2) * t327;
t524 = pkin(3) * t316;
t523 = t328 * pkin(6);
t321 = t331 * pkin(6);
t521 = -pkin(6) - t332;
t386 = Icges(6,5) * t326 + Icges(6,6) * t329;
t215 = Icges(6,3) * t315 - t316 * t386;
t505 = Icges(6,4) * t326;
t389 = Icges(6,2) * t329 + t505;
t216 = Icges(6,6) * t315 - t316 * t389;
t504 = Icges(6,4) * t329;
t395 = Icges(6,1) * t326 + t504;
t217 = Icges(6,5) * t315 - t316 * t395;
t492 = t217 * t326;
t100 = t215 * t315 + (-t216 * t329 - t492) * t316;
t484 = t315 * t322;
t138 = t389 * t484 + (Icges(6,6) * t322 + (Icges(6,2) * t326 - t504) * qJD(5)) * t316;
t137 = t386 * t484 + (Icges(6,3) * t322 + (-Icges(6,5) * t329 + Icges(6,6) * t326) * qJD(5)) * t316;
t450 = qJD(5) * t316;
t478 = t322 * t329;
t482 = t316 * t322;
t399 = t315 * t137 + t215 * t482 + t484 * t492 + (t315 * t478 + t326 * t450) * t216;
t139 = t395 * t484 + (Icges(6,5) * t322 + (-Icges(6,1) * t329 + t505) * qJD(5)) * t316;
t493 = t139 * t326;
t520 = ((-t493 + (-qJD(5) * t217 - t138) * t329) * t316 + t399) * t315 + t100 * t482;
t519 = rSges(3,1) * t330;
t518 = rSges(4,1) * t316;
t517 = rSges(3,2) * t327;
t516 = rSges(5,2) * t315;
t515 = rSges(3,3) * t331;
t173 = Icges(6,5) * t266 + Icges(6,6) * t265 + Icges(6,3) * t480;
t175 = Icges(6,4) * t266 + Icges(6,2) * t265 + Icges(6,6) * t480;
t177 = Icges(6,1) * t266 + Icges(6,4) * t265 + Icges(6,5) * t480;
t381 = t175 * t329 + t177 * t326;
t415 = qJD(1) * t315 + qJD(5);
t440 = t316 * t477;
t339 = -t328 * t415 + t440;
t416 = qJD(5) * t315 + qJD(1);
t366 = t416 * t326;
t153 = t329 * t339 - t331 * t366;
t367 = t329 * t416;
t154 = t326 * t339 + t331 * t367;
t82 = Icges(6,5) * t154 + Icges(6,6) * t153 - Icges(6,3) * t555;
t84 = Icges(6,4) * t154 + Icges(6,2) * t153 - Icges(6,6) * t555;
t86 = Icges(6,1) * t154 + Icges(6,4) * t153 - Icges(6,5) * t555;
t22 = (t322 * t381 + t82) * t315 + (t173 * t322 - t326 * t86 - t329 * t84 + (t175 * t326 - t177 * t329) * qJD(5)) * t316;
t514 = t22 * t328;
t174 = Icges(6,5) * t268 + Icges(6,6) * t267 + Icges(6,3) * t481;
t176 = Icges(6,4) * t268 + Icges(6,2) * t267 + Icges(6,6) * t481;
t178 = Icges(6,1) * t268 + Icges(6,4) * t267 + Icges(6,5) * t481;
t380 = t176 * t329 + t178 * t326;
t365 = t415 * t331;
t151 = t329 * t365 + (t316 * t478 - t366) * t328;
t479 = t322 * t328;
t441 = t316 * t479;
t152 = t328 * t367 + (t365 + t441) * t326;
t455 = qJD(1) * t331;
t431 = t316 * t455;
t443 = t315 * t479;
t347 = t431 - t443;
t81 = Icges(6,5) * t152 + Icges(6,6) * t151 + Icges(6,3) * t347;
t83 = Icges(6,4) * t152 + Icges(6,2) * t151 + Icges(6,6) * t347;
t85 = Icges(6,1) * t152 + Icges(6,4) * t151 + Icges(6,5) * t347;
t23 = (t322 * t380 + t81) * t315 + (t174 * t322 - t326 * t85 - t329 * t83 + (t176 * t326 - t178 * t329) * qJD(5)) * t316;
t513 = t23 * t331;
t319 = t328 * rSges(5,1);
t318 = t328 * rSges(3,3);
t317 = t328 * rSges(4,3);
t512 = -rSges(5,3) - qJ(4);
t494 = qJ(4) * t315;
t410 = -rSges(4,2) * t315 + t518;
t249 = t410 * t322;
t491 = t249 * t328;
t490 = t252 * t330;
t489 = t253 * t330;
t488 = t254 * t327;
t487 = t255 * t327;
t483 = t315 * t331;
t472 = t331 * t332;
t471 = rSges(6,1) * t154 + rSges(6,2) * t153;
t221 = t321 + t472 - t545;
t299 = t331 * t310;
t222 = -pkin(1) * t331 + t328 * t521 + t299;
t470 = t221 * t328 + t222 * t331;
t210 = qJ(4) * t484 + (pkin(3) * t322 - qJD(4)) * t316;
t406 = -rSges(5,2) * t316 + rSges(5,3) * t315;
t469 = -t322 * t406 - t210;
t236 = -rSges(4,3) * t331 + t328 * t410;
t237 = rSges(4,1) * t480 - rSges(4,2) * t483 + t317;
t155 = t236 * t328 + t237 * t331;
t238 = -rSges(5,2) * t480 + rSges(5,3) * t483 + t319;
t261 = pkin(3) * t480 + qJ(4) * t483;
t468 = -t238 - t261;
t260 = (t494 + t524) * t328;
t467 = t260 * t328 + t261 * t331;
t279 = pkin(3) * t315 - qJ(4) * t316;
t263 = t279 * t456;
t405 = rSges(5,3) * t316 + t516;
t466 = -t405 * t456 + t263;
t465 = -t279 + t405;
t451 = qJD(4) * t315;
t464 = qJ(4) * t440 + t331 * t451;
t433 = t315 * t456;
t463 = rSges(4,2) * t433 + rSges(4,3) * t455;
t462 = t554 * t328;
t461 = t331 * t519 + t318;
t460 = t323 + t324;
t459 = qJD(1) * t224;
t458 = qJD(1) * t233;
t251 = Icges(3,3) * t328 + t331 * t388;
t457 = qJD(1) * t251;
t452 = qJD(2) * t330;
t449 = -pkin(3) - t528;
t448 = t331 * t517;
t446 = pkin(2) * t452;
t74 = t173 * t315 - t316 * t381;
t93 = t215 * t480 + t216 * t265 + t217 * t266;
t445 = -t93 / 0.2e1 - t74 / 0.2e1;
t75 = t174 * t315 - t316 * t380;
t94 = t215 * t481 + t216 * t267 + t217 * t268;
t444 = t94 / 0.2e1 + t75 / 0.2e1;
t289 = pkin(3) * t443;
t439 = t328 * (pkin(3) * t431 + t328 * t451 - t289 + (t315 * t455 + t441) * qJ(4)) + t331 * (-pkin(3) * t555 - qJ(4) * t433 + t464) + t260 * t455;
t359 = t281 * t322;
t438 = t328 * (-t328 * t359 + (t331 * t410 + t317) * qJD(1)) + t331 * (-rSges(4,1) * t555 - rSges(4,2) * t440 + t463) + t236 * t455;
t437 = -t261 + t543;
t436 = t328 * ((-t331 * t522 - t523) * qJD(1) - t462) + t331 * (-t331 * t447 + (t331 * t521 + t545) * qJD(1)) + t221 * t455;
t407 = rSges(6,1) * t326 + rSges(6,2) * t329;
t218 = rSges(6,3) * t315 - t316 * t407;
t202 = t218 * t456;
t435 = pkin(8) * t433 + t202 + t263;
t434 = t289 + t462;
t430 = t327 * t456;
t429 = -t484 / 0.2e1;
t428 = t456 / 0.2e1;
t427 = t455 / 0.2e1;
t426 = -t281 - t525;
t199 = t465 * t331;
t159 = qJD(1) * t226 - t328 * t486;
t424 = t227 * t322 + t159;
t161 = qJD(1) * t228 - t328 * t485;
t422 = t225 * t322 - t161;
t162 = qJD(1) * t229 + t328 * t351;
t421 = t322 * t540 - t162;
t164 = qJD(1) * t231 + t328 * t352;
t419 = -t322 * t541 - t164;
t417 = -t328 * t332 + t299;
t239 = -rSges(5,1) * t331 + t328 * t406;
t101 = t238 * t331 + t239 * t328 + t467;
t414 = rSges(5,1) * t455 + rSges(5,2) * t555 + rSges(5,3) * t440;
t413 = t465 - t525;
t412 = -pkin(8) * t315 - t218 - t279;
t411 = -t517 + t519;
t409 = rSges(6,1) * t152 + rSges(6,2) * t151;
t60 = t173 * t480 + t175 * t265 + t177 * t266;
t61 = t174 * t480 + t176 * t265 + t178 * t266;
t404 = t328 * t61 + t331 * t60;
t41 = t328 * t60 - t331 * t61;
t62 = t173 * t481 + t175 * t267 + t177 * t268;
t63 = t174 * t481 + t176 * t267 + t178 * t268;
t403 = t328 * t63 + t331 * t62;
t42 = t328 * t62 - t331 * t63;
t402 = t328 * t75 + t331 * t74;
t401 = t328 * t74 - t331 * t75;
t400 = t315 * t512 - t310;
t397 = Icges(3,1) * t327 + t509;
t393 = Icges(3,2) * t330 + t510;
t379 = t179 * t328 - t180 * t331;
t140 = t407 * t484 + (rSges(6,3) * t322 + (-rSges(6,1) * t329 + rSges(6,2) * t326) * qJD(5)) * t316;
t368 = -pkin(8) * t482 - t140 - t210;
t364 = -t446 + t469;
t363 = -pkin(1) - t411;
t362 = t417 + t261;
t191 = t413 * t331;
t150 = t412 * t331;
t361 = -t310 - t410;
t360 = t328 * (t405 * t479 + (t331 * t406 + t319) * qJD(1)) + t331 * (-rSges(5,3) * t433 + t414) + t239 * t455 + t439;
t272 = -t331 * pkin(4) + pkin(8) * t481;
t68 = t271 * t331 + t272 * t328 + t467 - t544;
t357 = t412 - t525;
t354 = t322 * t390;
t353 = t322 * t275;
t350 = qJD(2) * t397;
t349 = qJD(2) * t393;
t348 = qJD(2) * (-Icges(3,5) * t327 - Icges(3,6) * t330);
t106 = -t223 * t331 - t328 * t376;
t107 = -t224 * t331 - t547;
t112 = -t233 * t331 + t548;
t113 = t328 * t373 + t562;
t156 = -t331 * t353 - t542;
t157 = -t328 * t353 + t459;
t166 = t328 * t354 + t458;
t167 = qJD(1) * t539 + t331 * t354;
t18 = t153 * t175 + t154 * t177 - t173 * t555 + t265 * t84 + t266 * t86 + t480 * t82;
t19 = t153 * t176 + t154 * t178 - t174 * t555 + t265 * t83 + t266 * t85 + t480 * t81;
t9 = qJD(1) * t404 + t18 * t328 - t19 * t331;
t345 = t41 * t455 + t42 * t456 + ((-t106 - t113) * t456 + t565 * t455) * t331 + (t9 + ((t156 + t167) * t328 + (t547 - t548 - t565) * qJD(1)) * t328 + (t107 + t112) * t456 + t561 * t455 + ((-t157 - t166) * t328 + (t564 * t482 + t563 * t484 - t542) * t331 + (t557 * t328 + (-t161 + t164) * t331) * t316 + (t556 * t328 + (t159 - t162) * t331) * t315 + ((-t376 + t373 + t570) * t328 + t562 + t561) * qJD(1)) * t331) * t328;
t344 = t316 * t449 - t310 - t494;
t132 = t357 * t331;
t343 = (rSges(5,2) - pkin(3)) * t316 + t400;
t314 = pkin(4) * t455;
t87 = rSges(6,3) * t347 + t409;
t88 = -rSges(6,3) * t555 + t471;
t342 = t439 + (t180 + t272) * t455 + (-pkin(8) * t555 + t314 + t88) * t331 + (pkin(4) * t456 + pkin(8) * t347 + t87) * t328;
t341 = t368 - t446;
t340 = t344 * t328;
t27 = t315 * t93 + t316 * t404;
t28 = t315 * t94 + t316 * t403;
t16 = t151 * t175 + t152 * t177 + t173 * t347 + t267 * t84 + t268 * t86 + t481 * t82;
t17 = t151 * t176 + t152 * t178 + t174 * t347 + t267 * t83 + t268 * t85 + t481 * t81;
t33 = t137 * t481 + t138 * t267 + t139 * t268 + t151 * t216 + t152 * t217 + t215 * t347;
t3 = (-t322 * t403 + t33) * t315 + (-qJD(1) * t42 + t16 * t331 + t17 * t328 + t322 * t94) * t316;
t34 = t137 * t480 + t138 * t265 + t139 * t266 + t153 * t216 + t154 * t217 - t215 * t555;
t4 = (-t322 * t404 + t34) * t315 + (-qJD(1) * t41 + t18 * t331 + t19 * t328 + t322 * t93) * t316;
t8 = qJD(1) * t403 + t16 * t328 - t17 * t331;
t338 = t3 * t529 + t4 * t530 + t8 * t481 / 0.2e1 + t315 * (qJD(1) * t402 - t513 + t514) / 0.2e1 + t28 * t428 + t27 * t427 + t401 * t482 / 0.2e1 + t9 * t480 / 0.2e1 + (t316 * t427 + t328 * t429) * t42 + (t331 * t429 - t432 / 0.2e1) * t41;
t337 = rSges(3,2) * t430 + rSges(3,3) * t455 - t331 * t358;
t13 = (t331 * t157 + (t107 + t546) * qJD(1)) * t331 + (t106 * qJD(1) + (-t158 * t315 + t160 * t316 - t226 * t482 - t228 * t484 + t459) * t328 + (-t156 + t422 * t316 + t424 * t315 + (-t223 - t375) * qJD(1)) * t331) * t328;
t14 = (t331 * t166 + (t112 - t549) * qJD(1)) * t331 + (t113 * qJD(1) + (t163 * t315 - t165 * t316 + t229 * t482 + t231 * t484 + t458) * t328 + (-t167 - t419 * t316 + t421 * t315 + (t539 + t374) * qJD(1)) * t331) * t328;
t336 = (-t8 - t13 - t14) * t331 + t345;
t333 = t514 / 0.2e1 - t513 / 0.2e1 + (t315 * t557 - t316 * t556 - t328 * t553 + t331 * t559 + t34) * t530 + (t33 + t553 * t331 + t559 * t328 + (t421 + t424) * t316 + (t419 - t422) * t315) * t529 + (t315 * t563 + t316 * t564 + t558 * t328 - t560 * t331 + t75 + t94) * t428 + (t315 * t568 + t316 * t569 + t560 * t328 + t558 * t331 + t74 + t93) * t427;
t307 = pkin(2) * t430;
t288 = t411 * qJD(2);
t257 = -t448 + t461;
t256 = t328 * t411 - t515;
t220 = t426 * t331;
t219 = t426 * t328;
t207 = t523 + (pkin(1) - t517) * t331 + t461;
t206 = t328 * t363 + t321 + t515;
t198 = t465 * t328;
t193 = t237 + t417;
t192 = (rSges(4,3) - t332) * t331 + t361 * t328;
t190 = t413 * t328;
t185 = t328 * t348 + t457;
t184 = -qJD(1) * t250 + t331 * t348;
t170 = t552 + ((-rSges(3,3) - pkin(6)) * t328 + t363 * t331) * qJD(1);
t169 = (t321 + (-pkin(1) - t519) * t328) * qJD(1) + t337;
t149 = t412 * t328;
t142 = t238 + t362;
t141 = (rSges(5,1) - t332) * t331 + t343 * t328;
t134 = -t281 * t455 - t491 + (-t327 * t455 - t328 * t452) * pkin(2);
t133 = t281 * t456 + t307 + (-t249 - t446) * t331;
t131 = t357 * t328;
t121 = t251 * t328 - t331 * t371;
t120 = t250 * t328 - t550;
t119 = -t251 * t331 - t551;
t118 = -t250 * t331 - t328 * t372;
t117 = t281 * t479 + (t331 * t361 - t317) * qJD(1) + t462;
t116 = (-t310 - t518) * t456 + (-t359 - t554) * t331 + t463;
t115 = t179 * t315 - t218 * t480;
t114 = -t180 * t315 + t218 * t481;
t105 = qJD(1) * t199 + t328 * t469;
t104 = t331 * t469 + t466;
t103 = t362 - t543;
t102 = (pkin(4) - t332) * t331 + t340 + t408;
t98 = qJD(1) * t191 + t328 * t364;
t97 = t331 * t364 + t307 + t466;
t96 = t155 + t470;
t95 = t379 * t316;
t90 = (-t451 + (t316 * t512 - t516) * t322) * t328 + (t331 * t343 - t319) * qJD(1) + t434;
t89 = (-pkin(3) * t484 - t447) * t331 + (-t472 + (t400 - t524) * t328) * qJD(1) + t414 + t464;
t78 = t101 + t470;
t77 = qJD(1) * t150 + t328 * t368;
t76 = t331 * t368 + t435;
t73 = -t237 * t456 + t438;
t70 = qJD(1) * t132 + t328 * t341;
t69 = t331 * t341 + t307 + t435;
t51 = t68 + t470;
t50 = (-qJ(4) * t482 + (t322 * t528 - qJD(4)) * t315) * t328 + (t331 * t344 - t320) * qJD(1) - t409 + t434;
t49 = t314 + (t449 * t484 - t447) * t331 + (t340 - t472) * qJD(1) + t464 + t471;
t48 = (-t222 - t237) * t456 + t436 + t438;
t47 = (-t218 * t479 - t87) * t315 + (t140 * t328 - t180 * t322 + t218 * t455) * t316;
t46 = (t218 * t477 + t88) * t315 + (-t140 * t331 + t179 * t322 + t202) * t316;
t45 = t456 * t468 + t360;
t30 = (-t222 + t468) * t456 + t360 + t436;
t29 = t379 * t484 + (qJD(1) * t544 - t328 * t88 + t331 * t87) * t316;
t24 = t437 * t456 + t342;
t15 = (-t222 + t437) * t456 + t342 + t436;
t1 = [t399 + (t169 * t207 + t170 * t206) * t536 + (t116 * t193 + t117 * t192) * t535 + (t141 * t90 + t142 * t89) * t534 + (t102 * t50 + t103 * t49) * t533 - t329 * t217 * t450 + t572 * t484 + t571 * t482 + (-t393 + t398) * t453 + (t394 + t397) * t452 + t566 * t315 + (-t138 * t329 - t493 - t567) * t316; ((t489 / 0.2e1 + t487 / 0.2e1 - t207 * t527) * t331 + (t206 * t527 + t490 / 0.2e1 + t488 / 0.2e1) * t328) * qJD(1) + (-t372 * qJD(2) + (qJD(1) * t253 - t328 * t349) * t330 + (qJD(1) * t255 - t328 * t350) * t327) * t529 + t333 + (-t371 * qJD(2) + (-qJD(1) * t252 - t331 * t349) * t330 + (-qJD(1) * t254 - t331 * t350) * t327) * t530 + (t323 / 0.2e1 + t324 / 0.2e1) * t388 * qJD(2) + m(3) * ((-t169 * t328 - t170 * t331) * t296 + (-t206 * t331 - t207 * t328) * t288) + m(4) * (t116 * t219 + t117 * t220 + t133 * t192 + t134 * t193) + m(5) * (t141 * t97 + t142 * t98 + t190 * t89 + t191 * t90) + m(6) * (t102 * t69 + t103 * t70 + t131 * t49 + t132 * t50); t345 + ((t256 * t328 + t257 * t331) * ((qJD(1) * t256 + t337) * t331 + (-t552 + (-t257 - t448 + t318) * qJD(1)) * t328) + t460 * t296 * t288) * t536 - t331 * ((t331 * t185 + (t119 + t550) * qJD(1)) * t331 + (t118 * qJD(1) + (-t253 * t452 - t255 * t453 + t457) * t328 + (-t184 + (t488 + t490) * qJD(2) - t371 * qJD(1)) * t331) * t328) + t328 * ((t328 * t184 + (t120 + t551) * qJD(1)) * t328 + (t121 * qJD(1) + (t252 * t452 + t254 * t453) * t331 + (-t185 + (-t487 - t489) * qJD(2) + (t251 - t372) * qJD(1)) * t328) * t331) + (-t120 * t331 + t121 * t328) * t455 + (-t118 * t331 + t119 * t328) * t456 - t331 * t8 - t331 * t13 - t331 * t14 + (t131 * t70 + t132 * t69 + t15 * t51) * t533 + (t190 * t98 + t191 * t97 + t30 * t78) * t534 + (t133 * t220 + t134 * t219 + t48 * t96) * t535; t333 + m(4) * (-t192 * t331 - t193 * t328) * t249 + (-t116 * t328 - t117 * t331 + (t192 * t328 - t193 * t331) * qJD(1)) * t526 + m(5) * (t104 * t141 + t105 * t142 + t198 * t89 + t199 * t90) + m(6) * (t102 * t76 + t103 * t77 + t149 * t49 + t150 * t50); m(6) * (t131 * t77 + t132 * t76 + t149 * t70 + t15 * t68 + t150 * t69 + t24 * t51) + m(5) * (t101 * t30 + t104 * t191 + t105 * t190 + t198 * t98 + t199 * t97 + t45 * t78) + t336 + (-t133 * t331 - t134 * t328 + (-t219 * t331 + t220 * t328) * qJD(1)) * t526 + m(4) * (-t220 * t249 * t331 + t155 * t48 - t219 * t491 + t73 * t96); (t149 * t77 + t150 * t76 + t24 * t68) * t533 + (t101 * t45 + t104 * t199 + t105 * t198) * t534 + (t249 * t281 * t460 + t155 * t73) * t535 + t336; 0.2e1 * ((t102 * t331 + t103 * t328) * t531 + (t141 * t331 + t142 * t328) * t532) * t482 + 0.2e1 * ((-t102 * t456 + t103 * t455 + t328 * t49 + t331 * t50) * t531 + (-t141 * t456 + t142 * t455 + t328 * t89 + t331 * t90) * t532) * t315; 0.2e1 * ((t131 * t479 + t132 * t477 - t15) * t531 + (t190 * t479 + t191 * t477 - t30) * t532) * t316 + 0.2e1 * ((t131 * t455 - t132 * t456 + t322 * t51 + t328 * t70 + t331 * t69) * t531 + (t190 * t455 - t191 * t456 + t322 * t78 + t328 * t98 + t331 * t97) * t532) * t315; 0.2e1 * ((t149 * t479 + t150 * t477 - t24) * t531 + (t198 * t479 + t199 * t477 - t45) * t532) * t316 + 0.2e1 * ((t149 * t455 - t150 * t456 + t322 * t68 + t328 * t77 + t331 * t76) * t531 + (t101 * t322 + t104 * t331 + t105 * t328 + t198 * t455 - t199 * t456) * t532) * t315; 0.4e1 * (t532 + t531) * (-0.1e1 + t460) * t315 * t482; m(6) * (t102 * t47 + t103 * t46 + t114 * t50 + t115 * t49) + (-t328 * t444 + t331 * t445) * t484 + ((t34 / 0.2e1 + t22 / 0.2e1) * t331 + (t33 / 0.2e1 + t23 / 0.2e1) * t328 + (t328 * t445 + t331 * t444) * qJD(1)) * t316 + t520; m(6) * (t114 * t69 + t115 * t70 + t131 * t46 + t132 * t47 - t15 * t95 + t29 * t51) + t338; m(6) * (t114 * t76 + t115 * t77 + t149 * t46 + t150 * t47 - t24 * t95 + t29 * t68) + t338; m(6) * ((-t29 + (t114 * t331 + t115 * t328) * t322) * t316 + (-t322 * t95 + t328 * t46 + t331 * t47 + (-t114 * t328 + t115 * t331) * qJD(1)) * t315); (t114 * t47 + t115 * t46 - t29 * t95) * t533 + ((-t27 * t331 - t28 * t328 - t315 * t402) * t322 + t520) * t315 + (t331 * t4 + t328 * t3 + t402 * t482 + (t100 * t322 + t22 * t331 + t23 * t328) * t315 + (-t27 * t328 + t28 * t331 - t315 * t401) * qJD(1)) * t316;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

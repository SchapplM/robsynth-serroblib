% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR8_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR8_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:42
% EndTime: 2019-12-31 19:06:01
% DurationCPUTime: 13.27s
% Computational Cost: add. (14103->639), mult. (25632->820), div. (0->0), fcn. (26954->8), ass. (0->339)
t510 = sin(qJ(3));
t511 = sin(qJ(1));
t512 = cos(qJ(3));
t513 = cos(qJ(1));
t252 = -t511 * t510 - t513 * t512;
t253 = t513 * t510 - t511 * t512;
t296 = qJ(4) + qJ(5);
t285 = sin(t296);
t286 = cos(t296);
t368 = Icges(6,5) * t286 - Icges(6,6) * t285;
t130 = -Icges(6,3) * t253 + t252 * t368;
t497 = Icges(6,4) * t285;
t376 = Icges(6,1) * t286 - t497;
t139 = Icges(6,5) * t252 + t253 * t376;
t496 = Icges(6,4) * t286;
t372 = -Icges(6,2) * t285 + t496;
t135 = Icges(6,6) * t252 + t253 * t372;
t482 = t135 * t285;
t569 = -t139 * t286 + t482;
t574 = -t130 + t569;
t131 = Icges(6,3) * t252 + t253 * t368;
t138 = -Icges(6,5) * t253 + t252 * t376;
t134 = -Icges(6,6) * t253 + t252 * t372;
t481 = t134 * t285;
t570 = -t138 * t286 + t481;
t573 = t131 + t570;
t464 = t253 * t286;
t501 = t131 * t252 + t139 * t464;
t572 = t252 * t570 + t501 + (t130 - t482) * t253;
t469 = t252 * t286;
t457 = -t130 * t253 + t138 * t469;
t571 = t253 * t569 - (t131 + t481) * t252 + t457;
t298 = cos(qJ(4));
t297 = sin(qJ(4));
t499 = Icges(5,4) * t297;
t378 = Icges(5,1) * t298 - t499;
t154 = -Icges(5,5) * t253 + t252 * t378;
t498 = Icges(5,4) * t298;
t374 = -Icges(5,2) * t297 + t498;
t150 = -Icges(5,6) * t253 + t252 * t374;
t479 = t150 * t297;
t359 = t154 * t298 - t479;
t155 = Icges(5,5) * t252 + t253 * t378;
t151 = Icges(5,6) * t252 + t253 * t374;
t480 = t151 * t297;
t567 = -t155 * t298 + t480;
t370 = Icges(5,5) * t298 - Icges(5,6) * t297;
t146 = -Icges(5,3) * t253 + t252 * t370;
t147 = Icges(5,3) * t252 + t253 * t370;
t462 = t253 * t298;
t456 = -t147 * t252 - t155 * t462;
t565 = -(t146 - t480) * t253 + t456;
t360 = -t150 * t298 - t154 * t297;
t362 = -t151 * t298 - t155 * t297;
t427 = qJD(4) * t253;
t428 = qJD(4) * t252;
t295 = qJD(1) - qJD(3);
t515 = -t295 / 0.2e1;
t564 = ((-t360 * t252 + t362 * t253) * qJD(4) + t360 * t428 - t362 * t427) * t515;
t560 = t146 * t253;
t559 = t147 * t253;
t552 = t147 + t479;
t288 = t513 * qJ(2);
t421 = t511 * pkin(1);
t265 = t421 - t288;
t283 = qJD(2) * t511;
t410 = qJD(1) * t513;
t431 = qJ(2) * t410 + t283;
t551 = qJD(1) * t265 - t283 + t431;
t82 = t134 * t286 + t138 * t285;
t81 = t135 * t286 + t139 * t285;
t294 = qJD(4) + qJD(5);
t204 = t252 * t294;
t549 = -t204 / 0.2e1;
t205 = t253 * t294;
t548 = t205 / 0.2e1;
t465 = t253 * t285;
t500 = t130 * t252 + t138 * t464;
t46 = -t134 * t465 + t500;
t547 = t46 * t204;
t470 = t252 * t285;
t502 = -t131 * t253 + t139 * t469;
t47 = -t135 * t470 + t502;
t546 = t47 * t205;
t202 = t295 * t252;
t203 = t295 * t253;
t385 = rSges(5,1) * t298 - rSges(5,2) * t297;
t255 = t385 * qJD(4);
t301 = qJD(1) ^ 2;
t420 = t511 * pkin(2);
t409 = qJD(1) * t511;
t442 = (-pkin(1) * t409 + t283 + t431) * qJD(1);
t348 = -t301 * t420 + t442;
t444 = t203 * pkin(3) + t202 * pkin(7);
t333 = t295 * t444 + t348;
t384 = rSges(5,1) * t297 + rSges(5,2) * t298;
t425 = qJD(4) * t297;
t414 = t252 * t425;
t349 = t203 * t298 + t414;
t340 = t349 * rSges(5,1) + t202 * rSges(5,3);
t424 = qJD(4) * t298;
t413 = t252 * t424;
t350 = -t203 * t297 + t413;
t73 = rSges(5,2) * t350 + t340;
t37 = t295 * t73 + (t202 * t384 + t253 * t255) * qJD(4) + t333;
t116 = -t202 * pkin(3) + t203 * pkin(7);
t292 = t513 * pkin(2);
t284 = qJD(2) * t513;
t429 = t513 * pkin(1) + t511 * qJ(2);
t397 = (-qJD(1) * t429 + 0.2e1 * t284) * qJD(1);
t327 = -t292 * t301 + t397;
t412 = t253 * t425;
t351 = -t202 * t298 + t412;
t341 = t351 * rSges(5,1) + t203 * rSges(5,3);
t477 = t202 * t297;
t352 = t253 * t424 + t477;
t72 = rSges(5,2) * t352 + t341;
t38 = (-t116 - t72) * t295 + (-t203 * t384 + t252 * t255) * qJD(4) + t327;
t436 = -t253 * pkin(3) - t252 * pkin(7);
t193 = t295 * t436;
t331 = t283 + (-t420 - t265) * qJD(1);
t323 = -t193 + t331;
t426 = qJD(4) * t384;
t236 = t252 * rSges(5,3);
t440 = rSges(5,1) * t462 + t236;
t463 = t253 * t297;
t160 = rSges(5,2) * t463 - t440;
t535 = t295 * t160;
t76 = t252 * t426 + t323 - t535;
t234 = t253 * rSges(5,3);
t467 = t252 * t298;
t437 = -rSges(5,1) * t467 + t234;
t468 = t252 * t297;
t161 = rSges(5,2) * t468 + t437;
t415 = t292 + t429;
t330 = t415 * qJD(1) - t284;
t534 = t252 * pkin(3) - t253 * pkin(7);
t77 = t253 * t426 + (t161 - t534) * t295 + t330;
t545 = (t297 * (-t202 * t76 - t203 * t77 + t252 * t37 - t253 * t38) + (t252 * t77 - t253 * t76) * t424) * rSges(5,2);
t112 = t294 * t202;
t383 = rSges(6,1) * t286 - rSges(6,2) * t285;
t219 = t383 * t294;
t382 = rSges(6,1) * t285 + rSges(6,2) * t286;
t458 = t298 * qJD(4) ^ 2;
t226 = pkin(4) * t414;
t281 = pkin(4) * t298 + pkin(3);
t299 = -pkin(8) - pkin(7);
t416 = -t202 * t299 + t203 * t281 + t226;
t55 = t416 - t444;
t460 = t285 * t294;
t353 = t203 * t286 + t252 * t460;
t346 = t353 * rSges(6,1) + t202 * rSges(6,3);
t459 = t286 * t294;
t354 = -t203 * t285 + t252 * t459;
t63 = rSges(6,2) * t354 + t346;
t506 = t55 + t63;
t19 = t112 * t382 + t205 * t219 + t506 * t295 + (t202 * t425 + t253 * t458) * pkin(4) + t333;
t113 = t294 * t203;
t396 = pkin(4) * t412;
t345 = -t202 * t281 - t203 * t299 + t396;
t54 = t345 - t116;
t355 = -t202 * t286 + t253 * t460;
t347 = t355 * rSges(6,1) + t203 * rSges(6,3);
t356 = t202 * t285 + t253 * t459;
t62 = rSges(6,2) * t356 + t347;
t20 = -t113 * t382 + t204 * t219 + (-t203 * t425 + t252 * t458) * pkin(4) + (-t116 - t54 - t62) * t295 + t327;
t439 = t252 * t299 - t253 * t281;
t119 = -t436 + t439;
t235 = t252 * rSges(6,3);
t441 = rSges(6,1) * t464 + t235;
t144 = rSges(6,2) * t465 - t441;
t538 = t204 * t382 + t226;
t529 = -t295 * (t119 + t144) + t538;
t43 = t323 + t529;
t443 = t252 * t281 + t253 * t299;
t120 = t534 - t443;
t233 = t253 * rSges(6,3);
t438 = -rSges(6,1) * t469 + t233;
t145 = rSges(6,2) * t470 + t438;
t453 = t120 + t145;
t537 = t205 * t382 + t396;
t44 = (-t534 + t453) * t295 + t330 + t537;
t381 = t252 * t44 - t253 * t43;
t544 = (t285 * (t19 * t252 - t20 * t253 - t202 * t43 - t203 * t44) + t381 * t459) * rSges(6,2);
t387 = t253 * rSges(4,1) - t252 * rSges(4,2);
t541 = t295 * t387;
t369 = Icges(5,5) * t297 + Icges(5,6) * t298;
t179 = t369 * t252;
t178 = t369 * t253;
t375 = Icges(6,1) * t285 + t496;
t540 = -t372 - t375;
t338 = -t421 - t420;
t539 = pkin(2) * t409 + t338 * qJD(1) + t551;
t536 = 0.2e1 * qJD(4);
t373 = Icges(5,2) * t298 + t499;
t377 = Icges(5,1) * t297 + t498;
t357 = -t297 * t373 + t298 * t377;
t95 = t253 * t357 + t179;
t93 = t95 * t295;
t96 = t252 * t357 - t178;
t94 = t96 * t295;
t336 = t513 * rSges(3,1) + t511 * rSges(3,3);
t533 = t429 + t336;
t528 = t193 + t539;
t446 = t373 * t253 - t155;
t448 = -t377 * t253 - t151;
t526 = t297 * t446 + t298 * t448;
t371 = Icges(6,2) * t286 + t497;
t524 = t205 * (t371 * t252 - t138) - t204 * (t371 * t253 - t139) + t295 * t540;
t521 = t202 / 0.2e1;
t520 = t203 / 0.2e1;
t519 = -t205 / 0.2e1;
t518 = t204 / 0.2e1;
t517 = -t252 / 0.2e1;
t516 = t253 / 0.2e1;
t514 = t295 / 0.2e1;
t509 = -pkin(3) + t281;
t508 = pkin(7) + t299;
t507 = t202 * t144 + t253 * t62;
t503 = t77 * t384;
t483 = t119 * t202;
t367 = Icges(6,5) * t285 + Icges(6,6) * t286;
t472 = t367 * t295;
t461 = t370 * t295;
t455 = t146 * t252 + t154 * t462;
t452 = t155 * t467 - t559;
t451 = t154 * t467 - t560;
t447 = -t377 * t252 - t150;
t445 = t373 * t252 - t154;
t433 = -t373 + t378;
t432 = -t374 - t377;
t423 = -t513 / 0.2e1;
t422 = t511 / 0.2e1;
t419 = t511 * rSges(3,1);
t406 = t428 / 0.2e1;
t405 = -t427 / 0.2e1;
t404 = pkin(4) * t297 + t382;
t403 = -Icges(6,1) * t355 - Icges(6,4) * t356 - Icges(6,5) * t203 - t135 * t294;
t402 = -Icges(6,1) * t353 - Icges(6,4) * t354 - Icges(6,5) * t202 - t134 * t294;
t401 = Icges(6,4) * t355 + Icges(6,2) * t356 + Icges(6,6) * t203 - t139 * t294;
t400 = Icges(6,4) * t353 + Icges(6,2) * t354 + Icges(6,6) * t202 - t138 * t294;
t173 = t382 * t252;
t399 = t295 * t173 + t205 * t383;
t398 = t540 * t294;
t395 = t439 - t441;
t394 = -t438 + t443;
t393 = t436 - t440;
t392 = t534 - t437;
t115 = t203 * rSges(4,1) - t202 * rSges(4,2);
t114 = -t202 * rSges(4,1) - t203 * rSges(4,2);
t386 = rSges(4,1) * t252 + rSges(4,2) * t253;
t379 = -t252 * t76 - t253 * t77;
t358 = -t285 * t371 + t286 * t375;
t49 = -t151 * t463 - t456;
t50 = -t150 * t463 + t455;
t344 = (-t252 * t49 + t253 * t50) * qJD(4);
t51 = -t151 * t468 + t452;
t52 = -t150 * t468 + t451;
t343 = (-t252 * t51 + t253 * t52) * qJD(4);
t339 = t358 * t295;
t78 = (t160 * t253 + t161 * t252) * qJD(4);
t337 = -t419 - t421;
t329 = -t295 * t368 + (-t204 * t253 + t252 * t205) * t367;
t328 = t288 + t338;
t326 = t297 * t445 + t298 * t447;
t325 = t341 + t116;
t324 = -t340 - t444;
t321 = -t346 - t416;
t312 = t285 * t400 + t286 * t402;
t57 = Icges(6,5) * t353 + Icges(6,6) * t354 + Icges(6,3) * t202;
t10 = -t130 * t202 + t203 * t570 + t252 * t312 + t253 * t57;
t45 = -t135 * t465 + t501;
t90 = t252 * t367 + t253 * t358;
t88 = t90 * t295;
t21 = -t204 * t45 + t205 * t46 + t88;
t48 = -t134 * t470 + t457;
t91 = t252 * t358 - t253 * t367;
t89 = t91 * t295;
t22 = -t204 * t47 + t205 * t48 + t89;
t28 = t285 * t403 - t286 * t401;
t29 = t285 * t402 - t286 * t400;
t306 = (-t375 * t252 - t134) * t205 - (-t375 * t253 - t135) * t204 + (-t371 + t376) * t295;
t303 = t524 * t285 + t306 * t286;
t216 = t368 * t294;
t218 = t376 * t294;
t309 = (-t294 * t371 + t218) * t286 + t398 * t285;
t35 = t202 * t358 - t203 * t367 + t216 * t252 + t253 * t309;
t36 = -t202 * t367 - t203 * t358 - t216 * t253 + t252 * t309;
t313 = t285 * t401 + t286 * t403;
t56 = Icges(6,5) * t355 + Icges(6,6) * t356 + Icges(6,3) * t203;
t7 = -t131 * t203 - t202 * t569 - t252 * t56 + t253 * t313;
t8 = -t130 * t203 - t202 * t570 - t252 * t57 + t253 * t312;
t9 = -t131 * t202 + t203 * t569 + t252 * t313 + t253 * t56;
t320 = (t10 * t205 + t112 * t48 + t113 * t47 - t204 * t9 + t295 * t36) * t516 + (t252 * t303 + t253 * t329) * t519 + (-t252 * t329 + t253 * t303) * t518 + t22 * t521 + t21 * t520 + (t112 * t46 + t113 * t45 - t204 * t7 + t205 * t8 + t295 * t35) * t517 + (t306 * t285 - t524 * t286) * t515 + t112 * (-t252 * t47 + t253 * t48) / 0.2e1 + t113 * (-t252 * t45 + t253 * t46) / 0.2e1 + (t10 * t253 + t202 * t48 + t203 * t47 - t252 * t9) * t548 + (t202 * t46 + t203 * t45 - t252 * t7 + t253 * t8) * t549 + (t202 * t82 + t203 * t81 - t252 * t28 + t253 * t29) * t514;
t317 = (t297 * t432 + t298 * t433) * t295;
t172 = t382 * t253;
t34 = t144 * t205 + t145 * t204 + (t119 * t253 + t120 * t252) * qJD(4);
t314 = t34 * (t205 * t172 + t173 * t204) + t43 * (-t172 * t295 + t204 * t383);
t68 = Icges(5,4) * t351 + Icges(5,2) * t352 + Icges(5,6) * t203;
t70 = Icges(5,1) * t351 + Icges(5,4) * t352 + Icges(5,5) * t203;
t308 = qJD(4) * t362 + t297 * t68 - t298 * t70;
t69 = Icges(5,4) * t349 + Icges(5,2) * t350 + Icges(5,6) * t202;
t71 = Icges(5,1) * t349 + Icges(5,4) * t350 + Icges(5,5) * t202;
t307 = qJD(4) * t360 + t297 * t69 - t298 * t71;
t250 = t374 * qJD(4);
t251 = t378 * qJD(4);
t305 = -t250 * t297 + t251 * t298 + (-t297 * t377 - t298 * t373) * qJD(4);
t304 = t345 + t347;
t26 = t93 + t344;
t27 = t94 + t343;
t31 = -t567 * qJD(4) - t297 * t70 - t298 * t68;
t32 = qJD(4) * t359 - t297 * t71 - t298 * t69;
t249 = t370 * qJD(4);
t39 = t202 * t357 - t203 * t369 + t249 * t252 + t253 * t305;
t40 = -t202 * t369 - t203 * t357 - t249 * t253 + t252 * t305;
t302 = t26 * t427 / 0.2e1 - (t91 + t82) * t112 / 0.2e1 - (t90 + t81) * t113 / 0.2e1 + (t36 + t29) * t519 + (t35 + t28) * t518 + (t40 + t32) * t405 + (-t218 * t285 + t371 * t460 + t286 * t398 - t251 * t297 + t373 * t425 + (-qJD(4) * t377 - t250) * t298) * t295 + (t27 + t39 + t31) * t406 - ((t96 - t360) * t202 + (t95 - t362) * t203) * qJD(4) / 0.2e1;
t290 = t513 * rSges(3,3);
t282 = rSges(3,3) * t410;
t185 = t384 * t252;
t184 = t384 * t253;
t165 = -t301 * t336 + t397;
t164 = qJD(1) * (-rSges(3,1) * t409 + t282) + t442;
t159 = t253 * t385 + t236;
t158 = t252 * t385 - t234;
t143 = t253 * t383 + t235;
t142 = t252 * t383 - t233;
t129 = -t386 * t295 + t330;
t128 = t331 + t541;
t118 = -t252 * t508 + t253 * t509;
t117 = t252 * t509 + t253 * t508;
t107 = t253 * t144;
t86 = -t295 * t114 + t327;
t85 = t295 * t115 + t348;
t67 = Icges(5,5) * t349 + Icges(5,6) * t350 + Icges(5,3) * t202;
t66 = Icges(5,5) * t351 + Icges(5,6) * t352 + Icges(5,3) * t203;
t33 = (-t150 * t253 - t151 * t252) * t297 + t452 + t455;
t30 = (-t134 * t253 - t135 * t252) * t285 + t500 + t502;
t6 = t112 * t144 - t113 * t145 + t205 * t62 + t204 * t63 + (-t120 * t203 + t252 * t55 + t253 * t54 + t483) * qJD(4);
t1 = [-t302 + (-t93 + ((t51 + (t147 - t359) * t253) * t253 + (t552 * t252 - t451 + t52 - t560) * t252) * qJD(4)) * t405 + (t94 + ((t253 * t567 + t451 + t49) * t253 + (-t552 * t253 - t33 + t50) * t252) * qJD(4)) * t406 + (t547 - t30 * t204 + t89 + (t45 + t571) * t205) * t518 + t22 * t549 + (t21 + t546 + (t574 * t205 - t472) * t252 + (t573 * t205 - t339) * t253 + (t48 - t571) * t204) * t519 + t564 + (-t34 * ((t118 + t119) * t428 + (t143 + t144) * t204) + t20 * (t328 - t395) + t19 * (-t394 + t415) + t544 + (-t304 - t330) * t43 + (t43 - (t118 + t143) * t295 - t321 + t528 - t538) * t44) * m(6) + (-(t503 + t78 * (t159 + t160)) * t428 + t38 * (t328 - t393) + t37 * (-t392 + t415) + t545 + (-t325 - t330) * t76 + (-t295 * t159 - t324 + t528 + t76) * t77) * m(5) + (t86 * (t328 + t387) + t85 * (-t386 + t415) + (-t114 - t330) * t128 + (t115 + t128 + t539 - t541) * t129) * m(4) + (t165 * (t288 + t290 + t337) + t164 * t533 + (t282 + (t419 - t290 + t337) * qJD(1) + t551) * (t533 * qJD(1) - t284)) * m(3); 0.2e1 * (t19 * t423 + t20 * t422) * m(6) + 0.2e1 * (t37 * t423 + t38 * t422) * m(5) + 0.2e1 * (t422 * t86 + t423 * t85) * m(4) + 0.2e1 * (t164 * t423 + t165 * t422) * m(3); (-t94 + ((-t49 - t565) * t253 + (-t50 + (t146 - t567) * t252 - t559) * t252) * qJD(4)) * t406 + (t93 + ((t33 - t51) * t253 + (t252 * t359 - t52 + t565) * t252) * qJD(4)) * t405 + t302 + t21 * t548 + (t30 * t205 - t546 + t88 - (t48 + t572) * t204) * t519 + (t22 - t547 + (-t573 * t204 + t472) * t253 + (-t574 * t204 - t339) * t252 + (-t45 + t572) * t205) * t518 + t564 + (-t34 * ((t117 + t120) * t427 + (t142 + t145) * t205) + t20 * t395 + t19 * t394 - t544 + (-t193 + t321 + t529) * t44 + (-(-t117 - t142 - t534) * t295 + t304 - t537) * t43) * m(6) + (-t76 * (-t158 - t534) * t295 - (-t252 * t503 + (t76 * t384 + t78 * (t158 + t161)) * t253) * qJD(4) + t38 * t393 + t76 * t325 + t37 * t392 - t545 + (-t193 + t324 - t535) * t77) * m(5) + (t114 * t128 - t115 * t129 - t86 * t387 + t85 * t386 - (-t128 * t386 - t387 * t129) * t295) * m(4); (-t202 * t360 - t203 * t362 - t252 * t31 + t253 * t32) * t514 + ((t178 * t428 + t461) * t252 + (t317 + (t326 * t253 + (-t179 - t526) * t252) * qJD(4)) * t253) * t406 + ((t179 * t427 - t461) * t253 + (t317 + (-t526 * t252 + (-t178 + t326) * t253) * qJD(4)) * t252) * t405 + ((t297 * t433 - t298 * t432) * t295 + ((t252 * t446 - t253 * t445) * t298 + (-t252 * t448 + t253 * t447) * t297) * qJD(4)) * t515 + t320 + (t295 * t39 + (-(-t147 * t203 - t202 * t567 - t252 * t66 + t253 * t308) * t252 + (-t146 * t203 + t202 * t359 - t252 * t67 + t253 * t307) * t253 + t202 * t50 + t203 * t49) * t536) * t517 + (t295 * t40 + (-(-t147 * t202 + t203 * t567 + t252 * t308 + t253 * t66) * t252 + (-t146 * t202 - t203 * t359 + t252 * t307 + t253 * t67) * t253 + t202 * t52 + t203 * t51) * t536) * t516 + (t343 + t27) * t521 + (t344 + t26) * t520 + (t6 * t107 + t34 * (t483 + t507) + (t6 * t119 + t19 * t404 + t34 * t54) * t253 + (-t34 * t453 - t404 * t43) * t203 + (t20 * t404 + t43 * (pkin(4) * t424 + t219) + t6 * t453 + t34 * t506) * t252 - (t43 * t413 + (t381 * t295 + t34 * (t252 ^ 2 + t253 ^ 2) * qJD(4)) * t297) * pkin(4) - t314 + (pkin(4) * t477 + t202 * t382 + t219 * t253 - t399) * t44) * m(6) + (-(-t184 * t76 + t185 * t77) * t295 - (t78 * (t184 * t253 + t185 * t252) - t379 * t385) * qJD(4) + 0.2e1 * t78 * (t160 * t202 - t161 * t203 + t252 * t73 + t253 * t72) - t379 * t255 - (-t202 * t77 + t203 * t76 - t252 * t38 - t253 * t37) * t384) * m(5); t320 + (t6 * (t145 * t252 + t107) + t34 * (-t145 * t203 + t252 * t63 + t507) - (-t252 * t43 - t253 * t44) * t219 - (-t19 * t253 - t20 * t252 - t202 * t44 + t203 * t43) * t382 - t44 * t399 - t314) * m(6);];
tauc = t1(:);

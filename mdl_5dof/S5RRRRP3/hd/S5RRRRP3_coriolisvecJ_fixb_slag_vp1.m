% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:07
% EndTime: 2019-12-31 21:49:23
% DurationCPUTime: 11.58s
% Computational Cost: add. (15929->468), mult. (11518->558), div. (0->0), fcn. (8862->8), ass. (0->281)
t259 = cos(qJ(4));
t257 = sin(qJ(4));
t430 = Icges(5,4) * t257;
t214 = Icges(5,2) * t259 + t430;
t250 = Icges(6,5) * t257;
t332 = Icges(6,3) * t259 - t250;
t526 = t214 + t332;
t429 = Icges(6,5) * t259;
t216 = Icges(6,1) * t257 - t429;
t251 = Icges(5,4) * t259;
t218 = Icges(5,1) * t257 + t251;
t531 = t216 + t218;
t256 = qJ(1) + qJ(2);
t252 = qJ(3) + t256;
t245 = sin(t252);
t211 = Icges(5,5) * t259 - Icges(5,6) * t257;
t246 = cos(t252);
t307 = t211 * t246;
t123 = Icges(5,3) * t245 + t307;
t213 = Icges(6,4) * t259 + Icges(6,6) * t257;
t308 = t213 * t246;
t125 = Icges(6,2) * t245 + t308;
t530 = t123 + t125;
t334 = Icges(6,1) * t259 + t250;
t128 = -Icges(6,4) * t246 + t245 * t334;
t413 = t245 * t257;
t202 = Icges(5,4) * t413;
t412 = t245 * t259;
t130 = Icges(5,1) * t412 - Icges(5,5) * t246 - t202;
t529 = t128 + t130;
t310 = t334 * t246;
t129 = Icges(6,4) * t245 + t310;
t219 = Icges(5,1) * t259 - t430;
t311 = t219 * t246;
t131 = Icges(5,5) * t245 + t311;
t528 = t129 + t131;
t209 = Icges(6,3) * t257 + t429;
t333 = -Icges(5,2) * t257 + t251;
t527 = t209 - t333;
t524 = t219 + t334;
t514 = -t526 * t257 + t531 * t259;
t409 = t246 * t259;
t201 = Icges(6,5) * t409;
t410 = t246 * t257;
t121 = Icges(6,6) * t245 + Icges(6,3) * t410 + t201;
t523 = t121 * t410 + t530 * t245 + t528 * t409;
t124 = -Icges(6,2) * t246 + t213 * t245;
t113 = t245 * t124;
t120 = -Icges(6,6) * t246 + t209 * t245;
t122 = Icges(5,5) * t412 - Icges(5,6) * t413 - Icges(5,3) * t246;
t522 = -t120 * t410 - t245 * t122 - t529 * t409 - t113;
t520 = t527 * qJD(4);
t519 = t524 * qJD(4);
t210 = Icges(5,5) * t257 + Icges(5,6) * t259;
t212 = Icges(6,4) * t257 - Icges(6,6) * t259;
t518 = t210 + t212;
t517 = -t211 - t213;
t255 = qJD(1) + qJD(2);
t247 = qJD(3) + t255;
t516 = (-Icges(5,6) + Icges(6,6)) * t247 + t526 * qJD(4);
t515 = (Icges(6,4) + Icges(5,5)) * t247 - t531 * qJD(4);
t126 = Icges(5,4) * t412 - Icges(5,2) * t413 - Icges(5,6) * t246;
t424 = t126 * t257;
t327 = -t130 * t259 + t424;
t425 = t124 * t246;
t330 = t120 * t257 + t128 * t259;
t467 = t245 * t330;
t49 = -t425 + t467;
t513 = -t122 * t246 - t245 * t327 + t49;
t485 = -t126 * t410 - t522;
t309 = t333 * t246;
t127 = Icges(5,6) * t245 + t309;
t484 = -t127 * t410 + t523;
t416 = t212 * t246;
t419 = t210 * t246;
t512 = t245 * t514 - t416 - t419;
t417 = t212 * t245;
t420 = t210 * t245;
t511 = t246 * t514 + t417 + t420;
t342 = -t121 * t413 + t125 * t246 - t129 * t412;
t102 = t131 * t412;
t348 = t123 * t246 - t102;
t52 = -t127 * t413 - t348;
t510 = -t342 + t52;
t226 = rSges(6,1) * t259 + rSges(6,3) * t257;
t509 = pkin(4) * t259 + qJ(5) * t257 + t226;
t492 = rSges(6,1) + pkin(4);
t414 = t245 * t247;
t508 = t516 * t246 - t527 * t414;
t411 = t246 * t247;
t507 = t209 * t411 + t516 * t245 - t247 * t309;
t506 = t515 * t246 - t524 * t414;
t505 = (-t310 - t311) * t247 - t515 * t245;
t483 = rSges(6,3) + qJ(5);
t504 = (t121 - t127) * t259 - t528 * t257;
t481 = (-t120 + t126) * t259 + t529 * t257;
t503 = t519 * t259 + t520 * t257 + t518 * t247 + (-t257 * t531 - t526 * t259) * qJD(4);
t502 = (Icges(6,2) + Icges(5,3)) * t247 - t518 * qJD(4);
t423 = t127 * t257;
t501 = -t121 * t257 - t259 * t528 + t423;
t500 = t327 - t330;
t170 = rSges(4,1) * t245 + rSges(4,2) * t246;
t154 = t247 * t170;
t258 = sin(qJ(1));
t438 = pkin(1) * qJD(1);
t368 = t258 * t438;
t248 = sin(t256);
t407 = t248 * t255;
t370 = pkin(2) * t407;
t109 = -t368 - t154 - t370;
t499 = t517 * qJD(4) + t514 * t247;
t498 = t511 * t247;
t238 = t246 * rSges(6,2);
t395 = t509 * t245 - t238;
t497 = t395 * t247;
t241 = t246 * pkin(8);
t172 = pkin(3) * t245 - t241;
t190 = pkin(8) * t411;
t496 = t247 * t172 + t190;
t495 = (t484 * t245 - t485 * t246) * qJD(4);
t494 = (t510 * t245 - t513 * t246) * qJD(4);
t493 = t512 * t247;
t491 = t493 + t494;
t490 = t495 + t498;
t489 = t500 * qJD(4) + t505 * t257 + t507 * t259;
t488 = -t501 * qJD(4) + t506 * t257 - t508 * t259;
t487 = -t499 * t245 + t503 * t246;
t486 = t503 * t245 + t499 * t246;
t382 = t492 * t257 - t483 * t259;
t349 = t246 * t382;
t479 = t245 * rSges(6,2) + pkin(4) * t409;
t478 = t504 * qJD(4) + t530 * t247 + t508 * t257 + t506 * t259;
t477 = t505 * t259 - t507 * t257 + (-t122 - t124) * t247 + t481 * qJD(4);
t387 = rSges(5,2) * t413 + t246 * rSges(5,3);
t133 = rSges(5,1) * t412 - t387;
t116 = t247 * t133;
t376 = qJD(4) * t259;
t359 = t246 * t376;
t408 = t247 * t257;
t367 = t245 * t408;
t312 = rSges(5,3) * t411 + (-t359 + t367) * rSges(5,2);
t377 = qJD(4) * t257;
t360 = t246 * t377;
t476 = -rSges(5,1) * t360 + t116 + t312 + t496;
t475 = t425 + t523;
t375 = qJD(5) * t257;
t197 = t246 * t375;
t462 = rSges(6,2) * t411 + t483 * t359 + t197;
t474 = t462 + t496 + t497;
t361 = t245 * t376;
t473 = t246 * t408 + t361;
t249 = cos(t256);
t176 = rSges(3,1) * t248 + rSges(3,2) * t249;
t421 = t176 * t255;
t140 = -t368 - t421;
t472 = (t307 + t308 + t500) * t247 + t502 * t245;
t471 = t502 * t246 + t501 * t247 + t517 * t414;
t470 = 0.2e1 * qJD(4);
t469 = t245 / 0.2e1;
t374 = qJD(5) * t259;
t394 = rSges(6,1) * t409 + t483 * t410 + t479;
t48 = -t374 + (t395 * t245 + t394 * t246) * qJD(4);
t468 = qJD(4) * t48;
t173 = t246 * pkin(3) + t245 * pkin(8);
t465 = t173 + t394;
t362 = t245 * t377;
t388 = t492 * t362;
t464 = rSges(5,1) * t409 + t245 * rSges(5,3);
t139 = t173 * t247;
t260 = cos(qJ(1));
t253 = t260 * pkin(1);
t261 = qJD(1) ^ 2;
t372 = t261 * t253;
t442 = pkin(2) * t255 ^ 2;
t313 = -t249 * t442 - t372;
t391 = -qJD(4) * t509 + t374;
t341 = t374 + t391;
t378 = qJD(4) * t246;
t358 = t245 * t375;
t440 = rSges(6,3) * t361 + t358 + t473 * qJ(5) - t388 + (t226 * t246 + t479) * t247;
t19 = t341 * t378 + (-t139 + (qJD(4) * t382 - t375) * t245 - t440) * t247 + t313;
t320 = -qJD(4) * t349 + t197;
t290 = t320 - t370;
t272 = t290 - t368;
t46 = (-t172 - t395) * t247 + t272;
t436 = t247 * t46;
t379 = qJD(4) * t245;
t321 = -t382 * t379 + t358;
t406 = t249 * t255;
t371 = pkin(2) * t406;
t291 = t321 + t371;
t369 = t260 * t438;
t47 = t247 * t465 + t291 + t369;
t262 = (-t47 * t492 * t377 + (-t257 * t483 - t259 * t492 - pkin(3)) * t436) * t246 + (-t19 * pkin(3) + (-t46 * qJD(5) - t19 * t483) * t257 + (-qJD(4) * t46 * t483 - t19 * t492) * t259 + (t46 * (-rSges(6,2) - pkin(8)) + t47 * (-pkin(3) - t509)) * t247) * t245;
t463 = t436 * t465 + t262;
t460 = -t370 + t476;
t223 = rSges(5,1) * t257 + rSges(5,2) * t259;
t168 = t223 * t379;
t135 = -rSges(5,2) * t410 + t464;
t300 = t135 + t173;
t336 = -t247 * t300 + t168;
t453 = -t370 + t474;
t439 = rSges(5,1) * t259;
t355 = -pkin(3) - t439;
t363 = t223 * t378;
t303 = -t363 - t370;
t280 = t303 - t368;
t57 = (-t133 - t172) * t247 + t280;
t437 = t246 * t57;
t305 = t369 + t371;
t58 = t305 - t336;
t264 = (t355 * t437 + (t57 * (-rSges(5,3) - pkin(8)) + t58 * t355) * t245) * t247;
t364 = rSges(5,1) * t362 + t473 * rSges(5,2);
t452 = t264 + (-t336 + t364) * t57;
t397 = -Icges(5,2) * t412 + t130 - t202;
t401 = t218 * t245 + t126;
t451 = -t257 * t397 - t259 * t401;
t399 = -t332 * t245 + t128;
t403 = -t216 * t245 + t120;
t450 = -t257 * t399 + t259 * t403;
t446 = t247 / 0.2e1;
t444 = pkin(1) * t258;
t443 = pkin(2) * t248;
t301 = -t247 * t412 - t360;
t441 = t492 * t301 - t483 * t367 + t462;
t171 = t246 * rSges(4,1) - rSges(4,2) * t245;
t422 = t171 * t247;
t418 = t211 * t247;
t415 = t213 * t247;
t402 = -Icges(6,1) * t410 + t121 + t201;
t400 = -t218 * t246 - t127;
t398 = -t332 * t246 + t129;
t396 = -t214 * t246 + t131;
t393 = t382 * t245;
t386 = -t332 + t334;
t385 = t209 - t216;
t384 = -t214 + t219;
t383 = t218 + t333;
t380 = t238 + t241;
t373 = t261 * t444;
t354 = -t379 / 0.2e1;
t351 = t378 / 0.2e1;
t177 = t249 * rSges(3,1) - rSges(3,2) * t248;
t347 = -t122 + t423;
t344 = t380 - t443;
t165 = rSges(3,1) * t406 - rSges(3,2) * t407;
t138 = rSges(4,1) * t411 - rSges(4,2) * t414;
t244 = pkin(2) * t249;
t343 = t171 + t244;
t339 = -rSges(5,2) * t257 + t439;
t338 = -t245 * t58 - t437;
t314 = -t248 * t442 - t373;
t63 = (t133 * t245 + t135 * t246) * qJD(4);
t302 = -t170 - t443;
t299 = t247 * (-pkin(3) * t414 + t190) + t314;
t298 = -t138 - t371;
t297 = t244 + t465;
t296 = -t257 * t398 + t259 * t402;
t295 = -t257 * t396 + t259 * t400;
t294 = t245 * t355 + t241 + t387;
t292 = t244 + t300;
t289 = (t257 * t385 + t259 * t386) * t247;
t288 = (-t257 * t383 + t259 * t384) * t247;
t279 = t294 - t443;
t263 = (((t52 - t102 + (t123 + t424) * t246 + t522) * t246 + (t49 - t467 + t475) * t245) * qJD(4) + t498) * t351 + (t514 * qJD(4) + t519 * t257 - t520 * t259) * t247 + (t487 + t488) * t379 / 0.2e1 + (((t246 * t347 - t475 + t484) * t246 + (t245 * t347 - t113 + t342 + t348 + t485) * t245) * qJD(4) + t491 - t493) * t354 - (t486 - t489 + t490) * t378 / 0.2e1 + ((t481 + t512) * t245 + (-t504 + t511) * t246) * qJD(4) * t446;
t188 = t339 * qJD(4);
t161 = t223 * t246;
t157 = t223 * t245;
t141 = t177 * t255 + t369;
t118 = -t165 * t255 - t372;
t117 = -t255 * t421 - t373;
t110 = t305 + t422;
t99 = -t138 * t247 + t313;
t98 = -t154 * t247 + t314;
t85 = t464 * t247 - t364;
t83 = rSges(5,1) * t301 + t312;
t42 = -t188 * t378 + (-t139 - t85 + t168) * t247 + t313;
t41 = t247 * t83 + (-t188 * t245 - t223 * t411) * qJD(4) + t299;
t18 = (t197 + t441) * t247 + (t245 * t341 - t247 * t349) * qJD(4) + t299;
t5 = (t375 + (t441 + t497) * t246 + (-t394 * t247 + t440) * t245) * qJD(4);
t1 = [t263 + m(3) * (t118 * (-t176 - t444) + t117 * (t177 + t253) + (-t165 - t369 + t141) * t140) + (t19 * (t344 - t444) + t46 * (-t305 + t388) + t18 * (t253 + t297) + t262 + (-t368 + t46 - t272 + t453) * t47) * m(6) + (t42 * (t279 - t444) + t57 * (-t305 + t364) + t41 * (t253 + t292) + t264 + (-t280 + t57 - t368 + t460) * t58) * m(5) + m(4) * (t99 * (t302 - t444) + t98 * (t253 + t343) + (t298 - t369 + t110) * t109); t263 + (t18 * t297 + t19 * t344 + (-t290 + t453) * t47 + (-t371 + t388 + t291) * t46 + t463) * m(6) + (t42 * t279 + t41 * t292 + (-t303 + t460) * t58 + t452) * m(5) + (t99 * t302 + t98 * t343 + (t298 + t371 + t422) * t109) * m(4) + (-(-t140 * t177 - t141 * t176) * t255 + t117 * t177 - t118 * t176 - t140 * t165 - t141 * t421) * m(3); t263 + (t18 * t465 + t19 * t380 + (-t320 + t474) * t47 + (t388 + t321) * t46 + t463) * m(6) + (t42 * t294 + t41 * t300 + (t363 + t476) * t58 + t452) * m(5) + (-(-t109 * t171 - t110 * t170) * t247 - t109 * t138 - t110 * t154 - t170 * t99 + t171 * t98) * m(4); -(((t383 - t385) * t259 + (t384 + t386) * t257) * t247 + (((-t397 - t399) * t246 + (t396 + t398) * t245) * t259 + ((t401 - t403) * t246 + (t400 + t402) * t245) * t257) * qJD(4)) * t247 / 0.2e1 + ((-t247 * t504 + t489) * t246 + (t481 * t247 + t488) * t245) * t446 + ((-t379 * t416 + t415) * t245 + (t289 + (-t450 * t246 + (t417 + t296) * t245) * qJD(4)) * t246 + (-t379 * t419 + t418) * t245 + (t288 + (-t451 * t246 + (t420 + t295) * t245) * qJD(4)) * t246) * t354 + ((-t378 * t417 - t415) * t246 + (t289 + (t296 * t245 + (t416 - t450) * t246) * qJD(4)) * t245 + (-t378 * t420 - t418) * t246 + (t288 + (t295 * t245 + (t419 - t451) * t246) * qJD(4)) * t245) * t351 + (-(t48 * t257 + (t245 * t47 + t246 * t46) * t259) * qJD(5) - (-t349 * t47 + t393 * t46) * t247 - ((-t349 * t48 - t46 * t509) * t246 + (-t48 * t393 - t47 * t509) * t245) * qJD(4) + (-t19 * t382 + t46 * t391 + t5 * t394 + t48 * t441 + (-t47 * t382 + t48 * t395) * t247) * t246 + (-t18 * t382 + t47 * t391 + t5 * t395 + t48 * t440 + (t382 * t46 - t394 * t48) * t247) * t245) * m(6) + (0.2e1 * t63 * ((t83 + t116) * t246 + (-t135 * t247 + t85) * t245) + t338 * t188 + ((-t247 * t58 - t42) * t246 + (t247 * t57 - t41) * t245) * t223 - (t157 * t57 - t161 * t58) * t247 - (t63 * (-t157 * t245 - t161 * t246) + t338 * t339) * qJD(4)) * m(5) + (t487 * t247 + ((t477 * t246 + t484 * t247) * t246 + (t471 * t245 + t485 * t247 + (-t472 + t478) * t246) * t245) * t470) * t469 - (t486 * t247 + ((t472 * t246 + t510 * t247) * t246 + (t478 * t245 + t513 * t247 + (-t471 + t477) * t246) * t245) * t470) * t246 / 0.2e1 + (t491 + t494) * t414 / 0.2e1 + (t490 + t495) * t411 / 0.2e1; (-t259 * t5 + 0.2e1 * (t468 / 0.2e1 + t18 * t469 + t19 * t246 / 0.2e1 - (t245 ^ 2 + t246 ^ 2) * t468 / 0.2e1) * t257) * m(6);];
tauc = t1(:);

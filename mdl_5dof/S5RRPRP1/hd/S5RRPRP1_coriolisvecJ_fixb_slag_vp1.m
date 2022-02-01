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
% m [6x1]
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:19:33
% EndTime: 2022-01-20 10:19:52
% DurationCPUTime: 11.74s
% Computational Cost: add. (12262->438), mult. (9083->521), div. (0->0), fcn. (6938->8), ass. (0->278)
t257 = sin(qJ(4));
t259 = cos(qJ(4));
t207 = Icges(6,5) * t259 - Icges(6,6) * t257;
t209 = Icges(5,5) * t259 - Icges(5,6) * t257;
t513 = t207 + t209;
t530 = Icges(5,5) + Icges(6,5);
t529 = Icges(5,6) + Icges(6,6);
t528 = Icges(5,3) + Icges(6,3);
t255 = qJ(1) + qJ(2);
t247 = pkin(8) + t255;
t244 = cos(t247);
t243 = sin(t247);
t407 = t243 * t259;
t408 = t243 * t257;
t126 = Icges(6,4) * t407 - Icges(6,2) * t408 - Icges(6,6) * t244;
t128 = Icges(5,4) * t407 - Icges(5,2) * t408 - Icges(5,6) * t244;
t518 = t126 + t128;
t200 = Icges(6,4) * t408;
t130 = Icges(6,1) * t407 - Icges(6,5) * t244 - t200;
t201 = Icges(5,4) * t408;
t132 = Icges(5,1) * t407 - Icges(5,5) * t244 - t201;
t525 = t130 + t132;
t426 = Icges(6,4) * t257;
t215 = Icges(6,1) * t259 - t426;
t302 = t215 * t244;
t131 = Icges(6,5) * t243 + t302;
t427 = Icges(5,4) * t257;
t217 = Icges(5,1) * t259 - t427;
t303 = t217 * t244;
t133 = Icges(5,5) * t243 + t303;
t516 = t131 + t133;
t210 = Icges(6,2) * t259 + t426;
t212 = Icges(5,2) * t259 + t427;
t527 = t210 + t212;
t250 = Icges(6,4) * t259;
t214 = Icges(6,1) * t257 + t250;
t251 = Icges(5,4) * t259;
t216 = Icges(5,1) * t257 + t251;
t523 = -t214 - t216;
t526 = t513 * t244;
t498 = t528 * t244 - t530 * t407 + t529 * t408;
t497 = t528 * t243 + t526;
t522 = t215 + t217;
t326 = -Icges(6,2) * t257 + t250;
t327 = -Icges(5,2) * t257 + t251;
t521 = t326 + t327;
t520 = t518 * t257;
t519 = t516 * t407;
t300 = t326 * t244;
t127 = Icges(6,6) * t243 + t300;
t301 = t327 * t244;
t129 = Icges(5,6) * t243 + t301;
t517 = t127 + t129;
t509 = t527 * t257 + t523 * t259;
t493 = -t525 * t259 + t520;
t515 = t521 * qJD(4);
t514 = t522 * qJD(4);
t206 = Icges(6,5) * t257 + Icges(6,6) * t259;
t208 = Icges(5,5) * t257 + Icges(5,6) * t259;
t512 = t208 + t206;
t254 = qJD(1) + qJD(2);
t511 = -qJD(4) * t527 + t529 * t254;
t510 = t523 * qJD(4) + t530 * t254;
t508 = -t497 * t244 + t519;
t404 = t244 * t259;
t507 = t498 * t243 - t525 * t404;
t467 = t497 * t243 + t516 * t404;
t506 = -t493 * t243 + t498 * t244;
t479 = -t517 * t408 + t508;
t405 = t244 * t257;
t478 = -t518 * t405 - t507;
t477 = -t517 * t405 + t467;
t412 = t208 * t244;
t415 = t206 * t244;
t505 = -t509 * t243 - t412 - t415;
t413 = t208 * t243;
t416 = t206 * t243;
t504 = -t509 * t244 + t413 + t416;
t503 = t517 * t257;
t476 = t525 * t257 + t518 * t259;
t475 = t516 * t257 + t517 * t259;
t410 = t243 * t254;
t502 = t511 * t244 - t521 * t410;
t501 = (t300 + t301) * t254 + t511 * t243;
t500 = t510 * t244 - t522 * t410;
t499 = (-t302 - t303) * t254 - t510 * t243;
t496 = t514 * t259 - t515 * t257 + t512 * t254 + (t523 * t257 - t259 * t527) * qJD(4);
t495 = -t512 * qJD(4) + t528 * t254;
t494 = -t516 * t259 + t503;
t256 = -qJ(5) - pkin(7);
t219 = t244 * t256;
t439 = pkin(4) * t259;
t246 = pkin(3) + t439;
t492 = -rSges(6,1) * t407 - t243 * t246 - t219;
t491 = t513 * qJD(4) + t509 * t254;
t490 = t504 * t254;
t240 = t244 * pkin(7);
t166 = pkin(3) * t243 - t240;
t460 = -rSges(6,2) * t408 - t244 * rSges(6,3);
t397 = -t166 + t460 - t492;
t489 = t397 * t254;
t488 = (t477 * t243 - t478 * t244) * qJD(4);
t487 = (t479 * t243 - t506 * t244) * qJD(4);
t486 = t505 * t254;
t485 = t486 + t487;
t484 = t488 + t490;
t483 = t493 * qJD(4) + t499 * t257 - t501 * t259;
t482 = -t494 * qJD(4) + t500 * t257 + t502 * t259;
t481 = t491 * t243 + t496 * t244;
t480 = t496 * t243 - t491 * t244;
t315 = rSges(6,1) * t404 + t243 * rSges(6,3);
t474 = -t244 * t246 - t315;
t249 = cos(t255);
t245 = pkin(2) * t249;
t473 = -t244 * rSges(4,1) - t245;
t472 = -t475 * qJD(4) + t497 * t254 - t502 * t257 + t500 * t259;
t471 = t476 * qJD(4) + t498 * t254 + t501 * t257 + t499 * t259;
t248 = sin(t255);
t441 = pkin(2) * t248;
t470 = -t441 + t492;
t469 = t498 + t503;
t376 = qJD(4) * t259;
t402 = t254 * t257;
t468 = t243 * t376 + t244 * t402;
t258 = sin(qJ(1));
t434 = pkin(1) * qJD(1);
t369 = t258 * t434;
t173 = rSges(3,1) * t248 + rSges(3,2) * t249;
t417 = t173 * t254;
t139 = -t369 - t417;
t466 = (t493 + t526) * t254 + t495 * t243;
t465 = t495 * t244 + t494 * t254 - t513 * t410;
t464 = 0.2e1 * qJD(4);
t409 = t243 * t256;
t463 = -rSges(6,2) * t405 - t409 - t474;
t205 = rSges(5,1) * t404;
t459 = t243 * rSges(5,3) + t205;
t137 = -rSges(5,2) * t405 + t459;
t239 = t243 * pkin(7);
t241 = t244 * pkin(3);
t167 = t241 + t239;
t347 = t167 + t245;
t462 = t137 + t347;
t461 = -rSges(4,2) * t243 - t473;
t377 = qJD(4) * t257;
t361 = t243 * t377;
t227 = qJD(5) * t244;
t385 = pkin(4) * t361 + t227;
t313 = rSges(6,1) * t361 + t468 * rSges(6,2) + t254 * t409 + t385;
t384 = rSges(5,2) * t408 + t244 * rSges(5,3);
t135 = rSges(5,1) * t407 - t384;
t119 = t254 * t135;
t161 = t254 * t166;
t406 = t244 * t254;
t195 = pkin(7) * t406;
t358 = t244 * t376;
t368 = t243 * t402;
t304 = rSges(5,3) * t406 + (-t358 + t368) * rSges(5,2);
t359 = t244 * t377;
t458 = -rSges(5,1) * t359 + t119 + t161 + t195 + t304;
t221 = rSges(5,1) * t257 + rSges(5,2) * t259;
t379 = qJD(4) * t243;
t163 = t221 * t379;
t453 = t254 * t462 - t163;
t433 = t259 * rSges(6,2);
t220 = rSges(6,1) * t257 + t433;
t396 = -t167 + t463;
t450 = t254 * (t347 + t396) - t220 * t379 - t385;
t387 = rSges(6,2) * t368 + rSges(6,3) * t406;
t449 = t161 + t387 + t489;
t389 = -Icges(5,2) * t407 + t132 - t201;
t393 = t216 * t243 + t128;
t448 = -t257 * t389 - t259 * t393;
t391 = -Icges(6,2) * t407 + t130 - t200;
t395 = t214 * t243 + t126;
t447 = -t257 * t391 - t259 * t395;
t253 = t254 ^ 2;
t446 = t243 / 0.2e1;
t445 = -t244 / 0.2e1;
t443 = t254 / 0.2e1;
t442 = pkin(1) * t258;
t440 = pkin(2) * t253;
t260 = cos(qJ(1));
t252 = t260 * pkin(1);
t438 = pkin(3) - t246;
t295 = -t254 * t407 - t359;
t226 = qJD(5) * t243;
t317 = -pkin(4) * t359 + t226;
t437 = -t195 + (t243 * t438 - t219) * t254 + t317 + rSges(6,1) * t295 - rSges(6,2) * t358 + t387;
t436 = -t313 + (-t244 * t438 - t239 + t315) * t254;
t435 = rSges(5,1) * t259;
t242 = t249 * rSges(3,1);
t414 = t207 * t254;
t411 = t209 * t254;
t403 = t248 * t254;
t394 = -t214 * t244 - t127;
t392 = -t216 * t244 - t129;
t390 = -t210 * t244 + t131;
t388 = -t212 * t244 + t133;
t383 = -t210 + t215;
t382 = t214 + t326;
t381 = -t212 + t217;
t380 = t216 + t327;
t378 = qJD(4) * t244;
t370 = t260 * t434;
t61 = t370 + t453;
t375 = t61 * t441;
t374 = pkin(2) * t403;
t262 = qJD(1) ^ 2;
t373 = t262 * t442;
t372 = t262 * t252;
t364 = rSges(5,1) * t361 + t468 * rSges(5,2);
t362 = t221 * t378;
t355 = -pkin(3) - t435;
t354 = -t379 / 0.2e1;
t351 = t378 / 0.2e1;
t164 = rSges(4,1) * t243 + rSges(4,2) * t244;
t296 = -t164 - t441;
t348 = -t166 - t441;
t346 = pkin(4) * t257 + t220;
t223 = rSges(6,1) * t259 - rSges(6,2) * t257;
t345 = -t223 - t439;
t174 = -rSges(3,2) * t248 + t242;
t337 = -pkin(4) * t405 - t220 * t244;
t160 = -rSges(3,2) * t403 + t242 * t254;
t335 = t226 - t369;
t184 = t223 * qJD(4);
t334 = -pkin(4) * t376 - t184;
t331 = -rSges(5,2) * t257 + t435;
t297 = -t362 - t369;
t60 = (-t135 + t348) * t254 + t297;
t330 = -t243 * t61 - t244 * t60;
t321 = t135 * t243 + t137 * t244;
t314 = t346 * t378;
t312 = (-t439 * qJD(4) - t184) * qJD(4);
t306 = -t248 * t440 - t373;
t305 = -t249 * t440 - t372;
t294 = t254 * (-pkin(3) * t410 + t195) + t306;
t293 = -t257 * t390 + t259 * t394;
t292 = -t257 * t388 + t259 * t392;
t290 = -t460 + t470;
t289 = (-t257 * t382 + t259 * t383) * t254;
t288 = (-t257 * t380 + t259 * t381) * t254;
t281 = t245 + t463;
t280 = t243 * t355 + t240 + t384 - t441;
t279 = -t314 + t335;
t92 = rSges(5,1) * t295 + t304;
t94 = t459 * t254 - t364;
t278 = (t92 + t119) * t244 + (-t137 * t254 + t94) * t243;
t116 = t254 * t296 - t369;
t117 = t254 * t461 + t370;
t266 = (t116 * t473 + t117 * t296) * t254;
t265 = ((t467 * t243 + ((t497 + t520) * t244 + t479 + t507 - t519) * t244) * qJD(4) + t490) * t351 + (-t509 * qJD(4) + t514 * t257 + t515 * t259) * t254 + (((t469 * t244 - t467 + t477) * t244 + (t469 * t243 + t478 - t508) * t243) * qJD(4) + t485 - t486) * t354 + (t481 + t482) * t379 / 0.2e1 - (t480 - t483 + t484) * t378 / 0.2e1 + ((t476 + t505) * t243 + (t475 + t504) * t244) * qJD(4) * t443;
t264 = (t60 * (-t205 - t241 - t245) - t375 + (t60 * (-rSges(5,3) - pkin(7)) + t61 * t355) * t243) * t254;
t50 = (t348 - t397) * t254 + t279;
t51 = t370 + t450;
t263 = (t50 * (-t245 + t474) + t51 * t470) * t254 + t51 * (-t433 + (-rSges(6,1) - pkin(4)) * t257) * t378;
t192 = rSges(4,2) * t410;
t185 = t331 * qJD(4);
t158 = t254 * t164;
t157 = t221 * t244;
t155 = t221 * t243;
t154 = t220 * t243;
t141 = t167 * t254;
t140 = t174 * t254 + t370;
t121 = -t160 * t254 - t372;
t120 = -t254 * t417 - t373;
t102 = -t254 * (rSges(4,1) * t406 - t192) + t305;
t101 = -t164 * t253 + t306;
t66 = qJD(4) * t321 + qJD(3);
t49 = -t185 * t378 + (-t141 - t94 + t163) * t254 + t305;
t48 = t254 * t92 + (-t185 * t243 - t221 * t406) * qJD(4) + t294;
t42 = qJD(3) + (t397 * t243 + t396 * t244) * qJD(4);
t25 = t278 * qJD(4);
t24 = t312 * t244 + (t346 * t379 - t141 + t227 - t436) * t254 + t305;
t23 = t312 * t243 + (-t314 + t226 + t437) * t254 + t294;
t5 = ((t437 + t489) * t244 + (-t396 * t254 + t436) * t243) * qJD(4);
t1 = [m(3) * (t121 * (-t173 - t442) + t120 * (t174 + t252) + (-t160 - t370 + t140) * t139) + t265 + (t24 * (t290 - t442) + t50 * (t313 - t370) + t23 * (t252 + t281) + t263 + (-t279 + t50 + t374 + t335 + t449) * t51) * m(6) + (t49 * (t280 - t442) + t60 * (t364 - t370) + t48 * (t252 + t462) + t264 + (-t297 + t60 + t374 - t369 + t458) * t61) * m(5) + (t102 * (t296 - t442) + t116 * (t192 - t370) + t101 * (t252 + t461) + t266 + (t116 + t158 + t374) * t117) * m(4); t265 + (t23 * t281 + t24 * t290 + t263 + (t220 * t378 + t254 * t441 + t226 - t317 + t449) * t51 + (t313 + t450) * t50) * m(6) + (t375 * t254 + t49 * t280 + t48 * t462 + t264 + (t362 + t458) * t61 + (t364 + t453) * t60) * m(5) + (t117 * t158 - (-t116 * t461 - t117 * t441) * t254 + t101 * t461 + t102 * t296 + t116 * t192 + t266) * m(4) + (-(-t139 * t174 - t140 * t173) * t254 + t120 * t174 - t121 * t173 - t139 * t160 - t140 * t417) * m(3); m(5) * t25 + m(6) * t5; -(((t380 + t382) * t259 + (t381 + t383) * t257) * t254 + (((-t389 - t391) * t244 + (t388 + t390) * t243) * t259 + ((t393 + t395) * t244 + (t392 + t394) * t243) * t257) * qJD(4)) * t254 / 0.2e1 + ((t475 * t254 + t483) * t244 + (t476 * t254 + t482) * t243) * t443 + ((-t379 * t412 + t411) * t243 + (t288 + (-t448 * t244 + (t413 + t292) * t243) * qJD(4)) * t244 + (-t379 * t415 + t414) * t243 + (t289 + (-t447 * t244 + (t416 + t293) * t243) * qJD(4)) * t244) * t354 + ((-t378 * t416 - t414) * t244 + (t289 + (t293 * t243 + (t415 - t447) * t244) * qJD(4)) * t243 + (-t378 * t413 - t411) * t244 + (t288 + (t292 * t243 + (t412 - t448) * t244) * qJD(4)) * t243) * t351 + (t25 * t321 + t66 * t278 + t330 * t185 + ((-t254 * t61 - t49) * t244 + (t254 * t60 - t48) * t243) * t221 - (t155 * t60 - t157 * t61) * t254 - (t66 * (-t155 * t243 - t157 * t244) + t330 * t331) * qJD(4)) * m(5) + (t481 * t254 + ((t471 * t244 + t477 * t254) * t244 + (t465 * t243 + t478 * t254 + (-t466 + t472) * t244) * t243) * t464) * t446 + (t480 * t254 + ((t466 * t244 + t479 * t254) * t244 + (t472 * t243 + t506 * t254 + (-t465 + t471) * t244) * t243) * t464) * t445 + ((-t23 * t346 + t51 * t334 + t5 * t397 + t42 * t436 + (t50 * t220 - t396 * t42) * t254) * t243 + (-t24 * t346 + t50 * t334 + t5 * t396 + t42 * t437 + (-t346 * t51 + t397 * t42) * t254) * t244 - (t50 * t154 + t337 * t51) * t254 - ((t337 * t42 + t345 * t50) * t244 + (t51 * t345 + (-pkin(4) * t408 - t154) * t42) * t243) * qJD(4)) * m(6) + (t485 + t487) * t410 / 0.2e1 + (t484 + t488) * t406 / 0.2e1; 0.2e1 * (t23 * t445 + t24 * t446) * m(6);];
tauc = t1(:);

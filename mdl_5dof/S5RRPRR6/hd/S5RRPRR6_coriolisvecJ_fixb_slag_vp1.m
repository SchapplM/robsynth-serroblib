% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:50
% EndTime: 2022-01-20 11:17:11
% DurationCPUTime: 12.67s
% Computational Cost: add. (24245->670), mult. (23385->928), div. (0->0), fcn. (22651->10), ass. (0->351)
t332 = qJ(1) + qJ(2);
t324 = sin(t332);
t326 = cos(t332);
t337 = cos(qJ(4));
t334 = cos(pkin(9));
t335 = sin(qJ(4));
t463 = t334 * t335;
t264 = t324 * t463 + t326 * t337;
t462 = t334 * t337;
t467 = t326 * t335;
t265 = t324 * t462 - t467;
t333 = sin(pkin(9));
t475 = t324 * t333;
t168 = Icges(5,5) * t265 - Icges(5,6) * t264 + Icges(5,3) * t475;
t247 = Icges(5,4) * t265;
t171 = -Icges(5,2) * t264 + Icges(5,6) * t475 + t247;
t246 = Icges(5,4) * t264;
t175 = -Icges(5,1) * t265 - Icges(5,5) * t475 + t246;
t79 = -t168 * t334 - (t171 * t335 + t175 * t337) * t333;
t331 = qJ(4) + qJ(5);
t323 = sin(t331);
t325 = cos(t331);
t472 = t325 * t326;
t474 = t324 * t334;
t241 = t323 * t474 + t472;
t471 = t326 * t323;
t242 = t325 * t474 - t471;
t150 = Icges(6,5) * t242 - Icges(6,6) * t241 + Icges(6,3) * t475;
t231 = Icges(6,4) * t242;
t153 = -Icges(6,2) * t241 + Icges(6,6) * t475 + t231;
t230 = Icges(6,4) * t241;
t157 = -Icges(6,1) * t242 - Icges(6,5) * t475 + t230;
t74 = -t150 * t334 - (t153 * t323 + t157 * t325) * t333;
t317 = t326 * qJ(3);
t285 = pkin(2) * t324 - t317;
t330 = qJD(1) + qJD(2);
t275 = t330 * t285;
t314 = qJD(3) * t324;
t470 = t326 * t330;
t299 = qJ(3) * t470;
t448 = t299 + t314;
t551 = t448 - t314 + t275;
t235 = -Icges(6,3) * t334 + (Icges(6,5) * t325 - Icges(6,6) * t323) * t333;
t486 = Icges(6,4) * t325;
t236 = -Icges(6,6) * t334 + (-Icges(6,2) * t323 + t486) * t333;
t487 = Icges(6,4) * t323;
t237 = -Icges(6,5) * t334 + (Icges(6,1) * t325 - t487) * t333;
t468 = t326 * t334;
t477 = t324 * t325;
t243 = -t323 * t468 + t477;
t478 = t323 * t324;
t244 = t325 * t468 + t478;
t469 = t326 * t333;
t101 = t235 * t469 + t236 * t243 + t237 * t244;
t329 = qJD(4) + qJD(5);
t413 = t329 * t333;
t259 = t324 * t413;
t260 = t326 * t413;
t529 = -t329 * t334 + t330;
t66 = t150 * t469 + t243 * t153 - t157 * t244;
t152 = Icges(6,5) * t244 + Icges(6,6) * t243 + Icges(6,3) * t469;
t488 = Icges(6,4) * t244;
t155 = Icges(6,2) * t243 + Icges(6,6) * t469 + t488;
t232 = Icges(6,4) * t243;
t158 = Icges(6,1) * t244 + Icges(6,5) * t469 + t232;
t67 = t152 * t469 + t243 * t155 + t244 * t158;
t25 = t101 * t529 + t259 * t66 + t260 * t67;
t549 = -t171 * t264 - t175 * t265;
t266 = t324 * t337 - t326 * t463;
t473 = t324 * t335;
t267 = t326 * t462 + t473;
t548 = t266 * t171 - t175 * t267;
t64 = t150 * t475 - t153 * t241 - t157 * t242;
t542 = t168 * t469;
t339 = -pkin(8) - pkin(7);
t464 = t333 * t339;
t502 = pkin(4) * qJD(4);
t541 = t330 * t464 + t337 * t502;
t336 = sin(qJ(1));
t503 = pkin(1) * qJD(1);
t439 = t336 * t503;
t286 = rSges(3,1) * t324 + rSges(3,2) * t326;
t479 = t286 * t330;
t239 = -t439 - t479;
t506 = pkin(7) * t333;
t508 = pkin(3) * t334;
t400 = t506 + t508;
t269 = t400 * t324;
t245 = t330 * t269;
t540 = t245 + t551;
t394 = rSges(6,1) * t242 - rSges(6,2) * t241;
t161 = rSges(6,3) * t475 + t394;
t322 = pkin(4) * t337 + pkin(3);
t443 = pkin(4) * t467;
t403 = -t322 * t474 + t443;
t504 = pkin(7) + t339;
t177 = (t333 * t504 + t508) * t324 + t403;
t505 = pkin(3) - t322;
t229 = -t333 * t505 + t334 * t504;
t444 = qJD(4) * t333;
t428 = t324 * t444;
t208 = t229 * t428;
t238 = -rSges(6,3) * t334 + (rSges(6,1) * t325 - rSges(6,2) * t323) * t333;
t312 = -qJD(4) * t334 + t330;
t539 = -t161 * t529 + t312 * t177 + t238 * t259 + t208;
t396 = rSges(5,1) * t265 - rSges(5,2) * t264;
t179 = rSges(5,3) * t475 + t396;
t258 = -rSges(5,3) * t334 + (rSges(5,1) * t337 - rSges(5,2) * t335) * t333;
t538 = -t179 * t312 + t258 * t428;
t369 = -t508 - pkin(2) + (-rSges(5,3) - pkin(7)) * t333;
t194 = -qJD(4) * t265 + t266 * t330;
t195 = -qJD(4) * t264 + t267 * t330;
t397 = rSges(5,1) * t195 + rSges(5,2) * t194;
t465 = t330 * t333;
t432 = t326 * t465;
t114 = rSges(5,3) * t432 + t397;
t427 = t326 * t444;
t217 = t258 * t427;
t374 = (-rSges(5,1) * t335 - rSges(5,2) * t337) * t333;
t274 = qJD(4) * t374;
t338 = cos(qJ(1));
t327 = t338 * pkin(1);
t341 = qJD(1) ^ 2;
t441 = t341 * t327;
t445 = qJD(3) * t330;
t404 = t326 * t445 - t441;
t287 = t326 * pkin(2) + t324 * qJ(3);
t315 = qJD(3) * t326;
t219 = t287 * t330 - t315;
t454 = -t400 * t470 - t219;
t62 = t274 * t428 - t114 * t312 + (t217 + t454) * t330 + t404;
t494 = t62 * t324;
t402 = t314 - t439;
t361 = (-t269 - t285) * t330 + t402;
t86 = t361 + t538;
t438 = t338 * t503;
t401 = -t315 + t438;
t270 = pkin(3) * t468 + pkin(7) * t469;
t450 = t270 + t287;
t532 = t330 * t450;
t360 = t401 + t532;
t181 = t267 * rSges(5,1) + t266 * rSges(5,2) + rSges(5,3) * t469;
t526 = t181 * t312 - t217;
t87 = t360 + t526;
t537 = (t86 * t369 * t326 + (-t86 * qJ(3) + t87 * (-rSges(5,3) * t333 - pkin(2) - t400)) * t324) * t330 + t369 * t494;
t382 = -rSges(6,3) * t333 - t322 * t334 - pkin(2);
t311 = pkin(4) * t473;
t411 = t502 * t463;
t430 = -t324 * t411 - t541 * t326;
t522 = t334 * t505 + t506;
t118 = (-t522 * t326 + t311) * t330 + t430;
t209 = t229 * t427;
t221 = t330 * t260;
t268 = (-rSges(6,1) * t323 - rSges(6,2) * t325) * t333;
t225 = t329 * t268;
t507 = pkin(4) * t335;
t414 = t333 ^ 2 * qJD(4) ^ 2 * t507;
t391 = t529 * t325;
t416 = t330 * t334 - t329;
t166 = t324 * t391 - t416 * t471;
t392 = t323 * t529;
t167 = t324 * t392 + t416 * t472;
t395 = rSges(6,1) * t167 + rSges(6,2) * t166;
t99 = rSges(6,3) * t432 + t395;
t44 = -t324 * t414 - t118 * t312 + t221 * t238 + t225 * t259 - t529 * t99 + (t209 + t454) * t330 + t404;
t496 = t44 * t324;
t59 = t361 + t539;
t163 = t244 * rSges(6,1) + t243 * rSges(6,2) + rSges(6,3) * t469;
t383 = t322 * t468 - t326 * t464 + t311;
t178 = t383 - t270;
t521 = t163 * t529 + t178 * t312 - t260 * t238 - t209;
t60 = t360 + t521;
t536 = (t59 * t382 * t326 + (t59 * (-qJ(3) - t507) + t60 * t382) * t324) * t330 + (-pkin(2) + (-rSges(6,3) + t339) * t333) * t496;
t72 = t542 + t548;
t534 = t72 - t542;
t530 = -rSges(4,1) * t468 - t324 * rSges(4,3);
t367 = -rSges(4,2) * t469 + t287 - t530;
t533 = t330 * t367;
t440 = rSges(4,1) * t474;
t447 = rSges(4,2) * t475 + t326 * rSges(4,3);
t227 = t440 - t447;
t433 = t324 * t465;
t449 = rSges(4,2) * t433 + rSges(4,3) * t470;
t531 = t330 * t227 + t449;
t362 = t541 * t324 - t326 * t411 + t330 * t443;
t164 = t326 * t391 + t416 * t478;
t165 = t326 * t392 - t416 * t477;
t461 = t165 * rSges(6,1) + t164 * rSges(6,2);
t528 = t362 + t461 - t539 + t540;
t527 = t299 + t275;
t65 = t152 * t475 - t241 * t155 + t242 * t158;
t170 = Icges(5,5) * t267 + Icges(5,6) * t266 + Icges(5,3) * t469;
t491 = Icges(5,4) * t267;
t173 = Icges(5,2) * t266 + Icges(5,6) * t469 + t491;
t248 = Icges(5,4) * t266;
t176 = Icges(5,1) * t267 + Icges(5,5) * t469 + t248;
t71 = t170 * t475 - t264 * t173 + t265 * t176;
t192 = -qJD(4) * t267 + t264 * t330;
t193 = qJD(4) * t266 - t265 * t330;
t455 = t193 * rSges(5,1) + t192 * rSges(5,2);
t523 = t455 - t538;
t520 = t324 * (-Icges(5,2) * t265 - t175 - t246) + t326 * (-Icges(5,2) * t267 + t176 + t248);
t262 = (-Icges(6,2) * t325 - t487) * t333;
t519 = t259 * (-Icges(6,2) * t242 - t157 - t230) + t260 * (-Icges(6,2) * t244 + t158 + t232) + t529 * (t237 + t262);
t220 = t329 * t433;
t518 = -t220 / 0.2e1;
t517 = t221 / 0.2e1;
t516 = -t259 / 0.2e1;
t515 = t259 / 0.2e1;
t513 = t260 / 0.2e1;
t511 = -t326 / 0.2e1;
t510 = -t334 / 0.2e1;
t509 = pkin(1) * t336;
t70 = t168 * t475 + t549;
t499 = t324 * t70;
t93 = Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t432;
t95 = Icges(6,4) * t167 + Icges(6,2) * t166 + Icges(6,6) * t432;
t97 = Icges(6,1) * t167 + Icges(6,4) * t166 + Icges(6,5) * t432;
t34 = -t334 * t93 + ((-t153 * t329 + t97) * t325 + (t157 * t329 - t95) * t323) * t333;
t498 = t34 * t259;
t92 = Icges(6,5) * t165 + Icges(6,6) * t164 - Icges(6,3) * t433;
t94 = Icges(6,4) * t165 + Icges(6,2) * t164 - Icges(6,6) * t433;
t96 = Icges(6,1) * t165 + Icges(6,4) * t164 - Icges(6,5) * t433;
t35 = -t334 * t92 + ((-t155 * t329 + t96) * t325 + (-t158 * t329 - t94) * t323) * t333;
t497 = t35 * t260;
t495 = t60 * t225;
t493 = t74 * t221;
t75 = -t152 * t334 + (-t155 * t323 + t158 * t325) * t333;
t492 = t75 * t220;
t490 = Icges(5,4) * t335;
t489 = Icges(5,4) * t337;
t253 = -Icges(5,3) * t334 + (Icges(5,5) * t337 - Icges(5,6) * t335) * t333;
t254 = -Icges(5,6) * t334 + (-Icges(5,2) * t335 + t489) * t333;
t255 = -Icges(5,5) * t334 + (Icges(5,1) * t337 - t490) * t333;
t115 = t253 * t475 - t254 * t264 + t255 * t265;
t484 = t115 * t312;
t476 = t324 * t330;
t458 = -t163 - t178;
t278 = (-Icges(5,1) * t335 - t489) * t333;
t452 = -t254 + t278;
t277 = (-Icges(5,2) * t337 - t490) * t333;
t451 = t255 + t277;
t442 = t341 * t509;
t435 = t225 * t475 + t238 * t432 + t334 * t99;
t211 = t238 * t475;
t73 = t170 * t469 + t266 * t173 + t267 * t176;
t426 = t475 / 0.2e1;
t425 = t469 / 0.2e1;
t424 = -rSges(4,1) * t334 - pkin(2);
t423 = -t444 / 0.2e1;
t422 = t444 / 0.2e1;
t288 = t326 * rSges(3,1) - rSges(3,2) * t324;
t98 = -rSges(6,3) * t433 + t461;
t421 = -t161 * t330 - t98;
t117 = t522 * t476 + t362;
t420 = t177 * t330 - t117;
t188 = -rSges(6,1) * t241 - rSges(6,2) * t242;
t189 = rSges(6,1) * t243 - rSges(6,2) * t244;
t419 = t260 * t188 - t189 * t259;
t418 = t189 * t529 - t260 * t268;
t417 = -t188 * t529 + t259 * t268;
t410 = -t433 / 0.2e1;
t409 = t330 * t425;
t408 = t324 * t423;
t407 = t324 * t422;
t406 = t326 * t423;
t405 = t326 * t422;
t257 = rSges(3,1) * t470 - rSges(3,2) * t476;
t393 = t324 * t86 - t326 * t87;
t390 = t315 - t532;
t389 = (-Icges(5,5) * t264 - Icges(5,6) * t265) * t324 + (Icges(5,5) * t266 - Icges(5,6) * t267) * t326;
t388 = t330 * t408;
t387 = t330 * t405;
t386 = t330 * (-pkin(2) * t476 + t448) + t324 * t445 - t442;
t381 = t315 - t397;
t380 = t317 - t396;
t375 = -t330 ^ 2 * t269 + t386;
t373 = (t326 * t71 + t499) * t333;
t372 = (t324 * t72 + t326 * t73) * t333;
t370 = t333 * t389;
t263 = (-Icges(6,1) * t323 - t486) * t333;
t276 = (-Icges(5,5) * t335 - Icges(5,6) * t337) * t333;
t261 = (-Icges(6,5) * t323 - Icges(6,6) * t325) * t333;
t368 = t181 + t450;
t366 = (-Icges(6,5) * t241 - Icges(6,6) * t242) * t259 + (Icges(6,5) * t243 - Icges(6,6) * t244) * t260 + t261 * t529;
t105 = (t179 * t326 - t181 * t324) * t444;
t365 = (Icges(5,1) * t266 - t173 - t491) * t326 + (-Icges(5,1) * t264 - t171 - t247) * t324;
t363 = t324 * t424 + t317 + t447;
t15 = -t161 * t220 - t163 * t221 - t259 * t98 + t260 * t99 + ((-t178 * t330 + t118) * t326 + t420 * t324) * t444;
t54 = t161 * t260 - t163 * t259 + (-t177 * t326 - t178 * t324) * t444;
t358 = t44 * (t334 * t161 + t211) + (t15 * t161 + t54 * t99) * t469;
t100 = t235 * t475 - t236 * t241 + t237 * t242;
t357 = t315 - t395 - t430;
t356 = t333 * t366;
t355 = t317 - t394 + t403;
t16 = t153 * t164 - t157 * t165 + t243 * t95 + t244 * t97 + (-t150 * t476 + t326 * t93) * t333;
t17 = t155 * t164 + t158 * t165 + t243 * t94 + t244 * t96 + (-t152 * t476 + t326 * t92) * t333;
t18 = t153 * t166 - t157 * t167 - t241 * t95 + t242 * t97 + (t150 * t470 + t324 * t93) * t333;
t19 = t155 * t166 + t158 * t167 - t241 * t94 + t242 * t96 + (t152 * t470 + t324 * t92) * t333;
t348 = (Icges(6,1) * t243 - t155 - t488) * t260 + (-Icges(6,1) * t241 - t153 - t231) * t259 + (-t236 + t263) * t529;
t222 = t329 * t261;
t223 = t329 * t262;
t224 = t329 * t263;
t50 = t164 * t236 + t165 * t237 + t223 * t243 + t224 * t244 + (t222 * t326 - t235 * t476) * t333;
t51 = t166 * t236 + t167 * t237 - t223 * t241 + t224 * t242 + (t222 * t324 + t235 * t470) * t333;
t81 = -t222 * t334 + ((-t236 * t329 + t224) * t325 + (-t237 * t329 - t223) * t323) * t333;
t78 = t81 * t529;
t354 = (t18 * t259 + t19 * t260 - t220 * t65 + t221 * t64 + t51 * t529) * t426 + (-t241 * t519 + t242 * t348 + t324 * t356) * t516 - (t519 * t243 + t348 * t244 + t326 * t356) * t260 / 0.2e1 - (-t366 * t334 + (-t323 * t519 + t325 * t348) * t333) * t529 / 0.2e1 + (t16 * t259 + t17 * t260 - t220 * t67 + t221 * t66 + t50 * t529) * t425 + t25 * t410 + (t100 * t529 + t259 * t64 + t260 * t65) * t409 + (-t101 * t334 + (t324 * t66 + t326 * t67) * t333) * t518 + (-t100 * t334 + (t324 * t64 + t326 * t65) * t333) * t517 + (-t334 * t51 + ((t330 * t64 + t19) * t326 + (-t330 * t65 + t18) * t324) * t333) * t515 + (-t334 * t50 + ((t330 * t66 + t17) * t326 + (-t330 * t67 + t16) * t324) * t333) * t513 + (-t492 + t493 + t497 + t78 + t498) * t510 + t529 * (-t334 * t81 + ((t330 * t74 + t35) * t326 + (-t330 * t75 + t34) * t324) * t333) / 0.2e1;
t353 = t383 + t163 + t287;
t107 = Icges(5,5) * t193 + Icges(5,6) * t192 - Icges(5,3) * t433;
t108 = Icges(5,5) * t195 + Icges(5,6) * t194 + Icges(5,3) * t432;
t109 = Icges(5,4) * t193 + Icges(5,2) * t192 - Icges(5,6) * t433;
t110 = Icges(5,4) * t195 + Icges(5,2) * t194 + Icges(5,6) * t432;
t111 = Icges(5,1) * t193 + Icges(5,4) * t192 - Icges(5,5) * t433;
t112 = Icges(5,1) * t195 + Icges(5,4) * t194 + Icges(5,5) * t432;
t351 = ((t330 * t72 + t109 * t266 + t111 * t267 + t173 * t192 + t176 * t193 + (t107 * t326 - t170 * t476) * t333) * t326 + (-t330 * t73 + t110 * t266 + t112 * t267 + t171 * t192 - t175 * t193 + (t108 * t326 - t168 * t476) * t333) * t324) * t333;
t350 = ((t330 * t70 - t109 * t264 + t111 * t265 + t173 * t194 + t176 * t195 + (t107 * t324 + t170 * t470) * t333) * t326 + (-t330 * t71 - t110 * t264 + t112 * t265 + t171 * t194 - t175 * t195 + (t108 * t324 + t168 * t470) * t333) * t324) * t333;
t45 = -t108 * t334 + (-t110 * t335 + t112 * t337 + (-t171 * t337 + t175 * t335) * qJD(4)) * t333;
t46 = -t107 * t334 + (-t109 * t335 + t111 * t337 + (-t173 * t337 - t176 * t335) * qJD(4)) * t333;
t80 = -t170 * t334 + (-t173 * t335 + t176 * t337) * t333;
t349 = ((t330 * t79 + t46) * t326 + (-t330 * t80 + t45) * t324) * t333;
t159 = (-t227 - t285) * t330 + t402;
t160 = t401 + t533;
t345 = (t159 * t424 * t326 + (t159 * (-rSges(4,3) - qJ(3)) + t160 * t424) * t324) * t330;
t116 = t253 * t469 + t254 * t266 + t255 * t267;
t106 = t116 * t312;
t39 = qJD(4) * t373 + t484;
t40 = qJD(4) * t372 + t106;
t271 = qJD(4) * t276;
t272 = qJD(4) * t277;
t273 = qJD(4) * t278;
t57 = t192 * t254 + t193 * t255 + t266 * t272 + t267 * t273 + (-t253 * t476 + t271 * t326) * t333;
t58 = t194 * t254 + t195 * t255 - t264 * t272 + t265 * t273 + (t253 * t470 + t271 * t324) * t333;
t102 = -t271 * t334 + (-t272 * t335 + t273 * t337 + (-t254 * t337 - t255 * t335) * qJD(4)) * t333;
t91 = t102 * t312;
t344 = (t106 + ((-t549 + t70 + t73) * t326 + t534 * t324) * t444) * t408 + t25 * t516 + t91 + t78 + t498 / 0.2e1 + t497 / 0.2e1 + t50 * t513 + t100 * t517 + t101 * t518 - t492 / 0.2e1 + t493 / 0.2e1 + (t51 + t25) * t515 + (t46 + t57) * t405 + (t116 + t80) * t388 + (t115 + t79) * t387 + (t58 + t40 + t45) * t407 + (-t484 + ((t534 - t548 - t71) * t326 - t499) * t444 + t39) * t406;
t292 = rSges(4,2) * t432;
t250 = t266 * pkin(4);
t249 = t264 * pkin(4);
t240 = t288 * t330 + t438;
t214 = -t257 * t330 - t441;
t213 = -t330 * t479 - t442;
t206 = t330 * t211;
t204 = rSges(5,1) * t266 - rSges(5,2) * t267;
t203 = -rSges(5,1) * t264 - rSges(5,2) * t265;
t120 = (t530 * t330 - t219 + t292) * t330 + t404;
t119 = t330 * (-t330 * t440 + t449) + t386;
t113 = -rSges(5,3) * t433 + t455;
t61 = t113 * t312 + (t258 * t476 - t274 * t326) * t444 + t375;
t43 = t117 * t312 + t208 * t330 + t220 * t238 - t225 * t260 + t326 * t414 + t529 * t98 + t375;
t1 = [m(3) * (t214 * (-t286 - t509) + t213 * (t288 + t327) + (-t257 - t438 + t240) * t239) + t344 + (t44 * (t355 - t509) + t59 * (t357 - t438) + t43 * (t327 + t353) + (t59 + t528) * t60 + t536) * m(6) + (t62 * (t380 - t509) + t86 * (t381 - t438) + t61 * (t327 + t368) + (t245 + t86 + t523 + t527) * t87 + t537) * m(5) + (t120 * (t363 - t509) + t159 * (t292 - t401) + t119 * (t327 + t367) + t345 + (t159 + t527 + t531) * t160) * m(4); t344 + (t43 * t353 + t44 * t355 + t528 * t60 + (-t390 + t357 + t521) * t59 + t536) * m(6) + (t61 * t368 + t62 * t380 + (t523 + t540) * t87 + (-t390 + t381 + t526) * t86 + t537) * m(5) + (t119 * t367 + t120 * t363 + t345 + (t531 + t551) * t160 + (t292 + t533) * t159) * m(4) + (-(-t239 * t288 - t240 * t286) * t330 + t213 * t288 - t214 * t286 - t239 * t257 - t240 * t479) * m(3); 0.2e1 * (t496 / 0.2e1 + t43 * t511) * m(6) + 0.2e1 * (t494 / 0.2e1 + t61 * t511) * m(5) + 0.2e1 * (t119 * t511 + t120 * t324 / 0.2e1) * m(4); (qJD(4) * t351 + t312 * t57) * t425 + (qJD(4) * t350 + t312 * t58) * t426 + t40 * t410 + t39 * t409 + t354 + ((t266 * t451 + t267 * t452 + t276 * t469) * t312 + (t520 * t266 + t365 * t267 + t326 * t370) * t444) * t406 + ((-t264 * t451 + t265 * t452 + t276 * t475) * t312 + (-t264 * t520 + t265 * t365 + t324 * t370) * t444) * t408 + (-t334 * t57 + t351) * t405 + (-t334 * t58 + t350) * t407 + (-t116 * t334 + t372) * t388 + (-t115 * t334 + t373) * t387 + (qJD(4) * t349 + t91) * t510 + t312 * (-t102 * t334 + t349) / 0.2e1 - t312 * (-t276 * t334 * t312 + ((-t335 * t451 + t337 * t452) * t312 + ((-t335 * t520 + t337 * t365) * t333 - t389 * t334) * qJD(4)) * t333) / 0.2e1 + (t59 * t435 + t60 * t206 + (-t44 * t177 + t59 * t118 + t43 * t458 + t60 * (-t117 - t98)) * t334 + ((t43 * (-t229 - t238) - t495 - t15 * t177 + t54 * t118 + (t59 * t229 + t458 * t54) * t330) * t326 + (t15 * t458 + t54 * (t420 + t421) + (t330 * t60 + t44) * t229) * t324) * t333 + t358 - t59 * (t249 * t312 + t417) - t60 * (t250 * t312 + t418) - t54 * ((-t249 * t326 - t250 * t324) * t444 + t419)) * m(6) + (-(-t203 * t86 + t204 * t87) * t312 - (t105 * (t203 * t326 - t204 * t324) + t393 * t374) * t444 + (-t113 * t87 + t114 * t86 + t179 * t62 - t181 * t61) * t334 + (0.2e1 * t105 * ((-t181 * t330 + t114) * t326 + (-t179 * t330 - t113) * t324) + t393 * t274 + ((t330 * t86 - t61) * t326 + (t330 * t87 + t62) * t324) * t258) * t333) * m(5); t354 + (-t43 * t163 * t334 + ((-t163 * t330 * t54 - t238 * t43 - t495) * t326 + (-t15 * t163 + t421 * t54) * t324) * t333 + t358 - t54 * t419 + (-t334 * t98 + t206 - t418) * t60 + (t435 - t417) * t59) * m(6);];
tauc = t1(:);

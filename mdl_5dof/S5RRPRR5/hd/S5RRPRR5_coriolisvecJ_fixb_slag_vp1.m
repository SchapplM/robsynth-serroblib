% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR5
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:31
% EndTime: 2020-01-03 12:03:45
% DurationCPUTime: 8.95s
% Computational Cost: add. (18232->593), mult. (12208->752), div. (0->0), fcn. (9444->10), ass. (0->350)
t487 = rSges(4,2) * sin(pkin(9));
t539 = pkin(2) - t487;
t287 = qJ(1) + qJ(2);
t277 = sin(t287);
t284 = pkin(9) + qJ(4);
t275 = cos(t284);
t263 = Icges(5,4) * t275;
t274 = sin(t284);
t344 = -Icges(5,2) * t274 + t263;
t278 = cos(t287);
t286 = qJD(1) + qJD(2);
t441 = t278 * t286;
t476 = Icges(5,4) * t274;
t215 = Icges(5,2) * t275 + t476;
t512 = -Icges(5,6) * t286 + qJD(4) * t215;
t111 = -t277 * t512 + t344 * t441;
t218 = Icges(5,1) * t275 - t476;
t329 = t218 * t286;
t523 = Icges(5,1) * t274 + t263;
t509 = -Icges(5,5) * t286 + qJD(4) * t523;
t113 = -t277 * t509 + t278 * t329;
t444 = t275 * t278;
t446 = t274 * t278;
t154 = Icges(5,4) * t444 - Icges(5,2) * t446 + Icges(5,6) * t277;
t239 = Icges(5,4) * t446;
t156 = Icges(5,1) * t444 + Icges(5,5) * t277 - t239;
t340 = t154 * t275 + t156 * t274;
t155 = -Icges(5,5) * t278 + t218 * t277;
t457 = t155 * t275;
t327 = t344 * t277;
t153 = -Icges(5,6) * t278 + t327;
t459 = t153 * t274;
t341 = t457 - t459;
t538 = -qJD(4) * t341 - t111 * t275 - t113 * t274 + t286 * t340;
t522 = t277 * rSges(3,1) + t278 * rSges(3,2);
t191 = t522 * t286;
t291 = sin(qJ(1));
t483 = pkin(1) * qJD(1);
t393 = t291 * t483;
t174 = t393 + t191;
t276 = qJ(5) + t284;
t268 = sin(t276);
t269 = cos(t276);
t475 = Icges(6,4) * t268;
t205 = Icges(6,1) * t269 - t475;
t328 = t205 * t277;
t142 = -Icges(6,5) * t278 + t328;
t450 = t268 * t278;
t232 = Icges(6,4) * t450;
t448 = t269 * t278;
t143 = Icges(6,1) * t448 + Icges(6,5) * t277 - t232;
t202 = Icges(6,2) * t269 + t475;
t285 = qJD(4) + qJD(5);
t221 = t277 * t285;
t442 = t278 * t285;
t262 = Icges(6,4) * t269;
t343 = -Icges(6,2) * t268 + t262;
t524 = Icges(6,1) * t268 + t262;
t532 = t524 + t343;
t302 = t221 * (-Icges(6,2) * t448 + t143 - t232) - t442 * (-t202 * t277 + t142) + t286 * t532;
t326 = t343 * t277;
t140 = -Icges(6,6) * t278 + t326;
t141 = Icges(6,4) * t448 - Icges(6,2) * t450 + Icges(6,6) * t277;
t503 = t221 * (t278 * t524 + t141) - t442 * (t277 * t524 + t140) + t286 * (t202 - t205);
t537 = t302 * t268 + t269 * t503;
t311 = t277 * (-Icges(5,2) * t444 + t156 - t239) - t278 * (-t215 * t277 + t155);
t507 = -t278 * t153 + t154 * t277;
t536 = -t311 * t274 - t275 * t507;
t535 = -rSges(4,3) - qJ(3);
t417 = t523 + t344;
t418 = t215 - t218;
t533 = (t274 * t417 + t275 * t418) * t286;
t531 = 0.2e1 * qJD(4);
t530 = rSges(5,2) * t275;
t201 = Icges(6,5) * t269 - Icges(6,6) * t268;
t138 = -Icges(6,3) * t278 + t201 * t277;
t461 = t141 * t268;
t342 = -t143 * t269 + t461;
t334 = -t138 + t342;
t529 = t442 * t334;
t289 = cos(pkin(9));
t491 = rSges(4,1) * t289;
t257 = t277 * t491;
t335 = t539 * t277 + t257;
t528 = t286 * (t278 * t535 + t335);
t466 = qJ(3) * t277;
t225 = pkin(2) * t278 + t466;
t270 = t289 * pkin(3) + pkin(2);
t290 = -pkin(7) - qJ(3);
t412 = t278 * t270 - t277 * t290;
t137 = t225 - t412;
t195 = t286 * t225;
t525 = -t286 * t137 + t195;
t258 = t278 * t491;
t409 = t277 * rSges(4,3) + t258;
t170 = -t278 * t487 + t409;
t407 = qJD(3) * t278;
t443 = t277 * t286;
t411 = pkin(2) * t441 + qJ(3) * t443;
t414 = rSges(4,3) * t443 + t286 * t258;
t521 = -t286 * t170 - t195 + t407 + t411 + t414;
t243 = rSges(5,2) * t446;
t400 = rSges(5,1) * t444;
t160 = rSges(5,3) * t277 - t243 + t400;
t135 = t286 * t160;
t210 = t270 * t441;
t219 = rSges(5,1) * t274 + t530;
t404 = qJD(4) * t277;
t333 = -t219 * t404 - t407;
t419 = rSges(5,3) * t443 + t286 * t400;
t520 = t210 + t419 - t135 - t333 - t525;
t492 = pkin(4) * t275;
t229 = t270 + t492;
t194 = t278 * t229;
t283 = -pkin(8) + t290;
t121 = t277 * t283 - t194 + t412;
t118 = t286 * t121;
t235 = rSges(6,2) * t450;
t398 = rSges(6,1) * t448;
t145 = rSges(6,3) * t277 - t235 + t398;
t134 = t286 * t145;
t188 = t229 * t441;
t484 = rSges(6,2) * t269;
t489 = rSges(6,1) * t268;
t206 = t484 + t489;
t382 = t274 * t404;
t323 = -pkin(4) * t382 - t407;
t308 = -t206 * t221 + t323;
t421 = rSges(6,3) * t443 + t286 * t398;
t519 = t188 + t421 + t118 - t134 - t308 - t525;
t454 = t202 * t285;
t518 = -Icges(6,6) * t286 + t454;
t171 = t206 * t277;
t172 = t206 * t278;
t485 = rSges(6,2) * t268;
t488 = rSges(6,1) * t269;
t207 = -t485 + t488;
t261 = t278 * t290;
t413 = t277 * t270 + t261;
t426 = t277 * t229 + t278 * t283;
t120 = -t413 + t426;
t449 = t269 * t277;
t234 = rSges(6,1) * t449;
t451 = t268 * t277;
t144 = -rSges(6,2) * t451 - rSges(6,3) * t278 + t234;
t57 = t144 * t221 + t145 * t442 + (t120 * t277 - t121 * t278) * qJD(4);
t264 = qJD(3) * t277;
t403 = qJD(4) * t278;
t381 = t274 * t403;
t361 = pkin(4) * t381;
t332 = t206 * t442 - t264 + t361;
t360 = t120 + t144 + t413;
t58 = t286 * t360 + t332 + t393;
t292 = cos(qJ(1));
t273 = t292 * t483;
t416 = t225 - t137;
t439 = -t121 + t145;
t59 = t273 + (t416 + t439) * t286 + t308;
t517 = -(-t286 * t171 + t207 * t442) * t58 - t57 * (-t171 * t221 - t172 * t442) - t59 * (-t172 * t286 - t221 * t207);
t200 = Icges(6,5) * t268 + Icges(6,6) * t269;
t516 = -Icges(6,3) * t286 + t200 * t285;
t515 = -Icges(6,5) * t286 + t285 * t524;
t103 = t153 * t275 + t155 * t274;
t214 = Icges(5,5) * t275 - Icges(5,6) * t274;
t325 = t214 * t277;
t151 = -Icges(5,3) * t278 + t325;
t514 = qJD(4) * t103 + t111 * t274 - t113 * t275 - t151 * t286;
t213 = Icges(5,5) * t274 + Icges(5,6) * t275;
t513 = -Icges(5,3) * t286 + qJD(4) * t213;
t197 = t344 * qJD(4);
t198 = t218 * qJD(4);
t511 = qJD(4) * (t215 * t275 + t274 * t523) + t197 * t274 - t198 * t275 - t213 * t286;
t110 = t278 * t512 + t286 * t327;
t112 = t277 * t329 + t278 * t509;
t152 = Icges(5,5) * t444 - Icges(5,6) * t446 + Icges(5,3) * t277;
t510 = qJD(4) * t340 - t110 * t274 + t112 * t275 - t152 * t286;
t362 = -t205 * t285 + t454;
t363 = t532 * t285;
t506 = -t200 * t286 + t268 * t363 + t269 * t362;
t139 = Icges(6,5) * t448 - Icges(6,6) * t450 + Icges(6,3) * t277;
t367 = t143 * t285 - t278 * t518 - t286 * t326;
t369 = t141 * t285 + t278 * t515 + t286 * t328;
t505 = -t139 * t286 + t268 * t367 + t269 * t369;
t368 = t142 * t285 - t277 * t518 + t343 * t441;
t370 = t140 * t285 - t205 * t441 + t277 * t515;
t504 = -t138 * t286 + t268 * t368 + t269 * t370;
t189 = t286 * t221;
t502 = t189 / 0.2e1;
t190 = t286 * t442;
t501 = -t190 / 0.2e1;
t500 = -t221 / 0.2e1;
t499 = t221 / 0.2e1;
t498 = t442 / 0.2e1;
t497 = -t442 / 0.2e1;
t496 = -t277 / 0.2e1;
t495 = -t278 / 0.2e1;
t494 = -t286 / 0.2e1;
t493 = t286 / 0.2e1;
t281 = t291 * pkin(1);
t282 = t292 * pkin(1);
t490 = rSges(5,1) * t275;
t486 = rSges(5,2) * t274;
t482 = t286 * t59;
t68 = t273 + (t160 + t416) * t286 + t333;
t481 = t286 * t68;
t456 = t200 * t277;
t88 = -t202 * t450 + t448 * t524 + t456;
t480 = t88 * t286;
t408 = t283 - t290;
t100 = t361 + (t408 * t278 + (t229 - t270) * t277) * t286;
t394 = t285 * t484;
t395 = t286 * t485;
t422 = rSges(6,3) * t441 + t277 * t395;
t97 = t278 * t394 + (t268 * t442 + t269 * t443) * rSges(6,1) - t422;
t479 = -t100 - t97;
t453 = t213 * t277;
t106 = -t215 * t446 + t444 * t523 + t453;
t465 = t106 * t286;
t462 = t140 * t268;
t460 = t142 * t269;
t458 = t154 * t274;
t164 = t200 * t278;
t324 = t201 * t286;
t181 = t213 * t278;
t452 = t214 * t286;
t447 = t274 * t277;
t445 = t275 * t277;
t440 = t286 * t290;
t246 = qJ(3) * t441;
t410 = t246 + t264;
t162 = pkin(2) * t443 - t410;
t438 = -t246 - (t261 + (-pkin(2) + t270) * t277) * t286 - t162;
t294 = qJD(1) ^ 2;
t272 = t294 * t282;
t432 = t286 * (-t407 + t411) + t272;
t425 = t194 - t235;
t396 = t286 * t486;
t420 = rSges(5,3) * t441 + t277 * t396;
t397 = t286 * t487;
t415 = rSges(4,3) * t441 + t277 * t397;
t406 = qJD(4) * t274;
t405 = qJD(4) * t275;
t402 = qJD(4) ^ 2 * t492;
t401 = t294 * t281;
t399 = t285 * t489;
t392 = pkin(4) * t406;
t98 = -t277 * t399 + (-t221 * t269 - t268 * t441) * rSges(6,2) + t421;
t391 = t144 * t441 + (-t134 + t98) * t277;
t390 = t286 * (-t277 * t440 + t210 - t411) + t432;
t242 = rSges(5,1) * t445;
t159 = -rSges(5,2) * t447 - rSges(5,3) * t278 + t242;
t389 = t159 + t413;
t387 = t234 + t426;
t385 = t242 + t413;
t384 = -t243 + t412;
t383 = t219 * t403;
t380 = t275 * t404;
t378 = t443 / 0.2e1;
t377 = -t441 / 0.2e1;
t375 = t404 / 0.2e1;
t373 = t403 / 0.2e1;
t371 = pkin(4) * t274 + t206;
t366 = -t139 - t462;
t365 = -t139 + t460;
t226 = rSges(3,1) * t278 - rSges(3,2) * t277;
t175 = t226 * t286 + t273;
t356 = t286 * t264 - t401;
t192 = rSges(3,1) * t441 - rSges(3,2) * t443;
t354 = -t264 + t393;
t353 = -qJD(3) - t397;
t352 = t409 + t466;
t351 = t264 - t383;
t348 = -t486 + t490;
t67 = t286 * t389 - t351 + t393;
t347 = -t277 * t68 + t278 * t67;
t90 = -t141 * t269 - t143 * t268;
t339 = -t156 * t275 + t458;
t338 = -t202 * t268 + t269 * t524;
t336 = -t215 * t274 + t275 * t523;
t129 = t155 * t445;
t70 = -t151 * t278 - t153 * t447 + t129;
t130 = t156 * t445;
t71 = t152 * t278 + t154 * t447 - t130;
t331 = (-t277 * t71 - t278 * t70) * qJD(4);
t131 = t153 * t446;
t72 = -t151 * t277 - t155 * t444 + t131;
t73 = t277 * t152 - t278 * t339;
t330 = (-t277 * t73 - t278 * t72) * qJD(4);
t102 = (t159 * t277 + t160 * t278) * qJD(4);
t320 = t164 * t221 - t442 * t456 - t324;
t318 = -t277 * t324 - t278 * t516 + t286 * t342;
t317 = -t278 * t324 + t516 * t277 + (t460 - t462) * t286;
t316 = -t278 * t513 + (-t325 + t339) * t286;
t315 = -t214 * t441 + t277 * t513 + t286 * t341;
t314 = -t201 * t285 + t286 * t338;
t313 = -t214 * qJD(4) + t286 * t336;
t13 = t317 * t277 + t278 * t504;
t14 = t318 * t277 - t278 * t505;
t15 = -t277 * t504 + t317 * t278;
t16 = t277 * t505 + t318 * t278;
t123 = t142 * t449;
t63 = -t138 * t278 - t140 * t451 + t123;
t124 = t143 * t449;
t64 = t139 * t278 + t141 * t451 - t124;
t87 = t277 * t338 - t164;
t83 = t87 * t286;
t28 = -t221 * t64 - t442 * t63 + t83;
t125 = t140 * t450;
t65 = -t138 * t277 - t142 * t448 + t125;
t66 = t139 * t277 - t278 * t342;
t29 = -t221 * t66 - t442 * t65 - t480;
t42 = t314 * t277 + t278 * t506;
t43 = -t277 * t506 + t314 * t278;
t44 = -t268 * t370 + t269 * t368;
t45 = t268 * t369 - t269 * t367;
t89 = t140 * t269 + t142 * t268;
t309 = (-t13 * t442 - t14 * t221 + t189 * t65 - t190 * t66 + t286 * t42) * t496 + (t320 * t277 + t537 * t278) * t499 + (-t537 * t277 + t320 * t278) * t498 + (-t15 * t442 - t16 * t221 + t189 * t63 - t190 * t64 + t286 * t43) * t495 + (-t268 * t503 + t269 * t302) * t494 + t28 * t378 + t29 * t377 + ((-t286 * t66 - t13) * t278 + (t286 * t65 - t14) * t277) * t500 + (-t277 * t64 - t278 * t63) * t502 + (-t277 * t66 - t278 * t65) * t501 + ((-t286 * t64 - t15) * t278 + (t286 * t63 - t16) * t277) * t497 + ((-t286 * t90 - t44) * t278 + (t286 * t89 - t45) * t277) * t493;
t307 = (-pkin(2) - t491) * t443 + t410 + t415;
t306 = -t392 - t394 - t399;
t116 = t354 + t528;
t79 = (-t257 * t286 - t162 + t415) * t286 + t356;
t80 = (t278 * t353 + t414) * t286 + t432;
t298 = (t116 * t353 + t80 * t535 + t79 * t539) * t278;
t105 = t277 * t336 - t181;
t99 = t105 * t286;
t34 = t99 + t331;
t35 = t330 - t465;
t49 = qJD(4) * t339 + t110 * t275 + t112 * t274;
t54 = t313 * t277 + t278 * t511;
t55 = -t277 * t511 + t313 * t278;
t297 = (t83 - (t277 * t366 + t123 + t66) * t442 + (t124 - t125 + t65 + (t138 - t461) * t277) * t221 + (t221 * t365 - t529) * t278) * t499 + (t99 + ((t72 + t130 - t131 + (t151 - t458) * t277) * t277 + (-t129 - t73 + (t151 - t339) * t278 + (t457 + t459) * t277) * t278) * qJD(4)) * t375 + t90 * t501 + t190 * t88 / 0.2e1 + (t89 + t87) * t502 + (t480 - (-t125 + t64) * t442 + (-t123 + t63) * t221 + (-t221 * t334 - t365 * t442) * t278 + (-t221 * t366 + t529) * t277 + t29) * t498 + (t44 + t43) * t497 + (t35 + t465 + ((t131 - t71 + (t152 - t457) * t278) * t278 + (-t129 + t70 + (t152 + t459) * t277) * t277) * qJD(4)) * t373 + (qJD(4) * t336 + t197 * t275 + t198 * t274 - t268 * t362 + t269 * t363) * t286 + (t45 + t42 + t28) * t500 - (t49 + t54 + t34) * t404 / 0.2e1 - (t55 - t538) * t403 / 0.2e1 + (t278 * t106 + (t103 + t105) * t277) * qJD(4) * t493;
t114 = t403 * t530 + (t275 * t443 + t381) * rSges(5,1) - t420;
t199 = t348 * qJD(4);
t52 = -t199 * t404 + (-t114 - t383 + t438) * t286 + t356;
t115 = -rSges(5,1) * t382 + (-t274 * t441 - t380) * rSges(5,2) + t419;
t53 = t199 * t403 + (t115 + t333) * t286 + t390;
t296 = (t52 * t490 + t68 * (-rSges(5,1) * t406 - rSges(5,2) * t405 - t440) - t53 * rSges(5,3) + t67 * (-qJD(3) - t396)) * t278 + (-t53 * t486 + t52 * rSges(5,3) - t67 * t219 * qJD(4) + (t68 * (-t270 - t490) - t67 * t290) * t286) * t277;
t179 = t207 * t285;
t30 = -t277 * t402 - t179 * t221 - t190 * t206 + (-t361 + t438 + t479) * t286 + t356;
t101 = t188 - t210 + (-t286 * t408 - t392) * t277;
t31 = t278 * t402 + t179 * t442 - t189 * t206 + (t101 + t323 + t98) * t286 + t390;
t295 = (t30 * (rSges(6,3) - t283) - t31 * t485 + t58 * t306 + (t59 * (-t229 - t488) - t58 * t283) * t286) * t277 + (t30 * t488 + t59 * (-t283 * t286 + t306) - t31 * rSges(6,3) + t58 * (-qJD(3) - t395)) * t278;
t187 = t219 * t278;
t186 = t219 * t277;
t150 = t192 * t286 + t272;
t149 = -t191 * t286 - t401;
t132 = t277 * t144;
t117 = -t407 + t273 + (t170 + t225) * t286;
t12 = t144 * t190 - t145 * t189 + t221 * t98 - t442 * t97 + ((t120 * t286 - t100) * t278 + (t101 + t118) * t277) * qJD(4);
t1 = [m(3) * (t149 * (t226 + t282) + t150 * (t281 + t522) + (-t175 + t192 + t273) * t174) + t297 + (t30 * (t282 + t425) + t59 * (-t354 + t422) + t31 * (t281 + t387) + t295 + (t59 + t519) * t58) * m(6) + (t52 * (t282 + t384) + t68 * (-t354 + t420) + t53 * (t281 + t385) + t296 + (t68 + t520) * t67) * m(5) + (t79 * (t282 + t352) + t117 * (t307 - t393) + t80 * (t281 + t335) + t298 + (t117 + t521) * t116) * m(4); t297 + (t30 * t425 + t31 * t387 + t360 * t482 + t295 + (t264 + t422 + t332) * t59 + t519 * t58) * m(6) + (t52 * t384 + t53 * t385 + t389 * t481 + t296 + (t264 + t420 - t351) * t68 + t520 * t67) * m(5) + (t80 * t335 + t79 * t352 + t298 + (-t264 + t307 + t528) * t117 + t521 * t116) * m(4) + (t149 * t226 + t150 * t522 + t174 * t192 - t175 * t191 - (t174 * t226 - t175 * t522) * t286) * m(3); m(4) * (-t277 * t80 - t278 * t79) + m(5) * (-t277 * t53 - t278 * t52) + m(6) * (-t277 * t31 - t278 * t30); t309 + (t538 * t278 + (t103 * t286 - t49) * t277) * t493 + ((t181 * t404 - t452) * t277 + (t533 + (-t277 * t453 - t536) * qJD(4)) * t278) * t375 + ((-t403 * t453 - t452) * t278 + (-t533 + (t278 * t181 + t536) * qJD(4)) * t277) * t373 + ((-t274 * t418 + t275 * t417) * t286 + (-t274 * t507 + t275 * t311) * qJD(4)) * t494 + (t286 * t54 + ((-t315 * t277 - t278 * t514 - t286 * t73) * t278 + (-t316 * t277 + t278 * t510 + t286 * t72) * t277) * t531) * t496 + (t286 * t55 + ((t277 * t514 - t315 * t278 - t286 * t71) * t278 + (-t277 * t510 - t316 * t278 + t286 * t70) * t277) * t531) * t495 + (t331 + t34) * t378 + (t330 + t35) * t377 + (t12 * t132 + t57 * t391 + (t12 * t120 + t57 * t101 - t30 * t371 + t59 * (-pkin(4) * t405 - t179) + (t57 * t121 - t371 * t58) * t286) * t277 + (t12 * t439 + t57 * t479 + t31 * t371 + t58 * t179 + (t57 * t120 - t371 * t59) * t286) * t278 - (-t59 * t380 + ((-t277 * t58 - t278 * t59) * t286 + t57 * (-t277 ^ 2 - t278 ^ 2) * qJD(4)) * t274) * pkin(4) + t517) * m(6) + (-(-t186 * t67 - t187 * t68) * t286 - (t102 * (-t186 * t277 - t187 * t278) + t347 * t348) * qJD(4) + 0.2e1 * t102 * ((t159 * t286 - t114) * t278 + (t115 - t135) * t277) + t347 * t199 + ((t53 - t481) * t278 + (-t286 * t67 - t52) * t277) * t219) * m(5); t309 + (t12 * (t145 * t278 + t132) + t57 * (-t278 * t97 + t391) + (-t277 * t59 + t278 * t58) * t179 + ((t31 - t482) * t278 + (-t286 * t58 - t30) * t277) * t206 + t517) * m(6);];
tauc = t1(:);

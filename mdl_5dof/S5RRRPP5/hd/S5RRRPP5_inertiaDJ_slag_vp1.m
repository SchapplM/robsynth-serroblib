% Calculate time derivative of joint inertia matrix for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:36
% EndTime: 2019-12-31 20:57:57
% DurationCPUTime: 12.27s
% Computational Cost: add. (9596->617), mult. (12506->807), div. (0->0), fcn. (9738->6), ass. (0->340)
t265 = qJ(2) + qJ(3);
t257 = cos(t265);
t256 = sin(t265);
t424 = Icges(5,5) * t256;
t211 = -Icges(5,3) * t257 + t424;
t428 = Icges(6,4) * t256;
t213 = -Icges(6,2) * t257 + t428;
t431 = Icges(4,4) * t256;
t215 = Icges(4,2) * t257 + t431;
t497 = t211 + t213 - t215;
t427 = Icges(6,4) * t257;
t216 = Icges(6,1) * t256 - t427;
t423 = Icges(5,5) * t257;
t217 = Icges(5,1) * t256 - t423;
t430 = Icges(4,4) * t257;
t218 = Icges(4,1) * t256 + t430;
t496 = t217 + t216 + t218;
t450 = rSges(6,1) + pkin(4);
t267 = sin(qJ(1));
t269 = cos(qJ(1));
t317 = Icges(6,5) * t257 + Icges(6,6) * t256;
t154 = -Icges(6,3) * t267 + t269 * t317;
t319 = Icges(4,5) * t257 - Icges(4,6) * t256;
t158 = Icges(4,3) * t267 + t269 * t319;
t322 = Icges(5,4) * t257 + Icges(5,6) * t256;
t162 = Icges(5,2) * t267 + t269 * t322;
t495 = -t154 + t158 + t162;
t318 = Icges(5,3) * t256 + t423;
t156 = Icges(5,6) * t267 + t269 * t318;
t321 = Icges(6,2) * t256 + t427;
t160 = -Icges(6,6) * t267 + t269 * t321;
t323 = -Icges(4,2) * t256 + t430;
t164 = Icges(4,6) * t267 + t269 * t323;
t494 = t156 + t160 - t164;
t326 = Icges(6,1) * t257 + t428;
t166 = -Icges(6,5) * t267 + t269 * t326;
t327 = Icges(5,1) * t257 + t424;
t168 = Icges(5,4) * t267 + t269 * t327;
t328 = Icges(4,1) * t257 - t431;
t170 = Icges(4,5) * t267 + t269 * t328;
t493 = t166 + t168 + t170;
t262 = qJD(2) + qJD(3);
t492 = (-t318 - t321 + t323) * t262;
t491 = (t326 + t327 + t328) * t262;
t210 = Icges(6,5) * t256 - Icges(6,6) * t257;
t212 = Icges(4,5) * t256 + Icges(4,6) * t257;
t214 = Icges(5,4) * t256 - Icges(5,6) * t257;
t481 = t210 - t212 - t214;
t479 = t497 * t256 + t496 * t257;
t399 = t257 * t269;
t401 = t256 * t269;
t490 = rSges(6,2) * t401 + t450 * t399;
t332 = rSges(6,1) * t257 + rSges(6,2) * t256;
t172 = t269 * rSges(6,3) + t267 * t332;
t489 = t267 * t257 * pkin(4) + t269 * qJ(5) + t172;
t438 = t267 * rSges(6,3);
t488 = -t267 * qJ(5) - t438 + t490;
t397 = t262 * t269;
t361 = t257 * t397;
t487 = rSges(6,2) * t361 - qJD(5) * t267;
t270 = -pkin(7) - pkin(6);
t266 = sin(qJ(2));
t375 = qJD(2) * t266;
t366 = pkin(2) * t375;
t486 = qJD(1) * t270 + t366;
t153 = Icges(6,3) * t269 + t267 * t317;
t157 = -Icges(4,3) * t269 + t267 * t319;
t161 = -Icges(5,2) * t269 + t267 * t322;
t155 = -Icges(5,6) * t269 + t267 * t318;
t167 = -Icges(5,4) * t269 + t267 * t327;
t314 = t155 * t256 + t167 * t257;
t466 = t269 * t314;
t159 = Icges(6,6) * t269 + t267 * t321;
t165 = Icges(6,5) * t269 + t267 * t326;
t312 = t159 * t256 + t165 * t257;
t467 = t269 * t312;
t163 = -Icges(4,6) * t269 + t267 * t323;
t169 = -Icges(4,5) * t269 + t267 * t328;
t310 = t163 * t256 - t169 * t257;
t468 = t269 * t310;
t485 = -t466 - t467 + t468 + (t153 - t157 - t161) * t267;
t309 = t164 * t256 - t170 * t257;
t311 = t160 * t256 + t166 * t257;
t313 = t156 * t256 + t168 * t257;
t484 = (-t309 + t311 + t313) * t269 + t495 * t267;
t483 = -t155 - t159 + t163;
t482 = t165 + t167 + t169;
t403 = t218 * t262;
t404 = t217 * t262;
t405 = t216 * t262;
t406 = t215 * t262;
t407 = t213 * t262;
t408 = t211 * t262;
t480 = (t408 - t406 + t407 + t491) * t257 + (-t404 - t403 - t405 - t492) * t256 - t481 * qJD(1);
t105 = -qJD(1) * t165 - t269 * t405;
t107 = -qJD(1) * t167 - t269 * t404;
t109 = -qJD(1) * t169 - t269 * t403;
t478 = t494 * t262 + t105 + t107 + t109;
t103 = -qJD(1) * t163 - t269 * t406;
t95 = -qJD(1) * t155 - t269 * t408;
t99 = -qJD(1) * t159 - t269 * t407;
t477 = -t493 * t262 - t103 + t95 + t99;
t476 = (t317 - t319 - t322) * t262 + t479 * qJD(1);
t268 = cos(qJ(2));
t236 = rSges(3,1) * t266 + rSges(3,2) * t268;
t297 = qJD(2) * t236;
t475 = t267 * t297;
t432 = Icges(3,4) * t268;
t325 = -Icges(3,2) * t266 + t432;
t195 = Icges(3,6) * t267 + t269 * t325;
t433 = Icges(3,4) * t266;
t330 = Icges(3,1) * t268 - t433;
t197 = Icges(3,5) * t267 + t269 * t330;
t307 = t195 * t266 - t197 * t268;
t474 = t267 * t307;
t473 = t267 * t309;
t472 = t267 * t311;
t471 = t267 * t313;
t252 = pkin(2) * t268 + pkin(1);
t445 = pkin(1) - t252;
t470 = t267 * t445;
t194 = -Icges(3,6) * t269 + t267 * t325;
t196 = -Icges(3,5) * t269 + t267 * t330;
t308 = t194 * t266 - t196 * t268;
t469 = t269 * t308;
t465 = qJD(1) * t153;
t464 = qJD(1) * t157;
t463 = qJD(1) * t161;
t320 = Icges(3,5) * t268 - Icges(3,6) * t266;
t192 = -Icges(3,3) * t269 + t267 * t320;
t459 = 2 * m(3);
t458 = 2 * m(4);
t457 = 2 * m(5);
t456 = 2 * m(6);
t263 = t267 ^ 2;
t264 = t269 ^ 2;
t455 = m(5) / 0.2e1;
t454 = m(6) / 0.2e1;
t453 = t267 / 0.2e1;
t452 = -t269 / 0.2e1;
t451 = -rSges(5,1) - pkin(3);
t449 = m(3) * t236;
t222 = rSges(4,1) * t256 + rSges(4,2) * t257;
t448 = m(4) * t222;
t447 = pkin(2) * t266;
t446 = t267 * pkin(6);
t261 = t269 * pkin(6);
t444 = -pkin(6) - t270;
t443 = rSges(3,1) * t268;
t442 = rSges(4,1) * t257;
t441 = rSges(5,1) * t256;
t440 = rSges(3,2) * t266;
t439 = rSges(3,3) * t269;
t260 = t267 * rSges(5,2);
t259 = t267 * rSges(3,3);
t258 = t267 * rSges(4,3);
t436 = -rSges(6,2) - qJ(4);
t435 = -rSges(5,3) - qJ(4);
t434 = rSges(6,3) + qJ(5);
t334 = -rSges(4,2) * t256 + t442;
t191 = t334 * t262;
t413 = t191 * t267;
t412 = t194 * t268;
t411 = t195 * t268;
t410 = t196 * t266;
t409 = t197 * t266;
t402 = t256 * t262;
t400 = t257 * t262;
t398 = t262 * t267;
t396 = t269 * t270;
t151 = t261 + t396 - t470;
t239 = t269 * t252;
t152 = -t269 * pkin(1) + t267 * t444 + t239;
t395 = t267 * t151 + t269 * t152;
t142 = qJ(4) * t402 + (pkin(3) * t262 - qJD(4)) * t257;
t333 = rSges(5,1) * t257 + rSges(5,3) * t256;
t394 = -t333 * t262 - t142;
t174 = -t269 * rSges(4,3) + t267 * t334;
t177 = rSges(4,1) * t399 - rSges(4,2) * t401 + t258;
t92 = t267 * t174 + t269 * t177;
t176 = rSges(5,1) * t399 + rSges(5,3) * t401 + t260;
t202 = pkin(3) * t399 + qJ(4) * t401;
t393 = -t176 - t202;
t201 = (pkin(3) * t257 + qJ(4) * t256) * t267;
t392 = t267 * t201 + t269 * t202;
t219 = pkin(3) * t256 - qJ(4) * t257;
t378 = qJD(1) * t267;
t203 = t219 * t378;
t221 = -rSges(5,3) * t257 + t441;
t391 = t221 * t378 + t203;
t390 = -t219 - t221;
t373 = qJD(4) * t256;
t389 = qJ(4) * t361 + t269 * t373;
t377 = qJD(1) * t269;
t388 = rSges(5,2) * t377 + rSges(5,3) * t361;
t356 = t256 * t378;
t387 = rSges(4,2) * t356 + rSges(4,3) * t377;
t385 = t486 * t267;
t384 = t269 * t443 + t259;
t383 = t263 + t264;
t382 = qJD(1) * t154;
t381 = qJD(1) * t158;
t380 = qJD(1) * t162;
t193 = Icges(3,3) * t267 + t269 * t320;
t379 = qJD(1) * t193;
t374 = qJD(2) * t268;
t371 = qJD(5) * t269;
t370 = -pkin(3) - t450;
t362 = t256 * t397;
t282 = -t257 * t378 - t362;
t298 = t222 * t262;
t369 = t174 * t377 + t267 * (-t267 * t298 + (t269 * t334 + t258) * qJD(1)) + t269 * (rSges(4,1) * t282 - rSges(4,2) * t361 + t387);
t363 = t256 * t398;
t230 = pkin(3) * t363;
t355 = t257 * t377;
t368 = t201 * t377 + t267 * (pkin(3) * t355 + t267 * t373 - t230 + (t256 * t377 + t257 * t398) * qJ(4)) + t269 * (pkin(3) * t282 - qJ(4) * t356 + t389);
t367 = t269 * t440;
t365 = pkin(2) * t374;
t364 = -t270 - t434;
t360 = t267 * ((-t269 * t445 - t446) * qJD(1) - t385) + t269 * (-t269 * t366 + (t269 * t444 + t470) * qJD(1)) + t151 * t377;
t359 = -t202 - t488;
t220 = rSges(6,1) * t256 - rSges(6,2) * t257;
t358 = pkin(4) * t356 + t220 * t378 + t203;
t357 = t230 + t385;
t354 = t266 * t378;
t351 = -t222 - t447;
t96 = qJD(1) * t156 - t267 * t408;
t349 = t167 * t262 - t96;
t133 = t390 * t269;
t100 = qJD(1) * t160 - t267 * t407;
t347 = t165 * t262 - t100;
t104 = qJD(1) * t164 - t267 * t406;
t345 = t169 * t262 + t104;
t106 = qJD(1) * t166 - t267 * t405;
t343 = t159 * t262 + t106;
t108 = qJD(1) * t168 - t267 * t404;
t341 = t155 * t262 + t108;
t110 = qJD(1) * t170 - t267 * t403;
t339 = t163 * t262 - t110;
t338 = -t267 * t270 + t239;
t173 = -t269 * rSges(5,2) + t267 * t333;
t45 = t267 * t173 + t269 * t176 + t392;
t337 = t390 - t447;
t336 = -pkin(4) * t256 - t219 - t220;
t335 = -t440 + t443;
t278 = t256 * t436 + t257 * t370 - t252;
t66 = t267 * t278 + t269 * t364;
t67 = t267 * t364 + t202 + t239 + t490;
t331 = t267 * t67 + t269 * t66;
t329 = Icges(3,1) * t266 + t432;
t324 = Icges(3,2) * t268 + t433;
t303 = -pkin(4) * t400 - t332 * t262 - t142;
t302 = t173 * t377 + t267 * (-t221 * t398 + (t269 * t333 + t260) * qJD(1)) + t269 * (rSges(5,1) * t282 - rSges(5,3) * t356 + t388) + t368;
t301 = -t365 + t394;
t300 = -pkin(1) - t335;
t127 = t337 * t269;
t117 = t336 * t269;
t299 = -t252 - t334;
t37 = t489 * t267 + t488 * t269 + t392;
t296 = t336 - t447;
t291 = t262 * t214;
t289 = t262 * t212;
t287 = t262 * t210;
t286 = qJD(2) * t329;
t285 = qJD(2) * t324;
t284 = qJD(2) * (-Icges(3,5) * t266 - Icges(3,6) * t268);
t101 = -t269 * t291 - t463;
t102 = -t267 * t291 + t380;
t48 = t153 * t269 + t267 * t312;
t49 = t154 * t269 + t472;
t50 = -t161 * t269 + t267 * t314;
t51 = -t162 * t269 + t471;
t52 = -t157 * t269 - t267 * t310;
t53 = -t158 * t269 - t473;
t93 = -t269 * t287 - t465;
t94 = -t267 * t287 + t382;
t97 = -t269 * t289 - t464;
t98 = -t267 * t289 + t381;
t283 = (((-t93 + t101 + t97) * t267 + (-t472 - t471 + t473 - t485) * qJD(1)) * t267 + (t49 + t51 + t53) * t378 + t484 * t377) * t267 + ((-t48 - t50 - t52) * t378 + t485 * t377 + (t484 * qJD(1) + (-t463 - t464 + t465 + t482 * t402 + t483 * t400 + (-t106 - t108 - t110) * t257 + (-t100 + t104 - t96) * t256) * t269 + (-t102 + t94 - t98 + t478 * t257 + t477 * t256 + (t312 + t314 - t310 + t495) * qJD(1)) * t267) * t267) * t269;
t114 = t296 * t269;
t281 = t368 + t489 * t377 + (-rSges(6,1) * t362 + pkin(4) * t282 - qJ(5) * t377 - qJD(1) * t172 + t487) * t269 + (-qJ(5) * t378 + t371 + (t355 - t363) * pkin(4) - t220 * t398 + (t269 * t332 - t438) * qJD(1)) * t267;
t280 = t256 * t435 + t257 * t451 - t252;
t279 = t303 - t365;
t277 = rSges(3,2) * t354 + rSges(3,3) * t377 - t269 * t297;
t276 = t280 * t267;
t4 = (-t269 * t94 + (t49 - t467) * qJD(1)) * t269 + (t48 * qJD(1) + (t105 * t257 + t160 * t400 - t166 * t402 + t256 * t99 - t382) * t267 + (t93 - t343 * t257 + t347 * t256 + (t153 + t311) * qJD(1)) * t269) * t267;
t5 = (t269 * t102 + (t51 - t466) * qJD(1)) * t269 + (t50 * qJD(1) + (t107 * t257 + t156 * t400 - t168 * t402 + t256 * t95 + t380) * t267 + (-t101 - t341 * t257 + t349 * t256 + (-t161 + t313) * qJD(1)) * t269) * t267;
t6 = (t269 * t98 + (t53 + t468) * qJD(1)) * t269 + (t52 * qJD(1) + (-t103 * t256 + t109 * t257 - t164 * t400 - t170 * t402 + t381) * t267 + (-t97 + t339 * t257 + t345 * t256 + (-t157 - t309) * qJD(1)) * t269) * t267;
t275 = (-t4 - t5 - t6) * t269 + t283;
t271 = (t478 * t256 - t477 * t257 - t476 * t267 + t480 * t269) * t453 + (t476 * t269 + t480 * t267 + (t345 + t347 + t349) * t257 + (-t339 + t341 + t343) * t256) * t452 + (t482 * t256 + t483 * t257 + t479 * t267 + t481 * t269) * t378 / 0.2e1 + (t493 * t256 - t494 * t257 - t481 * t267 + t479 * t269) * t377 / 0.2e1;
t249 = pkin(2) * t354;
t229 = t335 * qJD(2);
t199 = -t367 + t384;
t198 = t267 * t335 - t439;
t150 = t351 * t269;
t149 = t351 * t267;
t139 = t446 + (pkin(1) - t440) * t269 + t384;
t138 = t267 * t300 + t261 + t439;
t132 = t390 * t267;
t129 = t177 + t338;
t128 = (rSges(4,3) - t270) * t269 + t299 * t267;
t126 = t337 * t267;
t121 = t267 * t284 + t379;
t120 = -qJD(1) * t192 + t269 * t284;
t116 = t336 * t267;
t113 = t296 * t267;
t112 = t475 + ((-rSges(3,3) - pkin(6)) * t267 + t300 * t269) * qJD(1);
t111 = (t261 + (-pkin(1) - t443) * t267) * qJD(1) + t277;
t85 = t338 - t393;
t84 = (rSges(5,2) - t270) * t269 + t276;
t81 = -t222 * t377 - t413 + (-t266 * t377 - t267 * t374) * pkin(2);
t80 = t222 * t378 + t249 + (-t191 - t365) * t269;
t65 = t193 * t267 - t307 * t269;
t64 = t192 * t267 - t469;
t63 = -t193 * t269 - t474;
t62 = -t192 * t269 - t267 * t308;
t61 = t222 * t398 + (t269 * t299 - t258) * qJD(1) + t385;
t60 = (-t252 - t442) * t378 + (-t298 - t486) * t269 + t387;
t47 = qJD(1) * t133 + t267 * t394;
t46 = t269 * t394 + t391;
t44 = qJD(1) * t127 + t267 * t301;
t43 = t269 * t301 + t249 + t391;
t42 = t92 + t395;
t41 = qJD(1) * t117 + t267 * t303;
t40 = t269 * t303 + t358;
t39 = qJD(1) * t114 + t267 * t279;
t38 = t269 * t279 + t249 + t358;
t36 = (-t373 + (t257 * t435 + t441) * t262) * t267 + (t269 * t280 - t260) * qJD(1) + t357;
t35 = (t402 * t451 - t366) * t269 + (t276 - t396) * qJD(1) + t388 + t389;
t34 = t45 + t395;
t33 = -t371 + (-t373 + (t256 * t450 + t257 * t436) * t262) * t267 + (t267 * t434 + t269 * t278) * qJD(1) + t357;
t32 = (t370 * t402 - t366) * t269 + t66 * qJD(1) + t389 + t487;
t31 = -t177 * t378 + t369;
t24 = t37 + t395;
t11 = (-t152 - t177) * t378 + t360 + t369;
t10 = t378 * t393 + t302;
t9 = t359 * t378 + t281;
t8 = (-t152 + t393) * t378 + t302 + t360;
t7 = (-t152 + t359) * t378 + t281 + t360;
t1 = [(t32 * t67 + t33 * t66) * t456 + (t35 * t85 + t36 * t84) * t457 + (t128 * t61 + t129 * t60) * t458 + (t111 * t139 + t112 * t138) * t459 + (t330 - t324) * t375 + (t329 + t325) * t374 + t497 * t402 + t496 * t400 + t492 * t257 + t491 * t256; t271 + m(3) * ((-t111 * t267 - t112 * t269) * t236 + (-t138 * t269 - t139 * t267) * t229) + m(6) * (t113 * t32 + t114 * t33 + t38 * t66 + t39 * t67) + m(5) * (t126 * t35 + t127 * t36 + t43 * t84 + t44 * t85) + m(4) * (t128 * t80 + t129 * t81 + t149 * t60 + t150 * t61) + (-qJD(2) * t308 + (qJD(1) * t195 - t267 * t285) * t268 + (qJD(1) * t197 - t267 * t286) * t266) * t452 + (-qJD(2) * t307 + (-qJD(1) * t194 - t269 * t285) * t268 + (-qJD(1) * t196 - t269 * t286) * t266) * t453 + (t263 / 0.2e1 + t264 / 0.2e1) * t320 * qJD(2) + ((t411 / 0.2e1 + t409 / 0.2e1 - t139 * t449) * t269 + (t138 * t449 + t412 / 0.2e1 + t410 / 0.2e1) * t267) * qJD(1); ((t198 * t267 + t199 * t269) * ((qJD(1) * t198 + t277) * t269 + (-t475 + (-t199 - t367 + t259) * qJD(1)) * t267) + t383 * t236 * t229) * t459 + t267 * ((t267 * t120 + (t64 + t474) * qJD(1)) * t267 + (t65 * qJD(1) + (t194 * t374 + t196 * t375) * t269 + (-t121 + (-t409 - t411) * qJD(2) + (t193 - t308) * qJD(1)) * t267) * t269) - t269 * ((t269 * t121 + (t63 + t469) * qJD(1)) * t269 + (t62 * qJD(1) + (-t195 * t374 - t197 * t375 + t379) * t267 + (-t120 + (t410 + t412) * qJD(2) - t307 * qJD(1)) * t269) * t267) - t269 * t4 - t269 * t5 - t269 * t6 + (t267 * t63 - t269 * t62) * t378 + (t267 * t65 - t269 * t64) * t377 + t283 + (t113 * t39 + t114 * t38 + t24 * t7) * t456 + (t126 * t44 + t127 * t43 + t34 * t8) * t457 + (t11 * t42 + t149 * t81 + t150 * t80) * t458; m(6) * (t116 * t32 + t117 * t33 + t40 * t66 + t41 * t67) + m(5) * (t132 * t35 + t133 * t36 + t46 * t84 + t47 * t85) + (-t267 * t60 - t269 * t61 + (t128 * t267 - t129 * t269) * qJD(1)) * t448 + m(4) * (-t128 * t269 - t129 * t267) * t191 + t271; m(4) * (-t150 * t191 * t269 + t11 * t92 - t149 * t413 + t31 * t42) + m(6) * (t113 * t41 + t114 * t40 + t116 * t39 + t117 * t38 + t24 * t9 + t37 * t7) + m(5) * (t10 * t34 + t126 * t47 + t127 * t46 + t132 * t44 + t133 * t43 + t45 * t8) + (-t267 * t81 - t269 * t80 + (-t149 * t269 + t150 * t267) * qJD(1)) * t448 + t275; (t116 * t41 + t117 * t40 + t37 * t9) * t456 + (t10 * t45 + t132 * t47 + t133 * t46) * t457 + (t191 * t222 * t383 + t31 * t92) * t458 + t275; 0.2e1 * (t331 * t454 + (t267 * t85 + t269 * t84) * t455) * t400 + 0.2e1 * ((t267 * t32 + t269 * t33 + t377 * t67 - t378 * t66) * t454 + (t267 * t35 + t269 * t36 + t377 * t85 - t378 * t84) * t455) * t256; 0.2e1 * ((t113 * t398 + t114 * t397 - t7) * t454 + (t126 * t398 + t127 * t397 - t8) * t455) * t257 + 0.2e1 * ((t113 * t377 - t114 * t378 + t24 * t262 + t267 * t39 + t269 * t38) * t454 + (t126 * t377 - t127 * t378 + t262 * t34 + t267 * t44 + t269 * t43) * t455) * t256; 0.2e1 * ((t116 * t398 + t117 * t397 - t9) * t454 + (t132 * t398 + t133 * t397 - t10) * t455) * t257 + 0.2e1 * ((t116 * t377 - t117 * t378 + t262 * t37 + t267 * t41 + t269 * t40) * t454 + (t132 * t377 - t133 * t378 + t262 * t45 + t267 * t47 + t269 * t46) * t455) * t256; 0.4e1 * (t455 + t454) * (-0.1e1 + t383) * t256 * t400; m(6) * (-qJD(1) * t331 - t267 * t33 + t269 * t32); m(6) * (-t267 * t38 + t269 * t39 + (-t113 * t267 - t114 * t269) * qJD(1)); m(6) * (-t267 * t40 + t269 * t41 + (-t116 * t267 - t117 * t269) * qJD(1)); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

% Calculate time derivative of joint inertia matrix for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:38:00
% EndTime: 2019-03-09 02:38:27
% DurationCPUTime: 15.22s
% Computational Cost: add. (23381->760), mult. (19008->1053), div. (0->0), fcn. (17463->12), ass. (0->375)
t497 = Icges(4,3) + Icges(5,3);
t259 = qJ(3) + pkin(10);
t252 = sin(t259);
t255 = cos(t259);
t265 = sin(qJ(3));
t267 = cos(qJ(3));
t496 = Icges(4,5) * t267 + Icges(5,5) * t255 - Icges(4,6) * t265 - Icges(5,6) * t252;
t260 = qJ(1) + pkin(9);
t253 = sin(t260);
t256 = cos(t260);
t492 = t253 * t497 + t496 * t256;
t262 = cos(pkin(11));
t261 = sin(pkin(11));
t408 = t256 * t261;
t199 = t253 * t262 - t255 * t408;
t407 = t256 * t262;
t414 = t253 * t261;
t200 = t255 * t407 + t414;
t416 = t252 * t256;
t108 = Icges(6,5) * t200 + Icges(6,6) * t199 + Icges(6,3) * t416;
t434 = Icges(4,4) * t267;
t320 = -Icges(4,2) * t265 + t434;
t187 = Icges(4,6) * t253 + t256 * t320;
t435 = Icges(4,4) * t265;
t326 = Icges(4,1) * t267 - t435;
t189 = Icges(4,5) * t253 + t256 * t326;
t303 = t187 * t265 - t189 * t267;
t432 = Icges(5,4) * t255;
t318 = -Icges(5,2) * t252 + t432;
t161 = Icges(5,6) * t253 + t256 * t318;
t433 = Icges(5,4) * t252;
t324 = Icges(5,1) * t255 - t433;
t164 = Icges(5,5) * t253 + t256 * t324;
t305 = t161 * t252 - t164 * t255;
t474 = t303 + t305;
t495 = -t108 * t416 + t474 * t256;
t493 = t496 * t253 - t256 * t497;
t491 = (-Icges(4,5) * t265 - Icges(5,5) * t252 - Icges(4,6) * t267 - Icges(5,6) * t255) * qJD(3);
t186 = -Icges(4,6) * t256 + t253 * t320;
t188 = -Icges(4,5) * t256 + t253 * t326;
t304 = t186 * t265 - t188 * t267;
t160 = -Icges(5,6) * t256 + t253 * t318;
t163 = -Icges(5,5) * t256 + t253 * t324;
t306 = t160 * t252 - t163 * t255;
t490 = t304 + t306;
t459 = t253 / 0.2e1;
t489 = -t256 / 0.2e1;
t488 = -qJD(1) / 0.2e1;
t487 = t492 * qJD(1);
t410 = t255 * t262;
t198 = t253 * t410 - t408;
t411 = t255 * t261;
t292 = t253 * t411 + t407;
t417 = t252 * t253;
t109 = Icges(6,4) * t198 - Icges(6,2) * t292 + Icges(6,6) * t417;
t110 = Icges(6,4) * t200 + Icges(6,2) * t199 + Icges(6,6) * t416;
t111 = Icges(6,1) * t198 - Icges(6,4) * t292 + Icges(6,5) * t417;
t112 = Icges(6,1) * t200 + Icges(6,4) * t199 + Icges(6,5) * t416;
t285 = t303 * t253;
t286 = t304 * t256;
t287 = t305 * t253;
t288 = t306 * t256;
t370 = t108 * t417;
t107 = Icges(6,5) * t198 - Icges(6,6) * t292 + Icges(6,3) * t417;
t371 = t107 * t416;
t486 = t109 * t199 - t110 * t292 + t111 * t200 + t112 * t198 + t253 * t493 - t256 * t492 - t285 - t286 - t287 - t288 + t370 + t371;
t248 = pkin(3) * t267 + pkin(2);
t228 = t256 * t248;
t263 = -qJ(4) - pkin(7);
t444 = -pkin(7) - t263;
t152 = -pkin(2) * t256 + t253 * t444 + t228;
t442 = rSges(5,1) * t255;
t344 = -rSges(5,2) * t252 + t442;
t168 = -rSges(5,3) * t256 + t253 * t344;
t412 = t255 * t256;
t244 = t253 * rSges(5,3);
t476 = -rSges(5,2) * t416 + t244;
t169 = rSges(5,1) * t412 + t476;
t212 = rSges(5,1) * t252 + rSges(5,2) * t255;
t289 = qJD(3) * t212;
t246 = t256 * pkin(7);
t406 = t256 * t263;
t446 = pkin(2) - t248;
t482 = t253 * t446;
t151 = t246 + t406 - t482;
t242 = qJD(4) * t253;
t384 = qJD(3) * t265;
t376 = pkin(3) * t384;
t392 = qJD(1) * t253;
t365 = qJD(4) * t256 + t253 * t376 + t263 * t392;
t390 = qJD(1) * t256;
t448 = pkin(7) * t253;
t367 = t253 * ((-t256 * t446 - t448) * qJD(1) - t365) + t256 * (-t256 * t376 + t242 + (t256 * t444 + t482) * qJD(1)) + t151 * t390;
t364 = t252 * t392;
t398 = rSges(5,2) * t364 + rSges(5,3) * t390;
t24 = (qJD(1) * t168 - t256 * t289 + t398) * t256 + (-t253 * t289 + (-t152 - t169 + t476) * qJD(1)) * t253 + t367;
t466 = 2 * m(5);
t485 = t24 * t466;
t441 = rSges(4,2) * t265;
t443 = rSges(4,1) * t267;
t345 = -t441 + t443;
t440 = rSges(4,3) * t256;
t190 = t253 * t345 - t440;
t378 = t256 * t441;
t245 = t253 * rSges(4,3);
t396 = t256 * t443 + t245;
t191 = -t378 + t396;
t237 = rSges(4,1) * t265 + rSges(4,2) * t267;
t290 = qJD(3) * t237;
t389 = qJD(1) * t265;
t363 = t253 * t389;
t271 = rSges(4,2) * t363 + rSges(4,3) * t390 - t256 * t290;
t42 = (qJD(1) * t190 + t271) * t256 + (-t253 * t290 + (-t191 - t378 + t245) * qJD(1)) * t253;
t467 = 2 * m(4);
t484 = t42 * t467;
t247 = pkin(5) * t262 + pkin(4);
t445 = pkin(4) - t247;
t483 = t252 * t445;
t264 = -pkin(8) - qJ(5);
t405 = qJ(5) + t264;
t481 = t255 * t405;
t258 = pkin(11) + qJ(6);
t251 = sin(t258);
t254 = cos(t258);
t182 = -t251 * t412 + t253 * t254;
t183 = t251 * t253 + t254 * t412;
t100 = t183 * rSges(7,1) + t182 * rSges(7,2) + rSges(7,3) * t416;
t238 = pkin(5) * t414;
t477 = -t264 * t416 + t238;
t480 = t247 * t412 + t100 + t477;
t479 = -qJD(1) * t493 + t256 * t491;
t478 = -t253 * t491 - t487;
t452 = sin(qJ(1)) * pkin(1);
t475 = t246 - t452;
t471 = t110 * t199 + t112 * t200 + t253 * t492 - t495;
t372 = t107 * t417;
t470 = t109 * t292 - t111 * t198 + t253 * t490 + t256 * t493 - t372;
t391 = qJD(1) * t255;
t352 = -qJD(6) + t391;
t385 = qJD(3) * t256;
t358 = t252 * t385;
t468 = t253 * t352 + t358;
t465 = 2 * m(6);
t464 = 2 * m(7);
t249 = t253 ^ 2;
t250 = t256 ^ 2;
t463 = m(6) / 0.2e1;
t462 = m(7) / 0.2e1;
t316 = Icges(6,4) * t262 - Icges(6,2) * t261;
t177 = -Icges(6,6) * t255 + t252 * t316;
t461 = t177 / 0.2e1;
t322 = Icges(6,1) * t262 - Icges(6,4) * t261;
t178 = -Icges(6,5) * t255 + t252 * t322;
t460 = t178 / 0.2e1;
t458 = -t255 / 0.2e1;
t456 = t256 / 0.2e1;
t455 = -t261 / 0.2e1;
t454 = t262 / 0.2e1;
t453 = m(4) * t237;
t451 = pkin(3) * t265;
t450 = pkin(4) * t252;
t449 = pkin(4) * t255;
t257 = cos(qJ(1)) * pkin(1);
t447 = qJD(1) / 0.2e1;
t439 = rSges(7,3) * t252;
t438 = -rSges(6,3) - qJ(5);
t437 = -rSges(7,3) + t264;
t272 = -t252 * t405 - t255 * t445;
t380 = pkin(5) * t408;
t413 = t254 * t256;
t415 = t253 * t255;
t180 = -t251 * t415 - t413;
t181 = -t251 * t256 + t254 * t415;
t340 = -rSges(7,1) * t181 - rSges(7,2) * t180;
t99 = rSges(7,3) * t417 - t340;
t436 = t253 * t272 - t380 + t99;
t431 = Icges(7,4) * t251;
t430 = Icges(7,4) * t254;
t321 = Icges(7,1) * t254 - t431;
t162 = -Icges(7,5) * t255 + t252 * t321;
t423 = t162 * t254;
t422 = t186 * t267;
t421 = t187 * t267;
t420 = t188 * t265;
t419 = t189 * t265;
t418 = t247 * t252;
t409 = t255 * t264;
t227 = pkin(4) * t412;
t196 = qJ(5) * t416 + t227;
t404 = -t196 + t480;
t403 = t253 * t151 + t256 * t152;
t339 = rSges(7,1) * t254 - rSges(7,2) * t251;
t167 = -rSges(7,3) * t255 + t252 * t339;
t402 = t167 + t481 - t483;
t401 = -t152 - t196;
t211 = -qJ(5) * t255 + t450;
t232 = pkin(3) * t363;
t400 = t211 * t392 + t232;
t399 = qJD(1) * t380 + t264 * t364;
t382 = qJD(5) * t252;
t222 = t256 * t382;
t397 = t222 + t242;
t395 = t249 + t250;
t388 = qJD(3) * t252;
t387 = qJD(3) * t253;
t386 = qJD(3) * t255;
t383 = qJD(3) * t267;
t381 = qJD(6) * t252;
t359 = t255 * t385;
t353 = -qJD(6) * t255 + qJD(1);
t300 = t353 * t254;
t84 = t251 * t468 + t256 * t300;
t301 = t353 * t251;
t85 = -t254 * t468 + t256 * t301;
t379 = t85 * rSges(7,1) + t84 * rSges(7,2) + rSges(7,3) * t359;
t375 = pkin(3) * t383;
t95 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t417;
t97 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t417;
t336 = -t251 * t95 + t254 * t97;
t93 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t417;
t38 = t252 * t336 - t255 * t93;
t311 = Icges(7,5) * t254 - Icges(7,6) * t251;
t156 = -Icges(7,3) * t255 + t252 * t311;
t315 = -Icges(7,2) * t251 + t430;
t159 = -Icges(7,6) * t255 + t252 * t315;
t52 = t156 * t417 + t159 * t180 + t162 * t181;
t374 = t38 / 0.2e1 + t52 / 0.2e1;
t96 = Icges(7,4) * t183 + Icges(7,2) * t182 + Icges(7,6) * t416;
t98 = Icges(7,1) * t183 + Icges(7,4) * t182 + Icges(7,5) * t416;
t335 = -t251 * t96 + t254 * t98;
t94 = Icges(7,5) * t183 + Icges(7,6) * t182 + Icges(7,3) * t416;
t39 = t252 * t335 - t255 * t94;
t53 = t156 * t416 + t159 * t182 + t162 * t183;
t373 = t39 / 0.2e1 + t53 / 0.2e1;
t141 = qJD(1) * t292 + t261 * t358;
t142 = -qJD(1) * t198 - t262 * t358;
t366 = t142 * rSges(6,1) + t141 * rSges(6,2) + rSges(6,3) * t359;
t125 = t200 * rSges(6,1) + t199 * rSges(6,2) + rSges(6,3) * t416;
t362 = t159 * t386;
t361 = t252 * t387;
t360 = t253 * t386;
t357 = t386 / 0.2e1;
t356 = -t211 - t451;
t355 = -t212 - t451;
t354 = -t247 * t255 - t248;
t337 = qJ(5) * t252 + t449;
t195 = t337 * t253;
t351 = t253 * t195 + t256 * t196 + t403;
t341 = rSges(6,1) * t262 - rSges(6,2) * t261;
t179 = -rSges(6,3) * t255 + t252 * t341;
t350 = -t179 + t356;
t349 = -qJD(3) * t337 + qJD(5) * t255 - t375;
t86 = t253 * t300 + (-t256 * t352 + t361) * t251;
t87 = t352 * t413 + (-t254 * t388 + t301) * t253;
t348 = t87 * rSges(7,1) + t86 * rSges(7,2);
t347 = -t253 * t263 + t228 + t257;
t171 = t355 * t256;
t143 = qJD(1) * t199 + t261 * t361;
t144 = qJD(1) * t200 - t262 * t361;
t343 = -t144 * rSges(6,1) - t143 * rSges(6,2);
t342 = -rSges(6,1) * t198 + rSges(6,2) * t292;
t338 = -t406 - t452;
t27 = t180 * t95 + t181 * t97 + t417 * t93;
t28 = t180 * t96 + t181 * t98 + t417 * t94;
t16 = t253 * t28 - t256 * t27;
t334 = t253 * t27 + t256 * t28;
t29 = t182 * t95 + t183 * t97 + t416 * t93;
t30 = t182 * t96 + t183 * t98 + t416 * t94;
t17 = t253 * t30 - t256 * t29;
t333 = t253 * t29 + t256 * t30;
t332 = t253 * t39 - t256 * t38;
t331 = t253 * t38 + t256 * t39;
t274 = t252 * t437 + t354;
t56 = -t452 + (pkin(5) * t261 - t263) * t256 + t274 * t253 + t340;
t57 = t347 + t480;
t330 = t253 * t57 + t256 * t56;
t64 = t167 * t417 + t255 * t99;
t65 = -t100 * t255 - t167 * t416;
t329 = t253 * t65 + t256 * t64;
t277 = t252 * t438 - t248 - t449;
t269 = t253 * t277 + t338;
t72 = t269 + t342;
t73 = t347 + t125 + t196;
t328 = t253 * t73 + t256 * t72;
t327 = -t100 * t253 + t256 * t99;
t325 = Icges(4,1) * t265 + t434;
t323 = Icges(5,1) * t252 + t432;
t319 = Icges(4,2) * t267 + t435;
t317 = Icges(5,2) * t255 + t433;
t312 = Icges(6,5) * t262 - Icges(6,6) * t261;
t302 = t356 - t402;
t299 = -(rSges(6,3) * t252 + t255 * t341) * qJD(3) + t349;
t298 = -pkin(2) - t345;
t213 = qJ(5) * t359;
t217 = pkin(4) * t361;
t276 = t252 * t390 + t360;
t296 = t195 * t390 + t253 * (qJ(5) * t276 + qJD(1) * t227 + t253 * t382 - t217) + t256 * (-qJ(5) * t364 + t213 + t222 + (-t253 * t391 - t358) * pkin(4)) + t367;
t114 = t350 * t256;
t295 = -t248 - t344;
t126 = (-rSges(7,1) * t251 - rSges(7,2) * t254) * t381 + (t255 * t339 + t439) * qJD(3);
t291 = -t272 * qJD(3) - t126 + t349;
t284 = qJD(3) * t325;
t283 = qJD(3) * t323;
t282 = qJD(3) * t319;
t281 = qJD(3) * t317;
t77 = t302 * t256;
t275 = t359 - t364;
t115 = (-Icges(7,5) * t251 - Icges(7,6) * t254) * t381 + (Icges(7,3) * t252 + t255 * t311) * qJD(3);
t121 = (-Icges(7,1) * t251 - t430) * t381 + (Icges(7,5) * t252 + t255 * t321) * qJD(3);
t270 = -t255 * t115 + t156 * t388 + t386 * t423 + (t121 * t252 - t159 * t381) * t254;
t221 = t345 * qJD(3);
t205 = t344 * qJD(3);
t170 = t355 * t253;
t155 = (Icges(6,5) * t252 + t255 * t322) * qJD(3);
t154 = (Icges(6,6) * t252 + t255 * t316) * qJD(3);
t137 = t448 + t257 + (pkin(2) - t441) * t256 + t396;
t136 = t253 * t298 + t440 + t475;
t128 = t169 + t347;
t127 = -t452 + (rSges(5,3) - t263) * t256 + t295 * t253;
t124 = rSges(6,3) * t417 - t342;
t118 = (-Icges(7,2) * t254 - t431) * t381 + (Icges(7,6) * t252 + t255 * t315) * qJD(3);
t113 = t350 * t253;
t102 = -t212 * t390 - t205 * t253 + (-t253 * t383 - t256 * t389) * pkin(3);
t101 = t212 * t392 + t232 + (-t205 - t375) * t256;
t89 = t237 * t387 + (-t257 + (-rSges(4,3) - pkin(7)) * t253 + t298 * t256) * qJD(1);
t88 = ((-pkin(2) - t443) * t253 + t475) * qJD(1) + t271;
t76 = t302 * t253;
t75 = t212 * t387 + (t256 * t295 - t244 - t257) * qJD(1) + t365;
t74 = t242 + qJD(3) * t171 + ((-t248 - t442) * t253 + t338) * qJD(1) + t398;
t71 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t276;
t70 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t275;
t69 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t276;
t68 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t275;
t59 = -t156 * t255 + (-t159 * t251 + t423) * t252;
t58 = t59 * t388;
t55 = qJD(1) * t114 + t253 * t299;
t54 = t179 * t392 + t256 * t299 + t400;
t51 = t327 * t252;
t50 = rSges(7,3) * t276 + t348;
t49 = -rSges(7,3) * t364 + t379;
t48 = Icges(7,1) * t87 + Icges(7,4) * t86 + Icges(7,5) * t276;
t47 = Icges(7,1) * t85 + Icges(7,4) * t84 + Icges(7,5) * t275;
t46 = Icges(7,4) * t87 + Icges(7,2) * t86 + Icges(7,6) * t276;
t45 = Icges(7,4) * t85 + Icges(7,2) * t84 + Icges(7,6) * t275;
t44 = Icges(7,5) * t87 + Icges(7,6) * t86 + Icges(7,3) * t276;
t43 = Icges(7,5) * t85 + Icges(7,6) * t84 + Icges(7,3) * t275;
t41 = t217 + (t386 * t438 - t382) * t253 + (t256 * t277 - t257) * qJD(1) + t343 + t365;
t40 = t213 + (-t450 - t451) * t385 + t269 * qJD(1) + t366 + t397;
t33 = qJD(1) * t77 + t253 * t291;
t32 = t256 * t291 + t392 * t402 + t400;
t31 = t124 * t253 + t125 * t256 + t351;
t26 = (-t382 + (t255 * t437 + t418) * qJD(3)) * t253 + (t256 * t274 - t238 - t257) * qJD(1) - t348 + t365;
t25 = (-t409 - t418 - t451) * t385 + ((t354 - t439) * t253 + t338) * qJD(1) + t379 + t397 + t399;
t23 = (t167 * t387 + t50) * t255 + (-qJD(3) * t99 + t126 * t253 + t167 * t390) * t252;
t22 = (-t167 * t385 - t49) * t255 + (qJD(3) * t100 - t126 * t256 + t167 * t392) * t252;
t21 = (-t362 + (-qJD(6) * t162 - t118) * t252) * t251 + t270;
t20 = t253 * t436 + t256 * t404 + t351;
t19 = t115 * t417 + t118 * t180 + t121 * t181 + t156 * t276 + t159 * t86 + t162 * t87;
t18 = t115 * t416 + t118 * t182 + t121 * t183 + t156 * t275 + t159 * t84 + t162 * t85;
t15 = t327 * t386 + (-t253 * t49 + t256 * t50 + (-t100 * t256 - t253 * t99) * qJD(1)) * t252;
t14 = t252 * t333 - t255 * t53;
t13 = t252 * t334 - t255 * t52;
t12 = t253 * (rSges(6,3) * t360 - t343) + t256 * t366 + (t256 * t124 + (-t125 + t401) * t253) * qJD(1) + t296;
t11 = (qJD(3) * t335 - t43) * t255 + (qJD(3) * t94 - t251 * t45 + t254 * t47 + (-t251 * t98 - t254 * t96) * qJD(6)) * t252;
t10 = (qJD(3) * t336 - t44) * t255 + (qJD(3) * t93 - t251 * t46 + t254 * t48 + (-t251 * t97 - t254 * t95) * qJD(6)) * t252;
t9 = t94 * t360 + t180 * t45 + t181 * t47 + t86 * t96 + t87 * t98 + (t253 * t43 + t390 * t94) * t252;
t8 = t93 * t360 + t180 * t46 + t181 * t48 + t86 * t95 + t87 * t97 + (t253 * t44 + t390 * t93) * t252;
t7 = t94 * t359 + t182 * t45 + t183 * t47 + t84 * t96 + t85 * t98 + (t256 * t43 - t392 * t94) * t252;
t6 = t93 * t359 + t182 * t46 + t183 * t48 + t84 * t95 + t85 * t97 + (t256 * t44 - t392 * t93) * t252;
t5 = (-t213 + t49 + t399) * t256 + (t217 + t50) * t253 + (t250 * (-t409 + t483) + (-t418 - t481) * t249) * qJD(3) + (t436 * t256 + (t401 - t404 + t477) * t253) * qJD(1) + t296;
t4 = qJD(1) * t334 + t253 * t9 - t256 * t8;
t3 = qJD(1) * t333 + t253 * t7 - t256 * t6;
t2 = (qJD(3) * t334 - t19) * t255 + (-qJD(1) * t16 + qJD(3) * t52 + t253 * t8 + t256 * t9) * t252;
t1 = (qJD(3) * t333 - t18) * t255 + (-qJD(1) * t17 + qJD(3) * t53 + t253 * t6 + t256 * t7) * t252;
t34 = [t270 + (t25 * t57 + t26 * t56) * t464 + (t40 * t73 + t41 * t72) * t465 + (t127 * t75 + t128 * t74) * t466 + (t136 * t89 + t137 * t88) * t467 + (t326 - t319) * t384 + (t320 + t325) * t383 + (-t154 * t261 + t155 * t262) * t252 + (-Icges(6,3) * t255 + t252 * t312 - t317 + t324) * t388 + (-t118 * t252 - t162 * t381 - t362) * t251 + (-Icges(6,3) * t252 - t177 * t261 + t178 * t262 - t255 * t312 + t318 + t323) * t386; 0; 0; m(7) * (t25 * t76 + t26 * t77 + t32 * t56 + t33 * t57) + m(6) * (t113 * t40 + t114 * t41 + t54 * t72 + t55 * t73) + m(5) * (t101 * t127 + t102 * t128 + t170 * t74 + t171 * t75) + m(4) * ((-t253 * t88 - t256 * t89) * t237 + (-t136 * t256 - t137 * t253) * t221) + ((t161 * t488 + t281 * t459 + Icges(6,5) * t144 / 0.2e1 + Icges(6,6) * t143 / 0.2e1 + Icges(6,3) * t276 / 0.2e1) * t256 + (-Icges(6,5) * t142 / 0.2e1 - Icges(6,6) * t141 / 0.2e1 - Icges(6,3) * t275 / 0.2e1 + t160 * t488 + t281 * t489) * t253) * t255 + ((t421 / 0.2e1 + t419 / 0.2e1 + t199 * t461 + t200 * t460 - t137 * t453 + (t161 / 0.2e1 - t108 / 0.2e1) * t255 + (t164 / 0.2e1 + t110 * t455 + t112 * t454) * t252 + t373) * t256 + (t422 / 0.2e1 + t420 / 0.2e1 - t292 * t461 + t198 * t460 + t136 * t453 + (t160 / 0.2e1 - t107 / 0.2e1) * t255 + (t163 / 0.2e1 + t109 * t455 + t111 * t454) * t252 + t374) * t253) * qJD(1) + ((-qJD(1) * t186 - t256 * t282) * t267 + (-qJD(1) * t188 - t256 * t284) * t265 + t141 * t177 + t142 * t178 + t154 * t199 + t155 * t200 + t11 + t18 + (-qJD(1) * t163 - t256 * t283 - t261 * t68 + t262 * t70) * t252) * t459 + ((qJD(1) * t187 - t253 * t282) * t267 + (qJD(1) * t189 - t253 * t284) * t265 + t143 * t177 + t144 * t178 - t154 * t292 + t155 * t198 + t19 + t10 + (qJD(1) * t164 - t253 * t283 - t261 * t69 + t262 * t71) * t252) * t489 + (t286 / 0.2e1 - t285 / 0.2e1 + t288 / 0.2e1 - t287 / 0.2e1 + t496 * (t249 / 0.2e1 + t250 / 0.2e1) + (t108 * t252 - t110 * t411 + t112 * t410) * t459 + (t107 * t252 - t109 * t411 + t111 * t410) * t489) * qJD(3); m(4) * t42 + m(5) * t24 + m(6) * t12 + m(7) * t5; t17 * t390 + (t20 * t5 + t32 * t77 + t33 * t76) * t464 + t16 * t392 + (t113 * t55 + t114 * t54 + t31 * t12) * t465 + (t171 * t101 + t170 * t102 + t24 * t403) * t466 + t395 * t237 * t221 * t467 + (t168 * t485 + t190 * t484 + t3 + (t141 * t110 + t142 * t112 + t199 * t68 + t200 * t70 + t253 * t479) * t253 + t471 * t390 + (t253 * t474 - t370 + t486) * t392) * t253 + (-t4 + t169 * t485 + t191 * t484 + (t143 * t109 + t144 * t111 + t198 * t71 + t256 * t478 - t292 * t69) * t256 + t470 * t392 + (-t141 * t109 - t142 * t111 - t199 * t69 - t200 * t71 - t143 * t110 - t144 * t112 + t292 * t68 - t198 * t70 + (t160 * t386 + t163 * t388 + t186 * t383 + t188 * t384 + t479) * t256 + (t161 * t386 + t164 * t388 + t187 * t383 + t189 * t384 + t478 - t487) * t253 + ((-t160 * t255 - t163 * t252 - t420 - t422) * t256 + (-t161 * t255 - t164 * t252 - t419 - t421) * t253) * qJD(3) + (t372 + (-t490 + t492) * t253 + t470 + t471 + t495) * qJD(1)) * t253 + (-t256 * t490 + t371 - t486) * t390) * t256; m(7) * (qJD(1) * t330 - t25 * t256 + t253 * t26) + m(6) * (qJD(1) * t328 + t253 * t41 - t256 * t40) + m(5) * (t253 * t75 - t256 * t74 + (t127 * t256 + t128 * t253) * qJD(1)); 0; m(7) * (t253 * t32 - t256 * t33 + (t253 * t76 + t256 * t77) * qJD(1)) + m(6) * (t253 * t54 - t256 * t55 + (t113 * t253 + t114 * t256) * qJD(1)) + m(5) * (t101 * t253 - t102 * t256 + (t170 * t253 + t171 * t256) * qJD(1)); 0; 0.2e1 * (t328 * t463 + t330 * t462) * t386 + 0.2e1 * ((t25 * t253 + t256 * t26 + t390 * t57 - t392 * t56) * t462 + (t253 * t40 + t256 * t41 + t390 * t73 - t392 * t72) * t463) * t252; (m(6) + m(7)) * t388; 0.2e1 * ((t385 * t77 + t387 * t76 - t5) * t462 + (t113 * t387 + t114 * t385 - t12) * t463) * t255 + 0.2e1 * ((qJD(3) * t20 + t253 * t33 + t256 * t32 + t390 * t76 - t392 * t77) * t462 + (qJD(3) * t31 + t113 * t390 - t114 * t392 + t253 * t55 + t256 * t54) * t463) * t252; 0; 0.4e1 * (t463 + t462) * (-0.1e1 + t395) * t252 * t386; m(7) * (t22 * t57 + t23 * t56 + t25 * t65 + t26 * t64) + t58 + (-t21 + (t253 * t374 + t256 * t373) * qJD(3)) * t255 + ((t11 / 0.2e1 + t18 / 0.2e1) * t256 + (t10 / 0.2e1 + t19 / 0.2e1) * t253 + (-t253 * t373 + t256 * t374) * qJD(1)) * t252; m(7) * t15; m(7) * (t15 * t20 + t22 * t76 + t23 * t77 + t32 * t64 + t33 * t65 + t51 * t5) + (t14 * t447 + (qJD(1) * t39 - t10) * t458 + t17 * t357 - t2 / 0.2e1) * t256 + (t1 / 0.2e1 + (qJD(1) * t38 + t11) * t458 + t16 * t357 + t13 * t447) * t253 + (qJD(3) * t332 / 0.2e1 + t4 * t459 + t3 * t456 + (t16 * t456 - t253 * t17 / 0.2e1) * qJD(1)) * t252; m(7) * (qJD(1) * t329 - t22 * t256 + t23 * t253); m(7) * ((qJD(3) * t329 - t15) * t255 + (qJD(3) * t51 + t22 * t253 + t23 * t256 + (-t253 * t64 + t256 * t65) * qJD(1)) * t252); (t51 * t15 + t22 * t65 + t23 * t64) * t464 + (t21 * t255 - t58 + (t253 * t13 + t256 * t14 - t255 * t331) * qJD(3)) * t255 + (t256 * t1 + t253 * t2 - t255 * (t10 * t253 + t11 * t256) + (t252 * t331 - t255 * t59) * qJD(3) + (t256 * t13 - t253 * t14 + t255 * t332) * qJD(1)) * t252;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t34(1) t34(2) t34(4) t34(7) t34(11) t34(16); t34(2) t34(3) t34(5) t34(8) t34(12) t34(17); t34(4) t34(5) t34(6) t34(9) t34(13) t34(18); t34(7) t34(8) t34(9) t34(10) t34(14) t34(19); t34(11) t34(12) t34(13) t34(14) t34(15) t34(20); t34(16) t34(17) t34(18) t34(19) t34(20) t34(21);];
Mq  = res;

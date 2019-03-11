% Calculate time derivative of joint inertia matrix for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR8_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:47
% EndTime: 2019-03-09 02:59:12
% DurationCPUTime: 18.64s
% Computational Cost: add. (6930->765), mult. (18474->1046), div. (0->0), fcn. (16585->6), ass. (0->369)
t264 = cos(qJ(3));
t261 = sin(qJ(3));
t436 = Icges(5,5) * t261;
t444 = Icges(4,4) * t261;
t496 = t436 - t444 + (Icges(4,1) + Icges(5,1)) * t264;
t435 = Icges(5,5) * t264;
t440 = Icges(6,4) * t264;
t495 = t435 - t440 + (Icges(6,1) + Icges(5,3)) * t261;
t494 = t496 * qJD(3);
t493 = t495 * qJD(3);
t262 = sin(qJ(1));
t492 = -t262 / 0.2e1;
t462 = t262 / 0.2e1;
t265 = cos(qJ(1));
t392 = qJD(3) * t265;
t491 = -t392 / 0.2e1;
t490 = t392 / 0.2e1;
t489 = -qJD(1) / 0.2e1;
t488 = qJD(1) / 0.2e1;
t393 = qJD(3) * t264;
t364 = t262 * t393;
t396 = qJD(1) * t265;
t270 = t261 * t396 + t364;
t365 = t261 * t392;
t397 = qJD(1) * t264;
t487 = t262 * t397 + t365;
t260 = sin(qJ(6));
t263 = cos(qJ(6));
t438 = Icges(7,4) * t263;
t308 = -Icges(7,2) * t260 + t438;
t141 = Icges(7,6) * t264 + t261 * t308;
t439 = Icges(7,4) * t260;
t315 = Icges(7,1) * t263 - t439;
t148 = Icges(7,5) * t264 + t261 * t315;
t486 = -t141 * t260 + t148 * t263;
t485 = -qJD(3) / 0.2e1;
t421 = t261 * t262;
t238 = pkin(8) * t421;
t417 = t263 * t265;
t418 = t262 * t264;
t170 = t260 * t418 - t417;
t171 = -t260 * t265 - t263 * t418;
t97 = t171 * rSges(7,1) + t170 * rSges(7,2) + rSges(7,3) * t421;
t484 = -t238 - t97;
t394 = qJD(3) * t262;
t366 = t261 * t394;
t353 = qJD(6) + t397;
t354 = qJD(6) * t264 + qJD(1);
t419 = t262 * t263;
t86 = t354 * t419 + (t265 * t353 - t366) * t260;
t395 = qJD(3) * t261;
t87 = -t353 * t417 + (t260 * t354 + t263 * t395) * t262;
t45 = t87 * rSges(7,1) + t86 * rSges(7,2) + rSges(7,3) * t270;
t483 = -pkin(5) * t366 - pkin(8) * t270 - t45;
t346 = rSges(4,1) * t261 + rSges(4,2) * t264;
t286 = t265 * t346;
t443 = Icges(4,4) * t264;
t320 = Icges(4,1) * t261 + t443;
t474 = -Icges(4,5) * t262 + t265 * t320;
t313 = Icges(4,2) * t264 + t444;
t476 = -Icges(4,6) * t262 + t265 * t313;
t289 = -t261 * t474 - t264 * t476;
t278 = t289 * t262;
t146 = Icges(4,6) * t265 + t262 * t313;
t153 = Icges(4,5) * t265 + t262 * t320;
t290 = t146 * t264 + t153 * t261;
t279 = t290 * t265;
t310 = Icges(6,2) * t261 + t440;
t142 = -Icges(6,6) * t265 - t262 * t310;
t441 = Icges(6,4) * t261;
t317 = Icges(6,1) * t264 + t441;
t149 = -Icges(6,5) * t265 - t262 * t317;
t292 = t142 * t261 + t149 * t264;
t281 = t292 * t265;
t304 = -Icges(5,3) * t264 + t436;
t137 = Icges(5,6) * t265 + t262 * t304;
t318 = Icges(5,1) * t261 - t435;
t151 = Icges(5,4) * t265 + t262 * t318;
t294 = t137 * t264 - t151 * t261;
t283 = t294 * t265;
t406 = rSges(5,1) * t421 + t265 * rSges(5,2);
t482 = -t141 * t263 - t148 * t260;
t205 = rSges(6,1) * t261 - rSges(6,2) * t264;
t481 = qJD(3) * t205;
t480 = t262 * t353 + t365;
t303 = Icges(6,5) * t264 + Icges(6,6) * t261;
t136 = -Icges(6,3) * t262 + t265 * t303;
t479 = -Icges(5,6) * t262 + t265 * t304;
t306 = Icges(4,5) * t261 + Icges(4,6) * t264;
t478 = -Icges(4,3) * t262 + t265 * t306;
t143 = -Icges(6,6) * t262 + t265 * t310;
t311 = Icges(5,4) * t261 - Icges(5,6) * t264;
t477 = -Icges(5,2) * t262 + t265 * t311;
t150 = -Icges(6,5) * t262 + t265 * t317;
t475 = -Icges(5,4) * t262 + t265 * t318;
t342 = rSges(7,1) * t263 - rSges(7,2) * t260;
t390 = qJD(6) * t261;
t120 = (-rSges(7,1) * t260 - rSges(7,2) * t263) * t390 + (-rSges(7,3) * t261 + t264 * t342) * qJD(3);
t155 = rSges(7,3) * t264 + t261 * t342;
t398 = qJD(1) * t262;
t363 = t264 * t392;
t370 = t261 * t398;
t269 = -t363 + t370;
t288 = t265 * t354;
t84 = t260 * t480 - t263 * t288;
t85 = -t260 * t288 - t263 * t480;
t349 = t85 * rSges(7,1) + t84 * rSges(7,2);
t44 = rSges(7,3) * t269 + t349;
t416 = t264 * t265;
t172 = -t260 * t416 - t419;
t173 = -t262 * t260 + t263 * t416;
t343 = -t173 * rSges(7,1) - t172 * rSges(7,2);
t420 = t261 * t265;
t98 = -rSges(7,3) * t420 - t343;
t22 = (-t155 * t392 - t44) * t264 + (qJD(3) * t98 - t120 * t265 + t155 * t398) * t261;
t23 = (-t155 * t394 + t45) * t264 + (-qJD(3) * t97 - t262 * t120 - t155 * t396) * t261;
t69 = -t155 * t421 + t264 * t97;
t70 = -t155 * t420 - t264 * t98;
t473 = qJD(1) * (t262 * t69 + t265 * t70) + t22 * t262 - t23 * t265;
t472 = 2 * m(4);
t471 = 2 * m(5);
t470 = 2 * m(6);
t469 = 2 * m(7);
t258 = t262 ^ 2;
t259 = t265 ^ 2;
t468 = m(5) / 0.2e1;
t467 = m(6) / 0.2e1;
t466 = m(7) / 0.2e1;
t90 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t421;
t92 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t421;
t94 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t421;
t27 = t170 * t92 + t171 * t94 + t421 * t90;
t91 = Icges(7,5) * t173 + Icges(7,6) * t172 - Icges(7,3) * t420;
t93 = Icges(7,4) * t173 + Icges(7,2) * t172 - Icges(7,6) * t420;
t95 = Icges(7,1) * t173 + Icges(7,4) * t172 - Icges(7,5) * t420;
t28 = t170 * t93 + t171 * t95 + t421 * t91;
t18 = t28 * t262 + t265 * t27;
t465 = t18 / 0.2e1;
t464 = -pkin(1) - pkin(7);
t463 = -pkin(3) - pkin(4);
t461 = t264 / 0.2e1;
t460 = rSges(3,2) - pkin(1);
t459 = rSges(7,3) + pkin(8);
t453 = rSges(4,2) * t261;
t208 = rSges(4,1) * t264 - t453;
t458 = m(4) * t208;
t457 = pkin(5) * t261;
t456 = pkin(5) * t264;
t455 = rSges(5,1) * t261;
t454 = rSges(6,1) * t264;
t452 = rSges(6,2) * t261;
t451 = rSges(5,3) * t264;
t450 = rSges(6,3) * t265;
t449 = t262 * rSges(5,2);
t448 = t262 * rSges(4,3);
t254 = t265 * rSges(4,3);
t372 = pkin(4) * t363 + qJ(5) * t396 + qJD(5) * t262;
t376 = pkin(3) * t363 + qJ(4) * t487;
t391 = qJD(4) * t264;
t88 = t265 * (pkin(3) * t370 + t265 * t391 - t376);
t447 = t265 * (pkin(4) * t370 - t372) + t88;
t129 = pkin(4) * t270 + qJ(5) * t398 - qJD(5) * t265;
t377 = pkin(3) * t270 + qJ(4) * t366;
t96 = (-qJ(4) * t396 - qJD(4) * t262) * t264 + t377;
t446 = -t129 - t96;
t348 = -pkin(8) * t261 + t456;
t445 = t348 * t265 + t98;
t427 = qJ(5) * t265;
t106 = (-Icges(7,2) * t263 - t439) * t390 + (-Icges(7,6) * t261 + t264 * t308) * qJD(3);
t422 = t260 * t106;
t415 = -t348 * qJD(3) - t120;
t165 = qJD(4) * t261 + (-pkin(3) * t261 + qJ(4) * t264) * qJD(3);
t206 = pkin(3) * t264 + qJ(4) * t261;
t414 = t262 * t165 + t206 * t396;
t209 = pkin(8) * t264 + t457;
t413 = t155 + t209;
t240 = pkin(3) * t421;
t174 = -qJ(4) * t418 + t240;
t412 = rSges(5,3) * t418 - t174 - t406;
t243 = pkin(3) * t420;
t175 = qJ(4) * t416 - t243;
t164 = t265 * t175;
t194 = -pkin(4) * t420 - t262 * qJ(5);
t411 = t265 * t194 + t164;
t239 = pkin(4) * t421;
t410 = -t174 - t239 + t427;
t409 = -t175 - t194;
t180 = t262 * t206;
t408 = pkin(4) * t418 + t180;
t407 = rSges(6,1) * t366 + rSges(6,3) * t398;
t253 = t265 * qJ(2);
t405 = t243 + t253;
t404 = qJ(2) * t396 + qJD(2) * t262;
t403 = t265 * pkin(1) + t262 * qJ(2);
t402 = t258 + t259;
t135 = -Icges(6,3) * t265 - t262 * t303;
t401 = qJD(1) * t135;
t139 = Icges(4,3) * t265 + t262 * t306;
t400 = qJD(1) * t139;
t144 = Icges(5,2) * t265 + t262 * t311;
t399 = qJD(1) * t144;
t389 = -rSges(5,2) + t464;
t388 = -rSges(4,3) + t464;
t387 = rSges(6,3) + t464;
t386 = pkin(4) * t395;
t339 = -t260 * t92 + t263 * t94;
t39 = Icges(7,5) * t87 + Icges(7,6) * t86 + Icges(7,3) * t270;
t41 = Icges(7,4) * t87 + Icges(7,2) * t86 + Icges(7,6) * t270;
t43 = Icges(7,1) * t87 + Icges(7,4) * t86 + Icges(7,5) * t270;
t10 = (qJD(3) * t339 + t39) * t264 + (-qJD(3) * t90 - t260 * t41 + t263 * t43 + (-t260 * t94 - t263 * t92) * qJD(6)) * t261;
t113 = (-Icges(7,1) * t260 - t438) * t390 + (-Icges(7,5) * t261 + t264 * t315) * qJD(3);
t301 = Icges(7,5) * t263 - Icges(7,6) * t260;
t134 = Icges(7,3) * t264 + t261 * t301;
t99 = (-Icges(7,5) * t260 - Icges(7,6) * t263) * t390 + (-Icges(7,3) * t261 + t264 * t301) * qJD(3);
t17 = t170 * t106 + t171 * t113 + t134 * t270 + t86 * t141 + t87 * t148 + t421 * t99;
t385 = t10 / 0.2e1 + t17 / 0.2e1;
t338 = -t260 * t93 + t263 * t95;
t38 = Icges(7,5) * t85 + Icges(7,6) * t84 + Icges(7,3) * t269;
t40 = Icges(7,4) * t85 + Icges(7,2) * t84 + Icges(7,6) * t269;
t42 = Icges(7,1) * t85 + Icges(7,4) * t84 + Icges(7,5) * t269;
t11 = (qJD(3) * t338 + t38) * t264 + (-qJD(3) * t91 - t260 * t40 + t263 * t42 + (-t260 * t95 - t263 * t93) * qJD(6)) * t261;
t16 = t172 * t106 + t173 * t113 + t134 * t269 + t84 * t141 + t85 * t148 - t420 * t99;
t384 = -t11 / 0.2e1 - t16 / 0.2e1;
t33 = t261 * t339 + t264 * t90;
t48 = t134 * t421 + t141 * t170 + t148 * t171;
t383 = t33 / 0.2e1 + t48 / 0.2e1;
t34 = t261 * t338 + t264 * t91;
t49 = -t134 * t420 + t172 * t141 + t173 * t148;
t382 = t34 / 0.2e1 + t49 / 0.2e1;
t381 = t464 * t262;
t367 = t264 * t396;
t380 = pkin(4) * t367 + t414;
t344 = t452 + t454;
t379 = t262 * t344 + t410 + t450;
t178 = t206 * t398;
t378 = pkin(4) * t487 + t178;
t375 = -rSges(5,1) * t270 - rSges(5,3) * t366;
t374 = rSges(4,1) * t270 + rSges(4,2) * t367;
t157 = rSges(4,1) * t421 + rSges(4,2) * t418 + t254;
t371 = t265 * pkin(7) + t403;
t362 = t262 * t391;
t361 = (-pkin(5) - qJ(4)) * t264;
t360 = -pkin(4) * t264 - t206;
t359 = (-rSges(5,3) - qJ(4)) * t264;
t191 = t346 * qJD(3);
t358 = t191 * t402;
t357 = qJD(1) * t413;
t356 = t389 * t262;
t355 = pkin(5) * t418 + t410 + t484;
t250 = qJD(2) * t265;
t352 = t250 + t376;
t351 = t240 + t371;
t350 = t405 - t194;
t347 = (t303 + t306 + t311) * t485;
t207 = rSges(5,1) * t264 + rSges(5,3) * t261;
t345 = t451 - t455;
t323 = t377 + t404;
t266 = t129 + t323;
t24 = (t265 * t361 + t381) * qJD(1) + t266 - t362 - t483;
t271 = t352 + t372;
t25 = (-t391 + (t264 * t459 + t457) * qJD(3)) * t265 + (t464 * t265 + (t456 - qJ(2) + (-t459 + t463) * t261) * t262) * qJD(1) + t271 - t349;
t340 = -t24 * t265 + t262 * t25;
t337 = t262 * t27 - t265 * t28;
t29 = t172 * t92 + t173 * t94 - t420 * t90;
t30 = t172 * t93 + t173 * t95 - t420 * t91;
t19 = t30 * t262 + t265 * t29;
t336 = t262 * t29 - t265 * t30;
t31 = t262 * t357 + (-t165 + t415) * t265 + t378;
t32 = t265 * t357 + (-t386 - t415) * t262 + t380;
t335 = t32 * t262 - t265 * t31;
t334 = t34 * t262 + t265 * t33;
t333 = t262 * t33 - t265 * t34;
t285 = -t452 + (-rSges(6,1) - qJ(4)) * t264;
t267 = t285 * t265;
t35 = (-rSges(6,2) * qJD(3) - qJD(4)) * t418 + (t381 + t267) * qJD(1) + t266 + t407;
t36 = (-t391 + t481) * t265 + (t387 * t265 + (t454 - qJ(2) + (rSges(6,2) + t463) * t261) * t262) * qJD(1) + t271;
t332 = t262 * t36 - t265 * t35;
t46 = -t362 + (t265 * t359 + t356) * qJD(1) + t323 - t375;
t47 = (qJD(3) * t207 - t391) * t265 + (t389 * t265 + (t451 - qJ(2) + (-rSges(5,1) - pkin(3)) * t261) * t262) * qJD(1) + t352;
t331 = t262 * t47 - t265 * t46;
t51 = t381 + (t261 * t459 + t361) * t265 + t343 + t350;
t322 = t239 + t351;
t52 = t262 * t361 + t322 - t427 - t484;
t330 = t262 * t51 - t265 * t52;
t192 = t344 * qJD(3);
t54 = t205 * t398 + (-t165 - t192) * t265 + t378;
t55 = t205 * t396 + (t192 - t386) * t262 + t380;
t329 = t55 * t262 - t265 * t54;
t328 = t262 * t70 - t265 * t69;
t190 = t345 * qJD(3);
t71 = t207 * t398 + t178 + (-t165 - t190) * t265;
t72 = t262 * t190 + t207 * t396 + t414;
t326 = t72 * t262 - t265 * t71;
t75 = t262 * t387 + t267 + t350;
t76 = (-rSges(6,3) - qJ(5)) * t265 + t285 * t262 + t322;
t325 = t262 * t75 - t265 * t76;
t324 = t262 * t98 + t265 * t97;
t314 = -Icges(4,2) * t261 + t443;
t312 = Icges(5,4) * t264 + Icges(5,6) * t261;
t309 = -Icges(6,2) * t264 + t441;
t307 = Icges(4,5) * t264 - Icges(4,6) * t261;
t302 = Icges(6,5) * t261 - Icges(6,6) * t264;
t293 = t261 * t475 - t264 * t479;
t291 = t143 * t261 + t150 * t264;
t287 = rSges(3,3) * t265 + t262 * t460;
t284 = (t468 + t467 + t466) * t395;
t282 = t293 * t262;
t280 = t291 * t262;
t274 = qJD(3) * t314;
t273 = qJD(3) * t309;
t268 = t261 * t263 * t113 + t264 * t99 + t393 * t486;
t163 = -rSges(3,2) * t265 + t262 * rSges(3,3) + t403;
t162 = t253 + t287;
t161 = -t262 * rSges(6,3) + t265 * t344;
t160 = t448 - t286;
t159 = t265 * t345 + t449;
t132 = (-t206 - t207) * t265;
t131 = t207 * t262 + t180;
t127 = t250 + (t460 * t265 + (-rSges(3,3) - qJ(2)) * t262) * qJD(1);
t126 = qJD(1) * t287 + t404;
t124 = t371 + t157;
t123 = t262 * t388 + t253 + t286;
t122 = (-t205 + t360) * t265;
t121 = t205 * t262 + t408;
t110 = qJD(1) * t477 + t312 * t394;
t109 = -t312 * t392 + t399;
t105 = qJD(1) * t478 + t307 * t394;
t104 = -t307 * t392 + t400;
t101 = -qJD(1) * t136 + t302 * t394;
t100 = -t302 * t392 + t401;
t83 = t262 * t359 + t351 + t406;
t82 = (t359 + t455) * t265 + t356 + t405;
t78 = (t360 - t413) * t265;
t77 = t262 * t413 + t408;
t74 = t250 + t208 * t392 + (t388 * t265 + (-qJ(2) - t346) * t262) * qJD(1);
t73 = (-rSges(4,2) * t395 + qJD(1) * t388) * t262 + t374 + t404;
t68 = -t262 * t136 + t265 * t291;
t67 = -t262 * t135 + t281;
t66 = -t262 * t478 - t265 * t289;
t65 = t262 * t139 - t279;
t64 = -t262 * t477 + t265 * t293;
t63 = t262 * t144 + t283;
t62 = -t136 * t265 - t280;
t61 = -t135 * t265 - t262 * t292;
t60 = -t265 * t478 + t278;
t59 = t139 * t265 + t262 * t290;
t58 = -t265 * t477 - t282;
t57 = t144 * t265 - t262 * t294;
t56 = t134 * t264 + t261 * t486;
t53 = t159 * t265 + t262 * t412 + t164;
t50 = t324 * t261;
t37 = t161 * t265 + t262 * t379 + t411;
t26 = t262 * t355 + t265 * t445 + t411;
t21 = t88 + (-t96 + (-t159 - t175 + t449) * qJD(1) + t375) * t262 + (-t207 * t392 + (t412 + t406) * qJD(1)) * t265;
t20 = ((-qJD(3) * t134 + qJD(6) * t482 - t422) * t261 + t268) * t264;
t15 = -t259 * t481 + (rSges(6,2) * t364 - t407 + t446) * t262 + ((-t161 + t409) * t262 + (t379 - t450) * t265) * qJD(1) + t447;
t14 = t324 * t393 + (t262 * t44 + t265 * t45 + (-t262 * t97 + t265 * t98) * qJD(1)) * t261;
t13 = t261 * t336 + t49 * t264;
t12 = t261 * t337 + t48 * t264;
t9 = t91 * t364 + t170 * t40 + t171 * t42 + t86 * t93 + t87 * t95 + (t262 * t38 + t396 * t91) * t261;
t8 = t90 * t364 + t170 * t41 + t171 * t43 + t86 * t92 + t87 * t94 + (t262 * t39 + t396 * t90) * t261;
t7 = -t91 * t363 + t172 * t40 + t173 * t42 + t84 * t93 + t85 * t95 + (-t265 * t38 + t398 * t91) * t261;
t6 = -t90 * t363 + t172 * t41 + t173 * t43 + t84 * t92 + t85 * t94 + (-t265 * t39 + t398 * t90) * t261;
t5 = ((t409 - t445) * qJD(1) + t446 + t483) * t262 + (t44 - t209 * t392 + (t355 + t238) * qJD(1)) * t265 + t447;
t4 = -qJD(1) * t337 + t9 * t262 + t265 * t8;
t3 = -qJD(1) * t336 + t7 * t262 + t265 * t6;
t2 = (qJD(3) * t337 + t17) * t264 + (qJD(1) * t18 - qJD(3) * t48 + t262 * t8 - t265 * t9) * t261;
t1 = (qJD(3) * t336 + t16) * t264 + (qJD(1) * t19 - qJD(3) * t49 + t262 * t6 - t265 * t7) * t261;
t79 = [t268 + 0.2e1 * m(3) * (t126 * t163 + t127 * t162) + (t123 * t74 + t124 * t73) * t472 + (t46 * t83 + t47 * t82) * t471 + (t35 * t76 + t36 * t75) * t470 + (t24 * t52 + t25 * t51) * t469 + t482 * t390 + (-t320 - t318 - t310) * t393 + (-t274 + t493) * t264 + (t317 + t313 - t304 - t134) * t395 + (t273 - t422 - t494) * t261; m(7) * ((t262 * t52 + t265 * t51) * qJD(1) + t340) + m(5) * ((t262 * t83 + t265 * t82) * qJD(1) + t331) + m(6) * ((t262 * t76 + t265 * t75) * qJD(1) + t332) + m(4) * (t262 * t74 - t265 * t73 + (t123 * t265 + t124 * t262) * qJD(1)) + m(3) * (-t126 * t265 + t262 * t127 + (t162 * t265 + t163 * t262) * qJD(1)); 0; m(5) * (t131 * t47 + t132 * t46 + t71 * t83 + t72 * t82) + m(6) * (t121 * t36 + t122 * t35 + t54 * t76 + t55 * t75) + m(7) * (t24 * t78 + t25 * t77 + t31 * t52 + t32 * t51) + (m(4) * (t124 * t191 - t208 * t73) + t347 * t265 + t385 + (t273 * t492 + (t143 + t474 + t475) * t488 + t494 * t462) * t264) * t265 + (m(4) * (-t123 * t191 + t208 * t74) + t347 * t262 + (t142 * t489 + t309 * t490 + t496 * t491 + (t151 + t153) * t488) * t264 - t384) * t262 + ((t274 * t492 + t479 * t488 + (t150 + t476) * t489 + t493 * t462) * t265 + (t146 * t489 + t314 * t490 + t495 * t491 + (t137 + t149) * t488) * t262) * t261 + (t281 / 0.2e1 + t283 / 0.2e1 - t279 / 0.2e1 - t278 / 0.2e1 + (t291 + t293) * t462) * qJD(3) + ((t124 * t458 + (t142 / 0.2e1 - t151 / 0.2e1 - t153 / 0.2e1) * t264 + (-t149 / 0.2e1 - t137 / 0.2e1 + t146 / 0.2e1) * t261 - t383) * t262 + (t123 * t458 + (-t474 / 0.2e1 - t143 / 0.2e1 - t475 / 0.2e1) * t264 + (t476 / 0.2e1 + t150 / 0.2e1 - t479 / 0.2e1) * t261 + t382) * t265) * qJD(1); m(5) * ((t131 * t265 + t132 * t262) * qJD(1) + t326) + m(6) * ((t121 * t265 + t122 * t262) * qJD(1) + t329) + m(7) * ((t262 * t78 + t265 * t77) * qJD(1) + t335) - m(4) * t358; (t26 * t5 + t31 * t78 + t32 * t77) * t469 + t265 * t4 + t262 * t3 + (t131 * t72 + t132 * t71 + t21 * t53) * t471 + (t121 * t55 + t122 * t54 + t15 * t37) * t470 + ((-t262 * t157 + t160 * t265) * (-t262 * t374 + (-t208 * t259 + t258 * t453) * qJD(3) + ((-t157 + t254) * t265 + (-t160 + t286 + t448) * t262) * qJD(1)) - t208 * t358) * t472 + t265 * ((-t265 * t101 + (t62 - t281) * qJD(1)) * t265 + (-t61 * qJD(1) + (-t143 * t393 + t150 * t395) * t262 + (-t100 + (-t142 * t264 + t149 * t261) * qJD(3) + (t135 - t291) * qJD(1)) * t265) * t262) + t262 * ((-t262 * t100 + (-t67 - t280) * qJD(1)) * t262 + (t68 * qJD(1) + (t142 * t393 - t149 * t395 - t401) * t265 + (-t101 + (t143 * t264 - t150 * t261) * qJD(3) - t292 * qJD(1)) * t262) * t265) + t265 * ((t265 * t105 + (t60 + t279) * qJD(1)) * t265 + (-t59 * qJD(1) + (-t393 * t474 + t395 * t476) * t262 + (t104 + (-t146 * t261 + t153 * t264) * qJD(3) + (-t139 + t289) * qJD(1)) * t265) * t262) + t262 * ((t262 * t104 + (-t65 + t278) * qJD(1)) * t262 + (t66 * qJD(1) + (t146 * t395 - t153 * t393 + t400) * t265 + (t105 + (-t261 * t476 + t264 * t474) * qJD(3) + t290 * qJD(1)) * t262) * t265) + t262 * ((t262 * t109 + (-t63 - t282) * qJD(1)) * t262 + (t64 * qJD(1) + (-t137 * t395 - t151 * t393 + t399) * t265 + (t110 + (t261 * t479 + t264 * t475) * qJD(3) - t294 * qJD(1)) * t262) * t265) + t265 * ((t265 * t110 + (t58 - t283) * qJD(1)) * t265 + (-t57 * qJD(1) + (-t393 * t475 - t395 * t479) * t262 + (t109 + (t137 * t261 + t151 * t264) * qJD(3) + (-t144 - t293) * qJD(1)) * t265) * t262) + (-t18 + (-t57 - t59 - t61) * t265 + (-t58 - t60 - t62) * t262) * t398 + (t19 + (t63 + t65 + t67) * t265 + (t64 + t66 + t68) * t262) * t396; 0.2e1 * (t330 * t466 + (t262 * t82 - t265 * t83) * t468 + t325 * t467) * t395 + 0.2e1 * ((-t396 * t51 - t398 * t52 - t340) * t466 + (-t396 * t82 - t398 * t83 - t331) * t468 + (-t396 * t75 - t398 * t76 - t332) * t467) * t264; 0.2e1 * t402 * t284; 0.2e1 * ((-t392 * t78 + t394 * t77 + t5) * t466 + (t131 * t394 - t132 * t392 + t21) * t468 + (t121 * t394 - t122 * t392 + t15) * t467) * t261 + 0.2e1 * ((qJD(3) * t26 - t396 * t77 - t398 * t78 - t335) * t466 + (qJD(3) * t53 - t131 * t396 - t132 * t398 - t326) * t468 + (qJD(3) * t37 - t121 * t396 - t122 * t398 - t329) * t467) * t264; 0.4e1 * (0.1e1 - t402) * t264 * t284; m(7) * (qJD(1) * t330 - t262 * t24 - t25 * t265) + m(6) * (qJD(1) * t325 - t262 * t35 - t265 * t36); 0; m(7) * (-t262 * t31 - t265 * t32 + (t262 * t77 - t265 * t78) * qJD(1)) + m(6) * (-t262 * t54 - t265 * t55 + (t121 * t262 - t122 * t265) * qJD(1)); 0; 0; t20 + m(7) * (t22 * t51 + t23 * t52 + t24 * t69 + t25 * t70) + (t262 * t383 - t265 * t382) * t393 + (-qJD(3) * t56 + t384 * t265 + t385 * t262 + (t262 * t382 + t265 * t383) * qJD(1)) * t261; m(7) * t473; m(7) * (t14 * t26 + t22 * t77 + t23 * t78 + t31 * t69 + t32 * t70 + t5 * t50) + (t2 / 0.2e1 + (qJD(1) * t34 + t10) * t461 + t13 * t488 - t19 * t393 / 0.2e1) * t265 + (t12 * t489 + (-qJD(1) * t33 + t11) * t461 + t1 / 0.2e1 + t393 * t465) * t262 + (t334 * t485 + t4 * t462 - t265 * t3 / 0.2e1 + (t19 * t462 + t265 * t465) * qJD(1)) * t261; m(7) * ((qJD(3) * t328 + t14) * t261 + (qJD(3) * t50 - t473) * t264); m(7) * (qJD(1) * t328 - t22 * t265 - t23 * t262); (t14 * t50 + t22 * t70 + t23 * t69) * t469 + (t20 + (t262 * t12 - t265 * t13 + t264 * t333) * qJD(3)) * t264 + (t262 * t2 - t265 * t1 + t264 * (t10 * t262 - t11 * t265) + (-t261 * t333 - 0.2e1 * t56 * t264) * qJD(3) + (t265 * t12 + t262 * t13 + t264 * t334) * qJD(1)) * t261;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t79(1) t79(2) t79(4) t79(7) t79(11) t79(16); t79(2) t79(3) t79(5) t79(8) t79(12) t79(17); t79(4) t79(5) t79(6) t79(9) t79(13) t79(18); t79(7) t79(8) t79(9) t79(10) t79(14) t79(19); t79(11) t79(12) t79(13) t79(14) t79(15) t79(20); t79(16) t79(17) t79(18) t79(19) t79(20) t79(21);];
Mq  = res;

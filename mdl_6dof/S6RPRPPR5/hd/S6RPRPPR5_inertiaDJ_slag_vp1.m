% Calculate time derivative of joint inertia matrix for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:11
% EndTime: 2019-03-09 02:50:33
% DurationCPUTime: 17.30s
% Computational Cost: add. (18324->797), mult. (20749->1135), div. (0->0), fcn. (18963->10), ass. (0->402)
t270 = pkin(9) + qJ(3);
t262 = cos(t270);
t260 = sin(t270);
t453 = Icges(5,6) * t260;
t461 = Icges(4,4) * t260;
t521 = t453 + t461 + (Icges(4,2) + Icges(5,3)) * t262;
t452 = Icges(5,6) * t262;
t460 = Icges(4,4) * t262;
t520 = t452 + t460 + (Icges(4,1) + Icges(5,2)) * t260;
t519 = -Icges(6,3) / 0.2e1;
t518 = t521 * qJD(3);
t517 = t520 * qJD(3);
t273 = sin(pkin(10));
t279 = sin(qJ(1));
t434 = t279 * t273;
t275 = cos(pkin(10));
t280 = cos(qJ(1));
t438 = t275 * t280;
t201 = t260 * t438 - t434;
t413 = qJD(3) * t279;
t384 = t262 * t413;
t149 = qJD(1) * t201 + t275 * t384;
t516 = -t149 / 0.2e1;
t439 = t273 * t280;
t394 = t260 * t439;
t433 = t279 * t275;
t202 = t394 + t433;
t150 = qJD(1) * t202 + t273 * t384;
t515 = -t150 / 0.2e1;
t514 = t279 / 0.2e1;
t513 = -t280 / 0.2e1;
t512 = -qJD(1) / 0.2e1;
t472 = qJD(1) / 0.2e1;
t412 = qJD(3) * t280;
t385 = t260 * t412;
t417 = qJD(1) * t279;
t511 = t262 * t417 + t385;
t269 = pkin(10) + qJ(6);
t261 = cos(t269);
t259 = sin(t269);
t458 = Icges(7,4) * t259;
t338 = Icges(7,2) * t261 + t458;
t157 = Icges(7,6) * t260 - t262 * t338;
t457 = Icges(7,4) * t261;
t343 = Icges(7,1) * t259 + t457;
t158 = Icges(7,5) * t260 - t262 * t343;
t510 = t157 * t261 + t158 * t259;
t420 = t279 ^ 2 + t280 ^ 2;
t509 = -0.1e1 + t420;
t508 = qJD(3) / 0.2e1;
t332 = -Icges(5,3) * t260 + t452;
t175 = Icges(5,5) * t279 - t280 * t332;
t334 = Icges(5,2) * t262 - t453;
t177 = Icges(5,4) * t279 - t280 * t334;
t321 = t175 * t260 - t177 * t262;
t507 = t279 * t321;
t342 = -Icges(4,2) * t260 + t460;
t172 = Icges(4,6) * t279 + t280 * t342;
t346 = Icges(4,1) * t262 - t461;
t174 = Icges(4,5) * t279 + t280 * t346;
t322 = t172 * t260 - t174 * t262;
t506 = t279 * t322;
t497 = Icges(5,4) * t280 + t279 * t334;
t498 = Icges(5,5) * t280 + t279 * t332;
t320 = -t260 * t498 + t262 * t497;
t505 = t280 * t320;
t171 = -Icges(4,6) * t280 + t279 * t342;
t173 = -Icges(4,5) * t280 + t279 * t346;
t323 = t171 * t260 - t173 * t262;
t504 = t280 * t323;
t436 = t279 * t259;
t443 = t261 * t280;
t191 = t260 * t443 - t436;
t435 = t279 * t261;
t444 = t260 * t280;
t192 = t259 * t444 + t435;
t440 = t262 * t280;
t105 = t192 * rSges(7,1) + t191 * rSges(7,2) + rSges(7,3) * t440;
t252 = pkin(5) * t275 + pkin(4);
t244 = t279 * t252;
t277 = -pkin(8) - qJ(5);
t393 = t277 * t440;
t503 = pkin(5) * t394 + t105 + t244 - t393;
t266 = t279 * rSges(5,1);
t502 = -rSges(5,2) * t440 + t266;
t267 = t279 * pkin(4);
t215 = qJ(5) * t440 + t267;
t265 = t279 * rSges(4,3);
t501 = -rSges(4,2) * t444 + t265;
t278 = -pkin(7) - qJ(2);
t276 = cos(pkin(9));
t253 = pkin(2) * t276 + pkin(1);
t469 = rSges(4,1) * t262;
t364 = -rSges(4,2) * t260 + t469;
t312 = -t253 - t364;
t142 = (rSges(4,3) - t278) * t280 + t312 * t279;
t186 = rSges(4,1) * t440 + t501;
t377 = t280 * t253 - t279 * t278;
t143 = t186 + t377;
t500 = t142 * t280 + t143 * t279;
t193 = t259 * t280 + t260 * t435;
t194 = t260 * t436 - t443;
t360 = -t194 * rSges(7,1) - t193 * rSges(7,2);
t441 = t262 * t279;
t106 = rSges(7,3) * t441 - t360;
t328 = t105 * t279 - t106 * t280;
t415 = qJD(3) * t260;
t386 = t260 * t413;
t416 = qJD(1) * t280;
t387 = t262 * t416;
t300 = -t386 + t387;
t375 = qJD(1) * t260 + qJD(6);
t288 = t280 * t375 + t384;
t376 = qJD(6) * t260 + qJD(1);
t319 = t259 * t376;
t89 = t261 * t288 - t279 * t319;
t318 = t261 * t376;
t90 = t259 * t288 + t279 * t318;
t369 = t90 * rSges(7,1) + t89 * rSges(7,2);
t48 = rSges(7,3) * t300 + t369;
t383 = t262 * t412;
t287 = -t279 * t375 + t383;
t91 = t261 * t287 - t280 * t319;
t92 = t259 * t287 + t280 * t318;
t470 = t92 * rSges(7,1) + t91 * rSges(7,2);
t49 = -rSges(7,3) * t511 + t470;
t14 = t328 * t415 + (-t279 * t49 + t280 * t48 + (-t105 * t280 - t279 * t106) * qJD(1)) * t262;
t359 = rSges(7,1) * t259 + rSges(7,2) * t261;
t160 = rSges(7,3) * t260 - t262 * t359;
t61 = -t106 * t260 + t160 * t441;
t62 = t260 * t105 - t160 * t440;
t350 = t279 * t62 + t280 * t61;
t499 = qJD(3) * t350 - t14;
t337 = Icges(4,5) * t262 - Icges(4,6) * t260;
t169 = -Icges(4,3) * t280 + t279 * t337;
t340 = Icges(5,4) * t262 - Icges(5,5) * t260;
t496 = Icges(5,1) * t280 + t279 * t340;
t316 = rSges(3,1) * t276 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t465 = rSges(3,3) + qJ(2);
t163 = t279 * t465 + t280 * t316;
t495 = 2 * m(4);
t494 = 2 * m(5);
t493 = 2 * m(6);
t492 = 2 * m(7);
t491 = m(5) / 0.2e1;
t490 = -m(6) / 0.2e1;
t489 = m(6) / 0.2e1;
t488 = -m(7) / 0.2e1;
t487 = m(7) / 0.2e1;
t339 = Icges(6,4) * t273 + Icges(6,2) * t275;
t165 = Icges(6,6) * t260 - t262 * t339;
t486 = t165 / 0.2e1;
t344 = Icges(6,1) * t273 + Icges(6,4) * t275;
t166 = Icges(6,5) * t260 - t262 * t344;
t485 = t166 / 0.2e1;
t484 = t201 / 0.2e1;
t483 = t202 / 0.2e1;
t482 = t260 / 0.2e1;
t481 = -t273 / 0.2e1;
t480 = t273 / 0.2e1;
t479 = -t275 / 0.2e1;
t478 = t275 / 0.2e1;
t477 = t280 / 0.2e1;
t476 = -rSges(7,3) - pkin(3);
t225 = rSges(4,1) * t260 + rSges(4,2) * t262;
t475 = m(4) * t225;
t474 = pkin(3) * t262;
t473 = pkin(5) * t273;
t335 = Icges(7,5) * t259 + Icges(7,6) * t261;
t156 = Icges(7,3) * t260 - t262 * t335;
t409 = qJD(6) * t262;
t414 = qJD(3) * t262;
t93 = (-Icges(7,5) * t261 + Icges(7,6) * t259) * t409 + (Icges(7,3) * t262 + t260 * t335) * qJD(3);
t365 = t259 * t157 * t409 + t156 * t414 + t260 * t93 + t415 * t510;
t95 = (-Icges(7,1) * t261 + t458) * t409 + (Icges(7,5) * t262 + t260 * t343) * qJD(3);
t466 = t259 * t95;
t60 = t156 * t260 - t262 * t510;
t94 = (Icges(7,2) * t259 - t457) * t409 + (Icges(7,6) * t262 + t260 * t338) * qJD(3);
t471 = ((-t466 + (-qJD(6) * t158 - t94) * t261) * t262 + t365) * t260 + t60 * t414;
t468 = rSges(5,2) * t260;
t467 = rSges(6,3) * t260;
t464 = -rSges(5,3) - qJ(4);
t463 = rSges(7,3) - t277;
t450 = qJ(4) * t260;
t449 = qJ(4) * t262;
t442 = t262 * t273;
t437 = t278 * t280;
t432 = -qJ(5) - t277;
t431 = -t215 + t503;
t268 = t280 * pkin(4);
t301 = t260 * t473 + t262 * t432;
t430 = -t252 * t280 + t279 * t301 + t106 + t268;
t203 = t260 * t433 + t439;
t151 = -qJD(1) * t203 + t275 * t383;
t204 = t260 * t434 - t438;
t370 = t273 * t383;
t152 = -qJD(1) * t204 + t370;
t429 = t152 * rSges(6,1) + t151 * rSges(6,2);
t405 = pkin(5) * t442;
t428 = t260 * t432 + t160 - t405;
t356 = t450 + t474;
t182 = qJD(3) * t356 - qJD(4) * t262;
t358 = -rSges(5,2) * t262 + rSges(5,3) * t260;
t427 = -t358 * qJD(3) - t182;
t197 = t356 * t279;
t198 = pkin(3) * t440 + qJ(4) * t444;
t426 = t279 * t197 + t280 * t198;
t425 = -t198 - t215;
t223 = pkin(3) * t260 - t449;
t200 = t223 * t417;
t389 = t260 * t417;
t424 = qJ(5) * t389 + t200;
t357 = rSges(5,3) * t262 + t468;
t423 = -t223 + t357;
t411 = qJD(4) * t260;
t422 = qJ(4) * t383 + t280 * t411;
t264 = qJD(2) * t280;
t421 = t278 * t417 + t264;
t170 = Icges(4,3) * t279 + t280 * t337;
t419 = qJD(1) * t170;
t179 = Icges(5,1) * t279 - t280 * t340;
t418 = qJD(1) * t179;
t410 = qJD(5) * t262;
t408 = t489 + t487;
t407 = -rSges(6,3) - pkin(3) - qJ(5);
t101 = Icges(7,1) * t192 + Icges(7,4) * t191 + Icges(7,5) * t440;
t99 = Icges(7,4) * t192 + Icges(7,2) * t191 + Icges(7,6) * t440;
t347 = t101 * t259 + t261 * t99;
t43 = Icges(7,5) * t92 + Icges(7,6) * t91 - Icges(7,3) * t511;
t45 = Icges(7,4) * t92 + Icges(7,2) * t91 - Icges(7,6) * t511;
t47 = Icges(7,1) * t92 + Icges(7,4) * t91 - Icges(7,5) * t511;
t97 = Icges(7,5) * t192 + Icges(7,6) * t191 + Icges(7,3) * t440;
t10 = (qJD(3) * t347 + t43) * t260 + (qJD(3) * t97 - t259 * t47 - t261 * t45 + (-t101 * t261 + t259 * t99) * qJD(6)) * t262;
t17 = -t156 * t511 + t91 * t157 + t92 * t158 + t191 * t94 + t192 * t95 + t440 * t93;
t402 = t10 / 0.2e1 + t17 / 0.2e1;
t100 = Icges(7,4) * t194 + Icges(7,2) * t193 + Icges(7,6) * t441;
t102 = Icges(7,1) * t194 + Icges(7,4) * t193 + Icges(7,5) * t441;
t330 = t100 * t261 + t102 * t259;
t42 = Icges(7,5) * t90 + Icges(7,6) * t89 + Icges(7,3) * t300;
t44 = Icges(7,4) * t90 + Icges(7,2) * t89 + Icges(7,6) * t300;
t46 = Icges(7,1) * t90 + Icges(7,4) * t89 + Icges(7,5) * t300;
t98 = Icges(7,5) * t194 + Icges(7,6) * t193 + Icges(7,3) * t441;
t11 = (qJD(3) * t330 + t42) * t260 + (qJD(3) * t98 - t259 * t46 - t261 * t44 + (t100 * t259 - t102 * t261) * qJD(6)) * t262;
t16 = t156 * t300 + t89 * t157 + t90 * t158 + t193 * t94 + t194 * t95 + t441 * t93;
t401 = t11 / 0.2e1 + t16 / 0.2e1;
t35 = t260 * t97 - t262 * t347;
t50 = t156 * t440 + t191 * t157 + t192 * t158;
t400 = -t35 / 0.2e1 - t50 / 0.2e1;
t36 = t260 * t98 - t262 * t330;
t51 = t156 * t441 + t157 * t193 + t158 * t194;
t399 = t36 / 0.2e1 + t51 / 0.2e1;
t112 = Icges(6,5) * t202 + Icges(6,6) * t201 + Icges(6,3) * t440;
t398 = t112 * t441;
t397 = t112 * t440;
t113 = Icges(6,5) * t204 + Icges(6,6) * t203 + Icges(6,3) * t441;
t396 = t113 * t441;
t395 = t113 * t440;
t239 = pkin(3) * t386;
t392 = t279 * (pkin(3) * t387 + t279 * t411 - t239 + (t260 * t416 + t384) * qJ(4)) + t280 * (-pkin(3) * t511 - qJ(4) * t389 + t422) + t197 * t416;
t130 = t202 * rSges(6,1) + t201 * rSges(6,2) + rSges(6,3) * t440;
t263 = qJD(2) * t279;
t391 = t263 + t422;
t390 = t239 + t421;
t381 = -t340 * qJD(3) / 0.2e1 + t337 * t508;
t380 = -t415 / 0.2e1;
t379 = -qJ(4) - t473;
t146 = t423 * t280;
t378 = -qJ(5) * t260 - t223;
t216 = qJ(5) * t441 - t268;
t374 = t280 * t215 + t279 * t216 + t426;
t373 = pkin(5) * t370 + t252 * t416 + t277 * t511;
t240 = t280 * t410;
t372 = t240 + t391;
t371 = rSges(5,1) * t416 + rSges(5,2) * t511 + rSges(5,3) * t383;
t361 = rSges(6,1) * t273 + rSges(6,2) * t275;
t167 = -t262 * t361 + t467;
t368 = -t167 + t378;
t367 = t112 / 0.2e1 + t174 / 0.2e1 - t177 / 0.2e1;
t366 = t173 / 0.2e1 + t497 / 0.2e1 + t113 / 0.2e1;
t363 = -t150 * rSges(6,1) - t149 * rSges(6,2);
t362 = -t204 * rSges(6,1) - t203 * rSges(6,2);
t29 = t192 * t101 + t191 * t99 + t440 * t97;
t30 = t191 * t100 + t192 * t102 + t440 * t98;
t355 = t279 * t30 + t280 * t29;
t18 = t29 * t279 - t280 * t30;
t31 = t101 * t194 + t193 * t99 + t441 * t97;
t32 = t100 * t193 + t102 * t194 + t441 * t98;
t354 = t279 * t32 + t280 * t31;
t19 = t31 * t279 - t280 * t32;
t353 = t279 * t36 + t280 * t35;
t352 = t35 * t279 - t280 * t36;
t309 = t260 * t379 - t253;
t285 = (-pkin(3) - t463) * t262 + t309;
t58 = (t252 - t278) * t280 + t285 * t279 + t360;
t313 = t377 + t198;
t59 = t313 + t503;
t351 = t279 * t59 + t280 * t58;
t295 = t262 * t407 - t253 - t450;
t284 = t279 * t295 - t437;
t64 = t268 + t284 + t362;
t65 = t313 + t130 + t215;
t349 = t279 * t65 + t280 * t64;
t348 = t260 * t464 - t253;
t336 = Icges(6,5) * t273 + Icges(6,6) * t275;
t292 = (rSges(5,2) - pkin(3)) * t262 + t348;
t103 = (rSges(5,1) - t278) * t280 + t292 * t279;
t187 = rSges(5,3) * t444 + t502;
t104 = t187 + t313;
t329 = t103 * t280 + t104 * t279;
t317 = -t410 - t411;
t314 = t378 - t428;
t231 = qJ(5) * t386;
t256 = pkin(4) * t416;
t311 = t279 * (qJD(1) * t215 + t279 * t410 - t231) + t280 * (-qJ(5) * t511 + t240 + t256) + t216 * t416 + t392;
t108 = t368 * t280;
t308 = qJD(3) * t225;
t305 = qJD(3) * (Icges(5,4) * t260 + Icges(5,5) * t262);
t304 = qJD(3) * (-Icges(4,5) * t260 - Icges(4,6) * t262);
t81 = t314 * t280;
t298 = -qJ(5) * t414 - qJD(5) * t260 - t182;
t5 = (qJ(5) * t385 - t256 + t373 + t49) * t280 + (t231 + t48 + (t260 * t277 + t405) * t413) * t279 + (t430 * t280 + (-t393 + (-pkin(4) + t252) * t279 + t425 - t431) * t279) * qJD(1) + t311;
t80 = t314 * t279;
t296 = t412 * t81 + t413 * t80 - t5;
t107 = t368 * t279;
t131 = rSges(6,3) * t441 - t362;
t15 = t279 * (-rSges(6,3) * t386 - t363) + t280 * (-rSges(6,3) * t385 + t429) + (t280 * t131 + (-t130 + t425) * t279) * qJD(1) + t311;
t294 = t107 * t413 + t108 * t412 - t15;
t293 = -(rSges(6,3) * t262 + t260 * t361) * qJD(3) + t298;
t96 = (-rSges(7,1) * t261 + rSges(7,2) * t259) * t409 + (rSges(7,3) * t262 + t260 * t359) * qJD(3);
t291 = -t301 * qJD(3) + t298 - t96;
t114 = Icges(6,4) * t202 + Icges(6,2) * t201 + Icges(6,6) * t440;
t116 = Icges(6,1) * t202 + Icges(6,4) * t201 + Icges(6,5) * t440;
t290 = t114 * t478 + t116 * t480 - t172 / 0.2e1 + t175 / 0.2e1;
t115 = Icges(6,4) * t204 + Icges(6,2) * t203 + Icges(6,6) * t441;
t117 = Icges(6,1) * t204 + Icges(6,4) * t203 + Icges(6,5) * t441;
t289 = t115 * t479 + t117 * t481 + t498 / 0.2e1 + t171 / 0.2e1;
t286 = rSges(4,2) * t389 + rSges(4,3) * t416 - t280 * t308;
t162 = -t279 * t316 + t280 * t465;
t21 = (t160 * t412 + t49) * t260 + (qJD(3) * t105 + t160 * t417 - t280 * t96) * t262;
t22 = (-t160 * t413 - t48) * t260 + (-qJD(3) * t106 + t160 * t416 + t279 * t96) * t262;
t54 = t328 * t262;
t283 = -qJD(3) * t54 + t21 * t279 + t22 * t280 + (-t279 * t61 + t280 * t62) * qJD(1);
t25 = t476 * t385 + (-t437 + (t262 * t476 + t309) * t279) * qJD(1) + t372 + t373 + t470;
t26 = ((t260 * t463 + t262 * t379) * qJD(3) + t317) * t279 + (t280 * t285 - t244) * qJD(1) - t369 + t390;
t33 = qJD(1) * t284 + t385 * t407 + t256 + t372 + t429;
t34 = t231 + ((-t449 + t467) * qJD(3) + t317) * t279 + (t280 * t295 - t267) * qJD(1) + t363 + t390;
t282 = (t25 * t279 + t26 * t280 + t416 * t59 - t417 * t58) * t487 + (t279 * t33 + t280 * t34 + t416 * t65 - t417 * t64) * t489;
t24 = t279 * t430 + t280 * t431 + t374;
t27 = t280 * t291 + t417 * t428 + t424;
t28 = qJD(1) * t81 + t279 * t291;
t41 = t130 * t280 + t279 * t131 + t374;
t52 = t167 * t417 + t280 * t293 + t424;
t53 = qJD(1) * t108 + t279 * t293;
t281 = (qJD(3) * t24 + t27 * t280 + t279 * t28 + t416 * t80 - t417 * t81) * t487 + (qJD(3) * t41 + t107 * t416 - t108 * t417 + t279 * t53 + t280 * t52) * t489;
t214 = t364 * qJD(3);
t188 = -rSges(5,1) * t280 + t279 * t358;
t185 = -rSges(4,3) * t280 + t279 * t364;
t155 = (Icges(6,5) * t262 + t260 * t344) * qJD(3);
t154 = (Icges(6,6) * t262 + t260 * t339) * qJD(3);
t145 = t423 * t279;
t141 = -qJD(1) * t163 + t264;
t140 = qJD(1) * t162 + t263;
t137 = t509 * t260 * t414;
t129 = qJD(1) * t496 + t280 * t305;
t128 = t279 * t305 + t418;
t119 = t279 * t304 + t419;
t118 = -qJD(1) * t169 + t280 * t304;
t85 = t225 * t413 + (t280 * t312 - t265) * qJD(1) + t421;
t84 = t263 + (-t437 + (-t253 - t469) * t279) * qJD(1) + t286;
t83 = qJD(1) * t146 + t279 * t427;
t82 = t280 * t427 - t357 * t417 + t200;
t79 = t279 * t320 + t280 * t496;
t78 = -t179 * t280 + t507;
t77 = -t279 * t496 + t505;
t76 = t279 * t179 + t280 * t321;
t75 = t279 * t170 - t280 * t322;
t74 = t279 * t169 - t504;
t73 = -t170 * t280 - t506;
t72 = -t169 * t280 - t279 * t323;
t71 = Icges(6,1) * t152 + Icges(6,4) * t151 - Icges(6,5) * t511;
t70 = Icges(6,1) * t150 + Icges(6,4) * t149 + Icges(6,5) * t300;
t69 = Icges(6,4) * t152 + Icges(6,2) * t151 - Icges(6,6) * t511;
t68 = Icges(6,4) * t150 + Icges(6,2) * t149 + Icges(6,6) * t300;
t63 = t187 * t280 + t279 * t188 + t426;
t57 = (-t411 + (t262 * t464 - t468) * qJD(3)) * t279 + (t280 * t292 - t266) * qJD(1) + t390;
t56 = -pkin(3) * t385 + (-t437 + (t348 - t474) * t279) * qJD(1) + t371 + t391;
t40 = t115 * t203 + t117 * t204 + t396;
t39 = t114 * t203 + t116 * t204 + t398;
t38 = t201 * t115 + t202 * t117 + t395;
t37 = t201 * t114 + t202 * t116 + t397;
t23 = (qJD(1) * t188 + t371) * t280 + (t357 * t413 + (-t187 - t198 + t502) * qJD(1)) * t279 + t392;
t13 = t51 * t260 + t262 * t354;
t12 = t50 * t260 + t262 * t355;
t9 = -t98 * t385 + t91 * t100 + t92 * t102 + t191 * t44 + t192 * t46 + (t280 * t42 - t417 * t98) * t262;
t8 = -t97 * t385 + t92 * t101 + t191 * t45 + t192 * t47 + t91 * t99 + (t280 * t43 - t417 * t97) * t262;
t7 = -t98 * t386 + t89 * t100 + t90 * t102 + t193 * t44 + t194 * t46 + (t279 * t42 + t416 * t98) * t262;
t6 = -t97 * t386 + t90 * t101 + t193 * t45 + t194 * t47 + t89 * t99 + (t279 * t43 + t416 * t97) * t262;
t4 = qJD(1) * t355 + t8 * t279 - t280 * t9;
t3 = qJD(1) * t354 + t6 * t279 - t280 * t7;
t2 = (-qJD(3) * t355 + t17) * t260 + (-qJD(1) * t18 + qJD(3) * t50 + t279 * t9 + t280 * t8) * t262;
t1 = (-qJD(3) * t354 + t16) * t260 + (-qJD(1) * t19 + qJD(3) * t51 + t279 * t7 + t280 * t6) * t262;
t20 = [(t25 * t59 + t26 * t58) * t492 + (t33 * t65 + t34 * t64) * t493 + (t103 * t57 + t104 * t56) * t494 + (t142 * t85 + t143 * t84) * t495 + 0.2e1 * m(3) * (t140 * t163 + t141 * t162) + t365 - t261 * t158 * t409 - t155 * t442 + (-t154 * t275 - t261 * t94 - t466) * t262 + (Icges(6,3) * t260 - t262 * t336 + t332 + t342 + t520) * t414 + (Icges(6,3) * t262 + t165 * t275 + t166 * t273 + t260 * t336 + t334 + t346 - t521) * t415; m(7) * (qJD(1) * t351 - t25 * t280 + t279 * t26) + m(5) * (qJD(1) * t329 + t279 * t57 - t280 * t56) + m(6) * (qJD(1) * t349 + t279 * t34 - t280 * t33) + m(4) * (qJD(1) * t500 + t279 * t85 - t280 * t84) + m(3) * (-t140 * t280 + t279 * t141 + (t162 * t280 + t163 * t279) * qJD(1)); 0; (t165 * t516 + t166 * t515 - t203 * t154 / 0.2e1 - t204 * t155 / 0.2e1 + t381 * t280 - t401) * t280 + (t151 * t486 + t152 * t485 + t154 * t484 + t155 * t483 + t279 * t381 + t402) * t279 + m(4) * ((-t279 * t84 - t280 * t85) * t225 - t500 * t214) + m(7) * (t25 * t80 + t26 * t81 + t27 * t58 + t28 * t59) + m(6) * (t107 * t33 + t108 * t34 + t52 * t64 + t53 * t65) + m(5) * (t103 * t82 + t104 * t83 + t145 * t56 + t146 * t57) + ((Icges(6,5) * t515 + Icges(6,6) * t516 + t174 * t512 + t177 * t472 + t300 * t519 + t517 * t514) * t280 + (Icges(6,5) * t152 / 0.2e1 + Icges(6,6) * t151 / 0.2e1 + t511 * t519 + t517 * t513 + (t173 + t497) * t512) * t279) * t260 + ((t172 * t512 + t175 * t472 + t68 * t478 + t70 * t480 + t514 * t518) * t280 + (t69 * t479 + t71 * t481 + t518 * t513 + (t171 + t498) * t512) * t279) * t262 + ((t279 * t367 - t280 * t366) * t262 + (t279 * t290 + t280 * t289) * t260) * qJD(3) + ((-t143 * t475 + t165 * t484 + t166 * t483 + t260 * t367 - t262 * t290 - t400) * t280 + (t142 * t475 + t203 * t486 + t204 * t485 + t260 * t366 + t262 * t289 + t399) * t279) * qJD(1); m(5) * (t82 * t279 - t280 * t83 + (t145 * t279 + t146 * t280) * qJD(1)) + m(6) * (t52 * t279 - t280 * t53 + (t107 * t279 + t108 * t280) * qJD(1)) + m(7) * (t27 * t279 - t28 * t280 + (t279 * t80 + t280 * t81) * qJD(1)); (t24 * t5 + t27 * t81 + t28 * t80) * t492 - t280 * t3 + t279 * t4 + (t145 * t83 + t146 * t82 + t63 * t23) * t494 + (t107 * t53 + t108 * t52 + t41 * t15) * t493 + ((t279 * t185 + t186 * t280) * ((qJD(1) * t185 + t286) * t280 + (-t279 * t308 + (-t186 + t501) * qJD(1)) * t279) + t420 * t225 * t214) * t495 - t280 * ((t280 * t119 + (t73 + t504) * qJD(1)) * t280 + (t72 * qJD(1) + (-t172 * t414 - t174 * t415 + t419) * t279 + (-t118 + (t171 * t262 + t173 * t260) * qJD(3) - t322 * qJD(1)) * t280) * t279) + t279 * ((t279 * t129 + (t77 - t507) * qJD(1)) * t279 + (t76 * qJD(1) + (t414 * t498 + t415 * t497) * t280 + (-t128 + (t175 * t262 + t177 * t260) * qJD(3) + (t179 + t320) * qJD(1)) * t279) * t280) - t280 * ((t280 * t128 + (t78 - t505) * qJD(1)) * t280 + (t79 * qJD(1) + (t175 * t414 + t177 * t415 + t418) * t279 + (-t129 + (t260 * t497 + t262 * t498) * qJD(3) + t321 * qJD(1)) * t280) * t279) - t280 * ((-t149 * t115 - t150 * t117 - t203 * t68 - t204 * t70 + (t39 - t395) * qJD(1)) * t280 + (t149 * t114 + t150 * t116 + t203 * t69 + t204 * t71 + (t40 + t397) * qJD(1)) * t279) + t279 * ((t279 * t118 + (t74 + t506) * qJD(1)) * t279 + (t75 * qJD(1) + (t171 * t414 + t173 * t415) * t280 + (-t119 + (-t172 * t262 - t174 * t260) * qJD(3) + (t170 - t323) * qJD(1)) * t279) * t280) + t279 * ((t151 * t114 + t152 * t116 + t201 * t69 + t202 * t71 + (t38 - t398) * qJD(1)) * t279 + (-t151 * t115 - t152 * t117 - t201 * t68 - t202 * t70 + (t37 + t396) * qJD(1)) * t280) + (t19 + (-t40 - t72 - t79) * t280 + (t39 + t73 + t78) * t279) * t417 + (t18 + (-t38 - t74 - t77) * t280 + (t37 + t75 + t76) * t279) * t416; 0.2e1 * (t329 * t491 + t349 * t489 + t351 * t487) * t414 + 0.2e1 * ((-t103 * t417 + t104 * t416 + t279 * t56 + t280 * t57) * t491 + t282) * t260; 0; 0.2e1 * (t296 * t487 + (t145 * t413 + t146 * t412 - t23) * t491 + t294 * t489) * t262 + 0.2e1 * ((qJD(3) * t63 + t145 * t416 - t146 * t417 + t279 * t83 + t280 * t82) * t491 + t281) * t260; 0.4e1 * (t491 + t408) * t137; 0.2e1 * (t349 * t490 + t351 * t488) * t415 + 0.2e1 * t282 * t262; 0; 0.2e1 * (t294 * t490 + t296 * t488) * t260 + 0.2e1 * t281 * t262; 0.2e1 * t408 * (-t260 ^ 2 + t262 ^ 2) * t509 * qJD(3); -0.4e1 * t408 * t137; m(7) * (t21 * t59 + t22 * t58 + t25 * t62 + t26 * t61) + (-t279 * t399 + t280 * t400) * t415 + (t402 * t280 + t401 * t279 + (t279 * t400 + t280 * t399) * qJD(1)) * t262 + t471; m(7) * (qJD(1) * t350 - t21 * t280 + t22 * t279); m(7) * (t14 * t24 + t21 * t80 + t22 * t81 + t61 * t27 + t62 * t28 - t54 * t5) + (t18 * t380 - t1 / 0.2e1 + t12 * t472 + (qJD(1) * t35 - t11) * t482) * t280 + (t19 * t380 + t13 * t472 + t2 / 0.2e1 + (qJD(1) * t36 + t10) * t482) * t279 + (t3 * t514 + t4 * t477 + t352 * t508 + (t19 * t477 - t279 * t18 / 0.2e1) * qJD(1)) * t262; m(7) * (t283 * t260 + t262 * t499); m(7) * (-t260 * t499 + t283 * t262); (-t14 * t54 + t21 * t62 + t22 * t61) * t492 + ((-t280 * t12 - t279 * t13 - t260 * t353) * qJD(3) + t471) * t260 + (t280 * t2 + t279 * t1 + t260 * (t10 * t280 + t11 * t279) + (t60 * t260 + t262 * t353) * qJD(3) + (-t279 * t12 + t280 * t13 - t260 * t352) * qJD(1)) * t262;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;

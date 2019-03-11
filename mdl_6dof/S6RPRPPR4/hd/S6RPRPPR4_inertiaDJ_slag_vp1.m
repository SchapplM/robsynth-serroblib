% Calculate time derivative of joint inertia matrix for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:47
% EndTime: 2019-03-09 02:47:11
% DurationCPUTime: 15.41s
% Computational Cost: add. (21261->836), mult. (33119->1147), div. (0->0), fcn. (34597->10), ass. (0->368)
t473 = Icges(6,2) + Icges(5,3);
t476 = t473 / 0.2e1;
t474 = Icges(6,4) + Icges(5,5);
t475 = -Icges(5,6) + Icges(6,6);
t472 = Icges(6,6) / 0.2e1 - Icges(5,6) / 0.2e1;
t279 = sin(qJ(1));
t442 = t279 / 0.2e1;
t281 = cos(qJ(1));
t441 = -t281 / 0.2e1;
t465 = -qJD(1) / 0.2e1;
t270 = pkin(9) + qJ(3);
t265 = sin(t270);
t266 = cos(t270);
t273 = sin(pkin(10));
t275 = cos(pkin(10));
t316 = Icges(6,5) * t275 + Icges(6,3) * t273;
t320 = Icges(5,4) * t275 - Icges(5,2) * t273;
t464 = ((t316 - t320) * t266 + t475 * t265) * qJD(3);
t323 = Icges(6,1) * t275 + Icges(6,5) * t273;
t324 = Icges(5,1) * t275 - Icges(5,4) * t273;
t463 = ((t323 + t324) * t266 + t474 * t265) * qJD(3);
t195 = -Icges(6,6) * t266 + t265 * t316;
t198 = -Icges(5,6) * t266 + t265 * t320;
t462 = t198 - t195;
t199 = -Icges(6,4) * t266 + t265 * t323;
t200 = -Icges(5,5) * t266 + t265 * t324;
t461 = t200 + t199;
t445 = m(7) / 0.2e1;
t446 = m(6) / 0.2e1;
t384 = t446 + t445;
t460 = 0.2e1 * t384;
t416 = t265 * t281;
t278 = sin(qJ(6));
t280 = cos(qJ(6));
t308 = t273 * t278 + t275 * t280;
t215 = t308 * t265;
t425 = Icges(4,4) * t266;
t322 = -Icges(4,2) * t265 + t425;
t207 = Icges(4,6) * t279 + t281 * t322;
t426 = Icges(4,4) * t265;
t326 = Icges(4,1) * t266 - t426;
t209 = Icges(4,5) * t279 + t281 * t326;
t311 = t207 * t265 - t209 * t266;
t298 = t311 * t279;
t206 = -Icges(4,6) * t281 + t279 * t322;
t208 = -Icges(4,5) * t281 + t279 * t326;
t312 = t206 * t265 - t208 * t266;
t299 = t312 * t281;
t409 = t279 * t273;
t411 = t275 * t281;
t231 = t266 * t409 + t411;
t408 = t279 * t275;
t412 = t273 * t281;
t232 = t266 * t408 - t412;
t417 = t265 * t279;
t131 = Icges(6,5) * t232 + Icges(6,6) * t417 + Icges(6,3) * t231;
t137 = Icges(5,4) * t232 - Icges(5,2) * t231 + Icges(5,6) * t417;
t458 = t131 - t137;
t233 = t266 * t412 - t408;
t234 = t266 * t411 + t409;
t132 = Icges(6,5) * t234 + Icges(6,6) * t416 + Icges(6,3) * t233;
t138 = Icges(5,4) * t234 - Icges(5,2) * t233 + Icges(5,6) * t416;
t457 = t132 - t138;
t139 = Icges(6,1) * t232 + Icges(6,4) * t417 + Icges(6,5) * t231;
t141 = Icges(5,1) * t232 - Icges(5,4) * t231 + Icges(5,5) * t417;
t456 = t139 + t141;
t140 = Icges(6,1) * t234 + Icges(6,4) * t416 + Icges(6,5) * t233;
t142 = Icges(5,1) * t234 - Icges(5,4) * t233 + Icges(5,5) * t416;
t455 = t140 + t142;
t454 = t273 * t475 + t474 * t275;
t269 = t279 * rSges(4,3);
t453 = -rSges(4,2) * t416 + t269;
t277 = -pkin(7) - qJ(2);
t276 = cos(pkin(9));
t263 = pkin(2) * t276 + pkin(1);
t434 = rSges(4,1) * t266;
t343 = -rSges(4,2) * t265 + t434;
t304 = -t263 - t343;
t169 = (rSges(4,3) - t277) * t281 + t304 * t279;
t413 = t266 * t281;
t213 = rSges(4,1) * t413 + t453;
t349 = t281 * t263 - t279 * t277;
t170 = t213 + t349;
t452 = t169 * t281 + t170 * t279;
t318 = Icges(4,5) * t266 - Icges(4,6) * t265;
t204 = -Icges(4,3) * t281 + t279 * t318;
t307 = rSges(3,1) * t276 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t430 = rSges(3,3) + qJ(2);
t194 = t279 * t430 + t281 * t307;
t451 = 2 * m(4);
t450 = 2 * m(5);
t449 = 0.2e1 * m(6);
t448 = 0.2e1 * m(7);
t271 = t279 ^ 2;
t272 = t281 ^ 2;
t447 = m(5) / 0.2e1;
t444 = -pkin(4) - pkin(5);
t443 = t266 / 0.2e1;
t439 = -rSges(6,1) - pkin(4);
t438 = -rSges(7,3) - pkin(8);
t437 = pkin(3) * t265;
t436 = pkin(3) * t266;
t435 = qJD(1) / 0.2e1;
t387 = qJD(3) * t279;
t359 = t265 * t387;
t391 = qJD(1) * t281;
t361 = t266 * t391;
t392 = qJD(1) * t279;
t180 = -t275 * t392 + (-t359 + t361) * t273;
t433 = t180 * rSges(6,3);
t432 = t231 * rSges(6,3);
t431 = -rSges(6,2) - qJ(4);
t429 = -rSges(5,3) - qJ(4);
t160 = t231 * t280 - t232 * t278;
t161 = t231 * t278 + t232 * t280;
t338 = -t161 * rSges(7,1) - t160 * rSges(7,2);
t97 = -rSges(7,3) * t417 - t338;
t428 = t232 * pkin(5) - pkin(8) * t417 + t97;
t225 = t234 * pkin(5);
t382 = pkin(8) * t416;
t162 = t233 * t280 - t234 * t278;
t163 = t233 * t278 + t234 * t280;
t406 = t163 * rSges(7,1) + t162 * rSges(7,2);
t98 = -rSges(7,3) * t416 + t406;
t427 = t225 - t382 + t98;
t419 = t265 * t273;
t418 = t265 * t275;
t415 = t266 * t273;
t414 = t266 * t275;
t410 = t277 * t281;
t309 = t273 * t280 - t275 * t278;
t214 = t309 * t265;
t126 = rSges(7,1) * t215 + rSges(7,2) * t214 + rSges(7,3) * t266;
t407 = pkin(5) * t418 + pkin(8) * t266 + t126;
t166 = t234 * pkin(4) + t233 * qJ(5);
t229 = pkin(3) * t413 + qJ(4) * t416;
t405 = -t166 - t229;
t404 = -t180 * qJ(5) - t231 * qJD(5);
t386 = qJD(3) * t281;
t358 = t265 * t386;
t179 = -qJD(1) * t232 - t275 * t358;
t362 = t265 * t392;
t403 = t179 * pkin(5) + pkin(8) * t362;
t337 = qJ(4) * t265 + t436;
t211 = qJD(3) * t337 - qJD(4) * t266;
t336 = pkin(4) * t275 + qJ(5) * t273;
t389 = qJD(3) * t266;
t402 = -qJD(5) * t419 - t336 * t389 - t211;
t340 = rSges(5,1) * t275 - rSges(5,2) * t273;
t401 = -(rSges(5,3) * t265 + t266 * t340) * qJD(3) - t211;
t202 = -rSges(5,3) * t266 + t265 * t340;
t243 = -qJ(4) * t266 + t437;
t400 = -t202 - t243;
t220 = t336 * t265;
t230 = t243 * t392;
t399 = t220 * t392 + t230;
t228 = t337 * t279;
t398 = t279 * t228 + t281 * t229;
t397 = -t220 - t243;
t354 = t266 * t386;
t385 = qJD(4) * t265;
t396 = qJ(4) * t354 + t281 * t385;
t268 = qJD(2) * t281;
t395 = t277 * t392 + t268;
t394 = t271 + t272;
t205 = Icges(4,3) * t279 + t281 * t318;
t393 = qJD(1) * t205;
t390 = qJD(3) * t265;
t388 = qJD(3) * t273;
t383 = -qJ(4) - t438;
t178 = qJD(1) * t231 + t273 * t358;
t79 = -qJD(6) * t163 - t178 * t280 - t179 * t278;
t80 = qJD(6) * t162 - t178 * t278 + t179 * t280;
t381 = t80 * rSges(7,1) + t79 * rSges(7,2) + rSges(7,3) * t362;
t91 = Icges(7,5) * t161 + Icges(7,6) * t160 - Icges(7,3) * t417;
t93 = Icges(7,4) * t161 + Icges(7,2) * t160 - Icges(7,6) * t417;
t95 = Icges(7,1) * t161 + Icges(7,4) * t160 - Icges(7,5) * t417;
t31 = t214 * t93 + t215 * t95 + t266 * t91;
t123 = Icges(7,5) * t215 + Icges(7,6) * t214 + Icges(7,3) * t266;
t124 = Icges(7,4) * t215 + Icges(7,2) * t214 + Icges(7,6) * t266;
t125 = Icges(7,1) * t215 + Icges(7,4) * t214 + Icges(7,5) * t266;
t35 = -t123 * t417 + t124 * t160 + t125 * t161;
t379 = -t31 / 0.2e1 - t35 / 0.2e1;
t92 = Icges(7,5) * t163 + Icges(7,6) * t162 - Icges(7,3) * t416;
t94 = Icges(7,4) * t163 + Icges(7,2) * t162 - Icges(7,6) * t416;
t96 = Icges(7,1) * t163 + Icges(7,4) * t162 - Icges(7,5) * t416;
t32 = t214 * t94 + t215 * t96 + t266 * t92;
t36 = -t123 * t416 + t162 * t124 + t163 * t125;
t378 = t36 / 0.2e1 + t32 / 0.2e1;
t133 = Icges(5,5) * t232 - Icges(5,6) * t231 + Icges(5,3) * t417;
t377 = t133 * t417;
t376 = t133 * t416;
t134 = Icges(5,5) * t234 - Icges(5,6) * t233 + Icges(5,3) * t416;
t375 = t134 * t417;
t374 = t134 * t416;
t135 = Icges(6,4) * t232 + Icges(6,2) * t417 + Icges(6,6) * t231;
t373 = t135 * t417;
t372 = t135 * t416;
t136 = Icges(6,4) * t234 + Icges(6,2) * t416 + Icges(6,6) * t233;
t371 = t136 * t417;
t370 = t136 * t416;
t251 = pkin(3) * t359;
t355 = t266 * t387;
t292 = t265 * t391 + t355;
t369 = t279 * (pkin(3) * t361 + qJ(4) * t292 + t279 * t385 - t251) + t281 * (-qJ(4) * t362 + (-t266 * t392 - t358) * pkin(3) + t396) + t228 * t391;
t368 = t179 * rSges(5,1) + t178 * rSges(5,2) + rSges(5,3) * t354;
t367 = t179 * rSges(6,1) + rSges(6,2) * t354 - t178 * rSges(6,3);
t339 = rSges(6,1) * t275 + rSges(6,3) * t273;
t366 = -(rSges(6,2) * t265 + t266 * t339) * qJD(3) + t402;
t201 = -rSges(6,2) * t266 + t265 * t339;
t365 = -t201 + t397;
t151 = t234 * rSges(6,1) + rSges(6,2) * t416 + t233 * rSges(6,3);
t152 = t234 * rSges(5,1) - t233 * rSges(5,2) + rSges(5,3) * t416;
t267 = qJD(2) * t279;
t364 = t267 + t396;
t363 = t251 + t395;
t360 = t265 * t389;
t353 = t195 / 0.2e1 - t198 / 0.2e1;
t352 = t199 / 0.2e1 + t200 / 0.2e1;
t351 = -t389 / 0.2e1;
t350 = -t263 - t436;
t154 = t400 * t281;
t348 = pkin(3) * t358;
t129 = -qJD(6) * t215 + t309 * t389;
t130 = qJD(6) * t214 + t308 * t389;
t86 = rSges(7,1) * t130 + rSges(7,2) * t129 - rSges(7,3) * t390;
t347 = -(pkin(5) * t414 - pkin(8) * t265) * qJD(3) - t86 + t402;
t346 = t397 - t407;
t219 = t231 * qJ(5);
t165 = t232 * pkin(4) + t219;
t345 = t279 * t165 + t281 * t166 + t398;
t122 = t365 * t281;
t181 = qJD(1) * t234 - t275 * t359;
t81 = -qJD(6) * t161 + t180 * t280 - t181 * t278;
t82 = qJD(6) * t160 + t180 * t278 + t181 * t280;
t344 = t82 * rSges(7,1) + t81 * rSges(7,2);
t244 = rSges(4,1) * t265 + rSges(4,2) * t266;
t342 = -t181 * rSges(5,1) + t180 * rSges(5,2);
t341 = -t232 * rSges(5,1) + t231 * rSges(5,2);
t27 = t160 * t93 + t161 * t95 - t417 * t91;
t28 = t160 * t94 + t161 * t96 - t417 * t92;
t17 = -t27 * t281 + t28 * t279;
t335 = t27 * t279 + t28 * t281;
t29 = t162 * t93 + t163 * t95 - t416 * t91;
t30 = t162 * t94 + t163 * t96 - t416 * t92;
t18 = t30 * t279 - t281 * t29;
t334 = t279 * t29 + t281 * t30;
t333 = t32 * t279 - t281 * t31;
t332 = -t279 * t31 - t281 * t32;
t288 = t265 * t383 + t350;
t51 = t232 * t444 + t279 * t288 - t219 + t338 - t410;
t306 = t349 + t229;
t286 = t166 + t306;
t52 = t416 * t438 + t225 + t286 + t406;
t331 = t279 * t52 + t281 * t51;
t63 = -t126 * t417 - t266 * t97;
t64 = t126 * t416 + t266 * t98;
t330 = t279 * t64 + t281 * t63;
t294 = t265 * t431 + t350;
t283 = t279 * t294 - t410;
t70 = t232 * t439 - t219 + t283 - t432;
t71 = t286 + t151;
t329 = t279 * t71 + t281 * t70;
t328 = t279 * t98 - t281 * t97;
t327 = t363 + t404;
t325 = Icges(4,1) * t265 + t425;
t321 = Icges(4,2) * t266 + t426;
t293 = t265 * t429 + t350;
t282 = t279 * t293 - t410;
t117 = t282 + t341;
t118 = t306 + t152;
t315 = t117 * t281 + t118 * t279;
t310 = t231 * t279 + t233 * t281;
t75 = t346 * t281;
t301 = t179 * pkin(4) - t178 * qJ(5) + t233 * qJD(5);
t305 = t281 * t301 + t165 * t391 + t279 * (t181 * pkin(4) - t404) + t369;
t300 = qJD(3) * t244;
t297 = qJD(3) * t325;
t296 = qJD(3) * t321;
t295 = qJD(3) * (-Icges(4,5) * t265 - Icges(4,6) * t266);
t291 = -t354 + t362;
t83 = Icges(7,5) * t130 + Icges(7,6) * t129 - Icges(7,3) * t390;
t84 = Icges(7,4) * t130 + Icges(7,2) * t129 - Icges(7,6) * t390;
t85 = Icges(7,1) * t130 + Icges(7,4) * t129 - Icges(7,5) * t390;
t287 = -t123 * t390 + t129 * t124 + t130 * t125 + t214 * t84 + t215 * t85 + t266 * t83;
t285 = t301 + t364;
t284 = rSges(4,2) * t362 + rSges(4,3) * t391 - t281 * t300;
t193 = -t279 * t307 + t281 * t430;
t238 = t343 * qJD(3);
t212 = -rSges(4,3) * t281 + t279 * t343;
t168 = -qJD(1) * t194 + t268;
t167 = qJD(1) * t193 + t267;
t153 = t400 * t279;
t150 = rSges(5,3) * t417 - t341;
t149 = t232 * rSges(6,1) + rSges(6,2) * t417 + t432;
t144 = t279 * t295 + t393;
t143 = -qJD(1) * t204 + t281 * t295;
t121 = t365 * t279;
t120 = t244 * t387 + (t281 * t304 - t269) * qJD(1) + t395;
t119 = t267 + (-t410 + (-t263 - t434) * t279) * qJD(1) + t284;
t116 = t279 * t205 - t281 * t311;
t115 = t279 * t204 - t299;
t114 = -t205 * t281 - t298;
t113 = -t204 * t281 - t279 * t312;
t112 = Icges(5,1) * t181 - Icges(5,4) * t180 + Icges(5,5) * t292;
t111 = Icges(5,1) * t179 + Icges(5,4) * t178 - Icges(5,5) * t291;
t110 = Icges(6,1) * t181 + Icges(6,4) * t292 + Icges(6,5) * t180;
t109 = Icges(6,1) * t179 - Icges(6,4) * t291 - Icges(6,5) * t178;
t108 = Icges(5,4) * t181 - Icges(5,2) * t180 + Icges(5,6) * t292;
t107 = Icges(5,4) * t179 + Icges(5,2) * t178 - Icges(5,6) * t291;
t102 = Icges(6,5) * t181 + Icges(6,6) * t292 + Icges(6,3) * t180;
t101 = Icges(6,5) * t179 - Icges(6,6) * t291 - Icges(6,3) * t178;
t88 = qJD(1) * t154 + t279 * t401;
t87 = t202 * t392 + t281 * t401 + t230;
t74 = t346 * t279;
t67 = t279 * t150 + t152 * t281 + t398;
t66 = (t389 * t429 - t385) * t279 + t293 * t391 + t342 + t363;
t65 = qJD(1) * t282 - t348 + t364 + t368;
t62 = qJD(1) * t122 + t279 * t366;
t61 = t201 * t392 + t281 * t366 + t399;
t60 = -t233 * t138 + t234 * t142 + t374;
t59 = -t233 * t137 + t234 * t141 + t376;
t58 = t233 * t132 + t234 * t140 + t370;
t57 = t233 * t131 + t234 * t139 + t372;
t56 = -t138 * t231 + t142 * t232 + t375;
t55 = -t137 * t231 + t141 * t232 + t377;
t54 = t132 * t231 + t140 * t232 + t371;
t53 = t131 * t231 + t139 * t232 + t373;
t49 = t328 * t265;
t48 = t123 * t266 + t124 * t214 + t125 * t215;
t47 = t279 * t149 + t151 * t281 + t345;
t46 = -t433 + t439 * t181 + (t389 * t431 - t385) * t279 + t294 * t391 + t327;
t45 = qJD(1) * t283 + t285 - t348 + t367;
t44 = -rSges(7,3) * t292 + t344;
t43 = -rSges(7,3) * t354 + t381;
t42 = Icges(7,1) * t82 + Icges(7,4) * t81 - Icges(7,5) * t292;
t41 = Icges(7,1) * t80 + Icges(7,4) * t79 + Icges(7,5) * t291;
t40 = Icges(7,4) * t82 + Icges(7,2) * t81 - Icges(7,6) * t292;
t39 = Icges(7,4) * t80 + Icges(7,2) * t79 + Icges(7,6) * t291;
t38 = Icges(7,5) * t82 + Icges(7,6) * t81 - Icges(7,3) * t292;
t37 = Icges(7,5) * t80 + Icges(7,6) * t79 + Icges(7,3) * t291;
t34 = qJD(1) * t75 + t279 * t347;
t33 = t281 * t347 + t392 * t407 + t399;
t26 = t279 * t428 + t281 * t427 + t345;
t25 = t279 * (rSges(5,3) * t355 - t342) + t281 * t368 + (t281 * t150 + (-t152 - t229) * t279) * qJD(1) + t369;
t24 = t444 * t181 + (t383 * t389 - t385) * t279 + t288 * t391 + t327 - t344;
t23 = (t266 * t438 - t437) * t386 + (-t410 + (-t263 - t337) * t279) * qJD(1) + t285 + t381 + t403;
t22 = (-t126 * t387 - t44) * t266 + (qJD(3) * t97 - t126 * t391 - t279 * t86) * t265;
t21 = (t126 * t386 + t43) * t266 + (-qJD(3) * t98 - t126 * t392 + t281 * t86) * t265;
t20 = t287 * t266;
t19 = t279 * (t181 * rSges(6,1) + rSges(6,2) * t355 + t433) + t281 * t367 + (t281 * t149 + (-t151 + t405) * t279) * qJD(1) + t305;
t16 = -t123 * t292 + t81 * t124 + t82 * t125 + t160 * t84 + t161 * t85 - t417 * t83;
t15 = t123 * t291 + t79 * t124 + t80 * t125 + t162 * t84 + t163 * t85 - t416 * t83;
t14 = t328 * t389 + (t279 * t43 - t281 * t44 + (t279 * t97 + t281 * t98) * qJD(1)) * t265;
t13 = -t265 * t334 + t36 * t266;
t12 = -t265 * t335 + t35 * t266;
t11 = t129 * t94 + t130 * t96 + t214 * t39 + t215 * t41 + t266 * t37 - t390 * t92;
t10 = t129 * t93 + t130 * t95 + t214 * t40 + t215 * t42 + t266 * t38 - t390 * t91;
t9 = (-pkin(8) * t354 + t403 + t43) * t281 + (t181 * pkin(5) - pkin(8) * t355 + t44) * t279 + (t428 * t281 + (-t382 + t405 - t427) * t279) * qJD(1) + t305;
t8 = -t92 * t355 + t160 * t39 + t161 * t41 + t81 * t94 + t82 * t96 + (-t279 * t37 - t391 * t92) * t265;
t7 = -t91 * t355 + t160 * t40 + t161 * t42 + t81 * t93 + t82 * t95 + (-t279 * t38 - t391 * t91) * t265;
t6 = -t92 * t354 + t162 * t39 + t163 * t41 + t79 * t94 + t80 * t96 + (-t281 * t37 + t392 * t92) * t265;
t5 = -t91 * t354 + t162 * t40 + t163 * t42 + t79 * t93 + t80 * t95 + (-t281 * t38 + t392 * t91) * t265;
t4 = qJD(1) * t335 + t8 * t279 - t281 * t7;
t3 = qJD(1) * t334 + t6 * t279 - t281 * t5;
t2 = (-qJD(3) * t335 + t16) * t266 + (qJD(1) * t17 - qJD(3) * t35 - t279 * t7 - t281 * t8) * t265;
t1 = (-qJD(3) * t334 + t15) * t266 + (qJD(1) * t18 - qJD(3) * t36 - t279 * t5 - t281 * t6) * t265;
t50 = [t287 + (t23 * t52 + t24 * t51) * t448 + (t45 * t71 + t46 * t70) * t449 + (t117 * t66 + t118 * t65) * t450 + (t119 * t170 + t120 * t169) * t451 + 0.2e1 * m(3) * (t167 * t194 + t168 * t193) + t464 * t419 + t463 * t418 - t462 * t266 * t388 + (t265 * t454 - t266 * t473 - t321 + t326) * t390 + (-t265 * t473 - t454 * t266 + t461 * t275 + t322 + t325) * t389; m(7) * (qJD(1) * t331 - t23 * t281 + t279 * t24) + m(6) * (qJD(1) * t329 + t279 * t46 - t281 * t45) + m(5) * (qJD(1) * t315 + t279 * t66 - t281 * t65) + m(4) * (qJD(1) * t452 - t119 * t281 + t279 * t120) + m(3) * (-t167 * t281 + t279 * t168 + (t193 * t281 + t194 * t279) * qJD(1)); 0; m(4) * ((-t119 * t279 - t120 * t281) * t244 - t452 * t238) + m(7) * (t23 * t74 + t24 * t75 + t33 * t51 + t34 * t52) + m(6) * (t121 * t45 + t122 * t46 + t61 * t70 + t62 * t71) + m(5) * (t117 * t87 + t118 * t88 + t153 * t65 + t154 * t66) + ((t207 * t465 + t296 * t442 + t474 * t181 / 0.2e1 + t292 * t476 + t472 * t180) * t281 + (t206 * t465 + t296 * t441 - t474 * t179 / 0.2e1 + t291 * t476 + t472 * t178) * t279) * t266 + ((t233 * t353 + t234 * t352 + t378) * t281 + (t231 * t353 + t232 * t352 - t379) * t279 + m(4) * (t169 * t279 - t170 * t281) * t244 + ((t207 / 0.2e1 - t136 / 0.2e1 - t134 / 0.2e1) * t281 + (t206 / 0.2e1 - t135 / 0.2e1 - t133 / 0.2e1) * t279) * t266 + (t273 * t457 + t275 * t455 + t209) * t416 / 0.2e1) * qJD(1) + ((t272 / 0.2e1 + t271 / 0.2e1) * t318 + t299 / 0.2e1 - t298 / 0.2e1) * qJD(3) + (t11 + t15 + (t455 * t414 + t457 * t415) * qJD(3) + (-t281 * t297 + (t109 + t111) * t275 + (t101 - t107) * t273 + (t134 + t136) * qJD(3) + (t273 * t458 + t275 * t456) * qJD(1)) * t265 + t463 * t234 + t464 * t233 + t461 * t179 + t462 * t178) * t442 + (t16 + t10 + (qJD(1) * t209 - t279 * t297 + (t110 + t112) * t275 + (t102 - t108) * t273) * t265 + (t458 * t415 + t456 * t414 + (t133 + t135) * t265) * qJD(3) + t463 * t232 + t464 * t231 + t461 * t181 - t462 * t180) * t441; m(5) * (t87 * t279 - t281 * t88 + (t153 * t279 + t154 * t281) * qJD(1)) + m(6) * (t61 * t279 - t281 * t62 + (t121 * t279 + t122 * t281) * qJD(1)) + m(7) * (t33 * t279 - t281 * t34 + (t279 * t74 + t281 * t75) * qJD(1)); -t281 * t4 + t279 * t3 + (t26 * t9 + t33 * t75 + t34 * t74) * t448 + (t121 * t62 + t122 * t61 + t19 * t47) * t449 + (t153 * t88 + t154 * t87 + t25 * t67) * t450 + t279 * ((t233 * t101 + t234 * t109 - t178 * t132 + t179 * t140 + (t57 - t371) * qJD(1)) * t279 + (-t233 * t102 - t234 * t110 + t178 * t131 - t179 * t139 + (t58 + t373) * qJD(1)) * t281) + ((t279 * t212 + t213 * t281) * ((qJD(1) * t212 + t284) * t281 + (-t279 * t300 + (-t213 + t453) * qJD(1)) * t279) + t394 * t244 * t238) * t451 - t281 * ((-t231 * t102 - t232 * t110 - t180 * t131 - t181 * t139 + (t54 - t372) * qJD(1)) * t281 + (t231 * t101 + t232 * t109 + t180 * t132 + t181 * t140 + (t53 + t370) * qJD(1)) * t279) + t279 * ((t279 * t143 + (t115 + t298) * qJD(1)) * t279 + (t116 * qJD(1) + (t206 * t389 + t208 * t390) * t281 + (-t144 + (-t207 * t266 - t209 * t265) * qJD(3) + (t205 - t312) * qJD(1)) * t279) * t281) - t281 * ((t144 * t281 + (t114 + t299) * qJD(1)) * t281 + (t113 * qJD(1) + (-t207 * t389 - t209 * t390 + t393) * t279 + (-t143 + (t206 * t266 + t208 * t265) * qJD(3) - t311 * qJD(1)) * t281) * t279) - t281 * ((t231 * t108 - t232 * t112 + t180 * t137 - t181 * t141 + (t56 - t376) * qJD(1)) * t281 + (-t231 * t107 + t232 * t111 - t180 * t138 + t181 * t142 + (t55 + t374) * qJD(1)) * t279) + t279 * ((-t233 * t107 + t234 * t111 + t178 * t138 + t179 * t142 + (t59 - t375) * qJD(1)) * t279 + (t233 * t108 - t234 * t112 - t178 * t137 - t179 * t141 + (t60 + t377) * qJD(1)) * t281) + (t17 + (-t113 - t53 - t55) * t281 + (t114 + t54 + t56) * t279) * t392 + (t18 + (-t115 - t57 - t59) * t281 + (t116 + t58 + t60) * t279) * t391; 0.2e1 * (t315 * t447 + t329 * t446 + t331 * t445) * t389 + 0.2e1 * ((t23 * t279 + t24 * t281 + t391 * t52 - t392 * t51) * t445 + (t279 * t45 + t281 * t46 + t391 * t71 - t392 * t70) * t446 + (-t117 * t392 + t118 * t391 + t279 * t65 + t281 * t66) * t447) * t265; 0; 0.2e1 * ((t386 * t75 + t387 * t74 - t9) * t445 + (t121 * t387 + t122 * t386 - t19) * t446 + (t153 * t387 + t154 * t386 - t25) * t447) * t266 + 0.2e1 * ((qJD(3) * t26 + t279 * t34 + t281 * t33 + t391 * t74 - t392 * t75) * t445 + (qJD(3) * t47 + t121 * t391 - t122 * t392 + t279 * t62 + t281 * t61) * t446 + (qJD(3) * t67 + t153 * t391 - t154 * t392 + t279 * t88 + t281 * t87) * t447) * t265; 0.4e1 * (t447 + t384) * (-0.1e1 + t394) * t360; m(7) * (-t178 * t51 + t180 * t52 + t23 * t231 + t233 * t24) + m(6) * (-t178 * t70 + t180 * t71 + t231 * t45 + t233 * t46); (qJD(1) * t310 - t178 * t279 - t180 * t281) * t460; m(7) * (-t178 * t75 + t180 * t74 + t231 * t34 + t233 * t33 + (t26 * t389 + t265 * t9) * t273) + m(6) * (t121 * t180 - t122 * t178 + t231 * t62 + t233 * t61 + (t19 * t265 + t389 * t47) * t273); ((t310 - t415) * t389 + (t265 * t388 - t178 * t281 + t180 * t279 + (t231 * t281 - t233 * t279) * qJD(1)) * t265) * t460; 0.4e1 * t384 * (t273 ^ 2 * t360 - t178 * t233 + t180 * t231); m(7) * (t21 * t52 + t22 * t51 + t23 * t64 + t24 * t63) + t20 + (t279 * t379 - t281 * t378) * t389 + (-qJD(3) * t48 + (-t15 / 0.2e1 - t11 / 0.2e1) * t281 + (-t10 / 0.2e1 - t16 / 0.2e1) * t279 + (t279 * t378 + t281 * t379) * qJD(1)) * t265; m(7) * (qJD(1) * t330 - t21 * t281 + t22 * t279); m(7) * (t14 * t26 + t21 * t74 + t22 * t75 + t33 * t63 + t34 * t64 + t49 * t9) + (t18 * t351 + (qJD(1) * t32 - t10) * t443 - t2 / 0.2e1 + t13 * t435) * t281 + (t17 * t351 + (qJD(1) * t31 + t11) * t443 + t12 * t435 + t1 / 0.2e1) * t279 + (-t279 * t4 / 0.2e1 + t3 * t441 - qJD(3) * t333 / 0.2e1 + (t17 * t441 + t18 * t442) * qJD(1)) * t265; m(7) * ((qJD(3) * t330 - t14) * t266 + (qJD(3) * t49 + t21 * t279 + t22 * t281 + (-t279 * t63 + t281 * t64) * qJD(1)) * t265); m(7) * (-t178 * t63 + t180 * t64 + t21 * t231 + t22 * t233 + (t14 * t265 + t389 * t49) * t273); (t14 * t49 + t21 * t64 + t22 * t63) * t448 + (t20 + (-t279 * t12 - t281 * t13 + t266 * t332) * qJD(3)) * t266 + (-t281 * t1 - t279 * t2 + t266 * (-t10 * t279 - t11 * t281) + (-t265 * t332 - 0.2e1 * t48 * t266) * qJD(3) + (-t281 * t12 + t279 * t13 + t266 * t333) * qJD(1)) * t265;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t50(1) t50(2) t50(4) t50(7) t50(11) t50(16); t50(2) t50(3) t50(5) t50(8) t50(12) t50(17); t50(4) t50(5) t50(6) t50(9) t50(13) t50(18); t50(7) t50(8) t50(9) t50(10) t50(14) t50(19); t50(11) t50(12) t50(13) t50(14) t50(15) t50(20); t50(16) t50(17) t50(18) t50(19) t50(20) t50(21);];
Mq  = res;

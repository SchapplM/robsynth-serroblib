% Calculate time derivative of joint inertia matrix for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR7_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:33:09
% EndTime: 2019-03-09 02:33:23
% DurationCPUTime: 7.92s
% Computational Cost: add. (24834->684), mult. (24674->951), div. (0->0), fcn. (22995->10), ass. (0->361)
t280 = sin(qJ(1));
t282 = cos(qJ(1));
t272 = pkin(10) + qJ(4);
t262 = qJ(5) + t272;
t257 = cos(t262);
t256 = sin(t262);
t435 = Icges(6,4) * t256;
t319 = Icges(6,2) * t257 + t435;
t179 = Icges(6,6) * t282 + t280 * t319;
t434 = Icges(6,4) * t257;
t323 = Icges(6,1) * t256 + t434;
t181 = Icges(6,5) * t282 + t280 * t323;
t307 = t179 * t257 + t181 * t256;
t483 = t280 * t307;
t273 = qJD(4) + qJD(5);
t407 = t273 * t280;
t366 = t257 * t407;
t380 = qJD(1) * t282;
t482 = t256 * t380 + t366;
t454 = rSges(7,3) + pkin(9);
t356 = t454 * t257;
t450 = pkin(5) * t256;
t481 = -t356 + t450;
t281 = cos(qJ(6));
t279 = sin(qJ(6));
t432 = Icges(7,4) * t281;
t318 = -Icges(7,2) * t279 + t432;
t413 = t256 * t273;
t433 = Icges(7,4) * t279;
t111 = -t318 * t413 + (Icges(7,6) * t273 + (-Icges(7,2) * t281 - t433) * qJD(6)) * t257;
t169 = Icges(7,6) * t256 + t257 * t318;
t322 = Icges(7,1) * t281 - t433;
t170 = Icges(7,5) * t256 + t257 * t322;
t480 = -t111 * t279 + (-t169 * t281 - t170 * t279) * qJD(6);
t399 = t281 * t282;
t401 = t280 * t279;
t207 = -t256 * t401 + t399;
t400 = t280 * t281;
t404 = t279 * t282;
t208 = t256 * t400 + t404;
t410 = t257 * t280;
t138 = Icges(7,5) * t208 + Icges(7,6) * t207 - Icges(7,3) * t410;
t140 = Icges(7,4) * t208 + Icges(7,2) * t207 - Icges(7,6) * t410;
t142 = Icges(7,1) * t208 + Icges(7,4) * t207 - Icges(7,5) * t410;
t54 = -t138 * t410 + t140 * t207 + t142 * t208;
t209 = t256 * t404 + t400;
t210 = -t256 * t399 + t401;
t409 = t257 * t282;
t139 = Icges(7,5) * t210 + Icges(7,6) * t209 + Icges(7,3) * t409;
t141 = Icges(7,4) * t210 + Icges(7,2) * t209 + Icges(7,6) * t409;
t143 = Icges(7,1) * t210 + Icges(7,4) * t209 + Icges(7,5) * t409;
t55 = -t139 * t410 + t141 * t207 + t143 * t208;
t39 = t55 * t280 + t282 * t54;
t315 = Icges(6,5) * t256 + Icges(6,6) * t257;
t177 = Icges(6,3) * t282 + t280 * t315;
t89 = t177 * t282 + t483;
t466 = -Icges(6,5) * t280 + t282 * t323;
t468 = -Icges(6,6) * t280 + t282 * t319;
t306 = -t256 * t466 - t257 * t468;
t296 = t306 * t280;
t470 = -Icges(6,3) * t280 + t282 * t315;
t90 = -t282 * t470 + t296;
t479 = -t90 * t280 - t282 * t89 - t39;
t220 = rSges(6,1) * t257 - rSges(6,2) * t256;
t405 = t273 * t282;
t478 = t220 * t405;
t261 = cos(t272);
t260 = sin(t272);
t437 = Icges(5,4) * t260;
t320 = Icges(5,2) * t261 + t437;
t189 = Icges(5,6) * t282 + t280 * t320;
t436 = Icges(5,4) * t261;
t324 = Icges(5,1) * t260 + t436;
t191 = Icges(5,5) * t282 + t280 * t324;
t305 = t189 * t261 + t191 * t260;
t477 = t282 * t305;
t476 = t282 * t307;
t333 = rSges(6,1) * t256 + rSges(6,2) * t257;
t475 = t282 * t333;
t446 = rSges(5,1) * t260;
t334 = rSges(5,2) * t261 + t446;
t297 = t282 * t334;
t330 = rSges(7,1) * t281 - rSges(7,2) * t279;
t113 = -t330 * t413 + (rSges(7,3) * t273 + (-rSges(7,1) * t279 - rSges(7,2) * t281) * qJD(6)) * t257;
t171 = rSges(7,3) * t256 + t257 * t330;
t474 = t280 * t113 + t171 * t380;
t276 = sin(pkin(10));
t452 = pkin(3) * t276;
t254 = t282 * t452;
t278 = -pkin(7) - qJ(3);
t388 = t280 * t278 + t254;
t445 = rSges(5,2) * t260;
t228 = rSges(5,1) * t261 - t445;
t253 = t278 * t380;
t265 = qJD(2) * t282;
t340 = -qJD(3) * t280 + t265;
t376 = qJD(4) * t282;
t456 = -rSges(5,3) - pkin(1);
t100 = t253 + t228 * t376 + (t456 * t282 + (-qJ(2) - t334 - t452) * t280) * qJD(1) + t340;
t377 = qJD(4) * t280;
t352 = t260 * t377;
t387 = qJ(2) * t380 + qJD(2) * t280;
t357 = qJD(3) * t282 + t387;
t351 = t261 * t377;
t353 = t261 * t380;
t358 = rSges(5,1) * t351 + rSges(5,2) * t353 + t380 * t446;
t99 = -rSges(5,2) * t352 + (t254 + (t278 + t456) * t280) * qJD(1) + t357 + t358;
t473 = t280 * t100 - t282 * t99;
t267 = t282 * qJ(2);
t156 = t280 * t456 + t267 + t297 + t388;
t269 = t282 * rSges(5,3);
t408 = t261 * t280;
t196 = rSges(5,2) * t408 + t280 * t446 + t269;
t373 = t280 * t452;
t270 = t282 * pkin(1);
t386 = t280 * qJ(2) + t270;
t157 = -t278 * t282 + t196 + t373 + t386;
t472 = -t156 * t280 + t157 * t282;
t338 = qJD(1) * t256 + qJD(6);
t365 = t257 * t405;
t471 = t280 * t338 - t365;
t316 = Icges(5,5) * t260 + Icges(5,6) * t261;
t469 = -Icges(5,3) * t280 + t282 * t316;
t467 = -Icges(5,6) * t280 + t282 * t320;
t465 = -Icges(5,5) * t280 + t282 * t324;
t194 = t319 * t273;
t195 = t323 * t273;
t217 = Icges(6,5) * t257 - Icges(6,6) * t256;
t218 = -Icges(6,2) * t256 + t434;
t219 = Icges(6,1) * t257 - t435;
t463 = t256 * (t218 * t273 + t195) - t257 * (t219 * t273 - t194) + qJD(1) * t217;
t462 = 2 * m(5);
t461 = 2 * m(6);
t460 = 2 * m(7);
t274 = t280 ^ 2;
t275 = t282 ^ 2;
t459 = t280 / 0.2e1;
t458 = t282 / 0.2e1;
t457 = rSges(3,2) - pkin(1);
t455 = -rSges(6,3) - pkin(1);
t453 = m(5) * t228;
t451 = pkin(4) * t261;
t449 = pkin(5) * t257;
t448 = t280 * pkin(1);
t314 = Icges(7,5) * t281 - Icges(7,6) * t279;
t110 = -t314 * t413 + (Icges(7,3) * t273 + (-Icges(7,5) * t279 - Icges(7,6) * t281) * qJD(6)) * t257;
t112 = -t322 * t413 + (Icges(7,5) * t273 + (-Icges(7,1) * t279 - t432) * qJD(6)) * t257;
t168 = Icges(7,3) * t256 + t257 * t314;
t406 = t273 * t281;
t411 = t257 * t273;
t422 = t169 * t279;
t289 = t257 * t281 * t112 + t168 * t411 + t413 * t422 + (-t170 * t406 + t110) * t256;
t84 = t168 * t256 + (t170 * t281 - t422) * t257;
t447 = (t257 * t480 + t289) * t256 + t84 * t411;
t312 = t140 * t279 - t142 * t281;
t339 = qJD(6) * t256 + qJD(1);
t129 = -t339 * t400 + (-t282 * t338 - t366) * t279;
t130 = t338 * t399 + (t257 * t406 - t279 * t339) * t280;
t354 = t257 * t380;
t368 = t256 * t407;
t292 = -t354 + t368;
t70 = Icges(7,5) * t130 + Icges(7,6) * t129 + Icges(7,3) * t292;
t72 = Icges(7,4) * t130 + Icges(7,2) * t129 + Icges(7,6) * t292;
t74 = Icges(7,1) * t130 + Icges(7,4) * t129 + Icges(7,5) * t292;
t20 = (t273 * t312 + t70) * t256 + (t138 * t273 - t279 * t72 + t281 * t74 + (-t140 * t281 - t142 * t279) * qJD(6)) * t257;
t444 = t20 * t282;
t311 = t141 * t279 - t143 * t281;
t302 = t282 * t339;
t127 = -t279 * t471 + t281 * t302;
t128 = t279 * t302 + t281 * t471;
t367 = t256 * t405;
t381 = qJD(1) * t280;
t291 = -t257 * t381 - t367;
t69 = Icges(7,5) * t128 + Icges(7,6) * t127 + Icges(7,3) * t291;
t71 = Icges(7,4) * t128 + Icges(7,2) * t127 + Icges(7,6) * t291;
t73 = Icges(7,1) * t128 + Icges(7,4) * t127 + Icges(7,5) * t291;
t21 = (t273 * t311 + t69) * t256 + (t139 * t273 - t279 * t71 + t281 * t73 + (-t141 * t281 - t143 * t279) * qJD(6)) * t257;
t443 = t21 * t280;
t442 = t280 * rSges(5,3);
t268 = t282 * rSges(6,3);
t440 = rSges(4,3) + qJ(3);
t332 = t128 * rSges(7,1) + t127 * rSges(7,2);
t75 = rSges(7,3) * t291 + t332;
t439 = (t75 + t291 * pkin(9) + (t256 * t381 - t365) * pkin(5)) * t282;
t359 = pkin(5) * t482 + pkin(9) * t368;
t364 = t130 * rSges(7,1) + t129 * rSges(7,2) + rSges(7,3) * t368;
t76 = -rSges(7,3) * t354 + t364;
t438 = pkin(9) * t354 - t359 - t76;
t250 = pkin(4) * t408;
t172 = t220 * t280 + t250;
t419 = t172 * t280;
t418 = t189 * t260;
t417 = t467 * t260;
t416 = t191 * t261;
t415 = t465 * t261;
t198 = t333 * t273;
t414 = t198 * t282;
t412 = t256 * t280;
t167 = t280 * t171;
t402 = t280 * t217;
t331 = -t210 * rSges(7,1) - t209 * rSges(7,2);
t146 = rSges(7,3) * t409 - t331;
t134 = t282 * t146;
t235 = pkin(4) * t260 + t452;
t398 = -qJ(2) - t235;
t337 = pkin(9) * t257 - t450;
t201 = t337 * t282;
t397 = t282 * t201 + t134;
t392 = t208 * rSges(7,1) + t207 * rSges(7,2);
t145 = -rSges(7,3) * t410 + t392;
t247 = pkin(5) * t412;
t396 = pkin(9) * t410 - t145 - t247;
t395 = -t146 - t201;
t223 = t280 * t235;
t271 = -pkin(8) + t278;
t162 = -t373 + t223 + (-t271 + t278) * t282;
t184 = rSges(6,1) * t412 + rSges(6,2) * t410 + t268;
t394 = -t162 - t184;
t222 = pkin(9) * t256 + t449;
t393 = -t171 - t222;
t136 = t280 * t222 + t167;
t391 = -t282 * t235 - t280 * t271;
t371 = pkin(4) * t376;
t390 = qJD(1) * t250 + t260 * t371;
t389 = -t261 * t371 - t271 * t380;
t385 = t274 + t275;
t384 = qJD(1) * t177;
t187 = Icges(5,3) * t282 + t280 * t316;
t383 = qJD(1) * t187;
t379 = qJD(4) * t260;
t378 = qJD(4) * t261;
t374 = -pkin(1) - t440;
t372 = rSges(6,2) * t413;
t62 = t139 * t256 - t257 * t311;
t80 = t168 * t409 + t209 * t169 + t210 * t170;
t370 = -t62 / 0.2e1 - t80 / 0.2e1;
t61 = t138 * t256 - t257 * t312;
t79 = -t168 * t410 + t169 * t207 + t170 * t208;
t369 = t79 / 0.2e1 + t61 / 0.2e1;
t363 = -t162 + t396;
t362 = pkin(4) * t351 + t235 * t380 + t271 * t381;
t361 = t267 - t391;
t360 = rSges(6,1) * t482 + rSges(6,2) * t354;
t163 = t171 * t381;
t350 = -t381 / 0.2e1;
t349 = t380 / 0.2e1;
t348 = t198 * t385;
t216 = t334 * qJD(4);
t347 = t216 * t385;
t122 = qJD(1) * t179 - t218 * t405;
t346 = t273 * t466 - t122;
t123 = qJD(1) * t468 + t218 * t407;
t345 = t181 * t273 + t123;
t124 = qJD(1) * t181 - t219 * t405;
t344 = -t273 * t468 - t124;
t125 = qJD(1) * t466 + t219 * t407;
t343 = -t179 * t273 + t125;
t199 = t337 * t273;
t67 = t280 * t199 + t222 * t380 + t474;
t16 = t129 * t140 + t130 * t142 + t138 * t292 + t207 * t72 + t208 * t74 - t410 * t70;
t17 = t129 * t141 + t130 * t143 + t139 * t292 + t207 * t71 + t208 * t73 - t410 * t69;
t329 = t280 * t54 - t282 * t55;
t10 = -qJD(1) * t329 + t16 * t282 + t17 * t280;
t120 = -t217 * t405 + t384;
t121 = qJD(1) * t470 + t273 * t402;
t56 = t138 * t409 + t209 * t140 + t210 * t142;
t57 = t139 * t409 + t209 * t141 + t210 * t143;
t40 = t57 * t280 + t282 * t56;
t14 = t127 * t140 + t128 * t142 + t138 * t291 + t209 * t72 + t210 * t74 + t409 * t70;
t15 = t127 * t141 + t128 * t143 + t139 * t291 + t209 * t71 + t210 * t73 + t409 * t69;
t328 = t280 * t56 - t282 * t57;
t9 = -qJD(1) * t328 + t14 * t282 + t15 * t280;
t91 = t280 * t177 - t476;
t92 = -t280 * t470 - t282 * t306;
t336 = t40 * t380 + (t91 * t380 + t10 + (t121 * t282 + (t90 + t476) * qJD(1)) * t282) * t282 + (t9 + t92 * t380 + (t280 * t120 + (-t91 + t296) * qJD(1)) * t280 + ((-t411 * t466 + t413 * t468 + t121) * t280 + (t179 * t413 - t181 * t411 + t120 + t384) * t282 + ((t122 + t346) * t280 + (-t123 + t345) * t282) * t257 + ((t124 + t344) * t280 + (-t125 + t343) * t282) * t256 + (t92 - t89 + t483 + (-t177 + t306) * t282) * qJD(1)) * t282) * t280;
t335 = rSges(4,1) * t276 + rSges(4,2) * cos(pkin(10));
t327 = t62 * t280 + t282 * t61;
t326 = t280 * t61 - t282 * t62;
t325 = Icges(5,1) * t261 - t437;
t321 = -Icges(5,2) * t260 + t436;
t317 = Icges(5,5) * t261 - Icges(5,6) * t260;
t131 = t280 * t455 + t361 + t475;
t301 = -t271 * t282 + t223 + t386;
t132 = t301 + t184;
t313 = t131 * t280 - t132 * t282;
t310 = t145 * t282 + t146 * t280;
t304 = -t260 * t465 - t261 * t467;
t303 = t218 * t257 + t219 * t256;
t300 = t340 - t389;
t299 = t357 + t362;
t298 = rSges(3,3) * t282 + t280 * t457;
t66 = t163 + t222 * t381 + (-t113 - t199) * t282;
t295 = t304 * t280;
t294 = qJD(4) * t325;
t293 = qJD(4) * t321;
t290 = qJD(1) * t303 - t315 * t273;
t25 = t79 * t256 - t257 * t329;
t26 = t80 * t256 - t257 * t328;
t31 = t110 * t409 + t209 * t111 + t210 * t112 + t127 * t169 + t128 * t170 + t168 * t291;
t3 = (t273 * t328 + t31) * t256 + (-qJD(1) * t40 - t14 * t280 + t15 * t282 + t273 * t80) * t257;
t32 = -t110 * t410 + t207 * t111 + t208 * t112 + t129 * t169 + t130 * t170 + t168 * t292;
t4 = (t273 * t329 + t32) * t256 + (-qJD(1) * t39 - t16 * t280 + t17 * t282 + t273 * t79) * t257;
t288 = t3 * t459 + t4 * t458 + t256 * (-qJD(1) * t326 + t443 + t444) / 0.2e1 + t26 * t349 - t40 * t367 / 0.2e1 + t327 * t411 / 0.2e1 - t10 * t410 / 0.2e1 + t9 * t409 / 0.2e1 + (t368 / 0.2e1 - t354 / 0.2e1) * t39 + (t257 * t40 + t25) * t350;
t287 = t381 * t479 + t336;
t286 = t280 * t374 + t282 * t335;
t85 = (qJD(1) * t455 - t372) * t280 + t299 + t360;
t86 = t478 + (t455 * t282 + (-t333 + t398) * t280) * qJD(1) + t300;
t285 = m(6) * (t280 * t86 - t282 * t85 + (t131 * t282 + t132 * t280) * qJD(1));
t107 = t220 * t381 + t390 + t414;
t244 = pkin(4) * t353;
t108 = t220 * t380 + t244 + (-pkin(4) * t379 - t198) * t280;
t173 = (-t220 - t451) * t282;
t284 = m(6) * (-t107 * t282 + t108 * t280 + (t172 * t282 + t173 * t280) * qJD(1));
t283 = t444 / 0.2e1 + t443 / 0.2e1 + (t256 * t346 - t257 * t344 + t290 * t280 + t282 * t463 + t31) * t459 + (-t256 * t345 + t257 * t343 - t280 * t463 + t290 * t282 + t32) * t458 + (-t179 * t256 + t181 * t257 + t217 * t282 + t280 * t303 + t61 + t79) * t350 + (t256 * t468 - t257 * t466 - t282 * t303 + t402 + t62 + t80) * t349;
t212 = -rSges(3,2) * t282 + t280 * rSges(3,3) + t386;
t211 = t267 + t298;
t197 = t442 - t297;
t185 = t280 * rSges(6,3) - t475;
t176 = t265 + (t457 * t282 + (-rSges(3,3) - qJ(2)) * t280) * qJD(1);
t175 = qJD(1) * t298 + t387;
t174 = t282 * t185;
t166 = t280 * t335 + t282 * t440 + t386;
t165 = t267 + t286;
t161 = t388 + t391;
t160 = t282 * t161;
t151 = qJD(1) * t469 + t317 * t377;
t150 = -t317 * t376 + t383;
t149 = (t374 * t282 + (-qJ(2) - t335) * t280) * qJD(1) + t340;
t148 = qJD(1) * t286 + t357;
t144 = -qJD(1) * t388 + t362;
t137 = t393 * t282;
t133 = t282 * (t253 + (t235 - t452) * t381 + t389);
t126 = (-rSges(6,3) * qJD(1) - t372) * t280 + t360;
t119 = -t184 * t280 + t174;
t116 = t282 * (-t478 + (t280 * t333 + t268) * qJD(1));
t115 = (t393 - t451) * t282;
t114 = t250 + t136;
t98 = -t280 * t469 - t282 * t304;
t97 = t280 * t187 - t477;
t96 = -t282 * t469 + t295;
t95 = t187 * t282 + t280 * t305;
t94 = -t256 * t146 + t171 * t409;
t93 = t145 * t256 + t167 * t257;
t88 = -t280 * t356 + t247 + t301 + t392;
t87 = t282 * t481 + t331 + t361 - t448;
t83 = t310 * t257;
t81 = t280 * t394 + t160 + t174;
t65 = t280 * t396 + t397;
t64 = -pkin(4) * t352 + t244 + t67;
t63 = t66 + t390;
t60 = -t280 * t126 + t116 + (-t184 * t282 - t185 * t280) * qJD(1);
t47 = t280 * t363 + t160 + t397;
t46 = (t256 * t454 + t449) * t405 + (-t270 + (t398 - t481) * t280) * qJD(1) + t300 - t332;
t45 = (-t282 * t356 - t448) * qJD(1) + t299 + t359 + t364;
t44 = (-t167 * t273 + t76) * t256 + (t145 * t273 + t474) * t257;
t43 = (-t171 * t405 - t75) * t256 + (t113 * t282 - t146 * t273 - t163) * t257;
t42 = t116 + t133 + (-t126 - t144) * t280 + (t394 * t282 + (-t161 - t185) * t280) * qJD(1);
t28 = t310 * t413 + (-t280 * t75 - t282 * t76 + (t145 * t280 - t134) * qJD(1)) * t257;
t27 = t438 * t280 + (t280 * t395 + t282 * t396) * qJD(1) + t439;
t22 = t133 + (-t144 + t438) * t280 + (t363 * t282 + (-t161 + t395) * t280) * qJD(1) + t439;
t1 = [t320 * t379 + t289 - t324 * t378 + (t45 * t88 + t46 * t87) * t460 + (t131 * t86 + t132 * t85) * t461 + (t100 * t156 + t157 * t99) * t462 + 0.2e1 * m(3) * (t175 * t212 + t176 * t211) + 0.2e1 * m(4) * (t148 * t166 + t149 * t165) + t256 * t194 - t260 * t294 - t261 * t293 - t219 * t413 - t218 * t411 + (-t195 + t480) * t257; m(7) * (t280 * t46 - t282 * t45 + (t280 * t88 + t282 * t87) * qJD(1)) + t285 + m(5) * ((t156 * t282 + t157 * t280) * qJD(1) + t473) + m(3) * (-t175 * t282 + t280 * t176 + (t211 * t282 + t212 * t280) * qJD(1)) + m(4) * (-t148 * t282 + t280 * t149 + (t165 * t282 + t166 * t280) * qJD(1)); 0; m(7) * (t280 * t45 + t282 * t46 + (-t280 * t87 + t282 * t88) * qJD(1)) + m(6) * (-qJD(1) * t313 + t280 * t85 + t282 * t86) + m(5) * (qJD(1) * t472 + t100 * t282 + t280 * t99) + m(4) * (t280 * t148 + t149 * t282 + (-t165 * t280 + t166 * t282) * qJD(1)); 0; 0; m(5) * (t216 * t472 + t228 * t473) - (t275 / 0.2e1 + t274 / 0.2e1) * t316 * qJD(4) + (-qJD(4) * t305 - (qJD(1) * t467 + t280 * t293) * t260 + (qJD(1) * t465 + t280 * t294) * t261) * t458 + (-qJD(4) * t304 - (qJD(1) * t189 - t321 * t376) * t260 + (qJD(1) * t191 - t325 * t376) * t261) * t459 + m(7) * (t114 * t46 + t115 * t45 + t63 * t88 + t64 * t87) + m(6) * (t107 * t132 + t108 * t131 + t172 * t86 + t173 * t85) + ((t157 * t453 + t418 / 0.2e1 - t416 / 0.2e1) * t280 + (t156 * t453 + t417 / 0.2e1 - t415 / 0.2e1) * t282) * qJD(1) + t283; t284 + m(7) * (t64 * t280 - t282 * t63 + (t114 * t282 + t115 * t280) * qJD(1)) - m(5) * t347; m(6) * (t107 * t280 + t108 * t282 + (t173 * t282 - t419) * qJD(1)) + m(7) * (t63 * t280 + t282 * t64 + (-t114 * t280 + t115 * t282) * qJD(1)); (t114 * t64 + t115 * t63 + t22 * t47) * t460 + (t107 * t173 + t108 * t172 + t42 * t81) * t461 + ((-t280 * t196 + t197 * t282) * (-t280 * t358 + (-t228 * t275 + t274 * t445) * qJD(4) + ((-t196 + t269) * t282 + (-t197 + t297 + t442) * t280) * qJD(1)) - t228 * t347) * t462 + t282 * ((t151 * t282 + (t96 + t477) * qJD(1)) * t282 + (-t95 * qJD(1) + (-t378 * t465 + t379 * t467) * t280 + (t150 + (t416 - t418) * qJD(4) + (-t187 + t304) * qJD(1)) * t282) * t280) + (t98 * t280 + t282 * t97) * t380 + t280 * ((t280 * t150 + (-t97 + t295) * qJD(1)) * t280 + (t98 * qJD(1) + (t189 * t379 - t191 * t378 + t383) * t282 + (t151 + (t415 - t417) * qJD(4) + t305 * qJD(1)) * t280) * t282) + t336 + (-t96 * t280 - t282 * t95 + t479) * t381; m(7) * (t136 * t46 + t137 * t45 + t66 * t88 + t67 * t87) + t220 * t285 - m(6) * t313 * t198 + t283; m(7) * (t67 * t280 - t66 * t282 + (t136 * t282 + t137 * t280) * qJD(1)) - m(6) * t348; m(7) * (t66 * t280 + t282 * t67 + (-t136 * t280 + t137 * t282) * qJD(1)); m(7) * (t114 * t67 + t115 * t66 + t136 * t64 + t137 * t63 + t22 * t65 + t27 * t47) + m(6) * (t119 * t42 + t173 * t414 - t198 * t419 + t60 * t81) + t220 * t284 + t287; (t119 * t60 - t220 * t348) * t461 + (t136 * t67 + t137 * t66 + t27 * t65) * t460 + t287; m(7) * (t43 * t87 + t44 * t88 + t45 * t93 + t46 * t94) + (t280 * t369 + t282 * t370) * t413 + ((t21 / 0.2e1 + t31 / 0.2e1) * t282 + (-t32 / 0.2e1 - t20 / 0.2e1) * t280 + (t280 * t370 - t282 * t369) * qJD(1)) * t257 + t447; m(7) * (t43 * t280 - t282 * t44 + (t280 * t93 + t282 * t94) * qJD(1)); m(7) * (t44 * t280 + t282 * t43 + (-t280 * t94 + t282 * t93) * qJD(1)); t288 + m(7) * (t114 * t43 + t115 * t44 - t22 * t83 + t28 * t47 + t63 * t93 + t64 * t94); t288 + m(7) * (t136 * t43 + t137 * t44 - t27 * t83 + t28 * t65 + t66 * t93 + t67 * t94); (-t28 * t83 + t43 * t94 + t44 * t93) * t460 + ((t280 * t25 + t256 * t326 - t282 * t26) * t273 + t447) * t256 + (-t280 * t4 + t282 * t3 - t326 * t411 + (-t20 * t280 + t21 * t282 + t273 * t84) * t256 + (-t282 * t25 - t256 * t327 - t280 * t26) * qJD(1)) * t257;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

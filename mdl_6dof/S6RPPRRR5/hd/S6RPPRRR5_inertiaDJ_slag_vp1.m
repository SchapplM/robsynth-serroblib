% Calculate time derivative of joint inertia matrix for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:16
% EndTime: 2019-03-09 02:28:31
% DurationCPUTime: 9.38s
% Computational Cost: add. (17301->679), mult. (24536->964), div. (0->0), fcn. (22851->8), ass. (0->353)
t262 = sin(qJ(1));
t265 = cos(qJ(1));
t259 = qJ(4) + qJ(5);
t250 = sin(t259);
t263 = cos(qJ(6));
t385 = t263 * t265;
t260 = sin(qJ(6));
t389 = t260 * t262;
t206 = -t250 * t389 + t385;
t386 = t262 * t263;
t388 = t260 * t265;
t207 = t250 * t386 + t388;
t251 = cos(t259);
t394 = t251 * t262;
t132 = Icges(7,5) * t207 + Icges(7,6) * t206 - Icges(7,3) * t394;
t134 = Icges(7,4) * t207 + Icges(7,2) * t206 - Icges(7,6) * t394;
t136 = Icges(7,1) * t207 + Icges(7,4) * t206 - Icges(7,5) * t394;
t208 = -t250 * t388 - t386;
t209 = t250 * t385 - t389;
t393 = t251 * t265;
t49 = -t132 * t393 + t134 * t208 + t136 * t209;
t133 = Icges(7,5) * t209 + Icges(7,6) * t208 - Icges(7,3) * t393;
t135 = Icges(7,4) * t209 + Icges(7,2) * t208 - Icges(7,6) * t393;
t137 = Icges(7,1) * t209 + Icges(7,4) * t208 - Icges(7,5) * t393;
t50 = -t133 * t393 + t135 * t208 + t137 * t209;
t322 = t262 * t50 - t265 * t49;
t309 = Icges(6,5) * t250 + Icges(6,6) * t251;
t170 = Icges(6,3) * t265 + t262 * t309;
t419 = Icges(6,4) * t250;
t312 = Icges(6,2) * t251 + t419;
t172 = Icges(6,6) * t265 + t262 * t312;
t418 = Icges(6,4) * t251;
t315 = Icges(6,1) * t250 + t418;
t174 = Icges(6,5) * t265 + t262 * t315;
t301 = t172 * t251 + t174 * t250;
t285 = t301 * t265;
t88 = -t170 * t262 + t285;
t171 = -Icges(6,3) * t262 + t265 * t309;
t173 = -Icges(6,6) * t262 + t265 * t312;
t175 = -Icges(6,5) * t262 + t265 * t315;
t300 = t173 * t251 + t175 * t250;
t89 = -t171 * t262 + t265 * t300;
t461 = t262 * t89 - t265 * t88 + t322;
t47 = -t132 * t394 + t134 * t206 + t136 * t207;
t48 = -t133 * t394 + t135 * t206 + t137 * t207;
t324 = t262 * t48 - t265 * t47;
t86 = t170 * t265 + t262 * t301;
t284 = t300 * t262;
t87 = t171 * t265 + t284;
t458 = t262 * t87 - t265 * t86 + t324;
t376 = t209 * rSges(7,1) + t208 * rSges(7,2);
t139 = -rSges(7,3) * t393 + t376;
t396 = t250 * t265;
t233 = pkin(5) * t396;
t460 = -pkin(9) * t393 + t139 + t233;
t366 = qJD(1) * t262;
t347 = t251 * t366;
t256 = qJD(4) + qJD(5);
t390 = t256 * t265;
t355 = t250 * t390;
t276 = t347 + t355;
t416 = Icges(7,4) * t263;
t311 = -Icges(7,2) * t260 + t416;
t397 = t250 * t256;
t417 = Icges(7,4) * t260;
t109 = -t311 * t397 + (Icges(7,6) * t256 + (-Icges(7,2) * t263 - t417) * qJD(6)) * t251;
t165 = Icges(7,6) * t250 + t251 * t311;
t314 = Icges(7,1) * t263 - t417;
t166 = Icges(7,5) * t250 + t251 * t314;
t459 = -t109 * t260 + (-t165 * t263 - t166 * t260) * qJD(6);
t353 = t251 * t390;
t334 = -qJD(6) * t250 - qJD(1);
t295 = t263 * t334;
t333 = qJD(1) * t250 + qJD(6);
t450 = t262 * t333 - t353;
t114 = t260 * t450 + t265 * t295;
t294 = t334 * t260;
t115 = -t263 * t450 + t265 * t294;
t70 = t115 * rSges(7,1) + t114 * rSges(7,2) + rSges(7,3) * t276;
t457 = -pkin(5) * t353 - pkin(9) * t276 - t70;
t325 = rSges(7,1) * t263 - rSges(7,2) * t260;
t111 = -t325 * t397 + (rSges(7,3) * t256 + (-rSges(7,1) * t260 - rSges(7,2) * t263) * qJD(6)) * t251;
t167 = rSges(7,3) * t250 + t251 * t325;
t365 = qJD(1) * t265;
t456 = -t262 * t111 - t167 * t365;
t257 = t262 ^ 2;
t258 = t265 ^ 2;
t369 = t257 + t258;
t254 = t265 * qJ(2);
t261 = sin(qJ(4));
t264 = cos(qJ(4));
t329 = rSges(5,1) * t261 + rSges(5,2) * t264;
t432 = -pkin(1) - qJ(3);
t289 = -t329 + t432;
t441 = -rSges(5,3) - pkin(7);
t271 = t262 * t289 + t265 * t441;
t153 = t254 + t271;
t370 = t265 * pkin(1) + t262 * qJ(2);
t349 = t265 * qJ(3) + t370;
t384 = t264 * t265;
t387 = t261 * t265;
t374 = rSges(5,1) * t387 + rSges(5,2) * t384;
t154 = t262 * t441 + t349 + t374;
t455 = t153 * t265 + t154 * t262;
t361 = qJD(4) * t265;
t345 = t264 * t361;
t234 = rSges(5,1) * t345;
t346 = t261 * t361;
t372 = qJ(2) * t365 + qJD(2) * t262;
t350 = qJD(3) * t265 + t372;
t96 = -rSges(5,2) * t346 + qJD(1) * t271 + t234 + t350;
t430 = rSges(5,2) * t261;
t230 = rSges(5,1) * t264 - t430;
t246 = pkin(7) * t366;
t249 = qJD(2) * t265;
t97 = t246 + t249 + (-t230 * qJD(4) - qJD(3)) * t262 + ((rSges(5,3) - qJ(2)) * t262 + t289 * t265) * qJD(1);
t454 = t262 * t96 + t265 * t97;
t453 = qJD(1) * t170;
t328 = rSges(6,1) * t250 + rSges(6,2) * t251;
t428 = rSges(6,3) * t265;
t178 = t262 * t328 + t428;
t452 = qJD(1) * t178;
t310 = Icges(5,5) * t261 + Icges(5,6) * t264;
t189 = Icges(5,3) * t265 + t262 * t310;
t451 = qJD(1) * t189;
t421 = Icges(5,4) * t261;
t313 = Icges(5,2) * t264 + t421;
t191 = Icges(5,6) * t265 + t262 * t313;
t420 = Icges(5,4) * t264;
t316 = Icges(5,1) * t261 + t420;
t193 = Icges(5,5) * t265 + t262 * t316;
t213 = -Icges(6,2) * t250 + t418;
t214 = Icges(6,1) * t251 - t419;
t296 = t213 * t251 + t214 * t250;
t448 = qJD(1) * t296 - t309 * t256;
t447 = 2 * m(5);
t446 = 2 * m(6);
t445 = 2 * m(7);
t444 = -t262 / 0.2e1;
t443 = t265 / 0.2e1;
t442 = rSges(3,2) - pkin(1);
t440 = rSges(7,3) + pkin(9);
t212 = Icges(6,5) * t251 - Icges(6,6) * t250;
t286 = t256 * t212;
t119 = t265 * t286 - t453;
t368 = qJD(1) * t171;
t120 = t262 * t286 + t368;
t287 = t213 * t256;
t121 = -qJD(1) * t172 + t265 * t287;
t288 = t214 * t256;
t123 = -qJD(1) * t174 + t265 * t288;
t124 = qJD(1) * t175 + t262 * t288;
t335 = -t172 * t256 + t124;
t122 = qJD(1) * t173 + t262 * t287;
t337 = t174 * t256 + t122;
t395 = t251 * t256;
t392 = t256 * t262;
t354 = t251 * t392;
t116 = t262 * t295 + (-t265 * t333 - t354) * t260;
t391 = t256 * t263;
t117 = t333 * t385 + (t251 * t391 + t294) * t262;
t277 = t250 * t392 - t251 * t365;
t65 = Icges(7,5) * t117 + Icges(7,6) * t116 + Icges(7,3) * t277;
t67 = Icges(7,4) * t117 + Icges(7,2) * t116 + Icges(7,6) * t277;
t69 = Icges(7,1) * t117 + Icges(7,4) * t116 + Icges(7,5) * t277;
t15 = t116 * t134 + t117 * t136 + t132 * t277 + t206 * t67 + t207 * t69 - t394 * t65;
t64 = Icges(7,5) * t115 + Icges(7,6) * t114 + Icges(7,3) * t276;
t66 = Icges(7,4) * t115 + Icges(7,2) * t114 + Icges(7,6) * t276;
t68 = Icges(7,1) * t115 + Icges(7,4) * t114 + Icges(7,5) * t276;
t16 = t116 * t135 + t117 * t137 + t133 * t277 + t206 * t66 + t207 * t68 - t394 * t64;
t323 = t262 * t47 + t265 * t48;
t9 = -qJD(1) * t323 + t15 * t265 - t16 * t262;
t439 = (t9 + (-t86 * qJD(1) + (-t121 * t251 - t123 * t250 + t173 * t397 - t175 * t395 + t368) * t262) * t262 + ((-t87 + t285) * qJD(1) + (-t119 + t337 * t251 + t335 * t250 + (-t170 - t300) * qJD(1)) * t262 + t265 * t120) * t265) * t265;
t336 = t173 * t256 - t123;
t338 = -t175 * t256 - t121;
t12 = (t119 * t262 + (-t88 + t284) * qJD(1)) * t262 + (-t89 * qJD(1) + (t122 * t251 + t124 * t250 - t172 * t397 + t174 * t395 - t453) * t265 + (-t120 + t338 * t251 + t336 * t250 + (t171 - t301) * qJD(1)) * t262) * t265;
t13 = t114 * t134 + t115 * t136 + t132 * t276 + t208 * t67 + t209 * t69 - t393 * t65;
t14 = t114 * t135 + t115 * t137 + t133 * t276 + t208 * t66 + t209 * t68 - t393 * t64;
t321 = t262 * t49 + t265 * t50;
t8 = -qJD(1) * t321 + t13 * t265 - t14 * t262;
t438 = -t12 - t8;
t437 = m(5) * t230;
t436 = pkin(4) * t261;
t435 = pkin(4) * t264;
t434 = pkin(5) * t250;
t433 = pkin(5) * t251;
t308 = Icges(7,5) * t263 - Icges(7,6) * t260;
t108 = -t308 * t397 + (Icges(7,3) * t256 + (-Icges(7,5) * t260 - Icges(7,6) * t263) * qJD(6)) * t251;
t110 = -t314 * t397 + (Icges(7,5) * t256 + (-Icges(7,1) * t260 - t416) * qJD(6)) * t251;
t164 = Icges(7,3) * t250 + t251 * t308;
t406 = t165 * t260;
t275 = t251 * t263 * t110 + t164 * t395 + t397 * t406 + (-t166 * t391 + t108) * t250;
t82 = t164 * t250 + (t166 * t263 - t406) * t251;
t431 = (t251 * t459 + t275) * t250 + t82 * t395;
t429 = rSges(6,3) * t262;
t306 = t134 * t260 - t136 * t263;
t19 = (t256 * t306 + t65) * t250 + (t132 * t256 - t260 * t67 + t263 * t69 + (-t134 * t263 - t136 * t260) * qJD(6)) * t251;
t427 = t19 * t265;
t305 = t135 * t260 - t137 * t263;
t20 = (t256 * t305 + t64) * t250 + (t133 * t256 - t260 * t66 + t263 * t68 + (-t135 * t263 - t137 * t260) * qJD(6)) * t251;
t426 = t20 * t262;
t423 = t366 * t434 + t457;
t327 = t117 * rSges(7,1) + t116 * rSges(7,2);
t71 = rSges(7,3) * t277 + t327;
t422 = -t277 * pkin(9) - (t250 * t365 + t354) * pkin(5) - t71;
t187 = t328 * t256;
t403 = t187 * t265;
t402 = t191 * t261;
t192 = -Icges(5,6) * t262 + t265 * t313;
t401 = t192 * t261;
t400 = t193 * t264;
t194 = -Icges(5,5) * t262 + t265 * t316;
t399 = t194 * t264;
t215 = rSges(6,1) * t251 - rSges(6,2) * t250;
t398 = t215 * t265;
t161 = t262 * t167;
t105 = t265 * t111;
t162 = t265 * t167;
t266 = -pkin(8) - pkin(7);
t383 = -qJ(2) - t266;
t330 = -pkin(9) * t251 + t434;
t188 = t330 * t256;
t382 = -t265 * t188 + t105;
t381 = t460 * t366;
t326 = -rSges(7,1) * t207 - rSges(7,2) * t206;
t138 = -rSges(7,3) * t394 - t326;
t380 = -t330 * t262 - t138;
t216 = pkin(9) * t250 + t433;
t378 = -t167 - t216;
t245 = t265 * t266;
t292 = pkin(7) * t265 - t262 * t436;
t205 = -t245 - t292;
t377 = -t178 - t205;
t142 = t262 * t216 + t161;
t143 = t265 * t216 + t162;
t364 = qJD(1) * t266;
t375 = pkin(4) * t345 + t265 * t364;
t373 = pkin(4) * t387 + t262 * t266;
t371 = t245 + t254;
t190 = -Icges(5,3) * t262 + t265 * t310;
t367 = qJD(1) * t190;
t363 = qJD(4) * t261;
t362 = qJD(4) * t264;
t359 = -rSges(4,3) + t432;
t242 = pkin(4) * t384;
t358 = pkin(4) * t363;
t59 = t132 * t250 - t251 * t306;
t76 = -t164 * t394 + t165 * t206 + t166 * t207;
t357 = t76 / 0.2e1 + t59 / 0.2e1;
t60 = t133 * t250 - t251 * t305;
t77 = -t164 * t393 + t165 * t208 + t166 * t209;
t356 = t77 / 0.2e1 + t60 / 0.2e1;
t352 = -t205 + t380;
t348 = t215 * t365;
t343 = t461 * t265;
t342 = -t366 / 0.2e1;
t341 = -t365 / 0.2e1;
t340 = t187 * t369;
t222 = t329 * qJD(4);
t339 = t222 * t369;
t73 = -t262 * t188 + t216 * t365 - t456;
t332 = t432 - t436;
t179 = rSges(6,1) * t396 + rSges(6,2) * t393 - t429;
t331 = -pkin(4) * t362 - qJD(3);
t320 = t262 * t60 - t265 * t59;
t319 = t262 * t59 + t265 * t60;
t318 = t350 + t375;
t317 = t349 + t373;
t274 = -t328 + t332;
t270 = t262 * t274 - t428;
t130 = t270 + t371;
t131 = t179 + t317;
t307 = t130 * t265 + t131 * t262;
t304 = t138 * t265 - t139 * t262;
t299 = t191 * t264 + t193 * t261;
t298 = t192 * t264 + t194 * t261;
t293 = rSges(6,1) * t353 - rSges(6,2) * t355;
t291 = rSges(3,3) * t265 + t262 * t442;
t290 = t332 - t434;
t283 = t299 * t265;
t282 = t298 * t262;
t281 = qJD(4) * (Icges(5,1) * t264 - t421);
t280 = qJD(4) * (-Icges(5,2) * t261 + t420);
t279 = qJD(4) * (Icges(5,5) * t264 - Icges(5,6) * t261);
t278 = rSges(4,2) * t265 + t262 * t359;
t25 = t250 * t76 - t251 * t323;
t26 = t250 * t77 - t251 * t321;
t30 = -t108 * t393 + t109 * t208 + t110 * t209 + t114 * t165 + t115 * t166 + t164 * t276;
t3 = (t256 * t321 + t30) * t250 + (qJD(1) * t322 - t13 * t262 - t14 * t265 + t256 * t77) * t251;
t31 = -t108 * t394 + t109 * t206 + t110 * t207 + t116 * t165 + t117 * t166 + t164 * t277;
t4 = (t256 * t323 + t31) * t250 + (qJD(1) * t324 - t15 * t262 - t16 * t265 + t256 * t76) * t251;
t273 = t3 * t444 + t4 * t443 + t250 * (-qJD(1) * t319 - t426 + t427) / 0.2e1 + t25 * t342 - t322 * t347 / 0.2e1 - t320 * t395 / 0.2e1 - t9 * t394 / 0.2e1 - t8 * t393 / 0.2e1 + (-t262 * t324 - t265 * t322) * t397 / 0.2e1 + (-t251 * t324 + t26) * t341;
t272 = t251 * t440 + t290;
t80 = qJD(1) * t270 + t293 + t318;
t81 = t249 + (-t215 * t256 + t331) * t262 + ((rSges(6,3) + t383) * t262 + t274 * t265) * qJD(1);
t269 = m(6) * (t262 * t80 + t265 * t81 + (-t130 * t262 + t131 * t265) * qJD(1));
t185 = t312 * t256;
t186 = t315 * t256;
t267 = -qJD(1) * t212 + (-t185 + t288) * t251 + (-t186 - t287) * t250;
t268 = t427 / 0.2e1 - t426 / 0.2e1 + (t250 * t338 - t251 * t336 - t262 * t448 + t267 * t265 + t30) * t444 + (-t250 * t337 + t251 * t335 + t267 * t262 + t265 * t448 + t31) * t443 + (-t172 * t250 + t174 * t251 + t212 * t265 + t262 * t296 + t59 + t76) * t342 + (-t173 * t250 + t175 * t251 - t212 * t262 + t265 * t296 + t60 + t77) * t341;
t240 = t262 * t435;
t236 = qJD(1) * t242;
t204 = pkin(7) * t262 + t373;
t200 = -rSges(3,2) * t265 + rSges(3,3) * t262 + t370;
t199 = t254 + t291;
t196 = -rSges(5,3) * t262 + t374;
t195 = t265 * rSges(5,3) + t262 * t329;
t183 = t204 * t366;
t182 = rSges(4,2) * t262 + rSges(4,3) * t265 + t349;
t181 = t254 + t278;
t169 = t242 + t398;
t168 = t215 * t262 + t240;
t163 = t179 * t366;
t160 = t249 + (t442 * t265 + (-rSges(3,3) - qJ(2)) * t262) * qJD(1);
t159 = qJD(1) * t291 + t372;
t157 = qJD(1) * t292 + t375;
t156 = t262 * t364 + t246 + (t261 * t365 + t262 * t362) * pkin(4);
t151 = -qJD(3) * t262 + t249 + ((-rSges(4,2) - qJ(2)) * t262 + t359 * t265) * qJD(1);
t150 = qJD(1) * t278 + t350;
t145 = t262 * t279 + t367;
t144 = t265 * t279 - t451;
t129 = t242 + t143;
t128 = t240 + t142;
t127 = t215 * t392 + (t265 * t328 - t429) * qJD(1);
t126 = t293 - t452;
t118 = -t178 * t262 - t179 * t265;
t107 = t348 + t236 + (-t187 - t358) * t262;
t106 = -t215 * t366 - t403 + (-t264 * t366 - t346) * pkin(4);
t95 = -t190 * t262 + t265 * t298;
t94 = -t189 * t262 + t283;
t93 = t190 * t265 + t282;
t92 = t189 * t265 + t262 * t299;
t91 = t139 * t250 + t162 * t251;
t90 = -t138 * t250 - t161 * t251;
t85 = -t393 * t440 + t233 + t317 + t376;
t84 = t262 * t272 + t326 + t371;
t83 = (-t179 - t204) * t265 + t377 * t262;
t78 = t304 * t251;
t72 = t366 * t378 + t382;
t63 = t262 * t380 - t265 * t460;
t62 = -t262 * t358 + t236 + t73;
t61 = -pkin(4) * t346 + (t378 - t435) * t366 + t382;
t58 = -t127 * t262 + t163 + (-t126 - t452) * t265;
t53 = (-t204 - t460) * t265 + t352 * t262;
t44 = t249 + ((-t250 * t440 - t433) * t256 + t331) * t262 + (t262 * t383 + t265 * t272) * qJD(1) - t327;
t43 = t290 * t366 + t318 - t457;
t42 = t163 + t183 + (-t127 - t156) * t262 + (qJD(1) * t377 - t126 - t157) * t265;
t41 = (t161 * t256 - t71) * t250 + (-t138 * t256 + t456) * t251;
t40 = (-t162 * t256 + t70) * t250 + (t139 * t256 - t167 * t366 + t105) * t251;
t27 = t304 * t397 + (t262 * t70 - t265 * t71 + (t138 * t262 + t139 * t265) * qJD(1)) * t251;
t24 = t422 * t262 + (qJD(1) * t380 + t423) * t265 + t381;
t21 = t183 + (-t156 + t422) * t262 + (qJD(1) * t352 - t157 + t423) * t265 + t381;
t1 = [-t316 * t362 + t313 * t363 - t261 * t281 - t264 * t280 + (t153 * t97 + t154 * t96) * t447 + 0.2e1 * m(4) * (t150 * t182 + t151 * t181) + 0.2e1 * m(3) * (t159 * t200 + t160 * t199) + t275 + (t43 * t85 + t44 * t84) * t445 + (t130 * t81 + t131 * t80) * t446 + t250 * t185 - t214 * t397 - t213 * t395 + (-t186 + t459) * t251; m(7) * (t262 * t44 - t265 * t43 + (t262 * t85 + t265 * t84) * qJD(1)) + m(6) * (qJD(1) * t307 + t262 * t81 - t265 * t80) + m(5) * (qJD(1) * t455 + t262 * t97 - t265 * t96) + m(4) * (-t150 * t265 + t151 * t262 + (t181 * t265 + t182 * t262) * qJD(1)) + m(3) * (-t159 * t265 + t160 * t262 + (t199 * t265 + t200 * t262) * qJD(1)); 0; m(7) * (t262 * t43 + t265 * t44 + (-t262 * t84 + t265 * t85) * qJD(1)) + t269 + m(5) * ((-t153 * t262 + t154 * t265) * qJD(1) + t454) + m(4) * (t150 * t262 + t151 * t265 + (-t181 * t262 + t182 * t265) * qJD(1)); 0; 0; -(t258 / 0.2e1 + t257 / 0.2e1) * t310 * qJD(4) + m(7) * (t128 * t43 + t129 * t44 + t61 * t84 + t62 * t85) + m(6) * (t106 * t130 + t107 * t131 + t168 * t80 + t169 * t81) + ((t401 / 0.2e1 - t399 / 0.2e1 + t154 * t437) * t265 + (-t153 * t437 + t402 / 0.2e1 - t400 / 0.2e1) * t262) * qJD(1) + (-qJD(4) * t299 - (qJD(1) * t192 + t262 * t280) * t261 + (qJD(1) * t194 + t262 * t281) * t264) * t443 + (-qJD(4) * t298 - (-qJD(1) * t191 + t265 * t280) * t261 + (-qJD(1) * t193 + t265 * t281) * t264) * t444 + m(5) * (-t222 * t455 + t230 * t454) + t268; m(6) * (t106 * t262 - t107 * t265 + (t168 * t262 + t169 * t265) * qJD(1)) + m(7) * (t262 * t61 - t265 * t62 + (t128 * t262 + t129 * t265) * qJD(1)); m(6) * (t106 * t265 + t107 * t262 + (t168 * t265 - t169 * t262) * qJD(1)) + m(7) * (t262 * t62 + t265 * t61 + (t128 * t265 - t129 * t262) * qJD(1)) - m(5) * t339; -t262 * t8 + (t128 * t62 + t129 * t61 + t53 * t21) * t445 - t262 * t12 + (t106 * t169 + t107 * t168 + t42 * t83) * t446 + (-t230 * t339 + (-t265 * t234 + (-t230 * t257 + t258 * t430) * qJD(4) + (rSges(5,3) * t369 - t265 * t195 + t262 * t196) * qJD(1)) * (-t195 * t262 - t196 * t265)) * t447 + t265 * ((t265 * t145 + (-t93 + t283) * qJD(1)) * t265 + (-t92 * qJD(1) + (t192 * t363 - t194 * t362 + t367) * t262 + (-t144 + (t400 - t402) * qJD(4) + (-t189 - t298) * qJD(1)) * t265) * t262) - t262 * ((t262 * t144 + (-t94 + t282) * qJD(1)) * t262 + (-t95 * qJD(1) + (-t191 * t363 + t193 * t362 - t451) * t265 + (-t145 + (-t399 + t401) * qJD(4) + (t190 - t299) * qJD(1)) * t262) * t265) + t439 + (t262 * t93 - t265 * t92 + t458) * t366 + (t262 * t95 - t265 * t94 + t461) * t365; m(7) * (t142 * t43 + t143 * t44 + t72 * t84 + t73 * t85) + t268 - m(6) * t307 * t187 + t215 * t269; m(7) * (t262 * t72 - t265 * t73 + (t142 * t262 + t143 * t265) * qJD(1)); m(7) * (t73 * t262 + t72 * t265 + (t142 * t265 - t143 * t262) * qJD(1)) - m(6) * t340; qJD(1) * t343 + m(7) * (t128 * t73 + t129 * t72 + t142 * t62 + t143 * t61 + t63 * t21 + t24 * t53) + m(6) * (t106 * t398 + t118 * t42 + t168 * t348 - t169 * t403 + t58 * t83) + (m(6) * (t107 * t215 - t168 * t187) + t438 + (-m(6) * t169 * t215 + t458) * qJD(1)) * t262 + t439; t438 * t262 + (t118 * t58 - t215 * t340) * t446 + (t142 * t73 + t143 * t72 + t63 * t24) * t445 + (t262 * t458 + t343) * qJD(1) + t439; m(7) * (t40 * t85 + t41 * t84 + t43 * t91 + t44 * t90) + (t262 * t357 + t265 * t356) * t397 + ((-t30 / 0.2e1 - t20 / 0.2e1) * t265 + (-t31 / 0.2e1 - t19 / 0.2e1) * t262 + (t262 * t356 - t265 * t357) * qJD(1)) * t251 + t431; m(7) * (t262 * t41 - t265 * t40 + (t262 * t91 + t265 * t90) * qJD(1)); m(7) * (t262 * t40 + t265 * t41 + (-t262 * t90 + t265 * t91) * qJD(1)); m(7) * (t128 * t40 + t129 * t41 - t21 * t78 + t27 * t53 + t61 * t90 + t62 * t91) + t273; m(7) * (t142 * t40 + t143 * t41 - t24 * t78 + t27 * t63 + t72 * t90 + t73 * t91) + t273; (-t27 * t78 + t40 * t91 + t41 * t90) * t445 + ((t262 * t25 + t250 * t319 + t265 * t26) * t256 + t431) * t250 + (-t265 * t3 - t262 * t4 - t319 * t395 + (-t19 * t262 - t20 * t265 + t256 * t82) * t250 + (-t265 * t25 + t250 * t320 + t262 * t26) * qJD(1)) * t251;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

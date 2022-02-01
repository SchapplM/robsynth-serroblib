% Calculate time derivative of joint inertia matrix for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:40
% EndTime: 2022-01-20 12:07:58
% DurationCPUTime: 6.64s
% Computational Cost: add. (19037->564), mult. (12740->778), div. (0->0), fcn. (9780->10), ass. (0->335)
t261 = qJ(1) + qJ(2);
t251 = sin(t261);
t253 = cos(t261);
t264 = cos(qJ(3));
t362 = qJD(3) * t264;
t258 = qJD(1) + qJD(2);
t262 = sin(qJ(3));
t379 = t258 * t262;
t444 = -t251 * t362 - t253 * t379;
t260 = qJ(3) + qJ(4);
t250 = sin(t260);
t382 = t253 * t258;
t252 = cos(t260);
t257 = qJD(3) + qJD(4);
t384 = t252 * t257;
t443 = -t250 * t382 - t251 * t384;
t411 = rSges(4,2) * t262;
t414 = rSges(4,1) * t264;
t442 = -t411 + t414;
t406 = Icges(4,4) * t264;
t313 = -Icges(4,2) * t262 + t406;
t291 = t313 * t253;
t163 = Icges(4,6) * t251 + t291;
t407 = Icges(4,4) * t262;
t316 = Icges(4,1) * t264 - t407;
t294 = t316 * t253;
t165 = Icges(4,5) * t251 + t294;
t300 = t163 * t262 - t165 * t264;
t441 = t251 * t300;
t404 = Icges(5,4) * t252;
t312 = -Icges(5,2) * t250 + t404;
t290 = t312 * t253;
t144 = Icges(5,6) * t251 + t290;
t405 = Icges(5,4) * t250;
t315 = Icges(5,1) * t252 - t405;
t293 = t315 * t253;
t146 = Icges(5,5) * t251 + t293;
t304 = t144 * t250 - t146 * t252;
t440 = t251 * t304;
t254 = qJ(5) + t260;
t244 = sin(t254);
t245 = cos(t254);
t402 = Icges(6,4) * t245;
t311 = -Icges(6,2) * t244 + t402;
t289 = t311 * t253;
t136 = Icges(6,6) * t251 + t289;
t403 = Icges(6,4) * t244;
t314 = Icges(6,1) * t245 - t403;
t292 = t314 * t253;
t138 = Icges(6,5) * t251 + t292;
t306 = t136 * t244 - t138 * t245;
t439 = t251 * t306;
t162 = -Icges(4,6) * t253 + t251 * t313;
t164 = -Icges(4,5) * t253 + t251 * t316;
t302 = t162 * t262 - t164 * t264;
t438 = t253 * t302;
t143 = -Icges(5,6) * t253 + t251 * t312;
t145 = -Icges(5,5) * t253 + t251 * t315;
t305 = t143 * t250 - t145 * t252;
t437 = t253 * t305;
t135 = -Icges(6,6) * t253 + t251 * t311;
t137 = -Icges(6,5) * t253 + t251 * t314;
t307 = t135 * t244 - t137 * t245;
t436 = t253 * t307;
t246 = t264 * pkin(3) + pkin(2);
t214 = pkin(4) * t252 + t246;
t369 = t214 - t246;
t435 = t369 * t251;
t234 = t251 * rSges(6,3);
t412 = rSges(6,1) * t245;
t434 = t253 * t412 + t234;
t235 = t251 * rSges(5,3);
t413 = rSges(5,1) * t252;
t433 = t253 * t413 + t235;
t249 = qJD(5) + t257;
t183 = Icges(6,2) * t245 + t403;
t184 = Icges(6,1) * t244 + t402;
t299 = t183 * t244 - t184 * t245;
t308 = Icges(6,5) * t245 - Icges(6,6) * t244;
t432 = t308 * t249 + t258 * t299;
t198 = Icges(5,2) * t252 + t405;
t199 = Icges(5,1) * t250 + t404;
t298 = t198 * t250 - t199 * t252;
t309 = Icges(5,5) * t252 - Icges(5,6) * t250;
t431 = t309 * t257 + t258 * t298;
t223 = Icges(4,2) * t264 + t407;
t224 = Icges(4,1) * t262 + t406;
t297 = t223 * t262 - t224 * t264;
t310 = Icges(4,5) * t264 - Icges(4,6) * t262;
t430 = t310 * qJD(3) + t258 * t297;
t266 = -pkin(8) - pkin(7);
t380 = t253 * t266;
t415 = pkin(2) - t246;
t429 = t251 * t415 - t380;
t428 = 2 * m(3);
t427 = 2 * m(4);
t426 = 2 * m(5);
t425 = 2 * m(6);
t424 = t251 / 0.2e1;
t423 = -t253 / 0.2e1;
t208 = t442 * qJD(3);
t422 = m(4) * t208;
t230 = rSges(4,1) * t262 + rSges(4,2) * t264;
t421 = m(4) * t230;
t410 = rSges(5,2) * t250;
t174 = (-t410 + t413) * t257;
t420 = m(5) * t174;
t409 = rSges(5,2) * t252;
t200 = rSges(5,1) * t250 + t409;
t419 = m(5) * t200;
t418 = pkin(3) * t262;
t241 = t251 * pkin(7);
t133 = -Icges(6,3) * t253 + t251 * t308;
t33 = -t133 * t253 - t251 * t307;
t395 = t183 * t249;
t279 = Icges(6,6) * t258 - t395;
t70 = t251 * t279 + t258 * t289;
t331 = t137 * t249 + t70;
t394 = t184 * t249;
t281 = Icges(6,5) * t258 - t394;
t72 = t251 * t281 + t258 * t292;
t333 = -t135 * t249 + t72;
t286 = t308 * t253;
t134 = Icges(6,3) * t251 + t286;
t34 = -t134 * t253 - t439;
t388 = t245 * t249;
t391 = t244 * t249;
t399 = t134 * t258;
t400 = t133 * t258;
t182 = Icges(6,5) * t244 + Icges(6,6) * t245;
t277 = Icges(6,3) * t258 - t182 * t249;
t386 = t251 * t258;
t67 = t253 * t277 - t308 * t386;
t68 = t251 * t277 + t258 * t286;
t69 = t253 * t279 - t311 * t386;
t71 = t253 * t281 - t314 * t386;
t2 = (t253 * t68 + (t34 + t436) * t258) * t253 + (t33 * t258 + (-t136 * t388 - t138 * t391 - t244 * t69 + t245 * t71 + t399) * t251 + (-t400 - t67 + (t138 * t258 - t333) * t245 + (-t136 * t258 + t331) * t244) * t253) * t251;
t417 = t253 * t2;
t263 = sin(qJ(1));
t416 = t263 * pkin(1);
t408 = pkin(1) * qJD(1);
t236 = t251 * rSges(4,3);
t188 = rSges(6,1) * t244 + rSges(6,2) * t245;
t336 = -pkin(4) * t250 - t188;
t295 = t336 - t418;
t115 = t295 * t253;
t401 = t115 * t258;
t141 = -Icges(5,3) * t253 + t251 * t309;
t398 = t141 * t258;
t287 = t309 * t253;
t142 = Icges(5,3) * t251 + t287;
t397 = t142 * t258;
t396 = t174 * t251;
t393 = t198 * t257;
t392 = t199 * t257;
t390 = t244 * t251;
t389 = t244 * t253;
t387 = t250 * t257;
t385 = t251 * t264;
t383 = t253 * t257;
t259 = pkin(9) - t266;
t381 = t253 * t259;
t378 = t258 * t266;
t233 = t259 * t251;
t324 = -t253 * t246 + t251 * t266;
t373 = t253 * t214 + t233;
t111 = t324 + t373;
t357 = rSges(6,2) * t389;
t140 = -t357 + t434;
t377 = -t111 - t140;
t242 = t253 * pkin(7);
t131 = t242 - t429;
t365 = -t253 * pkin(2) - t241;
t132 = -t324 + t365;
t376 = t251 * t131 + t253 * t132;
t210 = rSges(6,2) * t390;
t370 = t253 * rSges(6,3) + t210;
t139 = t251 * t412 - t370;
t75 = t251 * t139 + t253 * t140;
t367 = t253 * rSges(5,3) + t251 * t410;
t147 = t251 * t413 - t367;
t148 = -t253 * t410 + t433;
t86 = t251 * t147 + t253 * t148;
t351 = t250 * t386;
t375 = pkin(4) * t351 + t188 * t386;
t363 = qJD(3) * t262;
t354 = pkin(3) * t363;
t185 = -pkin(4) * t387 - t354;
t374 = t253 * t185 + t258 * t381;
t372 = rSges(5,2) * t351 + rSges(5,3) * t382;
t348 = t251 * t379;
t371 = rSges(4,2) * t348 + rSges(4,3) * t382;
t343 = t251 * t363;
t368 = pkin(3) * t343 + t251 * t378;
t366 = t253 * rSges(4,3) + t251 * t411;
t364 = t251 ^ 2 + t253 ^ 2;
t330 = t138 * t249 + t69;
t332 = -t136 * t249 + t71;
t35 = t133 * t251 - t436;
t36 = t134 * t251 - t253 * t306;
t361 = t251 * ((t251 * t67 + (t35 + t439) * t258) * t251 + (t36 * t258 + (t135 * t388 + t137 * t391 + t244 * t70 - t245 * t72 - t400) * t253 + (t399 - t68 + (t137 * t258 + t332) * t245 + (-t135 * t258 - t330) * t244) * t251) * t253) + (t251 * t34 - t253 * t33) * t386 + (t251 * t36 - t253 * t35) * t382;
t229 = pkin(7) * t382;
t341 = t253 * t363;
t323 = pkin(3) * t341;
t360 = t131 * t382 + t251 * ((-t253 * t415 - t241) * t258 - t368) + t253 * (t429 * t258 - t229 - t323);
t356 = rSges(6,2) * t388;
t269 = -t253 * t356 + t258 * t210 + rSges(6,3) * t382 + (-t245 * t386 - t249 * t389) * rSges(6,1);
t346 = -t249 * rSges(6,1) * t390 - t251 * t356 - t258 * t357;
t359 = t139 * t382 + t251 * (t434 * t258 + t346) + t253 * t269;
t345 = -t251 * rSges(5,1) * t387 + t443 * rSges(5,2);
t358 = t147 * t382 + t251 * (t433 * t258 + t345) + t253 * (-t383 * t409 + (-t250 * t383 - t252 * t386) * rSges(5,1) + t372);
t265 = cos(qJ(1));
t355 = t265 * t408;
t353 = pkin(3) * t362;
t352 = t263 * t408;
t344 = -rSges(4,1) * t343 + t444 * rSges(4,2);
t340 = t386 / 0.2e1;
t339 = t382 / 0.2e1;
t338 = -pkin(2) - t414;
t337 = -t200 - t418;
t335 = -t246 - t413;
t334 = -t214 - t412;
t202 = t253 * rSges(3,1) - rSges(3,2) * t251;
t282 = Icges(5,5) * t258 - t392;
t85 = t251 * t282 + t258 * t293;
t329 = -t143 * t257 + t85;
t84 = t253 * t282 - t315 * t386;
t328 = -t144 * t257 + t84;
t280 = Icges(5,6) * t258 - t393;
t83 = t251 * t280 + t258 * t290;
t327 = t145 * t257 + t83;
t82 = t253 * t280 - t312 * t386;
t326 = t146 * t257 + t82;
t325 = t185 * t251 + t258 * t233;
t110 = (-t259 - t266) * t253 + t435;
t25 = t251 * t110 + t253 * t111 + t75;
t159 = (-rSges(6,2) * t244 + t412) * t249;
t322 = -pkin(4) * t384 - t159;
t176 = -rSges(3,1) * t382 + rSges(3,2) * t386;
t321 = t361 - t417;
t37 = -t141 * t253 - t251 * t305;
t38 = -t142 * t253 - t440;
t39 = t141 * t251 - t437;
t40 = t142 * t251 - t253 * t304;
t197 = Icges(5,5) * t250 + Icges(5,6) * t252;
t278 = Icges(5,3) * t258 - t197 * t257;
t80 = t253 * t278 - t309 * t386;
t81 = t251 * t278 + t258 * t287;
t320 = (t251 * t38 - t253 * t37) * t386 + (t251 * t40 - t253 * t39) * t382 + t251 * ((t251 * t80 + (t39 + t440) * t258) * t251 + (t40 * t258 + (t143 * t384 + t145 * t387 + t250 * t83 - t252 * t85 - t398) * t253 + (t397 - t81 + (t145 * t258 + t328) * t252 + (-t143 * t258 - t326) * t250) * t251) * t253) + t361;
t201 = -rSges(3,1) * t251 - rSges(3,2) * t253;
t318 = t335 * t251;
t222 = Icges(4,5) * t262 + Icges(4,6) * t264;
t303 = t162 * t264 + t164 * t262;
t301 = t163 * t264 + t165 * t262;
t296 = t110 * t382 + t251 * (t369 * t382 + t325 + t368) + t253 * (t323 + (t380 - t435) * t258 + t374) + t359;
t167 = t442 * t253 + t236;
t175 = t201 * t258;
t288 = t310 * t253;
t155 = t311 * t249;
t156 = t314 * t249;
t272 = t182 * t258 + (t156 - t395) * t245 + (-t155 - t394) * t244;
t285 = (t244 * t332 + t245 * t330 + t432 * t251 + t272 * t253) * t424 + (t244 * t333 + t245 * t331 + t272 * t251 - t432 * t253) * t423 + (t135 * t245 + t137 * t244 - t182 * t253 - t251 * t299) * t340 + (t136 * t245 + t138 * t244 + t182 * t251 - t253 * t299) * t339;
t120 = t167 - t365;
t102 = t140 + t373;
t284 = t322 - t353;
t4 = (t253 * t81 + (t38 + t437) * t258) * t253 + (t37 * t258 + (-t144 * t384 - t146 * t387 - t250 * t82 + t252 * t84 + t397) * t251 + (-t398 - t80 + (t146 * t258 - t329) * t252 + (-t144 * t258 + t327) * t250) * t253) * t251;
t283 = (-t4 - t2) * t253 + t320;
t119 = t251 * t338 + t242 + t366;
t109 = t148 - t324;
t276 = Icges(4,5) * t258 - qJD(3) * t224;
t275 = Icges(4,6) * t258 - qJD(3) * t223;
t274 = Icges(4,3) * t258 - qJD(3) * t222;
t101 = t251 * t334 + t370 + t381;
t108 = t318 + t367 - t380;
t172 = t312 * t257;
t173 = t315 * t257;
t271 = t197 * t258 + (t173 - t393) * t252 + (-t172 - t392) * t250;
t273 = t285 + (t250 * t328 + t431 * t251 + t252 * t326 + t271 * t253) * t424 + (t250 * t329 + t271 * t251 + t252 * t327 - t431 * t253) * t423 + (t143 * t252 + t145 * t250 - t197 * t253 - t251 * t298) * t340 + (t144 * t252 + t146 * t250 + t197 * t251 - t253 * t298) * t339;
t74 = (t338 * t253 + (-rSges(4,3) - pkin(7)) * t251) * t258 - t344;
t204 = t313 * qJD(3);
t205 = t316 * qJD(3);
t270 = -t204 * t262 + t205 * t264 + t222 * t258 + (-t223 * t264 - t224 * t262) * qJD(3);
t46 = (t253 * t335 - t235) * t258 - t345 + t368;
t32 = (t253 * t334 - t234) * t258 - t325 - t346;
t268 = t245 * t155 + t244 * t156 + t252 * t172 + t250 * t173 - t183 * t391 + t184 * t388 - t198 * t387 + t199 * t384 + t264 * t204 + t262 * t205 - t223 * t363 + t224 * t362;
t73 = -rSges(4,2) * t253 * t362 - pkin(2) * t386 + t229 + (-t258 * t385 - t341) * rSges(4,1) + t371;
t267 = t273 + (-qJD(3) * t300 + t430 * t251 + t270 * t253 + t262 * (t253 * t276 - t316 * t386) + t264 * (t253 * t275 - t313 * t386)) * t424 + (-qJD(3) * t302 + (t251 * t276 + t258 * t294) * t262 + t270 * t251 - t430 * t253 + t264 * (t251 * t275 + t258 * t291)) * t423 + (-t222 * t253 - t251 * t297 + t303) * t340 + (t222 * t251 - t253 * t297 + t301) * t339;
t31 = -t214 * t386 + t269 + t374;
t45 = t258 * t318 + (-t200 * t257 - t354 - t378) * t253 + t372;
t256 = t265 * pkin(1);
t209 = pkin(3) * t348;
t178 = t202 + t256;
t177 = t201 - t416;
t166 = rSges(4,1) * t385 - t366;
t161 = Icges(4,3) * t251 + t288;
t160 = -Icges(4,3) * t253 + t251 * t310;
t153 = t176 - t355;
t152 = t175 - t352;
t150 = t337 * t253;
t149 = t337 * t251;
t128 = t336 * t253;
t127 = t336 * t251;
t114 = t295 * t251;
t113 = t120 + t256;
t112 = t119 - t416;
t105 = t109 + t256;
t104 = t108 - t416;
t96 = t251 * t274 + t258 * t288;
t95 = t253 * t274 - t310 * t386;
t92 = t102 + t256;
t91 = t101 - t416;
t79 = t444 * pkin(3) - t200 * t382 - t396;
t78 = t200 * t386 + t209 + (-t174 - t353) * t253;
t60 = t74 - t355;
t59 = t73 - t352;
t58 = t443 * pkin(4) - t159 * t251 - t188 * t382;
t57 = t253 * t322 + t375;
t50 = t161 * t251 - t253 * t300;
t49 = t160 * t251 - t438;
t48 = -t161 * t253 - t441;
t47 = -t160 * t253 - t251 * t302;
t44 = t251 * t284 + t401;
t43 = t253 * t284 + t209 + t375;
t42 = t46 - t355;
t41 = t45 - t352;
t30 = t32 - t355;
t29 = t31 - t352;
t28 = t86 + t376;
t22 = -t148 * t386 + t358;
t17 = -t140 * t386 + t359;
t14 = t25 + t376;
t7 = (-t132 - t148) * t386 + t358 + t360;
t6 = t377 * t386 + t296;
t5 = (-t132 + t377) * t386 + t296 + t360;
t1 = [t268 + (t29 * t92 + t30 * t91) * t425 + (t104 * t42 + t105 * t41) * t426 + (t112 * t60 + t113 * t59) * t427 + (t152 * t178 + t153 * t177) * t428; t268 + m(6) * (t101 * t30 + t102 * t29 + t31 * t92 + t32 * t91) + m(5) * (t104 * t46 + t105 * t45 + t108 * t42 + t109 * t41) + m(4) * (t112 * t74 + t113 * t73 + t119 * t60 + t120 * t59) + m(3) * (t152 * t202 + t153 * t201 + t175 * t178 + t176 * t177); t268 + (t101 * t32 + t102 * t31) * t425 + (t108 * t46 + t109 * t45) * t426 + (t119 * t74 + t120 * t73) * t427 + (t175 * t202 + t176 * t201) * t428; m(5) * (t104 * t78 + t105 * t79 + t149 * t41 + t150 * t42) + m(6) * (t114 * t29 + t115 * t30 + t43 * t91 + t44 * t92) + ((-t113 * t258 - t60) * t253 + (t112 * t258 - t59) * t251) * t421 + (-t112 * t253 - t113 * t251) * t422 + t267; ((-t120 * t258 - t74) * t253 + (t119 * t258 - t73) * t251) * t421 + (-t119 * t253 - t120 * t251) * t422 + m(5) * (t108 * t78 + t109 * t79 + t149 * t45 + t150 * t46) + m(6) * (t101 * t43 + t102 * t44 + t114 * t31 + t115 * t32) + t267; (t114 * t44 + t115 * t43 + t14 * t5) * t425 - t417 - t253 * t4 + (t149 * t79 + t150 * t78 + t28 * t7) * t426 + ((t166 * t251 + t167 * t253) * (((-t167 + t236) * t258 + t344) * t251 + (-qJD(3) * t230 * t253 + t258 * t166 + t371) * t253) + t364 * t230 * t208) * t427 + (t251 * t48 - t253 * t47) * t386 - t253 * ((t253 * t96 + (t48 + t438) * t258) * t253 + (t47 * t258 + (-t163 * t362 - t165 * t363) * t251 + (t303 * qJD(3) - t258 * t300 - t95) * t253) * t251) + (t251 * t50 - t253 * t49) * t382 + t251 * ((t251 * t95 + (t49 + t441) * t258) * t251 + (t50 * t258 + (t162 * t362 + t164 * t363) * t253 + (-t301 * qJD(3) - t258 * t302 - t96) * t251) * t253) + t320; ((-t105 * t258 - t42) * t253 + (t104 * t258 - t41) * t251) * t419 + m(6) * (t127 * t29 + t128 * t30 + t57 * t91 + t58 * t92) + (-t104 * t253 - t105 * t251) * t420 + t273; ((-t109 * t258 - t46) * t253 + (t108 * t258 - t45) * t251) * t419 + (-t108 * t253 - t109 * t251) * t420 + m(6) * (t101 * t57 + t102 * t58 + t127 * t31 + t128 * t32) + t273; m(6) * (t114 * t58 + t115 * t57 + t127 * t44 + t128 * t43 + t14 * t6 + t25 * t5) + m(5) * (-t150 * t174 * t253 - t149 * t396 + t22 * t28 + t7 * t86) + ((-t149 * t258 - t78) * t253 + (t150 * t258 - t79) * t251) * t419 + t283; (t174 * t200 * t364 + t22 * t86) * t426 + (t127 * t58 + t128 * t57 + t25 * t6) * t425 + t283; m(6) * ((-t251 * t92 - t253 * t91) * t159 + ((-t258 * t92 - t30) * t253 + (t258 * t91 - t29) * t251) * t188) + t285; m(6) * ((-t101 * t253 - t102 * t251) * t159 + ((-t102 * t258 - t32) * t253 + (t101 * t258 - t31) * t251) * t188) + t285; m(6) * (t14 * t17 + t5 * t75 + (-t114 * t251 - t115 * t253) * t159 + ((-t114 * t258 - t43) * t253 + (-t44 + t401) * t251) * t188) + t321; m(6) * (t17 * t25 + t6 * t75 + (-t127 * t251 - t128 * t253) * t159 + ((-t127 * t258 - t57) * t253 + (t128 * t258 - t58) * t251) * t188) + t321; (t159 * t188 * t364 + t17 * t75) * t425 + t321;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

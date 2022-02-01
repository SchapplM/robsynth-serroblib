% Calculate time derivative of joint inertia matrix for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:48:45
% EndTime: 2022-01-20 11:49:01
% DurationCPUTime: 8.38s
% Computational Cost: add. (12695->517), mult. (10040->686), div. (0->0), fcn. (7628->8), ass. (0->314)
t256 = qJ(3) + qJ(4);
t249 = cos(t256);
t247 = sin(t256);
t397 = Icges(6,4) * t247;
t193 = Icges(6,2) * t249 + t397;
t399 = Icges(5,4) * t247;
t194 = Icges(5,2) * t249 + t399;
t458 = t193 + t194;
t396 = Icges(6,4) * t249;
t195 = Icges(6,1) * t247 + t396;
t398 = Icges(5,4) * t249;
t196 = Icges(5,1) * t247 + t398;
t457 = t195 + t196;
t257 = qJ(1) + qJ(2);
t248 = sin(t257);
t250 = cos(t257);
t308 = -Icges(6,2) * t247 + t396;
t286 = t308 * t250;
t134 = Icges(6,6) * t248 + t286;
t309 = -Icges(5,2) * t247 + t398;
t287 = t309 * t250;
t136 = Icges(5,6) * t248 + t287;
t456 = t134 + t136;
t311 = Icges(6,1) * t249 - t397;
t289 = t311 * t250;
t138 = Icges(6,5) * t248 + t289;
t312 = Icges(5,1) * t249 - t399;
t290 = t312 * t250;
t140 = Icges(5,5) * t248 + t290;
t455 = t138 + t140;
t254 = qJD(3) + qJD(4);
t454 = (t308 + t309) * t254;
t453 = (t311 + t312) * t254;
t191 = Icges(6,5) * t247 + Icges(6,6) * t249;
t192 = Icges(5,5) * t247 + Icges(5,6) * t249;
t443 = t191 + t192;
t452 = t458 * t247 - t457 * t249;
t374 = t254 * t194;
t375 = t254 * t193;
t451 = -t374 - t375;
t133 = -Icges(6,6) * t250 + t248 * t308;
t135 = -Icges(5,6) * t250 + t248 * t309;
t450 = t133 + t135;
t137 = -Icges(6,5) * t250 + t248 * t311;
t139 = -Icges(5,5) * t250 + t248 * t312;
t449 = t137 + t139;
t260 = cos(qJ(3));
t244 = t260 * pkin(3) + pkin(2);
t210 = pkin(4) * t249 + t244;
t262 = -pkin(8) - pkin(7);
t253 = qJ(5) - t262;
t386 = t247 * t250;
t234 = t248 * rSges(6,3);
t381 = t249 * t250;
t425 = rSges(6,1) * t381 + t234;
t100 = -rSges(6,2) * t386 + t250 * t210 + t253 * t248 + t425;
t305 = Icges(6,5) * t249 - Icges(6,6) * t247;
t129 = -Icges(6,3) * t250 + t248 * t305;
t306 = Icges(5,5) * t249 - Icges(5,6) * t247;
t131 = -Icges(5,3) * t250 + t248 * t306;
t304 = t133 * t247 - t137 * t249;
t427 = t250 * t304;
t302 = t135 * t247 - t139 * t249;
t428 = t250 * t302;
t448 = t427 + t428 + (-t129 - t131) * t248;
t283 = t305 * t250;
t130 = Icges(6,3) * t248 + t283;
t284 = t306 * t250;
t132 = Icges(5,3) * t248 + t284;
t301 = t136 * t247 - t140 * t249;
t303 = t134 * t247 - t138 * t249;
t447 = (-t301 - t303) * t250 + (t130 + t132) * t248;
t387 = t247 * t248;
t363 = rSges(6,2) * t387 + t250 * rSges(6,3);
t384 = t248 * t249;
t365 = t210 - t244;
t426 = t365 * t248;
t446 = (-t253 - t262) * t250 + t426 + rSges(6,1) * t384 - t363;
t321 = -t250 * t244 + t248 * t262;
t445 = t321 + t100;
t255 = qJD(1) + qJD(2);
t258 = sin(qJ(3));
t358 = qJD(3) * t258;
t353 = pkin(3) * t358;
t385 = t247 * t254;
t444 = -pkin(4) * t385 + t253 * t255 - t353;
t379 = t250 * t254;
t383 = t248 * t255;
t282 = -t247 * t379 - t249 * t383;
t345 = t249 * t379;
t349 = t247 * t383;
t378 = t250 * t255;
t442 = rSges(6,1) * t282 + rSges(6,3) * t378 + (-t345 + t349) * rSges(6,2) + qJD(5) * t248 + t444 * t250;
t388 = t196 * t254;
t389 = t195 * t254;
t441 = t443 * t255 + (t451 + t453) * t249 + (-t388 - t389 - t454) * t247;
t350 = t248 * t385;
t380 = t249 * t254;
t435 = -t247 * t378 - t248 * t380;
t439 = rSges(6,1) * t350 - t435 * rSges(6,2) + qJD(5) * t250 - t444 * t248;
t276 = Icges(6,6) * t255 - t375;
t76 = t250 * t276 - t308 * t383;
t277 = Icges(5,6) * t255 - t374;
t78 = t250 * t277 - t309 * t383;
t438 = t455 * t254 + t76 + t78;
t278 = Icges(6,5) * t255 - t389;
t80 = t250 * t278 - t311 * t383;
t279 = Icges(5,5) * t255 - t388;
t82 = t250 * t279 - t312 * t383;
t437 = -t456 * t254 + t80 + t82;
t357 = qJD(3) * t260;
t373 = t255 * t258;
t436 = -t248 * t357 - t250 * t373;
t403 = rSges(4,2) * t258;
t406 = rSges(4,1) * t260;
t434 = -t403 + t406;
t433 = t452 * t255 + (t305 + t306) * t254;
t400 = Icges(4,4) * t260;
t310 = -Icges(4,2) * t258 + t400;
t288 = t310 * t250;
t156 = Icges(4,6) * t248 + t288;
t401 = Icges(4,4) * t258;
t313 = Icges(4,1) * t260 - t401;
t291 = t313 * t250;
t158 = Icges(4,5) * t248 + t291;
t297 = t156 * t258 - t158 * t260;
t432 = t248 * t297;
t431 = t248 * t301;
t430 = t248 * t303;
t155 = -Icges(4,6) * t250 + t248 * t310;
t157 = -Icges(4,5) * t250 + t248 * t313;
t299 = t155 * t258 - t157 * t260;
t429 = t250 * t299;
t235 = t248 * rSges(5,3);
t424 = rSges(5,1) * t381 + t235;
t221 = Icges(4,2) * t260 + t401;
t222 = Icges(4,1) * t258 + t400;
t294 = t221 * t258 - t222 * t260;
t307 = Icges(4,5) * t260 - Icges(4,6) * t258;
t421 = t307 * qJD(3) + t255 * t294;
t377 = t250 * t262;
t407 = pkin(2) - t244;
t420 = t248 * t407 - t377;
t419 = 2 * m(3);
t418 = 2 * m(4);
t417 = 2 * m(5);
t416 = 2 * m(6);
t415 = t248 / 0.2e1;
t414 = -t250 / 0.2e1;
t208 = t434 * qJD(3);
t413 = m(4) * t208;
t228 = rSges(4,1) * t258 + rSges(4,2) * t260;
t412 = m(4) * t228;
t405 = rSges(5,1) * t249;
t172 = (-rSges(5,2) * t247 + t405) * t254;
t411 = m(5) * t172;
t198 = rSges(5,1) * t247 + rSges(5,2) * t249;
t410 = m(5) * t198;
t409 = pkin(3) * t258;
t241 = t248 * pkin(7);
t259 = sin(qJ(1));
t408 = t259 * pkin(1);
t404 = rSges(6,1) * t249;
t402 = pkin(1) * qJD(1);
t236 = t248 * rSges(4,3);
t197 = rSges(6,1) * t247 + rSges(6,2) * t249;
t332 = -pkin(4) * t247 - t197;
t292 = t332 - t409;
t114 = t292 * t250;
t395 = t114 * t255;
t394 = t129 * t255;
t393 = t130 * t255;
t392 = t131 * t255;
t391 = t132 * t255;
t390 = t172 * t248;
t382 = t248 * t260;
t372 = t255 * t262;
t242 = t250 * pkin(7);
t125 = t242 - t420;
t360 = -t250 * pkin(2) - t241;
t126 = -t321 + t360;
t370 = t248 * t125 + t250 * t126;
t362 = rSges(5,2) * t387 + t250 * rSges(5,3);
t142 = rSges(5,1) * t384 - t362;
t144 = -rSges(5,2) * t386 + t424;
t84 = t248 * t142 + t250 * t144;
t369 = pkin(4) * t349 + t197 * t383;
t367 = rSges(5,2) * t349 + rSges(5,3) * t378;
t346 = t248 * t373;
t366 = rSges(4,2) * t346 + rSges(4,3) * t378;
t339 = t248 * t358;
t364 = pkin(3) * t339 + t248 * t372;
t361 = t250 * rSges(4,3) + t248 * t403;
t359 = t248 ^ 2 + t250 ^ 2;
t227 = pkin(7) * t378;
t337 = t250 * t358;
t320 = pkin(3) * t337;
t356 = t125 * t378 + t248 * ((-t250 * t407 - t241) * t255 - t364) + t250 * (t255 * t420 - t227 - t320);
t341 = -rSges(5,1) * t350 + t435 * rSges(5,2);
t355 = t142 * t378 + t248 * (t255 * t424 + t341) + t250 * (rSges(5,1) * t282 - rSges(5,2) * t345 + t367);
t261 = cos(qJ(1));
t354 = t261 * t402;
t352 = pkin(3) * t357;
t351 = t259 * t402;
t340 = -rSges(4,1) * t339 + t436 * rSges(4,2);
t336 = t383 / 0.2e1;
t335 = t378 / 0.2e1;
t334 = -pkin(2) - t406;
t333 = -t198 - t409;
t331 = -t244 - t405;
t330 = -t210 - t404;
t200 = t250 * rSges(3,1) - rSges(3,2) * t248;
t81 = t248 * t278 + t255 * t289;
t329 = -t133 * t254 + t81;
t83 = t248 * t279 + t255 * t290;
t327 = -t135 * t254 + t83;
t77 = t248 * t276 + t255 * t286;
t325 = t137 * t254 + t77;
t79 = t248 * t277 + t255 * t287;
t323 = t139 * t254 + t79;
t24 = t446 * t248 + t445 * t250;
t171 = (-rSges(6,2) * t247 + t404) * t254;
t319 = -pkin(4) * t380 - t171;
t174 = -rSges(3,1) * t378 + rSges(3,2) * t383;
t32 = -t129 * t250 - t248 * t304;
t33 = -t130 * t250 - t430;
t34 = -t131 * t250 - t248 * t302;
t35 = -t132 * t250 - t431;
t274 = Icges(6,3) * t255 - t191 * t254;
t72 = t250 * t274 - t305 * t383;
t73 = t248 * t274 + t255 * t283;
t275 = Icges(5,3) * t255 - t192 * t254;
t74 = t250 * t275 - t306 * t383;
t75 = t248 * t275 + t255 * t284;
t318 = ((-t32 - t34) * t383 + t448 * t378) * t250 + (((t72 + t74) * t248 + (t430 + t431 - t448) * t255) * t248 + (t33 + t35) * t383 + t447 * t378 + ((t393 - t73 + t391 - t75) * t248 + (t450 * t380 + t449 * t385 - t392 - t394) * t250 + t447 * t255 + ((t449 * t255 + t437) * t248 + (-t81 - t83) * t250) * t249 + ((-t450 * t255 - t438) * t248 + (t77 + t79) * t250) * t247) * t250) * t248;
t199 = -rSges(3,1) * t248 - rSges(3,2) * t250;
t314 = t331 * t248;
t220 = Icges(4,5) * t258 + Icges(4,6) * t260;
t300 = t155 * t260 + t157 * t258;
t298 = t156 * t260 + t158 * t258;
t293 = t446 * t378 + (t320 + (t377 - t426) * t255 + t442) * t250 + (t255 * t425 + t365 * t378 + t364 - t439) * t248;
t160 = t434 * t250 + t236;
t173 = t199 * t255;
t285 = t307 * t250;
t118 = t160 - t360;
t281 = t319 - t352;
t3 = (t250 * t73 + (t33 + t427) * t255) * t250 + (t32 * t255 + (-t134 * t380 - t138 * t385 - t247 * t76 + t249 * t80 + t393) * t248 + (-t394 - t72 + (t138 * t255 - t329) * t249 + (-t134 * t255 + t325) * t247) * t250) * t248;
t4 = (t250 * t75 + (t35 + t428) * t255) * t250 + (t34 * t255 + (-t136 * t380 - t140 * t385 - t247 * t78 + t249 * t82 + t391) * t248 + (-t392 - t74 + (t140 * t255 - t327) * t249 + (-t136 * t255 + t323) * t247) * t250) * t248;
t280 = (-t4 - t3) * t250 + t318;
t117 = t248 * t334 + t242 + t361;
t107 = t144 - t321;
t273 = Icges(4,5) * t255 - qJD(3) * t222;
t272 = Icges(4,6) * t255 - qJD(3) * t221;
t271 = Icges(4,3) * t255 - qJD(3) * t220;
t99 = t248 * t330 + t250 * t253 + t363;
t106 = t314 + t362 - t377;
t202 = t310 * qJD(3);
t203 = t313 * qJD(3);
t270 = t260 * t202 + t258 * t203 - t221 * t358 + t222 * t357 + t453 * t247 + t454 * t249 + t457 * t380;
t269 = (t437 * t247 + t433 * t248 + t438 * t249 + t441 * t250) * t415 + (-t433 * t250 + (t323 + t325) * t249 + t441 * t248 + (t327 + t329) * t247) * t414 + (t449 * t247 - t248 * t452 + t450 * t249 - t443 * t250) * t336 + (t455 * t247 + t443 * t248 + t456 * t249 - t250 * t452) * t335;
t65 = (t334 * t250 + (-rSges(4,3) - pkin(7)) * t248) * t255 - t340;
t266 = -t202 * t258 + t203 * t260 + t220 * t255 + (-t221 * t260 - t222 * t258) * qJD(3);
t43 = (t250 * t331 - t235) * t255 - t341 + t364;
t31 = (t250 * t330 - t234) * t255 + t439;
t264 = -t458 * t385 + t270;
t64 = -rSges(4,2) * t250 * t357 - pkin(2) * t383 + t227 + (-t255 * t382 - t337) * rSges(4,1) + t366;
t263 = t269 + (-qJD(3) * t297 + t248 * t421 + t266 * t250 + t258 * (t250 * t273 - t313 * t383) + t260 * (t250 * t272 - t310 * t383)) * t415 + (-qJD(3) * t299 + t266 * t248 - t250 * t421 + t258 * (t248 * t273 + t255 * t291) + t260 * (t248 * t272 + t255 * t288)) * t414 + (-t220 * t250 - t248 * t294 + t300) * t336 + (t220 * t248 - t250 * t294 + t298) * t335;
t30 = -t210 * t383 + t442;
t42 = t255 * t314 + (-t198 * t254 - t353 - t372) * t250 + t367;
t252 = t261 * pkin(1);
t209 = pkin(3) * t346;
t176 = t200 + t252;
t175 = t199 - t408;
t159 = rSges(4,1) * t382 - t361;
t154 = Icges(4,3) * t248 + t285;
t153 = -Icges(4,3) * t250 + t248 * t307;
t148 = t174 - t354;
t147 = t173 - t351;
t146 = t333 * t250;
t145 = t333 * t248;
t128 = t332 * t250;
t127 = t332 * t248;
t113 = t292 * t248;
t111 = t118 + t252;
t110 = t117 - t408;
t103 = t107 + t252;
t102 = t106 - t408;
t94 = t248 * t271 + t255 * t285;
t93 = t250 * t271 - t307 * t383;
t92 = t100 + t252;
t91 = t99 - t408;
t71 = t436 * pkin(3) - t198 * t378 - t390;
t70 = t198 * t383 + t209 + (-t172 - t352) * t250;
t59 = t435 * pkin(4) - t171 * t248 - t197 * t378;
t58 = t250 * t319 + t369;
t53 = t65 - t354;
t52 = t64 - t351;
t51 = t154 * t248 - t297 * t250;
t50 = t153 * t248 - t429;
t49 = -t154 * t250 - t432;
t48 = -t153 * t250 - t299 * t248;
t47 = t248 * t281 + t395;
t46 = t250 * t281 + t209 + t369;
t41 = t43 - t354;
t40 = t42 - t351;
t29 = t31 - t354;
t28 = t30 - t351;
t27 = t84 + t370;
t21 = -t144 * t383 + t355;
t12 = t24 + t370;
t7 = (-t126 - t144) * t383 + t355 + t356;
t6 = -t383 * t445 + t293;
t5 = (-t126 - t445) * t383 + t293 + t356;
t1 = [(t28 * t92 + t29 * t91) * t416 + (t102 * t41 + t103 * t40) * t417 + (t110 * t53 + t111 * t52) * t418 + (t147 * t176 + t148 * t175) * t419 + t270 + t451 * t247; m(6) * (t100 * t28 + t29 * t99 + t30 * t92 + t31 * t91) + m(5) * (t102 * t43 + t103 * t42 + t106 * t41 + t107 * t40) + m(4) * (t110 * t65 + t111 * t64 + t117 * t53 + t118 * t52) + m(3) * (t147 * t200 + t148 * t199 + t173 * t176 + t174 * t175) + t264; (t100 * t30 + t31 * t99) * t416 + (t106 * t43 + t107 * t42) * t417 + (t117 * t65 + t118 * t64) * t418 + (t173 * t200 + t174 * t199) * t419 + t264; t263 + m(6) * (t113 * t28 + t114 * t29 + t46 * t91 + t47 * t92) + m(5) * (t102 * t70 + t103 * t71 + t145 * t40 + t146 * t41) + (-t110 * t250 - t111 * t248) * t413 + ((-t111 * t255 - t53) * t250 + (t110 * t255 - t52) * t248) * t412; t263 + m(6) * (t100 * t47 + t113 * t30 + t114 * t31 + t46 * t99) + m(5) * (t106 * t70 + t107 * t71 + t145 * t42 + t146 * t43) + (-t117 * t250 - t118 * t248) * t413 + ((-t118 * t255 - t65) * t250 + (t117 * t255 - t64) * t248) * t412; (t113 * t47 + t114 * t46 + t12 * t5) * t416 - t250 * t3 - t250 * t4 + (t145 * t71 + t146 * t70 + t27 * t7) * t417 + ((t159 * t248 + t160 * t250) * (((-t160 + t236) * t255 + t340) * t248 + (-qJD(3) * t228 * t250 + t255 * t159 + t366) * t250) + t359 * t228 * t208) * t418 + (t248 * t49 - t250 * t48) * t383 - t250 * ((t250 * t94 + (t49 + t429) * t255) * t250 + (t48 * t255 + (-t156 * t357 - t158 * t358) * t248 + (t300 * qJD(3) - t255 * t297 - t93) * t250) * t248) + (t248 * t51 - t250 * t50) * t378 + t248 * ((t248 * t93 + (t50 + t432) * t255) * t248 + (t51 * t255 + (t155 * t357 + t157 * t358) * t250 + (-t298 * qJD(3) - t299 * t255 - t94) * t248) * t250) + t318; t269 + (-t102 * t250 - t103 * t248) * t411 + ((-t103 * t255 - t41) * t250 + (t102 * t255 - t40) * t248) * t410 + m(6) * (t127 * t28 + t128 * t29 + t58 * t91 + t59 * t92); t269 + (-t106 * t250 - t107 * t248) * t411 + ((-t107 * t255 - t43) * t250 + (t106 * t255 - t42) * t248) * t410 + m(6) * (t100 * t59 + t127 * t30 + t128 * t31 + t58 * t99); m(6) * (t113 * t59 + t114 * t58 + t12 * t6 + t127 * t47 + t128 * t46 + t24 * t5) + m(5) * (-t146 * t172 * t250 - t145 * t390 + t21 * t27 + t7 * t84) + ((-t145 * t255 - t70) * t250 + (t146 * t255 - t71) * t248) * t410 + t280; (t172 * t198 * t359 + t21 * t84) * t417 + (t127 * t59 + t128 * t58 + t24 * t6) * t416 + t280; m(6) * ((t255 * t91 - t28) * t250 + (t255 * t92 + t29) * t248); m(6) * ((t255 * t99 - t30) * t250 + (t100 * t255 + t31) * t248); m(6) * ((-t47 + t395) * t250 + (t113 * t255 + t46) * t248); m(6) * ((t128 * t255 - t59) * t250 + (t127 * t255 + t58) * t248); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

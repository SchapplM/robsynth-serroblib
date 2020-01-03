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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:11:05
% EndTime: 2020-01-03 12:11:23
% DurationCPUTime: 10.03s
% Computational Cost: add. (12695->547), mult. (10040->725), div. (0->0), fcn. (7628->8), ass. (0->314)
t260 = qJ(3) + qJ(4);
t252 = cos(t260);
t250 = sin(t260);
t403 = Icges(6,4) * t250;
t194 = Icges(6,2) * t252 + t403;
t405 = Icges(5,4) * t250;
t195 = Icges(5,2) * t252 + t405;
t457 = t194 + t195;
t402 = Icges(6,4) * t252;
t196 = Icges(6,1) * t250 + t402;
t404 = Icges(5,4) * t252;
t197 = Icges(5,1) * t250 + t404;
t456 = t196 + t197;
t258 = qJD(3) + qJD(4);
t312 = -Icges(6,2) * t250 + t402;
t313 = -Icges(5,2) * t250 + t404;
t455 = (t312 + t313) * t258;
t315 = Icges(6,1) * t252 - t403;
t316 = Icges(5,1) * t252 - t405;
t454 = (t315 + t316) * t258;
t192 = Icges(6,5) * t250 + Icges(6,6) * t252;
t193 = Icges(5,5) * t250 + Icges(5,6) * t252;
t450 = -t192 - t193;
t448 = t457 * t250 - t456 * t252;
t264 = cos(qJ(3));
t246 = t264 * pkin(3) + pkin(2);
t208 = pkin(4) * t252 + t246;
t261 = qJ(1) + qJ(2);
t251 = sin(t261);
t253 = cos(t261);
t266 = -pkin(8) - pkin(7);
t257 = -qJ(5) + t266;
t377 = t251 * t252;
t380 = t250 * t251;
t99 = rSges(6,1) * t377 - rSges(6,2) * t380 + t251 * t208 + (-rSges(6,3) + t257) * t253;
t374 = t252 * t253;
t379 = t250 * t253;
t453 = rSges(6,1) * t374 - rSges(6,2) * t379 + t251 * rSges(6,3) + t253 * t208;
t235 = t253 * t266;
t357 = t251 * t246 + t235;
t452 = -t357 + t99;
t210 = t253 * t246;
t353 = t257 - t266;
t451 = t251 * t353 + t210 - t453;
t259 = qJD(1) + qJD(2);
t309 = Icges(6,5) * t252 - Icges(6,6) * t250;
t310 = Icges(5,5) * t252 - Icges(5,6) * t250;
t449 = -t448 * t259 + (-t309 - t310) * t258;
t262 = sin(qJ(3));
t352 = qJD(3) * t262;
t349 = pkin(3) * t352;
t378 = t250 * t258;
t182 = -pkin(4) * t378 - t349;
t371 = t253 * t259;
t341 = t252 * t371;
t376 = t251 * t259;
t447 = rSges(6,1) * t341 + rSges(6,3) * t376 + t251 * t182 + t208 * t371;
t381 = t197 * t258;
t382 = t196 * t258;
t383 = t195 * t258;
t384 = t194 * t258;
t446 = t450 * t259 + (t384 + t383 - t454) * t252 + (t382 + t381 + t455) * t250;
t199 = rSges(5,1) * t250 + rSges(5,2) * t252;
t416 = m(5) * t199;
t129 = -Icges(6,3) * t253 + t251 * t309;
t133 = -Icges(6,6) * t253 + t251 * t312;
t290 = t315 * t251;
t137 = -Icges(6,5) * t253 + t290;
t308 = t133 * t250 - t137 * t252;
t283 = t308 * t253;
t36 = -t129 * t251 + t283;
t130 = -Icges(6,3) * t251 - t253 * t309;
t134 = -Icges(6,6) * t251 - t253 * t312;
t138 = -Icges(6,5) * t251 - t253 * t315;
t307 = t134 * t250 - t138 * t252;
t37 = -t130 * t251 + t253 * t307;
t285 = t310 * t251;
t131 = -Icges(5,3) * t253 + t285;
t288 = t313 * t251;
t135 = -Icges(5,6) * t253 + t288;
t139 = -Icges(5,5) * t253 + t251 * t316;
t306 = t135 * t250 - t139 * t252;
t282 = t306 * t253;
t38 = -t131 * t251 + t282;
t132 = -Icges(5,3) * t251 - t253 * t310;
t136 = -Icges(5,6) * t251 - t253 * t313;
t140 = -Icges(5,5) * t251 - t253 * t316;
t305 = t136 * t250 - t140 * t252;
t39 = -t132 * t251 + t253 * t305;
t445 = (t36 + t38) * t253 + (t37 + t39) * t251;
t406 = Icges(4,4) * t264;
t314 = -Icges(4,2) * t262 + t406;
t157 = -Icges(4,6) * t251 - t253 * t314;
t407 = Icges(4,4) * t262;
t317 = Icges(4,1) * t264 - t407;
t159 = -Icges(4,5) * t251 - t253 * t317;
t301 = t157 * t262 - t159 * t264;
t444 = t251 * t301;
t443 = t251 * t305;
t442 = t251 * t307;
t441 = -t251 * rSges(3,1) - t253 * rSges(3,2);
t344 = t250 * t376;
t440 = rSges(6,2) * t344 + rSges(6,3) * t371 + qJD(5) * t251 + t253 * t182;
t439 = -Icges(5,5) * t259 + t381;
t438 = -Icges(6,5) * t259 + t382;
t437 = -Icges(5,6) * t259 + t383;
t436 = -Icges(6,6) * t259 + t384;
t435 = -Icges(6,3) * t259 + t192 * t258;
t434 = -Icges(5,3) * t259 + t193 * t258;
t217 = Icges(4,5) * t262 + Icges(4,6) * t264;
t433 = -Icges(4,3) * t259 + qJD(3) * t217;
t218 = Icges(4,2) * t264 + t407;
t432 = -Icges(4,6) * t259 + qJD(3) * t218;
t219 = Icges(4,1) * t262 + t406;
t431 = -Icges(4,5) * t259 + qJD(3) * t219;
t203 = t314 * qJD(3);
t204 = t317 * qJD(3);
t430 = t203 * t262 - t204 * t264 - t217 * t259 + (t218 * t264 + t219 * t262) * qJD(3);
t351 = qJD(3) * t264;
t373 = t252 * t258;
t268 = t264 * t203 + t262 * t204 - t218 * t352 + t219 * t351 + t454 * t250 + t455 * t252 + t456 * t373 - t457 * t378;
t427 = 2 * m(3);
t426 = 2 * m(4);
t425 = 0.2e1 * m(5);
t424 = 2 * m(6);
t32 = -t129 * t253 - t251 * t308;
t33 = -t130 * t253 - t442;
t287 = t312 * t259;
t75 = -t251 * t436 + t253 * t287;
t331 = t137 * t258 + t75;
t79 = -t251 * t438 + t315 * t371;
t335 = -t133 * t258 + t79;
t390 = t130 * t259;
t391 = t129 * t259;
t284 = t309 * t259;
t70 = t251 * t284 + t253 * t435;
t71 = -t251 * t435 + t253 * t284;
t74 = t251 * t287 + t253 * t436;
t78 = t253 * t438 + t259 * t290;
t3 = (t253 * t71 + (-t33 + t283) * t259) * t253 + (t32 * t259 + (t134 * t373 + t138 * t378 + t250 * t74 - t252 * t78 - t390) * t251 + (-t391 + t70 + (-t138 * t259 - t335) * t252 + (t134 * t259 + t331) * t250) * t253) * t251;
t77 = -t251 * t437 + t313 * t371;
t329 = t139 * t258 + t77;
t291 = t316 * t259;
t81 = -t251 * t439 + t253 * t291;
t333 = -t135 * t258 + t81;
t34 = -t131 * t253 - t251 * t306;
t35 = -t132 * t253 - t443;
t388 = t132 * t259;
t389 = t131 * t259;
t72 = t253 * t434 + t259 * t285;
t73 = -t251 * t434 + t310 * t371;
t76 = t253 * t437 + t259 * t288;
t80 = t251 * t291 + t253 * t439;
t4 = (t253 * t73 + (-t35 + t282) * t259) * t253 + (t34 * t259 + (t136 * t373 + t140 * t378 + t250 * t76 - t252 * t80 - t388) * t251 + (-t389 + t72 + (-t140 * t259 - t333) * t252 + (t136 * t259 + t329) * t250) * t253) * t251;
t423 = -t3 - t4;
t422 = ((-t32 - t34) * t253 + (-t33 - t35) * t251) * t376;
t421 = -t251 / 0.2e1;
t420 = -t253 / 0.2e1;
t410 = rSges(4,2) * t262;
t413 = rSges(4,1) * t264;
t207 = (-t410 + t413) * qJD(3);
t419 = m(4) * t207;
t230 = rSges(4,1) * t262 + rSges(4,2) * t264;
t418 = m(4) * t230;
t412 = rSges(5,1) * t252;
t174 = (-rSges(5,2) * t250 + t412) * t258;
t417 = m(5) * t174;
t415 = pkin(3) * t262;
t372 = t253 * t258;
t278 = t250 * t372 + t252 * t376;
t340 = t253 * t352;
t322 = pkin(3) * t340;
t342 = t252 * t372;
t414 = t322 - (t353 * t253 + (t208 - t246) * t251) * t259 - rSges(6,1) * t278 - rSges(6,2) * t342 + t440;
t411 = rSges(6,1) * t252;
t409 = rSges(6,2) * t250;
t408 = pkin(1) * qJD(1);
t198 = rSges(6,1) * t250 + rSges(6,2) * t252;
t336 = -pkin(4) * t250 - t198;
t127 = t336 * t251;
t392 = t127 * t259;
t387 = t174 * t251;
t375 = t251 * t264;
t370 = t253 * t262;
t369 = t452 * t251;
t173 = (-t409 + t411) * t258;
t367 = pkin(4) * t342 + t253 * t173;
t362 = rSges(5,2) * t344 + rSges(5,3) * t371;
t360 = rSges(5,1) * t341 + rSges(5,3) * t376;
t350 = t251 * t410;
t359 = rSges(4,3) * t371 + t259 * t350;
t232 = t253 * t413;
t358 = rSges(4,3) * t376 + t259 * t232;
t128 = pkin(4) * t379 + t253 * t198;
t356 = -pkin(2) * t371 - pkin(7) * t376;
t355 = t253 * pkin(2) + t251 * pkin(7);
t354 = t251 ^ 2 + t253 ^ 2;
t263 = sin(qJ(1));
t348 = t263 * t408;
t244 = t251 * pkin(2);
t125 = pkin(7) * t253 - t244 + t357;
t323 = t251 * t266 - t210;
t126 = t323 + t355;
t185 = t246 * t371;
t294 = -t259 * t266 - t349;
t347 = t126 * t376 + t125 * t371 + t251 * (t251 * t294 + t185 + t356);
t142 = rSges(5,1) * t377 - rSges(5,2) * t380 - rSges(5,3) * t253;
t144 = -rSges(5,1) * t374 + rSges(5,2) * t379 - t251 * rSges(5,3);
t343 = t250 * t371;
t279 = -t251 * t373 - t343;
t345 = t251 * t378;
t346 = t144 * t376 + t142 * t371 + t251 * (-rSges(5,1) * t345 + rSges(5,2) * t279 + t360);
t339 = t253 * t351;
t338 = t376 / 0.2e1;
t337 = -t371 / 0.2e1;
t201 = t253 * rSges(3,1) - rSges(3,2) * t251;
t334 = -t134 * t258 + t78;
t332 = -t136 * t258 + t80;
t330 = t138 * t258 + t74;
t328 = t140 * t258 + t76;
t176 = rSges(3,1) * t371 - rSges(3,2) * t376;
t321 = rSges(4,1) * t375 - t350;
t145 = (-t199 - t415) * t251;
t1 = (t251 * t70 + (t36 + t442) * t259) * t251 + (-t37 * t259 + (-t133 * t373 - t137 * t378 - t250 * t75 + t252 * t79 + t391) * t253 + (t390 + t71 + (-t137 * t259 + t334) * t252 + (t133 * t259 - t330) * t250) * t251) * t253;
t2 = (t251 * t72 + (t38 + t443) * t259) * t251 + (-t39 * t259 + (-t135 * t373 - t139 * t378 - t250 * t77 + t252 * t81 + t389) * t253 + (t388 + t73 + (-t139 * t259 + t332) * t252 + (t135 * t259 - t328) * t250) * t251) * t253;
t318 = (-t2 - t1) * t251 + t422;
t311 = Icges(4,5) * t264 - Icges(4,6) * t262;
t289 = t314 * t251;
t156 = -Icges(4,6) * t253 + t289;
t292 = t317 * t251;
t158 = -Icges(4,5) * t253 + t292;
t304 = t156 * t264 + t158 * t262;
t303 = t156 * t262 - t158 * t264;
t302 = t157 * t264 + t159 * t262;
t297 = t218 * t262 - t219 * t264;
t161 = rSges(4,2) * t370 - t251 * rSges(4,3) - t232;
t296 = t336 - t415;
t295 = t451 * t376 + t452 * t371 + (-rSges(6,1) * t345 + rSges(6,2) * t279 - qJD(5) * t253 - t185 + (-t259 * t353 + t349) * t251 + t447) * t251;
t175 = t441 * t259;
t293 = qJD(3) * t230;
t286 = t311 * t251;
t281 = t303 * t253;
t280 = t296 * t259;
t277 = -t251 * t351 - t259 * t370;
t121 = -t161 + t355;
t273 = -t311 * qJD(3) - t259 * t297;
t106 = t142 + t357;
t107 = -t144 - t323;
t100 = -t251 * t257 + t453;
t120 = t244 + (-rSges(4,3) - pkin(7)) * t253 + t321;
t272 = -t198 * t258 - t257 * t259;
t270 = -t199 * t258 + t294;
t269 = (t446 * t253 + (t328 + t330) * t252 + t449 * t251 + (t332 + t334) * t250) * t421 + (t449 * t253 + (t329 + t331) * t252 - t446 * t251 + (t333 + t335) * t250) * t420 + (t450 * t253 + (t133 + t135) * t252 - t448 * t251 + (t137 + t139) * t250) * t338 + (t448 * t253 + (t134 + t136) * t252 + t450 * t251 + (t138 + t140) * t250) * t337;
t65 = -rSges(4,1) * t251 * t352 + rSges(4,2) * t277 - t356 + t358;
t228 = pkin(7) * t371;
t64 = -rSges(4,2) * t339 - pkin(2) * t376 + t228 + (-t259 * t375 - t340) * rSges(4,1) + t359;
t267 = t269 + (-qJD(3) * t301 + t273 * t251 + t253 * t430 + t262 * (t253 * t431 + t259 * t292) + t264 * (t253 * t432 + t259 * t289)) * t421 + (-qJD(3) * t303 - t251 * t430 + t273 * t253 + t262 * (-t251 * t431 + t317 * t371) + t264 * (-t251 * t432 + t314 * t371)) * t420 + (-t217 * t253 - t251 * t297 + t304) * t338 + (-t217 * t251 + t253 * t297 + t302) * t337;
t30 = (-t208 - t411) * t376 + t272 * t253 + t440;
t31 = (-t259 * t409 - qJD(5)) * t253 + t272 * t251 + t447;
t43 = -rSges(5,2) * t343 + t251 * t270 + t185 + t360;
t42 = (-t246 - t412) * t376 + t270 * t253 + t362;
t265 = cos(qJ(1));
t256 = t265 * pkin(1);
t254 = t263 * pkin(1);
t247 = t265 * t408;
t234 = pkin(3) * t370;
t216 = pkin(3) * t339;
t179 = t201 + t256;
t178 = t254 - t441;
t160 = -rSges(4,3) * t253 + t321;
t155 = -Icges(4,3) * t251 - t253 * t311;
t154 = -Icges(4,3) * t253 + t286;
t148 = t176 + t247;
t147 = t175 - t348;
t146 = t199 * t253 + t234;
t124 = t251 * t142;
t122 = t251 * t125;
t115 = t234 + t128;
t114 = t296 * t251;
t111 = t121 + t256;
t110 = t254 + t120;
t104 = t107 + t256;
t103 = t106 + t254;
t98 = t322 + t228 + (t235 + (-pkin(2) + t246) * t251) * t259;
t93 = -t251 * t433 + t311 * t371;
t92 = t253 * t433 + t259 * t286;
t91 = t100 + t256;
t90 = t254 + t99;
t86 = rSges(5,1) * t278 + rSges(5,2) * t342 - t362;
t82 = -t144 * t253 + t124;
t69 = pkin(3) * t277 - t199 * t371 - t387;
t68 = t145 * t259 + t174 * t253 + t216;
t59 = pkin(4) * t279 - t173 * t251 - t198 * t371;
t58 = t367 + t392;
t53 = t247 + t65;
t52 = t64 - t348;
t50 = -t155 * t251 + t301 * t253;
t49 = -t154 * t251 + t281;
t48 = -t155 * t253 - t444;
t47 = -t154 * t253 - t303 * t251;
t46 = t253 * t280 + (-pkin(3) * t351 - pkin(4) * t373 - t173) * t251;
t45 = t251 * t280 + t216 + t367;
t41 = t247 + t43;
t40 = t42 - t348;
t29 = t247 + t31;
t28 = t30 - t348;
t27 = t122 + t124 + (-t126 - t144) * t253;
t24 = -t253 * t451 + t369;
t21 = -t253 * t86 + t346;
t10 = t122 + (-t126 - t451) * t253 + t369;
t7 = (-t86 - t98) * t253 + t346 + t347;
t6 = t253 * t414 + t295;
t5 = (-t98 + t414) * t253 + t295 + t347;
t8 = [(t147 * t179 + t148 * t178) * t427 + (t110 * t53 + t111 * t52) * t426 + (t103 * t41 + t104 * t40) * t425 + (t28 * t91 + t29 * t90) * t424 + t268; m(3) * (t147 * t201 - t148 * t441 + t175 * t179 + t176 * t178) + m(4) * (t110 * t65 + t111 * t64 + t120 * t53 + t121 * t52) + m(5) * (t103 * t43 + t104 * t42 + t106 * t41 + t107 * t40) + m(6) * (t100 * t28 + t29 * t99 + t30 * t91 + t31 * t90) + t268; (t100 * t30 + t31 * t99) * t424 + (t106 * t43 + t107 * t42) * t425 + (t120 * t65 + t121 * t64) * t426 + (t175 * t201 - t176 * t441) * t427 + t268; t267 + m(6) * (t114 * t28 + t115 * t29 + t45 * t90 + t46 * t91) + m(5) * (t103 * t68 + t104 * t69 + t145 * t40 + t146 * t41) + ((-t111 * t259 + t53) * t253 + (-t110 * t259 - t52) * t251) * t418 + (t110 * t253 - t111 * t251) * t419; m(6) * (t100 * t46 + t114 * t30 + t115 * t31 + t45 * t99) + m(5) * (t106 * t68 + t107 * t69 + t145 * t42 + t146 * t43) + t267 + ((-t121 * t259 + t65) * t253 + (-t120 * t259 - t64) * t251) * t418 + (t120 * t253 - t121 * t251) * t419; (t10 * t5 + t114 * t46 + t115 * t45) * t424 - t251 * t2 - t251 * t1 - t253 * t3 + (t145 * t69 + t146 * t68 + t27 * t7) * t425 - t253 * t4 - t251 * ((t251 * t92 + (t49 + t444) * t259) * t251 + (-t50 * t259 + (-t156 * t351 - t158 * t352) * t253 + (-t302 * qJD(3) + t259 * t303 + t93) * t251) * t253) + (-t251 * t48 - t253 * t47) * t376 - t253 * ((t253 * t93 + (-t48 + t281) * t259) * t253 + (t47 * t259 + (t157 * t351 + t159 * t352) * t251 + (t304 * qJD(3) + t259 * t301 + t92) * t253) * t251) + ((t160 * t251 - t161 * t253) * ((t259 * t160 - t253 * t293 + t359) * t253 + (-t251 * t293 + (t161 + (-t410 - t413) * t253) * t259 + t358) * t251) + t354 * t230 * t207) * t426 + t422 + (t251 * t50 + t253 * t49 + t445) * t371; m(6) * (t127 * t28 + t128 * t29 + t58 * t90 + t59 * t91) + t269 + ((-t104 * t259 + t41) * t253 + (-t103 * t259 - t40) * t251) * t416 + (t103 * t253 - t104 * t251) * t417; m(6) * (t100 * t59 + t127 * t30 + t128 * t31 + t58 * t99) + t269 + ((-t107 * t259 + t43) * t253 + (-t106 * t259 - t42) * t251) * t416 + (t106 * t253 - t107 * t251) * t417; m(6) * (t10 * t6 + t114 * t59 + t115 * t58 + t127 * t46 + t128 * t45 + t24 * t5) + m(5) * (-t145 * t387 + t21 * t27 + t7 * t82 + (-t146 * t376 - t251 * t69) * t199) + (m(5) * (t146 * t174 + t199 * t68) + t423 + (-t145 * t416 + t445) * t259) * t253 + t318; (t174 * t199 * t354 + t21 * t82) * t425 + (t127 * t59 + t128 * t58 + t24 * t6) * t424 + (t445 * t259 + t423) * t253 + t318; m(6) * ((-t259 * t90 - t28) * t253 + (t259 * t91 - t29) * t251); m(6) * ((-t259 * t99 - t30) * t253 + (t100 * t259 - t31) * t251); m(6) * ((-t115 * t259 - t46) * t253 + (t114 * t259 - t45) * t251); m(6) * ((-t128 * t259 - t59) * t253 + (-t58 + t392) * t251); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;

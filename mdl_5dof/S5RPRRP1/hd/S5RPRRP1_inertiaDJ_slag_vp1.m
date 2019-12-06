% Calculate time derivative of joint inertia matrix for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:15
% EndTime: 2019-12-05 17:59:36
% DurationCPUTime: 10.21s
% Computational Cost: add. (5749->412), mult. (7956->559), div. (0->0), fcn. (6126->6), ass. (0->234)
t434 = -Icges(5,4) - Icges(6,4);
t433 = Icges(5,1) + Icges(6,1);
t432 = Icges(5,2) + Icges(6,2);
t212 = qJ(3) + qJ(4);
t201 = cos(t212);
t431 = t434 * t201;
t200 = sin(t212);
t430 = t434 * t200;
t429 = Icges(5,5) + Icges(6,5);
t428 = Icges(5,6) + Icges(6,6);
t427 = t201 * t432 - t430;
t426 = t200 * t433 - t431;
t214 = sin(qJ(1));
t216 = cos(qJ(1));
t394 = t214 * t429 - t216 * t426;
t396 = -t214 * t428 + t216 * t427;
t425 = Icges(5,3) + Icges(6,3);
t397 = t214 * t427 + t216 * t428;
t395 = t214 * t426 + t216 * t429;
t424 = t200 * t429 + t201 * t428;
t420 = t200 * t395 + t201 * t397;
t423 = t420 * t214;
t415 = t200 * t432 + t431;
t422 = t201 * t433 + t430;
t421 = t200 * t394 - t201 * t396;
t418 = t214 * t424 + t216 * t425;
t411 = -t214 * t425 + t216 * t424;
t213 = sin(qJ(3));
t354 = pkin(3) * t213;
t162 = pkin(4) * t200 + t354;
t254 = rSges(6,1) * t200 + rSges(6,2) * t201;
t231 = t162 + t254;
t419 = t216 * t231;
t209 = qJD(3) + qJD(4);
t417 = t427 * t209;
t416 = t426 * t209;
t393 = -t200 * t428 + t201 * t429;
t413 = t421 * t214;
t391 = t200 * t422 - t201 * t415;
t410 = -t216 * t418 - t423;
t409 = t216 * t411 - t413;
t317 = t209 * t216;
t408 = qJD(1) * t397 + t317 * t415;
t318 = t209 * t214;
t407 = -qJD(1) * t396 + t318 * t415;
t406 = qJD(1) * t395 - t317 * t422;
t405 = qJD(1) * t394 - t318 * t422;
t404 = t418 * qJD(1);
t403 = t420 * t216;
t402 = t214 * t418 - t403;
t401 = -t214 * t411 - t216 * t421;
t400 = -t317 * t393 + t404;
t399 = qJD(1) * t411 + t318 * t393;
t217 = -pkin(7) - pkin(6);
t208 = -qJ(5) + t217;
t310 = t214 * t217 + t216 * t354;
t398 = t310 - t419 + (rSges(6,3) - t208) * t214;
t392 = qJD(1) * t391 - t209 * t424;
t390 = -t209 * t394 - t408;
t389 = -t209 * t395 + t407;
t388 = t209 * t396 + t406;
t387 = -t209 * t397 - t405;
t215 = cos(qJ(3));
t297 = qJD(3) * t215;
t300 = qJD(1) * t216;
t386 = t213 * t300 + t214 * t297;
t385 = t200 * t300 + t201 * t318;
t384 = (-t209 * t422 + t417) * t201 + (-t209 * t415 + t416) * t200 + t393 * qJD(1);
t383 = t214 * t409 + t216 * t410;
t154 = rSges(5,1) * t201 - rSges(5,2) * t200;
t382 = t154 * t317;
t342 = Icges(4,4) * t213;
t246 = Icges(4,2) * t215 + t342;
t130 = Icges(4,6) * t216 + t214 * t246;
t341 = Icges(4,4) * t215;
t250 = Icges(4,1) * t213 + t341;
t132 = Icges(4,5) * t216 + t214 * t250;
t235 = t130 * t215 + t132 * t213;
t381 = t216 * t235;
t255 = rSges(5,1) * t200 + rSges(5,2) * t201;
t378 = t216 * t255;
t256 = rSges(4,1) * t213 + rSges(4,2) * t215;
t229 = t216 * t256;
t204 = t216 * rSges(6,3);
t319 = t201 * t214;
t321 = t200 * t214;
t377 = -rSges(6,1) * t321 - rSges(6,2) * t319 - t214 * t162 - t204;
t320 = t201 * t209;
t144 = pkin(3) * t297 + pkin(4) * t320;
t279 = t201 * t300;
t301 = qJD(1) * t214;
t376 = -rSges(6,1) * t385 - rSges(6,2) * t279 - qJD(5) * t216 - t144 * t214 - t162 * t300 - t208 * t301;
t277 = t215 * t300;
t282 = rSges(4,1) * t386 + rSges(4,2) * t277;
t294 = -rSges(4,3) - pkin(1) - pkin(6);
t298 = qJD(3) * t213;
t309 = qJ(2) * t300 + qJD(2) * t214;
t54 = (-rSges(4,2) * t298 + qJD(1) * t294) * t214 + t282 + t309;
t348 = rSges(4,2) * t213;
t172 = rSges(4,1) * t215 - t348;
t199 = qJD(2) * t216;
t296 = qJD(3) * t216;
t55 = t199 + t172 * t296 + (t294 * t216 + (-qJ(2) - t256) * t214) * qJD(1);
t375 = t214 * t55 - t216 * t54;
t242 = Icges(4,5) * t213 + Icges(4,6) * t215;
t372 = -Icges(4,3) * t214 + t216 * t242;
t369 = -Icges(4,6) * t214 + t216 * t246;
t366 = -Icges(4,5) * t214 + t216 * t250;
t363 = 2 * m(4);
t362 = 2 * m(5);
t361 = 2 * m(6);
t210 = t214 ^ 2;
t211 = t216 ^ 2;
t360 = t214 / 0.2e1;
t359 = t216 / 0.2e1;
t358 = rSges(3,2) - pkin(1);
t357 = -rSges(5,3) - pkin(1);
t356 = -rSges(6,3) - pkin(1);
t355 = m(4) * t172;
t353 = pkin(3) * t215;
t352 = t216 * pkin(6);
t153 = rSges(6,1) * t201 - rSges(6,2) * t200;
t295 = qJD(5) * t214;
t288 = pkin(3) * t296;
t299 = qJD(1) * t217;
t311 = t215 * t288 + t216 * t299;
t314 = t216 * t208;
t351 = (-t144 * t216 - t153 * t317 + t295 + t311 + (-t314 + t204 + (-t354 + t231) * t214) * qJD(1)) * t216;
t281 = pkin(3) * t386 + t214 * t299;
t322 = t200 * t209;
t290 = rSges(6,2) * t322;
t350 = t281 - (-rSges(6,3) * qJD(1) - t290) * t214 + t376;
t349 = t398 * t216;
t347 = t214 * rSges(4,3);
t206 = t216 * rSges(4,3);
t205 = t216 * rSges(5,3);
t316 = t213 * t214;
t192 = pkin(3) * t316;
t344 = t192 - (-t208 + t217) * t216 + t377;
t127 = t255 * t209;
t327 = t127 * t216;
t326 = t130 * t213;
t325 = t213 * t369;
t324 = t132 * t215;
t323 = t215 * t366;
t315 = t214 * t215;
t116 = rSges(5,1) * t321 + rSges(5,2) * t319 + t205;
t142 = t192 + (-pkin(6) - t217) * t216;
t313 = -t116 - t142;
t177 = pkin(4) * t319;
t98 = t153 * t214 + t177;
t193 = pkin(3) * t315;
t312 = qJD(1) * t193 + t213 * t288;
t308 = pkin(1) * t216 + t214 * qJ(2);
t307 = t210 + t211;
t128 = Icges(4,3) * t216 + t214 * t242;
t304 = qJD(1) * t128;
t293 = t208 + t356;
t292 = pkin(4) * t322;
t291 = rSges(5,2) * t322;
t289 = pkin(3) * t298;
t287 = -t142 + t344;
t126 = t254 * t209;
t285 = pkin(4) * t279 - t126 * t214 + t153 * t300;
t283 = rSges(5,1) * t385 + rSges(5,2) * t279;
t135 = rSges(4,1) * t316 + rSges(4,2) * t315 + t206;
t273 = -pkin(4) * t201 - t153;
t264 = t307 * t127;
t160 = t256 * qJD(3);
t263 = t307 * t160;
t257 = ((t399 * t216 + (t403 - t409) * qJD(1)) * t216 + t402 * t300) * t216 + ((t400 * t214 + (-t402 + t413) * qJD(1)) * t214 + t401 * t300 + (t399 * t214 + (t400 + t404) * t216 + (t214 * t396 + t216 * t397) * t322 + (t214 * t394 - t216 * t395) * t320 + ((t390 + t408) * t214 + (-t389 + t407) * t216) * t201 + ((-t388 + t406) * t214 + (t387 + t405) * t216) * t200 + (t423 + (t421 - t418) * t216 + t401 + t410) * qJD(1)) * t216) * t214;
t251 = Icges(4,1) * t215 - t342;
t247 = -Icges(4,2) * t213 + t341;
t243 = Icges(4,5) * t215 - Icges(4,6) * t213;
t234 = -t213 * t366 - t215 * t369;
t48 = qJD(1) * t177 + t153 * t301 + (t126 + t292) * t216;
t230 = rSges(3,3) * t216 + t214 * t358;
t226 = t234 * t214;
t225 = qJD(3) * t251;
t224 = qJD(3) * t247;
t221 = t301 * t383 + t257;
t26 = (qJD(1) * t357 - t291) * t214 + t281 + t283 + t309;
t27 = t199 + t382 + (t357 * t216 + (-qJ(2) - t255 - t354) * t214) * qJD(1) + t311;
t203 = t216 * qJ(2);
t77 = t214 * t357 + t203 + t310 + t378;
t78 = -t216 * t217 + t116 + t192 + t308;
t220 = m(5) * (t214 * t27 - t216 * t26 + (t214 * t78 + t216 * t77) * qJD(1));
t219 = (t200 * t390 + t201 * t388 + t214 * t392 + t216 * t384) * t360 + (t200 * t389 + t201 * t387 - t214 * t384 + t216 * t392) * t359 - (-t200 * t397 + t201 * t395 + t214 * t391 + t216 * t393) * t301 / 0.2e1 + (t200 * t396 + t201 * t394 + t214 * t393 - t216 * t391) * t300 / 0.2e1;
t100 = t154 * t214 + t193;
t101 = (-t154 - t353) * t216;
t56 = t154 * t301 + t312 + t327;
t187 = pkin(3) * t277;
t57 = t154 * t300 + t187 + (-t127 - t289) * t214;
t218 = m(5) * (t57 * t214 - t56 * t216 + (t100 * t216 + t101 * t214) * qJD(1));
t141 = -pkin(6) * t214 - t310;
t138 = -rSges(3,2) * t216 + rSges(3,3) * t214 + t308;
t137 = t203 + t230;
t136 = t347 - t229;
t119 = t216 * t141;
t118 = rSges(5,3) * t214 - t378;
t99 = t273 * t216;
t97 = t216 * t118;
t95 = t199 + (t358 * t216 + (-rSges(3,3) - qJ(2)) * t214) * qJD(1);
t94 = qJD(1) * t230 + t309;
t93 = pkin(6) * t301 + t281;
t92 = t135 + t308 + t352;
t91 = t214 * t294 + t203 + t229;
t90 = (t273 - t353) * t216;
t89 = t193 + t98;
t88 = t216 * ((t192 - t352) * qJD(1) - t311);
t80 = qJD(3) * t214 * t243 + qJD(1) * t372;
t79 = -t243 * t296 + t304;
t76 = (-rSges(5,3) * qJD(1) - t291) * t214 + t283;
t62 = t308 - t314 - t377;
t61 = t214 * t293 + t203 + t419;
t60 = -t116 * t214 + t97;
t59 = t216 * (-t382 + (t214 * t255 + t205) * qJD(1));
t49 = -t214 * t292 + t285;
t42 = -t214 * t372 - t216 * t234;
t41 = t128 * t214 - t381;
t40 = -t216 * t372 + t226;
t39 = t128 * t216 + t214 * t235;
t37 = t187 + (-t289 - t292) * t214 + t285;
t36 = t48 + t312;
t25 = t214 * t313 + t119 + t97;
t24 = -t295 + t199 + (t153 * t209 + t144) * t216 + (t293 * t216 + (-qJ(2) - t231) * t214) * qJD(1);
t23 = (qJD(1) * t356 - t290) * t214 + t309 - t376;
t22 = t214 * t344 + t349;
t21 = -t214 * t76 + t59 + (-t116 * t216 - t118 * t214) * qJD(1);
t20 = t214 * t287 + t119 + t349;
t7 = t59 + t88 + (-t76 - t93) * t214 + (t313 * t216 + (-t118 - t141) * t214) * qJD(1);
t6 = t350 * t214 + (-t214 * t398 + t216 * t344) * qJD(1) + t351;
t5 = t88 + (-t93 + t350) * t214 + (t287 * t216 + (-t141 - t398) * t214) * qJD(1) + t351;
t1 = [(t23 * t62 + t24 * t61) * t361 + (t26 * t78 + t27 * t77) * t362 + (t54 * t92 + t55 * t91) * t363 - t213 * t225 - t250 * t297 - t215 * t224 + t246 * t298 + 0.2e1 * m(3) * (t137 * t95 + t138 * t94) - t422 * t322 + t415 * t320 - t416 * t201 + t417 * t200; m(6) * (t214 * t24 - t216 * t23 + (t214 * t62 + t216 * t61) * qJD(1)) + t220 + m(4) * ((t214 * t92 + t216 * t91) * qJD(1) + t375) + m(3) * (t214 * t95 - t216 * t94 + (t137 * t216 + t138 * t214) * qJD(1)); 0; t219 + ((t326 / 0.2e1 - t324 / 0.2e1 + t92 * t355) * t214 + (t91 * t355 + t325 / 0.2e1 - t323 / 0.2e1) * t216) * qJD(1) + m(4) * (t375 * t172 - (t214 * t91 - t216 * t92) * t160) + (-qJD(3) * t234 - t213 * (qJD(1) * t130 - t247 * t296) + t215 * (qJD(1) * t132 - t251 * t296)) * t360 + (-qJD(3) * t235 - t213 * (qJD(1) * t369 + t214 * t224) + t215 * (qJD(1) * t366 + t214 * t225)) * t359 + m(5) * (t100 * t27 + t101 * t26 + t56 * t78 + t57 * t77) + m(6) * (t23 * t90 + t24 * t89 + t36 * t62 + t37 * t61) - (t210 / 0.2e1 + t211 / 0.2e1) * t242 * qJD(3); t218 + m(6) * (t214 * t37 - t216 * t36 + (t214 * t90 + t216 * t89) * qJD(1)) - m(4) * t263; (t20 * t5 + t36 * t90 + t37 * t89) * t361 + (t100 * t57 + t101 * t56 + t25 * t7) * t362 + ((-t135 * t214 + t136 * t216) * (-t214 * t282 + (-t172 * t211 + t210 * t348) * qJD(3) + ((-t135 + t206) * t216 + (-t136 + t229 + t347) * t214) * qJD(1)) - t172 * t263) * t363 + t216 * ((t216 * t80 + (t40 + t381) * qJD(1)) * t216 + (-t39 * qJD(1) + (-t297 * t366 + t298 * t369) * t214 + (t79 + (t324 - t326) * qJD(3) + (-t128 + t234) * qJD(1)) * t216) * t214) + (t214 * t42 + t216 * t41) * t300 + t214 * ((t214 * t79 + (-t41 + t226) * qJD(1)) * t214 + (t42 * qJD(1) + (t130 * t298 - t132 * t297 + t304) * t216 + (t80 + (t323 - t325) * qJD(3) + t235 * qJD(1)) * t214) * t216) + t257 + (-t214 * t40 - t216 * t39 + t383) * t301; t219 + m(6) * (t23 * t99 + t24 * t98 + t48 * t62 + t49 * t61) + t154 * t220 - m(5) * (t214 * t77 - t216 * t78) * t127; m(6) * (t49 * t214 - t48 * t216 + (t214 * t99 + t216 * t98) * qJD(1)) - m(5) * t264; m(6) * (t20 * t6 + t22 * t5 + t36 * t99 + t37 * t98 + t48 * t90 + t49 * t89) + m(5) * (-t100 * t127 * t214 + t101 * t327 + t21 * t25 + t60 * t7) + t154 * t218 + t221; (-t154 * t264 + t21 * t60) * t362 + (t22 * t6 + t48 * t99 + t49 * t98) * t361 + t221; m(6) * (t214 * t23 + t216 * t24 + (-t214 * t61 + t216 * t62) * qJD(1)); 0; m(6) * (t214 * t36 + t216 * t37 + (-t214 * t89 + t216 * t90) * qJD(1)); m(6) * (t214 * t48 + t216 * t49 + (-t214 * t98 + t216 * t99) * qJD(1)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

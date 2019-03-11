% Calculate time derivative of joint inertia matrix for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR3_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:08
% EndTime: 2019-03-09 02:44:25
% DurationCPUTime: 13.93s
% Computational Cost: add. (13382->651), mult. (17992->888), div. (0->0), fcn. (16226->8), ass. (0->342)
t497 = Icges(5,4) + Icges(4,5);
t496 = -Icges(4,6) + Icges(5,6);
t248 = sin(qJ(3));
t251 = cos(qJ(3));
t306 = Icges(6,5) * t248 - Icges(6,6) * t251;
t493 = t496 * t248 + t251 * t497;
t495 = t306 - t493;
t494 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t421 = Icges(6,4) * t248;
t424 = Icges(4,4) * t248;
t492 = t421 + t424 + (Icges(6,1) + Icges(4,2)) * t251;
t415 = Icges(5,5) * t251;
t420 = Icges(6,4) * t251;
t423 = Icges(4,4) * t251;
t491 = -t415 + t420 + t423 + (Icges(4,1) + Icges(5,1) + Icges(6,2)) * t248;
t246 = qJ(1) + pkin(9);
t243 = sin(t246);
t244 = cos(t246);
t483 = t494 * t243 - t244 * t495;
t315 = -Icges(4,2) * t248 + t423;
t146 = Icges(4,6) * t243 + t244 * t315;
t322 = Icges(4,1) * t251 - t424;
t152 = Icges(4,5) * t243 + t244 * t322;
t293 = t146 * t248 - t152 * t251;
t311 = -Icges(6,2) * t251 + t421;
t142 = -Icges(6,6) * t243 + t244 * t311;
t317 = Icges(6,1) * t248 - t420;
t148 = -Icges(6,5) * t243 + t244 * t317;
t295 = t142 * t251 - t148 * t248;
t308 = Icges(5,3) * t248 + t415;
t138 = Icges(5,6) * t243 + t244 * t308;
t416 = Icges(5,5) * t248;
t320 = Icges(5,1) * t251 + t416;
t150 = Icges(5,4) * t243 + t244 * t320;
t297 = t138 * t248 + t150 * t251;
t455 = t293 + t295 - t297;
t490 = t455 * t244;
t137 = -Icges(5,6) * t244 + t243 * t308;
t147 = Icges(6,5) * t244 + t243 * t317;
t489 = t137 + t147;
t488 = t138 + t148;
t149 = -Icges(5,4) * t244 + t243 * t320;
t151 = -Icges(4,5) * t244 + t243 * t322;
t487 = t149 + t151;
t486 = t150 + t152;
t485 = t492 * qJD(3);
t484 = t243 * t495 + t494 * t244;
t482 = ((Icges(6,5) - t496) * t251 + (Icges(6,6) + t497) * t248) * qJD(3);
t145 = -Icges(4,6) * t244 + t243 * t315;
t294 = t145 * t248 - t151 * t251;
t141 = Icges(6,6) * t244 + t243 * t311;
t296 = t141 * t251 - t147 * t248;
t298 = t137 * t248 + t149 * t251;
t480 = t294 + t296 - t298;
t479 = -t243 / 0.2e1;
t478 = t243 / 0.2e1;
t477 = -t244 / 0.2e1;
t445 = t244 / 0.2e1;
t476 = -qJD(1) / 0.2e1;
t439 = qJD(1) / 0.2e1;
t375 = qJD(3) * t248;
t355 = t244 * t375;
t378 = qJD(1) * t251;
t358 = t243 * t378;
t475 = t355 + t358;
t247 = sin(qJ(6));
t250 = cos(qJ(6));
t418 = Icges(7,4) * t250;
t310 = -Icges(7,2) * t247 + t418;
t173 = Icges(7,6) * t248 - t251 * t310;
t419 = Icges(7,4) * t247;
t316 = Icges(7,1) * t250 - t419;
t174 = Icges(7,5) * t248 - t251 * t316;
t474 = -t173 * t247 + t174 * t250;
t274 = t293 * t243;
t275 = t294 * t244;
t276 = t295 * t243;
t277 = t296 * t244;
t278 = t297 * t243;
t279 = t298 * t244;
t473 = -t243 * t484 - t244 * t483 - t274 - t275 - t276 - t277 + t278 + t279;
t472 = t483 * qJD(1);
t471 = (t142 - t486) * t243 - (-t141 + t487) * t244;
t470 = (t146 - t488) * t243 + (t145 - t489) * t244;
t241 = t243 ^ 2;
t242 = t244 ^ 2;
t469 = qJD(3) / 0.2e1;
t435 = rSges(4,1) * t251;
t340 = -rSges(4,2) * t248 + t435;
t431 = rSges(4,3) * t244;
t156 = t243 * t340 - t431;
t397 = t244 * t251;
t399 = t244 * t248;
t467 = -rSges(4,2) * t399 + t243 * rSges(4,3);
t159 = rSges(4,1) * t397 + t467;
t219 = rSges(4,1) * t248 + rSges(4,2) * t251;
t280 = qJD(3) * t219;
t379 = qJD(1) * t248;
t359 = t243 * t379;
t380 = qJD(1) * t244;
t256 = rSges(4,2) * t359 + rSges(4,3) * t380 - t244 * t280;
t31 = (qJD(1) * t156 + t256) * t244 + (-t243 * t280 + (-t159 + t467) * qJD(1)) * t243;
t453 = 2 * m(4);
t468 = t31 * t453;
t245 = cos(qJ(1)) * pkin(1);
t466 = t243 * pkin(7) + t245;
t381 = qJD(1) * t243;
t465 = t243 * t480 - t244 * t484;
t462 = t243 * t483 - t490;
t461 = qJD(1) * t484 - t244 * t482;
t460 = t243 * t482 - t472;
t452 = 2 * m(5);
t451 = 2 * m(6);
t450 = 2 * m(7);
t449 = m(5) / 0.2e1;
t448 = m(6) / 0.2e1;
t447 = m(7) / 0.2e1;
t446 = -pkin(3) - pkin(4);
t444 = t248 / 0.2e1;
t443 = -rSges(5,1) - pkin(3);
t442 = rSges(7,3) + pkin(8);
t441 = m(4) * t219;
t440 = sin(qJ(1)) * pkin(1);
t438 = -pkin(5) - qJ(4);
t305 = Icges(7,5) * t250 - Icges(7,6) * t247;
t371 = qJD(6) * t251;
t123 = (Icges(7,5) * t247 + Icges(7,6) * t250) * t371 + (Icges(7,3) * t251 + t248 * t305) * qJD(3);
t124 = (Icges(7,2) * t250 + t419) * t371 + (Icges(7,6) * t251 + t248 * t310) * qJD(3);
t125 = (Icges(7,1) * t247 + t418) * t371 + (Icges(7,5) * t251 + t248 * t316) * qJD(3);
t172 = Icges(7,3) * t248 - t251 * t305;
t374 = qJD(3) * t251;
t254 = t248 * t123 + t172 * t374 + t474 * t375 + (-t125 * t251 + t173 * t371) * t250 + (t124 * t251 + t174 * t371) * t247;
t74 = t172 * t248 - t251 * t474;
t437 = t254 * t248 + t74 * t374;
t351 = -qJD(6) * t248 - qJD(1);
t259 = -t247 * t374 + t250 * t351;
t350 = qJD(6) + t379;
t403 = t243 * t247;
t87 = t244 * t259 + t350 * t403;
t258 = t247 * t351 + t250 * t374;
t402 = t243 * t250;
t88 = t244 * t258 - t350 * t402;
t436 = t88 * rSges(7,1) + t87 * rSges(7,2);
t434 = rSges(5,1) * t248;
t433 = rSges(5,2) * t244;
t432 = rSges(6,2) * t248;
t430 = rSges(6,3) * t244;
t237 = t243 * rSges(5,2);
t429 = -rSges(6,1) - qJ(4);
t428 = -rSges(5,3) - qJ(4);
t427 = -rSges(6,3) - qJ(5);
t344 = pkin(5) * t248 + pkin(8) * t251;
t396 = t247 * t248;
t398 = t244 * t250;
t168 = -t243 * t396 + t398;
t395 = t248 * t250;
t400 = t244 * t247;
t169 = t243 * t395 + t400;
t336 = -rSges(7,1) * t169 - rSges(7,2) * t168;
t401 = t243 * t251;
t97 = rSges(7,3) * t401 - t336;
t426 = t243 * t344 + t97;
t170 = -t244 * t396 - t402;
t171 = t244 * t395 - t403;
t98 = t171 * rSges(7,1) + t170 * rSges(7,2) + rSges(7,3) * t397;
t425 = pkin(5) * t399 + pkin(8) * t397 + t98;
t407 = qJ(5) * t243;
t406 = qJ(5) * t244;
t334 = pkin(3) * t251 + qJ(4) * t248;
t175 = t334 * t243;
t176 = pkin(3) * t397 + qJ(4) * t399;
t394 = t243 * t175 + t244 * t176;
t230 = pkin(4) * t397;
t182 = t230 - t407;
t393 = -t176 - t182;
t335 = rSges(7,1) * t250 - rSges(7,2) * t247;
t177 = rSges(7,3) * t248 - t251 * t335;
t345 = pkin(5) * t251 - pkin(8) * t248;
t392 = t177 - t345;
t180 = qJD(3) * t334 - qJD(4) * t251;
t339 = rSges(5,1) * t251 + rSges(5,3) * t248;
t391 = -t339 * qJD(3) - t180;
t217 = pkin(3) * t248 - qJ(4) * t251;
t183 = t217 * t381;
t390 = pkin(4) * t359 + t183;
t354 = t244 * t374;
t373 = qJD(4) * t248;
t389 = qJ(4) * t354 + t244 * t373;
t388 = rSges(5,2) * t380 + rSges(5,3) * t354;
t356 = t243 * t375;
t387 = -pkin(4) * t356 - qJ(5) * t381;
t218 = -rSges(5,3) * t251 + t434;
t386 = -t217 - t218;
t385 = t241 + t242;
t377 = qJD(3) * t243;
t376 = qJD(3) * t244;
t372 = qJD(5) * t243;
t204 = pkin(3) * t356;
t357 = t244 * t378;
t370 = t175 * t380 + t243 * (pkin(3) * t357 + t243 * t373 - t204 + (t243 * t374 + t244 * t379) * qJ(4)) + t244 * (-pkin(3) * t475 - qJ(4) * t359 + t389);
t368 = rSges(6,2) * t397;
t94 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t397;
t96 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t397;
t324 = t247 * t94 - t250 * t96;
t39 = Icges(7,5) * t88 + Icges(7,6) * t87 - Icges(7,3) * t475;
t41 = Icges(7,4) * t88 + Icges(7,2) * t87 - Icges(7,6) * t475;
t43 = Icges(7,1) * t88 + Icges(7,4) * t87 - Icges(7,5) * t475;
t92 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t397;
t11 = (-qJD(3) * t324 + t39) * t248 + (qJD(3) * t92 + t247 * t41 - t250 * t43 + (t247 * t96 + t250 * t94) * qJD(6)) * t251;
t18 = t123 * t397 + t124 * t170 + t125 * t171 - t172 * t475 + t173 * t87 + t174 * t88;
t367 = t18 / 0.2e1 + t11 / 0.2e1;
t93 = Icges(7,4) * t169 + Icges(7,2) * t168 + Icges(7,6) * t401;
t95 = Icges(7,1) * t169 + Icges(7,4) * t168 + Icges(7,5) * t401;
t325 = t247 * t93 - t250 * t95;
t263 = -t356 + t357;
t89 = t243 * t259 - t350 * t400;
t90 = t243 * t258 + t350 * t398;
t40 = Icges(7,5) * t90 + Icges(7,6) * t89 + Icges(7,3) * t263;
t42 = Icges(7,4) * t90 + Icges(7,2) * t89 + Icges(7,6) * t263;
t44 = Icges(7,1) * t90 + Icges(7,4) * t89 + Icges(7,5) * t263;
t91 = Icges(7,5) * t169 + Icges(7,6) * t168 + Icges(7,3) * t401;
t10 = (-qJD(3) * t325 + t40) * t248 + (qJD(3) * t91 + t247 * t42 - t250 * t44 + (t247 * t95 + t250 * t93) * qJD(6)) * t251;
t19 = t123 * t401 + t124 * t168 + t125 * t169 + t172 * t263 + t173 * t89 + t174 * t90;
t366 = t19 / 0.2e1 + t10 / 0.2e1;
t32 = t248 * t91 + t251 * t325;
t53 = t168 * t173 + t169 * t174 + t172 * t401;
t365 = t53 / 0.2e1 + t32 / 0.2e1;
t33 = t248 * t92 + t251 * t324;
t54 = t170 * t173 + t171 * t174 + t172 * t397;
t364 = -t54 / 0.2e1 - t33 / 0.2e1;
t235 = pkin(7) * t380;
t363 = t235 + t389;
t362 = rSges(6,1) * t354 + rSges(6,2) * t475;
t158 = rSges(5,1) * t397 + rSges(5,3) * t399 + t237;
t361 = t244 * pkin(2) + t466;
t360 = -t442 + t446;
t353 = -t375 / 0.2e1;
t352 = -pkin(4) * t248 - t217;
t131 = t386 * t244;
t181 = pkin(4) * t401 + t406;
t349 = t243 * t181 + t244 * t182 + t394;
t338 = rSges(6,1) * t251 + t432;
t348 = t338 + t352;
t347 = -pkin(4) * t374 - t180;
t346 = t90 * rSges(7,1) + t89 * rSges(7,2);
t343 = t248 * t429 - pkin(2);
t342 = -t306 * qJD(3) / 0.2e1 + t493 * t469;
t337 = rSges(6,1) * t248 - rSges(6,2) * t251;
t27 = t168 * t93 + t169 * t95 + t401 * t91;
t28 = t168 * t94 + t169 * t96 + t401 * t92;
t16 = t243 * t28 - t244 * t27;
t333 = t243 * t27 + t244 * t28;
t29 = t170 * t93 + t171 * t95 + t397 * t91;
t30 = t170 * t94 + t171 * t96 + t397 * t92;
t17 = t243 * t30 - t244 * t29;
t332 = t243 * t29 + t244 * t30;
t331 = t243 * t33 - t244 * t32;
t330 = t243 * t32 + t244 * t33;
t239 = t244 * pkin(7);
t257 = t248 * t438 + t251 * t360 - pkin(2);
t253 = t243 * t257 - t406 - t440;
t50 = t239 + t253 + t336;
t323 = t361 + t176;
t285 = t230 + t323;
t51 = t285 - t407 + t425;
t329 = t243 * t51 + t244 * t50;
t69 = t177 * t401 - t248 * t97;
t70 = -t177 * t397 + t248 * t98;
t328 = t243 * t70 + t244 * t69;
t260 = (rSges(6,2) + t446) * t251 + t343;
t282 = t244 * t427 - t440;
t72 = t243 * t260 + t239 + t282;
t223 = rSges(6,1) * t399;
t73 = t243 * t427 + t223 + t285 - t368;
t327 = t243 * t73 + t244 * t72;
t326 = t243 * t98 - t244 * t97;
t307 = -Icges(5,3) * t251 + t416;
t292 = t352 - t392;
t291 = -t337 * qJD(3) + t347;
t290 = -pkin(2) - t340;
t289 = t243 * ((pkin(4) * t378 + qJD(5)) * t244 + t387) + t244 * (-pkin(4) * t475 - qJ(5) * t380 - t372) + t181 * t380 + t370;
t288 = t363 - t372;
t287 = -qJD(5) * t244 + t204 - t387;
t128 = t348 * t244;
t286 = -rSges(6,3) * t243 - t368;
t284 = t251 * t39 - t375 * t92;
t283 = t251 * t40 - t375 * t91;
t126 = (rSges(7,1) * t247 + rSges(7,2) * t250) * t371 + (rSges(7,3) * t251 + t248 * t335) * qJD(3);
t281 = -t344 * qJD(3) - t126 + t347;
t266 = qJD(3) * t307;
t80 = t292 * t244;
t261 = t248 * t428 + t251 * t443 - pkin(2);
t255 = t243 * t261 - t440;
t205 = pkin(5) * t354;
t195 = t340 * qJD(3);
t157 = t223 + t286;
t155 = t243 * t339 - t433;
t154 = t243 * t337 + t430;
t130 = t386 * t243;
t127 = t348 * t243;
t120 = t159 + t361;
t119 = t243 * t290 + t239 + t431 - t440;
t82 = t323 + t158;
t81 = t239 + t255 + t433;
t79 = t292 * t243;
t78 = t219 * t377 + (-t245 + (-rSges(4,3) - pkin(7)) * t243 + t290 * t244) * qJD(1);
t77 = t235 + (-t440 + (-pkin(2) - t435) * t243) * qJD(1) + t256;
t76 = qJD(1) * t131 + t243 * t391;
t75 = t218 * t381 + t244 * t391 + t183;
t68 = qJD(1) * t128 + t243 * t291;
t67 = t244 * t291 - t338 * t381 + t390;
t52 = t155 * t243 + t158 * t244 + t394;
t49 = t326 * t251;
t48 = t204 + (-t373 + (t251 * t428 + t434) * qJD(3)) * t243 + (-t245 + (-rSges(5,2) - pkin(7)) * t243 + t261 * t244) * qJD(1);
t47 = qJD(1) * t255 + t355 * t443 + t363 + t388;
t46 = rSges(7,3) * t263 + t346;
t45 = -rSges(7,3) * t475 + t436;
t38 = t154 * t243 + t157 * t244 + t349;
t37 = (-t373 + (t251 * t429 - t432) * qJD(3)) * t243 + (-t245 + (rSges(6,3) - pkin(7)) * t243 + t260 * t244) * qJD(1) + t287;
t36 = t446 * t355 + ((t251 * t446 + t343) * t243 + t282) * qJD(1) + t288 + t362;
t35 = qJD(1) * t80 + t243 * t281;
t34 = t244 * t281 + t381 * t392 + t390;
t26 = t243 * t426 + t244 * t425 + t349;
t24 = (-t373 + (t248 * t442 + t251 * t438) * qJD(3)) * t243 + (t244 * t257 - t466) * qJD(1) + t287 - t346;
t23 = qJD(1) * t253 + t355 * t360 + t205 + t288 + t436;
t22 = (-t177 * t377 - t46) * t248 + (-qJD(3) * t97 + t126 * t243 + t177 * t380) * t251;
t21 = (t177 * t376 + t45) * t248 + (qJD(3) * t98 - t126 * t244 + t177 * t381) * t251;
t20 = t244 * t388 + (-t218 * t241 - t242 * t434) * qJD(3) + (t244 * t155 + (-t158 - t176 + t237) * t243) * qJD(1) + t370;
t15 = ((t154 - t430) * qJD(1) + t362) * t244 + (t338 * t377 + (-t157 + t286 + t393) * qJD(1)) * t243 + t289;
t14 = t326 * t375 + (-t243 * t45 + t244 * t46 + (-t243 * t97 - t244 * t98) * qJD(1)) * t251;
t13 = t248 * t54 + t251 * t332;
t12 = t248 * t53 + t251 * t333;
t9 = t168 * t41 + t169 * t43 + t243 * t284 + t357 * t92 + t89 * t94 + t90 * t96;
t8 = t168 * t42 + t169 * t44 + t243 * t283 + t357 * t91 + t89 * t93 + t90 * t95;
t7 = t170 * t41 + t171 * t43 + t244 * t284 - t358 * t92 + t87 * t94 + t88 * t96;
t6 = t170 * t42 + t171 * t44 + t244 * t283 - t358 * t91 + t87 * t93 + t88 * t95;
t5 = (-pkin(8) * t355 + t205 + t45) * t244 + (t345 * t377 + t46) * t243 + (t426 * t244 + (t393 - t425) * t243) * qJD(1) + t289;
t4 = qJD(1) * t333 + t243 * t9 - t244 * t8;
t3 = qJD(1) * t332 + t243 * t7 - t244 * t6;
t2 = (-qJD(3) * t333 + t19) * t248 + (-qJD(1) * t16 + qJD(3) * t53 + t243 * t8 + t244 * t9) * t251;
t1 = (-qJD(3) * t332 + t18) * t248 + (-t17 * qJD(1) + qJD(3) * t54 + t243 * t6 + t244 * t7) * t251;
t25 = [t254 + (t119 * t78 + t120 * t77) * t453 + (t47 * t82 + t48 * t81) * t452 + (t36 * t73 + t37 * t72) * t451 + (t23 * t51 + t24 * t50) * t450 + (t307 - t311 + t322 + t320 - t492) * t375 + (-t308 - t317 + t315 + t491) * t374; 0; 0; m(5) * (t130 * t47 + t131 * t48 + t75 * t81 + t76 * t82) + m(6) * (t127 * t36 + t128 * t37 + t67 * t72 + t68 * t73) + m(7) * (t23 * t79 + t24 * t80 + t34 * t50 + t35 * t51) + (m(4) * (-t119 * t195 - t219 * t78) + t342 * t244 + (t146 * t476 + t266 * t479 + t439 * t488 + t478 * t485) * t251 - t366) * t244 + (m(4) * (-t120 * t195 - t219 * t77) + t342 * t243 + (t145 * t476 + t266 * t445 + t439 * t489 + t477 * t485) * t251 + t367) * t243 + ((t141 * t243 + t142 * t244) * t439 + (t243 * t487 + t244 * t486) * t476 + (t243 * t477 + t244 * t478) * t491 * qJD(3)) * t248 + (t277 / 0.2e1 - t279 / 0.2e1 + t275 / 0.2e1 - t276 / 0.2e1 + t278 / 0.2e1 - t274 / 0.2e1) * qJD(3) + ((-t120 * t441 + (t146 / 0.2e1 - t148 / 0.2e1 - t138 / 0.2e1) * t251 + (t152 / 0.2e1 - t142 / 0.2e1 + t150 / 0.2e1) * t248 - t364) * t244 + (t119 * t441 + (-t147 / 0.2e1 - t137 / 0.2e1 + t145 / 0.2e1) * t251 + (-t141 / 0.2e1 + t149 / 0.2e1 + t151 / 0.2e1) * t248 + t365) * t243) * qJD(1); m(4) * t31 + m(5) * t20 + m(6) * t15 + m(7) * t5; (t26 * t5 + t34 * t80 + t35 * t79) * t450 + t17 * t380 + t16 * t381 + (t127 * t68 + t128 * t67 + t15 * t38) * t451 + (t130 * t76 + t131 * t75 + t20 * t52) * t452 + t385 * t219 * t195 * t453 + (t159 * t468 + t460 * t242 + t465 * t381 - t4 + (-t244 * t480 - t473) * t380) * t244 + (t3 + t156 * t468 + t461 * t241 + t462 * t380 + ((t460 - t472) * t243 + t461 * t244 - t471 * t375 + t470 * t374 + (t248 * t471 - t251 * t470) * qJD(3) + ((-t480 + t483) * t243 + t490 + t462 + t465) * qJD(1)) * t244 + (t243 * t455 + t473) * t381) * t243; 0.2e1 * (t329 * t447 + t327 * t448 + (t243 * t82 + t244 * t81) * t449) * t374 + 0.2e1 * ((t23 * t243 + t24 * t244 + t380 * t51 - t381 * t50) * t447 + (t243 * t36 + t244 * t37 + t380 * t73 - t381 * t72) * t448 + (t243 * t47 + t244 * t48 + t380 * t82 - t381 * t81) * t449) * t248; (m(5) + m(6) + m(7)) * t375; 0.2e1 * ((t376 * t80 + t377 * t79 - t5) * t447 + (t127 * t377 + t128 * t376 - t15) * t448 + (t130 * t377 + t131 * t376 - t20) * t449) * t251 + 0.2e1 * ((qJD(3) * t26 + t243 * t35 + t244 * t34 + t380 * t79 - t381 * t80) * t447 + (qJD(3) * t38 + t127 * t380 - t128 * t381 + t243 * t68 + t244 * t67) * t448 + (qJD(3) * t52 + t130 * t380 - t131 * t381 + t243 * t76 + t244 * t75) * t449) * t248; 0.4e1 * (t449 + t448 + t447) * (-0.1e1 + t385) * t248 * t374; m(7) * (-qJD(1) * t329 + t23 * t244 - t24 * t243) + m(6) * (-qJD(1) * t327 - t243 * t37 + t244 * t36); 0; m(7) * (-t243 * t34 + t244 * t35 + (-t243 * t79 - t244 * t80) * qJD(1)) + m(6) * (-t243 * t67 + t244 * t68 + (-t127 * t243 - t128 * t244) * qJD(1)); 0; 0; m(7) * (t21 * t51 + t22 * t50 + t23 * t70 + t24 * t69) + (-t243 * t365 + t244 * t364) * t375 + (t367 * t244 + t366 * t243 + (t243 * t364 + t244 * t365) * qJD(1)) * t251 + t437; m(7) * t14; m(7) * (t14 * t26 + t21 * t79 + t22 * t80 + t34 * t69 + t35 * t70 - t49 * t5) + (t17 * t353 + (qJD(1) * t33 - t10) * t444 - t2 / 0.2e1 + t13 * t439) * t244 + ((qJD(1) * t32 + t11) * t444 + t12 * t439 + t1 / 0.2e1 + t16 * t353) * t243 + (t3 * t445 + t331 * t469 + t4 * t478 + (t16 * t445 + t17 * t479) * qJD(1)) * t251; m(7) * ((qJD(3) * t328 - t14) * t251 + (-qJD(3) * t49 + t21 * t243 + t22 * t244 + (-t243 * t69 + t244 * t70) * qJD(1)) * t248); m(7) * (-qJD(1) * t328 + t21 * t244 - t22 * t243); (-t14 * t49 + t21 * t70 + t22 * t69) * t450 + ((-t243 * t12 - t244 * t13 - t248 * t330) * qJD(3) + t437) * t248 + (t244 * t1 + t243 * t2 + t248 * (t10 * t243 + t11 * t244) + (t248 * t74 + t251 * t330) * qJD(3) + (t244 * t12 - t243 * t13 - t248 * t331) * qJD(1)) * t251;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;

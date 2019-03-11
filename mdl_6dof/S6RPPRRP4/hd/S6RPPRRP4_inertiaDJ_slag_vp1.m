% Calculate time derivative of joint inertia matrix for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:11
% EndTime: 2019-03-09 02:05:34
% DurationCPUTime: 14.47s
% Computational Cost: add. (19952->693), mult. (49132->953), div. (0->0), fcn. (55882->8), ass. (0->327)
t386 = sin(pkin(9));
t387 = cos(pkin(9));
t402 = sin(qJ(1));
t403 = cos(qJ(1));
t236 = -t386 * t402 - t387 * t403;
t237 = t403 * t386 - t402 * t387;
t263 = cos(qJ(5));
t261 = sin(qJ(5));
t264 = cos(qJ(4));
t367 = t261 * t264;
t182 = t236 * t263 - t237 * t367;
t365 = t263 * t264;
t183 = -t236 * t261 - t237 * t365;
t262 = sin(qJ(4));
t370 = t237 * t262;
t122 = Icges(7,5) * t183 - Icges(7,6) * t370 + Icges(7,3) * t182;
t126 = Icges(7,4) * t183 - Icges(7,2) * t370 + Icges(7,6) * t182;
t130 = Icges(7,1) * t183 - Icges(7,4) * t370 + Icges(7,5) * t182;
t47 = t122 * t182 - t126 * t370 + t130 * t183;
t124 = Icges(6,5) * t183 - Icges(6,6) * t182 - Icges(6,3) * t370;
t128 = Icges(6,4) * t183 - Icges(6,2) * t182 - Icges(6,6) * t370;
t132 = Icges(6,1) * t183 - Icges(6,4) * t182 - Icges(6,5) * t370;
t49 = -t124 * t370 - t128 * t182 + t132 * t183;
t451 = t47 + t49;
t184 = -t236 * t367 - t237 * t263;
t185 = -t236 * t365 + t237 * t261;
t372 = t236 * t262;
t123 = Icges(7,5) * t185 - Icges(7,6) * t372 + Icges(7,3) * t184;
t127 = Icges(7,4) * t185 - Icges(7,2) * t372 + Icges(7,6) * t184;
t131 = Icges(7,1) * t185 - Icges(7,4) * t372 + Icges(7,5) * t184;
t48 = t123 * t182 - t127 * t370 + t131 * t183;
t125 = Icges(6,5) * t185 - Icges(6,6) * t184 - Icges(6,3) * t372;
t129 = Icges(6,4) * t185 - Icges(6,2) * t184 - Icges(6,6) * t372;
t133 = Icges(6,1) * t185 - Icges(6,4) * t184 - Icges(6,5) * t372;
t50 = -t125 * t370 - t129 * t182 + t133 * t183;
t450 = t48 + t50;
t51 = t122 * t184 - t126 * t372 + t130 * t185;
t53 = -t124 * t372 - t128 * t184 + t132 * t185;
t449 = t51 + t53;
t52 = t123 * t184 - t127 * t372 + t131 * t185;
t54 = -t125 * t372 - t129 * t184 + t133 * t185;
t448 = t52 + t54;
t231 = t237 * qJD(1);
t350 = qJD(4) * t262;
t335 = t261 * t350;
t230 = t236 * qJD(1);
t364 = t264 * t230;
t111 = qJD(5) * t183 - t231 * t263 + t237 * t335 - t261 * t364;
t337 = t237 * t350;
t286 = t337 - t364;
t347 = qJD(5) * t264;
t112 = (t237 * t347 + t231) * t261 + (-qJD(5) * t236 + t286) * t263;
t349 = qJD(4) * t264;
t336 = t237 * t349;
t374 = t230 * t262;
t287 = t336 + t374;
t67 = Icges(7,5) * t112 - Icges(7,6) * t287 + Icges(7,3) * t111;
t73 = Icges(6,4) * t112 - Icges(6,2) * t111 - Icges(6,6) * t287;
t468 = t67 - t73;
t363 = t264 * t231;
t113 = qJD(5) * t185 - t230 * t263 + t236 * t335 + t261 * t363;
t284 = t236 * t350 + t363;
t114 = (t236 * t347 + t230) * t261 + (qJD(5) * t237 + t284) * t263;
t285 = -t231 * t262 + t236 * t349;
t68 = Icges(7,5) * t114 - Icges(7,6) * t285 + Icges(7,3) * t113;
t74 = Icges(6,4) * t114 - Icges(6,2) * t113 - Icges(6,6) * t285;
t467 = t68 - t74;
t69 = Icges(6,5) * t112 - Icges(6,6) * t111 - Icges(6,3) * t287;
t71 = Icges(7,4) * t112 - Icges(7,2) * t287 + Icges(7,6) * t111;
t466 = t69 + t71;
t70 = Icges(6,5) * t114 - Icges(6,6) * t113 - Icges(6,3) * t285;
t72 = Icges(7,4) * t114 - Icges(7,2) * t285 + Icges(7,6) * t113;
t465 = -t70 - t72;
t75 = Icges(7,1) * t112 - Icges(7,4) * t287 + Icges(7,5) * t111;
t77 = Icges(6,1) * t112 - Icges(6,4) * t111 - Icges(6,5) * t287;
t464 = t75 + t77;
t76 = Icges(7,1) * t114 - Icges(7,4) * t285 + Icges(7,5) * t113;
t78 = Icges(6,1) * t114 - Icges(6,4) * t113 - Icges(6,5) * t285;
t463 = t76 + t78;
t462 = t122 - t128;
t461 = t123 - t129;
t460 = t124 + t126;
t459 = -t125 - t127;
t458 = t130 + t132;
t457 = t131 + t133;
t388 = rSges(7,3) + qJ(6);
t406 = rSges(7,1) + pkin(5);
t456 = t388 * t182 + t406 * t183;
t455 = t111 * t462 + t112 * t458 + t182 * t468 + t183 * t464 - t287 * t460 - t370 * t466;
t454 = t111 * t461 + t112 * t457 + t182 * t467 + t183 * t463 + t287 * t459 + t370 * t465;
t453 = -t113 * t462 - t114 * t458 - t184 * t468 - t185 * t464 + t285 * t460 + t372 * t466;
t452 = t113 * t461 + t114 * t457 + t184 * t467 + t185 * t463 + t285 * t459 + t372 * t465;
t380 = Icges(7,5) * t263;
t304 = -Icges(7,3) * t261 - t380;
t204 = Icges(7,6) * t264 + t262 * t304;
t307 = -Icges(7,4) * t263 - Icges(7,6) * t261;
t206 = Icges(7,2) * t264 + t262 * t307;
t381 = Icges(7,5) * t261;
t310 = -Icges(7,1) * t263 - t381;
t208 = Icges(7,4) * t264 + t262 * t310;
t96 = t182 * t204 + t183 * t208 - t206 * t370;
t305 = -Icges(6,5) * t263 + Icges(6,6) * t261;
t205 = Icges(6,3) * t264 + t262 * t305;
t382 = Icges(6,4) * t263;
t308 = Icges(6,2) * t261 - t382;
t207 = Icges(6,6) * t264 + t262 * t308;
t383 = Icges(6,4) * t261;
t311 = -Icges(6,1) * t263 + t383;
t209 = Icges(6,5) * t264 + t262 * t311;
t97 = -t182 * t207 + t183 * t209 - t205 * t370;
t447 = t96 + t97;
t98 = t184 * t204 + t185 * t208 - t206 * t372;
t99 = -t184 * t207 + t185 * t209 - t205 * t372;
t446 = t98 + t99;
t445 = -t209 - t208;
t444 = -t236 * t448 - t237 * t449;
t443 = -t450 * t236 - t237 * t451;
t385 = Icges(5,4) * t262;
t312 = -Icges(5,1) * t264 + t385;
t159 = Icges(5,5) * t237 + t236 * t312;
t376 = t159 * t264;
t384 = Icges(5,4) * t264;
t309 = Icges(5,2) * t262 - t384;
t157 = Icges(5,6) * t237 + t236 * t309;
t378 = t157 * t262;
t297 = -t376 + t378;
t158 = -Icges(5,5) * t236 + t237 * t312;
t377 = t158 * t264;
t156 = -Icges(5,6) * t236 + t237 * t309;
t379 = t156 * t262;
t298 = -t377 + t379;
t348 = qJD(5) * t262;
t169 = (-Icges(7,3) * t263 + t381) * t348 + (-Icges(7,6) * t262 + t264 * t304) * qJD(4);
t171 = (Icges(7,4) * t261 - Icges(7,6) * t263) * t348 + (-Icges(7,2) * t262 + t264 * t307) * qJD(4);
t173 = (Icges(7,1) * t261 - t380) * t348 + (-Icges(7,4) * t262 + t264 * t310) * qJD(4);
t37 = t111 * t204 + t112 * t208 + t182 * t169 - t171 * t370 + t183 * t173 - t206 * t287;
t170 = (Icges(6,5) * t261 + Icges(6,6) * t263) * t348 + (-Icges(6,3) * t262 + t264 * t305) * qJD(4);
t172 = (Icges(6,2) * t263 + t383) * t348 + (-Icges(6,6) * t262 + t264 * t308) * qJD(4);
t174 = (Icges(6,1) * t261 + t382) * t348 + (-Icges(6,5) * t262 + t264 * t311) * qJD(4);
t38 = -t111 * t207 + t112 * t209 - t170 * t370 - t182 * t172 + t183 * t174 - t205 * t287;
t442 = (-t443 * qJD(4) - t37 - t38) * t264 + (t447 * qJD(4) + t451 * t230 - t450 * t231 + t454 * t236 + t455 * t237) * t262;
t39 = t113 * t204 + t114 * t208 + t184 * t169 - t171 * t372 + t185 * t173 - t206 * t285;
t40 = -t113 * t207 + t114 * t209 - t170 * t372 - t184 * t172 + t185 * t174 - t205 * t285;
t441 = (t444 * qJD(4) + t39 + t40) * t264 + (-t446 * qJD(4) - t449 * t230 + t448 * t231 - t452 * t236 + t453 * t237) * t262;
t303 = -t122 * t261 - t130 * t263;
t19 = (qJD(4) * t303 + t71) * t264 + (-qJD(4) * t126 - t261 * t67 - t263 * t75 + (-t122 * t263 + t130 * t261) * qJD(5)) * t262;
t301 = t128 * t261 - t132 * t263;
t21 = (qJD(4) * t301 + t69) * t264 + (-qJD(4) * t124 + t261 * t73 - t263 * t77 + (t128 * t263 + t132 * t261) * qJD(5)) * t262;
t440 = -t19 - t21;
t302 = -t123 * t261 - t131 * t263;
t20 = (qJD(4) * t302 + t72) * t264 + (-qJD(4) * t127 - t261 * t68 - t263 * t76 + (-t123 * t263 + t131 * t261) * qJD(5)) * t262;
t300 = t129 * t261 - t133 * t263;
t22 = (qJD(4) * t300 + t70) * t264 + (-qJD(4) * t125 + t261 * t74 - t263 * t78 + (t129 * t263 + t133 * t261) * qJD(5)) * t262;
t439 = t20 + t22;
t438 = t443 * t262 + t447 * t264;
t437 = t444 * t262 + t446 * t264;
t57 = t126 * t264 + t262 * t303;
t59 = t124 * t264 + t262 * t301;
t395 = t57 + t59;
t58 = t127 * t264 + t262 * t302;
t60 = t125 * t264 + t262 * t300;
t394 = t58 + t60;
t419 = 2 * m(5);
t321 = -rSges(5,1) * t264 + rSges(5,2) * t262;
t389 = t236 * rSges(5,3);
t160 = t237 * t321 - t389;
t371 = t236 * t264;
t161 = -rSges(5,1) * t371 + rSges(5,2) * t372 + t237 * rSges(5,3);
t282 = rSges(5,1) * t284 + t230 * rSges(5,3);
t390 = t231 * rSges(5,3);
t55 = t230 * t160 + t237 * (rSges(5,1) * t286 + t390) - t231 * t161 + t236 * t282 + (t236 * t285 + t237 * t287) * rSges(5,2);
t436 = t419 * t55;
t248 = -t262 * rSges(5,1) - rSges(5,2) * t264;
t435 = t248 ^ 2 * t419;
t344 = rSges(7,2) * t370;
t362 = -t344 + t456;
t434 = -t206 - t205;
t352 = t403 * pkin(1) + t402 * qJ(2);
t340 = t403 * pkin(2) + t352;
t366 = t262 * t263;
t368 = t261 * t262;
t331 = t263 * t348;
t425 = -t261 * t349 - t331;
t433 = t172 * t368 - t174 * t366 - t207 * t425 - t445 * t261 * t348 + (t170 + t171) * t264;
t432 = -rSges(3,1) * t403 - rSges(3,3) * t402 - t352;
t272 = -t182 * qJD(6) - t111 * t388 - t112 * t406;
t431 = t450 * t230 + t451 * t231 - t455 * t236 + t454 * t237;
t430 = t448 * t230 + t449 * t231 + t453 * t236 + t452 * t237;
t306 = -Icges(5,5) * t264 + Icges(5,6) * t262;
t154 = -Icges(5,3) * t236 + t237 * t306;
t155 = Icges(5,3) * t237 + t236 * t306;
t427 = t298 * t236;
t428 = t297 * t237;
t429 = t237 * t154 - t236 * t155 + t427 + t428;
t317 = -rSges(7,1) * t263 - rSges(7,3) * t261;
t391 = rSges(7,2) * t264;
t355 = t391 + (-pkin(5) * t263 - qJ(6) * t261 + t317) * t262;
t329 = t355 * t236;
t242 = t321 * qJD(4);
t426 = t242 * t248 * t419;
t399 = pkin(4) * t264;
t339 = pkin(3) + t399;
t404 = rSges(6,3) + pkin(8);
t424 = t262 * t404 + t339;
t405 = rSges(7,2) + pkin(8);
t423 = t262 * t405 + t339;
t422 = t264 * t394 + t437;
t421 = t264 * t395 + t438;
t357 = (-pkin(5) * t349 - qJ(6) * t348) * t263 + (-qJ(6) * t349 + (pkin(5) * qJD(5) - qJD(6)) * t262) * t261 + (rSges(7,1) * t261 - rSges(7,3) * t263) * t348 + (-rSges(7,2) * t262 + t264 * t317) * qJD(4);
t420 = -t231 * t355 + t236 * t357;
t418 = 2 * m(6);
t417 = 2 * m(7);
t411 = t306 * qJD(4) / 0.2e1;
t410 = -t262 / 0.2e1;
t409 = t262 / 0.2e1;
t408 = -t264 / 0.2e1;
t407 = t264 / 0.2e1;
t401 = m(5) * t242;
t400 = m(5) * t248;
t295 = -t204 * t261 - t208 * t263;
t393 = (t295 * t349 + (-qJD(4) * t206 - t261 * t169 + (-qJD(5) * t204 - t173) * t263) * t262 + (-t205 * t262 - t209 * t365) * qJD(4) + t433) * t264;
t392 = -rSges(7,2) * t285 + t184 * qJD(6) + t113 * t388 + t114 * t406;
t318 = -rSges(6,1) * t263 + rSges(6,2) * t261;
t213 = rSges(6,3) * t264 + t262 * t318;
t375 = t213 * t236;
t187 = -pkin(4) * t371 - pkin(8) * t372;
t277 = pkin(4) * t284 - pkin(8) * t285;
t361 = -t231 * t187 + t236 * t277;
t360 = rSges(7,2) * t372 - t184 * t388 - t185 * t406;
t359 = t434 * t264 + (-t207 * t261 + t209 * t263 - t295) * t262;
t322 = -pkin(8) * t262 - t399;
t243 = t322 * qJD(4);
t249 = -t262 * pkin(4) + pkin(8) * t264;
t356 = t231 * t249 - t236 * t243;
t354 = t213 + t249;
t258 = t403 * qJ(2);
t353 = qJD(1) * t258 + qJD(2) * t402;
t351 = qJD(4) * t237;
t345 = t402 * pkin(1);
t31 = -t236 * t47 + t237 * t48;
t32 = -t236 * t49 + t237 * t50;
t343 = t31 / 0.2e1 + t32 / 0.2e1;
t33 = -t236 * t51 + t237 * t52;
t34 = -t236 * t53 + t237 * t54;
t342 = t33 / 0.2e1 + t34 / 0.2e1;
t138 = t185 * rSges(6,1) - t184 * rSges(6,2) - rSges(6,3) * t372;
t341 = t249 + t355;
t328 = t360 * qJD(4);
t83 = t114 * rSges(6,1) - t113 * rSges(6,2) - rSges(6,3) * t285;
t325 = -t355 + t391;
t320 = t112 * rSges(6,1) - t111 * rSges(6,2);
t319 = -t183 * rSges(6,1) + t182 * rSges(6,2);
t136 = -rSges(6,3) * t370 - t319;
t299 = -t136 * t236 + t138 * t237;
t176 = (rSges(6,1) * t261 + rSges(6,2) * t263) * t348 + (-rSges(6,3) * t262 + t264 * t318) * qJD(4);
t296 = t176 * t236 - t213 * t231;
t292 = pkin(3) - t321;
t291 = -t37 / 0.2e1 - t38 / 0.2e1 - t19 / 0.2e1 - t21 / 0.2e1;
t290 = -t39 / 0.2e1 - t40 / 0.2e1 - t20 / 0.2e1 - t22 / 0.2e1;
t289 = t59 / 0.2e1 + t96 / 0.2e1 + t97 / 0.2e1 + t57 / 0.2e1;
t288 = t98 / 0.2e1 + t99 / 0.2e1 + t58 / 0.2e1 + t60 / 0.2e1;
t283 = -t236 * pkin(3) + pkin(7) * t237 + t340;
t281 = -pkin(2) * t402 - t345;
t276 = t258 + t281;
t275 = t283 + t187;
t273 = t236 * pkin(7) + t276;
t271 = -rSges(3,1) * t402 + rSges(3,3) * t403 - t345;
t270 = qJD(1) * t281 + t353;
t256 = qJD(2) * t403;
t269 = -qJD(1) * t340 + t256;
t268 = t231 * pkin(3) + t230 * pkin(7) + t270;
t267 = -t231 * pkin(7) + t269;
t216 = pkin(4) * t337;
t266 = -t216 + t267;
t265 = t268 + t277;
t245 = -Icges(5,5) * t262 - Icges(5,6) * t264;
t218 = t258 + t271;
t196 = t236 * t249;
t189 = qJD(1) * t432 + t256;
t188 = qJD(1) * t271 + t353;
t186 = t322 * t237;
t167 = -rSges(4,1) * t236 - rSges(4,2) * t237 + t340;
t166 = t237 * rSges(4,1) - t236 * rSges(4,2) + t276;
t162 = t237 * t186;
t152 = -t196 - t375;
t151 = t354 * t237;
t150 = t230 * rSges(4,1) + t231 * rSges(4,2) + t269;
t149 = t231 * rSges(4,1) - t230 * rSges(4,2) + t270;
t144 = -pkin(4) * t364 - pkin(8) * t287 + t216;
t143 = -t196 - t329;
t142 = t341 * t237;
t140 = t161 + t283;
t139 = t237 * t292 + t273 + t389;
t116 = Icges(5,5) * t284 + Icges(5,6) * t285 + Icges(5,3) * t230;
t115 = Icges(5,5) * t286 + Icges(5,6) * t287 + Icges(5,3) * t231;
t103 = t138 * t264 + t213 * t372;
t102 = -t136 * t264 - t213 * t370;
t101 = -t296 + t356;
t100 = (-t176 - t243) * t237 - t354 * t230;
t95 = t230 * t292 + t248 * t351 + t267 - t390;
t94 = rSges(5,2) * t285 + t268 + t282;
t93 = t275 + t138;
t92 = t237 * t424 + t273 + t319;
t87 = t299 * t262;
t86 = t356 - t420;
t85 = (-t243 - t357) * t237 - t341 * t230;
t81 = -rSges(6,3) * t287 + t320;
t80 = t262 * t329 - t264 * t360;
t79 = -t264 * t362 - t355 * t370;
t66 = t275 - t360;
t65 = t237 * t423 + t273 - t456;
t56 = t136 * t237 + t162 + (t138 + t187) * t236;
t46 = t230 * t424 + t404 * t336 + t266 - t320;
t45 = t265 + t83;
t44 = (-t236 * t362 - t237 * t360) * t262;
t43 = (qJD(4) * t375 + t83) * t264 + (-qJD(4) * t138 + t296) * t262;
t42 = (-t213 * t351 - t81) * t264 + (qJD(4) * t136 - t176 * t237 - t213 * t230) * t262;
t41 = t162 + t362 * t237 + (t187 - t360) * t236;
t36 = t230 * t423 + t405 * t336 + t266 + t272;
t35 = t265 + t392;
t30 = (qJD(4) * t329 + t392) * t264 + (t328 + t420) * t262;
t29 = (t325 * t351 + t272) * t264 + (qJD(4) * t362 + t230 * t325 - t237 * t357) * t262;
t24 = t299 * t349 + (t136 * t231 + t138 * t230 - t236 * t81 + t237 * t83) * t262;
t23 = -t138 * t231 + t236 * t83 + (t144 + t81) * t237 - (-t136 - t186) * t230 + t361;
t10 = t392 * t236 + t360 * t231 + (-rSges(7,2) * t336 + t144 - t272) * t237 - (-t186 + t344 - t362) * t230 + t361;
t9 = t360 * t374 + (-t262 * t392 + t264 * t328) * t237 + (-t287 * rSges(7,2) - t272) * t372 + t362 * t285;
t1 = [-t173 * t366 - t169 * t368 + 0.2e1 * m(3) * (-t188 * t432 + t189 * t218) + 0.2e1 * m(4) * (t149 * t167 + t150 * t166) + (t139 * t95 + t140 * t94) * t419 + (t45 * t93 + t46 * t92) * t418 + (t35 * t66 + t36 * t65) * t417 + t425 * t204 + (-Icges(5,2) * t264 - t312 - t385 + t434) * t350 + (Icges(5,1) * t262 + t445 * t263 - t309 + t384) * t349 + t433; m(7) * (t402 * t36 - t403 * t35 + (t402 * t66 + t403 * t65) * qJD(1)) + m(6) * (t402 * t46 - t403 * t45 + (t402 * t93 + t403 * t92) * qJD(1)) + m(5) * (t402 * t95 - t403 * t94 + (t139 * t403 + t140 * t402) * qJD(1)) + m(3) * (t402 * t189 - t403 * t188 + (t218 * t403 - t402 * t432) * qJD(1)) + m(4) * (t402 * t150 - t403 * t149 + (t166 * t403 + t167 * t402) * qJD(1)); 0; 0; 0; 0; ((Icges(5,4) * t284 + Icges(5,2) * t285 + Icges(5,6) * t230) * t408 + (Icges(5,1) * t284 + Icges(5,4) * t285 + Icges(5,5) * t230) * t410 + t237 * t411 + (t378 / 0.2e1 - t376 / 0.2e1) * qJD(4) - t290) * t237 + ((Icges(5,4) * t286 + Icges(5,2) * t287 + Icges(5,6) * t231) * t407 + (Icges(5,1) * t286 + Icges(5,4) * t287 + Icges(5,5) * t231) * t409 + t236 * t411 + (-t379 / 0.2e1 + t377 / 0.2e1) * qJD(4) + t291) * t236 + (t156 * t408 + t158 * t410 - t236 * t245 + t289) * t231 - (t157 * t407 + t159 * t409 - t237 * t245 - t288) * t230 + m(6) * (t100 * t93 + t101 * t92 - t151 * t45 + t152 * t46) + m(7) * (-t142 * t35 + t143 * t36 + t65 * t86 + t66 * t85) + (-t139 * t236 - t140 * t237) * t401 + (t139 * t231 - t140 * t230 - t236 * t95 - t237 * t94) * t400; m(6) * (t101 * t402 - t100 * t403 + (-t151 * t402 + t152 * t403) * qJD(1)) + m(7) * (t86 * t402 - t85 * t403 + (-t142 * t402 + t143 * t403) * qJD(1)) + (-t236 * t402 + t237 * t403) * t401 + (t402 * t231 + t403 * t230 + (-t236 * t403 - t237 * t402) * qJD(1)) * t400; -m(5) * t55 - m(6) * t23 - m(7) * t10; (t41 * t10 - t142 * t85 + t143 * t86) * t417 + (-t100 * t151 + t101 * t152 + t56 * t23) * t418 + (t32 + t31) * t231 - (-t34 - t33) * t230 + (t160 * t436 + (t237 * t116 + t426) * t237 + (-t428 + t429) * t231 + t430 + (0.3e1 * t237 * t155 + 0.2e1 * t236 * t297 + t435) * t230) * t237 + (t161 * t436 + (-t236 * t115 + t426) * t236 + (-t237 * t115 + t236 * t116) * t237 - t431 + (t427 - t429 + (-t154 - t297) * t237) * t230 + (0.3e1 * t236 * t154 - t435 + (-t155 - t298) * t237) * t231) * t236; m(6) * (t102 * t46 + t103 * t45 + t42 * t92 + t43 * t93) + m(7) * (t29 * t65 + t30 * t66 + t35 * t80 + t36 * t79) + (-t236 * t288 - t237 * t289) * t349 + (qJD(4) * t359 - t230 * t289 + t231 * t288 + t236 * t290 + t237 * t291) * t262 + t393; m(6) * (t42 * t402 - t43 * t403 + (t102 * t403 + t103 * t402) * qJD(1)) + m(7) * (t29 * t402 - t30 * t403 + (t402 * t80 + t403 * t79) * qJD(1)); -m(6) * t24 + m(7) * t9; m(6) * (t100 * t103 + t101 * t102 - t151 * t43 + t152 * t42 + t23 * t87 + t24 * t56) + m(7) * (t44 * t10 - t142 * t30 + t143 * t29 - t9 * t41 + t79 * t86 + t80 * t85) + (-t236 * t342 - t237 * t343) * t349 + (-t343 * t230 + t342 * t231 - t431 * t237 / 0.2e1 - (-t236 * t395 + t237 * t394) * qJD(4) / 0.2e1) * t262 + t437 * t230 / 0.2e1 + t438 * t231 / 0.2e1 - (t262 * t430 - t442) * t236 / 0.2e1 + t441 * t237 / 0.2e1 + (t230 * t394 + t231 * t395 + t236 * t440 + t237 * t439) * t407; (t29 * t79 + t30 * t80 - t44 * t9) * t417 + (t102 * t42 + t103 * t43 + t24 * t87) * t418 + ((-t236 * t422 - t237 * t421) * qJD(4) + t393) * t264 + ((t264 * t440 + t442) * t237 + (-t264 * t439 - t441) * t236 + t422 * t231 - t421 * t230 + ((t236 * t394 + t237 * t395) * t262 + 0.2e1 * t359 * t264) * qJD(4)) * t262; m(7) * (t111 * t66 + t113 * t65 + t182 * t35 + t184 * t36); m(7) * (t113 * t402 - t111 * t403 + (t182 * t402 + t184 * t403) * qJD(1)); -t425 * m(7); m(7) * (-t41 * t331 - t111 * t142 + t113 * t143 + t182 * t85 + t184 * t86 + (-t10 * t262 - t349 * t41) * t261); m(7) * (-t44 * t331 + t111 * t80 + t113 * t79 + t182 * t30 + t184 * t29 + (t262 * t9 - t349 * t44) * t261); (t182 * t111 + t184 * t113 - t368 * t425) * t417;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;

% Calculate time derivative of joint inertia matrix for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_inertiaDJ_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:25:52
% EndTime: 2019-03-29 15:26:05
% DurationCPUTime: 9.89s
% Computational Cost: add. (30398->607), mult. (26809->882), div. (0->0), fcn. (25037->10), ass. (0->330)
t276 = qJD(3) + qJD(4);
t278 = qJ(3) + qJ(4);
t273 = cos(t278);
t271 = sin(t278);
t429 = Icges(5,4) * t271;
t339 = Icges(5,1) * t273 - t429;
t234 = Icges(5,2) * t273 + t429;
t402 = t276 * t234;
t467 = t339 * t276 - t402;
t280 = sin(qJ(5));
t283 = cos(qJ(5));
t332 = Icges(6,5) * t283 - Icges(6,6) * t280;
t409 = t273 * t276;
t159 = t332 * t409 + (Icges(6,3) * t276 + (-Icges(6,5) * t280 - Icges(6,6) * t283) * qJD(5)) * t271;
t427 = Icges(6,4) * t280;
t338 = Icges(6,1) * t283 - t427;
t426 = Icges(6,4) * t283;
t165 = t338 * t409 + (Icges(6,5) * t276 + (-Icges(6,1) * t280 - t426) * qJD(5)) * t271;
t198 = -Icges(6,3) * t273 + t271 * t332;
t335 = -Icges(6,2) * t280 + t426;
t201 = -Icges(6,6) * t273 + t271 * t335;
t204 = -Icges(6,5) * t273 + t271 * t338;
t400 = t276 * t283;
t401 = t276 * t280;
t415 = t271 * t276;
t466 = t271 * t283 * t165 + t198 * t415 + (-t201 * t401 + t204 * t400 - t159) * t273;
t279 = qJ(1) + qJ(2);
t274 = cos(t279);
t370 = t409 / 0.2e1;
t272 = sin(t279);
t277 = qJD(1) + qJD(2);
t413 = t272 * t277;
t384 = t271 * t413;
t465 = t274 * t370 - t384 / 0.2e1;
t406 = t274 * t277;
t369 = t406 / 0.2e1;
t464 = t271 * t369 + t272 * t370;
t428 = Icges(5,4) * t273;
t336 = -Icges(5,2) * t271 + t428;
t218 = t336 * t276;
t235 = Icges(5,1) * t271 + t428;
t281 = sin(qJ(3));
t284 = cos(qJ(3));
t430 = Icges(4,4) * t284;
t337 = -Icges(4,2) * t281 + t430;
t240 = t337 * qJD(3);
t431 = Icges(4,4) * t281;
t340 = Icges(4,1) * t284 - t431;
t241 = t340 * qJD(3);
t253 = Icges(4,2) * t284 + t431;
t254 = Icges(4,1) * t281 + t430;
t162 = t335 * t409 + (Icges(6,6) * t276 + (-Icges(6,2) * t283 - t427) * qJD(5)) * t271;
t293 = -t280 * t162 + (-t201 * t283 - t204 * t280) * qJD(5);
t392 = qJD(3) * t284;
t393 = qJD(3) * t281;
t286 = t273 * t218 + t235 * t409 + t284 * t240 + t281 * t241 - t253 * t393 + t254 * t392 + t466 + (t293 + t467) * t271;
t404 = t274 * t283;
t412 = t272 * t280;
t213 = -t273 * t412 - t404;
t405 = t274 * t280;
t411 = t272 * t283;
t214 = t273 * t411 - t405;
t417 = t271 * t272;
t151 = Icges(6,5) * t214 + Icges(6,6) * t213 + Icges(6,3) * t417;
t153 = Icges(6,4) * t214 + Icges(6,2) * t213 + Icges(6,6) * t417;
t155 = Icges(6,1) * t214 + Icges(6,4) * t213 + Icges(6,5) * t417;
t331 = -t153 * t280 + t155 * t283;
t354 = -qJD(5) * t273 + t277;
t302 = t271 * t401 + t283 * t354;
t408 = t273 * t277;
t355 = -qJD(5) + t408;
t139 = t272 * t302 - t355 * t405;
t301 = -t271 * t400 + t280 * t354;
t140 = t272 * t301 + t355 * t404;
t307 = t271 * t406 + t272 * t409;
t77 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t307;
t79 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t307;
t81 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t307;
t21 = (t276 * t331 - t77) * t273 + (t151 * t276 - t280 * t79 + t283 * t81 + (-t153 * t283 - t155 * t280) * qJD(5)) * t271;
t34 = t139 * t201 + t140 * t204 + t159 * t417 + t162 * t213 + t165 * t214 + t198 * t307;
t463 = t21 + t34;
t215 = -t273 * t405 + t411;
t216 = t273 * t404 + t412;
t416 = t271 * t274;
t152 = Icges(6,5) * t216 + Icges(6,6) * t215 + Icges(6,3) * t416;
t154 = Icges(6,4) * t216 + Icges(6,2) * t215 + Icges(6,6) * t416;
t156 = Icges(6,1) * t216 + Icges(6,4) * t215 + Icges(6,5) * t416;
t330 = -t154 * t280 + t156 * t283;
t137 = t274 * t302 + t355 * t412;
t138 = t274 * t301 - t355 * t411;
t407 = t274 * t276;
t380 = t273 * t407;
t306 = t380 - t384;
t76 = Icges(6,5) * t138 + Icges(6,6) * t137 + Icges(6,3) * t306;
t78 = Icges(6,4) * t138 + Icges(6,2) * t137 + Icges(6,6) * t306;
t80 = Icges(6,1) * t138 + Icges(6,4) * t137 + Icges(6,5) * t306;
t22 = (t276 * t330 - t76) * t273 + (t152 * t276 - t280 * t78 + t283 * t80 + (-t154 * t283 - t156 * t280) * qJD(5)) * t271;
t33 = t137 * t201 + t138 * t204 + t159 * t416 + t162 * t215 + t165 * t216 + t198 * t306;
t462 = t22 + t33;
t68 = -t151 * t273 + t271 * t331;
t94 = t198 * t417 + t201 * t213 + t204 * t214;
t461 = t68 + t94;
t69 = -t152 * t273 + t271 * t330;
t95 = t198 * t416 + t201 * t215 + t204 * t216;
t460 = t69 + t95;
t394 = t272 ^ 2 + t274 ^ 2;
t399 = t277 * t281;
t457 = -t272 * t392 - t274 * t399;
t311 = t337 * t274;
t203 = Icges(4,6) * t272 + t311;
t313 = t340 * t274;
t206 = Icges(4,5) * t272 + t313;
t321 = t203 * t281 - t206 * t284;
t456 = t272 * t321;
t310 = t336 * t274;
t185 = Icges(5,6) * t272 + t310;
t312 = t339 * t274;
t187 = Icges(5,5) * t272 + t312;
t325 = t185 * t271 - t187 * t273;
t455 = t272 * t325;
t202 = -Icges(4,6) * t274 + t272 * t337;
t205 = -Icges(4,5) * t274 + t272 * t340;
t323 = t202 * t281 - t205 * t284;
t454 = t274 * t323;
t184 = -Icges(5,6) * t274 + t272 * t336;
t186 = -Icges(5,5) * t274 + t272 * t339;
t326 = t184 * t271 - t186 * t273;
t453 = t274 * t326;
t264 = t272 * rSges(5,3);
t434 = rSges(5,1) * t273;
t452 = t274 * t434 + t264;
t403 = t274 * t284;
t317 = rSges(4,1) * t403 + t272 * rSges(4,3);
t83 = t140 * rSges(6,1) + t139 * rSges(6,2) + rSges(6,3) * t307;
t320 = t234 * t271 - t235 * t273;
t333 = Icges(5,5) * t273 - Icges(5,6) * t271;
t451 = t333 * t276 + t277 * t320;
t319 = t253 * t281 - t254 * t284;
t334 = Icges(4,5) * t284 - Icges(4,6) * t281;
t450 = t334 * qJD(3) + t277 * t319;
t375 = t272 * t393;
t170 = -rSges(4,1) * t375 + rSges(4,2) * t457 + t277 * t317;
t448 = 2 * m(3);
t447 = 2 * m(4);
t446 = 2 * m(5);
t445 = 2 * m(6);
t444 = t272 / 0.2e1;
t443 = -t274 / 0.2e1;
t433 = rSges(4,2) * t281;
t244 = (rSges(4,1) * t284 - t433) * qJD(3);
t442 = m(4) * t244;
t258 = rSges(4,1) * t281 + rSges(4,2) * t284;
t441 = m(4) * t258;
t346 = rSges(6,1) * t283 - rSges(6,2) * t280;
t168 = t346 * t409 + (rSges(6,3) * t276 + (-rSges(6,1) * t280 - rSges(6,2) * t283) * qJD(5)) * t271;
t440 = m(6) * t168;
t207 = -rSges(6,3) * t273 + t271 * t346;
t439 = m(6) * t207;
t282 = sin(qJ(1));
t438 = pkin(1) * t282;
t437 = pkin(2) * t281;
t436 = pkin(2) * t284;
t435 = rSges(5,1) * t271;
t432 = pkin(1) * qJD(1);
t425 = t168 * t272;
t424 = t168 * t274;
t182 = -Icges(5,3) * t274 + t272 * t333;
t423 = t182 * t277;
t308 = t333 * t274;
t183 = Icges(5,3) * t272 + t308;
t422 = t183 * t277;
t220 = (-rSges(5,2) * t271 + t434) * t276;
t419 = t220 * t272;
t418 = t235 * t276;
t414 = t272 * t276;
t410 = t272 * t284;
t347 = -rSges(6,1) * t214 - rSges(6,2) * t213;
t157 = rSges(6,3) * t417 - t347;
t158 = t216 * rSges(6,1) + t215 * rSges(6,2) + rSges(6,3) * t416;
t96 = t272 * t157 + t274 * t158;
t396 = rSges(5,2) * t417 + t274 * rSges(5,3);
t190 = t272 * t434 - t396;
t191 = -rSges(5,2) * t416 + t452;
t134 = t272 * t190 + t274 * t191;
t397 = rSges(5,2) * t384 + rSges(5,3) * t406;
t395 = t394 * t436;
t263 = pkin(2) * t403;
t390 = t282 * t432;
t285 = cos(qJ(1));
t389 = t285 * t432;
t388 = pkin(2) * t393;
t387 = pkin(2) * t392;
t385 = t207 * t413;
t381 = t272 * t399;
t377 = -rSges(5,2) * t307 - t414 * t435;
t373 = t417 / 0.2e1;
t372 = t416 / 0.2e1;
t371 = t413 / 0.2e1;
t368 = -t207 - t437;
t236 = rSges(5,2) * t273 + t435;
t367 = -t236 - t437;
t56 = t152 * t416 + t154 * t215 + t156 * t216;
t366 = -t137 * t153 - t138 * t155 - t151 * t306 - t215 * t79 - t216 * t81 + t277 * t56 - t416 * t77;
t55 = t151 * t416 + t153 * t215 + t155 * t216;
t365 = t137 * t154 + t138 * t156 + t152 * t306 + t215 * t78 + t216 * t80 + t277 * t55 + t416 * t76;
t54 = t152 * t417 + t154 * t213 + t156 * t214;
t364 = -t139 * t153 - t140 * t155 - t151 * t307 - t213 * t79 - t214 * t81 + t277 * t54 - t417 * t77;
t53 = t151 * t417 + t153 * t213 + t155 * t214;
t363 = t139 * t154 + t140 * t156 + t152 * t307 + t213 * t78 + t214 * t80 + t277 * t53 + t417 * t76;
t238 = t274 * rSges(3,1) - rSges(3,2) * t272;
t104 = -t198 * t273 + (-t201 * t280 + t204 * t283) * t271;
t362 = (t271 * t293 + t466) * t273 - t104 * t415;
t298 = Icges(5,6) * t277 - t402;
t128 = t274 * t298 - t336 * t413;
t361 = t187 * t276 + t128;
t129 = t272 * t298 + t277 * t310;
t360 = t186 * t276 + t129;
t299 = Icges(5,5) * t277 - t418;
t130 = t274 * t299 - t339 * t413;
t359 = -t185 * t276 + t130;
t131 = t272 * t299 + t277 * t312;
t358 = -t184 * t276 + t131;
t305 = -t274 * t393 - t277 * t410;
t169 = rSges(4,1) * t305 + rSges(4,3) * t406 + (-t274 * t392 + t381) * rSges(4,2);
t208 = rSges(4,1) * t410 - t274 * rSges(4,3) - t272 * t433;
t357 = t208 * t277 + t169;
t209 = -t274 * t433 + t317;
t356 = -t209 * t277 + t170;
t142 = t263 + t158;
t222 = -rSges(3,1) * t406 + rSges(3,2) * t413;
t349 = -t434 - t436;
t237 = -rSges(3,1) * t272 - rSges(3,2) * t274;
t345 = t272 * t53 + t274 * t54;
t344 = t272 * t55 + t274 * t56;
t343 = t272 * t69 - t274 * t68;
t342 = t272 * t68 + t274 * t69;
t100 = -t183 * t274 - t455;
t101 = t182 * t272 - t453;
t102 = t183 * t272 - t274 * t325;
t233 = Icges(5,5) * t271 + Icges(5,6) * t273;
t297 = Icges(5,3) * t277 - t233 * t276;
t126 = t274 * t297 - t333 * t413;
t127 = t272 * t297 + t277 * t308;
t42 = t272 * t54 - t274 * t53;
t43 = t272 * t56 - t274 * t55;
t8 = t272 * t365 + t274 * t366;
t99 = -t182 * t274 - t272 * t326;
t341 = t43 * t406 + t42 * t413 + (-t101 * t406 - t99 * t413) * t274 + (t8 + (t102 * t277 + (t129 * t271 - t131 * t273 + t184 * t409 + t186 * t415 - t423) * t274) * t274 + t100 * t413 + t102 * t406 + ((t101 + t455) * t277 + (t422 - t127 + (t186 * t277 + t359) * t273 + (-t184 * t277 - t361) * t271) * t274 + t272 * t126) * t272) * t272;
t252 = Icges(4,5) * t281 + Icges(4,6) * t284;
t329 = t157 * t274 - t158 * t272;
t324 = t202 * t284 + t205 * t281;
t322 = t203 * t284 + t206 * t281;
t318 = t394 * t388;
t315 = -t207 * t406 - t425;
t314 = t349 * t272;
t82 = t138 * rSges(6,1) + t137 * rSges(6,2) + rSges(6,3) * t306;
t39 = t157 * t406 - t158 * t413 + t272 * t83 + t274 * t82;
t221 = t237 * t277;
t176 = t191 + t263;
t59 = -t191 * t413 + t272 * (t277 * t452 + t377) + t274 * (-rSges(5,2) * t380 + (-t271 * t407 - t272 * t408) * rSges(5,1) + t397) + t190 * t406;
t309 = t334 * t274;
t304 = t457 * pkin(2);
t12 = (t274 * t127 + (t100 + t453) * t277) * t274 + (t99 * t277 + (-t128 * t271 + t130 * t273 - t185 * t409 - t187 * t415 + t422) * t272 + (-t423 - t126 + (t187 * t277 - t358) * t273 + (-t185 * t277 + t360) * t271) * t274) * t272;
t9 = t272 * t363 + t274 * t364;
t303 = (-t12 - t9) * t274 + t341;
t175 = t314 + t396;
t26 = t271 * t345 - t273 * t94;
t27 = t271 * t344 - t273 * t95;
t3 = (t276 * t344 - t33) * t273 + (-t272 * t366 + t274 * t365 + t276 * t95) * t271;
t4 = (t276 * t345 - t34) * t273 + (-t272 * t364 + t274 * t363 + t276 * t94) * t271;
t300 = t3 * t444 + t8 * t372 + t9 * t373 + t4 * t443 - t273 * ((t277 * t69 - t21) * t274 + (t277 * t68 + t22) * t272) / 0.2e1 + t26 * t371 + t27 * t369 + t343 * t415 / 0.2e1 + t465 * t43 + t464 * t42;
t296 = Icges(4,5) * t277 - qJD(3) * t254;
t295 = Icges(4,6) * t277 - qJD(3) * t253;
t294 = Icges(4,3) * t277 - qJD(3) * t252;
t141 = (-rSges(6,3) * t271 - t436) * t272 + t347;
t291 = t462 * t372 + t463 * t373 + t460 * t465 + t461 * t464 - t362;
t290 = t233 * t277 + t467 * t273 + (-t218 - t418) * t271;
t289 = (t271 * t359 + t272 * t451 + t273 * t361 + t290 * t274 + t462) * t444 + (t271 * t358 + t290 * t272 + t273 * t360 - t274 * t451 + t463) * t443 + (t184 * t273 + t186 * t271 - t233 * t274 - t272 * t320 + t461) * t371 + (t185 * t273 + t187 * t271 + t233 * t272 - t274 * t320 + t460) * t369;
t288 = -t240 * t281 + t241 * t284 + t252 * t277 + (-t253 * t284 - t254 * t281) * qJD(3);
t251 = pkin(2) * t375;
t114 = t251 + (t274 * t349 - t264) * t277 - t377;
t70 = pkin(2) * t305 + t82;
t71 = -t263 * t277 + t251 - t83;
t287 = t289 + (-qJD(3) * t321 + (t274 * t295 - t337 * t413) * t284 + (t274 * t296 - t340 * t413) * t281 + t450 * t272 + t288 * t274) * t444 + (-qJD(3) * t323 + (t272 * t295 + t277 * t311) * t284 + (t272 * t296 + t277 * t313) * t281 + t288 * t272 - t274 * t450) * t443 + (-t252 * t274 - t272 * t319 + t324) * t371 + (t252 * t272 - t274 * t319 + t322) * t369;
t113 = t277 * t314 + (-t236 * t276 - t388) * t274 + t397;
t275 = t285 * pkin(1);
t245 = pkin(2) * t381;
t224 = t238 + t275;
t223 = t237 - t438;
t200 = Icges(4,3) * t272 + t309;
t199 = -Icges(4,3) * t274 + t272 * t334;
t195 = t222 - t389;
t194 = t221 - t390;
t193 = t367 * t274;
t192 = t367 * t272;
t189 = t209 + t275;
t188 = -t208 - t438;
t174 = t368 * t274;
t173 = t368 * t272;
t172 = t176 + t275;
t171 = t175 - t438;
t161 = t272 * t294 + t277 * t309;
t160 = t274 * t294 - t334 * t413;
t148 = -t170 - t389;
t147 = t169 - t390;
t125 = t275 + t142;
t124 = t141 - t438;
t123 = -t236 * t406 + t304 - t419;
t122 = t236 * t413 + t245 + (-t220 - t387) * t274;
t115 = t395 + t134;
t112 = t114 - t389;
t111 = t113 - t390;
t110 = -t158 * t273 - t207 * t416;
t109 = t157 * t273 + t207 * t417;
t108 = t200 * t272 - t274 * t321;
t107 = t199 * t272 - t454;
t106 = -t200 * t274 - t456;
t105 = -t199 * t274 - t272 * t323;
t98 = t304 + t315;
t97 = t385 + t245 + (-t168 - t387) * t274;
t93 = t329 * t271;
t90 = t395 + t96;
t67 = t71 - t389;
t66 = t70 - t390;
t50 = -t318 + t59;
t46 = (t207 * t414 + t83) * t273 + (-t157 * t276 - t315) * t271;
t45 = (-t207 * t407 - t82) * t273 + (t158 * t276 + t385 - t424) * t271;
t30 = -t318 + t39;
t23 = t329 * t409 + ((-t158 * t277 + t83) * t274 + (-t157 * t277 - t82) * t272) * t271;
t1 = [(t124 * t67 + t125 * t66) * t445 + (t111 * t172 + t112 * t171) * t446 + (t147 * t189 + t148 * t188) * t447 + (t194 * t224 + t195 * t223) * t448 + t286; m(6) * (t124 * t71 + t125 * t70 + t141 * t67 + t142 * t66) + m(5) * (t111 * t176 + t112 * t175 + t113 * t172 + t114 * t171) + m(4) * (t147 * t209 - t148 * t208 + t169 * t189 - t170 * t188) + m(3) * (t194 * t238 + t195 * t237 + t221 * t224 + t222 * t223) + t286; (t141 * t71 + t142 * t70) * t445 + (t113 * t176 + t114 * t175) * t446 + (t169 * t209 + t170 * t208) * t447 + (t221 * t238 + t222 * t237) * t448 + t286; t287 + ((-t189 * t277 - t148) * t274 + (t188 * t277 - t147) * t272) * t441 + (-t188 * t274 - t189 * t272) * t442 + m(5) * (t111 * t192 + t112 * t193 + t122 * t171 + t123 * t172) + m(6) * (t124 * t97 + t125 * t98 + t173 * t66 + t174 * t67); (t208 * t274 - t209 * t272) * t442 + t287 + (-t272 * t357 + t274 * t356) * t441 + m(5) * (t113 * t192 + t114 * t193 + t122 * t175 + t123 * t176) + m(6) * (t141 * t97 + t142 * t98 + t173 * t70 + t174 * t71); (t173 * t98 + t174 * t97 + t30 * t90) * t445 - t274 * t9 - t274 * t12 + (t115 * t50 + t122 * t193 + t123 * t192) * t446 + (-t107 * t274 + t108 * t272) * t406 + t272 * ((t272 * t160 + (t107 + t456) * t277) * t272 + (t108 * t277 + (t202 * t392 + t205 * t393) * t274 + (-t322 * qJD(3) - t277 * t323 - t161) * t272) * t274) + (-t105 * t274 + t106 * t272) * t413 - t274 * ((t274 * t161 + (t106 + t454) * t277) * t274 + (t105 * t277 + (-t203 * t392 - t206 * t393) * t272 + (t324 * qJD(3) - t277 * t321 - t160) * t274) * t272) + ((t208 * t272 + t209 * t274) * (t272 * t356 + t274 * t357) + t394 * t258 * t244) * t447 + t341; (-t124 * t274 - t125 * t272) * t440 + ((-t125 * t277 - t67) * t274 + (t124 * t277 - t66) * t272) * t439 + m(5) * ((-t171 * t274 - t172 * t272) * t220 + ((-t172 * t277 - t112) * t274 + (t171 * t277 - t111) * t272) * t236) + t289; ((-t142 * t277 - t71) * t274 + (t141 * t277 - t70) * t272) * t439 + m(5) * ((-t175 * t274 - t176 * t272) * t220 + ((-t176 * t277 - t114) * t274 + (t175 * t277 - t113) * t272) * t236) + t289 + (-t141 * t274 - t142 * t272) * t440; m(6) * (-t173 * t425 - t174 * t424 + t30 * t96 + t39 * t90) + ((-t173 * t277 - t97) * t274 + (t174 * t277 - t98) * t272) * t439 + t303 + (-t193 * t220 * t274 + t115 * t59 + t134 * t50 - t192 * t419 + ((-t192 * t277 - t122) * t274 + (t193 * t277 - t123) * t272) * t236) * m(5); (t220 * t236 * t394 + t134 * t59) * t446 + (t168 * t207 * t394 + t39 * t96) * t445 + t303; m(6) * (t109 * t67 + t110 * t66 + t124 * t46 + t125 * t45) + t291; m(6) * (t109 * t71 + t110 * t70 + t141 * t46 + t142 * t45) + t291; t300 + m(6) * (t109 * t97 + t110 * t98 + t173 * t45 + t174 * t46 + t23 * t90 + t30 * t93); t300 + m(6) * (t23 * t96 + t39 * t93 + (-t109 * t274 - t110 * t272) * t168 + ((-t110 * t277 - t46) * t274 + (t109 * t277 - t45) * t272) * t207); (t109 * t46 + t110 * t45 + t23 * t93) * t445 + ((t272 * t26 + t274 * t27 - t273 * t342) * t276 + t362) * t273 + (t274 * t3 + t272 * t4 + t342 * t415 + (-t104 * t276 - t21 * t272 - t22 * t274) * t273 + (t274 * t26 - t272 * t27 + t273 * t343) * t277) * t271;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11); t1(2) t1(3) t1(5) t1(8) t1(12); t1(4) t1(5) t1(6) t1(9) t1(13); t1(7) t1(8) t1(9) t1(10) t1(14); t1(11) t1(12) t1(13) t1(14) t1(15);];
Mq  = res;

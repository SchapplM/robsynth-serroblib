% Calculate time derivative of joint inertia matrix for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:35
% EndTime: 2019-03-08 21:13:12
% DurationCPUTime: 24.37s
% Computational Cost: add. (63925->1158), mult. (182626->1596), div. (0->0), fcn. (212991->12), ass. (0->469)
t414 = sin(pkin(10));
t416 = cos(pkin(10));
t421 = cos(qJ(2));
t417 = cos(pkin(6));
t419 = sin(qJ(2));
t522 = t417 * t419;
t401 = t414 * t421 + t416 * t522;
t415 = sin(pkin(6));
t529 = sin(qJ(3));
t476 = t415 * t529;
t530 = cos(qJ(3));
t378 = t401 * t530 - t416 * t476;
t413 = sin(pkin(11));
t521 = t417 * t421;
t430 = -t414 * t419 + t416 * t521;
t528 = cos(pkin(11));
t335 = t378 * t413 + t430 * t528;
t336 = t378 * t528 - t413 * t430;
t477 = t415 * t530;
t424 = -t401 * t529 - t416 * t477;
t234 = Icges(6,5) * t336 - Icges(6,6) * t424 + Icges(6,3) * t335;
t240 = Icges(5,4) * t336 - Icges(5,2) * t335 - Icges(5,6) * t424;
t560 = t234 - t240;
t488 = t414 * t522;
t403 = t416 * t421 - t488;
t380 = t403 * t530 + t414 * t476;
t402 = t414 * t521 + t416 * t419;
t337 = t380 * t413 - t402 * t528;
t338 = t380 * t528 + t402 * t413;
t425 = -t403 * t529 + t414 * t477;
t235 = Icges(6,5) * t338 - Icges(6,6) * t425 + Icges(6,3) * t337;
t241 = Icges(5,4) * t338 - Icges(5,2) * t337 - Icges(5,6) * t425;
t559 = t235 - t241;
t242 = Icges(6,1) * t336 - Icges(6,4) * t424 + Icges(6,5) * t335;
t244 = Icges(5,1) * t336 - Icges(5,4) * t335 - Icges(5,5) * t424;
t558 = t242 + t244;
t243 = Icges(6,1) * t338 - Icges(6,4) * t425 + Icges(6,5) * t337;
t245 = Icges(5,1) * t338 - Icges(5,4) * t337 - Icges(5,5) * t425;
t557 = t243 + t245;
t236 = Icges(5,5) * t336 - Icges(5,6) * t335 - Icges(5,3) * t424;
t238 = Icges(6,4) * t336 - Icges(6,2) * t424 + Icges(6,6) * t335;
t294 = Icges(4,4) * t378 + Icges(4,2) * t424 - Icges(4,6) * t430;
t556 = -t236 - t238 + t294;
t237 = Icges(5,5) * t338 - Icges(5,6) * t337 - Icges(5,3) * t425;
t239 = Icges(6,4) * t338 - Icges(6,2) * t425 + Icges(6,6) * t337;
t295 = Icges(4,4) * t380 + Icges(4,2) * t425 + Icges(4,6) * t402;
t555 = -t237 - t239 + t295;
t292 = Icges(4,5) * t378 + Icges(4,6) * t424 - Icges(4,3) * t430;
t296 = Icges(4,1) * t378 + Icges(4,4) * t424 - Icges(4,5) * t430;
t552 = -t292 * t430 + t296 * t378 + t335 * t560 + t558 * t336 + t556 * t424;
t293 = Icges(4,5) * t380 + Icges(4,6) * t425 + Icges(4,3) * t402;
t297 = Icges(4,1) * t380 + Icges(4,4) * t425 + Icges(4,5) * t402;
t551 = -t293 * t430 + t297 * t378 + t335 * t559 + t336 * t557 + t424 * t555;
t550 = t292 * t402 + t296 * t380 + t337 * t560 + t558 * t338 + t556 * t425;
t549 = t293 * t402 + t297 * t380 + t337 * t559 + t338 * t557 + t425 * t555;
t405 = t417 * t529 + t419 * t477;
t473 = t415 * t528;
t375 = t405 * t413 + t421 * t473;
t523 = t415 * t421;
t376 = t405 * t528 - t413 * t523;
t404 = -t417 * t530 + t419 * t476;
t302 = Icges(6,5) * t376 + Icges(6,6) * t404 + Icges(6,3) * t375;
t304 = Icges(6,4) * t376 + Icges(6,2) * t404 + Icges(6,6) * t375;
t306 = Icges(6,1) * t376 + Icges(6,4) * t404 + Icges(6,5) * t375;
t144 = t302 * t335 - t304 * t424 + t306 * t336;
t303 = Icges(5,5) * t376 - Icges(5,6) * t375 + Icges(5,3) * t404;
t305 = Icges(5,4) * t376 - Icges(5,2) * t375 + Icges(5,6) * t404;
t307 = Icges(5,1) * t376 - Icges(5,4) * t375 + Icges(5,5) * t404;
t145 = -t303 * t424 - t305 * t335 + t307 * t336;
t355 = Icges(4,5) * t405 - Icges(4,6) * t404 - Icges(4,3) * t523;
t356 = Icges(4,4) * t405 - Icges(4,2) * t404 - Icges(4,6) * t523;
t357 = Icges(4,1) * t405 - Icges(4,4) * t404 - Icges(4,5) * t523;
t209 = -t355 * t430 + t356 * t424 + t357 * t378;
t554 = t144 + t145 + t209;
t146 = t302 * t337 - t304 * t425 + t306 * t338;
t147 = -t303 * t425 - t305 * t337 + t307 * t338;
t210 = t355 * t402 + t356 * t425 + t357 * t380;
t553 = t146 + t147 + t210;
t525 = t414 * t415;
t524 = t415 * t416;
t391 = t430 * qJD(2);
t392 = t401 * qJD(2);
t369 = pkin(2) * t391 + pkin(8) * t392;
t393 = t402 * qJD(2);
t492 = qJD(2) * t421;
t394 = -qJD(2) * t488 + t416 * t492;
t370 = -pkin(2) * t393 + pkin(8) * t394;
t496 = t369 * t525 + t370 * t524;
t474 = t415 * t492;
t382 = t405 * qJD(3) + t529 * t474;
t489 = m(6) / 0.2e1 + m(7) / 0.2e1;
t383 = -qJD(3) * t404 + t474 * t530;
t493 = qJD(2) * t419;
t475 = t415 * t493;
t316 = Icges(4,5) * t383 - Icges(4,6) * t382 + Icges(4,3) * t475;
t317 = Icges(4,4) * t383 - Icges(4,2) * t382 + Icges(4,6) * t475;
t318 = Icges(4,1) * t383 - Icges(4,4) * t382 + Icges(4,5) * t475;
t331 = qJD(3) * t378 + t391 * t529;
t332 = qJD(3) * t424 + t391 * t530;
t130 = -t316 * t430 + t317 * t424 + t318 * t378 - t331 * t356 + t332 * t357 + t355 * t392;
t418 = sin(qJ(6));
t420 = cos(qJ(6));
t325 = t375 * t420 - t376 * t418;
t326 = t375 * t418 + t376 * t420;
t220 = Icges(7,5) * t326 + Icges(7,6) * t325 - Icges(7,3) * t404;
t221 = Icges(7,4) * t326 + Icges(7,2) * t325 - Icges(7,6) * t404;
t222 = Icges(7,1) * t326 + Icges(7,4) * t325 - Icges(7,5) * t404;
t279 = t335 * t420 - t336 * t418;
t280 = t335 * t418 + t336 * t420;
t103 = t220 * t424 + t221 * t279 + t222 * t280;
t298 = t332 * t413 - t392 * t528;
t299 = t332 * t528 + t392 * t413;
t164 = -qJD(6) * t280 + t298 * t420 - t299 * t418;
t165 = qJD(6) * t279 + t298 * t418 + t299 * t420;
t107 = Icges(7,5) * t165 + Icges(7,6) * t164 - Icges(7,3) * t331;
t109 = Icges(7,4) * t165 + Icges(7,2) * t164 - Icges(7,6) * t331;
t111 = Icges(7,1) * t165 + Icges(7,4) * t164 - Icges(7,5) * t331;
t172 = Icges(7,5) * t280 + Icges(7,6) * t279 + Icges(7,3) * t424;
t174 = Icges(7,4) * t280 + Icges(7,2) * t279 + Icges(7,6) * t424;
t176 = Icges(7,1) * t280 + Icges(7,4) * t279 + Icges(7,5) * t424;
t21 = t107 * t424 + t109 * t279 + t111 * t280 + t164 * t174 + t165 * t176 - t172 * t331;
t282 = t337 * t418 + t338 * t420;
t334 = qJD(3) * t425 - t393 * t530;
t300 = t334 * t413 - t394 * t528;
t301 = t334 * t528 + t394 * t413;
t166 = -qJD(6) * t282 + t300 * t420 - t301 * t418;
t281 = t337 * t420 - t338 * t418;
t167 = qJD(6) * t281 + t300 * t418 + t301 * t420;
t333 = qJD(3) * t380 - t393 * t529;
t108 = Icges(7,5) * t167 + Icges(7,6) * t166 - Icges(7,3) * t333;
t110 = Icges(7,4) * t167 + Icges(7,2) * t166 - Icges(7,6) * t333;
t112 = Icges(7,1) * t167 + Icges(7,4) * t166 - Icges(7,5) * t333;
t173 = Icges(7,5) * t282 + Icges(7,6) * t281 + Icges(7,3) * t425;
t175 = Icges(7,4) * t282 + Icges(7,2) * t281 + Icges(7,6) * t425;
t177 = Icges(7,1) * t282 + Icges(7,4) * t281 + Icges(7,5) * t425;
t22 = t108 * t424 + t110 * t279 + t112 * t280 + t164 * t175 + t165 * t177 - t173 * t331;
t353 = t383 * t413 - t473 * t493;
t354 = t383 * t528 + t413 * t475;
t215 = -qJD(6) * t326 + t353 * t420 - t354 * t418;
t216 = qJD(6) * t325 + t353 * t418 + t354 * t420;
t153 = Icges(7,5) * t216 + Icges(7,6) * t215 - Icges(7,3) * t382;
t154 = Icges(7,4) * t216 + Icges(7,2) * t215 - Icges(7,6) * t382;
t155 = Icges(7,1) * t216 + Icges(7,4) * t215 - Icges(7,5) * t382;
t43 = t153 * t424 + t154 * t279 + t155 * t280 + t164 * t221 + t165 * t222 - t220 * t331;
t74 = t172 * t424 + t174 * t279 + t176 * t280;
t75 = t173 * t424 + t175 * t279 + t177 * t280;
t3 = -t21 * t430 + t22 * t402 + t74 * t392 + t75 * t394 + (t103 * t493 - t421 * t43) * t415;
t190 = Icges(6,5) * t299 + Icges(6,6) * t331 + Icges(6,3) * t298;
t194 = Icges(6,4) * t299 + Icges(6,2) * t331 + Icges(6,6) * t298;
t198 = Icges(6,1) * t299 + Icges(6,4) * t331 + Icges(6,5) * t298;
t58 = t190 * t335 - t194 * t424 + t198 * t336 + t234 * t298 + t238 * t331 + t242 * t299;
t191 = Icges(6,5) * t301 + Icges(6,6) * t333 + Icges(6,3) * t300;
t195 = Icges(6,4) * t301 + Icges(6,2) * t333 + Icges(6,6) * t300;
t199 = Icges(6,1) * t301 + Icges(6,4) * t333 + Icges(6,5) * t300;
t59 = t191 * t335 - t195 * t424 + t199 * t336 + t235 * t298 + t239 * t331 + t243 * t299;
t192 = Icges(5,5) * t299 - Icges(5,6) * t298 + Icges(5,3) * t331;
t196 = Icges(5,4) * t299 - Icges(5,2) * t298 + Icges(5,6) * t331;
t200 = Icges(5,1) * t299 - Icges(5,4) * t298 + Icges(5,5) * t331;
t60 = -t192 * t424 - t196 * t335 + t200 * t336 + t236 * t331 - t240 * t298 + t244 * t299;
t193 = Icges(5,5) * t301 - Icges(5,6) * t300 + Icges(5,3) * t333;
t197 = Icges(5,4) * t301 - Icges(5,2) * t300 + Icges(5,6) * t333;
t201 = Icges(5,1) * t301 - Icges(5,4) * t300 + Icges(5,5) * t333;
t61 = -t193 * t424 - t197 * t335 + t201 * t336 + t237 * t331 - t241 * t298 + t245 * t299;
t266 = Icges(6,5) * t354 + Icges(6,6) * t382 + Icges(6,3) * t353;
t268 = Icges(6,4) * t354 + Icges(6,2) * t382 + Icges(6,6) * t353;
t270 = Icges(6,1) * t354 + Icges(6,4) * t382 + Icges(6,5) * t353;
t82 = t266 * t335 - t268 * t424 + t270 * t336 + t298 * t302 + t299 * t306 + t304 * t331;
t267 = Icges(5,5) * t354 - Icges(5,6) * t353 + Icges(5,3) * t382;
t269 = Icges(5,4) * t354 - Icges(5,2) * t353 + Icges(5,6) * t382;
t271 = Icges(5,1) * t354 - Icges(5,4) * t353 + Icges(5,5) * t382;
t83 = -t267 * t424 - t269 * t335 + t271 * t336 - t298 * t305 + t299 * t307 + t303 * t331;
t253 = Icges(4,5) * t332 - Icges(4,6) * t331 + Icges(4,3) * t392;
t255 = Icges(4,4) * t332 - Icges(4,2) * t331 + Icges(4,6) * t392;
t257 = Icges(4,1) * t332 - Icges(4,4) * t331 + Icges(4,5) * t392;
t93 = -t253 * t430 + t255 * t424 + t257 * t378 + t292 * t392 - t294 * t331 + t296 * t332;
t254 = Icges(4,5) * t334 - Icges(4,6) * t333 + Icges(4,3) * t394;
t256 = Icges(4,4) * t334 - Icges(4,2) * t333 + Icges(4,6) * t394;
t258 = Icges(4,1) * t334 - Icges(4,4) * t333 + Icges(4,5) * t394;
t94 = -t254 * t430 + t256 * t424 + t258 * t378 + t293 * t392 - t295 * t331 + t297 * t332;
t548 = t3 + (-t60 - t58 - t93) * t430 + (t554 * t493 + (-t130 - t82 - t83) * t421) * t415 + (t61 + t59 + t94) * t402 + t551 * t394 + t552 * t392;
t131 = t316 * t402 + t317 * t425 + t318 * t380 - t333 * t356 + t334 * t357 + t355 * t394;
t104 = t220 * t425 + t221 * t281 + t222 * t282;
t23 = t107 * t425 + t109 * t281 + t111 * t282 + t166 * t174 + t167 * t176 - t172 * t333;
t24 = t108 * t425 + t110 * t281 + t112 * t282 + t166 * t175 + t167 * t177 - t173 * t333;
t44 = t153 * t425 + t154 * t281 + t155 * t282 + t166 * t221 + t167 * t222 - t220 * t333;
t76 = t172 * t425 + t174 * t281 + t176 * t282;
t77 = t173 * t425 + t175 * t281 + t177 * t282;
t4 = -t23 * t430 + t24 * t402 + t76 * t392 + t77 * t394 + (t104 * t493 - t421 * t44) * t415;
t62 = t190 * t337 - t194 * t425 + t198 * t338 + t234 * t300 + t238 * t333 + t242 * t301;
t63 = t191 * t337 - t195 * t425 + t199 * t338 + t235 * t300 + t239 * t333 + t243 * t301;
t64 = -t192 * t425 - t196 * t337 + t200 * t338 + t236 * t333 - t240 * t300 + t244 * t301;
t65 = -t193 * t425 - t197 * t337 + t201 * t338 + t237 * t333 - t241 * t300 + t245 * t301;
t84 = t266 * t337 - t268 * t425 + t270 * t338 + t300 * t302 + t301 * t306 + t304 * t333;
t85 = -t267 * t425 - t269 * t337 + t271 * t338 - t300 * t305 + t301 * t307 + t303 * t333;
t95 = t253 * t402 + t255 * t425 + t257 * t380 + t292 * t394 - t294 * t333 + t296 * t334;
t96 = t254 * t402 + t256 * t425 + t258 * t380 + t293 * t394 - t295 * t333 + t297 * t334;
t547 = t4 + (-t64 - t62 - t95) * t430 + (t553 * t493 + (-t131 - t84 - t85) * t421) * t415 + (t65 + t63 + t96) * t402 + t549 * t394 + t550 * t392;
t204 = rSges(6,1) * t301 + rSges(6,2) * t333 + rSges(6,3) * t300;
t205 = rSges(5,1) * t301 - rSges(5,2) * t300 + rSges(5,3) * t333;
t260 = rSges(4,1) * t334 - rSges(4,2) * t333 + rSges(4,3) * t394;
t114 = rSges(7,1) * t167 + rSges(7,2) * t166 - rSges(7,3) * t333;
t519 = pkin(5) * t301 - pkin(9) * t333 + t114;
t546 = -m(4) * t260 - m(5) * t205 - m(6) * t204 - m(7) * t519;
t545 = -0.2e1 * t392;
t544 = 0.2e1 * t430;
t543 = m(5) / 0.2e1;
t540 = -t331 / 0.2e1;
t539 = -t333 / 0.2e1;
t538 = t424 / 0.2e1;
t537 = t425 / 0.2e1;
t536 = -t382 / 0.2e1;
t535 = t392 / 0.2e1;
t534 = t394 / 0.2e1;
t533 = -t404 / 0.2e1;
t532 = t414 / 0.2e1;
t531 = -t416 / 0.2e1;
t527 = Icges(3,4) * t419;
t526 = Icges(3,4) * t421;
t113 = rSges(7,1) * t165 + rSges(7,2) * t164 - rSges(7,3) * t331;
t520 = pkin(5) * t299 - pkin(9) * t331 + t113;
t156 = rSges(7,1) * t216 + rSges(7,2) * t215 - rSges(7,3) * t382;
t518 = pkin(5) * t354 - pkin(9) * t382 + t156;
t178 = rSges(7,1) * t280 + rSges(7,2) * t279 + rSges(7,3) * t424;
t517 = pkin(5) * t336 + pkin(9) * t424 + t178;
t179 = rSges(7,1) * t282 + rSges(7,2) * t281 + rSges(7,3) * t425;
t516 = pkin(5) * t338 + pkin(9) * t425 + t179;
t233 = pkin(3) * t334 + qJ(4) * t333 - qJD(4) * t425;
t515 = -t205 - t233;
t207 = pkin(4) * t301 + qJ(5) * t300 + qJD(5) * t337;
t514 = -t207 - t233;
t232 = pkin(3) * t332 + qJ(4) * t331 - qJD(4) * t424;
t219 = t402 * t232;
t328 = pkin(3) * t378 - qJ(4) * t424;
t291 = t394 * t328;
t513 = t219 + t291;
t223 = rSges(7,1) * t326 + rSges(7,2) * t325 - rSges(7,3) * t404;
t512 = pkin(5) * t376 - pkin(9) * t404 + t223;
t345 = t417 * t370;
t511 = t417 * t233 + t345;
t510 = -t232 - t369;
t247 = rSges(5,1) * t336 - rSges(5,2) * t335 - rSges(5,3) * t424;
t509 = -t247 - t328;
t249 = rSges(5,1) * t338 - rSges(5,2) * t337 - rSges(5,3) * t425;
t329 = pkin(3) * t380 - qJ(4) * t425;
t508 = -t249 - t329;
t283 = pkin(4) * t336 + qJ(5) * t335;
t314 = t402 * t328;
t507 = t402 * t283 + t314;
t265 = pkin(4) * t354 + qJ(5) * t353 + qJD(5) * t375;
t312 = pkin(3) * t383 + qJ(4) * t382 + qJD(4) * t404;
t506 = -t265 - t312;
t273 = rSges(5,1) * t354 - rSges(5,2) * t353 + rSges(5,3) * t382;
t505 = -t273 - t312;
t284 = pkin(4) * t338 + qJ(5) * t337;
t320 = t329 * t475;
t504 = t284 * t475 + t320;
t503 = -t283 - t328;
t502 = -t284 - t329;
t311 = rSges(5,1) * t376 - rSges(5,2) * t375 + rSges(5,3) * t404;
t374 = pkin(3) * t405 + qJ(4) * t404;
t501 = -t311 - t374;
t500 = t328 * t523 - t374 * t430;
t373 = pkin(2) * t403 + pkin(8) * t402;
t371 = t417 * t373;
t499 = t417 * t329 + t371;
t327 = pkin(4) * t376 + qJ(5) * t375;
t498 = -t327 - t374;
t497 = 0.2e1 * t496;
t372 = pkin(2) * t401 - pkin(8) * t430;
t495 = t372 * t525 + t373 * t524;
t494 = qJD(2) * t415;
t487 = t417 * t207 + t511;
t486 = -t204 + t514;
t206 = pkin(4) * t299 + qJ(5) * t298 + qJD(5) * t335;
t485 = -t206 + t510;
t484 = t232 * t523 - t312 * t430 + t392 * t374;
t246 = rSges(6,1) * t336 - rSges(6,2) * t424 + rSges(6,3) * t335;
t483 = -t246 + t503;
t248 = rSges(6,1) * t338 - rSges(6,2) * t425 + rSges(6,3) * t337;
t482 = -t248 + t502;
t272 = rSges(6,1) * t354 + rSges(6,2) * t382 + rSges(6,3) * t353;
t481 = -t272 + t506;
t480 = t417 * t284 + t499;
t310 = rSges(6,1) * t376 + rSges(6,2) * t404 + rSges(6,3) * t375;
t479 = -t310 + t498;
t472 = 0.2e1 * m(4);
t471 = 0.2e1 * m(5);
t470 = 0.2e1 * m(6);
t468 = 0.2e1 * m(7);
t319 = rSges(4,1) * t383 - rSges(4,2) * t382 + rSges(4,3) * t475;
t399 = (pkin(2) * t421 + pkin(8) * t419) * t494;
t467 = (-t319 - t399) * t415;
t358 = t405 * rSges(4,1) - t404 * rSges(4,2) - rSges(4,3) * t523;
t406 = (pkin(2) * t419 - pkin(8) * t421) * t415;
t466 = (-t358 - t406) * t415;
t465 = t543 + t489;
t464 = t514 - t519;
t463 = t506 - t518;
t462 = t503 - t517;
t461 = t502 - t516;
t183 = t402 * t206;
t263 = t394 * t283;
t460 = t183 + t263 + t513;
t459 = t233 * t544 + t329 * t545 + 0.2e1 * t219 + 0.2e1 * t291;
t458 = t498 - t512;
t226 = t232 * t525;
t227 = t233 * t524;
t457 = 0.2e1 * t226 + 0.2e1 * t227 + t497;
t456 = t226 + t227 + t496;
t455 = t283 * t523 - t327 * t430 + t500;
t454 = t328 * t525 + t329 * t524 + t495;
t20 = t392 * t461 + t394 * t517 + t402 * t520 - t430 * t464 + t460;
t202 = rSges(6,1) * t299 + rSges(6,2) * t331 + rSges(6,3) * t298;
t49 = t202 * t402 + t246 * t394 + t392 * t482 - t430 * t486 + t460;
t452 = m(6) * t49 + m(7) * t20;
t432 = t206 * t523 - t265 * t430 + t392 * t327 + t484;
t37 = -t518 * t430 + t512 * t392 + (t421 * t520 + t462 * t493) * t415 + t432;
t54 = -t430 * t272 + t392 * t310 + (t202 * t421 + t483 * t493) * t415 + t432;
t451 = m(6) * t54 + m(7) * t37;
t38 = t463 * t402 + t458 * t394 + (t421 * t464 + t493 * t516) * t415 + t504;
t55 = t481 * t402 + t479 * t394 + (t248 * t493 + t421 * t486) * t415 + t504;
t450 = m(6) * t55 + m(7) * t38;
t186 = t206 * t525;
t187 = t207 * t524;
t433 = t186 + t187 + t456;
t45 = (t414 * t520 + t416 * t519) * t415 + t433;
t66 = (t202 * t414 + t204 * t416) * t415 + t433;
t449 = m(6) * t66 + m(7) * t45;
t429 = (-t399 + t463) * t415;
t52 = t414 * t429 + t417 * t519 + t487;
t436 = (-t399 + t481) * t415;
t78 = t204 * t417 + t414 * t436 + t487;
t448 = m(6) * t78 + m(7) * t52;
t53 = (t485 - t520) * t417 + t416 * t429;
t79 = (-t202 + t485) * t417 + t416 * t436;
t447 = m(6) * t79 + m(7) * t53;
t431 = t283 * t525 + t284 * t524 + t454;
t72 = (t414 * t517 + t416 * t516) * t415 + t431;
t99 = (t246 * t414 + t248 * t416) * t415 + t431;
t446 = m(6) * t99 + m(7) * t72;
t102 = t246 * t402 - t430 * t482 + t507;
t73 = t402 * t517 - t430 * t461 + t507;
t445 = m(6) * t102 + m(7) * t73;
t435 = (-t406 + t479) * t415;
t117 = t248 * t417 + t414 * t435 + t480;
t428 = (-t406 + t458) * t415;
t88 = t414 * t428 + t417 * t516 + t480;
t444 = m(6) * t117 + m(7) * t88;
t118 = (-t372 + t483) * t417 + t416 * t435;
t89 = (-t372 + t462) * t417 + t416 * t428;
t443 = m(6) * t118 + m(7) * t89;
t127 = t246 * t523 - t310 * t430 + t455;
t90 = -t430 * t512 + t517 * t523 + t455;
t442 = m(6) * t127 + m(7) * t90;
t128 = t402 * t479 + t482 * t523;
t91 = t402 * t458 + t461 * t523;
t441 = m(6) * t128 + m(7) * t91;
t440 = (-t399 + t505) * t415;
t439 = (-t406 + t501) * t415;
t438 = 0.2e1 * t489;
t51 = t113 * t425 - t114 * t424 - t178 * t333 + t179 * t331;
t203 = rSges(5,1) * t299 - rSges(5,2) * t298 + rSges(5,3) * t331;
t259 = rSges(4,1) * t332 - rSges(4,2) * t331 + rSges(4,3) * t392;
t422 = m(4) * t259 + m(5) * t203 + m(6) * t202 + m(7) * t520;
t398 = (rSges(3,1) * t421 - rSges(3,2) * t419) * t494;
t397 = (Icges(3,1) * t421 - t527) * t494;
t396 = (-Icges(3,2) * t419 + t526) * t494;
t395 = (Icges(3,5) * t421 - Icges(3,6) * t419) * t494;
t388 = t417 * rSges(3,3) + (rSges(3,1) * t419 + rSges(3,2) * t421) * t415;
t387 = Icges(3,5) * t417 + (Icges(3,1) * t419 + t526) * t415;
t386 = Icges(3,6) * t417 + (Icges(3,2) * t421 + t527) * t415;
t368 = -rSges(3,1) * t393 - rSges(3,2) * t394;
t367 = rSges(3,1) * t391 - rSges(3,2) * t392;
t366 = -Icges(3,1) * t393 - Icges(3,4) * t394;
t365 = Icges(3,1) * t391 - Icges(3,4) * t392;
t364 = -Icges(3,4) * t393 - Icges(3,2) * t394;
t363 = Icges(3,4) * t391 - Icges(3,2) * t392;
t362 = -Icges(3,5) * t393 - Icges(3,6) * t394;
t361 = Icges(3,5) * t391 - Icges(3,6) * t392;
t351 = rSges(3,1) * t403 - rSges(3,2) * t402 + rSges(3,3) * t525;
t350 = rSges(3,1) * t401 + rSges(3,2) * t430 - rSges(3,3) * t524;
t349 = Icges(3,1) * t403 - Icges(3,4) * t402 + Icges(3,5) * t525;
t348 = Icges(3,1) * t401 + Icges(3,4) * t430 - Icges(3,5) * t524;
t347 = Icges(3,4) * t403 - Icges(3,2) * t402 + Icges(3,6) * t525;
t346 = Icges(3,4) * t401 + Icges(3,2) * t430 - Icges(3,6) * t524;
t309 = rSges(4,1) * t380 + rSges(4,2) * t425 + rSges(4,3) * t402;
t308 = rSges(4,1) * t378 + rSges(4,2) * t424 - rSges(4,3) * t430;
t231 = -t309 * t523 - t402 * t358;
t230 = t308 * t523 - t358 * t430;
t214 = -t355 * t523 - t404 * t356 + t405 * t357;
t213 = t308 * t402 + t309 * t430;
t212 = (-t308 - t372) * t417 + t416 * t466;
t211 = t309 * t417 + t414 * t466 + t371;
t180 = (t308 * t414 + t309 * t416) * t415 + t495;
t171 = -t293 * t523 - t404 * t295 + t405 * t297;
t170 = -t292 * t523 - t404 * t294 + t405 * t296;
t169 = (-t259 - t369) * t417 + t416 * t467;
t168 = t260 * t417 + t414 * t467 + t345;
t163 = t303 * t404 - t305 * t375 + t307 * t376;
t162 = t302 * t375 + t304 * t404 + t306 * t376;
t152 = t402 * t501 + t508 * t523;
t151 = t247 * t523 - t311 * t430 + t500;
t150 = (t259 * t414 + t260 * t416) * t415 + t496;
t149 = -t402 * t319 - t394 * t358 + (-t260 * t421 + t309 * t493) * t415;
t148 = -t430 * t319 + t392 * t358 + (t259 * t421 - t308 * t493) * t415;
t143 = (-t372 + t509) * t417 + t416 * t439;
t142 = t249 * t417 + t414 * t439 + t499;
t141 = -t404 * t317 + t405 * t318 - t382 * t356 + t383 * t357 + (-t316 * t421 + t355 * t493) * t415;
t140 = -t179 * t404 - t223 * t425;
t139 = t178 * t404 + t223 * t424;
t138 = t247 * t402 - t430 * t508 + t314;
t137 = t259 * t402 + t260 * t430 + t308 * t394 - t309 * t392;
t136 = t237 * t404 - t241 * t375 + t245 * t376;
t135 = t236 * t404 - t240 * t375 + t244 * t376;
t134 = t235 * t375 + t239 * t404 + t243 * t376;
t133 = t234 * t375 + t238 * t404 + t242 * t376;
t132 = (t247 * t414 + t249 * t416) * t415 + t454;
t116 = t178 * t425 - t179 * t424;
t115 = -t220 * t404 + t221 * t325 + t222 * t326;
t106 = (-t203 + t510) * t417 + t416 * t440;
t105 = t205 * t417 + t414 * t440 + t511;
t101 = -t404 * t256 + t405 * t258 - t382 * t295 + t383 * t297 + (-t254 * t421 + t293 * t493) * t415;
t100 = -t404 * t255 + t405 * t257 - t382 * t294 + t383 * t296 + (-t253 * t421 + t292 * t493) * t415;
t98 = t267 * t404 - t269 * t375 + t271 * t376 + t303 * t382 - t305 * t353 + t307 * t354;
t97 = t266 * t375 + t268 * t404 + t270 * t376 + t302 * t353 + t304 * t382 + t306 * t354;
t92 = (t203 * t414 + t205 * t416) * t415 + t456;
t87 = -t173 * t404 + t175 * t325 + t177 * t326;
t86 = -t172 * t404 + t174 * t325 + t176 * t326;
t81 = t320 + t505 * t402 + t501 * t394 + (t249 * t493 + t421 * t515) * t415;
t80 = -t430 * t273 + t392 * t311 + (t203 * t421 + t493 * t509) * t415 + t484;
t71 = t203 * t402 + t247 * t394 + t392 * t508 - t430 * t515 + t513;
t70 = t193 * t404 - t197 * t375 + t201 * t376 + t237 * t382 - t241 * t353 + t245 * t354;
t69 = t192 * t404 - t196 * t375 + t200 * t376 + t236 * t382 - t240 * t353 + t244 * t354;
t68 = t191 * t375 + t195 * t404 + t199 * t376 + t235 * t353 + t239 * t382 + t243 * t354;
t67 = t190 * t375 + t194 * t404 + t198 * t376 + t234 * t353 + t238 * t382 + t242 * t354;
t57 = -t114 * t404 - t156 * t425 - t179 * t382 + t223 * t333;
t56 = t113 * t404 + t156 * t424 + t178 * t382 - t223 * t331;
t50 = t141 * t417 + (-t100 * t416 + t101 * t414) * t415;
t48 = -t153 * t404 + t154 * t325 + t155 * t326 + t215 * t221 + t216 * t222 - t220 * t382;
t47 = t131 * t417 + (t414 * t96 - t416 * t95) * t415;
t46 = t130 * t417 + (t414 * t94 - t416 * t93) * t415;
t42 = t115 * t417 + (t414 * t87 - t416 * t86) * t415;
t41 = -t115 * t523 + t87 * t402 - t430 * t86;
t40 = -t115 * t404 + t424 * t86 + t425 * t87;
t39 = -t100 * t430 + t101 * t402 + t170 * t392 + t171 * t394 + (-t141 * t421 + t214 * t493) * t415;
t36 = t104 * t417 + (t414 * t77 - t416 * t76) * t415;
t35 = t103 * t417 + (t414 * t75 - t416 * t74) * t415;
t34 = -t104 * t523 + t77 * t402 - t430 * t76;
t33 = -t103 * t523 + t75 * t402 - t430 * t74;
t32 = -t104 * t404 + t424 * t76 + t425 * t77;
t31 = -t103 * t404 + t424 * t74 + t425 * t75;
t28 = -t108 * t404 + t110 * t325 + t112 * t326 - t173 * t382 + t175 * t215 + t177 * t216;
t27 = -t107 * t404 + t109 * t325 + t111 * t326 - t172 * t382 + t174 * t215 + t176 * t216;
t26 = t417 * t98 + (t414 * t70 - t416 * t69) * t415;
t25 = t417 * t97 + (t414 * t68 - t416 * t67) * t415;
t19 = t417 * t85 + (t414 * t65 - t416 * t64) * t415;
t18 = t417 * t84 + (t414 * t63 - t416 * t62) * t415;
t17 = t417 * t83 + (t414 * t61 - t416 * t60) * t415;
t16 = t417 * t82 + (t414 * t59 - t416 * t58) * t415;
t15 = t135 * t392 + t136 * t394 - t69 * t430 + t70 * t402 + (t163 * t493 - t421 * t98) * t415;
t14 = t133 * t392 + t134 * t394 - t67 * t430 + t68 * t402 + (t162 * t493 - t421 * t97) * t415;
t9 = t417 * t48 + (-t27 * t416 + t28 * t414) * t415;
t8 = t417 * t44 + (-t23 * t416 + t24 * t414) * t415;
t7 = t417 * t43 + (-t21 * t416 + t22 * t414) * t415;
t6 = -t27 * t430 + t28 * t402 + t86 * t392 + t87 * t394 + (t115 * t493 - t421 * t48) * t415;
t5 = -t115 * t382 + t27 * t424 + t28 * t425 - t331 * t86 - t333 * t87 - t404 * t48;
t2 = -t104 * t382 + t23 * t424 + t24 * t425 - t331 * t76 - t333 * t77 - t404 * t44;
t1 = -t103 * t382 + t21 * t424 + t22 * t425 - t331 * t74 - t333 * t75 - t404 * t43;
t10 = [0; m(4) * t497 / 0.2e1 + t457 * t543 + (m(3) * t368 - t546) * t524 + (m(3) * t367 + t422) * t525 + t489 * (0.2e1 * t186 + 0.2e1 * t187 + t457); t8 * t525 - t7 * t524 - t46 * t524 - t17 * t524 + t47 * t525 + t19 * t525 + t18 * t525 - t16 * t524 - ((-t347 * t392 + t349 * t391 - t362 * t524 + t364 * t430 + t366 * t401) * t525 - (-t346 * t392 + t348 * t391 - t361 * t524 + t363 * t430 + t365 * t401) * t524 + (-t386 * t392 + t387 * t391 - t395 * t524 + t396 * t430 + t397 * t401) * t417) * t524 + ((-t347 * t394 - t349 * t393 + t362 * t525 - t364 * t402 + t366 * t403) * t525 - (-t346 * t394 - t348 * t393 + t361 * t525 - t363 * t402 + t365 * t403) * t524 + (-t386 * t394 - t387 * t393 + t395 * t525 - t396 * t402 + t397 * t403) * t417) * t525 + t417 * t9 + (t45 * t72 + t52 * t88 + t53 * t89) * t468 + (t117 * t78 + t118 * t79 + t66 * t99) * t470 + (t105 * t142 + t106 * t143 + t132 * t92) * t471 + t417 * t50 + t417 * t26 + t417 * t25 + (t150 * t180 + t168 * t211 + t169 * t212) * t472 + t417 * (t417 ^ 2 * t395 + (((t364 * t421 + t366 * t419) * t414 - (t363 * t421 + t365 * t419) * t416 + ((-t347 * t419 + t349 * t421) * t414 - (-t346 * t419 + t348 * t421) * t416) * qJD(2)) * t415 + (-t361 * t416 + t362 * t414 + t396 * t421 + t397 * t419 + (-t386 * t419 + t387 * t421) * qJD(2)) * t417) * t415) + 0.2e1 * m(3) * ((-t350 * t417 - t388 * t524) * (-t367 * t417 - t398 * t524) + (t351 * t417 - t388 * t525) * (t368 * t417 - t398 * t525) + (t350 * t414 + t351 * t416) * t415 ^ 2 * (t367 * t414 + t368 * t416)); t459 * t543 + t422 * t402 - t546 * t430 + (m(4) * t308 + m(5) * t247 + m(6) * t246 + m(7) * t517) * t394 + (-m(4) * t309 - m(5) * t249 - m(6) * t248 - m(7) * t516) * t392 + t489 * (t207 * t544 + t284 * t545 + 0.2e1 * t183 + 0.2e1 * t263 + t459); (t20 * t72 + t37 * t89 + t38 * t88 + t45 * t73 + t52 * t91 + t53 * t90) * m(7) + (t19 / 0.2e1 + t18 / 0.2e1 + t8 / 0.2e1 + t47 / 0.2e1) * t402 - (t17 / 0.2e1 + t16 / 0.2e1 + t7 / 0.2e1 + t46 / 0.2e1) * t430 + (t15 / 0.2e1 + t14 / 0.2e1 + t6 / 0.2e1 + t39 / 0.2e1 + (t210 / 0.2e1 + t147 / 0.2e1 + t146 / 0.2e1) * t394 + (t209 / 0.2e1 + t144 / 0.2e1 + t145 / 0.2e1) * t392) * t417 + (-t524 * t552 + t525 * t551 + t35) * t535 + (-t524 * t550 + t525 * t549 + t36) * t534 + m(4) * (t137 * t180 + t148 * t212 + t149 * t211 + t150 * t213 + t168 * t231 + t169 * t230) + m(5) * (t105 * t152 + t106 * t151 + t132 * t71 + t138 * t92 + t142 * t81 + t143 * t80) + m(6) * (t102 * t66 + t117 * t55 + t118 * t54 + t127 * t79 + t128 * t78 + t49 * t99) + ((-t50 / 0.2e1 - t25 / 0.2e1 - t26 / 0.2e1 - t9 / 0.2e1) * t421 + (t42 / 0.2e1 + (-t170 / 0.2e1 - t133 / 0.2e1 - t135 / 0.2e1) * t524 + (t171 / 0.2e1 + t134 / 0.2e1 + t136 / 0.2e1) * t525 + (t214 / 0.2e1 + t162 / 0.2e1 + t163 / 0.2e1) * t417) * t493 + t547 * t532 + t548 * t531) * t415; (t137 * t213 + t148 * t230 + t149 * t231) * t472 + (t138 * t71 + t151 * t80 + t152 * t81) * t471 + (t102 * t49 + t127 * t54 + t128 * t55) * t470 + (t20 * t73 + t37 * t90 + t38 * t91) * t468 + (-t14 - t15 - t39 - t6) * t523 + t547 * t402 - t548 * t430 + (t549 * t402 - t550 * t430 - t523 * t553 + t34) * t394 + (t551 * t402 - t552 * t430 - t523 * t554 + t33) * t392 + (t41 + (-t162 - t163 - t214) * t523 + (t134 + t136 + t171) * t402 - (t133 + t135 + t170) * t430) * t475; 0.2e1 * t465 * t382; (m(5) * t92 + t449) * t404 + (m(5) * t132 + t446) * t382 - (m(5) * t106 + t447) * t425 - (m(5) * t105 + t448) * t424 + (m(5) * t143 + t443) * t333 + (m(5) * t142 + t444) * t331; (m(5) * t71 + t452) * t404 + (m(5) * t138 + t445) * t382 - (m(5) * t80 + t451) * t425 - (m(5) * t81 + t450) * t424 + (m(5) * t151 + t442) * t333 + (m(5) * t152 + t441) * t331; 0.4e1 * t465 * (-t331 * t424 - t333 * t425 + t382 * t404); t353 * t438; t298 * t444 + t300 * t443 + t335 * t448 + t337 * t447 + t353 * t446 + t375 * t449; t298 * t441 + t300 * t442 + t335 * t450 + t337 * t451 + t353 * t445 + t375 * t452; (-t298 * t424 - t300 * t425 + t331 * t335 + t333 * t337 + t353 * t404 + t375 * t382) * t438; 0.4e1 * t489 * (t298 * t335 + t300 * t337 + t353 * t375); t51 * m(7); t42 * t536 + t9 * t533 + t36 * t539 + t8 * t537 + t35 * t540 + t7 * t538 + t417 * t5 / 0.2e1 + (t116 * t45 + t139 * t53 + t140 * t52 + t51 * t72 + t56 * t89 + t57 * t88) * m(7) + (t1 * t531 + t2 * t532) * t415; (t116 * t20 + t139 * t37 + t140 * t38 + t51 * t73 + t56 * t90 + t57 * t91) * m(7) + t33 * t540 + t3 * t538 + t34 * t539 + t4 * t537 + t41 * t536 + t6 * t533 + t32 * t534 + t402 * t2 / 0.2e1 + t31 * t535 - t430 * t1 / 0.2e1 + (t40 * t493 / 0.2e1 - t421 * t5 / 0.2e1) * t415; (t116 * t382 + t139 * t333 + t140 * t331 + t404 * t51 - t424 * t57 - t425 * t56) * m(7); (t116 * t353 + t139 * t300 + t140 * t298 + t335 * t57 + t337 * t56 + t375 * t51) * m(7); -t333 * t32 + t425 * t2 - t331 * t31 + t424 * t1 - t382 * t40 - t404 * t5 + (t116 * t51 + t139 * t56 + t140 * t57) * t468;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;

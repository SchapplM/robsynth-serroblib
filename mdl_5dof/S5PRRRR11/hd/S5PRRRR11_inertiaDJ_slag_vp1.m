% Calculate time derivative of joint inertia matrix for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR11_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR11_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR11_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:14
% EndTime: 2024-09-27 21:45:22
% DurationCPUTime: 4.86s
% Computational Cost: add. (54549->692), mult. (33047->969), div. (0->0), fcn. (26701->20), ass. (0->362)
t343 = pkin(10) + qJ(2);
t333 = sin(t343);
t334 = cos(t343);
t350 = sin(qJ(3));
t348 = cos(pkin(5));
t352 = cos(qJ(3));
t418 = t348 * t352;
t365 = t333 * t418 + t334 * t350;
t359 = t365 * qJD(3);
t450 = pkin(3) * t359;
t332 = pkin(3) * t352 + pkin(2);
t439 = pkin(2) - t332;
t449 = t439 * t333;
t347 = sin(pkin(5));
t427 = t334 * t347;
t314 = pkin(7) * t427;
t448 = -pkin(2) * t333 + t314;
t346 = qJ(3) + qJ(4);
t336 = pkin(5) + t346;
t323 = sin(t336);
t317 = t323 / 0.2e1;
t337 = pkin(5) - t346;
t324 = sin(t337);
t296 = t317 + t324 / 0.2e1;
t326 = cos(t337);
t318 = t326 / 0.2e1;
t325 = cos(t336);
t299 = t318 - t325 / 0.2e1;
t243 = Icges(5,5) * t299 + Icges(5,6) * t296 + Icges(5,3) * t348;
t244 = Icges(5,4) * t299 + Icges(5,2) * t296 + Icges(5,6) * t348;
t245 = Icges(5,1) * t299 + Icges(5,4) * t296 + Icges(5,5) * t348;
t298 = t318 + t325 / 0.2e1;
t338 = sin(t346);
t255 = t298 * t334 - t333 * t338;
t297 = t317 - t324 / 0.2e1;
t339 = cos(t346);
t256 = t297 * t334 + t333 * t339;
t104 = -t243 * t427 + t244 * t255 + t245 * t256;
t257 = -t298 * t333 - t334 * t338;
t258 = -t297 * t333 + t334 * t339;
t429 = t333 * t347;
t105 = t243 * t429 + t244 * t257 + t245 * t258;
t344 = qJD(3) + qJD(4);
t442 = t344 / 0.2e1;
t305 = t324 * t442;
t422 = t344 * t323;
t271 = t305 - t422 / 0.2e1;
t423 = t339 * t344;
t189 = -qJD(2) * t255 - t271 * t333 - t334 * t423;
t304 = t325 * t442;
t421 = t344 * t326;
t269 = t304 + t421 / 0.2e1;
t424 = t338 * t344;
t190 = -qJD(2) * t256 - t269 * t333 - t334 * t424;
t403 = qJD(2) * t347;
t389 = t334 * t403;
t113 = Icges(5,5) * t190 + Icges(5,6) * t189 + Icges(5,3) * t389;
t191 = qJD(2) * t257 + t271 * t334 - t333 * t423;
t192 = qJD(2) * t258 + t269 * t334 - t333 * t424;
t390 = t333 * t403;
t114 = Icges(5,5) * t192 + Icges(5,6) * t191 + Icges(5,3) * t390;
t115 = Icges(5,4) * t190 + Icges(5,2) * t189 + Icges(5,6) * t389;
t116 = Icges(5,4) * t192 + Icges(5,2) * t191 + Icges(5,6) * t390;
t117 = Icges(5,1) * t190 + Icges(5,4) * t189 + Icges(5,5) * t389;
t118 = Icges(5,1) * t192 + Icges(5,4) * t191 + Icges(5,5) * t390;
t183 = Icges(5,5) * t256 + Icges(5,6) * t255 - Icges(5,3) * t427;
t184 = Icges(5,5) * t258 + Icges(5,6) * t257 + Icges(5,3) * t429;
t185 = Icges(5,4) * t256 + Icges(5,2) * t255 - Icges(5,6) * t427;
t186 = Icges(5,4) * t258 + Icges(5,2) * t257 + Icges(5,6) * t429;
t187 = Icges(5,1) * t256 + Icges(5,4) * t255 - Icges(5,5) * t427;
t188 = Icges(5,1) * t258 + Icges(5,4) * t257 + Icges(5,5) * t429;
t268 = t304 - t421 / 0.2e1;
t270 = t305 + t422 / 0.2e1;
t22 = t114 * t348 + t116 * t296 + t118 * t299 + t185 * t268 + t187 * t270;
t23 = t113 * t348 + t115 * t296 + t117 * t299 + t186 * t268 + t188 * t270;
t224 = Icges(5,5) * t270 + Icges(5,6) * t268;
t225 = Icges(5,4) * t270 + Icges(5,2) * t268;
t226 = Icges(5,1) * t270 + Icges(5,4) * t268;
t404 = qJD(2) * t334;
t34 = t189 * t244 + t190 * t245 + t225 * t257 + t226 * t258 + (t224 * t333 + t243 * t404) * t347;
t405 = qJD(2) * t333;
t367 = t224 * t348 + t225 * t296 + t226 * t299 + t244 * t268 + t245 * t270;
t54 = t367 * t348;
t57 = -t183 * t427 + t185 * t255 + t187 * t256;
t58 = -t184 * t427 + t186 * t255 + t188 * t256;
t59 = t183 * t429 + t185 * t257 + t187 * t258;
t60 = t184 * t429 + t186 * t257 + t188 * t258;
t72 = t183 * t348 + t185 * t296 + t187 * t299;
t73 = t184 * t348 + t186 * t296 + t188 * t299;
t447 = (t104 * t348 + (t333 * t58 - t334 * t57) * t347) * t390 + (t105 * t348 + (t333 * t60 - t334 * t59) * t347) * t389 + (t34 * t348 + ((t257 * t115 + t258 * t117 + t189 * t186 + t190 * t188) * t333 + t60 * t404 - (t257 * t116 + t258 * t118 + t189 * t185 + t190 * t187) * t334 + t59 * t405 + ((t113 * t333 + t184 * t404) * t333 - (t114 * t333 + t183 * t404) * t334) * t347) * t347) * t429 + t348 * (t54 + (-t22 * t334 + t23 * t333 + (t333 * t72 + t334 * t73) * qJD(2)) * t347);
t303 = pkin(4) * t339 + t332;
t406 = t303 - t332;
t351 = cos(qJ(4));
t331 = pkin(4) * t351 + pkin(3);
t353 = pkin(7) + pkin(8);
t345 = pkin(9) + t353;
t349 = sin(qJ(4));
t399 = pkin(4) * t349 * t352;
t259 = -t345 * t347 + (t331 * t350 + t399) * t348;
t419 = t348 * t350;
t302 = pkin(3) * t419 - t347 * t353;
t408 = t259 - t302;
t163 = t333 * t406 + t334 * t408;
t446 = 2 * m(4);
t445 = 2 * m(5);
t444 = 2 * m(6);
t335 = qJD(5) + t344;
t443 = t335 / 0.2e1;
t440 = pkin(2) * t334;
t438 = -pkin(3) + t331;
t401 = qJD(3) * t350;
t355 = (-t349 * t401 + (-t349 * t350 + t351 * t352) * qJD(4)) * pkin(4);
t400 = qJD(3) * t352;
t228 = (t331 * t400 + t355) * t348;
t396 = pkin(3) * t401;
t295 = -pkin(4) * t424 - t396;
t379 = -t228 * t333 + t295 * t334;
t102 = -qJD(2) * t163 + t379 + t450;
t328 = -qJ(5) + t337;
t322 = cos(t328);
t316 = t322 / 0.2e1;
t327 = qJ(5) + t336;
t321 = cos(t327);
t293 = t316 + t321 / 0.2e1;
t340 = qJ(5) + t346;
t329 = sin(t340);
t249 = t293 * t334 - t329 * t333;
t320 = sin(t328);
t301 = t320 * t443;
t319 = sin(t327);
t426 = t335 * t319;
t265 = t301 - t426 / 0.2e1;
t330 = cos(t340);
t428 = t334 * t335;
t153 = -qJD(2) * t249 - t265 * t333 - t330 * t428;
t315 = t319 / 0.2e1;
t292 = t315 - t320 / 0.2e1;
t250 = t292 * t334 + t330 * t333;
t300 = t321 * t443;
t425 = t335 * t322;
t263 = t300 + t425 / 0.2e1;
t154 = -qJD(2) * t250 - t263 * t333 - t329 * t428;
t100 = rSges(6,1) * t154 + rSges(6,2) * t153 + rSges(6,3) * t389;
t92 = t348 * t100;
t437 = t102 * t348 + t92;
t436 = Icges(4,4) * t350;
t435 = Icges(4,4) * t352;
t227 = rSges(5,1) * t270 + rSges(5,2) * t268;
t434 = t227 * t333;
t433 = t302 * t334;
t432 = t303 * t333;
t431 = t332 * t333;
t430 = t333 * t335;
t290 = t334 * t303;
t306 = t334 * t332;
t420 = t347 * t350;
t251 = -t293 * t333 - t329 * t334;
t155 = qJD(2) * t251 + t265 * t334 - t330 * t430;
t252 = -t292 * t333 + t330 * t334;
t156 = qJD(2) * t252 + t263 * t334 - t329 * t430;
t370 = -rSges(6,1) * t156 - rSges(6,2) * t155;
t101 = rSges(6,3) * t390 - t370;
t395 = pkin(3) * t400;
t377 = t348 * t395;
t407 = t302 * t405 + t333 * t396;
t103 = (-qJD(2) * t259 + t295) * t333 + (qJD(2) * t406 + t228 - t377) * t334 + t407;
t417 = -t101 - t103;
t369 = -rSges(6,1) * t250 - rSges(6,2) * t249;
t173 = -rSges(6,3) * t427 - t369;
t174 = rSges(6,1) * t252 + rSges(6,2) * t251 + rSges(6,3) * t429;
t106 = t173 * t429 + t174 * t427;
t159 = t348 * t174;
t164 = -t333 * t408 + t290 - t306;
t416 = t164 * t348 + t159;
t415 = -t163 - t173;
t414 = -t164 - t174;
t371 = -rSges(5,1) * t256 - rSges(5,2) * t255;
t193 = -rSges(5,3) * t427 - t371;
t312 = rSges(5,3) * t429;
t194 = rSges(5,1) * t258 + rSges(5,2) * t257 + t312;
t112 = t193 * t429 + t194 * t427;
t262 = t300 - t425 / 0.2e1;
t264 = t301 + t426 / 0.2e1;
t207 = rSges(6,1) * t264 + rSges(6,2) * t262;
t413 = -t207 - (t400 * t438 + t355) * t347;
t229 = t314 + t433 - t449;
t385 = -pkin(7) * t347 - t302;
t230 = t333 * t385 + t306 - t440;
t412 = t229 * t429 + t230 * t427;
t291 = t315 + t320 / 0.2e1;
t294 = t316 - t321 / 0.2e1;
t242 = rSges(6,1) * t294 + rSges(6,2) * t291 + rSges(6,3) * t348;
t222 = t242 * t390;
t241 = (t345 - t353) * t348 + (t350 * t438 + t399) * t347;
t411 = t241 * t390 + t222;
t410 = -t241 - t242;
t246 = rSges(5,1) * t299 + rSges(5,2) * t296 + rSges(5,3) * t348;
t285 = pkin(3) * t420 + (-pkin(7) + t353) * t348;
t409 = -t246 - t285;
t402 = qJD(3) * t347;
t398 = pkin(7) * t429;
t397 = t100 * t427 + t101 * t429 + t173 * t389;
t121 = rSges(5,1) * t190 + rSges(5,2) * t189 + rSges(5,3) * t389;
t372 = -rSges(5,1) * t192 - rSges(5,2) * t191;
t122 = rSges(5,3) * t390 - t372;
t394 = t121 * t427 + t122 * t429 + t193 * t389;
t181 = -t450 + (t334 * t385 + t449) * qJD(2);
t361 = t334 * t377 - t407;
t182 = (-t334 * t439 - t398) * qJD(2) + t361;
t393 = t181 * t427 + t182 * t429 + t229 * t389;
t281 = -t333 * t350 + t334 * t418;
t364 = t333 * t419 - t334 * t352;
t234 = -qJD(2) * t281 + qJD(3) * t364;
t282 = t333 * t352 + t334 * t419;
t235 = -qJD(2) * t282 - t359;
t142 = rSges(4,1) * t235 + rSges(4,2) * t234 + rSges(4,3) * t389;
t392 = -t285 + t410;
t219 = -rSges(4,1) * t364 - rSges(4,2) * t365 + rSges(4,3) * t429;
t169 = Icges(6,4) * t250 + Icges(6,2) * t249 - Icges(6,6) * t427;
t171 = Icges(6,1) * t250 + Icges(6,4) * t249 - Icges(6,5) * t427;
t95 = Icges(6,5) * t156 + Icges(6,6) * t155 + Icges(6,3) * t390;
t97 = Icges(6,4) * t156 + Icges(6,2) * t155 + Icges(6,6) * t390;
t99 = Icges(6,1) * t156 + Icges(6,4) * t155 + Icges(6,5) * t390;
t13 = t169 * t262 + t171 * t264 + t291 * t97 + t294 * t99 + t348 * t95;
t170 = Icges(6,4) * t252 + Icges(6,2) * t251 + Icges(6,6) * t429;
t172 = Icges(6,1) * t252 + Icges(6,4) * t251 + Icges(6,5) * t429;
t94 = Icges(6,5) * t154 + Icges(6,6) * t153 + Icges(6,3) * t389;
t96 = Icges(6,4) * t154 + Icges(6,2) * t153 + Icges(6,6) * t389;
t98 = Icges(6,1) * t154 + Icges(6,4) * t153 + Icges(6,5) * t389;
t14 = t170 * t262 + t172 * t264 + t291 * t96 + t294 * t98 + t348 * t94;
t167 = Icges(6,5) * t250 + Icges(6,6) * t249 - Icges(6,3) * t427;
t168 = Icges(6,5) * t252 + Icges(6,6) * t251 + Icges(6,3) * t429;
t204 = Icges(6,5) * t264 + Icges(6,6) * t262;
t205 = Icges(6,4) * t264 + Icges(6,2) * t262;
t206 = Icges(6,1) * t264 + Icges(6,4) * t262;
t238 = Icges(6,5) * t294 + Icges(6,6) * t291 + Icges(6,3) * t348;
t239 = Icges(6,4) * t294 + Icges(6,2) * t291 + Icges(6,6) * t348;
t240 = Icges(6,1) * t294 + Icges(6,4) * t291 + Icges(6,5) * t348;
t29 = t153 * t239 + t154 * t240 + t205 * t251 + t206 * t252 + (t204 * t333 + t238 * t404) * t347;
t368 = t204 * t348 + t205 * t291 + t206 * t294 + t239 * t262 + t240 * t264;
t42 = t368 * t348;
t46 = -t167 * t427 + t169 * t249 + t171 * t250;
t47 = -t168 * t427 + t170 * t249 + t172 * t250;
t48 = t167 * t429 + t169 * t251 + t171 * t252;
t49 = t168 * t429 + t170 * t251 + t172 * t252;
t61 = t167 * t348 + t169 * t291 + t171 * t294;
t62 = t168 * t348 + t170 * t291 + t172 * t294;
t82 = -t238 * t427 + t239 * t249 + t240 * t250;
t83 = t238 * t429 + t239 * t251 + t240 * t252;
t391 = (t29 * t348 + ((t153 * t170 + t154 * t172 + t251 * t96 + t252 * t98) * t333 + t49 * t404 - (t153 * t169 + t154 * t171 + t251 * t97 + t252 * t99) * t334 + t48 * t405 + ((t168 * t404 + t333 * t94) * t333 - (t167 * t404 + t333 * t95) * t334) * t347) * t347) * t429 + (t348 * t83 + (t333 * t49 - t334 * t48) * t347) * t389 + t348 * (t42 + (-t13 * t334 + t14 * t333 + (t333 * t61 + t334 * t62) * qJD(2)) * t347) + (t348 * t82 + (t333 * t47 - t334 * t46) * t347) * t390;
t388 = t429 / 0.2e1;
t387 = -t427 / 0.2e1;
t386 = t403 / 0.2e1;
t384 = rSges(6,3) * t347 - t259;
t383 = t413 * t333;
t382 = t413 * t347;
t381 = t410 * t347;
t380 = t409 * t347;
t378 = t347 ^ 2 * t395;
t41 = t163 * t429 + t164 * t427 + t106;
t376 = t333 * t386;
t375 = t334 * t386;
t374 = t392 * t347;
t236 = -qJD(2) * t365 - qJD(3) * t282;
t237 = -qJD(2) * t364 + qJD(3) * t281;
t373 = -rSges(4,1) * t237 - rSges(4,2) * t236;
t366 = t102 * t427 + t103 * t429 + t163 * t389 + t397;
t30 = t155 * t239 + t156 * t240 + t205 * t249 + t206 * t250 + (-t204 * t334 + t238 * t405) * t347;
t2 = t30 * t348 + ((t155 * t170 + t156 * t172 + t249 * t96 + t250 * t98) * t333 + t47 * t404 - (t155 * t169 + t156 * t171 + t249 * t97 + t250 * t99) * t334 + t46 * t405 + ((t168 * t405 - t334 * t94) * t333 - (t167 * t405 - t334 * t95) * t334) * t347) * t347;
t363 = -t2 * t427 + t391;
t362 = t181 * t348 - t333 * t378;
t360 = t42 + (t14 + t29) * t388 + (t13 + t30) * t387 + (t61 + t82) * t376 + (t62 + t83) * t375;
t218 = rSges(4,1) * t282 + rSges(4,2) * t281 - rSges(4,3) * t427;
t278 = Icges(4,6) * t348 + (Icges(4,2) * t352 + t436) * t347;
t279 = Icges(4,5) * t348 + (Icges(4,1) * t350 + t435) * t347;
t286 = (Icges(4,5) * t352 - Icges(4,6) * t350) * t402;
t287 = (-Icges(4,2) * t350 + t435) * t402;
t288 = (Icges(4,1) * t352 - t436) * t402;
t357 = t348 * t286 + t288 * t420 + (-t278 * t401 + t279 * t400 + t287 * t352) * t347;
t35 = t191 * t244 + t192 * t245 + t225 * t255 + t226 * t256 + (-t224 * t334 + t243 * t405) * t347;
t4 = t35 * t348 + ((t255 * t115 + t256 * t117 + t191 * t186 + t192 * t188) * t333 + t58 * t404 - (t255 * t116 + t256 * t118 + t191 * t185 + t192 * t187) * t334 + t57 * t405 + ((-t113 * t334 + t184 * t405) * t333 - (-t114 * t334 + t183 * t405) * t334) * t347) * t347;
t356 = (-t2 - t4) * t427 + t391 + t447;
t354 = t360 + t54 + (t23 + t34) * t388 + (t22 + t35) * t387 + (t104 + t72) * t376 + (t105 + t73) * t375;
t289 = (rSges(4,1) * t352 - rSges(4,2) * t350) * t402;
t280 = rSges(4,3) * t348 + (rSges(4,1) * t350 + rSges(4,2) * t352) * t347;
t277 = Icges(4,3) * t348 + (Icges(4,5) * t350 + Icges(4,6) * t352) * t347;
t260 = t285 * t390;
t231 = t246 * t390;
t223 = t348 * t230;
t217 = -Icges(4,1) * t364 - Icges(4,4) * t365 + Icges(4,5) * t429;
t216 = Icges(4,1) * t282 + Icges(4,4) * t281 - Icges(4,5) * t427;
t215 = -Icges(4,4) * t364 - Icges(4,2) * t365 + Icges(4,6) * t429;
t214 = Icges(4,4) * t282 + Icges(4,2) * t281 - Icges(4,6) * t427;
t213 = -Icges(4,5) * t364 - Icges(4,6) * t365 + Icges(4,3) * t429;
t212 = Icges(4,5) * t282 + Icges(4,6) * t281 - Icges(4,3) * t427;
t198 = t219 + t398 + t440;
t197 = -t218 + t448;
t176 = t348 * t194;
t158 = -t218 * t348 - t280 * t427;
t157 = t219 * t348 - t280 * t429;
t143 = rSges(4,3) * t390 - t373;
t141 = Icges(4,1) * t237 + Icges(4,4) * t236 + Icges(4,5) * t390;
t140 = Icges(4,1) * t235 + Icges(4,4) * t234 + Icges(4,5) * t389;
t139 = Icges(4,4) * t237 + Icges(4,2) * t236 + Icges(4,6) * t390;
t138 = Icges(4,4) * t235 + Icges(4,2) * t234 + Icges(4,6) * t389;
t137 = Icges(4,5) * t237 + Icges(4,6) * t236 + Icges(4,3) * t390;
t136 = Icges(4,5) * t235 + Icges(4,6) * t234 + Icges(4,3) * t389;
t135 = -t302 * t333 + t194 + t306;
t134 = -t431 + (rSges(5,3) * t347 - t302) * t334 + t371;
t133 = t277 * t429 - t278 * t365 - t279 * t364;
t132 = -t277 * t427 + t278 * t281 + t279 * t282;
t131 = (-t440 + (-rSges(4,3) - pkin(7)) * t429) * qJD(2) + t373;
t130 = qJD(2) * t448 + t142;
t129 = -t193 * t348 - t246 * t427;
t128 = -t246 * t429 + t176;
t127 = t357 * t348;
t126 = -t259 * t333 + t174 + t290;
t125 = t334 * t384 + t369 - t432;
t124 = -t173 * t348 - t242 * t427;
t123 = -t242 * t429 + t159;
t120 = t142 * t348 + (-t280 * t404 - t289 * t333) * t347;
t119 = -t143 * t348 + (t280 * t405 - t289 * t334) * t347;
t111 = t348 * t121;
t108 = t213 * t348 + (t215 * t352 + t217 * t350) * t347;
t107 = t212 * t348 + (t214 * t352 + t216 * t350) * t347;
t85 = (-t193 - t229) * t348 + t334 * t380;
t84 = t333 * t380 + t176 + t223;
t81 = t213 * t429 - t215 * t365 - t217 * t364;
t80 = t212 * t429 - t214 * t365 - t216 * t364;
t79 = -t213 * t427 + t215 * t281 + t217 * t282;
t78 = -t212 * t427 + t214 * t281 + t216 * t282;
t75 = (-t312 - t306) * qJD(2) - t361 + t372;
t74 = (-t431 - t433) * qJD(2) - t450 + t121;
t71 = t236 * t278 + t237 * t279 + t281 * t287 + t282 * t288 + (t277 * t405 - t286 * t334) * t347;
t70 = t234 * t278 + t235 * t279 - t365 * t287 - t364 * t288 + (t277 * t404 + t286 * t333) * t347;
t67 = t111 + (-t246 * t404 - t434) * t347;
t66 = -t122 * t348 - t227 * t427 + t231;
t65 = t334 * t381 + t348 * t415;
t64 = t333 * t381 + t416;
t63 = t412 + t112;
t53 = t92 + (-t207 * t333 - t242 * t404) * t347;
t52 = -t101 * t348 - t207 * t427 + t222;
t51 = -t228 * t334 - t295 * t333 + (-t333 * t384 - t290) * qJD(2) + t370;
t50 = (-t259 * t334 - t432) * qJD(2) + t379 + t100;
t45 = (-t229 + t415) * t348 + t334 * t374;
t44 = t333 * t374 + t223 + t416;
t43 = (t142 * t334 + t143 * t333 + (t218 * t334 - t219 * t333) * qJD(2)) * t347;
t40 = t111 + (t404 * t409 - t434) * t347 + t362;
t39 = t231 + t260 + (-t122 - t182) * t348 + (-t227 * t347 - t378) * t334;
t38 = t136 * t348 + (t138 * t352 + t140 * t350 + (-t215 * t350 + t217 * t352) * qJD(3)) * t347;
t37 = t137 * t348 + (t139 * t352 + t141 * t350 + (-t214 * t350 + t216 * t352) * qJD(3)) * t347;
t36 = t41 + t412;
t31 = -t194 * t390 + t394;
t28 = -t174 * t390 + t397;
t25 = (t404 * t410 + t383) * t347 + t437;
t24 = t334 * t382 + t348 * t417 + t411;
t19 = (t392 * t404 + t383) * t347 + t362 + t437;
t18 = t260 + (-t182 + t417) * t348 + (-t378 + t382) * t334 + t411;
t17 = (-t194 - t230) * t390 + t393 + t394;
t8 = t390 * t414 + t366;
t7 = (-t230 + t414) * t390 + t366 + t393;
t1 = [0; 0; (t125 * t51 + t126 * t50) * t444 + (t134 * t75 + t135 * t74) * t445 + (t130 * t198 + t131 * t197) * t446 + t357 + t367 + t368; m(4) * t43 + m(5) * t17 + m(6) * t7; t354 + m(4) * (t119 * t197 + t120 * t198 + t130 * t157 + t131 * t158) + m(6) * (t125 * t18 + t126 * t19 + t44 * t50 + t45 * t51) + m(5) * (t134 * t39 + t135 * t40 + t74 * t84 + t75 * t85) + t127 + ((-t37 / 0.2e1 - t71 / 0.2e1) * t334 + (t38 / 0.2e1 + t70 / 0.2e1) * t333 + ((t108 / 0.2e1 + t133 / 0.2e1) * t334 + (t107 / 0.2e1 + t132 / 0.2e1) * t333) * qJD(2)) * t347; (t133 * t348 + (t333 * t81 - t334 * t80) * t347) * t389 + (t70 * t348 + ((-t138 * t365 - t140 * t364 + t234 * t215 + t235 * t217) * t333 + t81 * t404 - (-t139 * t365 - t141 * t364 + t234 * t214 + t235 * t216) * t334 + t80 * t405 + ((t136 * t333 + t213 * t404) * t333 - (t137 * t333 + t212 * t404) * t334) * t347) * t347) * t429 + (t132 * t348 + (t333 * t79 - t334 * t78) * t347) * t390 + (t18 * t45 + t19 * t44 + t36 * t7) * t444 + (t17 * t63 + t39 * t85 + t40 * t84) * t445 + t348 * (t127 + (t38 * t333 - t37 * t334 + (t107 * t333 + t108 * t334) * qJD(2)) * t347) + (t158 * t119 + t157 * t120 + (t218 * t333 + t219 * t334) * t43 * t347) * t446 + t363 + (-t4 - t71 * t348 - ((t281 * t138 + t282 * t140 + t236 * t215 + t237 * t217) * t333 + t79 * t404 - (t281 * t139 + t282 * t141 + t236 * t214 + t237 * t216) * t334 + t78 * t405 + ((-t136 * t334 + t213 * t405) * t333 - (-t137 * t334 + t212 * t405) * t334) * t347) * t347) * t427 + t447; m(5) * t31 + m(6) * t8; t354 + m(6) * (t125 * t24 + t126 * t25 + t50 * t64 + t51 * t65) + m(5) * (t128 * t74 + t129 * t75 + t134 * t66 + t135 * t67); m(6) * (t18 * t65 + t19 * t64 + t24 * t45 + t25 * t44 + t36 * t8 + t41 * t7) + m(5) * (t112 * t17 + t128 * t40 + t129 * t39 + t31 * t63 + t66 * t85 + t67 * t84) + t356; (t112 * t31 + t128 * t67 + t129 * t66) * t445 + (t24 * t65 + t25 * t64 + t41 * t8) * t444 + t356; m(6) * t28; m(6) * (t123 * t50 + t124 * t51 + t125 * t52 + t126 * t53) + t360; m(6) * (t106 * t7 + t123 * t19 + t124 * t18 + t28 * t36 + t44 * t53 + t45 * t52) + t363; m(6) * (t106 * t8 + t123 * t25 + t124 * t24 + t28 * t41 + t52 * t65 + t53 * t64) + t363; (t106 * t28 + t123 * t53 + t124 * t52) * t444 + t363;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

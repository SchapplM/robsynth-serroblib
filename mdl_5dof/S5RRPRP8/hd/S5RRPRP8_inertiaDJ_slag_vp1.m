% Calculate time derivative of joint inertia matrix for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP8_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP8_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:22
% EndTime: 2019-12-31 20:03:42
% DurationCPUTime: 12.53s
% Computational Cost: add. (8072->667), mult. (22082->908), div. (0->0), fcn. (20994->6), ass. (0->319)
t250 = sin(qJ(2));
t253 = cos(qJ(2));
t390 = Icges(4,5) * t253;
t394 = Icges(3,4) * t253;
t485 = -t390 + t394 + (Icges(3,1) + Icges(4,1)) * t250;
t251 = sin(qJ(1));
t252 = cos(qJ(4));
t236 = pkin(4) * t252 + pkin(3);
t407 = pkin(3) - t236;
t334 = t407 * t253;
t249 = sin(qJ(4));
t378 = t249 * t250;
t459 = -pkin(4) * t378 + t334;
t484 = t459 * t251;
t483 = Icges(5,1) + Icges(6,1);
t482 = Icges(5,4) + Icges(6,4);
t481 = Icges(5,5) + Icges(6,5);
t480 = Icges(5,2) + Icges(6,2);
t479 = Icges(5,6) + Icges(6,6);
t298 = t252 * t253 + t378;
t377 = t249 * t253;
t299 = -t250 * t252 + t377;
t460 = -t482 * t298 - t483 * t299;
t478 = -t460 / 0.2e1;
t461 = -t480 * t298 - t482 * t299;
t477 = t461 / 0.2e1;
t351 = qJD(2) * t253;
t426 = t299 * qJD(4);
t476 = (-t249 * t351 + t426) * pkin(4);
t248 = -qJ(5) - pkin(7);
t397 = -rSges(6,3) + t248;
t475 = t485 * qJD(2);
t439 = t251 / 0.2e1;
t254 = cos(qJ(1));
t474 = -t254 / 0.2e1;
t473 = t254 / 0.2e1;
t472 = -qJD(1) / 0.2e1;
t471 = qJD(1) / 0.2e1;
t144 = -qJD(2) * t299 + t426;
t274 = t298 * qJD(4);
t145 = qJD(2) * t298 - t274;
t470 = t479 * t144 + t481 * t145;
t469 = t480 * t144 + t482 * t145;
t468 = t482 * t144 + t483 * t145;
t153 = -rSges(5,1) * t299 - rSges(5,2) * t298;
t246 = t251 ^ 2;
t247 = t254 ^ 2;
t360 = t246 + t247;
t467 = t153 * t360;
t185 = t299 * t251;
t186 = t298 * t251;
t115 = Icges(6,4) * t186 - Icges(6,2) * t185 + Icges(6,6) * t254;
t117 = Icges(5,4) * t186 - Icges(5,2) * t185 + Icges(5,6) * t254;
t466 = t115 + t117;
t187 = t299 * t254;
t188 = t298 * t254;
t116 = Icges(6,4) * t188 - Icges(6,2) * t187 - Icges(6,6) * t251;
t118 = Icges(5,4) * t188 - Icges(5,2) * t187 - Icges(5,6) * t251;
t465 = t116 + t118;
t119 = Icges(6,1) * t186 - Icges(6,4) * t185 + Icges(6,5) * t254;
t121 = Icges(5,1) * t186 - Icges(5,4) * t185 + Icges(5,5) * t254;
t464 = t119 + t121;
t120 = Icges(6,1) * t188 - Icges(6,4) * t187 - Icges(6,5) * t251;
t122 = Icges(5,1) * t188 - Icges(5,4) * t187 - Icges(5,5) * t251;
t463 = t120 + t122;
t462 = -t479 * t298 - t481 * t299;
t445 = t187 / 0.2e1;
t444 = -t188 / 0.2e1;
t442 = t298 / 0.2e1;
t440 = t299 / 0.2e1;
t113 = Icges(5,5) * t186 - Icges(5,6) * t185 + Icges(5,3) * t254;
t268 = t113 * t254 - t185 * t117 + t186 * t121;
t111 = Icges(6,5) * t186 - Icges(6,6) * t185 + Icges(6,3) * t254;
t270 = t111 * t254 - t185 * t115 + t186 * t119;
t112 = Icges(6,5) * t188 - Icges(6,6) * t187 - Icges(6,3) * t251;
t382 = t112 * t254;
t38 = -t185 * t116 + t186 * t120 + t382;
t114 = Icges(5,5) * t188 - Icges(5,6) * t187 - Icges(5,3) * t251;
t380 = t114 * t254;
t39 = -t185 * t118 + t186 * t122 + t380;
t433 = (-t268 - t270) * t254 + (t38 + t39) * t251;
t267 = -t113 * t251 - t117 * t187 + t121 * t188;
t269 = -t111 * t251 - t115 * t187 + t119 * t188;
t383 = t112 * t251;
t40 = -t116 * t187 + t120 * t188 - t383;
t381 = t114 * t251;
t41 = -t118 * t187 + t122 * t188 - t381;
t432 = (-t267 - t269) * t254 + (t40 + t41) * t251;
t313 = -Icges(3,2) * t250 + t394;
t169 = Icges(3,6) * t251 + t254 * t313;
t395 = Icges(3,4) * t250;
t317 = Icges(3,1) * t253 - t395;
t173 = Icges(3,5) * t251 + t254 * t317;
t300 = t169 * t250 - t173 * t253;
t283 = t300 * t251;
t168 = -Icges(3,6) * t254 + t251 * t313;
t172 = -Icges(3,5) * t254 + t251 * t317;
t301 = t168 * t250 - t172 * t253;
t284 = t301 * t254;
t309 = Icges(4,3) * t250 + t390;
t163 = Icges(4,6) * t251 + t254 * t309;
t391 = Icges(4,5) * t250;
t315 = Icges(4,1) * t253 + t391;
t171 = Icges(4,4) * t251 + t254 * t315;
t302 = t163 * t250 + t171 * t253;
t285 = t302 * t251;
t162 = -Icges(4,6) * t254 + t251 * t309;
t170 = -Icges(4,4) * t254 + t251 * t315;
t303 = t162 * t250 + t170 * t253;
t286 = t303 * t254;
t321 = -t186 * rSges(6,1) + t185 * rSges(6,2);
t399 = rSges(6,3) * t254;
t374 = -t321 + t399 + (-pkin(7) - t248) * t254 - t484;
t431 = t374 * t254;
t376 = t250 * t254;
t346 = pkin(4) * t376;
t375 = t253 * t254;
t430 = t188 * rSges(6,1) - t187 * rSges(6,2) + t236 * t375 + t249 * t346 + t397 * t251;
t356 = qJD(1) * t251;
t100 = -t144 * t254 - t298 * t356;
t350 = qJD(2) * t254;
t336 = t253 * t350;
t355 = qJD(1) * t254;
t410 = pkin(4) * t249;
t99 = t145 * t254 + t299 * t356;
t429 = qJD(4) * t252 * t346 + t100 * rSges(6,1) + t99 * rSges(6,2) - qJD(5) * t251 + t248 * t355 + t336 * t410;
t428 = -rSges(3,2) * t376 + t251 * rSges(3,3);
t310 = Icges(3,5) * t253 - Icges(3,6) * t250;
t164 = -Icges(3,3) * t254 + t251 * t310;
t311 = Icges(4,4) * t253 + Icges(4,6) * t250;
t166 = -Icges(4,2) * t254 + t251 * t311;
t233 = pkin(3) * t375;
t207 = -pkin(7) * t251 + t233;
t373 = -t207 + t430;
t33 = -t251 * t374 - t254 * t373;
t101 = -qJD(1) * t187 + t145 * t251;
t102 = qJD(1) * t188 - t144 * t251;
t322 = -t102 * rSges(6,1) - t101 * rSges(6,2);
t347 = pkin(4) * t377;
t331 = qJD(4) * t347;
t353 = qJD(2) * t250;
t352 = qJD(2) * t251;
t338 = t250 * t352;
t362 = pkin(3) * t338 + pkin(7) * t356;
t409 = pkin(7) * t254;
t425 = -(-rSges(6,3) * t356 - t322 + (-t459 * qJD(1) + qJD(5)) * t254 + (qJD(1) * t248 - t236 * t353 - t476) * t251 + t362) * t251 - (-rSges(6,3) * t355 + (t353 * t407 - t331) * t254 + (t409 + t484) * qJD(1) + t429) * t254;
t424 = 2 * m(3);
t423 = 2 * m(4);
t422 = 2 * m(5);
t421 = 2 * m(6);
t420 = m(4) / 0.2e1;
t419 = m(5) / 0.2e1;
t418 = m(6) / 0.2e1;
t417 = -pkin(2) - pkin(3);
t260 = Icges(6,5) * t102 + Icges(6,6) * t101 - Icges(6,3) * t356;
t262 = Icges(6,4) * t102 + Icges(6,2) * t101 - Icges(6,6) * t356;
t264 = Icges(6,1) * t102 + Icges(6,4) * t101 - Icges(6,5) * t356;
t47 = Icges(6,5) * t100 + Icges(6,6) * t99 - Icges(6,3) * t355;
t49 = Icges(6,4) * t100 + Icges(6,2) * t99 - Icges(6,6) * t355;
t51 = Icges(6,1) * t100 + Icges(6,4) * t99 - Icges(6,5) * t355;
t1 = (t100 * t120 + t99 * t116 - t187 * t49 + t188 * t51 - t251 * t47) * t251 - (t100 * t119 - t111 * t355 + t99 * t115 - t187 * t262 + t188 * t264 - t251 * t260) * t254 + (t40 * t254 + (t269 - t382) * t251) * qJD(1);
t261 = Icges(5,5) * t102 + Icges(5,6) * t101 - Icges(5,3) * t356;
t263 = Icges(5,4) * t102 + Icges(5,2) * t101 - Icges(5,6) * t356;
t265 = Icges(5,1) * t102 + Icges(5,4) * t101 - Icges(5,5) * t356;
t48 = Icges(5,5) * t100 + Icges(5,6) * t99 - Icges(5,3) * t355;
t50 = Icges(5,4) * t100 + Icges(5,2) * t99 - Icges(5,6) * t355;
t52 = Icges(5,1) * t100 + Icges(5,4) * t99 - Icges(5,5) * t355;
t2 = (t100 * t122 + t99 * t118 - t187 * t50 + t188 * t52 - t251 * t48) * t251 - (t100 * t121 - t113 * t355 + t99 * t117 - t187 * t263 + t188 * t265 - t251 * t261) * t254 + (t41 * t254 + (t267 - t380) * t251) * qJD(1);
t416 = t2 + t1;
t3 = (t101 * t116 + t102 * t120 - t185 * t49 + t186 * t51 + t254 * t47) * t251 - (t101 * t115 + t102 * t119 - t111 * t356 - t185 * t262 + t186 * t264 + t254 * t260) * t254 + (t38 * t254 + (t270 - t383) * t251) * qJD(1);
t4 = (t101 * t118 + t102 * t122 - t185 * t50 + t186 * t52 + t254 * t48) * t251 - (t101 * t117 + t102 * t121 - t113 * t356 - t185 * t263 + t186 * t265 + t254 * t261) * t254 + (t39 * t254 + (t268 - t381) * t251) * qJD(1);
t415 = -t4 - t3;
t414 = -rSges(4,1) - pkin(2);
t413 = -rSges(5,3) - pkin(7);
t219 = rSges(3,1) * t250 + rSges(3,2) * t253;
t412 = m(3) * t219;
t411 = m(5) * t153;
t408 = -pkin(2) - t236;
t403 = t100 * rSges(5,1) + t99 * rSges(5,2);
t402 = rSges(4,1) * t250;
t401 = rSges(4,2) * t254;
t400 = rSges(3,3) * t254;
t242 = t251 * rSges(4,2);
t398 = -rSges(4,3) - qJ(3);
t396 = rSges(6,1) * t145 + rSges(6,2) * t144 - qJD(2) * t334 + (t249 * t353 - t274) * pkin(4);
t385 = qJ(3) * t250;
t384 = qJ(3) * t253;
t323 = -t186 * rSges(5,1) + t185 * rSges(5,2);
t124 = rSges(5,3) * t254 - t323;
t379 = t124 * t254;
t372 = -rSges(6,1) * t299 - rSges(6,2) * t298 - t250 * t407 - t347;
t371 = t188 * rSges(5,1) - t187 * rSges(5,2);
t320 = pkin(2) * t253 + t385;
t193 = t320 * t251;
t194 = pkin(2) * t375 + qJ(3) * t376;
t370 = t251 * t193 + t254 * t194;
t189 = qJD(2) * t320 - qJD(3) * t253;
t325 = rSges(4,1) * t253 + rSges(4,3) * t250;
t369 = -t325 * qJD(2) - t189;
t368 = -t194 - t207;
t217 = pkin(2) * t250 - t384;
t195 = t217 * t356;
t341 = t250 * t356;
t367 = pkin(3) * t341 + t195;
t218 = -rSges(4,3) * t253 + t402;
t366 = -t217 - t218;
t349 = qJD(3) * t250;
t365 = qJ(3) * t336 + t254 * t349;
t364 = rSges(4,2) * t355 + rSges(4,3) * t336;
t363 = rSges(3,2) * t341 + rSges(3,3) * t355;
t361 = t254 * pkin(1) + t251 * pkin(6);
t359 = qJD(1) * t153;
t165 = Icges(3,3) * t251 + t254 * t310;
t358 = qJD(1) * t165;
t167 = Icges(4,2) * t251 + t254 * t311;
t357 = qJD(1) * t167;
t225 = pkin(2) * t338;
t337 = t250 * t350;
t271 = -t253 * t356 - t337;
t340 = t253 * t355;
t344 = t251 * (pkin(2) * t340 + t251 * t349 - t225 + (t250 * t355 + t351 * t251) * qJ(3)) + t254 * (pkin(2) * t271 - qJ(3) * t341 + t365) + t193 * t355;
t240 = pkin(6) * t355;
t342 = t240 + t365;
t181 = rSges(4,1) * t375 + rSges(4,3) * t376 + t242;
t335 = (t310 + t311) * qJD(2) / 0.2e1;
t333 = -pkin(3) * t250 - t217;
t158 = t366 * t254;
t332 = qJD(1) * t372;
t206 = t251 * t253 * pkin(3) + t409;
t330 = t251 * t206 + t254 * t207 + t370;
t329 = t361 + t194;
t328 = -t153 + t333;
t327 = -pkin(3) * t351 - t189;
t326 = rSges(3,1) * t253 - rSges(3,2) * t250;
t324 = -t102 * rSges(5,1) - t101 * rSges(5,2);
t319 = -t251 * (-rSges(5,3) * t356 - t324) - t254 * (-rSges(5,3) * t355 + t403);
t244 = t254 * pkin(6);
t258 = -pkin(1) + t408 * t253 + (-qJ(3) - t410) * t250;
t257 = t258 * t251;
t68 = t254 * t397 + t244 + t257 + t321;
t69 = t329 + t430;
t318 = t251 * t69 + t254 * t68;
t312 = Icges(3,2) * t253 + t395;
t308 = -Icges(4,3) * t253 + t391;
t126 = -rSges(5,3) * t251 + t371;
t67 = -t251 * t124 - t126 * t254;
t94 = rSges(5,1) * t145 + rSges(5,2) * t144;
t297 = t327 - t94;
t296 = t333 - t372;
t182 = rSges(3,1) * t375 + t428;
t295 = -pkin(1) - t326;
t294 = -t461 * t99 / 0.2e1 + t100 * t478 - t465 * t144 / 0.2e1 - t463 * t145 / 0.2e1 + t469 * t445 + t468 * t444 + (t49 + t50) * t442 + (t51 + t52) * t440 + t470 * t439 + t462 * t355 / 0.2e1;
t293 = t101 * t477 + t460 * t102 / 0.2e1 + t466 * t144 / 0.2e1 + t464 * t145 / 0.2e1 - t469 * t185 / 0.2e1 + t468 * t186 / 0.2e1 - (t262 + t263) * t298 / 0.2e1 - (t264 + t265) * t299 / 0.2e1 + t470 * t473 - t462 * t356 / 0.2e1;
t292 = t185 * t477 + t186 * t478 + t464 * t440 + t466 * t442 + t462 * t474;
t291 = t462 * t439 + t463 * t440 + t465 * t442 + t460 * t444 + t461 * t445;
t106 = t328 * t254;
t289 = t251 * (pkin(3) * t340 - t362) + t254 * (pkin(3) * t271 - pkin(7) * t355) + t206 * t355 + t344;
t288 = t327 - t396;
t287 = qJD(2) * t219;
t280 = qJD(2) * t312;
t279 = qJD(2) * (-Icges(4,4) * t250 + Icges(4,6) * t253);
t278 = qJD(2) * (-Icges(3,5) * t250 - Icges(3,6) * t253);
t277 = qJD(2) * t308;
t76 = t296 * t254;
t273 = t253 * t417 - pkin(1) - t385;
t266 = t250 * t398 + t253 * t414 - pkin(1);
t259 = t266 * t251;
t256 = t251 * t273 + t254 * t413;
t205 = t326 * qJD(2);
t180 = t251 * t326 - t400;
t179 = t251 * t325 - t401;
t157 = t366 * t251;
t156 = t182 + t361;
t155 = t251 * t295 + t244 + t400;
t135 = t251 * t279 + t357;
t134 = -qJD(1) * t166 + t254 * t279;
t133 = t251 * t278 + t358;
t132 = -qJD(1) * t164 + t254 * t278;
t128 = t329 + t181;
t127 = t244 + t259 + t401;
t108 = t219 * t352 + ((-rSges(3,3) - pkin(6)) * t251 + t295 * t254) * qJD(1);
t107 = rSges(3,1) * t271 - rSges(3,2) * t336 - pkin(1) * t356 + t240 + t363;
t105 = t328 * t251;
t104 = t372 * t254;
t103 = t372 * t251;
t86 = qJD(1) * t158 + t251 * t369;
t85 = t218 * t356 + t254 * t369 + t195;
t84 = t251 * t165 - t300 * t254;
t83 = t251 * t164 - t284;
t82 = t251 * t167 + t302 * t254;
t81 = t251 * t166 + t286;
t80 = -t165 * t254 - t283;
t79 = -t164 * t254 - t251 * t301;
t78 = -t167 * t254 + t285;
t77 = -t166 * t254 + t251 * t303;
t75 = t296 * t251;
t74 = t251 * t413 + t233 + t329 + t371;
t73 = t244 + t256 + t323;
t72 = t251 * t179 + t181 * t254 + t370;
t71 = t225 + (-t349 + (t253 * t398 + t402) * qJD(2)) * t251 + ((-rSges(4,2) - pkin(6)) * t251 + t266 * t254) * qJD(1);
t70 = qJD(1) * t259 + t337 * t414 + t342 + t364;
t46 = qJD(1) * t106 + t251 * t297;
t45 = t153 * t356 + t254 * t297 + t367;
t44 = -t67 + t330;
t43 = t251 * t396 + t254 * t332;
t42 = t254 * t396 - t356 * t372;
t32 = t225 + (-qJ(3) * t351 - t349) * t251 + ((rSges(5,3) - pkin(6)) * t251 + t273 * t254) * qJD(1) + t324 + t362;
t31 = qJD(1) * t256 + t337 * t417 + t342 + t403;
t30 = qJD(1) * t76 + t251 * t288;
t29 = t251 * t332 + t254 * t288 + t367;
t28 = -t33 + t330;
t27 = t254 * t364 + (-t218 * t246 - t247 * t402) * qJD(2) + (t254 * t179 + (-t181 - t194 + t242) * t251) * qJD(1) + t344;
t26 = t225 + (qJD(1) * t258 - qJD(5)) * t254 + (-t349 + (t236 * t250 - t384) * qJD(2) + (-pkin(6) - t397) * qJD(1) + t476) * t251 + t322;
t25 = (t353 * t408 - t331) * t254 + (t257 - t399) * qJD(1) + t342 + t429;
t20 = (t126 * t251 - t379) * qJD(1) + t319;
t11 = (t379 + (-t126 + t368) * t251) * qJD(1) + t289 - t319;
t10 = (t251 * t373 - t431) * qJD(1) + t425;
t5 = (t431 + (t368 - t373) * t251) * qJD(1) + t289 - t425;
t6 = [(t25 * t69 + t26 * t68) * t421 + (t31 * t74 + t32 * t73) * t422 + (t127 * t71 + t128 * t70) * t423 + (t107 * t156 + t108 * t155) * t424 - t468 * t299 - t469 * t298 + t460 * t145 + t461 * t144 + (t317 - t312 + t315 + t308) * t353 + (t313 - t309 + t485) * t351; m(4) * (t127 * t85 + t128 * t86 + t157 * t70 + t158 * t71) + m(5) * (t105 * t31 + t106 * t32 + t45 * t73 + t46 * t74) + m(6) * (t25 * t75 + t26 * t76 + t29 * t68 + t30 * t69) + (m(3) * (-t108 * t219 - t155 * t205) + t335 * t254 + (t163 * t471 - t251 * t277 / 0.2e1 + t169 * t472 + t280 * t439) * t253 + ((t171 + t173) * t472 + t475 * t439) * t250 - t293) * t254 + (m(3) * (-t107 * t219 - t156 * t205) + t335 * t251 + (t162 * t471 + t168 * t472 + t277 * t473 + t280 * t474) * t253 + (t475 * t474 + (t170 + t172) * t472) * t250 - t294) * t251 + (-t286 / 0.2e1 + t284 / 0.2e1 + t285 / 0.2e1 - t283 / 0.2e1) * qJD(2) + ((-t156 * t412 + (-t163 / 0.2e1 + t169 / 0.2e1) * t253 + (t171 / 0.2e1 + t173 / 0.2e1) * t250 - t291) * t254 + (t155 * t412 + (-t162 / 0.2e1 + t168 / 0.2e1) * t253 + (t170 / 0.2e1 + t172 / 0.2e1) * t250 - t292) * t251) * qJD(1); (t28 * t5 + t29 * t76 + t30 * t75) * t421 + (t105 * t46 + t106 * t45 + t11 * t44) * t422 - t254 * t3 + t251 * t1 - t254 * t4 + t251 * t2 + (t157 * t86 + t158 * t85 + t27 * t72) * t423 + t251 * ((t251 * t134 + (t81 - t285) * qJD(1)) * t251 + (t82 * qJD(1) + (-t162 * t351 + t170 * t353) * t254 + (-t135 + (t163 * t253 - t171 * t250) * qJD(2) + (t167 + t303) * qJD(1)) * t251) * t254) - t254 * ((t135 * t254 + (t78 - t286) * qJD(1)) * t254 + (t77 * qJD(1) + (t163 * t351 - t171 * t353 + t357) * t251 + (-t134 + (-t162 * t253 + t170 * t250) * qJD(2) + t302 * qJD(1)) * t254) * t251) + t251 * ((t251 * t132 + (t83 + t283) * qJD(1)) * t251 + (t84 * qJD(1) + (t168 * t351 + t172 * t353) * t254 + (-t133 + (-t169 * t253 - t173 * t250) * qJD(2) + (t165 - t301) * qJD(1)) * t251) * t254) - t254 * ((t133 * t254 + (t80 + t284) * qJD(1)) * t254 + (t79 * qJD(1) + (-t169 * t351 - t173 * t353 + t358) * t251 + (-t132 + (t168 * t253 + t172 * t250) * qJD(2) - t300 * qJD(1)) * t254) * t251) + ((t251 * t180 + t182 * t254) * ((qJD(1) * t180 - t254 * t287 + t363) * t254 + (-t251 * t287 + (-t182 + t428) * qJD(1)) * t251) + t360 * t219 * t205) * t424 + ((-t77 - t79) * t254 + (t78 + t80) * t251 + t433) * t356 + ((-t81 - t83) * t254 + (t82 + t84) * t251 + t432) * t355; 0.2e1 * (t318 * t418 + (t251 * t74 + t254 * t73) * t419 + (t127 * t254 + t128 * t251) * t420) * t351 + 0.2e1 * ((t25 * t251 + t254 * t26 + t355 * t69 - t356 * t68) * t418 + (t251 * t31 + t254 * t32 + t355 * t74 - t356 * t73) * t419 + (-t127 * t356 + t128 * t355 + t251 * t70 + t254 * t71) * t420) * t250; 0.2e1 * ((t350 * t76 + t352 * t75 - t5) * t418 + (t105 * t352 + t106 * t350 - t11) * t419 + (t157 * t352 + t158 * t350 - t27) * t420) * t253 + 0.2e1 * ((qJD(2) * t28 + t251 * t30 + t254 * t29 + t355 * t75 - t356 * t76) * t418 + (qJD(2) * t44 + t105 * t355 - t106 * t356 + t251 * t46 + t254 * t45) * t419 + (qJD(2) * t72 + t157 * t355 - t158 * t356 + t251 * t86 + t254 * t85) * t420) * t250; 0.4e1 * (t420 + t419 + t418) * (-0.1e1 + t360) * t250 * t351; m(6) * (t103 * t25 + t104 * t26 + t42 * t68 + t43 * t69) + (m(5) * (t153 * t32 + t73 * t94) + (t411 * t74 + t291) * qJD(1) + t293) * t254 + (m(5) * (t153 * t31 + t74 * t94) + (-t411 * t73 + t292) * qJD(1) + t294) * t251; m(6) * (t10 * t28 + t103 * t30 + t104 * t29 + t33 * t5 + t42 * t76 + t43 * t75) + m(5) * (t67 * t11 + t20 * t44) + (m(5) * (t105 * t359 + t106 * t94 + t153 * t45) - t415 - t432 * qJD(1)) * t254 + (m(5) * (t105 * t94 - t106 * t359 + t153 * t46) - t416 - t433 * qJD(1)) * t251; 0.2e1 * ((qJD(2) * t467 - t20) * t419 + (t103 * t352 + t104 * t350 - t10) * t418) * t253 + 0.2e1 * ((qJD(2) * t67 + t360 * t94) * t419 + (qJD(2) * t33 + t103 * t355 - t104 * t356 + t251 * t43 + t254 * t42) * t418) * t250; t415 * t254 + t416 * t251 + (t20 * t67 + t94 * t467) * t422 + (t10 * t33 + t103 * t43 + t104 * t42) * t421 + (t251 * t433 + t254 * t432) * qJD(1); m(6) * (-qJD(1) * t318 + t25 * t254 - t251 * t26); m(6) * (-t251 * t29 + t254 * t30 + (-t251 * t75 - t254 * t76) * qJD(1)); 0; m(6) * (-t251 * t42 + t254 * t43 + (-t103 * t251 - t104 * t254) * qJD(1)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;

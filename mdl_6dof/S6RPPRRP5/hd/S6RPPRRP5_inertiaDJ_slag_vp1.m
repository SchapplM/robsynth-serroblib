% Calculate time derivative of joint inertia matrix for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:07:42
% EndTime: 2019-03-09 02:08:02
% DurationCPUTime: 13.23s
% Computational Cost: add. (10913->755), mult. (28165->1047), div. (0->0), fcn. (26518->6), ass. (0->356)
t244 = sin(qJ(5));
t245 = sin(qJ(4));
t248 = cos(qJ(4));
t247 = cos(qJ(5));
t399 = Icges(7,4) * t247;
t289 = -Icges(7,2) * t244 + t399;
t175 = Icges(7,6) * t245 + t248 * t289;
t401 = Icges(6,4) * t247;
t290 = -Icges(6,2) * t244 + t401;
t176 = Icges(6,6) * t245 + t248 * t290;
t458 = -t176 - t175;
t460 = t458 * t244;
t286 = Icges(7,5) * t247 - Icges(7,6) * t244;
t171 = Icges(7,3) * t245 + t248 * t286;
t287 = Icges(6,5) * t247 - Icges(6,6) * t244;
t172 = Icges(6,3) * t245 + t248 * t287;
t459 = t172 + t171;
t400 = Icges(7,4) * t244;
t292 = Icges(7,1) * t247 - t400;
t179 = Icges(7,5) * t245 + t248 * t292;
t402 = Icges(6,4) * t244;
t293 = Icges(6,1) * t247 - t402;
t180 = Icges(6,5) * t245 + t248 * t293;
t454 = -t180 - t179;
t356 = qJD(4) * t245;
t457 = qJD(1) * t248 / 0.2e1;
t246 = sin(qJ(1));
t249 = cos(qJ(1));
t381 = t244 * t249;
t232 = pkin(5) * t247 + pkin(4);
t416 = pkin(4) - t232;
t328 = t416 * t245;
t243 = -qJ(6) - pkin(8);
t383 = t243 * t248;
t452 = t328 - t383;
t456 = -pkin(5) * t381 + t452 * t246;
t405 = rSges(7,3) - t243;
t455 = (-t454 * t247 + t460) * t248 + t459 * t245;
t352 = qJD(5) * t248;
t140 = (-Icges(7,5) * t244 - Icges(7,6) * t247) * t352 + (Icges(7,3) * t248 - t245 * t286) * qJD(4);
t141 = (-Icges(6,5) * t244 - Icges(6,6) * t247) * t352 + (Icges(6,3) * t248 - t245 * t287) * qJD(4);
t148 = (-Icges(7,1) * t244 - t399) * t352 + (Icges(7,5) * t248 - t245 * t292) * qJD(4);
t149 = (-Icges(6,1) * t244 - t401) * t352 + (Icges(6,5) * t248 - t245 * t293) * qJD(4);
t354 = qJD(4) * t248;
t453 = t459 * t354 - t356 * t460 + (t140 + t141) * t245 + ((t148 + t149) * t248 + t454 * t356) * t247;
t353 = qJD(4) * t249;
t333 = t245 * t353;
t358 = qJD(1) * t246;
t255 = t248 * t358 + t333;
t144 = (-Icges(7,2) * t247 - t400) * t352 + (Icges(7,6) * t248 - t245 * t289) * qJD(4);
t145 = (-Icges(6,2) * t247 - t402) * t352 + (Icges(6,6) * t248 - t245 * t290) * qJD(4);
t451 = (-t145 - t144) * t244;
t450 = t458 * t247;
t377 = t246 * t247;
t380 = t245 * t249;
t196 = -t244 * t380 - t377;
t374 = t247 * t249;
t378 = t246 * t244;
t197 = t245 * t374 - t378;
t373 = t248 * t249;
t123 = Icges(7,5) * t197 + Icges(7,6) * t196 - Icges(7,3) * t373;
t127 = Icges(7,4) * t197 + Icges(7,2) * t196 - Icges(7,6) * t373;
t131 = Icges(7,1) * t197 + Icges(7,4) * t196 - Icges(7,5) * t373;
t282 = t127 * t244 - t131 * t247;
t62 = t123 * t245 - t248 * t282;
t125 = Icges(6,5) * t197 + Icges(6,6) * t196 - Icges(6,3) * t373;
t129 = Icges(6,4) * t197 + Icges(6,2) * t196 - Icges(6,6) * t373;
t133 = Icges(6,1) * t197 + Icges(6,4) * t196 - Icges(6,5) * t373;
t280 = t129 * t244 - t133 * t247;
t64 = t125 * t245 - t248 * t280;
t411 = t62 + t64;
t194 = -t245 * t378 + t374;
t195 = t245 * t377 + t381;
t376 = t246 * t248;
t122 = Icges(7,5) * t195 + Icges(7,6) * t194 - Icges(7,3) * t376;
t126 = Icges(7,4) * t195 + Icges(7,2) * t194 - Icges(7,6) * t376;
t130 = Icges(7,1) * t195 + Icges(7,4) * t194 - Icges(7,5) * t376;
t283 = t126 * t244 - t130 * t247;
t61 = t122 * t245 - t248 * t283;
t124 = Icges(6,5) * t195 + Icges(6,6) * t194 - Icges(6,3) * t376;
t128 = Icges(6,4) * t195 + Icges(6,2) * t194 - Icges(6,6) * t376;
t132 = Icges(6,1) * t195 + Icges(6,4) * t194 - Icges(6,5) * t376;
t281 = t128 * t244 - t132 * t247;
t63 = t124 * t245 - t248 * t281;
t412 = t61 + t63;
t449 = -t246 * t412 - t249 * t411;
t448 = t246 * t411 - t249 * t412;
t320 = qJD(1) * t245 + qJD(5);
t321 = qJD(5) * t245 + qJD(1);
t332 = t246 * t354;
t117 = -t321 * t377 + (-t249 * t320 - t332) * t244;
t273 = t321 * t244;
t118 = t320 * t374 + (t247 * t354 - t273) * t246;
t355 = qJD(4) * t246;
t335 = t245 * t355;
t357 = qJD(1) * t249;
t337 = t248 * t357;
t256 = t335 - t337;
t69 = Icges(7,5) * t118 + Icges(7,6) * t117 + Icges(7,3) * t256;
t73 = Icges(7,4) * t118 + Icges(7,2) * t117 + Icges(7,6) * t256;
t77 = Icges(7,1) * t118 + Icges(7,4) * t117 + Icges(7,5) * t256;
t19 = (qJD(4) * t283 + t69) * t245 + (qJD(4) * t122 - t244 * t73 + t247 * t77 + (-t126 * t247 - t130 * t244) * qJD(5)) * t248;
t71 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t256;
t75 = Icges(6,4) * t118 + Icges(6,2) * t117 + Icges(6,6) * t256;
t79 = Icges(6,1) * t118 + Icges(6,4) * t117 + Icges(6,5) * t256;
t21 = (qJD(4) * t281 + t71) * t245 + (qJD(4) * t124 - t244 * t75 + t247 * t79 + (-t128 * t247 - t132 * t244) * qJD(5)) * t248;
t447 = t19 + t21;
t272 = t321 * t249;
t331 = t248 * t353;
t438 = t246 * t320 - t331;
t115 = t244 * t438 - t247 * t272;
t116 = -t244 * t272 - t247 * t438;
t68 = Icges(7,5) * t116 + Icges(7,6) * t115 + Icges(7,3) * t255;
t72 = Icges(7,4) * t116 + Icges(7,2) * t115 + Icges(7,6) * t255;
t76 = Icges(7,1) * t116 + Icges(7,4) * t115 + Icges(7,5) * t255;
t20 = (qJD(4) * t282 + t68) * t245 + (qJD(4) * t123 - t244 * t72 + t247 * t76 + (-t127 * t247 - t131 * t244) * qJD(5)) * t248;
t70 = Icges(6,5) * t116 + Icges(6,6) * t115 + Icges(6,3) * t255;
t74 = Icges(6,4) * t116 + Icges(6,2) * t115 + Icges(6,6) * t255;
t78 = Icges(6,1) * t116 + Icges(6,4) * t115 + Icges(6,5) * t255;
t22 = (qJD(4) * t280 + t70) * t245 + (qJD(4) * t125 - t244 * t74 + t247 * t78 + (-t129 * t247 - t133 * t244) * qJD(5)) * t248;
t446 = -t20 - t22;
t53 = -t124 * t376 + t128 * t194 + t132 * t195;
t54 = -t125 * t376 + t129 * t194 + t133 * t195;
t303 = t246 * t53 + t249 * t54;
t51 = -t122 * t376 + t126 * t194 + t130 * t195;
t52 = -t123 * t376 + t127 * t194 + t131 * t195;
t305 = t246 * t51 + t249 * t52;
t84 = -t171 * t376 + t175 * t194 + t179 * t195;
t85 = -t172 * t376 + t176 * t194 + t180 * t195;
t445 = (-t303 - t305) * t248 + (t84 + t85) * t245;
t57 = -t124 * t373 + t196 * t128 + t197 * t132;
t58 = -t125 * t373 + t196 * t129 + t197 * t133;
t299 = t246 * t57 + t249 * t58;
t55 = -t122 * t373 + t196 * t126 + t197 * t130;
t56 = -t123 * t373 + t196 * t127 + t197 * t131;
t301 = t246 * t55 + t249 * t56;
t86 = -t171 * t373 + t196 * t175 + t197 * t179;
t87 = -t172 * t373 + t196 * t176 + t197 * t180;
t413 = (-t299 - t301) * t248 + (t86 + t87) * t245;
t406 = pkin(5) * qJD(5);
t346 = t247 * t406;
t351 = qJD(6) * t248;
t444 = t232 * t331 + (-t351 + (-qJD(4) * t243 - t244 * t406) * t245) * t249 - t246 * t346 + t116 * rSges(7,1) + t115 * rSges(7,2) + t255 * rSges(7,3);
t310 = rSges(7,1) * t247 - rSges(7,2) * t244;
t330 = t244 * t352;
t415 = -pkin(8) - t243;
t372 = -pkin(5) * t330 + qJD(6) * t245 + (-rSges(7,1) * t244 - rSges(7,2) * t247) * t352 + (-t245 * t310 + t328 + (rSges(7,3) + t415) * t248) * qJD(4);
t443 = t246 * t372;
t442 = t197 * rSges(7,1) + t196 * rSges(7,2) + t232 * t380 - t405 * t373;
t241 = t246 ^ 2;
t242 = t249 ^ 2;
t360 = t241 + t242;
t220 = rSges(5,1) * t331;
t316 = rSges(5,1) * t245 + rSges(5,2) * t248;
t417 = -pkin(1) - qJ(3);
t263 = -t316 + t417;
t426 = -rSges(5,3) - pkin(7);
t251 = t246 * t263 + t249 * t426;
t363 = qJ(2) * t357 + qJD(2) * t246;
t341 = qJD(3) * t249 + t363;
t103 = -rSges(5,2) * t333 + qJD(1) * t251 + t220 + t341;
t408 = rSges(5,2) * t245;
t215 = rSges(5,1) * t248 - t408;
t236 = qJD(2) * t249;
t362 = pkin(7) * t358 + t236;
t104 = (-t215 * qJD(4) - qJD(3)) * t246 + ((rSges(5,3) - qJ(2)) * t246 + t263 * t249) * qJD(1) + t362;
t441 = t246 * t103 + t104 * t249;
t239 = t249 * qJ(2);
t156 = t239 + t251;
t237 = t246 * qJ(2);
t361 = t249 * pkin(1) + t237;
t340 = t249 * qJ(3) + t361;
t364 = rSges(5,1) * t380 + rSges(5,2) * t373;
t157 = t246 * t426 + t340 + t364;
t440 = t156 * t249 + t157 * t246;
t288 = Icges(5,5) * t245 + Icges(5,6) * t248;
t173 = Icges(5,3) * t249 + t246 * t288;
t439 = qJD(1) * t173;
t404 = Icges(5,4) * t245;
t291 = Icges(5,2) * t248 + t404;
t177 = Icges(5,6) * t249 + t246 * t291;
t403 = Icges(5,4) * t248;
t294 = Icges(5,1) * t245 + t403;
t181 = Icges(5,5) * t249 + t246 * t294;
t422 = pkin(4) * t248;
t216 = pkin(8) * t245 + t422;
t207 = t246 * t216;
t166 = t245 * t415 - t248 * t416;
t368 = rSges(7,3) * t245 + t248 * t310 + t166;
t324 = t368 * t246;
t107 = t207 + t324;
t208 = t249 * t216;
t323 = t368 * t249;
t108 = t208 + t323;
t423 = pkin(4) * t245;
t210 = (pkin(8) * t248 - t423) * qJD(4);
t199 = t249 * t210;
t325 = t372 * t249;
t49 = t199 + t325 + (-t216 - t368) * t358;
t366 = t246 * t210 + t216 * t357;
t50 = qJD(1) * t323 + t366 + t443;
t437 = qJD(1) * (t107 * t249 - t108 * t246) + t50 * t246 + t249 * t49;
t421 = pkin(5) * t244;
t339 = -pkin(7) - t421;
t317 = t339 * t249;
t318 = -t232 * t245 + t417;
t43 = (t317 + (t318 - t383) * t246) * qJD(1) + t341 + t444;
t254 = t248 * t405 + t318;
t312 = t118 * rSges(7,1) + t117 * rSges(7,2);
t44 = (qJD(1) * t254 - t346) * t249 + (-qJ(2) * qJD(1) + t351 - qJD(3) + t321 * t421 + (-t232 * t248 - t245 * t405) * qJD(4)) * t246 - t312 + t362;
t311 = -t195 * rSges(7,1) - t194 * rSges(7,2);
t91 = t246 * t254 + t239 + t311 + t317;
t92 = t246 * t339 + t340 + t442;
t436 = qJD(1) * (t246 * t91 - t249 * t92) - t246 * t43 - t249 * t44;
t322 = qJD(1) * t368;
t230 = pkin(4) * t380;
t201 = -pkin(8) * t373 + t230;
t370 = -pkin(5) * t378 - t201 + t442;
t342 = pkin(4) * t331 + t255 * pkin(8);
t410 = t456 * qJD(1) - t342 + t444;
t23 = (-t353 * t368 + t410) * t245 + (qJD(4) * t370 - t246 * t322 + t325) * t248;
t229 = pkin(8) * t376;
t371 = -rSges(7,3) * t376 + t229 - t311 - t456;
t226 = pkin(8) * t337;
t409 = -t226 - (-qJD(1) * t452 + t346) * t249 - (-pkin(5) * t273 + qJD(4) * t166 - t351) * t246 - t256 * rSges(7,3) - t312;
t24 = (qJD(4) * t324 + t409) * t245 + (-qJD(4) * t371 - t249 * t322 - t443) * t248;
t59 = -t245 * t371 - t368 * t376;
t60 = t245 * t370 + t248 * t323;
t435 = -(t246 * t59 - t249 * t60) * qJD(1) + t23 * t246 + t24 * t249;
t434 = t246 * t370 - t249 * t371;
t433 = 2 * m(5);
t432 = 2 * m(6);
t431 = 2 * m(7);
t430 = t245 / 0.2e1;
t427 = rSges(3,2) - pkin(1);
t425 = -rSges(6,3) - pkin(8);
t424 = m(5) * t215;
t420 = pkin(7) * t249;
t407 = rSges(6,3) * t248;
t313 = rSges(6,1) * t247 - rSges(6,2) * t244;
t153 = (-rSges(6,1) * t244 - rSges(6,2) * t247) * t352 + (-t245 * t313 + t407) * qJD(4);
t394 = t153 * t249;
t391 = t177 * t245;
t390 = t177 * t248;
t178 = -Icges(5,6) * t246 + t249 * t291;
t389 = t178 * t245;
t388 = t178 * t248;
t387 = t181 * t245;
t386 = t181 * t248;
t182 = -Icges(5,5) * t246 + t249 * t294;
t385 = t182 * t245;
t384 = t182 * t248;
t314 = -t195 * rSges(6,1) - t194 * rSges(6,2);
t137 = -rSges(6,3) * t376 - t314;
t348 = t246 * t423;
t200 = -t229 + t348;
t369 = -t137 - t200;
t367 = t197 * rSges(6,1) + t196 * rSges(6,2);
t174 = -Icges(5,3) * t246 + t249 * t288;
t359 = qJD(1) * t174;
t350 = -rSges(4,3) + t417;
t347 = m(7) * t356;
t304 = t54 * t246 - t249 * t53;
t306 = t52 * t246 - t249 * t51;
t345 = -t306 / 0.2e1 - t304 / 0.2e1;
t300 = t58 * t246 - t249 * t57;
t302 = t56 * t246 - t249 * t55;
t344 = -t302 / 0.2e1 - t300 / 0.2e1;
t343 = -t200 - t371;
t327 = t455 * t354 + ((t451 + (t454 * t244 + t450) * qJD(5)) * t248 + t453) * t245;
t209 = t316 * qJD(4);
t326 = t209 * t360;
t81 = t116 * rSges(6,1) + t115 * rSges(6,2) + t255 * rSges(6,3);
t319 = t417 - t423;
t315 = t118 * rSges(6,1) + t117 * rSges(6,2);
t298 = t246 * t60 + t249 * t59;
t296 = t246 * t92 + t249 * t91;
t284 = t107 * t246 + t108 * t249;
t139 = -rSges(6,3) * t373 + t367;
t279 = t137 * t249 - t139 * t246;
t276 = t387 + t390;
t275 = t385 + t388;
t33 = t117 * t175 + t118 * t179 - t140 * t376 + t194 * t144 + t195 * t148 + t171 * t256;
t34 = t117 * t176 + t118 * t180 - t141 * t376 + t194 * t145 + t195 * t149 + t172 * t256;
t271 = -t19 / 0.2e1 - t33 / 0.2e1 - t34 / 0.2e1 - t21 / 0.2e1;
t31 = t115 * t175 + t116 * t179 - t140 * t373 + t196 * t144 + t197 * t148 + t171 * t255;
t32 = t115 * t176 + t116 * t180 - t141 * t373 + t196 * t145 + t197 * t149 + t172 * t255;
t270 = -t20 / 0.2e1 - t31 / 0.2e1 - t32 / 0.2e1 - t22 / 0.2e1;
t269 = t84 / 0.2e1 + t85 / 0.2e1 + t63 / 0.2e1 + t61 / 0.2e1;
t268 = t86 / 0.2e1 + t87 / 0.2e1 + t64 / 0.2e1 + t62 / 0.2e1;
t267 = rSges(3,3) * t249 + t246 * t427;
t185 = rSges(6,3) * t245 + t248 * t313;
t266 = -t246 * t153 - t185 * t357;
t264 = t319 + t407;
t262 = t276 * t249;
t261 = t275 * t246;
t259 = qJD(4) * (-Icges(5,2) * t245 + t403);
t258 = qJD(4) * (Icges(5,5) * t248 - Icges(5,6) * t245);
t257 = rSges(4,2) * t249 + t246 * t350;
t189 = -rSges(3,2) * t249 + t246 * rSges(3,3) + t361;
t188 = t239 + t267;
t187 = t201 * t358;
t186 = -rSges(5,3) * t246 + t364;
t183 = t249 * rSges(5,3) + t246 * t316;
t170 = t246 * rSges(4,2) + rSges(4,3) * t249 + t340;
t169 = t239 + t257;
t163 = t236 + (t427 * t249 + (-rSges(3,3) - qJ(2)) * t246) * qJD(1);
t162 = qJD(1) * t267 + t363;
t161 = pkin(8) * t335 - t226 + (t245 * t357 + t332) * pkin(4);
t160 = -qJD(1) * t348 + t342;
t159 = t185 * t249 + t208;
t158 = t185 * t246 + t207;
t155 = -qJD(3) * t246 + t236 + ((-rSges(4,2) - qJ(2)) * t246 + t350 * t249) * qJD(1);
t154 = qJD(1) * t257 + t341;
t143 = t246 * t258 + t359;
t142 = t249 * t258 - t439;
t106 = t245 * t139 + t185 * t373;
t105 = -t137 * t245 - t185 * t376;
t102 = -t246 * pkin(7) + t373 * t425 + t230 + t340 + t367;
t101 = t246 * t264 + t229 + t239 + t314 - t420;
t100 = -t246 * t174 + t249 * t275;
t99 = -t246 * t173 + t262;
t96 = t174 * t249 + t261;
t95 = t173 * t249 + t246 * t276;
t90 = t279 * t248;
t89 = -t266 + t366;
t88 = t394 + t199 + (-t185 - t216) * t358;
t83 = rSges(6,3) * t256 + t315;
t65 = (-t139 - t201) * t249 + t369 * t246;
t48 = t226 + (-qJD(3) + (t245 * t425 - t422) * qJD(4)) * t246 + (t249 * t264 - t237) * qJD(1) - t315 + t362;
t47 = (t246 * t319 - t420) * qJD(1) + t81 + t341 + t342;
t46 = t434 * t248;
t45 = (-t201 - t370) * t249 + t343 * t246;
t42 = (t185 * t355 - t83) * t245 + (-qJD(4) * t137 + t266) * t248;
t41 = (-t185 * t353 + t81) * t245 + (qJD(4) * t139 - t185 * t358 + t394) * t248;
t30 = t279 * t356 + (t246 * t81 - t249 * t83 + (t246 * t137 + t139 * t249) * qJD(1)) * t248;
t29 = t187 + (qJD(1) * t139 - t161 - t83) * t246 + (qJD(1) * t369 - t160 - t81) * t249;
t18 = t117 * t129 + t118 * t133 + t125 * t256 + t194 * t74 + t195 * t78 - t376 * t70;
t17 = t117 * t128 + t118 * t132 + t124 * t256 + t194 * t75 + t195 * t79 - t376 * t71;
t16 = t117 * t127 + t118 * t131 + t123 * t256 + t194 * t72 + t195 * t76 - t376 * t68;
t15 = t117 * t126 + t118 * t130 + t122 * t256 + t194 * t73 + t195 * t77 - t376 * t69;
t14 = t115 * t129 + t116 * t133 + t125 * t255 + t196 * t74 + t197 * t78 - t373 * t70;
t13 = t115 * t128 + t116 * t132 + t124 * t255 + t196 * t75 + t197 * t79 - t373 * t71;
t12 = t115 * t127 + t116 * t131 + t123 * t255 + t196 * t72 + t197 * t76 - t373 * t68;
t11 = t115 * t126 + t116 * t130 + t122 * t255 + t196 * t73 + t197 * t77 - t373 * t69;
t10 = t187 + (qJD(1) * t370 - t161 + t409) * t246 + (qJD(1) * t343 - t160 - t410) * t249;
t9 = -t434 * t356 + (t409 * t249 + t410 * t246 + (t246 * t371 + t249 * t370) * qJD(1)) * t248;
t8 = -qJD(1) * t303 + t17 * t249 - t18 * t246;
t7 = -qJD(1) * t305 + t15 * t249 - t16 * t246;
t6 = -qJD(1) * t299 + t13 * t249 - t14 * t246;
t5 = -qJD(1) * t301 + t11 * t249 - t12 * t246;
t4 = (qJD(4) * t303 + t34) * t245 + (qJD(1) * t304 + qJD(4) * t85 - t17 * t246 - t18 * t249) * t248;
t3 = (qJD(4) * t305 + t33) * t245 + (qJD(1) * t306 + qJD(4) * t84 - t15 * t246 - t16 * t249) * t248;
t2 = (qJD(4) * t299 + t32) * t245 + (qJD(1) * t300 + qJD(4) * t87 - t13 * t246 - t14 * t249) * t248;
t1 = (qJD(4) * t301 + t31) * t245 + (qJD(1) * t302 + qJD(4) * t86 - t11 * t246 - t12 * t249) * t248;
t25 = [-t294 * t354 + (t43 * t92 + t44 * t91) * t431 + (t101 * t48 + t102 * t47) * t432 + (t103 * t157 + t104 * t156) * t433 + 0.2e1 * m(4) * (t154 * t170 + t155 * t169) + 0.2e1 * m(3) * (t162 * t189 + t163 * t188) + t454 * t330 + t352 * t450 + (-t259 + t451) * t248 + t453 + (-Icges(5,1) * t248 + t291 + t404) * t356; m(6) * (t246 * t48 - t249 * t47 + (t101 * t249 + t102 * t246) * qJD(1)) + m(7) * (qJD(1) * t296 + t246 * t44 - t249 * t43) + m(5) * (qJD(1) * t440 - t103 * t249 + t246 * t104) + m(3) * (-t162 * t249 + t246 * t163 + (t188 * t249 + t189 * t246) * qJD(1)) + m(4) * (-t154 * t249 + t246 * t155 + (t169 * t249 + t170 * t246) * qJD(1)); 0; m(6) * (t246 * t47 + t249 * t48 + (-t101 * t246 + t102 * t249) * qJD(1)) - m(7) * t436 + m(5) * ((-t156 * t246 + t157 * t249) * qJD(1) + t441) + m(4) * (t246 * t154 + t155 * t249 + (-t169 * t246 + t170 * t249) * qJD(1)); 0; 0; (-(qJD(1) * t178 + t246 * t259) * t245 / 0.2e1 + t182 * t457 + (-t390 / 0.2e1 - t387 / 0.2e1) * qJD(4) - t271) * t249 + ((-qJD(1) * t177 + t249 * t259) * t430 + t181 * t457 + (t388 / 0.2e1 + t385 / 0.2e1) * qJD(4) + t270) * t246 + m(5) * (-t209 * t440 + t215 * t441) + m(7) * (t107 * t43 + t108 * t44 + t49 * t91 + t50 * t92) + m(6) * (t101 * t88 + t102 * t89 + t158 * t47 + t159 * t48) - (t242 / 0.2e1 + t241 / 0.2e1) * t288 * qJD(4) + ((t389 / 0.2e1 - t384 / 0.2e1 + t157 * t424 - t268) * t249 + (-t156 * t424 + t391 / 0.2e1 - t386 / 0.2e1 - t269) * t246) * qJD(1); m(6) * (t88 * t246 - t249 * t89 + (t158 * t246 + t159 * t249) * qJD(1)) + m(7) * (qJD(1) * t284 + t49 * t246 - t249 * t50); m(6) * (t89 * t246 + t249 * t88 + (t158 * t249 - t159 * t246) * qJD(1)) + m(7) * t437 - m(5) * t326; (t45 * t10 + t107 * t50 + t108 * t49) * t431 - t246 * t6 - t246 * t5 + t249 * t8 + (t158 * t89 + t159 * t88 + t65 * t29) * t432 + t249 * t7 + (-t215 * t326 + (-t249 * t220 + (-t215 * t241 + t242 * t408) * qJD(4) + (rSges(5,3) * t360 - t249 * t183 + t246 * t186) * qJD(1)) * (-t246 * t183 - t186 * t249)) * t433 + t249 * ((t249 * t143 + (-t96 + t262) * qJD(1)) * t249 + (-t95 * qJD(1) + (t178 * t356 - t182 * t354 + t359) * t246 + (-t142 + (t386 - t391) * qJD(4) + (-t173 - t275) * qJD(1)) * t249) * t246) - t246 * ((t246 * t142 + (-t99 + t261) * qJD(1)) * t246 + (-t100 * qJD(1) + (-t177 * t356 + t181 * t354 - t439) * t249 + (-t143 + (-t384 + t389) * qJD(4) + (t174 - t276) * qJD(1)) * t246) * t249) + (t96 * t246 - t249 * t95 + t304 + t306) * t358 + (t100 * t246 - t249 * t99 + t300 + t302) * t357; m(7) * (t23 * t92 + t24 * t91 + t60 * t43 + t59 * t44) + m(6) * (t101 * t42 + t102 * t41 + t105 * t48 + t106 * t47) + (t246 * t269 + t249 * t268) * t356 + (t270 * t249 + t271 * t246 + (t246 * t268 - t249 * t269) * qJD(1)) * t248 + t327; m(6) * (t42 * t246 - t249 * t41 + (t105 * t249 + t106 * t246) * qJD(1)) + m(7) * (qJD(1) * t298 - t23 * t249 + t24 * t246); m(6) * (t41 * t246 + t249 * t42 + (-t105 * t246 + t106 * t249) * qJD(1)) + m(7) * t435; m(7) * (t46 * t10 + t107 * t23 + t108 * t24 + t9 * t45 + t59 * t49 + t60 * t50) + m(6) * (t105 * t88 + t106 * t89 + t158 * t41 + t159 * t42 - t29 * t90 + t30 * t65) + (t4 / 0.2e1 + t3 / 0.2e1 + t344 * t356) * t249 + (-t2 / 0.2e1 - t1 / 0.2e1 + t345 * t356) * t246 + ((t246 * t344 - t249 * t345) * qJD(1) - (t7 + t8) * t246 / 0.2e1 - (t5 + t6) * t249 / 0.2e1 - t448 * qJD(4) / 0.2e1) * t248 + (t449 * qJD(1) + t446 * t246 + t447 * t249) * t430 - (t445 * t246 + t413 * t249) * qJD(1) / 0.2e1; (t23 * t60 + t24 * t59 + t46 * t9) * t431 + (t105 * t42 + t106 * t41 - t30 * t90) * t432 + (((t245 * t411 + t413) * t249 + (t245 * t412 + t445) * t246) * qJD(4) + t327) * t245 + ((-t1 - t2) * t249 + (-t3 - t4) * t246 + (-t447 * t246 + t446 * t249) * t245 + (t455 * t245 + t449 * t248) * qJD(4) + (t448 * t245 + t413 * t246 - t249 * t445) * qJD(1)) * t248; m(7) * (t248 * t436 + t296 * t356); 0; t360 * t347; m(7) * ((qJD(4) * t284 + t10) * t245 + (qJD(4) * t45 - t437) * t248); m(7) * ((qJD(4) * t298 + t9) * t245 + (qJD(4) * t46 - t435) * t248); 0.2e1 * (0.1e1 - t360) * t248 * t347;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;

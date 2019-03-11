% Calculate time derivative of joint inertia matrix for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:20
% EndTime: 2019-03-09 01:55:37
% DurationCPUTime: 11.96s
% Computational Cost: add. (11147->633), mult. (16378->903), div. (0->0), fcn. (14975->8), ass. (0->316)
t212 = pkin(9) + qJ(4);
t201 = sin(t212);
t202 = cos(t212);
t356 = Icges(6,6) * t202;
t364 = Icges(5,4) * t202;
t410 = t356 + t364 + (-Icges(5,2) - Icges(6,3)) * t201;
t357 = Icges(6,6) * t201;
t365 = Icges(5,4) * t201;
t409 = -t357 - t365 + (Icges(5,1) + Icges(6,2)) * t202;
t408 = t410 * qJD(4);
t407 = t409 * qJD(4);
t219 = sin(qJ(1));
t221 = cos(qJ(1));
t261 = Icges(5,2) * t202 + t365;
t125 = Icges(5,6) * t221 + t219 * t261;
t251 = Icges(6,3) * t202 + t357;
t130 = Icges(6,5) * t221 - t219 * t251;
t406 = -t125 + t130;
t264 = Icges(5,1) * t201 + t364;
t127 = Icges(5,5) * t221 + t219 * t264;
t253 = Icges(6,2) * t201 + t356;
t132 = Icges(6,4) * t221 - t219 * t253;
t405 = t127 - t132;
t131 = Icges(6,4) * t219 + t221 * t253;
t393 = -Icges(5,5) * t219 + t221 * t264;
t404 = t393 + t131;
t129 = Icges(6,5) * t219 + t221 * t251;
t394 = -Icges(5,6) * t219 + t221 * t261;
t403 = t394 + t129;
t329 = qJD(4) * t219;
t303 = t202 * t329;
t332 = qJD(1) * t221;
t228 = t201 * t332 + t303;
t402 = -qJD(4) / 0.2e1;
t218 = sin(qJ(6));
t220 = cos(qJ(6));
t294 = -qJD(1) * t202 - qJD(6);
t242 = t294 * t221;
t295 = qJD(6) * t202 + qJD(1);
t328 = qJD(4) * t220;
t79 = t220 * t242 + (t201 * t328 + t218 * t295) * t219;
t305 = t201 * t329;
t344 = t219 * t220;
t80 = -t295 * t344 + (t242 + t305) * t218;
t42 = t80 * rSges(7,1) + t79 * rSges(7,2) + rSges(7,3) * t228;
t401 = -pkin(8) * t228 - t42;
t244 = t130 * t202 + t132 * t201;
t400 = t221 * t244;
t247 = t125 * t202 + t127 * t201;
t399 = t221 * t247;
t288 = rSges(5,1) * t201 + rSges(5,2) * t202;
t237 = t221 * t288;
t346 = t218 * t221;
t152 = -t202 * t344 - t346;
t343 = t220 * t221;
t345 = t219 * t218;
t153 = -t202 * t345 + t343;
t350 = t201 * t219;
t105 = t153 * rSges(7,1) + t152 * rSges(7,2) + rSges(7,3) * t350;
t210 = t221 * pkin(5);
t398 = -pkin(8) * t350 - t105 - t210;
t215 = sin(pkin(9));
t377 = pkin(3) * t215;
t198 = t221 * t377;
t207 = t221 * qJ(2);
t217 = -pkin(7) - qJ(3);
t310 = t219 * t217 + t198 + t207;
t380 = -rSges(5,3) - pkin(1);
t106 = t219 * t380 + t237 + t310;
t208 = t221 * rSges(5,3);
t348 = t202 * t219;
t136 = rSges(5,1) * t350 + rSges(5,2) * t348 + t208;
t337 = t221 * pkin(1) + t219 * qJ(2);
t240 = -t217 * t221 + t219 * t377 + t337;
t107 = t240 + t136;
t397 = -t106 * t219 + t107 * t221;
t338 = qJ(2) * t332 + qJD(2) * t219;
t309 = qJD(3) * t221 + t338;
t333 = qJD(1) * t219;
t266 = qJD(1) * t198 + t217 * t333 + t309;
t306 = t202 * t332;
t311 = rSges(5,1) * t228 + rSges(5,2) * t306;
t331 = qJD(4) * t201;
t61 = (-rSges(5,2) * t331 + qJD(1) * t380) * t219 + t266 + t311;
t375 = rSges(5,2) * t201;
t172 = rSges(5,1) * t202 - t375;
t205 = qJD(2) * t221;
t296 = -qJD(3) * t219 + t205;
t291 = t217 * t332 + t296;
t299 = -qJ(2) - t377;
t327 = qJD(4) * t221;
t62 = t172 * t327 + (t380 * t221 + (-t288 + t299) * t219) * qJD(1) + t291;
t396 = t219 * t62 - t221 * t61;
t256 = Icges(5,5) * t201 + Icges(5,6) * t202;
t395 = -Icges(5,3) * t219 + t221 * t256;
t259 = Icges(6,4) * t201 + Icges(6,5) * t202;
t133 = Icges(6,1) * t219 + t221 * t259;
t150 = t202 * t343 - t345;
t151 = t202 * t346 + t344;
t287 = -t151 * rSges(7,1) - t150 * rSges(7,2);
t349 = t201 * t221;
t104 = -rSges(7,3) * t349 - t287;
t286 = rSges(7,1) * t218 + rSges(7,2) * t220;
t122 = rSges(7,3) * t202 + t201 * t286;
t302 = t202 * t327;
t308 = t201 * t333;
t227 = -t302 + t308;
t304 = t201 * t327;
t224 = t219 * t294 - t304;
t81 = t220 * t224 - t295 * t346;
t82 = t218 * t224 + t295 * t343;
t292 = t82 * rSges(7,1) + t81 * rSges(7,2);
t43 = rSges(7,3) * t227 + t292;
t324 = qJD(6) * t201;
t83 = (rSges(7,1) * t220 - rSges(7,2) * t218) * t324 + (-rSges(7,3) * t201 + t202 * t286) * qJD(4);
t21 = (-t122 * t327 - t43) * t202 + (qJD(4) * t104 + t122 * t333 - t221 * t83) * t201;
t238 = t122 * t332 + t219 * t83;
t22 = (-t122 * t329 + t42) * t202 + (-qJD(4) * t105 - t238) * t201;
t59 = t105 * t202 - t122 * t350;
t60 = -t202 * t104 - t122 * t349;
t392 = qJD(1) * (t219 * t59 + t221 * t60) + t21 * t219 - t22 * t221;
t391 = 2 * m(5);
t390 = 2 * m(6);
t389 = 2 * m(7);
t213 = t219 ^ 2;
t214 = t221 ^ 2;
t388 = m(6) / 0.2e1;
t387 = m(7) / 0.2e1;
t386 = -pkin(1) - pkin(5);
t385 = t202 / 0.2e1;
t384 = t219 / 0.2e1;
t383 = t221 / 0.2e1;
t382 = -rSges(6,1) - pkin(1);
t381 = rSges(3,2) - pkin(1);
t379 = rSges(7,3) + pkin(8);
t378 = m(5) * t172;
t376 = t219 * pkin(5);
t374 = rSges(6,2) * t201;
t373 = rSges(6,2) * t202;
t372 = rSges(6,3) * t202;
t371 = t219 * rSges(6,1);
t370 = t219 * rSges(5,3);
t209 = t221 * rSges(6,1);
t367 = rSges(4,3) + qJ(3);
t362 = Icges(7,4) * t218;
t361 = Icges(7,4) * t220;
t354 = qJ(5) * t202;
t263 = Icges(7,1) * t218 + t361;
t121 = Icges(7,5) * t202 + t201 * t263;
t351 = t121 * t218;
t347 = t202 * t221;
t322 = pkin(8) * t349;
t342 = t104 - t322 + t376;
t135 = qJD(5) * t201 + (-pkin(4) * t201 + t354) * qJD(4);
t170 = pkin(4) * t202 + qJ(5) * t201;
t341 = t219 * t135 + t170 * t332;
t192 = pkin(4) * t350;
t315 = qJ(5) * t348;
t142 = t192 - t315;
t284 = t372 + t374;
t340 = t219 * t284 - t142 - t209;
t336 = t213 + t214;
t123 = Icges(5,3) * t221 + t219 * t256;
t335 = qJD(1) * t123;
t134 = Icges(6,1) * t221 - t219 * t259;
t334 = qJD(1) * t134;
t330 = qJD(4) * t202;
t326 = qJD(5) * t202;
t325 = qJD(5) * t219;
t323 = -pkin(1) - t367;
t255 = Icges(7,5) * t218 + Icges(7,6) * t220;
t119 = Icges(7,3) * t202 + t201 * t255;
t258 = Icges(7,2) * t220 + t362;
t120 = Icges(7,6) * t202 + t201 * t258;
t76 = (Icges(7,5) * t220 - Icges(7,6) * t218) * t324 + (-Icges(7,3) * t201 + t202 * t255) * qJD(4);
t77 = (-Icges(7,2) * t218 + t361) * t324 + (-Icges(7,6) * t201 + t202 * t258) * qJD(4);
t78 = (Icges(7,1) * t220 - t362) * t324 + (-Icges(7,5) * t201 + t202 * t263) * qJD(4);
t16 = t119 * t227 + t81 * t120 + t82 * t121 + t150 * t77 + t151 * t78 - t349 * t76;
t100 = Icges(7,4) * t151 + Icges(7,2) * t150 - Icges(7,6) * t349;
t102 = Icges(7,1) * t151 + Icges(7,4) * t150 - Icges(7,5) * t349;
t250 = t100 * t220 + t102 * t218;
t37 = Icges(7,5) * t82 + Icges(7,6) * t81 + Icges(7,3) * t227;
t39 = Icges(7,4) * t82 + Icges(7,2) * t81 + Icges(7,6) * t227;
t41 = Icges(7,1) * t82 + Icges(7,4) * t81 + Icges(7,5) * t227;
t98 = Icges(7,5) * t151 + Icges(7,6) * t150 - Icges(7,3) * t349;
t9 = (qJD(4) * t250 + t37) * t202 + (-qJD(4) * t98 + t218 * t41 + t220 * t39 + (-t100 * t218 + t102 * t220) * qJD(6)) * t201;
t321 = -t9 / 0.2e1 - t16 / 0.2e1;
t101 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t350;
t103 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t350;
t249 = t101 * t220 + t103 * t218;
t36 = Icges(7,5) * t80 + Icges(7,6) * t79 + Icges(7,3) * t228;
t38 = Icges(7,4) * t80 + Icges(7,2) * t79 + Icges(7,6) * t228;
t40 = Icges(7,1) * t80 + Icges(7,4) * t79 + Icges(7,5) * t228;
t99 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t350;
t10 = (qJD(4) * t249 + t36) * t202 + (-qJD(4) * t99 + t218 * t40 + t220 * t38 + (-t101 * t218 + t103 * t220) * qJD(6)) * t201;
t15 = t119 * t228 + t79 * t120 + t80 * t121 + t152 * t77 + t153 * t78 + t350 * t76;
t320 = t10 / 0.2e1 + t15 / 0.2e1;
t30 = t201 * t250 + t202 * t98;
t44 = -t119 * t349 + t150 * t120 + t151 * t121;
t319 = t30 / 0.2e1 + t44 / 0.2e1;
t31 = t201 * t249 + t202 * t99;
t45 = t119 * t350 + t120 * t152 + t121 * t153;
t318 = t31 / 0.2e1 + t45 / 0.2e1;
t317 = t386 * t219;
t316 = qJ(5) * t347;
t314 = -t142 + t398;
t313 = pkin(4) * t228 + qJ(5) * t305;
t312 = pkin(4) * t302 + qJ(5) * t304 + qJD(1) * t315;
t301 = qJD(6) * t120 * t218;
t300 = t256 * t402 + t259 * qJD(4) / 0.2e1;
t298 = pkin(8) * t202 + t122;
t162 = t288 * qJD(4);
t297 = t162 * t336;
t193 = pkin(4) * t349;
t293 = t193 + t310;
t290 = t298 * t219;
t289 = rSges(4,1) * t215 + rSges(4,2) * cos(pkin(9));
t285 = -rSges(6,3) * t201 + t373;
t225 = t266 + t313;
t23 = t225 - t202 * t325 + (-t316 + t317) * qJD(1) - t401;
t226 = t291 + t312;
t24 = (qJD(4) * t379 - qJD(5)) * t347 + (t386 * t221 + ((-pkin(4) - t379) * t201 + t299) * t219) * qJD(1) + t226 - t292;
t278 = t219 * t24 - t221 * t23;
t25 = t150 * t100 + t151 * t102 - t349 * t98;
t26 = t150 * t101 + t151 * t103 - t349 * t99;
t277 = t219 * t26 - t221 * t25;
t17 = t25 * t219 + t221 * t26;
t27 = t100 * t152 + t102 * t153 + t350 * t98;
t28 = t101 * t152 + t103 * t153 + t350 * t99;
t276 = t219 * t28 - t221 * t27;
t18 = t27 * t219 + t221 * t28;
t275 = t219 * t31 - t221 * t30;
t274 = t30 * t219 + t221 * t31;
t146 = t170 * t333;
t32 = t146 + qJD(1) * t290 + (pkin(8) * t331 - t135 - t83) * t221;
t33 = (-t305 + t306) * pkin(8) + t238 + t341;
t273 = t33 * t219 - t221 * t32;
t177 = rSges(6,3) * t305;
t236 = -t374 + (-rSges(6,3) - qJ(5)) * t202;
t222 = t219 * t382 + t221 * t236;
t34 = t177 + (-rSges(6,2) * qJD(4) - qJD(5)) * t348 + t222 * qJD(1) + t225;
t35 = (-qJD(4) * t285 - t326) * t221 + (t382 * t221 + (t372 + (rSges(6,2) - pkin(4)) * t201 + t299) * t219) * qJD(1) + t226;
t272 = t219 * t35 - t221 * t34;
t47 = t317 + (t201 * t379 - t354) * t221 + t287 + t293;
t235 = t192 + t240;
t48 = t235 - t315 - t398;
t271 = t219 * t47 - t221 * t48;
t270 = t219 * t60 - t221 * t59;
t161 = t284 * qJD(4);
t63 = -t285 * t333 + t146 + (-t135 - t161) * t221;
t64 = t219 * t161 - t285 * t332 + t341;
t268 = t64 * t219 - t221 * t63;
t65 = t222 + t293;
t66 = t219 * t236 + t209 + t235;
t267 = t219 * t65 - t221 * t66;
t260 = Icges(6,4) * t202 - Icges(6,5) * t201;
t257 = Icges(5,5) * t202 - Icges(5,6) * t201;
t248 = t104 * t219 + t105 * t221;
t246 = -t201 * t393 - t202 * t394;
t245 = t129 * t202 + t131 * t201;
t243 = (t388 + t387) * t331;
t241 = t220 * t121 * t324 + t330 * t351 + (t120 * t328 + t76) * t202 + (t218 * t78 + t220 * t77) * t201;
t239 = rSges(3,3) * t221 + t219 * t381;
t234 = t246 * t219;
t233 = t245 * t219;
t223 = t219 * t323 + t221 * t289;
t154 = t219 * t170;
t149 = -rSges(3,2) * t221 + t219 * rSges(3,3) + t337;
t148 = t207 + t239;
t143 = -t193 + t316;
t140 = t221 * t143;
t138 = t221 * t284 + t371;
t137 = t370 - t237;
t117 = t205 + (t381 * t221 + (-rSges(3,3) - qJ(2)) * t219) * qJD(1);
t116 = qJD(1) * t239 + t338;
t115 = t219 * t289 + t221 * t367 + t337;
t114 = t207 + t223;
t113 = (-t170 + t285) * t221;
t112 = -t219 * t285 + t154;
t97 = t260 * t327 + t334;
t96 = -qJD(1) * t133 - t260 * t329;
t87 = t395 * qJD(1) + t257 * t329;
t86 = -t257 * t327 + t335;
t85 = (t323 * t221 + (-qJ(2) - t289) * t219) * qJD(1) + t296;
t84 = qJD(1) * t223 + t309;
t75 = (-qJ(5) * t332 - t325) * t202 + t313;
t72 = t221 * (pkin(4) * t308 + t221 * t326 - t312);
t71 = (-t170 - t298) * t221;
t70 = t154 + t290;
t58 = t134 * t221 - t219 * t244;
t57 = t133 * t221 - t233;
t56 = t219 * t134 + t400;
t55 = t219 * t133 + t221 * t245;
t54 = -t219 * t395 - t221 * t246;
t53 = t219 * t123 - t399;
t52 = -t221 * t395 + t234;
t51 = t123 * t221 + t219 * t247;
t50 = t138 * t221 + t219 * t340 + t140;
t49 = t119 * t202 + (t120 * t220 + t351) * t201;
t46 = t248 * t201;
t29 = t219 * t314 + t221 * t342 + t140;
t20 = t72 + (-t75 - t177) * t219 + (t213 * t373 + t214 * t285) * qJD(4) + ((t340 + t209) * t221 + (-t138 - t143 + t371) * t219) * qJD(1);
t19 = ((-qJD(4) * t119 - t301) * t201 + t241) * t202;
t14 = t248 * t330 + (t219 * t43 + t221 * t42 + (t104 * t221 - t219 * t105) * qJD(1)) * t201;
t13 = t201 * t276 + t45 * t202;
t12 = t201 * t277 + t44 * t202;
t11 = t72 + (-pkin(8) * t302 + t43) * t221 + (-t75 + t401) * t219 + ((t314 + t210) * t221 + (-t143 + t322 - t342 + t376) * t219) * qJD(1);
t8 = -t99 * t302 + t81 * t101 + t82 * t103 + t150 * t38 + t151 * t40 + (-t221 * t36 + t333 * t99) * t201;
t7 = -t98 * t302 + t81 * t100 + t82 * t102 + t150 * t39 + t151 * t41 + (-t221 * t37 + t333 * t98) * t201;
t6 = t99 * t303 + t79 * t101 + t80 * t103 + t152 * t38 + t153 * t40 + (t219 * t36 + t332 * t99) * t201;
t5 = t98 * t303 + t79 * t100 + t80 * t102 + t152 * t39 + t153 * t41 + (t219 * t37 + t332 * t98) * t201;
t4 = -qJD(1) * t277 + t7 * t219 + t221 * t8;
t3 = -qJD(1) * t276 + t5 * t219 + t221 * t6;
t2 = (qJD(4) * t277 + t16) * t202 + (qJD(1) * t17 - qJD(4) * t44 + t219 * t8 - t221 * t7) * t201;
t1 = (qJD(4) * t276 + t15) * t202 + (qJD(1) * t18 - qJD(4) * t45 + t219 * t6 - t221 * t5) * t201;
t67 = [t241 + (t23 * t48 + t24 * t47) * t389 + (t34 * t66 + t35 * t65) * t390 + (t106 * t62 + t107 * t61) * t391 + 0.2e1 * m(4) * (t114 * t85 + t115 * t84) + 0.2e1 * m(3) * (t116 * t149 + t117 * t148) + (-t264 - t253) * t330 - t408 * t202 + (t261 + t251 - t119) * t331 + (-t301 - t407) * t201; m(7) * ((t219 * t48 + t221 * t47) * qJD(1) + t278) + m(6) * ((t219 * t66 + t221 * t65) * qJD(1) + t272) + m(5) * ((t106 * t221 + t107 * t219) * qJD(1) + t396) + m(4) * (t219 * t85 - t221 * t84 + (t114 * t221 + t115 * t219) * qJD(1)) + m(3) * (-t116 * t221 + t219 * t117 + (t148 * t221 + t149 * t219) * qJD(1)); 0; m(7) * (-qJD(1) * t271 + t219 * t23 + t221 * t24) + m(6) * (-qJD(1) * t267 + t219 * t34 + t221 * t35) + m(5) * (t397 * qJD(1) + t219 * t61 + t221 * t62) + m(4) * (t219 * t84 + t221 * t85 + (-t114 * t219 + t115 * t221) * qJD(1)); 0; 0; (t300 * t221 + t320) * t221 + (t300 * t219 - t321) * t219 + m(5) * (t397 * t162 + t396 * t172) + m(6) * (t112 * t35 + t113 * t34 + t63 * t66 + t64 * t65) + m(7) * (t23 * t71 + t24 * t70 + t32 * t48 + t33 * t47) + ((t403 * qJD(4) - t327 * t409) * t384 + (t406 * qJD(4) + t407 * t219) * t383 + (t404 * t383 + t405 * t384) * qJD(1)) * t202 + ((t404 * qJD(4) + t327 * t410) * t384 + (-t405 * qJD(4) - t408 * t219) * t383 + (-t403 * t383 + t406 * t384) * qJD(1)) * t201 + ((t107 * t378 + (-t127 / 0.2e1 + t132 / 0.2e1) * t202 + (t125 / 0.2e1 - t130 / 0.2e1) * t201 - t318) * t219 + (t106 * t378 + (-t393 / 0.2e1 - t131 / 0.2e1) * t202 + (t394 / 0.2e1 + t129 / 0.2e1) * t201 + t319) * t221) * qJD(1); m(6) * ((t112 * t221 + t113 * t219) * qJD(1) + t268) + m(7) * ((t219 * t71 + t221 * t70) * qJD(1) + t273) - m(5) * t297; m(6) * (t63 * t219 + t221 * t64 + (-t112 * t219 + t113 * t221) * qJD(1)) + m(7) * (t32 * t219 + t221 * t33 + (-t219 * t70 + t221 * t71) * qJD(1)); (t11 * t29 + t32 * t71 + t33 * t70) * t389 + t219 * t4 + t221 * t3 + (t112 * t64 + t113 * t63 + t20 * t50) * t390 + ((-t219 * t136 + t137 * t221) * (-t219 * t311 + (-t172 * t214 + t213 * t375) * qJD(4) + ((-t136 + t208) * t221 + (-t137 + t237 + t370) * t219) * qJD(1)) - t172 * t297) * t391 + t219 * ((t219 * t97 + (-t56 - t233) * qJD(1)) * t219 + (t55 * qJD(1) + (-t130 * t331 + t132 * t330 + t334) * t221 + (t96 + (-t129 * t201 + t131 * t202) * qJD(4) - t244 * qJD(1)) * t219) * t221) + t219 * ((t219 * t86 + (-t53 + t234) * qJD(1)) * t219 + (t54 * qJD(1) + (t125 * t331 - t127 * t330 + t335) * t221 + (t87 + (-t201 * t394 + t202 * t393) * qJD(4) + t247 * qJD(1)) * t219) * t221) + t221 * ((t221 * t96 + (t57 - t400) * qJD(1)) * t221 + (-t58 * qJD(1) + (t129 * t331 - t131 * t330) * t219 + (t97 + (t130 * t201 - t132 * t202) * qJD(4) + (-t134 - t245) * qJD(1)) * t221) * t219) + t221 * ((t221 * t87 + (t52 + t399) * qJD(1)) * t221 + (-t51 * qJD(1) + (-t330 * t393 + t331 * t394) * t219 + (t86 + (-t125 * t201 + t127 * t202) * qJD(4) + (-t123 + t246) * qJD(1)) * t221) * t219) + (-t18 + (-t51 - t58) * t221 + (-t52 - t57) * t219) * t333 + (t17 + (t53 + t56) * t221 + (t54 + t55) * t219) * t332; 0.2e1 * (t267 * t388 + t271 * t387) * t331 + 0.2e1 * ((-t332 * t47 - t333 * t48 - t278) * t387 + (-t332 * t65 - t333 * t66 - t272) * t388) * t202; 0.2e1 * t336 * t243; 0; 0.2e1 * ((-t327 * t71 + t329 * t70 + t11) * t387 + (t112 * t329 - t113 * t327 + t20) * t388) * t201 + 0.2e1 * ((qJD(4) * t29 - t332 * t70 - t333 * t71 - t273) * t387 + (qJD(4) * t50 - t112 * t332 - t113 * t333 - t268) * t388) * t202; 0.4e1 * (0.1e1 - t336) * t202 * t243; m(7) * (t21 * t47 + t22 * t48 + t23 * t59 + t24 * t60) + t19 + (t219 * t318 - t221 * t319) * t330 + (-t49 * qJD(4) + t321 * t221 + t320 * t219 + (t219 * t319 + t221 * t318) * qJD(1)) * t201; m(7) * t392; m(7) * (-qJD(1) * t270 + t21 * t221 + t22 * t219); m(7) * (t11 * t46 + t14 * t29 + t21 * t70 + t22 * t71 + t32 * t59 + t33 * t60) + (qJD(1) * t12 / 0.2e1 - t17 * t330 / 0.2e1 + (qJD(1) * t30 + t10) * t385 + t1 / 0.2e1) * t221 + (t2 / 0.2e1 + t18 * t330 / 0.2e1 + (-qJD(1) * t31 + t9) * t385 - qJD(1) * t13 / 0.2e1) * t219 + (t3 * t384 - t221 * t4 / 0.2e1 + t274 * t402 + (t17 * t384 + t18 * t383) * qJD(1)) * t201; m(7) * ((qJD(4) * t270 + t14) * t201 + (qJD(4) * t46 - t392) * t202); (t14 * t46 + t21 * t60 + t22 * t59) * t389 + (t19 + (-t221 * t12 + t219 * t13 + t202 * t275) * qJD(4)) * t202 + (t219 * t1 - t221 * t2 + t202 * (t10 * t219 - t9 * t221) + (-t201 * t275 - 0.2e1 * t49 * t202) * qJD(4) + (t219 * t12 + t221 * t13 + t202 * t274) * qJD(1)) * t201;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t67(1) t67(2) t67(4) t67(7) t67(11) t67(16); t67(2) t67(3) t67(5) t67(8) t67(12) t67(17); t67(4) t67(5) t67(6) t67(9) t67(13) t67(18); t67(7) t67(8) t67(9) t67(10) t67(14) t67(19); t67(11) t67(12) t67(13) t67(14) t67(15) t67(20); t67(16) t67(17) t67(18) t67(19) t67(20) t67(21);];
Mq  = res;

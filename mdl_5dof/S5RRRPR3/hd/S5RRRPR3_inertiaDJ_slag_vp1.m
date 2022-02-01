% Calculate time derivative of joint inertia matrix for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:28
% EndTime: 2022-01-20 11:42:44
% DurationCPUTime: 8.47s
% Computational Cost: add. (10947->496), mult. (7952->683), div. (0->0), fcn. (5984->10), ass. (0->297)
t235 = qJ(3) + pkin(9);
t226 = cos(t235);
t225 = sin(t235);
t366 = Icges(5,4) * t225;
t169 = Icges(5,2) * t226 + t366;
t240 = cos(qJ(3));
t238 = sin(qJ(3));
t368 = Icges(4,4) * t238;
t197 = Icges(4,2) * t240 + t368;
t413 = -t169 * t225 - t197 * t238;
t365 = Icges(5,4) * t226;
t170 = Icges(5,1) * t225 + t365;
t367 = Icges(4,4) * t240;
t198 = Icges(4,1) * t238 + t367;
t412 = -t170 * t226 - t198 * t240 - t413;
t168 = Icges(5,5) * t225 + Icges(5,6) * t226;
t196 = Icges(4,5) * t238 + Icges(4,6) * t240;
t411 = t168 + t196;
t237 = -qJ(4) - pkin(7);
t232 = pkin(8) - t237;
t234 = qJD(1) + qJD(2);
t383 = pkin(4) * t225;
t384 = pkin(3) * t238;
t293 = -t383 - t384;
t410 = qJD(3) * t293 + t232 * t234;
t236 = qJ(1) + qJ(2);
t228 = sin(t236);
t229 = cos(t236);
t330 = qJD(3) * t240;
t346 = t234 * t238;
t408 = -t228 * t330 - t229 * t346;
t172 = rSges(5,1) * t225 + rSges(5,2) * t226;
t375 = rSges(4,2) * t238;
t380 = rSges(4,1) * t240;
t407 = -t375 + t380;
t282 = Icges(5,5) * t226 - Icges(5,6) * t225;
t283 = Icges(4,5) * t240 - Icges(4,6) * t238;
t406 = t412 * t234 + (t282 + t283) * qJD(3);
t285 = -Icges(5,2) * t225 + t365;
t152 = t285 * qJD(3);
t288 = Icges(5,1) * t226 - t366;
t153 = t288 * qJD(3);
t286 = -Icges(4,2) * t238 + t367;
t176 = t286 * qJD(3);
t289 = Icges(4,1) * t240 - t368;
t177 = t289 * qJD(3);
t405 = -t152 * t225 + t153 * t226 - t176 * t238 + t177 * t240 + (0.2e1 * t168 + t196) * t234 + (-0.2e1 * t169 * t226 - 0.2e1 * t170 * t225 - t197 * t240 - t198 * t238) * qJD(3);
t264 = t286 * t229;
t133 = Icges(4,6) * t228 + t264;
t267 = t289 * t229;
t135 = Icges(4,5) * t228 + t267;
t273 = t133 * t238 - t135 * t240;
t404 = t228 * t273;
t263 = t285 * t229;
t120 = Icges(5,6) * t228 + t263;
t266 = t288 * t229;
t122 = Icges(5,5) * t228 + t266;
t277 = t120 * t225 - t122 * t226;
t403 = t228 * t277;
t227 = qJ(5) + t235;
t217 = sin(t227);
t218 = cos(t227);
t363 = Icges(6,4) * t218;
t284 = -Icges(6,2) * t217 + t363;
t262 = t284 * t229;
t108 = Icges(6,6) * t228 + t262;
t364 = Icges(6,4) * t217;
t287 = Icges(6,1) * t218 - t364;
t265 = t287 * t229;
t110 = Icges(6,5) * t228 + t265;
t279 = t108 * t217 - t110 * t218;
t402 = t228 * t279;
t132 = -Icges(4,6) * t229 + t228 * t286;
t134 = -Icges(4,5) * t229 + t228 * t289;
t275 = t132 * t238 - t134 * t240;
t401 = t229 * t275;
t119 = -Icges(5,6) * t229 + t228 * t285;
t121 = -Icges(5,5) * t229 + t228 * t288;
t278 = t119 * t225 - t121 * t226;
t400 = t229 * t278;
t107 = -Icges(6,6) * t229 + t228 * t284;
t109 = -Icges(6,5) * t229 + t228 * t287;
t280 = t107 * t217 - t109 * t218;
t399 = t229 * t280;
t210 = t228 * rSges(6,3);
t377 = rSges(6,1) * t218;
t398 = t229 * t377 + t210;
t233 = qJD(3) + qJD(5);
t156 = Icges(6,2) * t218 + t364;
t157 = Icges(6,1) * t217 + t363;
t272 = t156 * t217 - t157 * t218;
t281 = Icges(6,5) * t218 - Icges(6,6) * t217;
t397 = t233 * t281 + t234 * t272;
t349 = t229 * t237;
t230 = t240 * pkin(3);
t222 = t230 + pkin(2);
t381 = pkin(2) - t222;
t394 = t228 * t381 - t349;
t393 = 2 * m(3);
t392 = 2 * m(4);
t391 = 2 * m(5);
t390 = 2 * m(6);
t389 = t228 / 0.2e1;
t388 = -t229 / 0.2e1;
t187 = t407 * qJD(3);
t387 = m(4) * t187;
t204 = rSges(4,1) * t238 + rSges(4,2) * t240;
t386 = m(4) * t204;
t239 = sin(qJ(1));
t385 = pkin(1) * t239;
t382 = pkin(4) * t226;
t219 = t228 * pkin(7);
t378 = rSges(5,1) * t226;
t376 = rSges(6,1) * t228;
t374 = rSges(5,2) * t225;
t372 = rSges(6,2) * t217;
t371 = pkin(1) * qJD(1);
t212 = t228 * rSges(4,3);
t211 = t228 * rSges(5,3);
t158 = rSges(6,1) * t217 + rSges(6,2) * t218;
t269 = -t158 + t293;
t94 = t269 * t229;
t370 = t234 * t94;
t220 = t229 * pkin(7);
t111 = t220 - t394;
t296 = -t222 * t229 + t228 * t237;
t337 = -pkin(2) * t229 - t219;
t112 = -t296 + t337;
t369 = t111 * t228 + t112 * t229;
t358 = t156 * t233;
t357 = t157 * t233;
t354 = t217 * t233;
t353 = t218 * t233;
t352 = t228 * t234;
t351 = t228 * t240;
t350 = t229 * t234;
t347 = t234 * t237;
t182 = t228 * t372;
t341 = rSges(6,3) * t229 + t182;
t113 = t218 * t376 - t341;
t326 = t229 * t372;
t114 = -t326 + t398;
t62 = t113 * t228 + t114 * t229;
t345 = t410 * t229;
t186 = t222 + t382;
t344 = t186 * t229 + t228 * t232;
t190 = t228 * t374;
t343 = rSges(5,3) * t350 + t190 * t234;
t322 = t228 * t346;
t342 = rSges(4,2) * t322 + rSges(4,3) * t350;
t331 = qJD(3) * t238;
t316 = t228 * t331;
t340 = pkin(3) * t316 + t228 * t347;
t339 = rSges(5,3) * t229 + t190;
t338 = rSges(4,3) * t229 + t228 * t375;
t336 = t228 ^ 2 + t229 ^ 2;
t335 = qJD(3) * t225;
t334 = qJD(3) * t226;
t332 = qJD(3) * t229;
t325 = rSges(6,2) * t353;
t244 = -t229 * t325 + t234 * t182 + rSges(6,3) * t350 + (-t218 * t352 - t229 * t354) * rSges(6,1);
t320 = -t228 * t325 - t234 * t326 - t354 * t376;
t329 = t228 * (t234 * t398 + t320) + t229 * t244 + t113 * t350;
t203 = pkin(7) * t350;
t208 = qJD(4) * t228;
t314 = t229 * t331;
t295 = pkin(3) * t314;
t209 = qJD(4) * t229;
t317 = -t209 - t340;
t328 = t228 * ((-t229 * t381 - t219) * t234 + t317) + t229 * (t234 * t394 - t203 + t208 - t295) + t111 * t350;
t327 = t229 * t374;
t324 = t239 * t371;
t241 = cos(qJ(1));
t323 = t241 * t371;
t319 = -qJD(3) * t172 * t228 - t234 * t327;
t318 = -rSges(4,1) * t316 + rSges(4,2) * t408;
t313 = t352 / 0.2e1;
t312 = t350 / 0.2e1;
t311 = -pkin(2) - t380;
t310 = -t172 - t384;
t22 = -t186 * t352 + t208 + t244 + t345;
t20 = t22 - t324;
t302 = -t186 - t377;
t79 = t228 * t302 + t229 * t232 + t341;
t75 = t79 - t385;
t309 = t234 * t75 - t20;
t297 = t410 * t228;
t23 = t209 + (t229 * t302 - t210) * t234 - t297 - t320;
t21 = t23 - t323;
t231 = t241 * pkin(1);
t80 = t114 + t344;
t76 = t231 + t80;
t308 = t234 * t76 + t21;
t307 = t234 * t79 - t22;
t306 = t234 * t80 + t23;
t188 = pkin(3) * t322;
t141 = (-t372 + t377) * t233;
t258 = -t141 + (-t230 - t382) * qJD(3);
t36 = t188 + (t158 + t383) * t352 + t258 * t229;
t93 = t269 * t228;
t305 = t234 * t93 + t36;
t37 = t228 * t258 + t370;
t304 = -t37 + t370;
t303 = -t222 - t378;
t174 = rSges(3,1) * t229 - rSges(3,2) * t228;
t257 = Icges(6,5) * t234 - t357;
t61 = t228 * t257 + t234 * t265;
t301 = -t107 * t233 + t61;
t60 = t229 * t257 - t287 * t352;
t300 = -t108 * t233 + t60;
t256 = Icges(6,6) * t234 - t358;
t59 = t228 * t256 + t234 * t262;
t299 = t109 * t233 + t59;
t58 = t229 * t256 - t284 * t352;
t298 = t110 * t233 + t58;
t146 = -rSges(3,1) * t350 + rSges(3,2) * t352;
t105 = -Icges(6,3) * t229 + t228 * t281;
t24 = -t105 * t229 - t228 * t280;
t259 = t281 * t229;
t106 = Icges(6,3) * t228 + t259;
t25 = -t106 * t229 - t402;
t26 = t105 * t228 - t399;
t27 = t106 * t228 - t229 * t279;
t155 = Icges(6,5) * t217 + Icges(6,6) * t218;
t255 = Icges(6,3) * t234 - t155 * t233;
t56 = t229 * t255 - t281 * t352;
t57 = t228 * t255 + t234 * t259;
t294 = -t229 * ((t229 * t57 + (t25 + t399) * t234) * t229 + (t24 * t234 + (-t108 * t353 - t110 * t354 - t217 * t58 + t218 * t60) * t228 + (-t56 + (t110 * t234 - t301) * t218 + (-t108 * t234 + t299) * t217) * t229) * t228) + t228 * ((t228 * t56 + (t26 + t402) * t234) * t228 + (t27 * t234 + (t107 * t353 + t109 * t354 + t217 * t59 - t218 * t61) * t229 + (-t57 + (t109 * t234 + t300) * t218 + (-t107 * t234 - t298) * t217) * t228) * t229) + (t228 * t25 - t229 * t24) * t352 + (t228 * t27 - t229 * t26) * t350;
t173 = -rSges(3,1) * t228 - rSges(3,2) * t229;
t290 = t303 * t228;
t276 = t132 * t240 + t134 * t238;
t274 = t133 * t240 + t135 * t238;
t137 = t229 * t407 + t212;
t124 = t229 * t378 + t211 - t327;
t139 = t284 * t233;
t140 = t287 * t233;
t247 = t155 * t234 + (t140 - t358) * t218 + (-t139 - t357) * t217;
t268 = (t217 * t300 + t218 * t298 + t228 * t397 + t229 * t247) * t389 + (t217 * t301 + t218 * t299 + t228 * t247 - t229 * t397) * t388 + (t107 * t218 + t109 * t217 - t155 * t229 - t228 * t272) * t313 + (t108 * t218 + t110 * t217 + t155 * t228 - t229 * t272) * t312;
t145 = t173 * t234;
t261 = t283 * t229;
t260 = t282 * t229;
t102 = t137 - t337;
t101 = t228 * t311 + t220 + t338;
t90 = t124 - t296;
t254 = Icges(4,5) * t234 - qJD(3) * t198;
t252 = Icges(4,6) * t234 - qJD(3) * t197;
t250 = Icges(4,3) * t234 - qJD(3) * t196;
t249 = Icges(5,3) * t234 - qJD(3) * t168;
t89 = t290 + t339 - t349;
t248 = t139 * t218 + t140 * t217 + t152 * t226 + t153 * t225 - t156 * t354 + t157 * t353 + t170 * t334 + t176 * t240 + t177 * t238 + t198 * t330;
t55 = (t311 * t229 + (-rSges(4,3) - pkin(7)) * t228) * t234 - t318;
t35 = (t229 * t303 - t211) * t234 - t317 - t319;
t243 = qJD(3) * t413 + t248;
t242 = t268 + (t238 * (t229 * t254 - t289 * t352) + t240 * (t229 * t252 - t286 * t352) + (-t225 * t288 - t226 * t285) * t352 + t405 * t229 + t406 * t228 + (-t273 - t277) * qJD(3)) * t389 + (t238 * (t228 * t254 + t234 * t267) + t240 * (t228 * t252 + t234 * t264) + (t225 * t266 + t226 * t263) * t234 - t406 * t229 + t405 * t228 + (-t275 - t278) * qJD(3)) * t388 + (t119 * t226 + t121 * t225 - t228 * t412 - t229 * t411 + t276) * t313 + (t120 * t226 + t122 * t225 + t228 * t411 - t229 * t412 + t274) * t312;
t54 = -rSges(4,2) * t229 * t330 - pkin(2) * t352 + t203 + (-t234 * t351 - t314) * rSges(4,1) + t342;
t34 = t208 + t234 * t290 + (qJD(3) * t310 - t347) * t229 + t343;
t154 = (-t374 + t378) * qJD(3);
t149 = t174 + t231;
t148 = t173 - t385;
t136 = rSges(4,1) * t351 - t338;
t131 = Icges(4,3) * t228 + t261;
t130 = -Icges(4,3) * t229 + t228 * t283;
t129 = t146 - t323;
t128 = t145 - t324;
t126 = t310 * t229;
t125 = t310 * t228;
t123 = t228 * t378 - t339;
t118 = Icges(5,3) * t228 + t260;
t117 = -Icges(5,3) * t229 + t228 * t282;
t96 = t102 + t231;
t95 = t101 - t385;
t92 = t296 + t344;
t91 = (-t232 - t237) * t229 + t228 * (t186 - t222);
t88 = t231 + t90;
t87 = t89 - t385;
t82 = t228 * t250 + t234 * t261;
t81 = t229 * t250 - t283 * t352;
t70 = t228 * t249 + t234 * t260;
t69 = t229 * t249 - t282 * t352;
t68 = pkin(3) * t408 - t154 * t228 - t172 * t350;
t67 = t172 * t352 + t188 + (-pkin(3) * t330 - t154) * t229;
t49 = t55 - t323;
t48 = t54 - t324;
t41 = t131 * t228 - t229 * t273;
t40 = t130 * t228 - t401;
t39 = -t131 * t229 - t404;
t38 = -t130 * t229 - t228 * t275;
t33 = t35 - t323;
t32 = t34 - t324;
t31 = t118 * t228 - t229 * t277;
t30 = t117 * t228 - t400;
t29 = -t118 * t229 - t403;
t28 = -t117 * t229 - t228 * t278;
t11 = -t114 * t352 + t329;
t10 = t228 * t91 + t229 * t92 + t369 + t62;
t3 = t228 * (t297 + t340) + t229 * (t295 + t345) + ((t91 + t349) * t229 + (-t112 - t114 - t92) * t228) * t234 + t328 + t329;
t1 = [t248 + (t20 * t76 + t75 * t21) * t390 + (t32 * t88 + t33 * t87) * t391 + (t48 * t96 + t49 * t95) * t392 + (t128 * t149 + t129 * t148) * t393 - t169 * t335 - t197 * t331; t243 + m(6) * (t20 * t80 + t21 * t79 + t22 * t76 + t23 * t75) + m(5) * (t32 * t90 + t33 * t89 + t34 * t88 + t35 * t87) + m(4) * (t101 * t49 + t102 * t48 + t54 * t96 + t55 * t95) + m(3) * (t128 * t174 + t129 * t173 + t145 * t149 + t146 * t148); (t22 * t80 + t23 * t79) * t390 + (t34 * t90 + t35 * t89) * t391 + (t101 * t55 + t102 * t54) * t392 + (t145 * t174 + t146 * t173) * t393 + t243; t242 + ((-t234 * t96 - t49) * t229 + (t234 * t95 - t48) * t228) * t386 + (-t228 * t96 - t229 * t95) * t387 + m(6) * (t20 * t93 + t21 * t94 + t36 * t75 + t37 * t76) + m(5) * (t125 * t32 + t126 * t33 + t67 * t87 + t68 * t88); m(6) * (t22 * t93 + t23 * t94 + t36 * t79 + t37 * t80) + m(5) * (t125 * t34 + t126 * t35 + t67 * t89 + t68 * t90) + t242 + ((-t102 * t234 - t55) * t229 + (t101 * t234 - t54) * t228) * t386 + (-t101 * t229 - t102 * t228) * t387; (t10 * t3 + t36 * t94 + t37 * t93) * t390 + (t126 * t67 + t125 * t68 + (t123 * t228 + t124 * t229 + t369) * ((t123 * t234 - t172 * t332 + t343) * t229 + ((-t112 - t124 + t211) * t234 + t319) * t228 + t328)) * t391 - t229 * ((t229 * t82 + (t39 + t401) * t234) * t229 + (t38 * t234 + (-t133 * t330 - t135 * t331) * t228 + (t276 * qJD(3) - t234 * t273 - t81) * t229) * t228) - t229 * ((t229 * t70 + (t29 + t400) * t234) * t229 + (t28 * t234 + (-t120 * t334 - t122 * t335) * t228 + (-t69 + (qJD(3) * t119 + t122 * t234) * t226 + (qJD(3) * t121 - t120 * t234) * t225) * t229) * t228) + t228 * ((t228 * t81 + (t40 + t404) * t234) * t228 + (t41 * t234 + (t132 * t330 + t134 * t331) * t229 + (-t274 * qJD(3) - t234 * t275 - t82) * t228) * t229) + ((t136 * t228 + t137 * t229) * (((-t137 + t212) * t234 + t318) * t228 + (t136 * t234 - t204 * t332 + t342) * t229) + t336 * t204 * t187) * t392 + t228 * ((t228 * t69 + (t30 + t403) * t234) * t228 + (t31 * t234 + (t119 * t334 + t121 * t335) * t229 + (-t70 + (-qJD(3) * t120 + t121 * t234) * t226 + (-qJD(3) * t122 - t119 * t234) * t225) * t228) * t229) + t294 + ((-t28 - t38) * t229 + (t29 + t39) * t228) * t352 + ((-t30 - t40) * t229 + (t31 + t41) * t228) * t350; m(6) * (t228 * t308 + t229 * t309) + m(5) * ((t234 * t87 - t32) * t229 + (t234 * t88 + t33) * t228); m(6) * (t228 * t306 + t229 * t307) + m(5) * ((t234 * t89 - t34) * t229 + (t234 * t90 + t35) * t228); m(6) * (t228 * t305 + t229 * t304) + m(5) * ((t126 * t234 - t68) * t229 + (t125 * t234 + t67) * t228); 0; m(6) * ((-t228 * t76 - t229 * t75) * t141 + (t228 * t309 - t229 * t308) * t158) + t268; m(6) * ((-t228 * t80 - t229 * t79) * t141 + (t228 * t307 - t229 * t306) * t158) + t268; m(6) * (t11 * t10 + t62 * t3 + (-t228 * t93 - t229 * t94) * t141 + (t228 * t304 - t229 * t305) * t158) + t294; 0; (t141 * t158 * t336 + t62 * t11) * t390 + t294;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;

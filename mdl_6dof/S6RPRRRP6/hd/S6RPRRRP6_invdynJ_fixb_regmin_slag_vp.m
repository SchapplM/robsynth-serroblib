% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:16:47
% EndTime: 2019-03-09 06:17:00
% DurationCPUTime: 5.05s
% Computational Cost: add. (8331->456), mult. (19782->580), div. (0->0), fcn. (15449->14), ass. (0->245)
t220 = sin(qJ(4));
t352 = pkin(8) + pkin(9);
t281 = qJD(4) * t352;
t217 = cos(pkin(10));
t224 = cos(qJ(3));
t299 = t224 * t217;
t192 = qJD(1) * t299;
t216 = sin(pkin(10));
t221 = sin(qJ(3));
t310 = t216 * t221;
t278 = qJD(1) * t310;
t153 = t192 - t278;
t320 = t153 * t220;
t165 = t216 * t224 + t217 * t221;
t154 = t165 * qJD(1);
t110 = pkin(3) * t154 - pkin(8) * t153;
t223 = cos(qJ(4));
t341 = pkin(7) + qJ(2);
t182 = t341 * t216;
t166 = qJD(1) * t182;
t183 = t341 * t217;
t167 = qJD(1) * t183;
t367 = -t166 * t224 - t221 * t167;
t330 = t220 * t110 + t223 * t367;
t383 = -pkin(9) * t320 + t220 * t281 + t330;
t96 = t223 * t110;
t382 = pkin(4) * t154 - t220 * t367 + t96 + (-pkin(9) * t153 + t281) * t223;
t124 = qJD(3) * t223 - t220 * t154;
t125 = qJD(3) * t220 + t154 * t223;
t219 = sin(qJ(5));
t351 = cos(qJ(5));
t246 = t219 * t124 + t351 * t125;
t353 = t246 ^ 2;
t73 = -t351 * t124 + t125 * t219;
t71 = t73 ^ 2;
t381 = -t71 + t353;
t144 = qJD(4) - t153;
t214 = pkin(10) + qJ(3);
t204 = sin(t214);
t205 = cos(t214);
t222 = sin(qJ(1));
t225 = cos(qJ(1));
t260 = g(1) * t225 + g(2) * t222;
t235 = -g(3) * t205 + t204 * t260;
t289 = qJD(1) * qJD(2);
t354 = t341 * qJDD(1) + t289;
t137 = t354 * t216;
t138 = t354 * t217;
t293 = qJD(3) * t224;
t294 = qJD(3) * t221;
t250 = -t224 * t137 - t221 * t138 + t166 * t294 - t167 * t293;
t326 = qJDD(3) * pkin(3);
t63 = -t250 - t326;
t380 = qJD(4) * pkin(8) * t144 - t235 + t63;
t292 = qJD(4) * t220;
t379 = t292 - t320;
t378 = t246 * t73;
t377 = t73 * qJ(6);
t279 = t351 * t223;
t307 = t219 * t220;
t245 = t279 - t307;
t361 = qJD(4) + qJD(5);
t273 = t351 * qJD(5);
t363 = t351 * qJD(4) + t273;
t329 = t245 * t153 - t363 * t223 + t361 * t307;
t280 = t351 * t220;
t169 = t219 * t223 + t280;
t118 = t361 * t169;
t328 = -t169 * t153 + t118;
t327 = qJDD(1) * pkin(1);
t366 = g(1) * t222 - g(2) * t225;
t252 = -qJDD(2) + t327 + t366;
t286 = t217 * qJDD(1);
t287 = t216 * qJDD(1);
t282 = qJD(3) * t192 + t221 * t286 + t224 * t287;
t238 = qJD(3) * t278 - t282;
t376 = qJD(3) * qJD(4) - t238;
t291 = qJD(4) * t223;
t276 = t165 * t291;
t164 = -t299 + t310;
t155 = t164 * qJD(3);
t305 = t220 * t155;
t375 = t276 - t305;
t140 = qJD(5) + t144;
t283 = t154 * t291 + t376 * t220;
t249 = t223 * qJDD(3) - t283;
t290 = qJD(5) * t219;
t69 = t220 * qJDD(3) - t154 * t292 + t376 * t223;
t26 = -t124 * t273 + t125 * t290 - t219 * t249 - t351 * t69;
t374 = t140 * t73 - t26;
t215 = qJ(4) + qJ(5);
t208 = cos(t215);
t312 = t208 * t222;
t207 = sin(t215);
t313 = t207 * t225;
t129 = -t205 * t312 + t313;
t311 = t208 * t225;
t314 = t207 * t222;
t131 = t205 * t311 + t314;
t115 = -t221 * t166 + t224 * t167;
t109 = qJD(3) * pkin(8) + t115;
t256 = -t137 * t221 + t138 * t224;
t62 = qJDD(3) * pkin(8) + qJD(3) * t367 + t256;
t156 = t165 * qJD(3);
t257 = t221 * t287 - t224 * t286;
t113 = qJD(1) * t156 + t257;
t196 = pkin(2) * t217 + pkin(1);
t176 = -qJDD(1) * t196 + qJDD(2);
t65 = t113 * pkin(3) + pkin(8) * t238 + t176;
t177 = -qJD(1) * t196 + qJD(2);
t84 = -pkin(3) * t153 - pkin(8) * t154 + t177;
t242 = -t109 * t292 + t220 * t65 + t223 * t62 + t84 * t291;
t10 = pkin(9) * t249 + t242;
t55 = -t109 * t220 + t223 * t84;
t43 = -pkin(9) * t125 + t55;
t33 = pkin(4) * t144 + t43;
t56 = t109 * t223 + t220 * t84;
t44 = pkin(9) * t124 + t56;
t107 = qJDD(4) + t113;
t61 = t223 * t65;
t7 = pkin(4) * t107 - pkin(9) * t69 - qJD(4) * t56 - t220 * t62 + t61;
t275 = -t351 * t10 - t219 * t7 - t33 * t273 + t44 * t290;
t343 = g(3) * t208;
t108 = -qJD(3) * pkin(3) - t367;
t70 = -pkin(4) * t124 + t108;
t373 = g(1) * t131 - g(2) * t129 + t204 * t343 + t70 * t73 + t275;
t36 = pkin(5) * t73 + qJD(6) + t70;
t372 = t246 * t36;
t370 = t379 * pkin(4) - t115;
t369 = qJ(6) * t246;
t368 = t204 * t366;
t122 = t182 * t224 + t221 * t183;
t184 = t352 * t220;
t185 = t352 * t223;
t297 = -t219 * t184 + t351 * t185;
t365 = -t297 * qJD(5) + t383 * t219 - t382 * t351;
t364 = t184 * t273 + t185 * t290 + t382 * t219 + t383 * t351;
t362 = qJ(2) * qJDD(1);
t128 = t205 * t314 + t311;
t130 = -t205 * t313 + t312;
t345 = g(3) * t204;
t360 = -g(1) * t130 + g(2) * t128 + t207 * t345;
t40 = t351 * t44;
t19 = t219 * t33 + t40;
t233 = -qJD(5) * t19 - t219 * t10 + t351 * t7;
t359 = -t70 * t246 + t233 + t360;
t27 = qJD(5) * t246 + t219 * t69 - t351 * t249;
t358 = t140 * t246 - t27;
t104 = qJDD(5) + t107;
t357 = t104 * t169 - t140 * t329;
t356 = t260 * t205 + t345;
t355 = -t245 * t26 - t246 * t328;
t209 = t223 * pkin(4);
t342 = pkin(3) + t209;
t38 = t219 * t44;
t18 = t351 * t33 - t38;
t12 = t18 - t369;
t11 = pkin(5) * t140 + t12;
t340 = -t12 + t11;
t339 = t351 * t43 - t38;
t338 = -t328 * qJ(6) + qJD(6) * t245 - t364;
t337 = -pkin(5) * t154 + t329 * qJ(6) - t169 * qJD(6) + t365;
t112 = pkin(3) * t164 - pkin(8) * t165 - t196;
t101 = t223 * t112;
t123 = -t182 * t221 + t183 * t224;
t318 = t165 * t223;
t51 = pkin(4) * t164 - pkin(9) * t318 - t123 * t220 + t101;
t116 = t223 * t123;
t298 = t220 * t112 + t116;
t319 = t165 * t220;
t58 = -pkin(9) * t319 + t298;
t335 = t219 * t51 + t351 * t58;
t334 = t154 * t73;
t333 = t154 * t246;
t331 = t69 * t220;
t324 = t124 * t144;
t323 = t124 * t154;
t322 = t125 * t144;
t321 = t125 * t154;
t178 = pkin(4) * t220 + pkin(5) * t207;
t316 = t178 * t205;
t306 = t220 * t107;
t304 = t220 * t222;
t303 = t220 * t225;
t302 = t222 * t223;
t93 = t223 * t107;
t301 = t223 * t155;
t300 = t223 * t225;
t296 = t178 + t341;
t179 = pkin(5) * t208 + t209;
t295 = t216 ^ 2 + t217 ^ 2;
t277 = t165 * t292;
t271 = -t219 * t43 - t40;
t269 = -t219 * t58 + t351 * t51;
t267 = -qJD(4) * t84 - t62;
t266 = t295 * qJD(1) ^ 2;
t265 = -t351 * t184 - t185 * t219;
t264 = t144 * t223;
t86 = t165 * qJD(2) - t182 * t294 + t183 * t293;
t263 = -t169 * t27 + t329 * t73;
t262 = 0.2e1 * t295;
t261 = t245 * t104 - t328 * t140;
t258 = -t109 * t291 + t61;
t87 = pkin(4) * t319 + t122;
t173 = pkin(3) + t179;
t211 = -qJ(6) - t352;
t254 = t173 * t205 - t204 * t211;
t251 = -t379 * t144 + t93;
t64 = t375 * pkin(4) + t86;
t248 = t196 + t254;
t85 = -t164 * qJD(2) - t122 * qJD(3);
t111 = pkin(3) * t156 + pkin(8) * t155;
t97 = t223 * t111;
t25 = pkin(9) * t301 + pkin(4) * t156 - t220 * t85 + t97 + (-t116 + (pkin(9) * t165 - t112) * t220) * qJD(4);
t241 = t220 * t111 + t112 * t291 - t123 * t292 + t223 * t85;
t29 = -pkin(9) * t375 + t241;
t247 = t219 * t25 + t51 * t273 + t351 * t29 - t290 * t58;
t243 = -t277 - t301;
t240 = -pkin(8) * t107 + t108 * t144;
t237 = t252 + t327;
t232 = t262 * t289 - t260;
t231 = -t335 * qJD(5) - t219 * t29 + t351 * t25;
t30 = -pkin(4) * t249 + t63;
t8 = t27 * pkin(5) + qJDD(6) + t30;
t201 = t351 * pkin(4) + pkin(5);
t148 = t205 * t300 + t304;
t147 = -t205 * t303 + t302;
t146 = -t205 * t302 + t303;
t145 = t205 * t304 + t300;
t99 = t245 * t165;
t98 = t169 * t165;
t92 = qJ(6) * t245 + t297;
t91 = -qJ(6) * t169 + t265;
t35 = -t155 * t280 - t219 * t277 - t290 * t319 + (-t155 * t219 + t363 * t165) * t223;
t34 = t118 * t165 + t155 * t279 - t219 * t305;
t21 = -qJ(6) * t98 + t335;
t20 = pkin(5) * t164 - qJ(6) * t99 + t269;
t15 = t339 - t369;
t14 = t271 + t377;
t13 = t19 - t377;
t4 = -qJ(6) * t35 - qJD(6) * t98 + t247;
t3 = t156 * pkin(5) + t34 * qJ(6) - t99 * qJD(6) + t231;
t2 = -qJ(6) * t27 - qJD(6) * t73 - t275;
t1 = t104 * pkin(5) + t26 * qJ(6) - qJD(6) * t246 + t233;
t5 = [qJDD(1), t366, t260, t237 * t217, -t237 * t216, t262 * t362 + t232, pkin(1) * t252 + (t295 * t362 + t232) * qJ(2), -t154 * t155 - t165 * t238, -t165 * t113 - t155 * t153 - t154 * t156 + t164 * t238, -qJD(3) * t155 + qJDD(3) * t165, -qJD(3) * t156 - qJDD(3) * t164, 0, -qJD(3) * t86 - t122 * qJDD(3) - t113 * t196 + t156 * t177 + t164 * t176 + t205 * t366, -t85 * qJD(3) - t123 * qJDD(3) - t177 * t155 + t176 * t165 + t196 * t238 - t368, t125 * t243 + t69 * t318 -(t124 * t223 - t125 * t220) * t155 + (t223 * t249 - t331 + (-t124 * t220 - t125 * t223) * qJD(4)) * t165, t125 * t156 + t144 * t243 + t164 * t69 + t165 * t93, t124 * t156 - t144 * t375 + t164 * t249 - t165 * t306, t107 * t164 + t144 * t156 (-t123 * t291 + t97) * t144 + t101 * t107 + t258 * t164 + t55 * t156 - t86 * t124 - t122 * t249 + t108 * t276 - g(1) * t146 - g(2) * t148 + ((-qJD(4) * t112 - t85) * t144 - t123 * t107 + t267 * t164 + t63 * t165 - t108 * t155) * t220, -g(1) * t145 - g(2) * t147 - t298 * t107 + t243 * t108 + t122 * t69 + t86 * t125 - t241 * t144 - t56 * t156 - t242 * t164 + t63 * t318, -t246 * t34 - t26 * t99, -t246 * t35 + t26 * t98 - t27 * t99 + t34 * t73, t104 * t99 - t140 * t34 + t156 * t246 - t164 * t26, -t104 * t98 - t140 * t35 - t156 * t73 - t164 * t27, t104 * t164 + t140 * t156, -g(1) * t129 - g(2) * t131 + t104 * t269 + t140 * t231 + t18 * t156 + t164 * t233 + t87 * t27 + t30 * t98 + t70 * t35 + t64 * t73, -g(1) * t128 - g(2) * t130 - t335 * t104 - t247 * t140 - t19 * t156 + t275 * t164 + t246 * t64 - t87 * t26 + t30 * t99 - t70 * t34, -t1 * t99 + t11 * t34 - t13 * t35 - t2 * t98 + t20 * t26 - t21 * t27 - t246 * t3 - t4 * t73 + t368, t2 * t21 + t13 * t4 + t1 * t20 + t11 * t3 + t8 * (pkin(5) * t98 + t87) + t36 * (pkin(5) * t35 + t64) + (-g(1) * t296 - g(2) * t248) * t225 + (g(1) * t248 - g(2) * t296) * t222; 0, 0, 0, -t286, t287, -t266, -qJ(2) * t266 - t252, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t154 + t257 (t153 - t278) * qJD(3) + t282, 0, 0, 0, 0, 0, t251 + t323, -t144 ^ 2 * t223 - t306 - t321, 0, 0, 0, 0, 0, t261 - t334, -t333 - t357, t263 - t355, t1 * t245 - t11 * t328 - t13 * t329 - t154 * t36 + t169 * t2 - t366; 0, 0, 0, 0, 0, 0, 0, -t154 * t153, -t153 ^ 2 + t154 ^ 2 (-t153 - t278) * qJD(3) + t282, -t257, qJDD(3), qJD(3) * t115 - t154 * t177 + t235 + t250, -t153 * t177 - t256 + t356, t125 * t264 + t331 (t69 + t324) * t223 + (t249 - t322) * t220, t144 * t264 + t306 - t321, t251 - t323, -t144 * t154, -pkin(3) * t283 - t96 * t144 - t55 * t154 + t115 * t124 + (t144 * t367 + t240) * t220 + (t326 - t380) * t223, -pkin(3) * t69 - t115 * t125 + t330 * t144 + t56 * t154 + t380 * t220 + t240 * t223, -t169 * t26 - t246 * t329, t263 + t355, -t333 + t357, t261 + t334, -t140 * t154, t104 * t265 - t18 * t154 - t30 * t245 - t342 * t27 - t205 * t343 + t370 * t73 + t328 * t70 + (g(1) * t311 + g(2) * t312) * t204 + t365 * t140, -t297 * t104 + t364 * t140 + t19 * t154 + t30 * t169 - t207 * t235 + t370 * t246 + t26 * t342 - t329 * t70, -t1 * t169 + t11 * t329 - t13 * t328 + t2 * t245 - t246 * t337 + t26 * t91 - t27 * t92 - t338 * t73 - t356, t2 * t92 + t1 * t91 + t8 * (-pkin(5) * t245 - t342) - g(3) * t254 + (pkin(5) * t328 + t370) * t36 + t338 * t13 + t337 * t11 + t260 * (t173 * t204 + t205 * t211); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125 * t124, -t124 ^ 2 + t125 ^ 2, t69 - t324, t249 + t322, t107, -g(1) * t147 + g(2) * t145 - t108 * t125 + t144 * t56 + (t267 + t345) * t220 + t258, g(1) * t148 - g(2) * t146 - t108 * t124 + t144 * t55 + t223 * t345 - t242, t378, t381, t374, t358, t104, -t271 * t140 + (t104 * t351 - t125 * t73 - t140 * t290) * pkin(4) + t359, t339 * t140 + (-t219 * t104 - t125 * t246 - t140 * t273) * pkin(4) + t373, -t11 * t73 + t13 * t246 + t14 * t246 + t15 * t73 + t201 * t26 + (-t219 * t27 + (t219 * t246 - t351 * t73) * qJD(5)) * pkin(4), t1 * t201 - t13 * t15 - t11 * t14 - pkin(5) * t372 - g(1) * (t179 * t222 - t225 * t316) - g(2) * (-t179 * t225 - t222 * t316) + t178 * t345 + (-t36 * t125 + t2 * t219 + (-t11 * t219 + t13 * t351) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t378, t381, t374, t358, t104, t19 * t140 + t359, t140 * t18 + t373, pkin(5) * t26 - t340 * t73, t340 * t13 + (t1 + t360 - t372) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71 - t353, t11 * t246 + t13 * t73 - t235 + t8;];
tau_reg  = t5;

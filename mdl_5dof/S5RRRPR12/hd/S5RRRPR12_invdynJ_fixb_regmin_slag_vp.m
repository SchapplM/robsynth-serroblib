% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:38
% EndTime: 2019-12-31 21:40:57
% DurationCPUTime: 8.07s
% Computational Cost: add. (7227->551), mult. (18224->793), div. (0->0), fcn. (14549->14), ass. (0->248)
t243 = cos(qJ(3));
t348 = cos(pkin(5));
t293 = t348 * qJD(1);
t270 = t293 + qJD(2);
t239 = sin(qJ(3));
t240 = sin(qJ(2));
t236 = sin(pkin(5));
t323 = qJD(1) * t236;
t306 = t240 * t323;
t287 = t239 * t306;
t152 = -t243 * t270 + t287;
t145 = qJD(5) + t152;
t154 = t239 * t270 + t243 * t306;
t244 = cos(qJ(2));
t322 = qJD(1) * t244;
t305 = t236 * t322;
t206 = -qJD(3) + t305;
t235 = sin(pkin(10));
t237 = cos(pkin(10));
t110 = t154 * t235 + t237 * t206;
t242 = cos(qJ(5));
t112 = t154 * t237 - t206 * t235;
t238 = sin(qJ(5));
t347 = t112 * t238;
t373 = -t242 * t110 - t347;
t375 = t373 * t145;
t286 = pkin(1) * t293;
t172 = pkin(7) * t305 + t240 * t286;
t374 = t239 * qJD(4) + t172 + t206 * (pkin(3) * t239 - qJ(4) * t243);
t271 = t110 * t238 - t112 * t242;
t372 = t145 * t271;
t313 = t240 * qJDD(1);
t298 = t236 * t313;
t315 = qJD(1) * qJD(2);
t299 = t236 * t315;
t370 = -t244 * t299 - t298;
t241 = sin(qJ(1));
t363 = cos(qJ(1));
t282 = t348 * t363;
t181 = t240 * t282 + t241 * t244;
t307 = t236 * t363;
t128 = t181 * t243 - t239 * t307;
t180 = t240 * t241 - t244 * t282;
t232 = pkin(10) + qJ(5);
t229 = sin(t232);
t230 = cos(t232);
t369 = t128 * t229 - t180 * t230;
t368 = t128 * t230 + t180 * t229;
t291 = t348 * qJDD(1);
t220 = t291 + qJDD(2);
t255 = qJD(3) * t270;
t320 = qJD(2) * t244;
t301 = t239 * t320;
t318 = qJD(3) * t243;
t76 = -t243 * t220 + (qJD(1) * (t240 * t318 + t301) + t239 * t313) * t236 + t239 * t255;
t231 = t236 ^ 2;
t312 = 0.2e1 * t231;
t319 = qJD(3) * t239;
t311 = pkin(8) * t319;
t169 = -pkin(7) * t306 + t244 * t286;
t280 = pkin(2) * t240 - pkin(8) * t244;
t170 = t280 * t323;
t326 = t243 * t169 + t239 * t170;
t87 = qJ(4) * t306 + t326;
t353 = -t374 * t237 + (t311 + t87) * t235;
t367 = t374 * t235 + t237 * t87;
t335 = t236 * t240;
t221 = pkin(7) * t335;
t294 = t244 * t348;
t366 = pkin(1) * t294 - t221;
t314 = qJDD(1) * t244;
t219 = t236 * t314;
t285 = t240 * t299;
t166 = qJDD(3) - t219 + t285;
t365 = -pkin(3) * t166 + qJDD(4);
t75 = -qJD(3) * t287 + t239 * t220 + (t255 - t370) * t243;
t56 = -t237 * t166 + t235 * t75;
t57 = t166 * t235 + t237 * t75;
t13 = -qJD(5) * t271 + t238 * t57 + t242 * t56;
t245 = qJD(1) ^ 2;
t361 = pkin(8) * t240;
t182 = t363 * t240 + t241 * t294;
t360 = g(1) * t182;
t359 = g(2) * t180;
t358 = pkin(9) + qJ(4);
t138 = pkin(8) * t270 + t172;
t268 = -pkin(2) * t244 - pkin(1) - t361;
t165 = t268 * t236;
t144 = qJD(1) * t165;
t267 = qJD(2) * t286;
t283 = pkin(1) * t291;
t309 = -pkin(7) * t219 - t240 * t283 - t244 * t267;
t250 = -pkin(7) * t285 - t309;
t96 = pkin(8) * t220 + t250;
t264 = t280 * qJD(2);
t99 = (qJD(1) * t264 + qJDD(1) * t268) * t236;
t260 = -t138 * t319 + t144 * t318 + t239 * t99 + t243 * t96;
t21 = qJ(4) * t166 - qJD(4) * t206 + t260;
t289 = t370 * pkin(7) - t240 * t267 + t244 * t283;
t97 = -pkin(2) * t220 - t289;
t23 = pkin(3) * t76 - qJ(4) * t75 - qJD(4) * t154 + t97;
t9 = t237 * t21 + t235 * t23;
t295 = t240 * t348;
t332 = t236 * t244;
t325 = pkin(1) * t295 + pkin(7) * t332;
t164 = t348 * pkin(8) + t325;
t171 = t236 * t264;
t173 = t366 * qJD(2);
t259 = -t164 * t319 + t165 * t318 + t239 * t171 + t243 * t173;
t321 = qJD(2) * t240;
t42 = (qJ(4) * t321 - qJD(4) * t244) * t236 + t259;
t333 = t236 * t243;
t179 = t348 * t239 + t240 * t333;
t124 = qJD(3) * t179 + t236 * t301;
t178 = t239 * t335 - t348 * t243;
t302 = t236 * t320;
t125 = -qJD(3) * t178 + t243 * t302;
t174 = t325 * qJD(2);
t48 = t124 * pkin(3) - t125 * qJ(4) - t179 * qJD(4) + t174;
t17 = t235 * t48 + t237 * t42;
t137 = -pkin(2) * t270 - t169;
t62 = t152 * pkin(3) - t154 * qJ(4) + t137;
t73 = t243 * t138 + t239 * t144;
t65 = -qJ(4) * t206 + t73;
t30 = t235 * t62 + t237 * t65;
t72 = -t239 * t138 + t144 * t243;
t98 = pkin(3) * t154 + qJ(4) * t152;
t39 = t235 * t98 + t237 * t72;
t163 = -t348 * pkin(2) - t366;
t84 = t178 * pkin(3) - t179 * qJ(4) + t163;
t327 = t243 * t164 + t239 * t165;
t85 = -qJ(4) * t332 + t327;
t37 = t235 * t84 + t237 * t85;
t357 = pkin(8) * qJD(3);
t356 = qJ(4) * t76;
t336 = t235 * t243;
t139 = -t237 * t306 + t305 * t336;
t329 = t243 * t244;
t140 = (t235 * t240 + t237 * t329) * t323;
t337 = t235 * t238;
t188 = -t242 * t237 + t337;
t189 = t235 * t242 + t237 * t238;
t317 = qJD(5) * t239;
t355 = t139 * t238 - t140 * t242 - t188 * t318 - t189 * t317;
t316 = qJD(5) * t242;
t331 = t237 * t239;
t354 = -t242 * t139 - t140 * t238 + t189 * t318 + t316 * t331 - t317 * t337;
t352 = -t237 * t311 - t367;
t351 = t145 * t188;
t350 = t145 * t189;
t308 = pkin(4) * t235 + pkin(8);
t155 = t239 * t169;
t88 = -pkin(3) * t306 - t170 * t243 + t155;
t349 = -pkin(4) * t139 + t308 * t318 - t88;
t346 = t152 * t206;
t345 = t152 * t235;
t344 = t154 * t206;
t341 = t206 * t239;
t340 = t229 * t243;
t339 = t230 * t243;
t338 = t231 * t245;
t334 = t236 * t241;
t330 = t237 * t243;
t64 = pkin(3) * t206 + qJD(4) - t72;
t328 = -qJD(4) + t64;
t266 = pkin(3) * t243 + qJ(4) * t239 + pkin(2);
t151 = pkin(8) * t330 - t235 * t266;
t233 = t240 ^ 2;
t324 = -t244 ^ 2 + t233;
t310 = t244 * t338;
t304 = t239 * t322;
t303 = t236 * t321;
t300 = t244 * t315;
t8 = -t21 * t235 + t237 * t23;
t16 = -t235 * t42 + t237 * t48;
t29 = -t235 * t65 + t237 * t62;
t36 = -t235 * t85 + t237 * t84;
t38 = -t235 * t72 + t237 * t98;
t292 = -t239 * t164 + t165 * t243;
t290 = t138 * t318 + t144 * t319 + t239 * t96 - t243 * t99;
t281 = t236 * t245 * t348;
t127 = t181 * t239 + t243 * t307;
t183 = -t241 * t295 + t363 * t244;
t131 = t183 * t239 - t241 * t333;
t279 = -g(1) * t127 + g(2) * t131;
t278 = g(1) * t183 + g(2) * t181;
t4 = pkin(4) * t76 - pkin(9) * t57 + t8;
t5 = -pkin(9) * t56 + t9;
t277 = t238 * t4 + t242 * t5;
t126 = -pkin(9) * t235 * t239 + t151;
t276 = pkin(4) * t236 * t304 - pkin(9) * t140 + qJD(5) * t126 - (pkin(4) * t239 - pkin(9) * t330) * qJD(3) - t353;
t186 = t237 * t266;
t115 = -pkin(9) * t331 - t186 + (-pkin(8) * t235 - pkin(4)) * t243;
t275 = -pkin(9) * t139 - qJD(5) * t115 - (-pkin(8) * t331 - pkin(9) * t336) * qJD(3) + t367;
t86 = pkin(3) * t332 - t292;
t15 = pkin(4) * t152 - pkin(9) * t112 + t29;
t22 = -pkin(9) * t110 + t30;
t6 = t15 * t242 - t22 * t238;
t7 = t15 * t238 + t22 * t242;
t123 = t179 * t237 - t235 * t332;
t25 = pkin(4) * t178 - pkin(9) * t123 + t36;
t122 = t179 * t235 + t237 * t332;
t31 = -pkin(9) * t122 + t37;
t273 = -t238 * t31 + t242 * t25;
t272 = t238 * t25 + t242 * t31;
t66 = t242 * t122 + t123 * t238;
t67 = -t122 * t238 + t123 * t242;
t269 = 0.2e1 * t293 + qJD(2);
t265 = -t164 * t318 - t165 * t319 + t171 * t243 - t239 * t173;
t201 = t358 * t235;
t263 = pkin(9) * t345 - qJD(4) * t237 + qJD(5) * t201 + t39;
t202 = t358 * t237;
t262 = pkin(9) * t152 * t237 + pkin(4) * t154 + qJD(4) * t235 + qJD(5) * t202 + t38;
t12 = -qJD(5) * t347 - t110 * t316 - t238 * t56 + t242 * t57;
t258 = t236 * (t291 + t220);
t257 = g(1) * t131 + g(2) * t127 + g(3) * t178;
t132 = t183 * t243 + t239 * t334;
t256 = -g(1) * t132 - g(2) * t128 - g(3) * t179;
t24 = t290 + t365;
t253 = -t24 + t257;
t252 = g(3) * t332 - t359 - t360;
t251 = -g(3) * t335 - t278;
t249 = -t252 - t97;
t248 = -pkin(8) * t166 - t137 * t206;
t47 = -pkin(3) * t303 - t265;
t2 = -qJD(5) * t7 - t238 * t5 + t242 * t4;
t246 = t257 - t290;
t228 = -pkin(4) * t237 - pkin(3);
t190 = t308 * t239;
t168 = t188 * t239;
t167 = t189 * t239;
t150 = -pkin(8) * t336 - t186;
t104 = t125 * t237 + t235 * t303;
t103 = t125 * t235 - t237 * t303;
t82 = t132 * t230 + t182 * t229;
t81 = -t132 * t229 + t182 * t230;
t69 = qJDD(5) + t76;
t58 = -pkin(4) * t345 + t73;
t55 = pkin(4) * t122 + t86;
t41 = pkin(4) * t110 + t64;
t32 = pkin(4) * t103 + t47;
t27 = qJD(5) * t67 + t242 * t103 + t238 * t104;
t26 = -qJD(5) * t66 - t238 * t103 + t242 * t104;
t14 = -pkin(9) * t103 + t17;
t11 = pkin(4) * t56 + t24;
t10 = pkin(4) * t124 - pkin(9) * t104 + t16;
t1 = t6 * qJD(5) + t277;
t3 = [qJDD(1), g(1) * t241 - g(2) * t363, g(1) * t363 + g(2) * t241, (qJDD(1) * t233 + 0.2e1 * t240 * t300) * t231, (t244 * t313 - t324 * t315) * t312, t240 * t258 + t269 * t302, t244 * t258 - t269 * t303, t220 * t348, -t174 * t270 - t221 * t220 + t289 * t348 + g(1) * t181 - g(2) * t183 + (t220 * t294 + (-t240 * t315 + t314) * t312) * pkin(1), -t173 * t270 - t325 * t220 - t250 * t348 - g(1) * t180 + g(2) * t182 + (-t300 - t313) * pkin(1) * t312, t125 * t154 + t179 * t75, -t124 * t154 - t125 * t152 - t178 * t75 - t179 * t76, -t125 * t206 + t179 * t166 + (t154 * t321 - t244 * t75) * t236, t124 * t206 - t178 * t166 + (-t152 * t321 + t244 * t76) * t236, (-t166 * t244 - t206 * t321) * t236, -t265 * t206 + t292 * t166 + t174 * t152 + t163 * t76 + t97 * t178 + t137 * t124 + g(1) * t128 - g(2) * t132 + (t244 * t290 + t321 * t72) * t236, t259 * t206 - t327 * t166 + t174 * t154 + t163 * t75 + t97 * t179 + t137 * t125 + (t244 * t260 - t321 * t73) * t236 + t279, t16 * t152 + t36 * t76 + t8 * t178 + t29 * t124 + t47 * t110 + t86 * t56 + t24 * t122 + t64 * t103 - g(1) * (-t128 * t237 - t180 * t235) - g(2) * (t132 * t237 + t182 * t235), -t17 * t152 - t37 * t76 - t9 * t178 - t30 * t124 + t47 * t112 + t86 * t57 + t24 * t123 + t64 * t104 - g(1) * (t128 * t235 - t180 * t237) - g(2) * (-t132 * t235 + t182 * t237), -t103 * t30 - t104 * t29 - t110 * t17 - t112 * t16 - t122 * t9 - t123 * t8 - t36 * t57 - t37 * t56 - t279, t9 * t37 + t30 * t17 + t8 * t36 + t29 * t16 + t24 * t86 + t64 * t47 - g(1) * (-t241 * pkin(1) - t181 * pkin(2) - pkin(3) * t128 + pkin(7) * t307 - t180 * pkin(8) - qJ(4) * t127) - g(2) * (t363 * pkin(1) + t183 * pkin(2) + t132 * pkin(3) + pkin(7) * t334 + t182 * pkin(8) + t131 * qJ(4)), t12 * t67 - t26 * t271, -t12 * t66 - t13 * t67 + t26 * t373 + t27 * t271, t12 * t178 - t124 * t271 + t145 * t26 + t67 * t69, t124 * t373 - t13 * t178 - t145 * t27 - t66 * t69, t124 * t145 + t178 * t69, (-qJD(5) * t272 + t242 * t10 - t238 * t14) * t145 + t273 * t69 + t2 * t178 + t6 * t124 - t32 * t373 + t55 * t13 + t11 * t66 + t41 * t27 + g(1) * t368 - g(2) * t82, -t32 * t271 + t55 * t12 + t11 * t67 + t41 * t26 - (qJD(5) * t273 + t238 * t10 + t242 * t14) * t145 - t272 * t69 - t1 * t178 - t7 * t124 - g(1) * t369 - g(2) * t81; 0, 0, 0, -t240 * t310, t324 * t338, -t244 * t281 + t298, t240 * t281 + t219, t220, pkin(1) * t240 * t338 + t172 * t270 - t252 + t289, pkin(1) * t310 + t169 * t270 + (pkin(7) * t315 + g(3)) * t335 + t278 + t309, t75 * t239 - t243 * t344, (t75 + t346) * t243 + (-t76 + t344) * t239, -t206 * t318 + t239 * t166 + (-t154 * t240 + t206 * t329) * t323, t206 * t319 + t243 * t166 + (t152 * t240 - t244 * t341) * t323, t206 * t306, -t72 * t306 - pkin(2) * t76 - t172 * t152 - t155 * t206 + t248 * t239 + ((t170 + t357) * t206 + t249) * t243, -pkin(2) * t75 - t326 * t206 + t73 * t306 - t172 * t154 + t248 * t243 + (-t206 * t357 - t249) * t239, -t88 * t110 - t64 * t139 + t150 * t76 + t353 * t152 + t251 * t235 + (-t8 + (pkin(8) * t110 + t235 * t64) * qJD(3) - t252 * t237) * t243 + (pkin(8) * t56 - t206 * t29 + t24 * t235) * t239, -t88 * t112 - t64 * t140 - t151 * t76 - t352 * t152 + t251 * t237 + (t9 + (pkin(8) * t112 + t237 * t64) * qJD(3) + t252 * t235) * t243 + (pkin(8) * t57 + t206 * t30 + t24 * t237) * t239, t30 * t139 + t29 * t140 - t150 * t57 - t151 * t56 - t353 * t112 - t352 * t110 + (-t235 * t30 - t237 * t29) * t318 + (-t235 * t9 - t237 * t8 - t252) * t239, t8 * t150 + t9 * t151 - t64 * t88 + t352 * t30 + t353 * t29 + t266 * t360 + t266 * t359 + (t239 * t24 + t318 * t64 - t278) * pkin(8) - g(3) * (t244 * t266 + t361) * t236, -t12 * t168 - t271 * t355, -t12 * t167 + t168 * t13 + t271 * t354 + t355 * t373, -t12 * t243 + t145 * t355 - t168 * t69 + t271 * t341, t13 * t243 - t145 * t354 - t167 * t69 - t341 * t373, -t145 * t341 - t69 * t243, (t115 * t242 - t126 * t238) * t69 - t2 * t243 + t6 * t319 + t190 * t13 + t11 * t167 - g(1) * (-t182 * t339 + t183 * t229) - g(2) * (-t180 * t339 + t181 * t229) - t349 * t373 + t354 * t41 + (t238 * t275 - t242 * t276) * t145 + (-t6 * t304 - g(3) * (t229 * t240 + t230 * t329)) * t236, t190 * t12 - t11 * t168 - (t115 * t238 + t126 * t242) * t69 + t1 * t243 - t7 * t319 - g(1) * (t182 * t340 + t183 * t230) - g(2) * (t180 * t340 + t181 * t230) - t349 * t271 + t355 * t41 + (t238 * t276 + t242 * t275) * t145 + (t7 * t304 - g(3) * (-t229 * t329 + t230 * t240)) * t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154 * t152, -t152 ^ 2 + t154 ^ 2, t75 - t346, -t344 - t76, t166, -t137 * t154 - t206 * t73 + t246, t137 * t152 - t206 * t72 - t256 - t260, -t235 * t356 - pkin(3) * t56 - t73 * t110 - t29 * t154 + (t235 * t328 - t38) * t152 + t253 * t237, -t237 * t356 - pkin(3) * t57 - t73 * t112 + t30 * t154 + (t237 * t328 + t39) * t152 - t253 * t235, t39 * t110 + t38 * t112 + (-qJ(4) * t56 - qJD(4) * t110 - t152 * t29 + t9) * t237 + (qJ(4) * t57 + qJD(4) * t112 - t152 * t30 - t8) * t235 + t256, -t29 * t38 - t30 * t39 - t64 * t73 + (-t235 * t29 + t237 * t30) * qJD(4) + t253 * pkin(3) + (-t8 * t235 + t9 * t237 + t256) * qJ(4), t12 * t189 + t271 * t351, -t12 * t188 - t189 * t13 + t271 * t350 - t351 * t373, -t145 * t351 + t154 * t271 + t189 * t69, -t145 * t350 - t154 * t373 - t188 * t69, -t145 * t154, (-t201 * t242 - t202 * t238) * t69 + t228 * t13 + t11 * t188 - t6 * t154 + t58 * t373 + t350 * t41 + (t238 * t263 - t242 * t262) * t145 + t257 * t230, t228 * t12 + t11 * t189 - (-t201 * t238 + t202 * t242) * t69 + t58 * t271 + t7 * t154 - t351 * t41 + (t238 * t262 + t242 * t263) * t145 - t257 * t229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112 * t152 + t56, -t110 * t152 + t57, -t110 ^ 2 - t112 ^ 2, t110 * t30 + t112 * t29 - t246 + t365, 0, 0, 0, 0, 0, t13 - t372, t12 + t375; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271 * t373, t271 ^ 2 - t373 ^ 2, t12 - t375, -t13 - t372, t69, t7 * t145 + t41 * t271 - g(1) * t81 + g(2) * t369 - g(3) * (-t179 * t229 - t230 * t332) + t2, -t41 * t373 + g(1) * t82 + g(2) * t368 - g(3) * (-t179 * t230 + t229 * t332) - t277 + (t145 - qJD(5)) * t6;];
tau_reg = t3;

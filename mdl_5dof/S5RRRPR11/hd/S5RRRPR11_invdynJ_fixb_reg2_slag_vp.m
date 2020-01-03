% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPR11
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR11_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:20
% EndTime: 2019-12-31 21:35:33
% DurationCPUTime: 6.58s
% Computational Cost: add. (5602->575), mult. (12248->716), div. (0->0), fcn. (8039->8), ass. (0->275)
t209 = cos(qJ(3));
t206 = sin(qJ(2));
t315 = qJD(1) * t206;
t293 = t209 * t315;
t205 = sin(qJ(3));
t313 = qJD(2) * t205;
t132 = t293 + t313;
t301 = t206 * qJDD(1);
t188 = pkin(6) * t301;
t210 = cos(qJ(2));
t302 = qJD(1) * qJD(2);
t284 = t210 * t302;
t264 = qJDD(2) * pkin(2) - pkin(6) * t284 - t188;
t303 = t209 * qJD(2);
t309 = qJD(3) * t206;
t400 = qJD(1) * t309 - qJDD(2);
t402 = -t284 - t301;
t61 = -qJD(3) * t303 + t400 * t205 + t402 * t209;
t219 = -qJ(4) * t61 + qJD(4) * t132 + t264;
t377 = pkin(3) + pkin(4);
t314 = qJD(1) * t210;
t62 = (qJD(2) * (qJD(3) + t314) + t301) * t205 + t400 * t209;
t10 = -t377 * t62 + t219;
t207 = sin(qJ(1));
t211 = cos(qJ(1));
t267 = g(1) * t211 + g(2) * t207;
t384 = t206 * t267;
t411 = t10 + t384;
t195 = t206 * pkin(7);
t198 = t210 * pkin(2);
t297 = -pkin(1) - t198;
t246 = t297 - t195;
t123 = t246 * qJD(1);
t190 = pkin(6) * t314;
t152 = qJD(2) * pkin(7) + t190;
t72 = t209 * t123 - t205 * t152;
t323 = qJD(4) - t72;
t204 = sin(qJ(5));
t331 = t206 * t209;
t208 = cos(qJ(5));
t334 = t205 * t208;
t104 = t204 * t331 - t206 * t334;
t326 = t209 * t211;
t328 = t207 * t210;
t115 = t205 * t328 + t326;
t324 = t211 * t205;
t327 = t209 * t210;
t116 = t207 * t327 - t324;
t388 = t115 * t208 - t116 * t204;
t294 = t205 * t315;
t130 = t294 - t303;
t151 = -qJD(2) * pkin(2) + pkin(6) * t315;
t237 = qJ(4) * t132 - t151;
t43 = -t377 * t130 + t237;
t329 = t207 * t209;
t117 = t210 * t324 - t329;
t325 = t210 * t211;
t118 = t205 * t207 + t209 * t325;
t59 = t117 * t208 - t118 * t204;
t68 = t130 * t204 + t132 * t208;
t410 = -g(1) * t59 - g(2) * t388 + g(3) * t104 - t43 * t68;
t401 = -t314 + qJD(3);
t338 = t132 * t401;
t340 = t130 * t401;
t409 = (t62 + t338) * t205 + (t61 + t340) * t209;
t250 = -t208 * t130 + t132 * t204;
t362 = t68 * t250;
t270 = pkin(2) * t206 - pkin(7) * t210;
t137 = t270 * qJD(1);
t119 = t205 * t137;
t184 = qJ(4) * t315;
t310 = qJD(3) * t205;
t333 = t205 * t210;
t376 = pkin(7) - pkin(8);
t408 = t376 * t310 + t119 + t184 + (-pkin(6) * t331 + pkin(8) * t333) * qJD(1);
t154 = t376 * t209;
t370 = pkin(6) * t205;
t296 = -pkin(3) - t370;
t222 = -pkin(8) * t327 + (-pkin(4) + t296) * t206;
t337 = t137 * t209;
t407 = -qJD(1) * t222 + qJD(3) * t154 + t337;
t404 = t62 - t338;
t403 = -pkin(8) * t132 + t323;
t399 = -t250 ^ 2 + t68 ^ 2;
t14 = qJD(5) * t68 - t204 * t61 - t208 * t62;
t159 = qJD(5) - t401;
t392 = t159 * t68 - t14;
t311 = qJD(2) * t210;
t349 = t209 * t62;
t350 = t205 * t61;
t397 = ((t130 * t205 - t132 * t209) * qJD(3) - t349 + t350) * t206 - (t130 * t209 + t132 * t205) * t311;
t396 = 0.2e1 * pkin(1);
t343 = qJ(4) * t205;
t236 = -t377 * t209 - t343;
t126 = pkin(2) - t236;
t32 = -t377 * t401 + t403;
t157 = t401 * qJ(4);
t73 = t123 * t205 + t152 * t209;
t48 = pkin(8) * t130 + t73;
t41 = t157 + t48;
t12 = t204 * t32 + t208 * t41;
t193 = t210 * qJDD(1);
t389 = -t206 * t302 + t193;
t128 = qJDD(3) - t389;
t112 = t389 * pkin(6) + qJDD(2) * pkin(7);
t308 = qJD(3) * t209;
t140 = t270 * qJD(2);
t80 = qJD(1) * t140 + qJDD(1) * t246;
t281 = t205 * t112 + t123 * t310 + t152 * t308 - t209 * t80;
t260 = qJDD(4) + t281;
t7 = pkin(8) * t61 - t377 * t128 + t260;
t121 = t128 * qJ(4);
t155 = t401 * qJD(4);
t23 = t209 * t112 + t123 * t308 - t152 * t310 + t205 * t80;
t15 = t121 + t155 + t23;
t9 = pkin(8) * t62 + t15;
t2 = -t12 * qJD(5) - t204 * t9 + t208 * t7;
t394 = t12 * t159 + t2;
t304 = qJD(5) * t208;
t305 = qJD(5) * t204;
t13 = -t130 * t304 + t132 * t305 - t204 * t62 + t208 * t61;
t393 = -t159 * t250 + t13;
t390 = qJD(4) * t205 + t190;
t247 = t204 * t209 - t334;
t134 = t204 * t205 + t208 * t209;
t386 = t134 * qJD(5) - t204 * t310 - t208 * t308;
t369 = pkin(7) * t128;
t58 = pkin(3) * t130 - t237;
t385 = -t401 * t58 + t369;
t232 = t128 * t209 - t310 * t401;
t383 = qJD(1) * (-t130 * t206 - t333 * t401) - t232;
t233 = t128 * t205 + t308 * t401;
t380 = t206 * (-qJD(2) * t130 - t233) + t210 * (-t313 * t401 + t62);
t105 = t134 * t206;
t238 = -t204 * t7 - t208 * t9 - t32 * t304 + t41 * t305;
t251 = t115 * t204 + t116 * t208;
t60 = t117 * t204 + t118 * t208;
t379 = g(1) * t60 + g(2) * t251 + g(3) * t105 + t250 * t43 + t238;
t378 = t132 ^ 2;
t214 = qJD(1) ^ 2;
t373 = pkin(1) * t214;
t372 = pkin(3) * t128;
t371 = pkin(3) * t205;
t368 = pkin(7) * t132;
t367 = g(1) * t207;
t363 = g(3) * t210;
t153 = t376 * t205;
t83 = t153 * t208 - t154 * t204;
t361 = qJD(5) * t83 + t407 * t204 - t408 * t208;
t84 = t153 * t204 + t154 * t208;
t360 = -qJD(5) * t84 + t408 * t204 + t407 * t208;
t230 = t210 * t134;
t359 = qJD(1) * t230 + t386;
t358 = t204 * t308 + t205 * t304 - t208 * t310 - t209 * t305 - t247 * t314;
t300 = t377 * t205;
t342 = qJ(4) * t209;
t235 = -t300 + t342;
t357 = t401 * t235 + t390;
t11 = -t204 * t41 + t208 * t32;
t356 = t11 * t159;
t53 = t157 + t73;
t352 = t401 * t53;
t351 = t401 * t73;
t141 = -qJ(4) * t204 - t208 * t377;
t348 = qJD(5) * t141 - t204 * t48 + t403 * t208;
t142 = qJ(4) * t208 - t204 * t377;
t347 = -qJD(5) * t142 - t403 * t204 - t208 * t48;
t258 = -t342 + t371;
t346 = t401 * t258 - t390;
t345 = pkin(6) * qJDD(1);
t344 = qJ(4) * t130;
t339 = t132 * t130;
t335 = t205 * t206;
t332 = t206 * t207;
t330 = t206 * t211;
t319 = t198 + t195;
t144 = -pkin(1) - t319;
t322 = t205 * t140 + t144 * t308;
t321 = (g(1) * t326 + g(2) * t329) * t206;
t173 = pkin(6) * t327;
t92 = t205 * t144 + t173;
t320 = g(1) * t332 - g(2) * t330;
t318 = t211 * pkin(1) + t207 * pkin(6);
t202 = t206 ^ 2;
t203 = t210 ^ 2;
t317 = t202 - t203;
t316 = t202 + t203;
t312 = qJD(2) * t206;
t306 = qJD(4) * t209;
t299 = t206 * t214 * t210;
t298 = g(1) * t325 + g(2) * t328 + g(3) * t206;
t292 = t205 * t311;
t291 = t210 * t303;
t290 = t205 * t309;
t289 = t401 * t315;
t286 = t130 ^ 2 - t378;
t280 = -t115 * pkin(3) + qJ(4) * t116;
t279 = -t117 * pkin(3) + qJ(4) * t118;
t171 = pkin(6) * t333;
t91 = t144 * t209 - t171;
t277 = t159 ^ 2;
t276 = -pkin(7) * t349 - t298;
t274 = pkin(3) * t327 + qJ(4) * t333 + t319;
t273 = pkin(2) * t325 + pkin(7) * t330 + t318;
t272 = t206 * t284;
t271 = t296 * t206;
t269 = -g(1) * t115 + g(2) * t117;
t268 = g(1) * t116 - g(2) * t118;
t266 = -g(2) * t211 + t367;
t199 = t211 * pkin(6);
t263 = -t116 * pkin(3) - qJ(4) * t115 + t199;
t85 = -qJ(4) * t210 + t92;
t262 = qJD(3) * t173 - t140 * t209 + t144 * t310;
t261 = (qJD(3) * t130 - t61) * pkin(7);
t259 = pkin(3) * t209 + t343;
t257 = pkin(6) * t130 + t151 * t205;
t256 = pkin(6) * t132 + t151 * t209;
t197 = t210 * pkin(3);
t63 = pkin(4) * t210 + t171 + t197 + (-pkin(8) * t206 - t144) * t209;
t71 = pkin(8) * t335 + t85;
t27 = -t204 * t71 + t208 * t63;
t28 = t204 * t63 + t208 * t71;
t51 = -pkin(3) * t401 + t323;
t255 = -t205 * t53 + t209 * t51;
t253 = -t205 * t73 - t209 * t72;
t252 = qJD(3) * t151 - t369;
t213 = qJD(2) ^ 2;
t243 = qJDD(2) * t210 - t206 * t213;
t88 = -pkin(6) * t293 + t119;
t242 = pkin(2) + t259;
t241 = -g(1) * t117 - g(2) * t115 - g(3) * t335;
t240 = -pkin(7) * qJD(3) * t401 - t363;
t16 = pkin(3) * t62 - t219;
t231 = -t16 + t240;
t229 = t264 + t240;
t226 = t246 * t367;
t225 = t118 * pkin(3) + qJ(4) * t117 + t273;
t221 = -t241 - t281;
t220 = t205 * t340 - t349;
t45 = (-t206 * t303 - t210 * t310) * pkin(6) + t322;
t218 = t132 * t58 + qJDD(4) - t221;
t217 = t62 * t335 + (t206 * t308 + t292) * t130;
t216 = g(1) * t118 + g(2) * t116 + g(3) * t331 + t401 * t72 - t23;
t185 = qJ(4) * t312;
t176 = pkin(7) * t325;
t172 = pkin(7) * t328;
t166 = qJ(4) * t331;
t122 = -qJDD(5) + t128;
t99 = -t166 + (pkin(6) + t371) * t206;
t87 = pkin(6) * t294 + t337;
t86 = t197 - t91;
t82 = t166 + (-pkin(6) - t300) * t206;
t81 = -t128 * t210 + t312 * t401;
t79 = pkin(3) * t132 + t344;
t78 = qJD(1) * t271 - t337;
t77 = t184 + t88;
t49 = -t377 * t132 - t344;
t46 = t312 * t370 - t262;
t44 = (qJD(3) * t259 - t306) * t206 + (pkin(6) + t258) * t311;
t42 = qJD(2) * t271 + t262;
t40 = -t61 + t340;
t39 = -qJD(4) * t210 + t185 + t45;
t36 = (-t132 * t206 - t327 * t401) * qJD(1) + t233;
t35 = qJD(2) * t230 + (qJD(3) - qJD(5)) * t206 * t247;
t34 = t204 * t291 + t386 * t206 - t208 * t292;
t33 = (qJD(3) * t236 + t306) * t206 + (-pkin(6) + t235) * t311;
t30 = t185 + (-pkin(6) * qJD(2) + pkin(8) * qJD(3)) * t331 + (-qJD(4) + (-pkin(6) * qJD(3) + pkin(8) * qJD(2)) * t205) * t210 + t322;
t29 = pkin(8) * t290 + qJD(2) * t222 + t262;
t26 = t209 * t338 - t350;
t25 = -t61 * t331 + (-t290 + t291) * t132;
t18 = t260 - t372;
t17 = (t303 * t401 + t61) * t210 + (qJD(2) * t132 + t232) * t206;
t4 = -qJD(5) * t28 - t204 * t30 + t208 * t29;
t3 = qJD(5) * t27 + t204 * t29 + t208 * t30;
t1 = [0, 0, 0, 0, 0, qJDD(1), t266, t267, 0, 0, qJDD(1) * t202 + 0.2e1 * t272, 0.2e1 * t206 * t193 - 0.2e1 * t317 * t302, qJDD(2) * t206 + t210 * t213, qJDD(1) * t203 - 0.2e1 * t272, t243, 0, (-0.2e1 * pkin(1) * t302 - pkin(6) * qJDD(2)) * t206 + (-pkin(6) * t213 + qJDD(1) * t396 + t266) * t210, -t243 * pkin(6) + t402 * t396 - t320, 0.2e1 * t316 * t345 - t267, -g(1) * (-pkin(1) * t207 + t199) - g(2) * t318 + (t316 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), t25, t397, t17, t217, t380, t81, t128 * t91 + t401 * t46 + (qJD(2) * t257 + t281) * t210 + (pkin(6) * t62 + qJD(2) * t72 + t151 * t308 - t205 * t264) * t206 + t268, -t128 * t92 - t401 * t45 + (qJD(2) * t256 + t23) * t210 + (-pkin(6) * t61 - qJD(2) * t73 - t151 * t310 - t209 * t264) * t206 + t269, -t130 * t45 - t132 * t46 + t61 * t91 - t62 * t92 + t253 * t311 + (-t205 * t23 + t209 * t281 + (t205 * t72 - t209 * t73) * qJD(3)) * t206 + t320, t23 * t92 + t73 * t45 - t281 * t91 + t72 * t46 - g(1) * t199 - g(2) * t273 - t226 + (t151 * t311 - t206 * t264) * pkin(6), t25, t17, -t397, t81, -t380, t217, -t128 * t86 + t130 * t44 - t401 * t42 + t62 * t99 + (t313 * t58 + t18) * t210 + (-qJD(2) * t51 + t16 * t205 + t308 * t58) * t206 + t268, -t130 * t39 + t132 * t42 - t61 * t86 - t62 * t85 + t255 * t311 + (-t15 * t205 + t18 * t209 + (-t205 * t51 - t209 * t53) * qJD(3)) * t206 + t320, t128 * t85 - t132 * t44 + t401 * t39 + t61 * t99 + (-t303 * t58 - t15) * t210 + (qJD(2) * t53 - t16 * t209 + t310 * t58) * t206 - t269, -g(1) * t263 - g(2) * t225 + t15 * t85 + t16 * t99 + t18 * t86 + t53 * t39 + t51 * t42 + t58 * t44 - t226, -t105 * t13 + t35 * t68, t104 * t13 - t105 * t14 - t250 * t35 - t34 * t68, -t105 * t122 - t13 * t210 + t159 * t35 - t312 * t68, t104 * t14 + t250 * t34, t104 * t122 - t14 * t210 - t159 * t34 + t250 * t312, -t122 * t210 - t159 * t312, g(1) * t251 - g(2) * t60 + t10 * t104 - t11 * t312 - t27 * t122 + t82 * t14 + t4 * t159 + t2 * t210 + t250 * t33 + t43 * t34, g(1) * t388 - g(2) * t59 + t10 * t105 + t12 * t312 + t28 * t122 - t82 * t13 - t3 * t159 + t210 * t238 + t33 * t68 + t43 * t35, t104 * t238 - t105 * t2 - t11 * t35 - t12 * t34 + t13 * t27 - t14 * t28 - t250 * t3 - t4 * t68 - t320, -t238 * t28 + t12 * t3 + t2 * t27 + t11 * t4 + t10 * t82 + t43 * t33 - g(1) * (-pkin(4) * t116 + t263) - g(2) * (pkin(4) * t118 - pkin(8) * t330 + t225) - (-t376 * t206 + t297) * t367; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, t317 * t214, t301, t299, t193, qJDD(2), -t363 - t188 + (t267 + t373) * t206, (-t345 + t373) * t210 + t298, 0, 0, t26, -t409, t36, t220, -t383, -t289, -pkin(2) * t62 - t401 * t87 + t252 * t205 + t229 * t209 + (-t206 * t72 - t210 * t257) * qJD(1) + t321, pkin(2) * t61 + t401 * t88 + t252 * t209 + (t206 * t73 - t210 * t256) * qJD(1) + (-t229 - t384) * t205, t130 * t88 + t132 * t87 + (t72 * t314 + t23 + (-t72 + t368) * qJD(3)) * t209 + (t281 + t261 - t351) * t205 + t276, t264 * pkin(2) - t73 * t88 - t72 * t87 - t151 * t190 - g(1) * (-pkin(2) * t330 + t176) - g(2) * (-pkin(2) * t332 + t172) - g(3) * t319 + (qJD(3) * t253 + t205 * t281 + t23 * t209) * pkin(7), t26, t36, t409, -t289, t383, t220, t346 * t130 - t385 * t205 + t231 * t209 - t242 * t62 + t51 * t315 + t401 * t78 + t321, t130 * t77 - t132 * t78 + (-t51 * t314 + t15 + (t51 + t368) * qJD(3)) * t209 + (t18 + t261 - t352) * t205 + t276, -t53 * t315 - t242 * t61 - t401 * t77 - t346 * t132 + t385 * t209 + (t231 + t384) * t205, -t53 * t77 - t51 * t78 - g(1) * t176 - g(2) * t172 - g(3) * t274 + t346 * t58 + (qJD(3) * t255 + t15 * t209 + t18 * t205) * pkin(7) + (-t16 + t384) * t242, t13 * t247 - t359 * t68, t13 * t134 + t14 * t247 + t250 * t359 - t358 * t68, t122 * t247 - t159 * t359 + t315 * t68, t134 * t14 + t250 * t358, t122 * t134 - t159 * t358 - t250 * t315, t159 * t315, -g(3) * t230 + t11 * t315 - t83 * t122 + t126 * t14 + t411 * t134 + t360 * t159 + t357 * t250 + t358 * t43, -t12 * t315 + t84 * t122 - t126 * t13 - t361 * t159 + t357 * t68 - t359 * t43 + (t363 - t411) * t247, t11 * t359 - t12 * t358 + t13 * t83 + t134 * t238 - t14 * t84 + t2 * t247 - t250 * t361 - t360 * t68 + t298, -t238 * t84 + t2 * t83 + t10 * t126 - g(1) * (-pkin(8) * t325 + t176) - g(2) * (-pkin(8) * t328 + t172) - g(3) * (pkin(4) * t327 + t274) + t357 * t43 + t361 * t12 + t360 * t11 + (g(3) * pkin(8) + t267 * t126) * t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t339, -t286, t40, -t339, -t404, t128, -t132 * t151 + t221 + t351, t130 * t151 + t216, 0, 0, t339, t40, t286, t128, t404, -t339, -t130 * t79 - t218 + t351 + 0.2e1 * t372, pkin(3) * t61 - qJ(4) * t62 + (t53 - t73) * t132 + (t51 - t323) * t130, -t130 * t58 + t132 * t79 + 0.2e1 * t121 + 0.2e1 * t155 - t216, t15 * qJ(4) - t18 * pkin(3) - t58 * t79 - t51 * t73 - g(1) * t279 - g(2) * t280 - g(3) * (-pkin(3) * t335 + t166) + t323 * t53, -t362, -t399, t393, t362, -t392, t122, -t141 * t122 + t159 * t347 - t250 * t49 - t2 - t410, t142 * t122 - t159 * t348 - t49 * t68 - t379, t13 * t141 - t14 * t142 + (t11 - t348) * t250 + (-t12 - t347) * t68, -t238 * t142 + t2 * t141 - t43 * t49 - g(1) * (-pkin(4) * t117 + t279) - g(2) * (-pkin(4) * t115 + t280) - g(3) * (-t206 * t300 + t166) + t348 * t12 + t347 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128 + t339, t40, -t401 ^ 2 - t378, t218 - t352 - t372, 0, 0, 0, 0, 0, 0, -t122 * t208 - t132 * t250 - t204 * t277, t122 * t204 - t132 * t68 - t208 * t277, t392 * t204 + t208 * t393, -t132 * t43 + t394 * t208 + (-t238 - t356) * t204 + t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t362, t399, -t393, -t362, t392, -t122, t394 + t410, t356 + t379, 0, 0;];
tau_reg = t1;

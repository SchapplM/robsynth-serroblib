% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:14
% EndTime: 2019-12-31 22:02:27
% DurationCPUTime: 6.11s
% Computational Cost: add. (7440->558), mult. (16983->702), div. (0->0), fcn. (11647->10), ass. (0->272)
t235 = sin(qJ(2));
t238 = cos(qJ(2));
t274 = pkin(2) * t235 - pkin(7) * t238;
t170 = t274 * qJD(1);
t237 = cos(qJ(3));
t234 = sin(qJ(3));
t316 = qJD(1) * t235;
t290 = t234 * t316;
t117 = pkin(6) * t290 + t237 * t170;
t324 = t237 * t238;
t261 = pkin(3) * t235 - pkin(8) * t324;
t375 = pkin(8) + pkin(7);
t298 = qJD(3) * t375;
t402 = -qJD(1) * t261 - t237 * t298 - t117;
t315 = qJD(1) * t238;
t207 = -qJD(3) + t315;
t195 = -qJD(4) + t207;
t312 = qJD(2) * t237;
t162 = -t290 + t312;
t295 = t237 * t316;
t314 = qJD(2) * t234;
t163 = t295 + t314;
t233 = sin(qJ(4));
t374 = cos(qJ(4));
t97 = -t374 * t162 + t163 * t233;
t347 = t195 * t97;
t309 = qJD(3) * t235;
t285 = qJD(1) * t309;
t305 = qJD(1) * qJD(2);
t286 = t238 * t305;
t303 = t235 * qJDD(1);
t391 = qJD(2) * qJD(3) + t286 + t303;
t280 = t391 * t234 + t237 * t285;
t254 = t237 * qJDD(2) - t280;
t288 = t374 * qJD(4);
t307 = qJD(4) * t233;
t89 = (-qJDD(2) + t285) * t234 - t391 * t237;
t36 = -t162 * t288 + t163 * t307 - t233 * t254 + t374 * t89;
t401 = -t36 - t347;
t258 = t233 * t162 + t374 * t163;
t376 = t258 ^ 2;
t95 = t97 ^ 2;
t400 = -t95 + t376;
t147 = t234 * t170;
t328 = t235 * t237;
t331 = t234 * t238;
t399 = t147 + (-pkin(6) * t328 - pkin(8) * t331) * qJD(1) + t234 * t298;
t217 = pkin(6) * t303;
t342 = qJDD(2) * pkin(2);
t138 = pkin(6) * t286 + t217 - t342;
t236 = sin(qJ(1));
t239 = cos(qJ(1));
t271 = g(1) * t239 + g(2) * t236;
t359 = g(3) * t238;
t250 = -t235 * t271 + t359;
t247 = -t138 - t250;
t398 = -qJD(3) * pkin(7) * t207 - t247;
t397 = pkin(3) * t307;
t396 = t374 * pkin(3);
t395 = t97 * qJ(5);
t394 = t97 * t258;
t296 = t234 * t315;
t297 = t374 * t237;
t333 = t233 * t234;
t382 = qJD(3) + qJD(4);
t384 = t374 * qJD(3) + t288;
t322 = -t233 * t296 - t237 * t384 + t297 * t315 + t333 * t382;
t165 = t233 * t237 + t374 * t234;
t112 = t382 * t165;
t321 = -t165 * t315 + t112;
t341 = t258 * t195;
t37 = qJD(4) * t258 - t233 * t89 - t374 * t254;
t393 = -t37 - t341;
t311 = qJD(2) * t238;
t294 = t234 * t311;
t308 = qJD(3) * t237;
t392 = t235 * t308 + t294;
t184 = -qJD(2) * pkin(2) + pkin(6) * t316;
t119 = -pkin(3) * t162 + t184;
t232 = qJ(3) + qJ(4);
t223 = sin(t232);
t224 = cos(t232);
t325 = t236 * t238;
t129 = t223 * t239 - t224 * t325;
t323 = t238 * t239;
t131 = t223 * t236 + t224 * t323;
t360 = g(3) * t235;
t222 = t238 * qJDD(1);
t385 = -t235 * t305 + t222;
t159 = qJDD(3) - t385;
t173 = t274 * qJD(2);
t275 = pkin(2) * t238 + pkin(7) * t235;
t179 = -pkin(1) - t275;
t113 = qJD(1) * t173 + qJDD(1) * t179;
t103 = t237 * t113;
t154 = t179 * qJD(1);
t219 = pkin(6) * t315;
t185 = qJD(2) * pkin(7) + t219;
t108 = t154 * t234 + t185 * t237;
t137 = pkin(6) * t385 + qJDD(2) * pkin(7);
t49 = -qJD(3) * t108 - t137 * t234 + t103;
t27 = pkin(3) * t159 + pkin(8) * t89 + t49;
t310 = qJD(3) * t234;
t48 = t234 * t113 + t237 * t137 + t154 * t308 - t185 * t310;
t34 = pkin(8) * t254 + t48;
t107 = t237 * t154 - t185 * t234;
t77 = -pkin(8) * t163 + t107;
t69 = -pkin(3) * t207 + t77;
t78 = pkin(8) * t162 + t108;
t5 = t233 * t27 + t69 * t288 - t78 * t307 + t374 * t34;
t248 = g(1) * t131 - g(2) * t129 + t224 * t360 - t5;
t390 = t119 * t97 + t248;
t186 = t375 * t234;
t187 = t375 * t237;
t353 = -t186 * t288 - t187 * t307 + t402 * t233 - t399 * t374;
t116 = -t233 * t186 + t374 * t187;
t352 = -qJD(4) * t116 + t399 * t233 + t374 * t402;
t388 = t107 * t207 + t48;
t387 = qJ(5) * t258;
t276 = -t219 + (-t296 + t310) * pkin(3);
t153 = qJDD(4) + t159;
t281 = pkin(3) * t288;
t372 = pkin(3) * t233;
t386 = -t153 * t372 + t195 * t281;
t209 = pkin(6) * t324;
t122 = t234 * t179 + t209;
t142 = t234 * t325 + t237 * t239;
t144 = -t234 * t323 + t236 * t237;
t383 = -g(1) * t144 + g(2) * t142;
t139 = t153 * pkin(4);
t346 = t36 * qJ(5);
t381 = -t258 * qJD(5) + t139 + t346;
t128 = t223 * t325 + t224 * t239;
t130 = -t223 * t323 + t224 * t236;
t380 = -g(1) * t130 + g(2) * t128 + t223 * t360;
t283 = -t233 * t34 + t374 * t27;
t74 = t374 * t78;
t39 = t233 * t69 + t74;
t6 = -qJD(4) * t39 + t283;
t244 = t6 + t380;
t379 = -t119 * t258 + t244;
t284 = -pkin(4) * t97 - qJD(5);
t63 = t119 - t284;
t378 = -t63 * t258 + t380;
t377 = -0.2e1 * pkin(1);
t373 = pkin(3) * t163;
t371 = pkin(3) * t234;
t370 = pkin(3) * t237;
t369 = pkin(6) * t234;
t366 = g(1) * t236;
t361 = g(2) * t239;
t72 = t233 * t78;
t38 = t374 * t69 - t72;
t25 = t38 - t387;
t22 = -pkin(4) * t195 + t25;
t357 = -t25 + t22;
t356 = -t97 * t281 - t37 * t372;
t164 = -t297 + t333;
t355 = -t321 * qJ(5) - qJD(5) * t164 + t353;
t354 = -pkin(4) * t316 + t322 * qJ(5) - t165 * qJD(5) + t352;
t46 = t374 * t77 - t72;
t351 = t321 * pkin(4) + t276;
t350 = qJ(5) * t37;
t345 = t89 * t234;
t161 = t237 * t179;
t105 = -pkin(8) * t328 + t161 + (-pkin(3) - t369) * t238;
t332 = t234 * t235;
t114 = -pkin(8) * t332 + t122;
t58 = t233 * t105 + t374 * t114;
t344 = pkin(6) * qJDD(1);
t343 = qJD(5) * t97;
t339 = t108 * t207;
t338 = t162 * t207;
t337 = t162 * t234;
t336 = t163 * t162;
t335 = t163 * t207;
t334 = t163 * t237;
t330 = t234 * t239;
t329 = t235 * t236;
t327 = t235 * t239;
t326 = t235 * t375;
t313 = qJD(2) * t235;
t320 = t237 * t173 + t313 * t369;
t174 = pkin(3) * t332 + t235 * pkin(6);
t319 = t239 * pkin(1) + t236 * pkin(6);
t230 = t235 ^ 2;
t231 = t238 ^ 2;
t318 = t230 - t231;
t317 = t230 + t231;
t300 = t233 * t332;
t242 = qJD(1) ^ 2;
t299 = t235 * t242 * t238;
t120 = t392 * pkin(3) + pkin(6) * t311;
t216 = pkin(2) + t370;
t293 = t234 * t309;
t291 = t195 * t316;
t178 = pkin(4) * t224 + t370;
t45 = -t233 * t77 - t74;
t57 = t374 * t105 - t114 * t233;
t115 = -t374 * t186 - t187 * t233;
t279 = t374 * t311;
t278 = t235 * t286;
t210 = g(1) * t329;
t277 = -g(2) * t327 + t210;
t273 = -g(1) * t128 - g(2) * t130;
t272 = -g(1) * t129 - g(2) * t131;
t270 = pkin(6) * t162 - t184 * t234;
t269 = pkin(6) * t163 + t184 * t237;
t268 = -pkin(7) * t159 + qJD(3) * t184;
t267 = -t107 * t237 - t108 * t234;
t169 = pkin(2) + t178;
t229 = -qJ(5) - t375;
t266 = t169 * t238 - t229 * t235;
t264 = t216 * t238 + t326;
t262 = (g(1) * t327 + g(2) * t329 - t359) * t224;
t260 = t271 * t223;
t259 = -pkin(6) * qJDD(2) + t305 * t377;
t257 = t159 * t234 - t207 * t308;
t256 = t159 * t237 + t207 * t310;
t56 = t261 * qJD(2) + (-t209 + (pkin(8) * t235 - t179) * t234) * qJD(3) + t320;
t75 = t234 * t173 + t179 * t308 + (-t235 * t312 - t238 * t310) * pkin(6);
t60 = -pkin(8) * t392 + t75;
t15 = t105 * t288 - t114 * t307 + t233 * t56 + t374 * t60;
t255 = pkin(1) * t242 + t271;
t241 = qJD(2) ^ 2;
t253 = pkin(6) * t241 + qJDD(1) * t377 + t361;
t251 = t254 * t237;
t249 = -t238 * t271 - t360;
t246 = t248 + t350;
t16 = -t58 * qJD(4) - t233 * t60 + t374 * t56;
t68 = -pkin(3) * t254 + t138;
t17 = t37 * pkin(4) + qJDD(5) + t68;
t227 = t239 * pkin(6);
t215 = pkin(4) + t396;
t193 = t223 * t359;
t177 = pkin(4) * t223 + t371;
t145 = t234 * t236 + t237 * t323;
t143 = -t236 * t324 + t330;
t135 = t235 * t297 - t300;
t134 = t165 * t235;
t125 = pkin(4) * t164 - t216;
t121 = -pkin(6) * t331 + t161;
t118 = -pkin(6) * t295 + t147;
t109 = pkin(4) * t134 + t174;
t106 = -t153 * t238 - t195 * t313;
t82 = -qJ(5) * t164 + t116;
t81 = -qJ(5) * t165 + t115;
t76 = -qJD(3) * t122 + t320;
t71 = pkin(4) * t258 + t373;
t62 = t234 * t279 - t233 * t293 - qJD(4) * t300 + (t233 * t311 + t235 * t384) * t237;
t61 = t112 * t235 + t233 * t294 - t237 * t279;
t51 = pkin(4) * t62 + t120;
t50 = -qJ(5) * t134 + t58;
t47 = -pkin(4) * t238 - qJ(5) * t135 + t57;
t31 = t46 - t387;
t30 = t45 + t395;
t29 = -t153 * t164 + t321 * t195 + t97 * t316;
t28 = t153 * t165 + t322 * t195 - t258 * t316;
t26 = t39 - t395;
t14 = t134 * t37 + t62 * t97;
t13 = -t135 * t36 - t258 * t61;
t12 = t37 * t164 + t321 * t97;
t11 = -t36 * t165 - t258 * t322;
t10 = -t134 * t153 + t195 * t62 + t238 * t37 - t97 * t313;
t9 = t135 * t153 + t195 * t61 + t238 * t36 + t258 * t313;
t8 = -qJ(5) * t62 - qJD(5) * t134 + t15;
t7 = pkin(4) * t313 + t61 * qJ(5) - t135 * qJD(5) + t16;
t4 = t134 * t36 - t135 * t37 - t258 * t62 + t61 * t97;
t3 = t164 * t36 - t165 * t37 - t258 * t321 + t322 * t97;
t2 = -t343 + t5 - t350;
t1 = t6 + t381;
t18 = [0, 0, 0, 0, 0, qJDD(1), -t361 + t366, t271, 0, 0, qJDD(1) * t230 + 0.2e1 * t278, 0.2e1 * t235 * t222 - 0.2e1 * t318 * t305, qJDD(2) * t235 + t238 * t241, qJDD(1) * t231 - 0.2e1 * t278, qJDD(2) * t238 - t235 * t241, 0, t259 * t235 + (-t253 + t366) * t238, t235 * t253 + t238 * t259 - t210, 0.2e1 * t317 * t344 - t271, -g(1) * (-pkin(1) * t236 + t227) - g(2) * t319 + (t317 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), -t89 * t328 + (t237 * t311 - t293) * t163, (t162 * t237 - t163 * t234) * t311 + (t251 + t345 + (-t334 - t337) * qJD(3)) * t235, (-t207 * t312 + t89) * t238 + (qJD(2) * t163 + t256) * t235, -t162 * t392 - t254 * t332, (t207 * t314 - t254) * t238 + (qJD(2) * t162 - t257) * t235, -t159 * t238 - t207 * t313, -g(1) * t143 - g(2) * t145 + t121 * t159 - t76 * t207 + (-qJD(2) * t270 - t49) * t238 + (-pkin(6) * t254 + t107 * qJD(2) + t138 * t234 + t184 * t308) * t235, -g(1) * t142 - g(2) * t144 - t122 * t159 + t207 * t75 + (qJD(2) * t269 + t48) * t238 + (-pkin(6) * t89 - qJD(2) * t108 + t138 * t237 - t184 * t310) * t235, -t76 * t163 + t121 * t89 + t75 * t162 + t122 * t254 + t210 + t267 * t311 + (-t361 - t48 * t234 - t49 * t237 + (t107 * t234 - t108 * t237) * qJD(3)) * t235, t48 * t122 + t108 * t75 + t49 * t121 + t107 * t76 - g(1) * t227 - g(2) * (t239 * t275 + t319) - t179 * t366 + (t138 * t235 + t184 * t311) * pkin(6), t13, t4, t9, t14, t10, t106, t119 * t62 + t120 * t97 + t134 * t68 + t153 * t57 - t16 * t195 + t174 * t37 - t238 * t6 + t313 * t38 + t272, -t119 * t61 + t120 * t258 + t135 * t68 + t15 * t195 - t153 * t58 - t174 * t36 + t238 * t5 - t313 * t39 + t273, -t134 * t5 - t135 * t6 - t15 * t97 - t16 * t258 + t36 * t57 - t37 * t58 + t38 * t61 - t39 * t62 + t277, t5 * t58 + t39 * t15 + t6 * t57 + t38 * t16 + t68 * t174 + t119 * t120 - g(1) * (pkin(3) * t330 + t227) - g(2) * (t216 * t323 + t239 * t326 + t319) + (-g(1) * (-pkin(1) - t264) - g(2) * t371) * t236, t13, t4, t9, t14, t10, t106, -t1 * t238 + t109 * t37 + t134 * t17 + t153 * t47 - t195 * t7 + t22 * t313 + t51 * t97 + t62 * t63 + t272, -t109 * t36 + t135 * t17 - t153 * t50 + t195 * t8 + t2 * t238 + t258 * t51 - t26 * t313 - t61 * t63 + t273, -t1 * t135 - t134 * t2 + t22 * t61 - t258 * t7 - t26 * t62 + t36 * t47 - t37 * t50 - t8 * t97 + t277, t2 * t50 + t26 * t8 + t1 * t47 + t22 * t7 + t17 * t109 + t63 * t51 - g(1) * (t177 * t239 + t227) - g(2) * (t169 * t323 - t229 * t327 + t319) + (-g(1) * (-pkin(1) - t266) - g(2) * t177) * t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, t318 * t242, t303, t299, t222, qJDD(2), t235 * t255 - t217 - t359, t360 + (t255 - t344) * t238, 0, 0, -t207 * t334 - t345, (-t89 - t338) * t237 + (t254 + t335) * t234, (-t163 * t235 + t207 * t324) * qJD(1) + t257, t207 * t337 + t251, (-t162 * t235 - t207 * t331) * qJD(1) + t256, t207 * t316, -pkin(2) * t280 + t117 * t207 + t268 * t234 + (-t107 * t235 + t238 * t270) * qJD(1) + (t342 - t398) * t237, pkin(2) * t89 - t118 * t207 + t268 * t237 + (t108 * t235 - t238 * t269) * qJD(1) + t398 * t234, t117 * t163 - t118 * t162 + ((qJD(3) * t163 + t254) * pkin(7) + t388) * t237 + (-t49 + t339 + (-qJD(3) * t162 - t89) * pkin(7)) * t234 + t249, -t184 * t219 - t107 * t117 - t108 * t118 + t247 * pkin(2) + (qJD(3) * t267 - t49 * t234 + t48 * t237 + t249) * pkin(7), t11, t3, t28, t12, t29, t291, t115 * t153 + t321 * t119 + t164 * t68 - t352 * t195 - t216 * t37 + t276 * t97 - t38 * t316 + t262, -t116 * t153 + t165 * t68 + t216 * t36 + t193 + t353 * t195 - t322 * t119 + t276 * t258 + (qJD(1) * t39 - t260) * t235, t115 * t36 - t116 * t37 - t164 * t5 - t165 * t6 - t258 * t352 - t321 * t39 + t322 * t38 - t353 * t97 + t249, -g(3) * t264 + t6 * t115 + t5 * t116 + t119 * t276 - t68 * t216 + t352 * t38 + t353 * t39 + t271 * (t216 * t235 - t238 * t375), t11, t3, t28, t12, t29, t291, t125 * t37 + t153 * t81 + t164 * t17 - t195 * t354 - t22 * t316 + t321 * t63 + t351 * t97 + t262, -t125 * t36 - t153 * t82 + t165 * t17 + t193 - t322 * t63 + t355 * t195 + t351 * t258 + (qJD(1) * t26 - t260) * t235, -t1 * t165 - t164 * t2 + t22 * t322 - t258 * t354 - t26 * t321 - t355 * t97 + t36 * t81 - t37 * t82 + t249, -g(3) * t266 + t1 * t81 + t17 * t125 + t2 * t82 + t22 * t354 + t26 * t355 + t351 * t63 + t271 * (t169 * t235 + t229 * t238); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t336, -t162 ^ 2 + t163 ^ 2, -t89 + t338, t336, t254 - t335, t159, -t185 * t308 - t339 - t163 * t184 + t103 + (-qJD(3) * t154 - t137 + t360) * t234 + t383, g(1) * t145 - g(2) * t143 + g(3) * t328 - t162 * t184 - t388, 0, 0, t394, t400, t401, -t394, t393, t153, t45 * t195 + (t374 * t153 - t163 * t97 + t195 * t307) * pkin(3) + t379, -t195 * t46 - t258 * t373 + t386 + t390, t36 * t396 + t356 + (t39 + t45 + t397) * t258 + (-t38 + t46) * t97, -t38 * t45 - t39 * t46 + (t5 * t233 + t6 * t374 - t119 * t163 + g(3) * t332 + (-t38 * t233 + t39 * t374) * qJD(4) + t383) * pkin(3), t394, t400, t401, -t394, t393, t153, t215 * t153 + t30 * t195 - t71 * t97 + (-t74 + (pkin(3) * t195 - t69) * t233) * qJD(4) + t283 + t378 + t381, -t195 * t31 - t258 * t71 + t63 * t97 + t246 + t343 + t386, t215 * t36 + t356 + (t26 + t30 + t397) * t258 + (-t22 + t31) * t97, t1 * t215 - t26 * t31 - t22 * t30 - t63 * t71 - g(1) * (-t177 * t323 + t178 * t236) - g(2) * (-t177 * t325 - t178 * t239) + t177 * t360 + (t2 * t233 + (-t22 * t233 + t374 * t26) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t394, t400, t401, -t394, t393, t153, -t39 * t195 + t379, -t195 * t38 + t390, 0, 0, t394, t400, t401, -t394, t393, t153, t346 - t26 * t195 + 0.2e1 * t139 + (t284 - t63) * t258 + t244, -pkin(4) * t376 - t195 * t25 + (qJD(5) + t63) * t97 + t246, t36 * pkin(4) - t357 * t97, t357 * t26 + (t1 + t378) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 - t341, -t36 + t347, -t95 - t376, t22 * t258 + t26 * t97 + t17 + t250;];
tau_reg = t18;

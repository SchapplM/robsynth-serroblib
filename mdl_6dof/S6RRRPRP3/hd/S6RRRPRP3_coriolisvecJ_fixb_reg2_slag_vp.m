% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:41:32
% EndTime: 2019-03-09 16:41:48
% DurationCPUTime: 6.46s
% Computational Cost: add. (15000->519), mult. (36662->648), div. (0->0), fcn. (27016->8), ass. (0->270)
t269 = sin(qJ(3));
t270 = sin(qJ(2));
t271 = cos(qJ(2));
t397 = cos(qJ(3));
t235 = t269 * t271 + t270 * t397;
t222 = qJD(1) * t235;
t263 = qJD(2) + qJD(3);
t266 = sin(pkin(10));
t267 = cos(pkin(10));
t196 = t222 * t267 + t263 * t266;
t268 = sin(qJ(5));
t396 = cos(qJ(5));
t205 = t266 * t222;
t404 = t263 * t267 - t205;
t297 = t396 * t404;
t138 = -t268 * t196 + t297;
t420 = t138 ^ 2;
t332 = t397 * t271;
t251 = qJD(1) * t332;
t342 = qJD(1) * t270;
t330 = t269 * t342;
t220 = -t251 + t330;
t215 = qJD(5) + t220;
t419 = t138 * t215;
t411 = t196 * t396 + t268 * t404;
t400 = t411 ^ 2;
t352 = t269 * t270;
t310 = t263 * t352;
t344 = t263 * t251;
t177 = qJD(1) * t310 - t344;
t331 = t396 * t267;
t369 = t177 * t266;
t69 = t268 * (qJD(5) * t196 - t369) - qJD(5) * t297 + t177 * t331;
t418 = t69 - t419;
t417 = t215 * t411;
t233 = t266 * t396 + t268 * t267;
t219 = t233 * qJD(5);
t415 = t233 * t220 + t219;
t293 = -t268 * t266 + t331;
t326 = qJD(5) * t396;
t340 = qJD(5) * t268;
t402 = -t266 * t340 + t267 * t326;
t414 = t293 * t220 + t402;
t398 = -pkin(8) - pkin(7);
t246 = t398 * t271;
t239 = qJD(1) * t246;
t224 = t397 * t239;
t245 = t398 * t270;
t237 = qJD(1) * t245;
t186 = t237 * t269 - t224;
t341 = qJD(3) * t269;
t405 = -pkin(2) * t341 + t186;
t413 = t196 * t220 - t369;
t355 = t267 * t177;
t412 = t220 * t404 - t355;
t391 = qJD(2) * pkin(2);
t336 = t270 * t391;
t410 = 0.2e1 * t336;
t339 = qJD(1) * qJD(2);
t409 = -0.2e1 * t339;
t179 = pkin(3) * t222 + qJ(4) * t220;
t162 = pkin(2) * t342 + t179;
t223 = t269 * t239;
t187 = t237 * t397 + t223;
t120 = t266 * t162 + t267 * t187;
t359 = t220 * t266;
t338 = pkin(9) * t359;
t102 = t338 + t120;
t255 = pkin(2) * t269 + qJ(4);
t226 = (-pkin(9) - t255) * t266;
t260 = t267 * pkin(9);
t227 = t255 * t267 + t260;
t174 = t268 * t226 + t227 * t396;
t327 = qJD(3) * t397;
t249 = pkin(2) * t327 + qJD(4);
t119 = t267 * t162 - t187 * t266;
t358 = t220 * t267;
t313 = t222 * pkin(4) + pkin(9) * t358;
t87 = t119 + t313;
t384 = qJD(5) * t174 - t268 * t102 + t233 * t249 + t396 * t87;
t228 = t237 + t391;
t183 = t228 * t397 + t223;
t123 = t266 * t179 + t267 * t183;
t105 = t338 + t123;
t242 = (-pkin(9) - qJ(4)) * t266;
t243 = qJ(4) * t267 + t260;
t194 = t268 * t242 + t243 * t396;
t122 = t267 * t179 - t183 * t266;
t90 = t122 + t313;
t378 = qJD(4) * t233 + qJD(5) * t194 - t268 * t105 + t396 * t90;
t191 = t263 * t235;
t178 = t191 * qJD(1);
t325 = t270 * t339;
t101 = pkin(2) * t325 + pkin(3) * t178 + qJ(4) * t177 - qJD(4) * t222;
t333 = qJD(2) * t398;
t316 = qJD(1) * t333;
t229 = t270 * t316;
t230 = t271 * t316;
t319 = -t228 * t327 - t397 * t229 - t269 * t230 - t239 * t341;
t121 = qJD(4) * t263 - t319;
t61 = t267 * t101 - t121 * t266;
t62 = t266 * t101 + t267 * t121;
t307 = -t61 * t266 + t62 * t267;
t70 = qJD(5) * t411 - t177 * t233;
t408 = t70 - t417;
t202 = pkin(4) * t359;
t312 = t202 - t405;
t234 = -t332 + t352;
t147 = t178 * t234;
t406 = t191 * t220 + t147;
t403 = t397 * t245 + t269 * t246;
t401 = t415 * pkin(5) - t414 * qJ(6) - qJD(6) * t233;
t47 = -t138 * t222 + t178 * t293 - t215 * t415;
t259 = -pkin(2) * t271 - pkin(1);
t182 = pkin(3) * t234 - qJ(4) * t235 + t259;
t199 = t269 * t245 - t246 * t397;
t130 = t267 * t182 - t199 * t266;
t103 = pkin(4) * t234 - t235 * t260 + t130;
t131 = t266 * t182 + t267 * t199;
t357 = t235 * t266;
t116 = -pkin(9) * t357 + t131;
t388 = t268 * t103 + t396 * t116;
t190 = -qJD(2) * t332 - t271 * t327 + t310;
t366 = t190 * t267;
t117 = pkin(3) * t191 + qJ(4) * t190 - qJD(4) * t235 + t336;
t238 = t270 * t333;
t240 = t271 * t333;
t139 = qJD(3) * t403 + t397 * t238 + t269 * t240;
t67 = t267 * t117 - t139 * t266;
t51 = pkin(4) * t191 + pkin(9) * t366 + t67;
t367 = t190 * t266;
t68 = t266 * t117 + t267 * t139;
t63 = pkin(9) * t367 + t68;
t14 = -qJD(5) * t388 - t268 * t63 + t396 * t51;
t10 = t138 * t414 - t233 * t70 - t293 * t69 - t411 * t415;
t399 = t220 ^ 2;
t395 = pkin(5) * t178;
t393 = t222 * pkin(5);
t392 = t267 * pkin(4);
t44 = t396 * t102 + t268 * t87;
t161 = -t263 * pkin(3) + qJD(4) - t183;
t127 = -pkin(4) * t404 + t161;
t56 = -pkin(5) * t138 - qJ(6) * t411 + t127;
t390 = t411 * t56;
t50 = t396 * t105 + t268 * t90;
t295 = t226 * t396 - t268 * t227;
t128 = qJD(5) * t295 + t249 * t293;
t212 = t222 * qJ(6);
t35 = t212 + t44;
t387 = t128 - t35;
t386 = t128 - t44;
t385 = t393 + t384;
t383 = t401 + t312;
t184 = t269 * t228 - t224;
t144 = -t202 + t184;
t382 = -t144 + t401;
t294 = t242 * t396 - t268 * t243;
t151 = qJD(4) * t293 + qJD(5) * t294;
t37 = t212 + t50;
t381 = t151 - t37;
t380 = t151 - t50;
t379 = t393 + t378;
t126 = t228 * t341 + t269 * t229 - t397 * t230 - t239 * t327;
t377 = t126 * t403;
t376 = t411 * t138;
t372 = t295 * t178;
t371 = t174 * t178;
t370 = t177 * t235;
t368 = t178 * qJ(6);
t364 = t294 * t178;
t363 = t194 * t178;
t362 = t196 * t266;
t361 = t215 * t222;
t360 = t220 * t222;
t273 = qJD(1) ^ 2;
t351 = t271 * t273;
t272 = qJD(2) ^ 2;
t350 = t272 * t270;
t349 = t272 * t271;
t244 = qJD(1) * t259;
t157 = t220 * pkin(3) - t222 * qJ(4) + t244;
t165 = t263 * qJ(4) + t184;
t107 = t267 * t157 - t165 * t266;
t75 = pkin(4) * t220 - pkin(9) * t196 + t107;
t108 = t266 * t157 + t267 * t165;
t83 = pkin(9) * t404 + t108;
t33 = -t268 * t83 + t396 * t75;
t348 = qJD(6) - t33;
t343 = t270 ^ 2 - t271 ^ 2;
t337 = t397 * pkin(2);
t334 = t270 * t351;
t256 = -pkin(3) - t392;
t328 = -t400 + t420;
t32 = pkin(4) * t178 + pkin(9) * t355 + t61;
t41 = pkin(9) * t369 + t62;
t324 = t268 * t41 - t396 * t32 + t83 * t326 + t75 * t340;
t323 = pkin(1) * t409;
t322 = t108 * t222 + t126 * t266;
t321 = t249 * t266 + t119;
t258 = -t337 - pkin(3);
t318 = t128 * t138 - t174 * t70 + t295 * t69;
t317 = t138 * t151 - t194 * t70 + t294 * t69;
t314 = t271 * t325;
t159 = pkin(4) * t357 - t403;
t168 = t233 * t235;
t89 = -t190 * t233 + t235 * t402;
t309 = -t138 * t89 + t168 * t70;
t91 = -pkin(4) * t369 + t126;
t306 = -t107 * t358 - t108 * t359 + t307;
t305 = t267 * t404;
t304 = -t107 * t222 - t126 * t267;
t303 = t126 * t235 + t177 * t403;
t302 = -t400 - t420;
t4 = t324 - t395;
t34 = t268 * t75 + t396 * t83;
t299 = t215 * t34 - t324;
t7 = t268 * t32 + t75 * t326 - t340 * t83 + t396 * t41;
t57 = t103 * t396 - t268 * t116;
t13 = t103 * t326 - t116 * t340 + t268 * t51 + t396 * t63;
t292 = -t138 * t415 - t293 * t70;
t17 = pkin(5) * t70 + qJ(6) * t69 - qJD(6) * t411 + t91;
t26 = -t215 * pkin(5) + t348;
t291 = -t17 * t293 + t26 * t222 + t415 * t56;
t27 = t215 * qJ(6) + t34;
t290 = -t17 * t233 - t27 * t222 - t414 * t56;
t289 = -t244 * t222 - t126;
t288 = t244 * t220 + t319;
t287 = t414 * t127 + t34 * t222 + t91 * t233;
t286 = t415 * t127 - t33 * t222 - t91 * t293;
t285 = -t161 * t190 + t303;
t284 = t177 * t234 - t178 * t235 + t190 * t220;
t176 = -pkin(5) * t293 - t233 * qJ(6) + t256;
t2 = qJD(6) * t215 + t368 + t7;
t283 = t2 * t293 + t4 * t233 + t414 * t26 - t415 * t27;
t282 = t324 * t233 + t7 * t293 - t414 * t33 - t415 * t34;
t281 = t305 - t362;
t280 = pkin(3) * t177 - qJ(4) * t178 + (-qJD(4) + t161) * t220;
t279 = -t258 * t177 - t255 * t178 + (t161 - t249) * t220;
t169 = t293 * t235;
t88 = t190 * t293 + t219 * t235;
t277 = -t138 * t88 + t168 * t69 - t169 * t70 - t411 * t89;
t276 = -t138 * t191 + t168 * t178 + t215 * t89 + t234 * t70;
t140 = qJD(3) * t199 + t269 * t238 - t397 * t240;
t106 = -pkin(4) * t367 + t140;
t274 = t70 + t417;
t262 = t267 ^ 2;
t261 = t266 ^ 2;
t241 = t258 - t392;
t163 = -t337 + t176;
t148 = t222 ^ 2 - t399;
t142 = t344 + (t220 - t330) * t263;
t110 = t412 * t266;
t109 = t413 * t267;
t104 = t191 * t215 + t147;
t85 = t267 * t178 - t222 * t404 - t266 * t399;
t84 = t178 * t266 - t196 * t222 + t267 * t399;
t81 = pkin(5) * t168 - qJ(6) * t169 + t159;
t79 = pkin(5) * t411 - qJ(6) * t138;
t66 = t281 * t220 + (t261 - t262) * t177;
t55 = -t234 * pkin(5) - t57;
t54 = qJ(6) * t234 + t388;
t46 = t233 * t178 + t215 * t414 - t222 * t411;
t42 = -t69 - t419;
t21 = -t69 * t233 + t411 * t414;
t20 = -t169 * t69 - t411 * t88;
t19 = t89 * pkin(5) + t88 * qJ(6) - t169 * qJD(6) + t106;
t18 = t169 * t178 + t191 * t411 - t215 * t88 - t234 * t69;
t12 = -t191 * pkin(5) - t14;
t11 = qJ(6) * t191 + qJD(6) * t234 + t13;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t314, t343 * t409, t349, -0.2e1 * t314, -t350, 0, -pkin(7) * t349 + t270 * t323, pkin(7) * t350 + t271 * t323, 0, 0, -t190 * t222 - t370, -t191 * t222 + t284, -t190 * t263, t406, -t191 * t263, 0, -t140 * t263 + t259 * t178 + t244 * t191 + (qJD(1) * t234 + t220) * t336, -t139 * t263 - t259 * t177 - t244 * t190 + t222 * t410, -t139 * t220 + t140 * t222 - t178 * t199 + t183 * t190 - t184 * t191 + t234 * t319 + t303, t184 * t139 - t183 * t140 - t199 * t319 + t244 * t410 - t377, -t196 * t366 - t262 * t370, -t190 * t281 + 0.2e1 * t355 * t357, t196 * t191 - t267 * t284, -t261 * t370 + t367 * t404, t191 * t404 + t266 * t284, t406, t107 * t191 + t130 * t178 - t140 * t404 + t67 * t220 + t61 * t234 + t266 * t285, -t108 * t191 - t131 * t178 + t140 * t196 - t68 * t220 - t62 * t234 + t267 * t285, -t67 * t196 - t68 * t205 + (t107 * t190 + t130 * t177 - t235 * t61 + t263 * t68) * t267 + (t108 * t190 + t131 * t177 - t235 * t62) * t266, t107 * t67 + t108 * t68 + t130 * t61 + t131 * t62 + t140 * t161 - t377, t20, t277, t18, t309, -t276, t104, -t106 * t138 + t127 * t89 + t14 * t215 + t159 * t70 + t168 * t91 + t178 * t57 + t191 * t33 - t234 * t324, t106 * t411 - t127 * t88 - t13 * t215 - t159 * t69 + t169 * t91 - t178 * t388 - t191 * t34 - t234 * t7, t13 * t138 - t14 * t411 - t168 * t7 + t169 * t324 + t33 * t88 - t34 * t89 - t388 * t70 + t57 * t69, t106 * t127 + t13 * t34 + t14 * t33 + t159 * t91 - t324 * t57 + t388 * t7, t20, t18, -t277, t104, t276, t309, -t12 * t215 - t138 * t19 + t168 * t17 - t178 * t55 - t191 * t26 - t234 * t4 + t56 * t89 + t70 * t81, t11 * t138 + t12 * t411 - t168 * t2 + t169 * t4 - t26 * t88 - t27 * t89 - t54 * t70 - t55 * t69, t11 * t215 - t169 * t17 + t178 * t54 - t19 * t411 + t191 * t27 + t2 * t234 + t56 * t88 + t69 * t81, t11 * t27 + t12 * t26 + t17 * t81 + t19 * t56 + t2 * t54 + t4 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t334, t343 * t273, 0, t334, 0, 0, t273 * pkin(1) * t270, pkin(1) * t351, 0, 0, t360, t148, t142, -t360, 0, 0, t186 * t263 + (-t220 * t342 - t263 * t341) * pkin(2) + t289, t187 * t263 + (-t222 * t342 - t263 * t327) * pkin(2) + t288 (t184 - t186) * t222 + (-t183 + t187) * t220 + (t397 * t177 - t178 * t269 + (-t220 * t397 + t222 * t269) * qJD(3)) * pkin(2), t183 * t186 - t184 * t187 + (-t244 * t342 - t397 * t126 - t319 * t269 + (-t183 * t269 + t184 * t397) * qJD(3)) * pkin(2), t109, t66, t84, -t110, t85, -t360, -t119 * t220 + t266 * t279 + t404 * t405 + t304, t120 * t220 - t196 * t405 + t267 * t279 + t322, -t120 * t404 + t196 * t321 + t249 * t305 + t306, t126 * t258 + t307 * t255 - t405 * t161 + (t249 * t267 - t120) * t108 - t321 * t107, t21, t10, t46, t292, t47, -t361, -t138 * t312 - t215 * t384 + t241 * t70 + t286 + t372, -t215 * t386 - t241 * t69 + t312 * t411 + t287 - t371, -t138 * t44 + t384 * t411 + t282 + t318, t127 * t312 + t7 * t174 + t91 * t241 - t295 * t324 - t33 * t384 + t34 * t386, t21, t46, -t10, -t361, -t47, t292, -t138 * t383 + t163 * t70 - t215 * t385 + t291 + t372, -t138 * t35 + t385 * t411 + t283 + t318, t163 * t69 + t215 * t387 - t383 * t411 + t290 + t371, t17 * t163 + t2 * t174 + t26 * t385 + t27 * t387 - t295 * t4 + t383 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t360, t148, t142, -t360, 0, 0, t184 * t263 + t289, t183 * t263 + t288, 0, 0, t109, t66, t84, -t110, t85, -t360, -t122 * t220 + t184 * t404 + t266 * t280 + t304, t123 * t220 - t184 * t196 + t267 * t280 + t322, -t123 * t404 + t122 * t196 + (t305 + t362) * qJD(4) + t306, -t126 * pkin(3) - t107 * t122 - t108 * t123 - t161 * t184 + (-t107 * t266 + t108 * t267) * qJD(4) + t307 * qJ(4), t21, t10, t46, t292, t47, -t361, t138 * t144 - t215 * t378 + t256 * t70 + t286 + t364, -t144 * t411 - t215 * t380 - t256 * t69 + t287 - t363, -t138 * t50 + t378 * t411 + t282 + t317, -t127 * t144 + t7 * t194 + t91 * t256 - t294 * t324 - t33 * t378 + t34 * t380, t21, t46, -t10, -t361, -t47, t292, -t138 * t382 + t176 * t70 - t215 * t379 + t291 + t364, -t138 * t37 + t379 * t411 + t283 + t317, t176 * t69 + t215 * t381 - t382 * t411 + t290 + t363, t17 * t176 + t2 * t194 + t26 * t379 + t27 * t381 - t294 * t4 + t382 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t413, t412, -t196 ^ 2 - t404 ^ 2, t107 * t196 - t108 * t404 + t126, 0, 0, 0, 0, 0, 0, t274, -t418, t302, -t138 * t34 + t33 * t411 + t91, 0, 0, 0, 0, 0, 0, t274, t302, t418, -t138 * t27 - t26 * t411 + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t376, -t328, t42, t376, -t408, t178, -t127 * t411 + t299, -t127 * t138 + t215 * t33 - t7, 0, 0, -t376, t42, t328, t178, t408, t376, t138 * t79 + t299 - t390 + 0.2e1 * t395, pkin(5) * t69 - t70 * qJ(6) + (t27 - t34) * t411 - (t26 - t348) * t138, 0.2e1 * t368 + t56 * t138 + t79 * t411 + (0.2e1 * qJD(6) - t33) * t215 + t7, -t4 * pkin(5) + t2 * qJ(6) - t26 * t34 + t27 * t348 - t56 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t222 * t263 - t376, t42, -t215 ^ 2 - t400, -t215 * t27 + t390 + t4;];
tauc_reg  = t1;
% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% tauc_reg [6x35]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRR11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:34:15
% EndTime: 2019-03-09 19:34:33
% DurationCPUTime: 7.76s
% Computational Cost: add. (8809->533), mult. (22606->729), div. (0->0), fcn. (17731->10), ass. (0->248)
t232 = sin(pkin(6));
t241 = cos(qJ(2));
t323 = qJD(1) * t241;
t301 = t232 * t323;
t202 = -qJD(3) + t301;
t198 = qJD(5) + t202;
t237 = sin(qJ(2));
t313 = qJD(1) * qJD(2);
t297 = t232 * t313;
t206 = t237 * t297;
t234 = sin(qJ(6));
t238 = cos(qJ(6));
t314 = qJD(6) * t238;
t240 = cos(qJ(3));
t279 = t241 * t297;
t236 = sin(qJ(3));
t233 = cos(pkin(6));
t324 = qJD(1) * t233;
t285 = qJD(2) + t324;
t325 = qJD(1) * t232;
t302 = t237 * t325;
t152 = t236 * t302 - t240 * t285;
t391 = qJD(3) * t152;
t111 = -t240 * t279 + t391;
t154 = t236 * t285 + t240 * t302;
t320 = qJD(3) * t154;
t112 = t236 * t279 + t320;
t235 = sin(qJ(5));
t239 = cos(qJ(5));
t316 = qJD(5) * t239;
t317 = qJD(5) * t235;
t41 = -t239 * t111 + t235 * t112 + t152 * t316 - t154 * t317;
t98 = t152 * t235 + t154 * t239;
t12 = t198 * t314 + t238 * t41 + (-qJD(6) * t98 - t206) * t234;
t77 = t198 * t234 + t238 * t98;
t13 = qJD(6) * t77 + t206 * t238 + t234 * t41;
t372 = -t239 * t152 + t154 * t235;
t383 = qJD(6) + t372;
t360 = t77 * t383;
t75 = -t238 * t198 + t234 * t98;
t361 = t75 * t383;
t401 = (t13 + t360) * t234 - (t12 - t361) * t238;
t318 = qJD(3) * t240;
t398 = -t240 * t301 + t318;
t281 = t236 * t301;
t319 = qJD(3) * t236;
t392 = t281 - t319;
t356 = t12 * t234;
t393 = t238 * t383;
t397 = t393 * t77 + t356;
t290 = -t111 * t235 - t239 * t112;
t42 = qJD(5) * t98 + t290;
t353 = t234 * t42;
t396 = t383 * t393 - t77 * t98 + t353;
t363 = pkin(9) - pkin(10);
t309 = pkin(1) * t324;
t172 = -pkin(8) * t302 + t241 * t309;
t260 = (pkin(2) * t237 - pkin(9) * t241) * t232;
t173 = qJD(1) * t260;
t329 = t240 * t172 + t236 * t173;
t90 = qJ(4) * t302 + t329;
t395 = pkin(10) * t281 + t363 * t319 + t90;
t205 = t363 * t240;
t288 = -t236 * t172 + t173 * t240;
t336 = t240 * t241;
t364 = pkin(3) + pkin(4);
t394 = qJD(3) * t205 - (-pkin(10) * t336 - t364 * t237) * t325 + t288;
t175 = pkin(8) * t301 + t237 * t309;
t128 = pkin(9) * t285 + t175;
t167 = (-pkin(2) * t241 - pkin(9) * t237 - pkin(1)) * t232;
t142 = qJD(1) * t167;
t83 = -t236 * t128 + t240 * t142;
t335 = qJD(4) - t83;
t390 = pkin(5) * t98;
t385 = -pkin(10) * t154 + t335;
t53 = t364 * t202 + t385;
t196 = t202 * qJ(4);
t84 = t240 * t128 + t236 * t142;
t64 = pkin(10) * t152 + t84;
t59 = -t196 + t64;
t17 = t235 * t53 + t239 * t59;
t15 = pkin(11) * t198 + t17;
t127 = -t285 * pkin(2) - t172;
t68 = t152 * pkin(3) - t154 * qJ(4) + t127;
t58 = -pkin(4) * t152 - t68;
t18 = pkin(5) * t372 - pkin(11) * t98 + t58;
t9 = t15 * t238 + t18 * t234;
t389 = t9 * t98;
t388 = t383 * t98;
t275 = t15 * t234 - t18 * t238;
t387 = t275 * t98;
t386 = t98 * t372;
t187 = t235 * t236 + t239 * t240;
t117 = (qJD(3) - qJD(5)) * t187;
t130 = t187 * t301;
t332 = t117 - t130;
t331 = t398 * t235 + t236 * t316 + t392 * t239 - t240 * t317;
t338 = t232 * t241;
t312 = pkin(8) * t338;
t362 = pkin(1) * t237;
t370 = (t233 * t362 + t312) * qJD(2);
t164 = qJD(1) * t370;
t45 = t112 * pkin(3) + t111 * qJ(4) - t154 * qJD(4) + t164;
t384 = t398 * qJ(4) + t236 * qJD(4) + t175;
t339 = t232 * t237;
t182 = t233 * t236 + t240 * t339;
t299 = qJD(2) * t338;
t121 = qJD(3) * t182 + t236 * t299;
t181 = -t233 * t240 + t236 * t339;
t122 = -qJD(3) * t181 + t240 * t299;
t56 = t121 * pkin(3) - t122 * qJ(4) - t182 * qJD(4) + t370;
t382 = t372 ^ 2 - t98 ^ 2;
t194 = t202 * qJD(4);
t201 = qJ(4) * t206;
t174 = qJD(2) * t260;
t162 = qJD(1) * t174;
t337 = t233 * t241;
t371 = pkin(1) * t337 - pkin(8) * t339;
t176 = t371 * qJD(2);
t163 = qJD(1) * t176;
t254 = t128 * t319 - t142 * t318 - t236 * t162 - t240 * t163;
t43 = -t194 + t201 - t254;
t28 = pkin(10) * t112 + t43;
t322 = qJD(2) * t237;
t300 = t232 * t322;
t267 = t364 * t300;
t283 = t128 * t318 + t142 * t319 - t240 * t162 + t236 * t163;
t29 = pkin(10) * t111 - qJD(1) * t267 + t283;
t295 = t235 * t28 - t239 * t29 + t59 * t316 + t53 * t317;
t381 = -t58 * t98 - t295;
t351 = t238 * t42;
t380 = t75 * t98 + t351;
t378 = t198 * t372 + t41;
t258 = t235 * t29 + t239 * t28 + t53 * t316 - t317 * t59;
t377 = t372 * t58 - t258;
t350 = t392 * t364 + t384;
t204 = t363 * t236;
t133 = t204 * t235 + t205 * t239;
t376 = -qJD(5) * t133 + t395 * t235 + t394 * t239;
t268 = pkin(3) * t206;
t46 = -t268 + t283;
t74 = -t196 + t84;
t375 = t202 * t74 + t46;
t264 = t204 * t239 - t205 * t235;
t374 = -qJD(5) * t264 - t394 * t235 + t395 * t239;
t373 = t237 * t241;
t327 = t239 * qJ(4) - t235 * t364;
t190 = -pkin(11) + t327;
t4 = pkin(5) * t206 + t295;
t99 = t154 * pkin(3) + t152 * qJ(4);
t73 = -pkin(4) * t154 - t99;
t369 = (-pkin(11) * t372 + qJD(6) * t190 - t390 + t73) * t383 - t4;
t368 = (t383 * pkin(11) + t390) * t383 + t4;
t200 = -t240 * pkin(3) - t236 * qJ(4) - pkin(2);
t184 = t240 * pkin(4) - t200;
t188 = -t235 * t240 + t236 * t239;
t100 = pkin(5) * t187 - pkin(11) * t188 + t184;
t367 = (-t331 * pkin(5) + t332 * pkin(11) + qJD(6) * t133 - t350) * t383 - t100 * t42;
t166 = t312 + (pkin(9) + t362) * t233;
t289 = -t236 * t166 + t167 * t240;
t87 = pkin(3) * t338 - t289;
t65 = pkin(4) * t338 - pkin(10) * t182 + t87;
t330 = t240 * t166 + t236 * t167;
t86 = -qJ(4) * t338 + t330;
t70 = pkin(10) * t181 + t86;
t270 = t235 * t65 + t239 * t70;
t261 = -t166 * t318 - t167 * t319 + t174 * t240 - t236 * t176;
t38 = -pkin(10) * t122 - t261 - t267;
t253 = -t166 * t319 + t167 * t318 + t236 * t174 + t240 * t176;
t49 = qJ(4) * t300 - qJD(4) * t338 + t253;
t39 = pkin(10) * t121 + t49;
t366 = -qJD(5) * t270 - t235 * t39 + t239 * t38;
t32 = -pkin(4) * t112 - t45;
t10 = pkin(5) * t42 - pkin(11) * t41 + t32;
t3 = -pkin(11) * t206 + t258;
t1 = -qJD(6) * t275 + t10 * t234 + t238 * t3;
t365 = t154 ^ 2;
t243 = qJD(1) ^ 2;
t358 = -pkin(5) * t302 - t376;
t355 = t154 * t68;
t352 = t234 * t383;
t266 = -qJ(4) * t235 - t239 * t364;
t349 = -qJD(5) * t266 + t235 * t64 - t385 * t239;
t348 = t327 * qJD(5) + t385 * t235 + t239 * t64;
t347 = t152 * t202;
t346 = t154 * t152;
t345 = t154 * t202;
t343 = t188 * t238;
t342 = t202 * t236;
t341 = t202 * t240;
t229 = t232 ^ 2;
t340 = t229 * t243;
t333 = t392 * pkin(3) + t384;
t326 = t237 ^ 2 - t241 ^ 2;
t321 = qJD(2) * t240;
t315 = qJD(6) * t234;
t311 = pkin(9) * t342;
t310 = pkin(9) * t341;
t308 = pkin(9) * t322;
t307 = pkin(9) * t321;
t165 = -t233 * pkin(2) - t371;
t303 = t229 * t323;
t298 = t229 * t313;
t292 = t383 * t198;
t287 = t198 ^ 2;
t286 = 0.2e1 * t298;
t284 = qJD(2) + 0.2e1 * t324;
t277 = -0.2e1 * pkin(1) * t298;
t31 = pkin(11) * t338 + t270;
t108 = t181 * t235 + t182 * t239;
t265 = t181 * t239 - t182 * t235;
t85 = t181 * pkin(3) - t182 * qJ(4) + t165;
t69 = -pkin(4) * t181 - t85;
t33 = -pkin(5) * t265 - pkin(11) * t108 + t69;
t274 = t234 * t33 + t238 * t31;
t273 = -t234 * t31 + t238 * t33;
t16 = -t235 * t59 + t239 * t53;
t271 = -t235 * t70 + t239 * t65;
t259 = -t108 * t234 + t238 * t338;
t89 = t108 * t238 + t234 * t338;
t257 = t235 * t38 + t239 * t39 + t65 * t316 - t317 * t70;
t256 = -t202 * t84 - t283;
t101 = t130 * t234 + t238 * t302;
t252 = t234 * t117 + t188 * t314 - t101;
t102 = t130 * t238 - t234 * t302;
t251 = t238 * t117 - t188 * t315 - t102;
t250 = pkin(1) * (-t233 * t313 + t340);
t14 = -pkin(5) * t198 - t16;
t247 = -pkin(11) * t42 + (t14 + t16) * t383;
t246 = -t202 * t83 + t254;
t2 = -qJD(6) * t9 + t238 * t10 - t234 * t3;
t245 = -t190 * t42 + (-t14 + t349) * t383;
t244 = -t133 * t42 + t4 * t188 + (-pkin(11) * t302 - qJD(6) * t100 + t374) * t383;
t44 = -pkin(4) * t121 - t56;
t189 = pkin(5) - t266;
t93 = -pkin(3) * t302 - t288;
t78 = -t111 - t347;
t72 = pkin(3) * t202 + t335;
t55 = -pkin(3) * t300 - t261;
t52 = qJD(5) * t265 + t121 * t235 + t122 * t239;
t51 = qJD(5) * t108 - t121 * t239 + t122 * t235;
t30 = -pkin(5) * t338 - t271;
t20 = qJD(6) * t259 - t234 * t300 + t238 * t52;
t19 = qJD(6) * t89 + t234 * t52 + t238 * t300;
t11 = pkin(5) * t51 - pkin(11) * t52 + t44;
t6 = pkin(5) * t300 - t366;
t5 = -pkin(11) * t300 + t257;
t7 = [0, 0, 0, t286 * t373, -t326 * t286, t284 * t299, -t284 * t300, 0, -t164 * t233 + t237 * t277 - t285 * t370, -t163 * t233 - t176 * t285 + t241 * t277, -t111 * t182 + t122 * t154, t111 * t181 - t112 * t182 - t121 * t154 - t122 * t152, -t122 * t202 + (t111 * t241 + (qJD(1) * t182 + t154) * t322) * t232, t121 * t202 + (t112 * t241 + (-qJD(1) * t181 - t152) * t322) * t232 (-t202 * t232 - t303) * t322, -t261 * t202 + t370 * t152 + t165 * t112 + t164 * t181 + t127 * t121 + (t283 * t241 + (qJD(1) * t289 + t83) * t322) * t232, t253 * t202 + t370 * t154 - t165 * t111 + t164 * t182 + t127 * t122 + (-t254 * t241 + (-qJD(1) * t330 - t84) * t322) * t232, t112 * t85 + t121 * t68 + t152 * t56 + t181 * t45 + t202 * t55 + (t241 * t46 + (-qJD(1) * t87 - t72) * t322) * t232, -t111 * t87 - t112 * t86 - t121 * t74 + t122 * t72 - t152 * t49 + t154 * t55 - t181 * t43 + t182 * t46, t111 * t85 - t122 * t68 - t154 * t56 - t182 * t45 - t202 * t49 + (-t241 * t43 + (qJD(1) * t86 + t74) * t322) * t232, t43 * t86 + t45 * t85 + t46 * t87 + t49 * t74 + t55 * t72 + t56 * t68, t108 * t41 + t52 * t98, -t108 * t42 + t265 * t41 - t372 * t52 - t51 * t98, t198 * t52 + (t241 * t41 + (-qJD(1) * t108 - t98) * t322) * t232, -t198 * t51 + (-t241 * t42 + (-qJD(1) * t265 + t372) * t322) * t232 (-t198 * t232 - t303) * t322, t366 * t198 + t44 * t372 + t69 * t42 - t32 * t265 + t58 * t51 + (-t295 * t241 + (-qJD(1) * t271 - t16) * t322) * t232, -t257 * t198 + t44 * t98 + t69 * t41 + t32 * t108 + t58 * t52 + (-t258 * t241 + (t270 * qJD(1) + t17) * t322) * t232, t12 * t89 + t20 * t77, t12 * t259 - t13 * t89 - t19 * t77 - t20 * t75, -t12 * t265 + t20 * t383 + t42 * t89 + t51 * t77, t13 * t265 - t19 * t383 + t259 * t42 - t51 * t75, -t265 * t42 + t383 * t51 (-qJD(6) * t274 + t11 * t238 - t234 * t5) * t383 + t273 * t42 - t2 * t265 - t275 * t51 + t6 * t75 + t30 * t13 - t4 * t259 + t14 * t19 -(qJD(6) * t273 + t11 * t234 + t238 * t5) * t383 - t274 * t42 + t1 * t265 - t9 * t51 + t6 * t77 + t30 * t12 + t4 * t89 + t14 * t20; 0, 0, 0, -t340 * t373, t326 * t340, -t232 * t243 * t337, t285 * t302 - t206, 0, -pkin(8) * t279 + t175 * t285 + t237 * t250, pkin(8) * t206 + t172 * t285 + t241 * t250, -t111 * t236 - t154 * t341 (-t111 + t347) * t240 + (-t112 + t345) * t236, -t202 * t318 + (t202 * t336 + (t236 * qJD(2) - t154) * t237) * t325, t202 * t319 + (-t241 * t342 + (t152 + t321) * t237) * t325, t202 * t302, -pkin(2) * t112 - t164 * t240 + t288 * t202 - t175 * t152 + (t127 * t236 + t310) * qJD(3) + (-t83 * t237 + (-t127 * t241 - t308) * t236) * t325, pkin(2) * t111 + t164 * t236 - t329 * t202 - t175 * t154 + (t127 * t240 - t311) * qJD(3) + (-t127 * t336 + (t84 - t307) * t237) * t325, t112 * t200 - t202 * t93 - t240 * t45 - t333 * t152 + (t236 * t68 + t310) * qJD(3) + (t237 * t72 + (-t241 * t68 - t308) * t236) * t325, t152 * t90 - t154 * t93 + (t43 - t202 * t72 + (-t112 + t320) * pkin(9)) * t240 + ((-t111 + t391) * pkin(9) + t375) * t236, t111 * t200 + t202 * t90 - t236 * t45 + t333 * t154 + (-t240 * t68 + t311) * qJD(3) + (t68 * t336 + (-t74 + t307) * t237) * t325, t200 * t45 - t72 * t93 - t74 * t90 - t333 * t68 + (t236 * t46 + t240 * t43 + (-t236 * t74 + t240 * t72) * qJD(3)) * pkin(9), t188 * t41 + t332 * t98, -t187 * t41 - t188 * t42 - t331 * t98 - t332 * t372, t332 * t198 + (-qJD(2) * t188 + t98) * t302, -t331 * t198 + (qJD(2) * t187 - t372) * t302, t198 * t302, t184 * t42 + t32 * t187 + t350 * t372 + t331 * t58 + t376 * t198 + (-qJD(2) * t264 + t16) * t302, t184 * t41 + t32 * t188 + t350 * t98 + t332 * t58 + t374 * t198 + (qJD(2) * t133 - t17) * t302, t12 * t343 + t251 * t77, t101 * t77 + t102 * t75 + (-t234 * t77 - t238 * t75) * t117 + (-t356 - t13 * t238 + (t234 * t75 - t238 * t77) * qJD(6)) * t188, t12 * t187 + t251 * t383 + t331 * t77 + t343 * t42, -t13 * t187 - t188 * t353 - t252 * t383 - t331 * t75, t187 * t42 + t331 * t383, -t13 * t264 + t252 * t14 + t2 * t187 + t244 * t234 - t238 * t367 - t275 * t331 + t358 * t75, -t1 * t187 - t12 * t264 + t251 * t14 + t234 * t367 + t244 * t238 - t331 * t9 + t358 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t346, -t152 ^ 2 + t365, t78, -t345 - t112, t206, -t127 * t154 + t256, t127 * t152 + t246, -t152 * t99 + t256 + 0.2e1 * t268 - t355, pkin(3) * t111 - qJ(4) * t112 + (t74 - t84) * t154 + (t72 - t335) * t152, -t152 * t68 + t154 * t99 - 0.2e1 * t194 + 0.2e1 * t201 - t246, -pkin(3) * t46 + qJ(4) * t43 + t335 * t74 - t68 * t99 - t72 * t84, -t386, t382, -t378, -t198 * t98 + t42, t206, -t348 * t198 - t266 * t206 - t372 * t73 - t381, t349 * t198 + t327 * t206 - t73 * t98 - t377, -t397, t401, -t396, t352 * t383 - t380, t388, t189 * t13 + t245 * t234 - t238 * t369 + t348 * t75 - t387, t189 * t12 + t234 * t369 + t245 * t238 + t348 * t77 - t389; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206 + t346, t78, -t202 ^ 2 - t365, t355 + t375, 0, 0, 0, 0, 0, -t154 * t372 - t206 * t239 - t235 * t287, -t154 * t98 + t206 * t235 - t239 * t287, 0, 0, 0, 0, 0, -t154 * t393 + (-t234 * t292 - t13) * t239 + (t198 * t75 - t314 * t383 - t353) * t235, t154 * t352 + (-t238 * t292 - t12) * t239 + (t198 * t77 + t315 * t383 - t351) * t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t386, -t382, t378, -t290 + (-qJD(5) + t198) * t98, -t206, t17 * t198 + t381, t16 * t198 + t377, t397, -t401, t396, -t234 * t383 ^ 2 + t380, -t388, -pkin(5) * t13 - t17 * t75 + t247 * t234 - t238 * t368 + t387, -pkin(5) * t12 - t17 * t77 + t234 * t368 + t247 * t238 + t389; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t75, -t75 ^ 2 + t77 ^ 2, t12 + t361, -t13 + t360, t42, -t14 * t77 + t383 * t9 + t2, t14 * t75 - t275 * t383 - t1;];
tauc_reg  = t7;
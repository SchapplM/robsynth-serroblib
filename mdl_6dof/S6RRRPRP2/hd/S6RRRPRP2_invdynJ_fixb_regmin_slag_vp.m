% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:37:21
% EndTime: 2019-03-09 16:37:37
% DurationCPUTime: 6.56s
% Computational Cost: add. (12869->551), mult. (30612->681), div. (0->0), fcn. (22955->14), ass. (0->303)
t237 = cos(qJ(5));
t329 = qJD(5) * t237;
t234 = sin(qJ(3));
t235 = sin(qJ(2));
t238 = cos(qJ(3));
t239 = cos(qJ(2));
t169 = t234 * t235 - t238 * t239;
t152 = t169 * qJD(1);
t334 = qJD(1) * t235;
t344 = t234 * t239;
t154 = -qJD(1) * t344 - t238 * t334;
t232 = sin(pkin(10));
t368 = cos(pkin(10));
t302 = -t368 * t152 + t154 * t232;
t357 = t302 * t237;
t413 = t329 - t357;
t233 = sin(qJ(5));
t330 = qJD(5) * t233;
t419 = -t302 * t233 + t330;
t226 = t239 * pkin(2);
t389 = pkin(1) + t226;
t292 = pkin(3) * t169 - t389;
t187 = qJD(1) * t389;
t418 = qJDD(1) * t389;
t406 = qJD(5) - t302;
t417 = t406 ^ 2;
t323 = qJD(2) + qJD(3);
t404 = pkin(7) + pkin(8);
t189 = t404 * t239;
t177 = qJD(1) * t189;
t159 = t238 * t177;
t188 = t404 * t235;
t175 = qJD(1) * t188;
t383 = qJD(2) * pkin(2);
t163 = -t175 + t383;
t277 = -t163 * t234 - t159;
t367 = qJ(4) * t152;
t104 = -t277 - t367;
t306 = t368 * t104;
t147 = t154 * qJ(4);
t155 = t234 * t177;
t301 = t238 * t163 - t155;
t103 = t147 + t301;
t96 = pkin(3) * t323 + t103;
t61 = t232 * t96 + t306;
t59 = pkin(9) * t323 + t61;
t131 = pkin(3) * t152 + qJD(4) - t187;
t266 = -t232 * t152 - t154 * t368;
t73 = -pkin(4) * t302 - pkin(9) * t266 + t131;
t34 = t233 * t73 + t237 * t59;
t26 = qJ(6) * t406 + t34;
t416 = t26 * t406;
t107 = t233 * t323 + t237 * t266;
t299 = t107 * t406;
t298 = t406 * t233;
t295 = qJD(1) * t323;
t279 = t235 * t295;
t324 = t239 * qJDD(1);
t327 = qJD(1) * qJD(2);
t310 = t239 * t327;
t325 = t235 * qJDD(1);
t265 = -t310 - t325;
t326 = qJD(1) * qJD(3);
t407 = -t239 * t326 + t265;
t293 = t234 * t324 - t238 * t407;
t255 = -t234 * t279 + t293;
t347 = t232 * t238;
t348 = t232 * t234;
t243 = t368 * t255 + (-t279 + t324) * t347 + t407 * t348;
t415 = qJD(5) * t323 + t243;
t231 = qJ(2) + qJ(3);
t223 = pkin(10) + t231;
t209 = sin(t223);
t240 = cos(qJ(1));
t353 = t209 * t240;
t236 = sin(qJ(1));
t396 = g(2) * t236;
t414 = g(1) * t353 + t209 * t396;
t286 = t237 * pkin(5) + t233 * qJ(6);
t289 = g(1) * t240 + t396;
t130 = qJDD(2) * pkin(2) + t265 * t404;
t127 = t238 * t130;
t311 = t235 * t327;
t264 = -t311 + t324;
t133 = t404 * t264;
t297 = qJD(3) * t163 + t133;
t332 = qJD(3) * t238;
t313 = t177 * t332;
t321 = qJDD(2) + qJDD(3);
t48 = -t313 + t127 - t293 * qJ(4) + t154 * qJD(4) + t321 * pkin(3) + (qJ(4) * t279 - t297) * t234;
t333 = qJD(3) * t234;
t150 = t177 * t333;
t53 = -t152 * qJD(4) - t150 + (qJ(4) * t407 + t130) * t234 + ((-t235 * t326 + t264) * qJ(4) + t297) * t238;
t16 = t232 * t48 + t368 * t53;
t14 = pkin(9) * t321 + t16;
t170 = t235 * t238 + t344;
t248 = t170 * t295;
t262 = t169 * qJDD(1);
t245 = -t262 - t248;
t207 = pkin(2) * t311;
t322 = qJDD(4) + t207;
t366 = qJDD(1) * pkin(1);
t71 = -t232 * t255 + t368 * t245;
t23 = -pkin(2) * t324 - pkin(3) * t245 - t71 * pkin(4) - pkin(9) * t243 + t322 - t366;
t269 = t237 * t14 + t233 * t23 + t73 * t329 - t330 * t59;
t68 = qJDD(5) - t71;
t384 = qJ(6) * t68;
t2 = qJD(6) * t406 + t269 + t384;
t304 = t233 * t14 - t237 * t23 + t59 * t329 + t73 * t330;
t403 = pkin(5) * t68;
t4 = qJDD(6) + t304 - t403;
t411 = t2 * t237 + t4 * t233;
t109 = t175 * t234 - t159 + t367;
t339 = -t238 * t175 - t155;
t110 = t147 + t339;
t305 = t368 * t234;
t385 = pkin(2) * qJD(3);
t371 = t368 * t109 - t110 * t232 + (t305 + t347) * t385;
t144 = (t238 * t368 - t348) * t385;
t79 = t232 * t109 + t110 * t368;
t370 = t144 - t79;
t410 = t419 * pkin(5) - qJ(6) * t413 - t233 * qJD(6);
t337 = -t234 * t188 + t238 * t189;
t210 = cos(t223);
t409 = -t210 * pkin(4) - t209 * pkin(9);
t408 = g(1) * t236 - g(2) * t240;
t405 = t107 ^ 2;
t402 = pkin(3) * t154;
t400 = pkin(3) * t232;
t399 = pkin(5) * t266;
t395 = g(3) * t209;
t394 = g(3) * t210;
t225 = cos(t231);
t393 = g(3) * t225;
t392 = g(3) * t233;
t15 = -t232 * t53 + t368 * t48;
t13 = -pkin(4) * t321 - t15;
t51 = -t233 * t321 - t237 * t415 + t266 * t330;
t52 = t233 * t415 - t237 * t321 + t266 * t329;
t5 = t52 * pkin(5) + t51 * qJ(6) - t107 * qJD(6) + t13;
t391 = t233 * t5;
t97 = t232 * t104;
t70 = t103 * t368 - t97;
t83 = pkin(4) * t266 - pkin(9) * t302 - t402;
t388 = t233 * t83 + t237 * t70;
t220 = pkin(2) * t334;
t80 = t220 + t83;
t387 = t233 * t80 + t237 * t79;
t125 = t169 * t368 + t170 * t232;
t126 = -t232 * t169 + t170 * t368;
t87 = pkin(4) * t125 - pkin(9) * t126 + t292;
t300 = -t238 * t188 - t189 * t234;
t117 = -qJ(4) * t170 + t300;
t118 = -qJ(4) * t169 + t337;
t89 = t232 * t117 + t118 * t368;
t386 = t233 * t87 + t237 * t89;
t382 = t406 * t34;
t105 = t233 * t266 - t237 * t323;
t60 = t368 * t96 - t97;
t58 = -pkin(4) * t323 - t60;
t36 = t105 * pkin(5) - t107 * qJ(6) + t58;
t381 = t302 * t36;
t380 = t302 * t58;
t218 = pkin(2) * t238 + pkin(3);
t146 = pkin(2) * t305 + t232 * t218;
t141 = pkin(9) + t146;
t379 = t141 * t68;
t211 = pkin(9) + t400;
t378 = t211 * t68;
t377 = t233 * t51;
t62 = t233 * t68;
t128 = t323 * t169;
t129 = t323 * t170;
t94 = -t128 * t368 - t232 * t129;
t376 = t233 * t94;
t375 = t237 * t52;
t63 = t237 * t68;
t374 = t237 * t94;
t373 = -t105 * t329 - t233 * t52;
t372 = t410 + t371;
t69 = t103 * t232 + t306;
t369 = -t69 + t410;
t365 = t105 * t266;
t364 = t105 * t302;
t363 = t105 * t233;
t362 = t107 * t105;
t361 = t107 * t266;
t360 = t107 * t237;
t359 = t406 * t266;
t356 = t126 * t237;
t355 = t144 * t237;
t354 = t154 * t152;
t352 = t210 * t240;
t224 = sin(t231);
t351 = t224 * t236;
t350 = t224 * t240;
t228 = -qJ(4) - t404;
t349 = t228 * t240;
t345 = t233 * t236;
t343 = t236 * t237;
t342 = t237 * t240;
t341 = t240 * t233;
t33 = -t233 * t59 + t237 * t73;
t340 = qJD(6) - t33;
t338 = t414 * t237;
t214 = pkin(3) * t225;
t336 = t214 + t226;
t229 = t235 ^ 2;
t335 = -t239 ^ 2 + t229;
t331 = qJD(5) * t211;
t222 = t235 * t383;
t35 = t36 * t330;
t55 = t58 * t330;
t318 = t210 * t392 - t233 * t414;
t316 = -t5 - t394;
t315 = qJD(2) * t404;
t314 = t368 * pkin(3);
t312 = -t13 - t394;
t308 = pkin(3) * t129 + t222;
t176 = t235 * t315;
t178 = t239 * t315;
t263 = -t238 * t176 - t234 * t178 - t188 * t332 - t189 * t333;
t76 = -qJ(4) * t129 - qJD(4) * t169 + t263;
t254 = -qJD(3) * t337 + t176 * t234 - t238 * t178;
t77 = qJ(4) * t128 - qJD(4) * t170 + t254;
t40 = t232 * t76 - t368 * t77;
t32 = t233 * t79 - t237 * t80 - t399;
t303 = t233 * t144 - t32;
t88 = -t368 * t117 + t118 * t232;
t296 = t238 * t323;
t136 = t210 * t345 + t342;
t138 = t210 * t341 - t343;
t291 = -g(1) * t136 + g(2) * t138;
t137 = t210 * t343 - t341;
t139 = t210 * t342 + t345;
t290 = g(1) * t137 - g(2) * t139;
t285 = pkin(5) * t233 - qJ(6) * t237;
t284 = t266 * t61 + t302 * t60;
t283 = -t379 - t380;
t282 = -t378 - t380;
t25 = -pkin(5) * t406 + t340;
t281 = -t233 * t26 + t237 * t25;
t280 = t233 * t25 + t237 * t26;
t145 = -pkin(2) * t348 + t218 * t368;
t278 = t210 * t286 + t214 - t409;
t276 = t25 * t266 + t338 + t35;
t275 = -t266 * t33 + t338 + t55;
t274 = -t406 * t419 + t63;
t273 = t406 * t413 + t62;
t272 = pkin(4) + t286;
t140 = -pkin(4) - t145;
t271 = -0.2e1 * pkin(1) * t327 - pkin(7) * qJDD(2);
t270 = -t126 * t330 + t374;
t41 = t232 * t77 + t368 * t76;
t93 = -t128 * t232 + t129 * t368;
t46 = pkin(4) * t93 - pkin(9) * t94 + t308;
t268 = t233 * t46 + t237 * t41 + t87 * t329 - t330 * t89;
t267 = -t141 * t330 + t355;
t261 = t13 * t233 + t266 * t34 + t58 * t329 + t318;
t260 = g(1) * t350 + g(2) * t351 - t187 * t154 + t127 - t393;
t259 = -t26 * t266 + t36 * t357 - t318 - t391;
t258 = g(1) * t138 + g(2) * t136 + t209 * t392 - t304;
t241 = qJD(2) ^ 2;
t257 = -pkin(7) * t241 + 0.2e1 * t366 + t408;
t242 = qJD(1) ^ 2;
t256 = pkin(1) * t242 - pkin(7) * qJDD(1) + t289;
t253 = qJD(5) * t281 + t411;
t252 = -t377 - t375 + (t360 + t363) * qJD(5);
t251 = t107 * t36 + qJDD(6) - t258;
t250 = -g(1) * t139 - g(2) * t137 - t237 * t395 + t269;
t249 = g(3) * t224 - t234 * t130 - t187 * t152 + t225 * t289 - t238 * t297 + t150;
t247 = t289 * t272 * t209;
t246 = -t210 * t289 + t25 * t413 - t26 * t419 - t395 + t411;
t244 = pkin(3) * t248 + qJDD(1) * t292 + t322;
t212 = -t314 - pkin(4);
t183 = pkin(9) * t352;
t181 = t236 * t210 * pkin(9);
t179 = -pkin(2) * t235 - pkin(3) * t224;
t174 = pkin(1) + t336;
t162 = -t314 - t272;
t161 = t240 * t174;
t149 = t207 - t418;
t132 = t140 - t286;
t115 = t266 * qJ(6);
t111 = -t152 ^ 2 + t154 ^ 2;
t92 = -t154 * t323 + t245;
t91 = t152 * t323 + t255;
t81 = pkin(5) * t107 + qJ(6) * t105;
t54 = t126 * t285 + t88;
t39 = -pkin(5) * t125 + t233 * t89 - t237 * t87;
t38 = qJ(6) * t125 + t386;
t31 = t115 + t387;
t30 = t233 * t70 - t237 * t83 - t399;
t29 = t115 + t388;
t28 = t105 * t406 - t51;
t19 = t237 * t299 - t377;
t18 = t273 - t361;
t17 = t274 + t365;
t9 = t285 * t94 + (qJD(5) * t286 - qJD(6) * t237) * t126 + t40;
t8 = (-t51 + t364) * t237 - t107 * t298 + t373;
t7 = -pkin(5) * t93 + qJD(5) * t386 + t233 * t41 - t237 * t46;
t6 = qJ(6) * t93 + qJD(6) * t125 + t268;
t1 = [qJDD(1), t408, t289, qJDD(1) * t229 + 0.2e1 * t235 * t310, 0.2e1 * t235 * t324 - 0.2e1 * t327 * t335, qJDD(2) * t235 + t239 * t241, qJDD(2) * t239 - t235 * t241, 0, t235 * t271 + t239 * t257, -t235 * t257 + t239 * t271, t154 * t128 + t170 * t255, t128 * t152 - t293 * t169 + t154 * t129 - t170 * t262 + (-t170 * t235 * t296 + (t169 * t235 - t170 * t239) * t234 * t323) * qJD(1), -t128 * t323 + t170 * t321, -t129 * t323 - t169 * t321, 0, t152 * t222 + t225 * t408 + t254 * t323 + t300 * t321 + (t149 - t418) * t169 - 0.2e1 * t187 * t129, -t389 * t293 + t149 * t170 + t187 * t128 - t263 * t323 - t337 * t321 + (t234 * t295 * t389 - t154 * t383) * t235 - t408 * t224, -t16 * t125 - t15 * t126 + t243 * t88 + t266 * t40 + t302 * t41 - t60 * t94 - t61 * t93 + t89 * t71 - t289, t16 * t89 + t61 * t41 - t15 * t88 - t60 * t40 + t244 * t292 + t131 * t308 - g(1) * (-t174 * t236 - t349) - g(2) * (-t228 * t236 + t161) t107 * t270 - t356 * t51 (-t105 * t237 - t107 * t233) * t94 + (t377 - t375 + (-t360 + t363) * qJD(5)) * t126, t107 * t93 - t125 * t51 + t270 * t406 + t356 * t68, -t126 * t62 - t105 * t93 - t125 * t52 + (-t126 * t329 - t376) * t406, t125 * t68 + t406 * t93, -t304 * t125 + t33 * t93 + t40 * t105 + t88 * t52 + ((-qJD(5) * t89 + t46) * t406 + t87 * t68 + t58 * qJD(5) * t126) * t237 + ((-qJD(5) * t87 - t41) * t406 - t89 * t68 + t13 * t126 + t58 * t94) * t233 + t290, -t268 * t406 - t386 * t68 - t269 * t125 - t34 * t93 + t40 * t107 - t88 * t51 + t58 * t374 + (t13 * t237 - t55) * t126 + t291, t36 * t376 + t105 * t9 - t406 * t7 - t125 * t4 - t25 * t93 - t39 * t68 + t52 * t54 + (t329 * t36 + t391) * t126 + t290, -t105 * t6 + t107 * t7 - t38 * t52 - t39 * t51 + t281 * t94 + t408 * t209 + (-qJD(5) * t280 - t2 * t233 + t237 * t4) * t126, -t36 * t374 - t107 * t9 + t406 * t6 + t125 * t2 + t26 * t93 + t38 * t68 + t51 * t54 + (-t237 * t5 + t35) * t126 - t291, t2 * t38 + t26 * t6 + t5 * t54 + t36 * t9 + t4 * t39 + t25 * t7 - g(1) * (-pkin(5) * t137 - qJ(6) * t136 - t349) - g(2) * (pkin(4) * t352 + pkin(5) * t139 + pkin(9) * t353 + qJ(6) * t138 + t161) + (-g(1) * (-t174 + t409) + g(2) * t228) * t236; 0, 0, 0, -t235 * t242 * t239, t335 * t242, t325, t324, qJDD(2), -g(3) * t239 + t235 * t256, g(3) * t235 + t239 * t256, -t354, t111, t91, t92, t321, -t313 + t159 * t323 + (-t175 * t323 - t297) * t234 + (-t152 * t334 + t238 * t321 - t323 * t333) * pkin(2) + t260, t339 * t323 + (-qJD(3) * t296 + t154 * t334 - t234 * t321) * pkin(2) + t249, -t145 * t243 + t146 * t71 + t266 * t371 + t302 * t370 + t284, t16 * t146 + t15 * t145 - t131 * (t220 - t402) - g(3) * t336 + t370 * t61 - t371 * t60 - t289 * t179, t19, t8, t18, t17, -t359, t140 * t52 + t312 * t237 + t283 * t233 + t371 * t105 + ((-qJD(5) * t141 - t80) * t237 - t370 * t233) * t406 + t275, -t140 * t51 + t283 * t237 + t371 * t107 + (-t267 + t387) * t406 + t261, t132 * t52 + t316 * t237 + (-t379 - t381) * t233 + t372 * t105 + (-t141 * t329 - t303) * t406 + t276, t303 * t107 + (t31 - t355) * t105 + t252 * t141 + t246, t132 * t51 + (-qJD(5) * t36 + t379) * t237 - t372 * t107 + (t267 - t31) * t406 + t259, t5 * t132 - t26 * t31 - t25 * t32 - g(1) * (t179 * t240 + t183) - g(2) * (t179 * t236 + t181) - g(3) * (t226 + t278) + t372 * t36 + t280 * t144 + t247 + t253 * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t354, t111, t91, t92, t321, -qJD(2) * t277 - t234 * t133 + t260, t301 * t323 + t249, -t243 * t314 - t266 * t69 - t302 * t70 + t400 * t71 + t284, t60 * t69 - t61 * t70 + (t131 * t154 + t15 * t368 + t16 * t232 + t224 * t289 - t393) * pkin(3), t19, t8, t18, t17, -t359, -t69 * t105 + t212 * t52 + (t406 * t70 + t282) * t233 + ((-t83 - t331) * t406 + t312) * t237 + t275, -t69 * t107 - t212 * t51 + t282 * t237 + (t211 * t330 + t388) * t406 + t261, t406 * t30 + t162 * t52 + (-t378 - t381) * t233 + t369 * t105 + (-t331 * t406 + t316) * t237 + t276, t105 * t29 - t107 * t30 + t211 * t252 + t246, t211 * t63 - t406 * t29 + t162 * t51 - t369 * t107 + (-t211 * t298 - t237 * t36) * qJD(5) + t259, t5 * t162 - t26 * t29 - t25 * t30 - g(1) * (-pkin(3) * t350 + t183) - g(2) * (-pkin(3) * t351 + t181) - g(3) * t278 + t369 * t36 + t253 * t211 + t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t266 ^ 2 - t302 ^ 2, t266 * t60 - t302 * t61 + t244 - t408, 0, 0, 0, 0, 0, t274 - t365, -t237 * t417 - t361 - t62, -t298 * t406 - t365 + t63 (t51 + t364) * t237 + t233 * t299 + t373, t273 + t361, -t266 * t36 + (-t4 + t416) * t237 + (t25 * t406 + t2) * t233 - t408; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t362, -t105 ^ 2 + t405, t28, -t52 + t299, t68, -t107 * t58 + t258 + t382, t105 * t58 + t33 * t406 - t250, -t105 * t81 - t251 + t382 + 0.2e1 * t403, pkin(5) * t51 - qJ(6) * t52 + (t26 - t34) * t107 + (t25 - t340) * t105, 0.2e1 * t384 - t105 * t36 + t107 * t81 + (0.2e1 * qJD(6) - t33) * t406 + t250, t2 * qJ(6) - t4 * pkin(5) - t36 * t81 - t25 * t34 - g(1) * (-pkin(5) * t138 + qJ(6) * t139) - g(2) * (-pkin(5) * t136 + qJ(6) * t137) + t340 * t26 + t285 * t395; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68 + t362, t28, -t405 - t417, t251 - t403 - t416;];
tau_reg  = t1;

% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRPR8
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 09:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRPR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:10:47
% EndTime: 2019-05-05 09:11:12
% DurationCPUTime: 11.26s
% Computational Cost: add. (29515->492), mult. (61778->710), div. (0->0), fcn. (49350->14), ass. (0->331)
t269 = cos(qJ(4));
t260 = sin(pkin(7));
t270 = cos(qJ(3));
t344 = qJD(2) * t270;
t332 = t260 * t344;
t252 = -qJD(4) + t332;
t250 = t252 ^ 2;
t265 = sin(qJ(4));
t262 = cos(pkin(7));
t345 = qJD(2) * t262;
t326 = qJD(3) + t345;
t266 = sin(qJ(3));
t346 = qJD(2) * t260;
t333 = t266 * t346;
t231 = t265 * t326 + t269 * t333;
t384 = t231 ^ 2;
t397 = -t384 - t250;
t341 = qJDD(2) * t260;
t239 = -qJD(3) * t333 + t270 * t341;
t234 = -qJDD(4) + t239;
t229 = t265 * t333 - t269 * t326;
t358 = t231 * t229;
t289 = t234 - t358;
t402 = t289 * t265;
t130 = t269 * t397 + t402;
t261 = sin(pkin(6));
t263 = cos(pkin(6));
t267 = sin(qJ(2));
t271 = cos(qJ(2));
t401 = t289 * t269;
t133 = t265 * t397 - t401;
t340 = t266 * qJDD(2);
t238 = (qJD(3) * t344 + t340) * t260;
t324 = qJDD(2) * t262 + qJDD(3);
t286 = t269 * t238 + t265 * t324;
t342 = qJD(4) - t252;
t154 = t229 * t342 - t286;
t302 = t133 * t266 - t154 * t270;
t76 = t260 * t130 + t262 * t302;
t89 = t133 * t270 + t154 * t266;
t446 = (t267 * t89 + t271 * t76) * t261 + t263 * (-t262 * t130 + t260 * t302);
t445 = pkin(2) * t76;
t442 = pkin(9) * t89;
t385 = t229 ^ 2;
t185 = -t250 - t385;
t288 = t234 + t358;
t403 = t269 * t288;
t123 = t185 * t265 - t403;
t439 = pkin(3) * t123;
t438 = pkin(10) * t123;
t404 = t265 * t288;
t125 = t185 * t269 + t404;
t437 = pkin(10) * t125;
t436 = t125 * t266;
t435 = t125 * t270;
t434 = t260 * t123;
t433 = t262 * t123;
t203 = t385 - t250;
t136 = t203 * t265 - t401;
t432 = t262 * t136;
t204 = -t384 + t250;
t135 = t204 * t269 - t404;
t184 = -t229 * qJD(4) + t286;
t359 = t229 * t252;
t155 = t359 - t184;
t431 = t262 * t135 - t260 * (t266 * (t204 * t265 + t403) - t270 * t155);
t383 = -2 * qJD(5);
t430 = pkin(3) * t130;
t429 = pkin(10) * t130;
t428 = pkin(10) * t133;
t425 = t266 * (t203 * t269 + t402);
t164 = t385 + t384;
t423 = pkin(3) * t164;
t420 = t155 * t265;
t419 = t155 * t269;
t418 = t164 * t266;
t417 = t164 * t270;
t349 = -g(3) + qJDD(1);
t373 = sin(pkin(12));
t374 = cos(pkin(12));
t293 = g(1) * t373 - g(2) * t374;
t405 = t263 * t293;
t280 = t261 * t349 + t405;
t294 = -g(1) * t374 - g(2) * t373;
t197 = t267 * t280 + t271 * t294;
t272 = qJD(2) ^ 2;
t188 = -t272 * pkin(2) + pkin(9) * t341 + t197;
t381 = pkin(3) * t270;
t322 = -pkin(10) * t266 - t381;
t237 = t322 * t346;
t323 = t326 ^ 2;
t196 = -t267 * t294 + t271 * t280;
t378 = pkin(9) * t260;
t274 = qJDD(2) * pkin(2) + t272 * t378 + t196;
t281 = -t261 * t293 + t263 * t349;
t407 = t260 * t281 + t262 * t274;
t347 = t407 * t270;
t112 = t266 * (t237 * t346 + t188) - t324 * pkin(3) - t323 * pkin(10) - t347;
t327 = -t265 * t238 + t269 * t324;
t183 = qJD(4) * t231 - t327;
t398 = t359 + t184;
t278 = t183 * pkin(4) - qJ(5) * t398 + t112;
t415 = t231 * t383 + t278;
t264 = sin(qJ(6));
t268 = cos(qJ(6));
t199 = -t268 * t229 - t252 * t264;
t201 = t229 * t264 - t252 * t268;
t158 = t201 * t199;
t180 = qJDD(6) + t184;
t400 = -t158 + t180;
t412 = t264 * t400;
t411 = t265 * t398;
t410 = t268 * t400;
t409 = t269 * t398;
t396 = t384 - t385;
t408 = t270 * t396;
t257 = t260 ^ 2;
t298 = qJD(2) * t326;
t406 = t257 * (-t262 * t272 + t298);
t296 = t260 * t298;
t244 = t270 * t296;
t399 = t238 + t244;
t243 = t266 * t296;
t211 = t239 - t243;
t189 = pkin(4) * t229 - qJ(5) * t231;
t120 = t270 * t188 + t266 * t407;
t113 = -t323 * pkin(3) + t324 * pkin(10) + t237 * t332 + t120;
t159 = t260 * t274 - t262 * t281;
t115 = -pkin(3) * t211 - pkin(10) * t399 - t159;
t73 = t265 * t113 - t269 * t115;
t57 = t234 * pkin(4) - t250 * qJ(5) + t231 * t189 + qJDD(5) + t73;
t42 = -pkin(5) * t155 + t288 * pkin(11) + t57;
t202 = pkin(5) * t231 + pkin(11) * t252;
t329 = -pkin(4) * t252 + t383;
t44 = -t385 * pkin(5) + t183 * pkin(11) + (-t202 + t329) * t231 + t278;
t26 = t264 * t44 - t268 * t42;
t27 = t264 * t42 + t268 * t44;
t12 = -t268 * t26 + t264 * t27;
t382 = pkin(4) + pkin(11);
t74 = t269 * t113 + t265 * t115;
t287 = -t250 * pkin(4) - t229 * t189 + t74;
t357 = t234 * qJ(5);
t43 = -t357 - t183 * pkin(5) - t385 * pkin(11) + (t383 - t202) * t252 + t287;
t394 = qJ(5) * t43 - t12 * t382;
t198 = t199 ^ 2;
t386 = t201 ^ 2;
t129 = -t198 - t386;
t128 = -t199 * qJD(6) + t264 * t183 - t268 * t234;
t223 = qJD(6) + t231;
t362 = t199 * t223;
t104 = t128 + t362;
t328 = t268 * t183 + t264 * t234;
t291 = (-qJD(6) + t223) * t201 + t328;
t65 = -t104 * t268 + t264 * t291;
t393 = qJ(5) * t129 - t382 * t65 - t12;
t297 = t128 - t362;
t375 = t268 * t43;
t220 = t223 ^ 2;
t331 = -t220 - t386;
t118 = t158 + t180;
t353 = t264 * t118;
t84 = t268 * t331 - t353;
t392 = qJ(5) * t297 - t382 * t84 + t375;
t100 = (qJD(6) + t223) * t201 - t328;
t376 = t264 * t43;
t138 = -t220 - t198;
t82 = t264 * t138 + t410;
t391 = qJ(5) * t100 - t382 * t82 + t376;
t284 = 0.2e1 * qJD(5) * t252 - t287;
t390 = -pkin(4) * t397 - qJ(5) * (t289 + t234) - t284;
t356 = t252 * t265;
t319 = -t269 * t183 - t229 * t356;
t334 = t270 * t358;
t355 = t252 * t269;
t389 = t262 * t319 + (t266 * (t183 * t265 - t229 * t355) + t334) * t260;
t320 = t265 * t184 - t231 * t355;
t388 = t262 * t320 + (t266 * (t269 * t184 + t231 * t356) - t334) * t260;
t295 = (t229 * t265 + t231 * t269) * t252;
t387 = t262 * t295 + (t270 * t234 + t266 * (t229 * t269 - t231 * t265) * t252) * t260;
t380 = pkin(4) * t265;
t379 = pkin(4) * t269;
t372 = t112 * t265;
t371 = t112 * t269;
t370 = t118 * t268;
t214 = t231 * t252;
t147 = t183 - t214;
t369 = t147 * t269;
t361 = t223 * t264;
t360 = t223 * t268;
t354 = t257 * t272;
t251 = t266 * t270 * t354;
t236 = t251 + t324;
t352 = t266 * t236;
t235 = -t251 + t324;
t350 = t270 * t235;
t339 = t265 * t158;
t338 = t269 * t158;
t337 = t201 * t361;
t258 = t266 ^ 2;
t336 = t258 * t354;
t259 = t270 ^ 2;
t335 = t259 * t354;
t330 = qJ(5) * t265 + pkin(3);
t40 = t265 * t73 + t269 * t74;
t54 = -t284 - t357;
t321 = -pkin(4) * t57 + qJ(5) * t54;
t148 = (-qJD(4) - t252) * t231 + t327;
t317 = pkin(4) * t155 + qJ(5) * t148;
t11 = t12 * t265 + t269 * t43;
t13 = t264 * t26 + t268 * t27;
t316 = t11 * t266 - t13 * t270;
t39 = t265 * t74 - t269 * t73;
t32 = t265 * t57 + t269 * t54;
t64 = t231 * t329 + t278;
t315 = t266 * t32 - t270 * t64;
t48 = t129 * t269 + t265 * t65;
t68 = t264 * t104 + t268 * t291;
t314 = t266 * t48 - t270 * t68;
t56 = t100 * t269 + t265 * t82;
t83 = t138 * t268 - t412;
t313 = t266 * t56 - t270 * t83;
t59 = t265 * t84 + t269 * t297;
t85 = -t264 * t331 - t370;
t312 = t266 * t59 - t270 * t85;
t309 = -t112 * t270 + t266 * t40;
t110 = t148 * t269 - t420;
t308 = t110 * t266 + t417;
t146 = t183 + t214;
t111 = -t146 * t269 - t420;
t307 = t111 * t266 + t417;
t119 = t266 * t188 - t347;
t306 = -t119 * t270 + t120 * t266;
t79 = t119 * t266 + t120 * t270;
t149 = -t231 * t342 + t327;
t305 = t149 * t270 + t436;
t304 = t147 * t270 - t436;
t209 = -t244 + t238;
t210 = t239 + t243;
t301 = -t209 * t270 + t210 * t266;
t221 = -t336 - t323;
t300 = t221 * t270 - t235 * t266;
t240 = -t323 - t335;
t299 = t236 * t270 + t240 * t266;
t277 = pkin(4) * t288 - qJ(5) * t185 + t57;
t242 = (-t258 - t259) * t354;
t241 = (t258 - t259) * t354;
t208 = (t340 + (0.2e1 * qJD(3) + t345) * t344) * t260;
t194 = t240 * t270 - t352;
t187 = -t221 * t266 - t350;
t167 = t220 - t386;
t166 = t198 - t220;
t165 = t199 * t360;
t160 = t209 * t266 + t210 * t270;
t157 = -t198 + t386;
t156 = t260 * t211 + t262 * t299;
t144 = -t260 * t208 + t262 * t300;
t139 = -t260 * t242 + t262 * t301;
t127 = -qJD(6) * t201 + t328;
t122 = t165 - t337;
t121 = (t199 * t264 + t201 * t268) * t223;
t109 = -t146 * t265 + t419;
t108 = -t147 * t265 + t409;
t107 = t149 * t265 + t409;
t106 = t148 * t265 + t419;
t98 = -t128 * t264 - t201 * t360;
t97 = -t127 * t264 + t165;
t96 = -t128 * t268 + t337;
t95 = -t127 * t268 - t199 * t361;
t94 = t121 * t269 + t180 * t265;
t93 = t167 * t264 - t410;
t92 = -t166 * t264 - t370;
t91 = -t166 * t268 + t353;
t90 = -t167 * t268 - t412;
t87 = -t147 * t266 - t435;
t86 = -t149 * t266 + t435;
t81 = t111 * t270 - t418;
t80 = t110 * t270 - t418;
t78 = t269 * t98 + t339;
t77 = t269 * t95 - t339;
t71 = t262 * t304 + t434;
t70 = t262 * t305 - t434;
t69 = t260 * t159 + t262 * t306;
t67 = t100 * t268 + t264 * t297;
t66 = t100 * t264 - t268 * t297;
t63 = pkin(3) * t154 + t372 - t428;
t62 = pkin(3) * t149 - t371 + t437;
t61 = t104 * t265 + t269 * t90;
t60 = t265 * t291 + t269 * t92;
t58 = t265 * t297 - t269 * t84;
t55 = t100 * t265 - t269 * t82;
t53 = -t260 * t109 + t262 * t307;
t52 = -t260 * t106 + t262 * t308;
t51 = t157 * t265 + t269 * t66;
t50 = (t147 - t214) * pkin(4) + t415;
t49 = pkin(4) * t214 - qJ(5) * t154 - t415;
t47 = t129 * t265 - t269 * t65;
t46 = qJ(5) * t164 + t57;
t45 = pkin(4) * t164 + t54;
t41 = t266 * t85 + t270 * t59;
t38 = t266 * t83 + t270 * t56;
t37 = pkin(5) * t65 - qJ(5) * t68;
t36 = t428 + t265 * t49 - (pkin(3) + t379) * t154;
t35 = t147 * t330 + t269 * t50 - t437;
t34 = -pkin(3) * t112 + pkin(10) * t40;
t33 = t112 * t266 + t270 * t40;
t31 = t265 * t54 - t269 * t57;
t30 = t266 * t68 + t270 * t48;
t29 = pkin(10) * t111 + t40 + t423;
t28 = -t260 * t58 + t262 * t312;
t24 = -t260 * t55 + t262 * t313;
t23 = pkin(5) * t297 - t382 * t85 - t376;
t22 = pkin(5) * t100 - t382 * t83 + t375;
t21 = pkin(10) * t110 + t265 * t46 + t269 * t45 + t423;
t20 = t266 * t64 + t270 * t32;
t19 = -t260 * t47 + t262 * t314;
t18 = -t260 * t39 + t262 * t309;
t17 = pkin(5) * t84 - qJ(5) * t85 - t27;
t16 = pkin(5) * t82 - qJ(5) * t83 - t26;
t15 = pkin(10) * t32 + (-t330 - t379) * t64;
t14 = -t260 * t31 + t262 * t315;
t10 = -t12 * t269 + t265 * t43;
t9 = pkin(5) * t129 - t382 * t68 - t13;
t8 = -pkin(3) * t85 + pkin(10) * t59 + t17 * t265 + t23 * t269;
t7 = -pkin(3) * t83 + pkin(10) * t56 + t16 * t265 + t22 * t269;
t6 = pkin(5) * t12 - qJ(5) * t13;
t5 = -pkin(3) * t68 + pkin(10) * t48 + t265 * t37 + t269 * t9;
t4 = pkin(5) * t43 - t13 * t382;
t3 = t11 * t270 + t13 * t266;
t2 = -t260 * t10 + t262 * t316;
t1 = -pkin(3) * t13 + pkin(10) * t11 + t265 * t6 + t269 * t4;
t25 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t349, 0, 0, 0, 0, 0, 0, (qJDD(2) * t271 - t267 * t272) * t261, (-qJDD(2) * t267 - t271 * t272) * t261, 0, t263 ^ 2 * t349 + (t271 * t196 + t267 * t197 - t405) * t261, 0, 0, 0, 0, 0, 0, t263 * (-t262 * t211 + t260 * t299) + (t156 * t271 + t194 * t267) * t261, t263 * (t262 * t208 + t260 * t300) + (t144 * t271 + t187 * t267) * t261, t263 * (t262 * t242 + t260 * t301) + (t139 * t271 + t160 * t267) * t261, t263 * (-t262 * t159 + t260 * t306) + (t267 * t79 + t271 * t69) * t261, 0, 0, 0, 0, 0, 0, t263 * (t260 * t305 + t433) + (t267 * t86 + t271 * t70) * t261, -t446, t263 * (t262 * t109 + t260 * t307) + (t267 * t81 + t271 * t53) * t261, t263 * (t260 * t309 + t262 * t39) + (t18 * t271 + t267 * t33) * t261, 0, 0, 0, 0, 0, 0, t263 * (t262 * t106 + t260 * t308) + (t267 * t80 + t271 * t52) * t261, t263 * (t260 * t304 - t433) + (t267 * t87 + t271 * t71) * t261, t446, t263 * (t260 * t315 + t262 * t31) + (t14 * t271 + t20 * t267) * t261, 0, 0, 0, 0, 0, 0, t263 * (t260 * t313 + t262 * t55) + (t24 * t271 + t267 * t38) * t261, t263 * (t260 * t312 + t262 * t58) + (t267 * t41 + t271 * t28) * t261, t263 * (t260 * t314 + t262 * t47) + (t19 * t271 + t267 * t30) * t261, t263 * (t262 * t10 + t260 * t316) + (t2 * t271 + t267 * t3) * t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t196, -t197, 0, 0, (t260 * t238 + t270 * t406) * t266, t262 * t241 + (t266 * t211 + t270 * t399) * t260, t262 * t209 + (t352 + t270 * (t323 - t336)) * t260, (t260 * t239 - t266 * t406) * t270, t262 * t210 + (t266 * (-t323 + t335) + t350) * t260, t262 * t324, pkin(2) * t156 - t262 * t119 + (pkin(9) * t194 + t159 * t270) * t260, pkin(2) * t144 - t262 * t120 + (pkin(9) * t187 - t159 * t266) * t260, pkin(2) * t139 + (pkin(9) * t160 + t79) * t260, pkin(2) * t69 + t378 * t79, t388, t262 * t107 + (t266 * (t149 * t269 - t411) - t408) * t260, t431, t389, t432 + (-t270 * t148 + t425) * t260, t387, pkin(2) * t70 + t262 * t62 + (t266 * (t372 - t438) + t270 * (t73 - t439) + pkin(9) * t86) * t260, -t445 + t262 * t63 + (t266 * (t371 - t429) + t270 * (t74 - t430) - t442) * t260, pkin(2) * t53 + t262 * t29 + (t266 * (-pkin(10) * t109 - t39) - t109 * t381 + pkin(9) * t81) * t260, pkin(2) * t18 + t262 * t34 + (pkin(9) * t33 + t322 * t39) * t260, t387, -t431, -t432 + (-t270 * t146 - t425) * t260, t388, t262 * t108 + (t266 * (-t369 - t411) - t408) * t260, t389, pkin(2) * t52 + t262 * t21 + (t266 * (-pkin(10) * t106 - t265 * t45 + t269 * t46) + t270 * (-pkin(3) * t106 - t317) + pkin(9) * t80) * t260, pkin(2) * t71 + t262 * t35 + (t266 * (qJ(5) * t369 - t265 * t50 + t438) + t270 * (-t277 + t439) + pkin(9) * t87) * t260, t445 + t262 * t36 + (t266 * (t154 * t380 + t269 * t49 + t429) + t270 * (-t390 + t430) + t442) * t260, pkin(2) * t14 + t262 * t15 + (t266 * (-pkin(10) * t31 + (-qJ(5) * t269 + t380) * t64) + t270 * (-pkin(3) * t31 - t321) + pkin(9) * t20) * t260, t262 * t78 + (t266 * (-t265 * t98 + t338) + t270 * t96) * t260, t262 * t51 + (t266 * (t157 * t269 - t265 * t66) + t270 * t67) * t260, t262 * t61 + (t266 * (t104 * t269 - t265 * t90) + t270 * t93) * t260, t262 * t77 + (t266 * (-t265 * t95 - t338) - t270 * t97) * t260, t262 * t60 + (t266 * (-t265 * t92 + t269 * t291) + t270 * t91) * t260, t262 * t94 + (t266 * (-t121 * t265 + t180 * t269) + t270 * t122) * t260, pkin(2) * t24 + t262 * t7 + (t266 * (-pkin(10) * t55 + t16 * t269 - t22 * t265) + t270 * (-pkin(3) * t55 - t391) + pkin(9) * t38) * t260, pkin(2) * t28 + t262 * t8 + (t266 * (-pkin(10) * t58 + t17 * t269 - t23 * t265) + t270 * (-pkin(3) * t58 - t392) + pkin(9) * t41) * t260, pkin(2) * t19 + t262 * t5 + (t266 * (-pkin(10) * t47 - t265 * t9 + t269 * t37) + t270 * (-pkin(3) * t47 - t393) + pkin(9) * t30) * t260, pkin(2) * t2 + t262 * t1 + (t266 * (-pkin(10) * t10 - t265 * t4 + t269 * t6) + t270 * (-pkin(3) * t10 - t394) + pkin(9) * t3) * t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t251, t241, t209, t251, t210, t324, -t119, -t120, 0, 0, t320, t107, t135, t319, t136, t295, t62, t63, t29, t34, t295, -t135, -t136, t320, t108, t319, t21, t35, t36, t15, t78, t51, t61, t77, t60, t94, t7, t8, t5, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t358, t396, -t155, -t358, t148, -t234, -t73, -t74, 0, 0, -t234, t155, t146, t358, t396, -t358, t317, t277, t390, t321, -t96, -t67, -t93, t97, -t91, -t122, t391, t392, t393, t394; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -t288, t397, t57, 0, 0, 0, 0, 0, 0, t82, t84, t65, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t157, t104, -t158, t291, t180, -t26, -t27, 0, 0;];
tauJ_reg  = t25;

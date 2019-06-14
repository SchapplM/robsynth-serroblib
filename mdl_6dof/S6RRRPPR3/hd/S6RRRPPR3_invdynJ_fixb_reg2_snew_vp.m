% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRRPPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:32:44
% EndTime: 2019-05-07 04:33:03
% DurationCPUTime: 8.83s
% Computational Cost: add. (17423->446), mult. (36907->541), div. (0->0), fcn. (25022->8), ass. (0->262)
t266 = sin(qJ(3));
t270 = cos(qJ(3));
t271 = cos(qJ(2));
t267 = sin(qJ(2));
t333 = qJD(1) * t267;
t224 = -t270 * t271 * qJD(1) + t266 * t333;
t254 = t267 * qJDD(1);
t327 = qJD(1) * qJD(2);
t319 = t271 * t327;
t235 = t254 + t319;
t255 = t271 * qJDD(1);
t320 = t267 * t327;
t236 = t255 - t320;
t153 = -t224 * qJD(3) + t270 * t235 + t266 * t236;
t261 = qJD(2) + qJD(3);
t343 = t261 * t224;
t390 = t153 - t343;
t259 = t261 ^ 2;
t337 = t267 * t270;
t226 = (t266 * t271 + t337) * qJD(1);
t380 = t226 ^ 2;
t388 = -t259 - t380;
t260 = qJDD(2) + qJDD(3);
t344 = t226 * t224;
t175 = -t344 - t260;
t423 = t175 * t266;
t130 = t270 * t388 + t423;
t422 = t175 * t270;
t133 = -t266 * t388 + t422;
t433 = pkin(7) * (t130 * t267 - t133 * t271);
t440 = -pkin(1) * t390 - t433;
t439 = pkin(2) * t130;
t381 = t224 ^ 2;
t387 = -t259 - t381;
t408 = t260 - t344;
t417 = t408 * t266;
t406 = t387 * t270 - t417;
t416 = t408 * t270;
t409 = t387 * t266 + t416;
t415 = pkin(7) * (t267 * t409 - t271 * t406);
t438 = pkin(8) * t130;
t437 = pkin(8) * t133;
t384 = -t381 - t380;
t403 = pkin(1) * t384;
t389 = t153 + t343;
t110 = t270 * t389;
t312 = t266 * t235 - t270 * t236;
t298 = -qJD(3) * t226 - t312;
t342 = t261 * t226;
t120 = t342 + t298;
t73 = t120 * t266 - t110;
t361 = t389 * t266;
t76 = t120 * t270 + t361;
t436 = -pkin(7) * (t267 * t73 - t271 * t76) - t403;
t435 = -pkin(2) * t390 + t437;
t432 = pkin(2) * t73;
t431 = pkin(8) * t73;
t420 = pkin(2) * t409;
t419 = pkin(8) * t406;
t418 = pkin(8) * t409;
t402 = pkin(2) * t384;
t430 = pkin(8) * t76 - t402;
t383 = t381 - t259;
t425 = t267 * (t270 * t383 + t423) - t271 * (-t266 * t383 + t422);
t424 = qJ(4) * t175;
t386 = t259 - t380;
t421 = t267 * (-t386 * t266 + t416) + t271 * (t270 * t386 + t417);
t378 = pkin(3) + pkin(4);
t412 = qJ(4) * t384;
t411 = t266 * t390;
t410 = t270 * t390;
t306 = t390 * qJ(4);
t264 = t271 ^ 2;
t273 = qJD(1) ^ 2;
t268 = sin(qJ(1));
t376 = cos(qJ(1));
t316 = t268 * g(1) - t376 * g(2);
t300 = qJDD(1) * pkin(1) + t316;
t301 = qJD(2) * pkin(2) - pkin(8) * t333;
t157 = t236 * pkin(2) + (pkin(8) * t264 + pkin(7)) * t273 - t301 * t333 + t300;
t281 = -pkin(3) * t342 + t157;
t276 = t306 + t281;
t382 = -pkin(4) * t298 + t381 * qJ(5) - qJDD(5);
t407 = t276 - t382;
t117 = t342 - t298;
t405 = 0.2e1 * t306;
t341 = t261 * t266;
t322 = t226 * t341;
t340 = t261 * t270;
t66 = t267 * (t270 * t153 - t322) + t271 * (t266 * t153 + t226 * t340);
t292 = (qJD(3) + t261) * t226 + t312;
t404 = t267 * (t292 * t270 + t411) - t271 * (-t292 * t266 + t410);
t328 = pkin(9) + t378;
t372 = pkin(5) + qJ(4);
t401 = qJ(5) * t389;
t265 = sin(qJ(6));
t269 = cos(qJ(6));
t191 = t224 * t265 + t269 * t261;
t193 = t224 * t269 - t261 * t265;
t143 = t193 * t191;
t150 = qJDD(6) + t153;
t391 = -t143 + t150;
t400 = t265 * t391;
t397 = t269 * t391;
t336 = t267 * t273;
t302 = g(1) * t376 + t268 * g(2);
t364 = qJDD(1) * pkin(7);
t229 = -t273 * pkin(1) - t302 + t364;
t339 = t267 * t229;
t140 = qJDD(2) * pkin(2) - t235 * pkin(8) - t339 + (pkin(2) * t336 + pkin(8) * t327 - g(3)) * t271;
t198 = -t267 * g(3) + t271 * t229;
t257 = t264 * t273;
t142 = -pkin(2) * t257 + t236 * pkin(8) - qJD(2) * t301 + t198;
t91 = -t270 * t140 + t266 * t142;
t296 = -t260 * pkin(3) - t259 * qJ(4) + qJDD(4) + t91;
t286 = -t260 * pkin(4) + t296 - t401;
t178 = pkin(3) * t224 - qJ(4) * t226;
t309 = pkin(4) * t224 - (2 * qJD(5)) + t178;
t40 = t226 * t309 + t286;
t385 = t380 - t381;
t186 = t191 ^ 2;
t187 = t193 ^ 2;
t219 = qJD(6) + t226;
t216 = t219 ^ 2;
t375 = pkin(3) * t266;
t374 = pkin(3) * t270;
t373 = t298 * pkin(3);
t331 = qJD(4) * t261;
t245 = 0.2e1 * t331;
t92 = t266 * t140 + t270 * t142;
t307 = t259 * pkin(3) - t260 * qJ(4) + t224 * t178 - t92;
t56 = t245 - t307;
t61 = t226 * t178 + t296;
t371 = -pkin(3) * t61 + qJ(4) * t56;
t179 = pkin(5) * t226 - pkin(9) * t224;
t199 = -pkin(4) * t261 - qJ(5) * t226;
t290 = pkin(4) * t381 - t261 * t199 + t307;
t284 = -qJ(5) * t298 + t245 - t290;
t33 = t260 * pkin(5) - t259 * pkin(9) + ((2 * qJD(5)) + t179) * t224 + t284;
t370 = t265 * t33;
t96 = t143 + t150;
t369 = t265 * t96;
t50 = t266 * t92 - t270 * t91;
t368 = t267 * t50;
t367 = t269 * t33;
t366 = t269 * t96;
t365 = qJ(4) * t387;
t360 = t157 * t266;
t359 = t157 * t270;
t347 = t219 * t191;
t346 = t219 * t265;
t345 = t219 * t269;
t241 = t271 * t336;
t338 = t267 * (qJDD(2) + t241);
t335 = t271 * (qJDD(2) - t241);
t119 = (-qJD(3) + t261) * t226 - t312;
t334 = -pkin(3) * t389 + qJ(4) * t119;
t332 = qJD(4) * t226;
t330 = qJD(5) * t224;
t329 = 0.2e1 * qJD(4) + t199;
t326 = -t187 - t216;
t325 = 0.2e1 * t332;
t324 = t266 * t143;
t323 = t270 * t143;
t100 = -t191 * qJD(6) - t265 * t260 - t269 * t298;
t318 = -qJ(4) * t266 - pkin(2);
t24 = -(-pkin(3) - pkin(9)) * t298 + t390 * pkin(5) + (-pkin(9) * t261 + t329) * t226 + t407;
t29 = -t259 * pkin(5) - t260 * pkin(9) + (-t179 + t309) * t226 + t286;
t15 = -t269 * t24 + t265 * t29;
t51 = t266 * t91 + t270 * t92;
t315 = t269 * t260 - t265 * t298;
t197 = t271 * g(3) + t339;
t313 = t267 * t197 + t271 * t198;
t42 = t284 + 0.2e1 * t330;
t311 = qJ(4) * t42 - t378 * t40;
t310 = -qJ(4) * t120 + t378 * t389;
t16 = t265 * t24 + t269 * t29;
t7 = -t269 * t15 + t265 * t16;
t8 = t265 * t15 + t269 * t16;
t304 = t100 - t347;
t297 = t271 * (-t224 * t266 - t226 * t270);
t294 = -t328 * t8 + t33 * t372;
t293 = (-qJD(6) + t219) * t193 - t315;
t291 = -pkin(3) * t388 - t424 + t56;
t109 = -t216 - t186;
t59 = t109 * t269 - t400;
t84 = (qJD(6) + t219) * t193 + t315;
t288 = -t328 * t59 + t372 * t84 + t367;
t65 = -t265 * t326 - t366;
t287 = t304 * t372 - t328 * t65 - t370;
t285 = pkin(3) * t408 + t365 - t61;
t102 = -t186 - t187;
t88 = -t100 - t347;
t49 = -t265 * t88 + t269 * t293;
t283 = t102 * t372 - t328 * t49 - t8;
t282 = t267 * (t224 * t340 - t266 * t298) + t271 * (t224 * t341 + t270 * t298);
t280 = t267 * t322 + (-t224 * t337 + t297) * t261;
t279 = t281 + t373;
t278 = -t378 * t388 + t42 - t424;
t277 = -t378 * t408 - t365 + t40;
t275 = t276 + t325;
t274 = t329 * t226 + t279 - t382;
t272 = qJD(2) ^ 2;
t263 = t267 ^ 2;
t256 = t263 * t273;
t237 = t255 - 0.2e1 * t320;
t234 = t254 + 0.2e1 * t319;
t228 = t273 * pkin(7) + t300;
t159 = -t187 + t216;
t158 = t186 - t216;
t141 = t187 - t186;
t99 = -qJD(6) * t193 - t315;
t98 = (t191 * t269 - t193 * t265) * t219;
t97 = (t191 * t265 + t193 * t269) * t219;
t93 = qJ(4) * t117 - qJ(5) * t408;
t80 = -t100 * t269 + t193 * t346;
t79 = -t100 * t265 - t193 * t345;
t78 = -t191 * t345 + t265 * t99;
t77 = -t191 * t346 - t269 * t99;
t75 = t119 * t270 + t361;
t72 = t119 * t266 - t110;
t70 = -t158 * t269 + t369;
t69 = t159 * t265 - t397;
t68 = -t158 * t265 - t366;
t67 = -t159 * t269 - t400;
t64 = t269 * t326 - t369;
t62 = qJ(5) * t175 + t378 * t390;
t58 = t265 * t109 + t397;
t53 = t61 - t412;
t52 = -pkin(3) * t384 + t56;
t48 = t265 * t304 + t269 * t84;
t47 = t265 * t293 + t269 * t88;
t46 = t265 * t84 - t269 * t304;
t44 = (-t292 + t298) * pkin(3) + t275;
t43 = t279 + t325 + t405;
t38 = t266 * t65 + t270 * t304;
t37 = t266 * t304 - t270 * t65;
t36 = t306 + t274;
t35 = t266 * t59 + t270 * t84;
t34 = t266 * t84 - t270 * t59;
t32 = t102 * t270 + t266 * t49;
t31 = t102 * t266 - t270 * t49;
t27 = t266 * t56 - t270 * t61;
t26 = -qJ(5) * t388 + t274 + t405;
t25 = -t40 + t401 + t412;
t22 = -0.2e1 * t330 - 0.2e1 * t331 + t378 * t384 + (t120 + t298) * qJ(5) + t290;
t21 = pkin(4) * t117 + qJ(5) * t387 - t226 * t199 - 0.2e1 * t332 + (t117 - t298) * pkin(3) - t407;
t20 = t266 * t40 + t270 * t42;
t19 = t266 * t42 - t270 * t40;
t18 = qJ(4) * t36 - qJ(5) * t40;
t17 = -qJ(5) * t49 + t372 * t47;
t13 = -qJ(5) * t42 + t36 * t378;
t12 = -qJ(5) * t304 + t328 * t64 - t367;
t11 = -qJ(5) * t84 + t328 * t58 - t370;
t10 = -qJ(5) * t65 + t372 * t64 - t16;
t9 = -qJ(5) * t59 + t372 * t58 - t15;
t5 = t266 * t8 + t270 * t33;
t4 = t266 * t33 - t270 * t8;
t3 = -qJ(5) * t102 + t328 * t47 + t7;
t2 = -qJ(5) * t8 + t372 * t7;
t1 = -qJ(5) * t33 + t328 * t7;
t6 = [0, 0, 0, 0, 0, qJDD(1), t316, t302, 0, 0, (t235 + t319) * t267, t234 * t271 + t237 * t267, t338 + t271 * (-t256 + t272), (t236 - t320) * t271, t267 * (t257 - t272) + t335, 0, t271 * t228 + pkin(1) * t237 + pkin(7) * (t271 * (-t257 - t272) - t338), -t267 * t228 - pkin(1) * t234 + pkin(7) * (-t335 - t267 * (-t256 - t272)), pkin(1) * (t256 + t257) + (t263 + t264) * t364 + t313, pkin(1) * t228 + pkin(7) * t313, t66, -t404, t421, t282, t425, t280, t267 * (-t360 - t418) + t271 * (-pkin(2) * t292 + t359 + t419) - pkin(1) * t292 - t415, t267 * (-t359 - t438) + t271 * (-t360 + t435) + t440, t267 * (-t50 - t431) + t271 * (t430 + t51) + t436, -pkin(8) * t368 + t271 * (pkin(2) * t157 + pkin(8) * t51) + pkin(1) * t157 + pkin(7) * (t271 * t51 - t368), t66, t421, t404, t280, -t425, t282, t267 * (-t266 * t44 - t418) + t271 * (t270 * t44 + t419) - t415 + (-qJ(4) * t337 + t271 * t318 - pkin(1)) * t292, t267 * (-pkin(8) * t72 - t266 * t52 + t270 * t53) + t271 * (pkin(8) * t75 + t266 * t53 + t270 * t52 - t402) - t403 + pkin(7) * (-t267 * t72 + t271 * t75), t267 * (t270 * t43 + t438) + t271 * (t266 * t43 - t437) + t433 + (-t267 * t375 + t271 * (pkin(2) + t374) + pkin(1)) * t390, (t267 * (qJ(4) * t270 - t375) + t271 * (-t318 + t374) + pkin(1)) * (t275 + t373) + (pkin(7) + pkin(8)) * (-t267 * t27 + t271 * (t266 * t61 + t270 * t56)), t282, t267 * (-t117 * t270 - t411) + t271 * (-t117 * t266 + t410), t425, t66, t421, (t267 * (-t224 * t270 + t226 * t266) + t297) * t261, t267 * (t26 * t270 - t266 * t62 + t438) + t271 * (t26 * t266 + t270 * t62 - t435) - t440, t267 * (-t21 * t266 + t270 * t93 + t418) + t271 * (pkin(2) * t117 + t21 * t270 + t266 * t93 - t419) + pkin(1) * t117 + t415, t267 * (-t22 * t266 + t25 * t270 + t431) + t271 * (t22 * t270 + t25 * t266 - t430) - t436, t267 * (-pkin(8) * t19 - t13 * t266 + t18 * t270) + t271 * (pkin(2) * t36 + pkin(8) * t20 + t13 * t270 + t18 * t266) + pkin(1) * t36 + pkin(7) * (-t19 * t267 + t20 * t271), t267 * (-t266 * t80 + t323) + t271 * (t270 * t80 + t324), t267 * (t141 * t270 - t266 * t48) + t271 * (t141 * t266 + t270 * t48), t267 * (-t266 * t69 - t270 * t88) + t271 * (-t266 * t88 + t270 * t69), t267 * (-t266 * t78 - t323) + t271 * (t270 * t78 - t324), t267 * (-t266 * t70 + t270 * t293) + t271 * (t266 * t293 + t270 * t70), t267 * (t150 * t270 - t266 * t98) + t271 * (t150 * t266 + t270 * t98), t267 * (-pkin(8) * t34 - t11 * t266 + t270 * t9) + t271 * (pkin(2) * t58 + pkin(8) * t35 + t11 * t270 + t266 * t9) + pkin(1) * t58 + pkin(7) * (-t267 * t34 + t271 * t35), t267 * (-pkin(8) * t37 + t10 * t270 - t12 * t266) + t271 * (pkin(2) * t64 + pkin(8) * t38 + t10 * t266 + t12 * t270) + pkin(1) * t64 + pkin(7) * (-t267 * t37 + t271 * t38), t267 * (-pkin(8) * t31 + t17 * t270 - t266 * t3) + t271 * (pkin(2) * t47 + pkin(8) * t32 + t17 * t266 + t270 * t3) + pkin(1) * t47 + pkin(7) * (-t267 * t31 + t271 * t32), t267 * (-pkin(8) * t4 - t1 * t266 + t2 * t270) + t271 * (pkin(2) * t7 + pkin(8) * t5 + t1 * t270 + t2 * t266) + pkin(1) * t7 + pkin(7) * (-t267 * t4 + t271 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, t256 - t257, t254, t241, t255, qJDD(2), -t197, -t198, 0, 0, t344, t385, t389, -t344, t119, t260, -t91 + t420, -t92 + t439, t432, pkin(2) * t50, t344, t389, -t385, t260, -t119, -t344, t285 + t420, pkin(2) * t72 + t334, t291 - t439, pkin(2) * t27 + t371, -t344, t385, t120, t344, t389, t260, t278 - t439, t277 - t420, t310 - t432, pkin(2) * t19 + t311, t79, t46, t67, t77, t68, t97, pkin(2) * t34 + t288, pkin(2) * t37 + t287, pkin(2) * t31 + t283, pkin(2) * t4 + t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t344, t385, t389, -t344, t119, t260, -t91, -t92, 0, 0, t344, t389, -t385, t260, -t119, -t344, t285, t334, t291, t371, -t344, t385, t120, t344, t389, t260, t278, t277, t310, t311, t79, t46, t67, t77, t68, t97, t288, t287, t283, t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t408, t389, t388, t61, 0, 0, 0, 0, 0, 0, t388, t408, -t389, t40, 0, 0, 0, 0, 0, 0, t59, t65, t49, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t390, t117, t384, t36, 0, 0, 0, 0, 0, 0, t58, t64, t47, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t141, -t88, -t143, t293, t150, -t15, -t16, 0, 0;];
tauJ_reg  = t6;

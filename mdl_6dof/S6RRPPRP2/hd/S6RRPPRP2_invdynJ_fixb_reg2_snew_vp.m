% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPPRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:11:44
% EndTime: 2019-05-06 09:12:03
% DurationCPUTime: 7.48s
% Computational Cost: add. (17150->432), mult. (40173->536), div. (0->0), fcn. (27491->8), ass. (0->259)
t258 = qJD(2) ^ 2;
t250 = sin(pkin(9));
t257 = cos(qJ(2));
t251 = cos(pkin(9));
t254 = sin(qJ(2));
t315 = t251 * t254;
t224 = (t250 * t257 + t315) * qJD(1);
t350 = t224 ^ 2;
t204 = t350 + t258;
t310 = qJD(1) * t254;
t222 = -t251 * t257 * qJD(1) + t250 * t310;
t317 = t224 * t222;
t374 = qJDD(2) + t317;
t387 = t374 * t250;
t135 = t204 * t251 + t387;
t386 = t374 * t251;
t137 = -t204 * t250 + t386;
t408 = pkin(7) * (t135 * t254 - t137 * t257);
t407 = pkin(2) * t135;
t406 = qJ(3) * t135;
t405 = qJ(3) * t137;
t238 = t254 * qJDD(1);
t302 = qJD(1) * qJD(2);
t294 = t257 * t302;
t231 = t238 + t294;
t239 = t257 * qJDD(1);
t295 = t254 * t302;
t232 = t239 - t295;
t184 = t231 * t250 - t251 * t232;
t309 = qJD(2) * t224;
t153 = t184 - t309;
t185 = t231 * t251 + t232 * t250;
t304 = t222 * qJD(2);
t357 = t185 + t304;
t393 = -t153 * t250 - t251 * t357;
t403 = pkin(2) * t393;
t402 = qJ(3) * t393;
t351 = t222 ^ 2;
t352 = -t351 - t350;
t392 = -t153 * t251 + t250 * t357;
t401 = -pkin(2) * t352 + qJ(3) * t392;
t205 = t350 - t258;
t168 = t317 - qJDD(2);
t396 = t168 * t251;
t397 = t168 * t250;
t400 = t254 * (-t205 * t250 + t396) + t257 * (t205 * t251 + t397);
t200 = t351 - t258;
t399 = t254 * (-t200 * t251 + t387) - t257 * (t200 * t250 + t386);
t398 = pkin(7) * (-t254 * t393 + t257 * t392) - pkin(1) * t352;
t253 = sin(qJ(5));
t256 = cos(qJ(5));
t194 = qJD(2) * t253 - t256 * t222;
t196 = qJD(2) * t256 + t222 * t253;
t151 = t196 * t194;
t179 = qJDD(5) + t185;
t365 = t151 - t179;
t389 = pkin(5) * t365;
t167 = -t258 - t351;
t122 = t167 * t250 - t396;
t125 = -t167 * t251 - t397;
t388 = pkin(7) * (t122 * t254 + t125 * t257);
t358 = t185 - t304;
t375 = t184 + t309;
t385 = t254 * (-t250 * t358 - t251 * t375) + t257 * (-t250 * t375 + t251 * t358);
t384 = pkin(2) * t122;
t383 = qJ(3) * t122;
t382 = qJ(3) * t125;
t361 = qJ(4) * t358;
t349 = pkin(3) + pkin(8);
t329 = t365 * t253;
t328 = t365 * t256;
t248 = t257 ^ 2;
t259 = qJD(1) ^ 2;
t273 = qJD(2) * pkin(2) - qJ(3) * t310;
t255 = sin(qJ(1));
t348 = cos(qJ(1));
t291 = g(1) * t255 - t348 * g(2);
t274 = qJDD(1) * pkin(1) + t291;
t148 = pkin(2) * t232 + (qJ(3) * t248 + pkin(7)) * t259 - t273 * t310 - qJDD(3) + t274;
t366 = pkin(3) * t309 - 0.2e1 * qJD(4) * t224 - t148;
t364 = 2 * qJD(3);
t192 = t194 ^ 2;
t217 = qJD(5) + t224;
t215 = t217 ^ 2;
t128 = -t215 - t192;
t85 = t128 * t253 - t328;
t86 = t128 * t256 + t329;
t363 = pkin(4) * t85 - qJ(4) * t86;
t193 = t196 ^ 2;
t139 = -t193 - t215;
t120 = t151 + t179;
t331 = t120 * t253;
t90 = t139 * t256 - t331;
t330 = t120 * t256;
t91 = -t139 * t253 - t330;
t362 = pkin(4) * t90 - qJ(4) * t91;
t132 = -qJD(5) * t194 + qJDD(2) * t256 + t184 * t253;
t164 = t217 * t194;
t360 = t132 - t164;
t308 = qJD(3) * t222;
t209 = -0.2e1 * t308;
t301 = qJD(4) * qJD(2);
t355 = t209 + 0.2e1 * t301;
t313 = t254 * t259;
t275 = g(1) * t348 + t255 * g(2);
t332 = qJDD(1) * pkin(7);
t228 = -t259 * pkin(1) - t275 + t332;
t316 = t228 * t254;
t143 = qJDD(2) * pkin(2) - qJ(3) * t231 - t316 + (pkin(2) * t313 + qJ(3) * t302 - g(3)) * t257;
t198 = -t254 * g(3) + t257 * t228;
t241 = t248 * t259;
t144 = -pkin(2) * t241 + t232 * qJ(3) - qJD(2) * t273 + t198;
t290 = -t251 * t143 + t144 * t250;
t270 = -qJDD(2) * pkin(3) - t258 * qJ(4) + qJDD(4) + t290;
t166 = pkin(3) * t222 - qJ(4) * t224;
t303 = t364 + t166;
t52 = -qJDD(2) * pkin(8) + t357 * pkin(4) + (pkin(8) * t222 + t303) * t224 + t270;
t199 = pkin(4) * t224 - qJD(2) * pkin(8);
t260 = -t361 + t366;
t66 = -pkin(4) * t351 + t184 * t349 - t199 * t224 + t260;
t34 = t253 * t66 - t256 * t52;
t354 = qJ(6) * t164 + 0.2e1 * qJD(6) * t196 + t34 + t389;
t287 = t253 * qJDD(2) - t256 * t184;
t131 = -t196 * qJD(5) - t287;
t160 = pkin(5) * t217 - qJ(6) * t196;
t35 = t253 * t52 + t256 * t66;
t26 = -t192 * pkin(5) + t131 * qJ(6) - 0.2e1 * qJD(6) * t194 - t217 * t160 + t35;
t353 = t350 - t351;
t347 = pkin(3) * t250;
t346 = pkin(3) * t251;
t127 = -t192 - t193;
t110 = t132 + t164;
t269 = (-qJD(5) + t217) * t196 - t287;
t73 = -t110 * t256 + t253 * t269;
t48 = t127 * t250 - t251 * t73;
t49 = t127 * t251 + t250 * t73;
t75 = t110 * t253 + t256 * t269;
t345 = pkin(7) * (-t254 * t48 + t257 * t49) - pkin(1) * t75;
t106 = (qJD(5) + t217) * t196 + t287;
t57 = t106 * t250 - t251 * t85;
t58 = t106 * t251 + t250 * t85;
t344 = pkin(7) * (-t254 * t57 + t257 * t58) - pkin(1) * t86;
t63 = t250 * t360 - t251 * t90;
t64 = t250 * t90 + t251 * t360;
t343 = pkin(7) * (-t254 * t63 + t257 * t64) - pkin(1) * t91;
t342 = qJ(3) * t48;
t341 = qJ(3) * t57;
t340 = qJ(3) * t63;
t25 = -qJ(6) * t132 - t354;
t337 = t25 * t253;
t336 = t25 * t256;
t311 = t250 * t143 + t251 * t144;
t280 = -t258 * pkin(3) + qJDD(2) * qJ(4) - t222 * t166 + t311;
t268 = -t184 * pkin(4) - pkin(8) * t351 + qJD(2) * t199 + t280;
t59 = t268 + t355;
t335 = t253 * t59;
t93 = t224 * t364 + t290;
t94 = t209 + t311;
t53 = t250 * t94 - t251 * t93;
t334 = t254 * t53;
t333 = t256 * t59;
t327 = t148 * t250;
t326 = t148 * t251;
t319 = t217 * t253;
t318 = t217 * t256;
t237 = t257 * t313;
t314 = t254 * (qJDD(2) + t237);
t312 = t257 * (qJDD(2) - t237);
t300 = t250 * t151;
t299 = t251 * t151;
t298 = -pkin(2) * t75 + qJ(3) * t49;
t297 = -pkin(2) * t86 + qJ(3) * t58;
t296 = -pkin(2) * t91 + qJ(3) * t64;
t293 = qJ(4) * t250 + pkin(2);
t40 = pkin(4) * t73 - qJ(4) * t75;
t54 = t250 * t93 + t251 * t94;
t197 = g(3) * t257 + t316;
t289 = t254 * t197 + t257 * t198;
t17 = t253 * t35 - t256 * t34;
t18 = t253 * t34 + t256 * t35;
t283 = pkin(2) * t57 + qJ(4) * t106 - t349 * t85;
t282 = pkin(2) * t63 + qJ(4) * t360 - t349 * t90;
t281 = pkin(2) * t48 + qJ(4) * t127 - t349 * t73;
t278 = pkin(4) * t106 - t349 * t86;
t277 = pkin(4) * t360 - t349 * t91;
t276 = pkin(4) * t127 - t349 * t75;
t271 = pkin(5) * t139 - t26;
t78 = t280 + t355;
t266 = t25 - t389;
t80 = t303 * t224 + t270;
t265 = t254 * (t251 * t185 - t250 * t309) + t257 * (t250 * t185 + t251 * t309);
t264 = t254 * (t184 * t250 + t251 * t304) + t257 * (-t251 * t184 + t250 * t304);
t263 = -t131 * pkin(5) - t192 * qJ(6) + t160 * t196 + qJDD(6) + t268;
t262 = (t254 * (-t222 * t251 + t224 * t250) + t257 * (-t222 * t250 - t224 * t251)) * qJD(2);
t41 = t263 + t355;
t261 = -pkin(3) * t184 - t366;
t247 = t254 ^ 2;
t240 = t247 * t259;
t233 = t239 - 0.2e1 * t295;
t230 = t238 + 0.2e1 * t294;
t227 = pkin(7) * t259 + t274;
t210 = 0.2e1 * t308;
t162 = -t193 + t215;
t161 = t192 - t215;
t146 = t193 - t192;
t117 = (-t194 * t256 + t196 * t253) * t217;
t116 = (t194 * t253 + t196 * t256) * t217;
t105 = pkin(5) * t110;
t102 = t132 * t256 - t196 * t319;
t101 = -t132 * t253 - t196 * t318;
t100 = -t131 * t253 + t194 * t318;
t99 = -t131 * t256 - t194 * t319;
t98 = t161 * t256 - t331;
t97 = -t162 * t253 - t328;
t96 = -t161 * t253 - t330;
t95 = -t162 * t256 + t329;
t81 = t261 + t361;
t79 = -pkin(5) * t360 - qJ(6) * t120;
t77 = (t375 + t184) * pkin(3) + t260;
t76 = t261 + 0.2e1 * t361;
t74 = -t106 * t256 - t253 * t360;
t72 = t106 * t253 - t256 * t360;
t68 = -qJ(4) * t352 + t80;
t67 = -pkin(3) * t352 + t78;
t60 = t254 * (-t116 * t250 + t179 * t251) + t257 * (t116 * t251 + t179 * t250);
t45 = t254 * (-t101 * t250 + t299) + t257 * (t101 * t251 + t300);
t44 = t254 * (-t250 * t99 - t299) + t257 * (t251 * t99 - t300);
t42 = t250 * t78 - t251 * t80;
t39 = -qJ(6) * t139 + t41;
t38 = t254 * (t110 * t251 - t250 * t95) + t257 * (t110 * t250 + t251 * t95);
t37 = t254 * (-t250 * t96 + t251 * t269) + t257 * (t250 * t269 + t251 * t96);
t32 = -t105 + t40;
t31 = t254 * (t146 * t251 - t250 * t72) + t257 * (t146 * t250 + t251 * t72);
t30 = -pkin(5) * t106 + qJ(6) * t128 + t210 - t263 - 0.2e1 * t301;
t28 = t277 - t335;
t27 = t278 + t333;
t24 = pkin(5) * t25;
t23 = -t35 + t362;
t22 = -t34 + t363;
t21 = (t110 + t132) * qJ(6) + t354;
t20 = -pkin(5) * t127 + qJ(6) * t269 + t26;
t19 = t271 + t362;
t16 = -t253 * t39 - t256 * t79 + t277;
t15 = -qJ(6) * t329 - t256 * t30 + t278;
t14 = -pkin(5) * t41 + qJ(6) * t26;
t13 = t266 + t363;
t12 = t17 * t250 + t251 * t59;
t11 = -t17 * t251 + t250 * t59;
t10 = t256 * t26 - t337;
t9 = t253 * t26 + t336;
t8 = -t18 + t276;
t7 = t250 * t9 + t251 * t41;
t6 = t250 * t41 - t251 * t9;
t5 = pkin(4) * t17 - qJ(4) * t18;
t4 = -t20 * t256 - t21 * t253 + t276;
t3 = pkin(4) * t59 - t18 * t349;
t2 = pkin(4) * t9 - qJ(4) * t10 + t24;
t1 = pkin(4) * t41 + qJ(6) * t337 - t10 * t349 - t14 * t256;
t29 = [0, 0, 0, 0, 0, qJDD(1), t291, t275, 0, 0, (t231 + t294) * t254, t230 * t257 + t233 * t254, t314 + t257 * (-t240 + t258), (t232 - t295) * t257, t254 * (t241 - t258) + t312, 0, t257 * t227 + pkin(1) * t233 + pkin(7) * (t257 * (-t241 - t258) - t314), -t254 * t227 - pkin(1) * t230 + pkin(7) * (-t312 - t254 * (-t240 - t258)), pkin(1) * (t240 + t241) + (t247 + t248) * t332 + t289, pkin(1) * t227 + pkin(7) * t289, t265, t385, -t400, t264, -t399, t262, t254 * (-t327 - t383) + t257 * (-pkin(2) * t375 + t326 - t382) - pkin(1) * t375 - t388, t254 * (-t326 + t406) + t257 * (-pkin(2) * t358 - t327 - t405) - pkin(1) * t358 + t408, t254 * (-t53 - t402) + t257 * (t401 + t54) + t398, -qJ(3) * t334 + t257 * (pkin(2) * t148 + qJ(3) * t54) + pkin(1) * t148 + pkin(7) * (t257 * t54 - t334), t262, t400, t399, t265, t385, t264, t254 * (-t250 * t67 + t251 * t68 - t402) + t257 * (t250 * t68 + t251 * t67 + t401) + t398, t254 * (-t250 * t77 + t383) + t257 * (t251 * t77 + t382) + t388 + (qJ(4) * t315 + t257 * t293 + pkin(1)) * t375, t254 * (t251 * t76 - t406) + t257 * (t250 * t76 + t405) - t408 + (-t254 * t347 + t257 * (pkin(2) + t346) + pkin(1)) * t358, (t254 * (qJ(4) * t251 - t347) + t257 * (t293 + t346) + pkin(1)) * t81 + (pkin(7) + qJ(3)) * (-t254 * t42 + t257 * (t250 * t80 + t251 * t78)), t45, t31, t38, t44, t37, t60, t254 * (t22 * t251 - t250 * t27 - t341) + t257 * (t22 * t250 + t251 * t27 + t297) + t344, t254 * (t23 * t251 - t250 * t28 - t340) + t257 * (t23 * t250 + t251 * t28 + t296) + t343, t254 * (-t250 * t8 + t251 * t40 - t342) + t257 * (t250 * t40 + t251 * t8 + t298) + t345, t254 * (-qJ(3) * t11 - t250 * t3 + t251 * t5) + t257 * (-pkin(2) * t18 + qJ(3) * t12 + t250 * t5 + t251 * t3) - pkin(1) * t18 + pkin(7) * (-t11 * t254 + t12 * t257), t45, t31, t38, t44, t37, t60, t254 * (t13 * t251 - t15 * t250 - t341) + t257 * (t13 * t250 + t15 * t251 + t297) + t344, t254 * (-t16 * t250 + t19 * t251 - t340) + t257 * (t16 * t251 + t19 * t250 + t296) + t343, t254 * (-t250 * t4 + t251 * t32 - t342) + t257 * (t250 * t32 + t251 * t4 + t298) + t345, t254 * (-qJ(3) * t6 - t1 * t250 + t2 * t251) + t257 * (-pkin(2) * t10 + qJ(3) * t7 + t1 * t251 + t2 * t250) - pkin(1) * t10 + pkin(7) * (-t254 * t6 + t257 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, t240 - t241, t238, t237, t239, qJDD(2), -t197, -t198, 0, 0, t317, t353, t357, -t317, -t153, qJDD(2), -t93 + t384, t210 - t311 - t407, t403, pkin(2) * t53, qJDD(2), -t357, t153, t317, t353, -t317, -pkin(3) * t357 - qJ(4) * t153 + t403, pkin(3) * t168 - qJ(4) * t167 - t384 + t80, pkin(3) * t204 + qJ(4) * t374 + t407 + t78, pkin(2) * t42 - pkin(3) * t80 + qJ(4) * t78, t102, t74, t97, t100, t98, t117, t283 + t335, t282 + t333, -t17 + t281, pkin(2) * t11 + qJ(4) * t59 - t17 * t349, t102, t74, t97, t100, t98, t117, qJ(6) * t328 - t253 * t30 + t283, -t253 * t79 + t256 * t39 + t282, -t20 * t253 + t21 * t256 + t281, pkin(2) * t6 + qJ(4) * t41 - qJ(6) * t336 - t14 * t253 - t349 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t375, t358, t352, -t148, 0, 0, 0, 0, 0, 0, t352, -t375, -t358, -t81, 0, 0, 0, 0, 0, 0, t86, t91, t75, t18, 0, 0, 0, 0, 0, 0, t86, t91, t75, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t357, -t168, -t204, t80, 0, 0, 0, 0, 0, 0, t85, t90, t73, t17, 0, 0, 0, 0, 0, 0, t85, t90, t73, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t146, t110, -t151, t269, t179, -t34, -t35, 0, 0, t151, t146, t110, -t151, t269, t179, t266, t271, -t105, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t360, t127, t41;];
tauJ_reg  = t29;

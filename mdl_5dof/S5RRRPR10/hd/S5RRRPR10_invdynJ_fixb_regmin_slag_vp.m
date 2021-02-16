% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR10
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
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:41:17
% EndTime: 2021-01-15 23:41:43
% DurationCPUTime: 8.12s
% Computational Cost: add. (7457->531), mult. (18918->740), div. (0->0), fcn. (15063->14), ass. (0->267)
t227 = cos(qJ(2));
t308 = qJD(1) * qJD(2);
t280 = t227 * t308;
t223 = sin(qJ(2));
t307 = qJDD(1) * t223;
t392 = t280 + t307;
t217 = sin(pkin(5));
t317 = qJD(1) * t227;
t283 = t217 * t317;
t391 = qJD(3) - t283;
t219 = cos(pkin(5));
t310 = t219 * qJD(1);
t199 = qJD(2) + t310;
t226 = cos(qJ(3));
t222 = sin(qJ(3));
t318 = qJD(1) * t217;
t288 = t223 * t318;
t267 = t222 * t288;
t128 = -t226 * t199 + t267;
t130 = t222 * t199 + t226 * t288;
t216 = sin(pkin(10));
t218 = cos(pkin(10));
t83 = t218 * t128 + t216 * t130;
t385 = qJD(5) + t83;
t221 = sin(qJ(5));
t225 = cos(qJ(5));
t256 = -t216 * t128 + t218 * t130;
t62 = t221 * t256 - t225 * t391;
t390 = t385 * t62;
t268 = t222 * t283;
t220 = qJ(4) + pkin(8);
t277 = qJD(3) * t220;
t301 = pkin(1) * t310;
t146 = -pkin(7) * t288 + t227 * t301;
t262 = pkin(2) * t223 - pkin(8) * t227;
t147 = t262 * t318;
t323 = t226 * t146 + t222 * t147;
t389 = qJ(4) * t268 + t226 * qJD(4) - t222 * t277 - t323;
t273 = -t222 * t146 + t226 * t147;
t327 = t226 * t227;
t388 = -t222 * qJD(4) - t226 * t277 - (pkin(3) * t223 - qJ(4) * t327) * t318 - t273;
t257 = -t221 * t391 - t225 * t256;
t387 = t257 * t385;
t314 = qJD(3) * t222;
t386 = -t268 + t314;
t339 = t217 * t227;
t369 = cos(qJ(1));
t290 = t369 * t223;
t224 = sin(qJ(1));
t330 = t224 * t227;
t241 = -t219 * t330 - t290;
t366 = g(1) * t241;
t252 = g(3) * t339 + t366;
t166 = t216 * t222 - t218 * t226;
t112 = t166 * t283;
t156 = t166 * qJD(3);
t384 = t112 - t156;
t167 = t216 * t226 + t218 * t222;
t321 = t391 * t167;
t383 = t392 * t217;
t275 = t225 * t385;
t305 = t219 * qJDD(1);
t198 = qJDD(2) + t305;
t313 = qJD(3) * t226;
t75 = -qJD(3) * t267 + t222 * t198 + t199 * t313 + t383 * t226;
t315 = qJD(2) * t227;
t285 = t222 * t315;
t76 = -t226 * t198 + t199 * t314 + (qJD(1) * (t223 * t313 + t285) + t222 * t307) * t217;
t48 = t216 * t75 + t218 * t76;
t47 = qJDD(5) + t48;
t349 = t221 * t47;
t382 = -t275 * t385 - t349;
t289 = t369 * t227;
t331 = t224 * t223;
t160 = t219 * t331 - t289;
t213 = qJ(3) + pkin(10);
t210 = sin(t213);
t211 = cos(t213);
t340 = t217 * t224;
t103 = t160 * t210 + t211 * t340;
t161 = -t219 * t290 - t330;
t291 = t217 * t369;
t239 = t161 * t210 - t211 * t291;
t341 = t217 * t223;
t381 = g(2) * t239 + g(3) * (-t210 * t341 + t219 * t211) + g(1) * t103;
t212 = t217 ^ 2;
t304 = 0.2e1 * t212;
t149 = pkin(7) * t283 + t223 * t301;
t117 = t199 * pkin(8) + t149;
t254 = -pkin(2) * t227 - pkin(8) * t223 - pkin(1);
t142 = t254 * t217;
t121 = qJD(1) * t142;
t74 = t226 * t117 + t222 * t121;
t59 = -t128 * qJ(4) + t74;
t350 = t216 * t59;
t73 = -t222 * t117 + t226 * t121;
t58 = -t130 * qJ(4) + t73;
t28 = t218 * t58 - t350;
t380 = t28 * t385;
t355 = -t216 * t389 + t388 * t218;
t353 = t216 * t388 + t389 * t218;
t312 = qJD(5) * t221;
t378 = -t225 * t47 + t312 * t385;
t264 = pkin(3) * t386 - t149;
t368 = pkin(1) * t223;
t320 = pkin(7) * t339 + t219 * t368;
t141 = t219 * pkin(8) + t320;
t324 = t226 * t141 + t222 * t142;
t306 = qJDD(1) * t227;
t197 = t217 * t306;
t281 = t223 * t308;
t265 = t217 * t281;
t144 = qJDD(3) - t197 + t265;
t377 = -t222 * t144 - t313 * t391;
t375 = t76 * pkin(3) + qJDD(4);
t162 = -t219 * t289 + t331;
t173 = t210 * t291;
t344 = t161 * t211;
t240 = t173 + t344;
t374 = t225 * t162 + t221 * t240;
t49 = -t216 * t76 + t218 * t75;
t21 = -qJD(5) * t257 - t225 * t144 + t221 * t49;
t209 = t226 * pkin(3) + pkin(2);
t101 = t166 * pkin(4) - t167 * pkin(9) - t209;
t185 = t220 * t226;
t282 = t222 * t220;
t115 = t218 * t185 - t216 * t282;
t364 = g(2) * t162;
t372 = t211 * (t364 - t366) - (-t321 * pkin(4) + t384 * pkin(9) + qJD(5) * t115 - t264) * t385 + t101 * t47;
t300 = pkin(1) * qJD(2) * t219;
t270 = qJD(1) * t300;
t298 = pkin(1) * t305;
t293 = -pkin(7) * t197 - t223 * t298 - t227 * t270;
t235 = -pkin(7) * t265 - t293;
t90 = t198 * pkin(8) + t235;
t246 = t262 * qJD(2);
t92 = (qJD(1) * t246 + qJDD(1) * t254) * t217;
t233 = -qJD(3) * t74 - t222 * t90 + t226 * t92;
t17 = t144 * pkin(3) - t75 * qJ(4) - t130 * qJD(4) + t233;
t245 = t117 * t314 - t121 * t313 - t222 * t92 - t226 * t90;
t19 = -t76 * qJ(4) - t128 * qJD(4) - t245;
t5 = t218 * t17 - t216 * t19;
t3 = -t144 * pkin(4) - t5;
t371 = (t130 * pkin(3) + pkin(4) * t256 + pkin(9) * t83) * t385 + t3 + t381;
t365 = g(2) * t160;
t363 = g(2) * t241;
t362 = g(2) * t224;
t361 = g(2) * t226;
t359 = g(3) * t217;
t358 = t198 * pkin(2);
t357 = t62 * t256;
t356 = t257 * t256;
t6 = t216 * t17 + t218 * t19;
t333 = t222 * t223;
t338 = t219 * t226;
t157 = t217 * t333 - t338;
t286 = t217 * t315;
t107 = -qJD(3) * t157 + t226 * t286;
t332 = t223 * t226;
t158 = t217 * t332 + t219 * t222;
t148 = t217 * t246;
t200 = pkin(7) * t341;
t337 = t219 * t227;
t150 = (pkin(1) * t337 - t200) * qJD(2);
t231 = -t324 * qJD(3) + t226 * t148 - t222 * t150;
t316 = qJD(2) * t223;
t287 = t217 * t316;
t33 = pkin(3) * t287 - t107 * qJ(4) - t158 * qJD(4) + t231;
t106 = qJD(3) * t158 + t217 * t285;
t244 = -t141 * t314 + t142 * t313 + t222 * t148 + t226 * t150;
t37 = -t106 * qJ(4) - t157 * qJD(4) + t244;
t14 = t216 * t33 + t218 * t37;
t54 = pkin(3) * t391 + t58;
t56 = t218 * t59;
t26 = t216 * t54 + t56;
t274 = -t222 * t141 + t226 * t142;
t61 = -pkin(3) * t339 - t158 * qJ(4) + t274;
t69 = -t157 * qJ(4) + t324;
t36 = t216 * t61 + t218 * t69;
t354 = pkin(4) * t288 - t355;
t311 = qJD(5) * t225;
t20 = t221 * t144 + t225 * t49 - t256 * t312 + t311 * t391;
t351 = t20 * t221;
t348 = t221 * t83;
t347 = t128 * t391;
t346 = t130 * t391;
t345 = t160 * t211;
t343 = t167 * t225;
t342 = t212 * qJD(1) ^ 2;
t336 = t221 * t227;
t334 = t222 * t217;
t328 = t225 * t227;
t322 = t221 * t162 - t225 * t173;
t151 = pkin(7) * t286 + t223 * t300;
t214 = t223 ^ 2;
t319 = -t227 ^ 2 + t214;
t309 = qJD(2) - t199;
t303 = g(3) * t341;
t296 = t227 * t342;
t295 = t217 * t336;
t294 = t217 * t328;
t269 = t383 * pkin(7) + t223 * t270 - t227 * t298;
t91 = t269 - t358;
t51 = t91 + t375;
t10 = t48 * pkin(4) - t49 * pkin(9) + t51;
t4 = t144 * pkin(9) + t6;
t292 = t225 * t10 - t221 * t4;
t272 = t199 + t310;
t271 = t198 + t305;
t87 = t106 * pkin(3) + t151;
t261 = g(1) * t160 + g(2) * t161;
t259 = -g(1) * t162 - t363;
t13 = -t216 * t37 + t218 * t33;
t25 = t218 * t54 - t350;
t35 = -t216 * t69 + t218 * t61;
t23 = pkin(9) * t391 + t26;
t116 = -t199 * pkin(2) - t146;
t81 = t128 * pkin(3) + qJD(4) + t116;
t29 = t83 * pkin(4) - pkin(9) * t256 + t81;
t9 = t221 * t29 + t225 * t23;
t8 = -t221 * t23 + t225 * t29;
t255 = t209 * t227 + t220 * t223;
t159 = -t209 * t223 + t220 * t227;
t253 = -t348 * t385 - t378;
t251 = g(1) * t369 + t362;
t250 = -g(1) * t224 + g(2) * t369;
t99 = -t216 * t157 + t218 * t158;
t78 = t221 * t99 + t294;
t140 = t200 + (-pkin(1) * t227 - pkin(2)) * t219;
t249 = t210 * t340 - t345;
t248 = -t160 * t226 + t224 * t334;
t247 = t217 * t226 + t219 * t333;
t96 = -t221 * t112 - t225 * t288;
t243 = -t156 * t221 + t167 * t311 - t96;
t97 = -t225 * t112 + t221 * t288;
t242 = -t156 * t225 - t167 * t312 - t97;
t236 = t252 - t364;
t100 = t157 * pkin(3) + t140;
t139 = t219 * t210 + t211 * t341;
t234 = g(1) * t249 + g(3) * t139;
t232 = -t221 * t249 - t225 * t241;
t230 = -t236 - t269;
t229 = -t115 * t47 + t3 * t167 + (pkin(9) * t288 - qJD(5) * t101 - t353) * t385 + t261;
t206 = -t218 * pkin(3) - pkin(4);
t205 = t216 * pkin(3) + pkin(9);
t154 = pkin(1) + t255;
t127 = t255 * t219;
t114 = t216 * t185 + t218 * t282;
t109 = -t161 * t222 + t226 * t291;
t105 = t217 * (t222 * pkin(3) + pkin(7)) + t159 * t219;
t98 = t218 * t157 + t216 * t158;
t79 = t225 * t99 - t295;
t68 = -t216 * t106 + t218 * t107;
t67 = t218 * t106 + t216 * t107;
t50 = t98 * pkin(4) - t99 * pkin(9) + t100;
t39 = -qJD(5) * t295 + t221 * t68 - t225 * t287 + t99 * t311;
t38 = -qJD(5) * t78 + t221 * t287 + t225 * t68;
t32 = -pkin(9) * t339 + t36;
t31 = pkin(4) * t339 - t35;
t27 = t216 * t58 + t56;
t24 = t67 * pkin(4) - t68 * pkin(9) + t87;
t22 = -pkin(4) * t391 - t25;
t12 = pkin(9) * t287 + t14;
t11 = -pkin(4) * t287 - t13;
t2 = -qJD(5) * t9 + t292;
t1 = qJD(5) * t8 + t221 * t10 + t225 * t4;
t7 = [qJDD(1), -t250, t251, (qJDD(1) * t214 + 0.2e1 * t223 * t280) * t212, (t223 * t306 - t319 * t308) * t304, (t223 * t271 + t272 * t315) * t217, (t227 * t271 - t272 * t316) * t217, t198 * t219, -t151 * t199 - t200 * t198 - t269 * t219 - g(1) * t161 + t365 + (t198 * t337 + (-t281 + t306) * t304) * pkin(1), -pkin(1) * t304 * t392 - t150 * t199 - t320 * t198 - t235 * t219 + t259, t130 * t107 + t75 * t158, -t130 * t106 - t107 * t128 - t75 * t157 - t158 * t76, t107 * t391 + t158 * t144 + (t130 * t316 - t227 * t75) * t217, -t106 * t391 - t157 * t144 + (-t128 * t316 + t227 * t76) * t217, (-t144 * t227 + t316 * t391) * t217, t231 * t391 + t274 * t144 - t233 * t339 + t73 * t287 + t151 * t128 + t140 * t76 + t91 * t157 + t116 * t106 - g(1) * ((-t219 * t332 + t334) * t369 - t224 * t327) - g(2) * t248, -t244 * t391 - t324 * t144 + t151 * t130 + t140 * t75 + t91 * t158 + t116 * t107 - g(1) * t109 - t222 * t365 + (-t224 * t361 - t227 * t245 - t74 * t316) * t217, t13 * t391 + t35 * t144 + t87 * t83 + t100 * t48 + t51 * t98 + t81 * t67 - g(1) * t240 + g(2) * t345 + (-t210 * t362 - t5 * t227 + t25 * t316) * t217, g(1) * t239 - g(2) * t103 + t100 * t49 - t14 * t391 - t36 * t144 + t51 * t99 + t81 * t68 + t87 * t256 + (t227 * t6 - t26 * t316) * t217, -t13 * t256 - t14 * t83 - t25 * t68 - t26 * t67 - t35 * t49 - t36 * t48 - t5 * t99 - t6 * t98 - t259, t6 * t36 + t26 * t14 + t5 * t35 + t25 * t13 + t51 * t100 + t81 * t87 - g(1) * (t105 * t369 - t154 * t224) - g(2) * (t105 * t224 + t154 * t369), t20 * t79 - t257 * t38, -t20 * t78 - t21 * t79 + t257 * t39 - t38 * t62, t20 * t98 - t257 * t67 + t38 * t385 + t47 * t79, -t21 * t98 - t385 * t39 - t47 * t78 - t62 * t67, t385 * t67 + t47 * t98, t2 * t98 + t8 * t67 + t11 * t62 + t31 * t21 + t3 * t78 + t22 * t39 + g(1) * t322 + ((-qJD(5) * t50 - t12) * t385 - t32 * t47 + t363) * t221 + ((-qJD(5) * t32 + t24) * t385 + t50 * t47 - g(1) * t344 - g(2) * t249) * t225, -(t225 * t12 + t221 * t24 + (-t221 * t32 + t225 * t50) * qJD(5)) * t385 - (t221 * t50 + t225 * t32) * t47 - t1 * t98 - t9 * t67 - t11 * t257 + t31 * t20 + t3 * t79 + t22 * t38 + g(1) * t374 - g(2) * t232; 0, 0, 0, -t223 * t296, t319 * t342, (t309 * t317 + t307) * t217, -t309 * t288 + t197, t198, t149 * t199 + t342 * t368 + t230, pkin(1) * t296 + t146 * t199 + (pkin(7) * t308 + g(3)) * t341 - t261 + t293, t75 * t222 + t226 * t346, (t75 - t347) * t226 + (-t346 - t76) * t222, (-t130 * t223 - t327 * t391) * t318 - t377, -t391 * t314 + t226 * t144 + (t222 * t227 * t391 + t128 * t223) * t318, -t391 * t288, -pkin(2) * t76 - t149 * t128 + t162 * t361 - t391 * t273 - t288 * t73 + t386 * t116 + t377 * pkin(8) + (-t91 - t252) * t226, -pkin(2) * t75 + t323 * t391 + t74 * t288 - t149 * t130 + (-pkin(8) * t144 + t116 * t391) * t226 + (pkin(8) * qJD(3) * t391 + t236 + t91) * t222, -t114 * t144 + t51 * t166 - t209 * t48 - t236 * t211 - t25 * t288 + t264 * t83 + t321 * t81 + t355 * t391, -t115 * t144 + t51 * t167 - t209 * t49 + t210 * t236 + t256 * t264 + t26 * t288 - t353 * t391 + t384 * t81, t114 * t49 - t115 * t48 - t6 * t166 - t5 * t167 - t25 * t384 - t256 * t355 - t26 * t321 - t353 * t83 + t261 - t303, t6 * t115 - t5 * t114 - t51 * t209 - g(1) * (-t127 * t224 + t159 * t369) - g(2) * (t127 * t369 + t159 * t224) + t264 * t81 + t353 * t26 + t355 * t25 - t255 * t359, t20 * t343 - t242 * t257, t97 * t62 - t257 * t96 - (t221 * t257 - t225 * t62) * t156 + (-t351 - t21 * t225 + (t221 * t62 + t225 * t257) * qJD(5)) * t167, t20 * t166 + t242 * t385 - t257 * t321 + t343 * t47, -t21 * t166 - t167 * t349 - t243 * t385 - t321 * t62, t47 * t166 + t321 * t385, t114 * t21 + t2 * t166 + t321 * t8 + t354 * t62 + t372 * t225 + t229 * t221 - (t211 * t328 + t221 * t223) * t359 + t243 * t22, -t1 * t166 + t114 * t20 - t321 * t9 - t354 * t257 - t372 * t221 + t229 * t225 - (-t211 * t336 + t223 * t225) * t359 + t242 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130 * t128, -t128 ^ 2 + t130 ^ 2, t75 + t347, -t76 + t346, t144, t74 * t391 - t116 * t130 - g(1) * (-t222 * t289 + t224 * t247) + g(2) * t109 + g(3) * t157 + t233, t73 * t391 + t116 * t128 + g(1) * t248 - g(2) * (t161 * t226 + t222 * t291) + g(3) * t158 + t245, t27 * t391 - t81 * t256 + (-t130 * t83 + t144 * t218) * pkin(3) + t5 - t381, t28 * t391 + t81 * t83 - g(2) * t240 + (-t130 * t256 - t144 * t216) * pkin(3) + t234 - t6, (-t216 * t48 - t218 * t49) * pkin(3) + (t26 - t27) * t256 + (-t25 + t28) * t83, t25 * t27 - t26 * t28 + (t250 * t247 - g(3) * t338 - t81 * t130 + t6 * t216 + t5 * t218 + (t227 * t251 + t303) * t222) * pkin(3), -t257 * t275 + t351, (t20 - t390) * t225 + (-t21 + t387) * t221, t356 - t382, t253 + t357, -t385 * t256, t221 * t380 + t206 * t21 - t27 * t62 - t8 * t256 + (t312 + t348) * t22 + (-t311 * t385 - t349) * t205 - t371 * t225, t225 * t380 + t206 * t20 + t27 * t257 + t9 * t256 + (t225 * t83 + t311) * t22 + t378 * t205 + t371 * t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256 * t391 + t48, -t391 * t83 + t49, -t256 ^ 2 - t83 ^ 2, t25 * t256 + t26 * t83 - t230 - t358 + t375, 0, 0, 0, 0, 0, t253 - t357, t356 + t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t257 * t62, t257 ^ 2 - t62 ^ 2, t20 + t390, -t21 - t387, t47, -t23 * t311 - t29 * t312 + t9 * t385 + t22 * t257 - g(1) * t232 - g(2) * t374 - g(3) * (-t139 * t221 - t294) + t292, t8 * t385 + t22 * t62 + g(2) * t322 + (qJD(5) * t23 - t10 - t252) * t221 + (-g(2) * t344 - qJD(5) * t29 + t234 - t4) * t225;];
tau_reg = t7;

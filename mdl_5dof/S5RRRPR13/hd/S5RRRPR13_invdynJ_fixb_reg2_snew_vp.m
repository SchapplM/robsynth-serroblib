% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPR13
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPR13_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:35
% EndTime: 2019-12-31 21:46:54
% DurationCPUTime: 7.24s
% Computational Cost: add. (19502->428), mult. (42616->589), div. (0->0), fcn. (32692->10), ass. (0->268)
t228 = cos(pkin(5));
t222 = qJD(1) * t228 + qJD(2);
t230 = sin(qJ(3));
t233 = cos(qJ(3));
t227 = sin(pkin(5));
t231 = sin(qJ(2));
t288 = qJD(1) * t231;
t277 = t227 * t288;
t194 = t222 * t230 + t233 * t277;
t234 = cos(qJ(2));
t289 = qJD(1) * t227;
t275 = qJD(2) * t289;
t285 = qJDD(1) * t227;
t202 = t231 * t285 + t234 * t275;
t270 = qJDD(1) * t228 + qJDD(2);
t150 = qJD(3) * t194 + t202 * t230 - t233 * t270;
t276 = t234 * t289;
t215 = -qJD(3) + t276;
t302 = t215 * t194;
t119 = t150 - t302;
t213 = t215 ^ 2;
t192 = -t233 * t222 + t230 * t277;
t332 = t192 ^ 2;
t152 = -t213 - t332;
t203 = -t231 * t275 + t234 * t285;
t197 = -qJDD(3) + t203;
t303 = t194 * t192;
t250 = t197 + t303;
t350 = t230 * t250;
t94 = t152 * t233 + t350;
t388 = pkin(1) * (-t119 * t234 + t231 * t94);
t349 = t233 * t250;
t92 = t152 * t230 - t349;
t387 = pkin(7) * (t119 * t231 + t234 * t94) - pkin(1) * t92;
t385 = pkin(2) * t92;
t384 = pkin(8) * t92;
t383 = pkin(8) * t94;
t331 = t194 ^ 2;
t343 = -t331 - t213;
t251 = t197 - t303;
t347 = t251 * t233;
t102 = t230 * t343 - t347;
t248 = t233 * t202 + t230 * t270;
t126 = (qJD(3) - t215) * t192 - t248;
t382 = pkin(1) * (t102 * t231 - t126 * t234);
t168 = -t331 + t213;
t104 = t168 * t233 - t350;
t381 = t228 * t104;
t348 = t251 * t230;
t99 = t233 * t343 + t348;
t380 = -pkin(7) * (t102 * t234 + t126 * t231) - pkin(1) * t99;
t167 = t332 - t213;
t105 = t167 * t230 - t347;
t345 = t150 + t302;
t377 = t228 * t105 + (t231 * (t167 * t233 + t348) + t234 * t345) * t227;
t329 = -2 * qJD(4);
t375 = pkin(2) * t99;
t374 = pkin(8) * t99;
t373 = pkin(8) * t102;
t372 = t231 * (t168 * t230 + t349);
t132 = t332 + t331;
t368 = pkin(2) * t132;
t367 = t119 * t230;
t313 = t119 * t233;
t362 = t132 * t231;
t361 = t132 * t234;
t151 = -t192 * qJD(3) + t248;
t175 = t192 * t215;
t124 = t175 + t151;
t269 = -t234 * pkin(2) - t231 * pkin(8);
t201 = t269 * t289;
t235 = qJD(1) ^ 2;
t326 = sin(qJ(1));
t327 = cos(qJ(1));
t255 = t327 * g(1) + t326 * g(2);
t198 = -t235 * pkin(1) + pkin(7) * t285 - t255;
t254 = t326 * g(1) - t327 * g(2);
t296 = t227 * t235;
t244 = qJDD(1) * pkin(1) + pkin(7) * t296 + t254;
t240 = t228 * t244;
t271 = t198 * t231 - t234 * t240;
t330 = t222 ^ 2;
t109 = (g(3) * t234 + t201 * t288) * t227 - t270 * pkin(2) - t330 * pkin(8) + t271;
t237 = t150 * pkin(3) - t124 * qJ(4) + t109;
t360 = t194 * t329 + t237;
t229 = sin(qJ(5));
t232 = cos(qJ(5));
t163 = -t232 * t192 - t215 * t229;
t165 = t192 * t229 - t215 * t232;
t128 = t165 * t163;
t147 = qJDD(5) + t151;
t346 = -t128 + t147;
t357 = t229 * t346;
t356 = t230 * t345;
t355 = t232 * t346;
t354 = t233 * t345;
t342 = t331 - t332;
t353 = t234 * t342;
t224 = t227 ^ 2;
t290 = qJD(1) * t222;
t351 = t224 * (-t228 * t235 + t290);
t125 = -t175 + t151;
t297 = t227 * t234;
t153 = g(3) * t297 + t271;
t298 = t227 * t231;
t236 = -g(3) * t298 + t231 * t240;
t154 = t234 * t198 + t236;
t344 = t231 * t153 + t234 * t154;
t155 = pkin(3) * t192 - qJ(4) * t194;
t110 = t270 * pkin(8) - t330 * pkin(2) + (t201 * t289 + t198) * t234 + t236;
t322 = t228 * g(3);
t111 = -t203 * pkin(2) - t202 * pkin(8) - t322 + ((pkin(2) * t231 - pkin(8) * t234) * t290 - t244) * t227;
t67 = t110 * t230 - t233 * t111;
t54 = t197 * pkin(3) - t213 * qJ(4) + t155 * t194 + qJDD(4) + t67;
t32 = t125 * pkin(4) + t250 * pkin(9) + t54;
t166 = pkin(4) * t194 + pkin(9) * t215;
t273 = -pkin(3) * t215 + t329;
t34 = -pkin(4) * t332 + pkin(9) * t150 + (-t166 + t273) * t194 + t237;
t19 = t229 * t34 - t232 * t32;
t20 = t229 * t32 + t232 * t34;
t10 = -t19 * t232 + t20 * t229;
t328 = pkin(3) + pkin(9);
t68 = t110 * t233 + t111 * t230;
t249 = -t213 * pkin(3) - t155 * t192 + t68;
t316 = qJ(4) * t197;
t33 = -t316 - t150 * pkin(4) - t332 * pkin(9) + (t329 - t166) * t215 + t249;
t340 = qJ(4) * t33 - t328 * t10;
t186 = qJD(5) + t194;
t272 = t232 * t150 + t197 * t229;
t252 = (-qJD(5) + t186) * t165 + t272;
t306 = t163 * t186;
t97 = -t163 * qJD(5) + t229 * t150 - t232 * t197;
t79 = t97 + t306;
t47 = t229 * t252 - t232 * t79;
t161 = t163 ^ 2;
t162 = t165 ^ 2;
t98 = -t161 - t162;
t339 = qJ(4) * t98 - t328 * t47 - t10;
t259 = t97 - t306;
t318 = t232 * t33;
t183 = t186 ^ 2;
t284 = -t162 - t183;
t89 = t128 + t147;
t319 = t229 * t89;
t60 = t232 * t284 - t319;
t338 = qJ(4) * t259 - t328 * t60 + t318;
t320 = t229 * t33;
t107 = -t183 - t161;
t58 = t229 * t107 + t355;
t75 = (qJD(5) + t186) * t165 - t272;
t337 = qJ(4) * t75 - t328 * t58 + t320;
t246 = 0.2e1 * qJD(4) * t215 - t249;
t336 = -pkin(3) * t343 - qJ(4) * (t251 + t197) - t246;
t301 = t215 * t230;
t266 = -t233 * t150 - t192 * t301;
t279 = t234 * t303;
t300 = t215 * t233;
t335 = t228 * t266 + (t231 * (t150 * t230 - t192 * t300) + t279) * t227;
t267 = t230 * t151 - t194 * t300;
t334 = t228 * t267 + (t231 * (t233 * t151 + t194 * t301) - t279) * t227;
t256 = (t192 * t230 + t194 * t233) * t215;
t333 = t228 * t256 + t197 * t297 + (t192 * t233 - t194 * t230) * t215 * t298;
t325 = pkin(3) * t230;
t324 = pkin(3) * t233;
t323 = t227 * pkin(7);
t317 = t232 * t89;
t315 = t109 * t230;
t314 = t109 * t233;
t305 = t186 * t229;
t304 = t186 * t232;
t299 = t224 * t235;
t214 = t234 * t231 * t299;
t200 = t214 + t270;
t294 = t231 * t200;
t199 = -t214 + t270;
t292 = t234 * t199;
t283 = t230 * t128;
t282 = t233 * t128;
t225 = t231 ^ 2;
t281 = t225 * t299;
t226 = t234 ^ 2;
t280 = t226 * t299;
t207 = t222 * t276;
t278 = t207 + t202;
t274 = qJ(4) * t230 + pkin(2);
t40 = t230 * t67 + t233 * t68;
t52 = -t246 - t316;
t268 = -pkin(3) * t54 + qJ(4) * t52;
t122 = (-qJD(3) - t215) * t192 + t248;
t265 = -pkin(3) * t122 - qJ(4) * t345;
t11 = t229 * t19 + t20 * t232;
t262 = t230 * t68 - t233 * t67;
t260 = -pkin(1) + t269;
t239 = pkin(3) * t250 - qJ(4) * t152 + t54;
t206 = t222 * t277;
t205 = (t225 - t226) * t299;
t204 = -t280 - t330;
t184 = -t281 - t330;
t178 = t227 * t244 + t322;
t174 = t203 - t206;
t173 = t203 + t206;
t172 = -t207 + t202;
t134 = -t162 + t183;
t133 = t161 - t183;
t127 = t162 - t161;
t96 = -qJD(5) * t165 + t272;
t91 = (-t163 * t232 + t165 * t229) * t186;
t90 = (t163 * t229 + t165 * t232) * t186;
t86 = t125 * t230 - t354;
t85 = t122 * t230 - t354;
t83 = -t126 * t233 - t367;
t82 = t124 * t233 - t367;
t81 = -t122 * t233 - t356;
t73 = -t165 * t305 + t232 * t97;
t72 = -t165 * t304 - t229 * t97;
t71 = -t163 * t304 + t229 * t96;
t70 = -t163 * t305 - t232 * t96;
t69 = t147 * t230 + t233 * t90;
t65 = t133 * t232 - t319;
t64 = -t134 * t229 + t355;
t63 = -t133 * t229 - t317;
t62 = -t134 * t232 - t357;
t61 = -t229 * t284 - t317;
t59 = t107 * t232 - t357;
t57 = t233 * t72 + t283;
t56 = t233 * t70 - t283;
t55 = t273 * t194 + t237;
t53 = pkin(2) * t126 + t315 - t373;
t51 = -pkin(2) * t119 - t314 + t383;
t50 = qJ(4) * t132 + t54;
t49 = t229 * t79 + t232 * t252;
t48 = -t229 * t259 - t232 * t75;
t46 = t229 * t75 - t232 * t259;
t45 = pkin(3) * t132 + t52;
t44 = (t119 - t302) * pkin(3) + t360;
t43 = pkin(3) * t302 - qJ(4) * t126 - t360;
t42 = t230 * t79 + t233 * t62;
t41 = t230 * t252 + t233 * t63;
t38 = t230 * t60 + t233 * t259;
t37 = t230 * t259 - t233 * t60;
t36 = t230 * t58 + t233 * t75;
t35 = t230 * t75 - t233 * t58;
t30 = t127 * t230 + t233 * t46;
t29 = t230 * t47 + t233 * t98;
t28 = t230 * t98 - t233 * t47;
t27 = -pkin(2) * t109 + pkin(8) * t40;
t26 = t230 * t54 + t233 * t52;
t25 = t230 * t52 - t233 * t54;
t24 = pkin(8) * t86 + t368 + t40;
t23 = t373 + t230 * t43 - (pkin(2) + t324) * t126;
t22 = t274 * t119 + t233 * t44 - t383;
t21 = pkin(4) * t47 - qJ(4) * t49;
t17 = pkin(8) * t85 + t230 * t50 + t233 * t45 + t368;
t16 = pkin(4) * t259 - t328 * t61 - t320;
t15 = pkin(4) * t75 - t328 * t59 + t318;
t14 = pkin(4) * t60 - qJ(4) * t61 - t20;
t13 = pkin(4) * t58 - qJ(4) * t59 - t19;
t12 = pkin(8) * t26 + (-t274 - t324) * t55;
t9 = t10 * t230 + t233 * t33;
t8 = -t10 * t233 + t230 * t33;
t7 = pkin(4) * t98 - t328 * t49 - t11;
t6 = -pkin(2) * t61 + pkin(8) * t38 + t14 * t230 + t16 * t233;
t5 = -pkin(2) * t59 + pkin(8) * t36 + t13 * t230 + t15 * t233;
t4 = pkin(4) * t10 - qJ(4) * t11;
t3 = pkin(4) * t33 - t328 * t11;
t2 = -pkin(2) * t49 + pkin(8) * t29 + t21 * t230 + t233 * t7;
t1 = -pkin(2) * t11 + pkin(8) * t9 + t230 * t4 + t233 * t3;
t18 = [0, 0, 0, 0, 0, qJDD(1), t254, t255, 0, 0, (t202 * t227 + t234 * t351) * t231, t228 * t205 + (t231 * t174 + t234 * t278) * t227, t228 * t172 + (t294 + t234 * (-t281 + t330)) * t227, (t203 * t227 - t231 * t351) * t234, t228 * t173 + (t231 * (t280 - t330) + t292) * t227, t228 * t270, (-t153 + pkin(1) * (t200 * t234 + t204 * t231)) * t228 + (t234 * t178 + pkin(1) * t174 + pkin(7) * (t204 * t234 - t294)) * t227, -t178 * t298 - t228 * t154 + pkin(1) * (-t227 * t278 + (t234 * t184 - t231 * t199) * t228) + (-t184 * t231 - t292) * t323, pkin(1) * ((-t234 * t172 + t231 * t173) * t228 - (-t225 - t226) * t224 * t296) + (t172 * t231 + t173 * t234) * t323 + t344 * t227, pkin(1) * (t227 * t178 + (-t153 * t234 + t154 * t231) * t228) + t344 * t323, t334, t228 * t82 + (t231 * (-t124 * t230 - t313) - t353) * t227, t381 + (-t234 * t125 - t372) * t227, t335, t377, t333, (t51 + t388) * t228 + (t231 * (t315 - t384) + t234 * (t67 - t385) + t387) * t227, (t53 - t382) * t228 + (t231 * (t314 - t374) + t234 * (t68 - t375) + t380) * t227, (t24 + pkin(1) * (t231 * t86 + t361)) * t228 + (-t231 * t262 + pkin(7) * (t234 * t86 - t362) + t260 * (-t125 * t233 - t356)) * t227, (t27 + pkin(1) * (-t109 * t234 + t231 * t40)) * t228 + (pkin(7) * (t231 * t109 + t234 * t40) + t260 * t262) * t227, t333, -t381 + (t234 * t122 + t372) * t227, -t377, t334, t228 * t83 + (t231 * (t126 * t230 - t313) - t353) * t227, t335, (t17 + pkin(1) * (t231 * t85 + t361)) * t228 + (t231 * (-pkin(8) * t81 - t230 * t45 + t233 * t50) + t234 * (-pkin(2) * t81 - t265) - pkin(1) * t81 + pkin(7) * (t234 * t85 - t362)) * t227, (t22 - t388) * t228 + (t231 * (qJ(4) * t313 - t230 * t44 + t384) + t234 * (-t239 + t385) - t387) * t227, (t23 + t382) * t228 + (t231 * (t126 * t325 + t233 * t43 + t374) + t234 * (-t336 + t375) - t380) * t227, (t12 + pkin(1) * (t231 * t26 - t234 * t55)) * t228 + (t231 * (-pkin(8) * t25 + (-qJ(4) * t233 + t325) * t55) + t234 * (-pkin(2) * t25 - t268) - pkin(1) * t25 + pkin(7) * (t231 * t55 + t234 * t26)) * t227, t228 * t57 + (t231 * (-t230 * t72 + t282) - t234 * t73) * t227, t228 * t30 + (t231 * (t127 * t233 - t230 * t46) - t234 * t48) * t227, t228 * t42 + (t231 * (-t230 * t62 + t233 * t79) - t234 * t64) * t227, t228 * t56 + (t231 * (-t230 * t70 - t282) + t234 * t71) * t227, t228 * t41 + (t231 * (-t230 * t63 + t233 * t252) - t234 * t65) * t227, t228 * t69 + (t231 * (t147 * t233 - t230 * t90) - t234 * t91) * t227, (t5 + pkin(1) * (t231 * t36 - t234 * t59)) * t228 + (t231 * (-pkin(8) * t35 + t13 * t233 - t15 * t230) + t234 * (-pkin(2) * t35 - t337) - pkin(1) * t35 + pkin(7) * (t231 * t59 + t234 * t36)) * t227, (t6 + pkin(1) * (t231 * t38 - t234 * t61)) * t228 + (t231 * (-pkin(8) * t37 + t14 * t233 - t16 * t230) + t234 * (-pkin(2) * t37 - t338) - pkin(1) * t37 + pkin(7) * (t231 * t61 + t234 * t38)) * t227, (t2 + pkin(1) * (t231 * t29 - t234 * t49)) * t228 + (t231 * (-pkin(8) * t28 + t21 * t233 - t230 * t7) + t234 * (-pkin(2) * t28 - t339) - pkin(1) * t28 + pkin(7) * (t231 * t49 + t234 * t29)) * t227, (t1 + pkin(1) * (-t11 * t234 + t231 * t9)) * t228 + (t231 * (-pkin(8) * t8 - t230 * t3 + t233 * t4) + t234 * (-pkin(2) * t8 - t340) - pkin(1) * t8 + pkin(7) * (t11 * t231 + t234 * t9)) * t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t214, t205, t172, t214, t173, t270, -t153, -t154, 0, 0, t267, t82, t104, t266, t105, t256, t51, t53, t24, t27, t256, -t104, -t105, t267, t83, t266, t17, t22, t23, t12, t57, t30, t42, t56, t41, t69, t5, t6, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t303, t342, t125, -t303, -t345, -t197, -t67, -t68, 0, 0, -t197, -t122, t345, t303, t342, -t303, t265, t239, t336, t268, t73, t48, t64, -t71, t65, t91, t337, t338, t339, t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, -t250, t343, t54, 0, 0, 0, 0, 0, 0, t58, t60, t47, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t127, t79, -t128, t252, t147, -t19, -t20, 0, 0;];
tauJ_reg = t18;

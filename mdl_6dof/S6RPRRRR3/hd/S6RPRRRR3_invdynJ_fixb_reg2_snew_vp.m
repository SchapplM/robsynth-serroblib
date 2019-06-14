% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:59:14
% EndTime: 2019-05-06 02:59:32
% DurationCPUTime: 7.48s
% Computational Cost: add. (53684->511), mult. (106501->721), div. (0->0), fcn. (75211->12), ass. (0->310)
t278 = sin(qJ(6));
t280 = sin(qJ(4));
t285 = cos(qJ(4));
t281 = sin(qJ(3));
t326 = qJD(1) * t281;
t238 = -t285 * qJD(3) + t280 * t326;
t239 = qJD(3) * t280 + t285 * t326;
t279 = sin(qJ(5));
t284 = cos(qJ(5));
t218 = t284 * t238 + t239 * t279;
t220 = -t238 * t279 + t239 * t284;
t283 = cos(qJ(6));
t183 = t283 * t218 + t220 * t278;
t185 = -t218 * t278 + t220 * t283;
t146 = t185 * t183;
t321 = qJD(1) * qJD(3);
t265 = t281 * t321;
t286 = cos(qJ(3));
t319 = t286 * qJDD(1);
t245 = -t265 + t319;
t237 = -qJDD(4) + t245;
t234 = -qJDD(5) + t237;
t229 = -qJDD(6) + t234;
t362 = -t146 - t229;
t368 = t278 * t362;
t189 = t220 * t218;
t360 = -t189 - t234;
t367 = t279 * t360;
t366 = t283 * t362;
t365 = t284 * t360;
t288 = qJD(1) ^ 2;
t337 = t239 * t238;
t291 = -t237 - t337;
t364 = t280 * t291;
t363 = t285 * t291;
t312 = t286 * t321;
t320 = t281 * qJDD(1);
t244 = t312 + t320;
t297 = -t280 * qJDD(3) - t285 * t244;
t214 = -qJD(4) * t238 - t297;
t298 = t285 * qJDD(3) - t280 * t244;
t289 = -qJD(4) * t239 + t298;
t159 = -t218 * qJD(5) + t284 * t214 + t279 * t289;
t308 = t279 * t214 - t284 * t289;
t294 = qJD(5) * t220 + t308;
t105 = -t183 * qJD(6) + t283 * t159 - t278 * t294;
t261 = qJD(1) * t286 - qJD(4);
t255 = -qJD(5) + t261;
t250 = -qJD(6) + t255;
t167 = t183 * t250;
t361 = t167 + t105;
t203 = t218 * t255;
t143 = -t203 + t159;
t359 = t203 + t159;
t228 = t238 * t261;
t194 = t214 - t228;
t309 = t278 * t159 + t283 * t294;
t84 = (qJD(6) + t250) * t185 + t309;
t140 = (qJD(5) + t255) * t220 + t308;
t190 = (qJD(4) + t261) * t239 - t298;
t181 = t183 ^ 2;
t182 = t185 ^ 2;
t358 = t218 ^ 2;
t217 = t220 ^ 2;
t357 = t238 ^ 2;
t236 = t239 ^ 2;
t248 = t250 ^ 2;
t254 = t255 ^ 2;
t259 = t261 ^ 2;
t356 = qJD(3) ^ 2;
t282 = sin(qJ(1));
t287 = cos(qJ(1));
t311 = t282 * g(1) - g(2) * t287;
t240 = qJDD(1) * pkin(1) + t311;
t302 = g(1) * t287 + g(2) * t282;
t241 = -pkin(1) * t288 - t302;
t274 = sin(pkin(11));
t275 = cos(pkin(11));
t307 = t275 * t240 - t274 * t241;
t210 = -qJDD(1) * pkin(2) - t288 * pkin(7) - t307;
t300 = -t245 + t265;
t301 = t244 + t312;
t176 = t300 * pkin(3) - t301 * pkin(8) + t210;
t327 = t274 * t240 + t275 * t241;
t211 = -pkin(2) * t288 + qJDD(1) * pkin(7) + t327;
t303 = -pkin(3) * t286 - pkin(8) * t281;
t306 = t288 * t303 + t211;
t328 = -g(3) + qJDD(2);
t310 = t281 * t328;
t187 = -t356 * pkin(3) + qJDD(3) * pkin(8) + t306 * t286 + t310;
t136 = -t285 * t176 + t280 * t187;
t114 = t291 * pkin(4) - t194 * pkin(9) - t136;
t137 = t280 * t176 + t285 * t187;
t225 = -pkin(4) * t261 - pkin(9) * t239;
t120 = -t357 * pkin(4) + pkin(9) * t289 + t261 * t225 + t137;
t72 = -t284 * t114 + t279 * t120;
t73 = t279 * t114 + t284 * t120;
t44 = t279 * t73 - t284 * t72;
t355 = pkin(4) * t44;
t101 = -t140 * t279 - t143 * t284;
t354 = pkin(4) * t101;
t57 = pkin(5) * t360 - t143 * pkin(10) - t72;
t200 = -pkin(5) * t255 - pkin(10) * t220;
t59 = -t358 * pkin(5) - t294 * pkin(10) + t255 * t200 + t73;
t33 = t278 * t59 - t283 * t57;
t34 = t278 * t57 + t283 * t59;
t17 = t278 * t34 - t283 * t33;
t353 = t17 * t279;
t352 = t17 * t284;
t264 = t286 * t328;
t186 = -qJDD(3) * pkin(3) - t356 * pkin(8) + t306 * t281 - t264;
t135 = -t289 * pkin(4) - t357 * pkin(9) + t239 * t225 + t186;
t89 = t294 * pkin(5) - t358 * pkin(10) + t220 * t200 + t135;
t351 = t278 * t89;
t350 = t280 * t44;
t349 = t283 * t89;
t348 = t285 * t44;
t127 = -t146 + t229;
t347 = t127 * t278;
t346 = t127 * t283;
t345 = t135 * t279;
t344 = t135 * t284;
t171 = -t189 + t234;
t343 = t171 * t279;
t342 = t171 * t284;
t341 = t186 * t280;
t340 = t186 * t285;
t206 = t237 - t337;
t339 = t206 * t280;
t338 = t206 * t285;
t336 = t250 * t278;
t335 = t250 * t283;
t334 = t255 * t279;
t333 = t255 * t284;
t332 = t261 * t280;
t331 = t261 * t285;
t260 = t286 * t288 * t281;
t253 = qJDD(3) + t260;
t330 = t281 * t253;
t252 = -t260 + qJDD(3);
t329 = t286 * t252;
t325 = qJD(4) - t261;
t16 = pkin(5) * t17;
t18 = t278 * t33 + t283 * t34;
t8 = t18 * t279 + t352;
t318 = pkin(4) * t8 + t16;
t317 = t286 * t146;
t316 = t286 * t189;
t315 = t286 * t337;
t314 = pkin(1) * t274 + pkin(7);
t87 = -t167 + t105;
t52 = -t278 * t84 - t283 * t87;
t54 = t278 * t87 - t283 * t84;
t28 = t279 * t54 + t284 * t52;
t50 = pkin(5) * t52;
t313 = pkin(4) * t28 + t50;
t45 = t279 * t72 + t284 * t73;
t99 = t136 * t280 + t285 * t137;
t198 = t211 * t281 - t264;
t199 = t286 * t211 + t310;
t160 = t281 * t198 + t286 * t199;
t138 = -t248 - t181;
t96 = t138 * t278 + t366;
t305 = pkin(5) * t96 - t33;
t197 = -t217 - t254;
t148 = t197 * t284 + t343;
t304 = pkin(4) * t148 - t73;
t299 = t136 * t285 - t137 * t280;
t177 = -t254 - t358;
t125 = t177 * t279 + t365;
t296 = pkin(4) * t125 - t72;
t162 = -t182 - t248;
t110 = t162 * t283 + t347;
t295 = pkin(5) * t110 - t34;
t97 = t138 * t283 - t368;
t60 = t279 * t97 + t284 * t96;
t293 = pkin(4) * t60 + t305;
t292 = -pkin(1) * t275 - pkin(2) + t303;
t111 = -t162 * t278 + t346;
t65 = t110 * t284 + t111 * t279;
t290 = pkin(4) * t65 + t295;
t271 = t286 ^ 2;
t270 = t281 ^ 2;
t268 = t271 * t288;
t266 = t270 * t288;
t258 = -t268 - t356;
t257 = -t266 - t356;
t249 = t266 + t268;
t247 = (t270 + t271) * qJDD(1);
t246 = -0.2e1 * t265 + t319;
t243 = 0.2e1 * t312 + t320;
t227 = -t236 + t259;
t226 = -t259 + t357;
t224 = -t257 * t281 - t329;
t223 = t258 * t286 - t330;
t222 = t236 - t357;
t221 = -t236 - t259;
t215 = -t259 - t357;
t205 = t236 + t357;
t202 = -t217 + t254;
t201 = -t254 + t358;
t195 = t325 * t238 + t297;
t193 = t214 + t228;
t191 = -t325 * t239 + t298;
t188 = t217 - t358;
t179 = -t221 * t280 + t338;
t178 = t221 * t285 + t339;
t175 = t215 * t285 - t364;
t174 = t215 * t280 + t363;
t166 = -t182 + t248;
t165 = t181 - t248;
t164 = (t218 * t284 - t220 * t279) * t255;
t163 = (t218 * t279 + t220 * t284) * t255;
t161 = -t217 - t358;
t156 = -t190 * t285 + t194 * t280;
t154 = t201 * t284 + t343;
t153 = -t202 * t279 + t365;
t152 = t201 * t279 - t342;
t151 = t202 * t284 + t367;
t150 = t179 * t286 - t195 * t281;
t149 = -t197 * t279 + t342;
t147 = t175 * t286 - t191 * t281;
t145 = t182 - t181;
t139 = (qJD(5) - t255) * t220 + t308;
t133 = t159 * t284 + t220 * t334;
t132 = t159 * t279 - t220 * t333;
t131 = -t218 * t333 + t279 * t294;
t130 = -t218 * t334 - t284 * t294;
t126 = t177 * t284 - t367;
t123 = (t183 * t283 - t185 * t278) * t250;
t122 = (t183 * t278 + t185 * t283) * t250;
t121 = -t181 - t182;
t118 = t165 * t283 + t347;
t117 = -t166 * t278 + t366;
t116 = t165 * t278 - t346;
t115 = t166 * t283 + t368;
t108 = -t148 * t280 + t149 * t285;
t107 = t148 * t285 + t149 * t280;
t106 = -pkin(9) * t148 + t344;
t104 = -qJD(6) * t185 - t309;
t103 = -t140 * t284 + t143 * t279;
t102 = -t139 * t284 - t279 * t359;
t100 = -t139 * t279 + t284 * t359;
t94 = -pkin(9) * t125 + t345;
t93 = -t125 * t280 + t126 * t285;
t92 = t125 * t285 + t126 * t280;
t91 = -t122 * t279 + t123 * t284;
t90 = t122 * t284 + t123 * t279;
t83 = (qJD(6) - t250) * t185 + t309;
t82 = t105 * t283 + t185 * t336;
t81 = t105 * t278 - t185 * t335;
t80 = -t104 * t278 - t183 * t335;
t79 = t104 * t283 - t183 * t336;
t77 = t108 * t286 + t281 * t359;
t76 = -pkin(4) * t359 + pkin(9) * t149 + t345;
t75 = -pkin(4) * t139 + pkin(9) * t126 - t344;
t74 = t139 * t281 + t286 * t93;
t70 = -t116 * t279 + t118 * t284;
t69 = -t115 * t279 + t117 * t284;
t68 = t116 * t284 + t118 * t279;
t67 = t115 * t284 + t117 * t279;
t66 = -t110 * t279 + t111 * t284;
t64 = -t101 * t280 + t103 * t285;
t63 = t101 * t285 + t103 * t280;
t62 = -pkin(10) * t110 + t349;
t61 = -t279 * t96 + t284 * t97;
t58 = -pkin(10) * t96 + t351;
t56 = t161 * t281 + t286 * t64;
t53 = -t278 * t361 - t283 * t83;
t51 = -t278 * t83 + t283 * t361;
t49 = -t279 * t81 + t284 * t82;
t48 = -t279 * t79 + t284 * t80;
t47 = t279 * t82 + t284 * t81;
t46 = t279 * t80 + t284 * t79;
t43 = -pkin(5) * t361 + pkin(10) * t111 + t351;
t42 = -t280 * t65 + t285 * t66;
t41 = t280 * t66 + t285 * t65;
t40 = -pkin(5) * t83 + pkin(10) * t97 - t349;
t39 = -pkin(4) * t135 + pkin(9) * t45;
t38 = -t280 * t60 + t285 * t61;
t37 = t280 * t61 + t285 * t60;
t36 = -pkin(9) * t101 - t44;
t35 = -pkin(4) * t161 + pkin(9) * t103 + t45;
t31 = t281 * t361 + t286 * t42;
t30 = -t279 * t52 + t284 * t54;
t29 = -t279 * t51 + t284 * t53;
t27 = t279 * t53 + t284 * t51;
t26 = t281 * t83 + t286 * t38;
t25 = t285 * t45 - t350;
t24 = t280 * t45 + t348;
t23 = t135 * t281 + t25 * t286;
t22 = -pkin(9) * t65 - t279 * t43 + t284 * t62;
t21 = -pkin(9) * t60 - t279 * t40 + t284 * t58;
t20 = -pkin(4) * t361 + pkin(9) * t66 + t279 * t62 + t284 * t43;
t19 = -pkin(4) * t83 + pkin(9) * t61 + t279 * t58 + t284 * t40;
t15 = -t28 * t280 + t285 * t30;
t14 = t28 * t285 + t280 * t30;
t13 = -pkin(5) * t89 + pkin(10) * t18;
t12 = t121 * t281 + t15 * t286;
t11 = -pkin(10) * t52 - t17;
t10 = -pkin(5) * t121 + pkin(10) * t54 + t18;
t9 = t18 * t284 - t353;
t7 = -pkin(9) * t28 - t10 * t279 + t11 * t284;
t6 = -pkin(4) * t121 + pkin(9) * t30 + t10 * t284 + t11 * t279;
t5 = -t280 * t8 + t285 * t9;
t4 = t280 * t9 + t285 * t8;
t3 = -pkin(9) * t8 - pkin(10) * t352 - t13 * t279;
t2 = t281 * t89 + t286 * t5;
t1 = -pkin(4) * t89 + pkin(9) * t9 - pkin(10) * t353 + t13 * t284;
t32 = [0, 0, 0, 0, 0, qJDD(1), t311, t302, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t275 - t274 * t288) + t307, pkin(1) * (-qJDD(1) * t274 - t275 * t288) - t327, 0, pkin(1) * (t274 * t327 + t275 * t307), t301 * t281, t243 * t286 + t246 * t281, t330 + t286 * (-t266 + t356), -t300 * t286, t281 * (t268 - t356) + t329, 0, -t286 * t210 + pkin(2) * t246 + pkin(7) * t223 + pkin(1) * (t223 * t274 + t246 * t275), t281 * t210 - pkin(2) * t243 + pkin(7) * t224 + pkin(1) * (t224 * t274 - t243 * t275), pkin(2) * t249 + pkin(7) * t247 + pkin(1) * (t247 * t274 + t249 * t275) + t160, -pkin(2) * t210 + pkin(7) * t160 + pkin(1) * (t160 * t274 - t210 * t275), t281 * (t214 * t285 + t239 * t332) - t315, t281 * (t191 * t285 - t193 * t280) - t286 * t222, t281 * (-t227 * t280 + t363) - t286 * t194, t281 * (-t238 * t331 - t280 * t289) + t315, t281 * (t226 * t285 + t339) + t286 * t190, t286 * t237 + t281 * (t238 * t285 - t239 * t280) * t261, t281 * (-pkin(8) * t174 + t341) + t286 * (-pkin(3) * t174 + t136) - pkin(2) * t174 + pkin(7) * t147 + pkin(1) * (t147 * t274 - t174 * t275), t281 * (-pkin(8) * t178 + t340) + t286 * (-pkin(3) * t178 + t137) - pkin(2) * t178 + pkin(7) * t150 + pkin(1) * (t150 * t274 - t178 * t275), t281 * t299 + t314 * (t156 * t286 - t205 * t281) + t292 * (-t190 * t280 - t194 * t285), t314 * (t186 * t281 + t286 * t99) - t292 * t299, t281 * (-t132 * t280 + t133 * t285) - t316, t281 * (-t100 * t280 + t102 * t285) - t286 * t188, t281 * (-t151 * t280 + t153 * t285) - t286 * t143, t281 * (-t130 * t280 + t131 * t285) + t316, t281 * (-t152 * t280 + t154 * t285) + t286 * t140, t281 * (-t163 * t280 + t164 * t285) + t286 * t234, t281 * (-pkin(8) * t92 - t280 * t75 + t285 * t94) + t286 * (-pkin(3) * t92 - t296) - pkin(2) * t92 + pkin(7) * t74 + pkin(1) * (t274 * t74 - t275 * t92), t281 * (-pkin(8) * t107 + t106 * t285 - t280 * t76) + t286 * (-pkin(3) * t107 - t304) - pkin(2) * t107 + pkin(7) * t77 + pkin(1) * (-t107 * t275 + t274 * t77), t281 * (-pkin(8) * t63 - t280 * t35 + t285 * t36) + t286 * (-pkin(3) * t63 - t354) - pkin(2) * t63 + pkin(7) * t56 + pkin(1) * (t274 * t56 - t275 * t63), t281 * (-pkin(8) * t24 - pkin(9) * t348 - t280 * t39) + t286 * (-pkin(3) * t24 - t355) - pkin(2) * t24 + pkin(7) * t23 + pkin(1) * (t23 * t274 - t24 * t275), t281 * (-t280 * t47 + t285 * t49) - t317, t281 * (-t27 * t280 + t285 * t29) - t286 * t145, t281 * (-t280 * t67 + t285 * t69) - t286 * t87, t281 * (-t280 * t46 + t285 * t48) + t317, t281 * (-t280 * t68 + t285 * t70) + t286 * t84, t281 * (-t280 * t90 + t285 * t91) + t286 * t229, t281 * (-pkin(8) * t37 - t19 * t280 + t21 * t285) + t286 * (-pkin(3) * t37 - t293) - pkin(2) * t37 + pkin(7) * t26 + pkin(1) * (t26 * t274 - t275 * t37), t281 * (-pkin(8) * t41 - t20 * t280 + t22 * t285) + t286 * (-pkin(3) * t41 - t290) - pkin(2) * t41 + pkin(7) * t31 + pkin(1) * (t274 * t31 - t275 * t41), t281 * (-pkin(8) * t14 - t280 * t6 + t285 * t7) + t286 * (-pkin(3) * t14 - t313) - pkin(2) * t14 + pkin(7) * t12 + pkin(1) * (t12 * t274 - t14 * t275), t281 * (-pkin(8) * t4 - t1 * t280 + t285 * t3) + t286 * (-pkin(3) * t4 - t318) - pkin(2) * t4 + pkin(7) * t2 + pkin(1) * (t2 * t274 - t275 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t328, 0, 0, 0, 0, 0, 0, t253 * t286 + t258 * t281, -t252 * t281 + t257 * t286, 0, -t198 * t286 + t199 * t281, 0, 0, 0, 0, 0, 0, t175 * t281 + t191 * t286, t179 * t281 + t195 * t286, t156 * t281 + t205 * t286, -t186 * t286 + t281 * t99, 0, 0, 0, 0, 0, 0, -t139 * t286 + t281 * t93, t108 * t281 - t286 * t359, -t161 * t286 + t281 * t64, -t135 * t286 + t25 * t281, 0, 0, 0, 0, 0, 0, t281 * t38 - t286 * t83, t281 * t42 - t286 * t361, -t121 * t286 + t15 * t281, t281 * t5 - t286 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t260, t266 - t268, t320, t260, t319, qJDD(3), -t198, -t199, 0, 0, t214 * t280 - t239 * t331, t191 * t280 + t193 * t285, t227 * t285 + t364, -t238 * t332 + t285 * t289, t226 * t280 - t338, (t238 * t280 + t239 * t285) * t261, pkin(3) * t191 + pkin(8) * t175 - t340, pkin(3) * t195 + pkin(8) * t179 + t341, pkin(3) * t205 + pkin(8) * t156 + t99, -pkin(3) * t186 + pkin(8) * t99, t132 * t285 + t133 * t280, t100 * t285 + t102 * t280, t151 * t285 + t153 * t280, t130 * t285 + t131 * t280, t152 * t285 + t154 * t280, t163 * t285 + t164 * t280, -pkin(3) * t139 + pkin(8) * t93 + t280 * t94 + t285 * t75, -pkin(3) * t359 + pkin(8) * t108 + t106 * t280 + t285 * t76, -pkin(3) * t161 + pkin(8) * t64 + t280 * t36 + t285 * t35, -pkin(3) * t135 + pkin(8) * t25 - pkin(9) * t350 + t285 * t39, t280 * t49 + t285 * t47, t27 * t285 + t280 * t29, t280 * t69 + t285 * t67, t280 * t48 + t285 * t46, t280 * t70 + t285 * t68, t280 * t91 + t285 * t90, -pkin(3) * t83 + pkin(8) * t38 + t19 * t285 + t21 * t280, -pkin(3) * t361 + pkin(8) * t42 + t20 * t285 + t22 * t280, -pkin(3) * t121 + pkin(8) * t15 + t280 * t7 + t285 * t6, -pkin(3) * t89 + pkin(8) * t5 + t1 * t285 + t280 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t337, t222, t194, -t337, -t190, -t237, -t136, -t137, 0, 0, t189, t188, t143, -t189, -t140, -t234, t296, t304, t354, t355, t146, t145, t87, -t146, -t84, -t229, t293, t290, t313, t318; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t188, t143, -t189, -t140, -t234, -t72, -t73, 0, 0, t146, t145, t87, -t146, -t84, -t229, t305, t295, t50, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t145, t87, -t146, -t84, -t229, -t33, -t34, 0, 0;];
tauJ_reg  = t32;

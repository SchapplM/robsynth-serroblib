% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 01:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:15:10
% EndTime: 2019-05-06 01:15:24
% DurationCPUTime: 5.06s
% Computational Cost: add. (21836->402), mult. (42908->522), div. (0->0), fcn. (29424->10), ass. (0->258)
t244 = sin(qJ(4));
t248 = cos(qJ(4));
t245 = sin(qJ(3));
t297 = qJD(1) * t245;
t207 = qJD(3) * t244 + t248 * t297;
t243 = sin(qJ(5));
t247 = cos(qJ(5));
t271 = qJD(3) * t248 - t244 * t297;
t187 = t207 * t243 - t247 * t271;
t189 = t207 * t247 + t243 * t271;
t150 = t189 * t187;
t292 = qJD(1) * qJD(3);
t231 = t245 * t292;
t249 = cos(qJ(3));
t290 = t249 * qJDD(1);
t213 = -t231 + t290;
t206 = -qJDD(4) + t213;
t203 = -qJDD(5) + t206;
t338 = t150 + t203;
t344 = pkin(5) * t338;
t278 = t249 * t292;
t291 = t245 * qJDD(1);
t212 = t278 + t291;
t181 = qJD(4) * t271 + t244 * qJDD(3) + t248 * t212;
t260 = t248 * qJDD(3) - t244 * t212;
t255 = -qJD(4) * t207 + t260;
t127 = -qJD(5) * t187 + t181 * t247 + t243 * t255;
t227 = qJD(1) * t249 - qJD(4);
t220 = -qJD(5) + t227;
t167 = t187 * t220;
t109 = -t167 + t127;
t339 = qJ(6) * t109;
t251 = qJD(1) ^ 2;
t246 = sin(qJ(1));
t250 = cos(qJ(1));
t276 = g(1) * t246 - g(2) * t250;
t208 = qJDD(1) * pkin(1) + t276;
t267 = g(1) * t250 + g(2) * t246;
t209 = -pkin(1) * t251 - t267;
t239 = sin(pkin(10));
t240 = cos(pkin(10));
t273 = t208 * t240 - t239 * t209;
t174 = -qJDD(1) * pkin(2) - pkin(7) * t251 - t273;
t264 = -t213 + t231;
t265 = t212 + t278;
t141 = pkin(3) * t264 - pkin(8) * t265 + t174;
t298 = t208 * t239 + t209 * t240;
t175 = -pkin(2) * t251 + qJDD(1) * pkin(7) + t298;
t268 = -pkin(3) * t249 - pkin(8) * t245;
t272 = t251 * t268 + t175;
t299 = -g(3) + qJDD(2);
t275 = t245 * t299;
t333 = qJD(3) ^ 2;
t147 = -pkin(3) * t333 + qJDD(3) * pkin(8) + t249 * t272 + t275;
t100 = -t141 * t248 + t244 * t147;
t197 = t271 * t227;
t155 = t197 + t181;
t261 = t271 * t207;
t337 = -t206 + t261;
t83 = pkin(4) * t337 - pkin(9) * t155 - t100;
t101 = t141 * t244 + t147 * t248;
t194 = -pkin(4) * t227 - pkin(9) * t207;
t269 = t271 ^ 2;
t86 = -pkin(4) * t269 + pkin(9) * t255 + t227 * t194 + t101;
t42 = t243 * t86 - t247 * t83;
t254 = 0.2e1 * qJD(6) * t189 + t339 + t344 + t42;
t253 = -t254 - t344;
t185 = t187 ^ 2;
t219 = t220 ^ 2;
t142 = -t219 - t185;
t304 = t247 * t338;
t92 = t142 * t243 - t304;
t91 = pkin(4) * t92;
t343 = t253 + t91;
t342 = t244 * t337;
t341 = t248 * t337;
t186 = t189 ^ 2;
t160 = -t186 - t219;
t136 = -t150 + t203;
t311 = t243 * t136;
t114 = t160 * t247 + t311;
t113 = pkin(4) * t114;
t274 = t243 * t181 - t247 * t255;
t126 = -qJD(5) * t189 - t274;
t163 = -pkin(5) * t220 - qJ(6) * t189;
t43 = t243 * t83 + t247 * t86;
t33 = -pkin(5) * t185 + qJ(6) * t126 - 0.2e1 * qJD(6) * t187 + t163 * t220 + t43;
t258 = pkin(5) * t160 - t33;
t340 = t113 + t258;
t230 = t249 * t299;
t146 = -qJDD(3) * pkin(3) - pkin(8) * t333 + t245 * t272 - t230;
t99 = -pkin(4) * t255 - pkin(9) * t269 + t207 * t194 + t146;
t310 = t243 * t338;
t58 = -pkin(5) * t126 - qJ(6) * t185 + t189 * t163 + qJDD(6) + t99;
t335 = t167 + t127;
t156 = t197 - t181;
t106 = (qJD(5) + t220) * t189 + t274;
t151 = (qJD(4) + t227) * t207 - t260;
t205 = t207 ^ 2;
t224 = t227 ^ 2;
t93 = t142 * t247 + t310;
t61 = t244 * t93 + t248 * t92;
t332 = pkin(3) * t61;
t305 = t247 * t136;
t115 = -t160 * t243 + t305;
t77 = t114 * t248 + t115 * t244;
t331 = pkin(3) * t77;
t19 = t243 * t43 - t247 * t42;
t330 = pkin(4) * t19;
t71 = -t106 * t243 - t109 * t247;
t73 = -t106 * t247 + t109 * t243;
t38 = t244 * t73 + t248 * t71;
t329 = pkin(8) * t38;
t328 = pkin(8) * t61;
t327 = pkin(8) * t77;
t326 = pkin(9) * t71;
t325 = pkin(9) * t92;
t324 = pkin(9) * t114;
t323 = t243 * t254;
t322 = t243 * t99;
t321 = t244 * t19;
t320 = t247 * t254;
t319 = t247 * t99;
t318 = t248 * t19;
t130 = -t185 - t186;
t39 = -t244 * t71 + t248 * t73;
t317 = -pkin(3) * t130 + pkin(8) * t39;
t105 = (qJD(5) - t220) * t189 + t274;
t62 = -t244 * t92 + t248 * t93;
t316 = -pkin(3) * t105 + pkin(8) * t62;
t78 = -t114 * t244 + t115 * t248;
t315 = -pkin(3) * t335 + pkin(8) * t78;
t313 = t220 * t243;
t312 = t220 * t247;
t309 = t244 * t146;
t170 = t206 + t261;
t308 = t244 * t170;
t307 = t244 * t207;
t226 = t249 * t251 * t245;
t218 = qJDD(3) + t226;
t306 = t245 * t218;
t303 = t248 * t146;
t302 = t248 * t170;
t301 = t248 * t207;
t217 = -t226 + qJDD(3);
t300 = t249 * t217;
t31 = t130 * t245 + t249 * t39;
t289 = pkin(1) * (t239 * t31 - t240 * t38) + pkin(7) * t31 - pkin(2) * t38;
t46 = t105 * t245 + t249 * t62;
t288 = pkin(1) * (t239 * t46 - t240 * t61) + pkin(7) * t46 - pkin(2) * t61;
t51 = t245 * t335 + t249 * t78;
t287 = pkin(1) * (t239 * t51 - t240 * t77) + pkin(7) * t51 - pkin(2) * t77;
t286 = t113 - t43;
t285 = t249 * t150;
t283 = pkin(1) * t239 + pkin(7);
t69 = pkin(4) * t71;
t282 = -pkin(3) * t38 - t69;
t11 = t243 * t33 - t320;
t27 = pkin(5) * t254;
t281 = pkin(4) * t11 - t27;
t280 = -pkin(4) * t105 + pkin(9) * t93;
t279 = -pkin(4) * t130 + pkin(9) * t73;
t277 = -pkin(4) * t335 + pkin(9) * t115;
t20 = t243 * t42 + t247 * t43;
t67 = t100 * t244 + t101 * t248;
t161 = t175 * t245 - t230;
t162 = t249 * t175 + t275;
t128 = t245 * t161 + t162 * t249;
t270 = t42 - t91;
t263 = t100 * t248 - t101 * t244;
t259 = t249 * t261;
t256 = -pkin(1) * t240 - pkin(2) + t268;
t236 = t249 ^ 2;
t235 = t245 ^ 2;
t233 = t236 * t251;
t232 = t235 * t251;
t223 = -t233 - t333;
t222 = -t232 - t333;
t216 = t232 + t233;
t215 = (t235 + t236) * qJDD(1);
t214 = -0.2e1 * t231 + t290;
t211 = 0.2e1 * t278 + t291;
t196 = -t205 + t224;
t195 = t269 - t224;
t193 = -t222 * t245 - t300;
t192 = t223 * t249 - t306;
t191 = t205 - t269;
t190 = -t205 - t224;
t183 = -t224 - t269;
t169 = t269 + t205;
t165 = -t186 + t219;
t164 = t185 - t219;
t152 = (-qJD(4) + t227) * t207 + t260;
t148 = t186 - t185;
t144 = -t190 * t244 + t302;
t143 = t190 * t248 + t308;
t140 = t183 * t248 - t342;
t139 = t183 * t244 + t341;
t132 = (t187 * t247 - t189 * t243) * t220;
t131 = (t187 * t243 + t189 * t247) * t220;
t122 = -t151 * t248 + t155 * t244;
t120 = t164 * t247 + t311;
t119 = -t165 * t243 - t304;
t118 = t164 * t243 - t305;
t117 = t165 * t247 - t310;
t116 = t144 * t249 - t156 * t245;
t111 = t140 * t249 - t152 * t245;
t102 = pkin(5) * t109;
t97 = t127 * t247 + t189 * t313;
t96 = t127 * t243 - t189 * t312;
t95 = -t126 * t243 - t187 * t312;
t94 = t126 * t247 - t187 * t313;
t88 = t131 * t248 + t132 * t244;
t87 = t245 * (-t131 * t244 + t132 * t248) + t249 * t203;
t85 = -pkin(5) * t335 + qJ(6) * t136;
t80 = t118 * t248 + t120 * t244;
t79 = t117 * t248 + t119 * t244;
t74 = t319 - t324;
t72 = -t105 * t247 - t243 * t335;
t70 = -t105 * t243 + t247 * t335;
t65 = t322 - t325;
t64 = t244 * t97 + t248 * t96;
t63 = t244 * t95 + t248 * t94;
t56 = t245 * (-t244 * t96 + t248 * t97) - t285;
t55 = t245 * (-t244 * t94 + t248 * t95) + t285;
t54 = t245 * (-t118 * t244 + t120 * t248) + t249 * t106;
t53 = t245 * (-t117 * t244 + t119 * t248) - t249 * t109;
t52 = -qJ(6) * t160 + t58;
t50 = t245 * t78 - t249 * t335;
t48 = t277 + t322;
t47 = t280 - t319;
t45 = -t105 * t249 + t245 * t62;
t40 = -pkin(5) * t105 + qJ(6) * t142 - t58;
t37 = t244 * t72 + t248 * t70;
t34 = t245 * (-t244 * t70 + t248 * t72) - t249 * t148;
t30 = -t130 * t249 + t245 * t39;
t26 = -t243 * t85 + t247 * t52 - t324;
t24 = qJ(6) * t304 - t243 * t40 - t325;
t23 = t254 + t339;
t22 = t243 * t52 + t247 * t85 + t277;
t21 = -pkin(5) * t130 - qJ(6) * t106 + t33;
t18 = qJ(6) * t310 + t247 * t40 + t280;
t17 = -pkin(5) * t58 + qJ(6) * t33;
t16 = -pkin(4) * t99 + pkin(9) * t20;
t15 = -t19 - t326;
t13 = t20 + t279;
t12 = t247 * t33 + t323;
t10 = t20 * t248 - t321;
t9 = t20 * t244 + t318;
t8 = t10 * t249 + t245 * t99;
t7 = -t21 * t243 + t23 * t247 - t326;
t6 = t21 * t247 + t23 * t243 + t279;
t5 = -t11 * t244 + t12 * t248;
t4 = t11 * t248 + t12 * t244;
t3 = -pkin(9) * t11 + qJ(6) * t320 - t17 * t243;
t2 = t245 * t58 + t249 * t5;
t1 = -pkin(4) * t58 + pkin(9) * t12 + qJ(6) * t323 + t17 * t247;
t14 = [0, 0, 0, 0, 0, qJDD(1), t276, t267, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t240 - t239 * t251) + t273, pkin(1) * (-qJDD(1) * t239 - t240 * t251) - t298, 0, pkin(1) * (t239 * t298 + t240 * t273), t265 * t245, t211 * t249 + t214 * t245, t306 + t249 * (-t232 + t333), -t264 * t249, t245 * (t233 - t333) + t300, 0, -t249 * t174 + pkin(2) * t214 + pkin(7) * t192 + pkin(1) * (t192 * t239 + t214 * t240), t245 * t174 - pkin(2) * t211 + pkin(7) * t193 + pkin(1) * (t193 * t239 - t211 * t240), pkin(2) * t216 + pkin(7) * t215 + pkin(1) * (t215 * t239 + t216 * t240) + t128, -pkin(2) * t174 + pkin(7) * t128 + pkin(1) * (t128 * t239 - t174 * t240), t245 * (t181 * t248 + t227 * t307) + t259, t245 * (t152 * t248 + t156 * t244) - t249 * t191, t245 * (-t196 * t244 + t341) - t249 * t155, t245 * (t197 * t248 - t244 * t255) - t259, t245 * (t195 * t248 + t308) + t249 * t151, t249 * t206 + t245 * (-t248 * t271 - t307) * t227, t245 * (-pkin(8) * t139 + t309) + t249 * (-pkin(3) * t139 + t100) - pkin(2) * t139 + pkin(7) * t111 + pkin(1) * (t111 * t239 - t139 * t240), t245 * (-pkin(8) * t143 + t303) + t249 * (-pkin(3) * t143 + t101) - pkin(2) * t143 + pkin(7) * t116 + pkin(1) * (t116 * t239 - t143 * t240), t245 * t263 + t283 * (t122 * t249 - t169 * t245) + t256 * (-t151 * t244 - t155 * t248), t283 * (t146 * t245 + t249 * t67) - t256 * t263, t56, t34, t53, t55, t54, t87, t245 * (-t244 * t47 + t248 * t65 - t328) + t249 * (t270 - t332) + t288, t245 * (-t244 * t48 + t248 * t74 - t327) + t249 * (-t286 - t331) + t287, t245 * (-t13 * t244 + t15 * t248 - t329) + t249 * t282 + t289, t245 * (-pkin(8) * t9 - pkin(9) * t318 - t16 * t244) + t249 * (-pkin(3) * t9 - t330) - pkin(2) * t9 + pkin(7) * t8 + pkin(1) * (t239 * t8 - t240 * t9), t56, t34, t53, t55, t54, t87, t245 * (-t18 * t244 + t24 * t248 - t328) + t249 * (-t332 - t343) + t288, t245 * (-t22 * t244 + t248 * t26 - t327) + t249 * (-t331 - t340) + t287, t245 * (-t244 * t6 + t248 * t7 - t329) + t249 * (t102 + t282) + t289, t245 * (-pkin(8) * t4 - t1 * t244 + t248 * t3) + t249 * (-pkin(3) * t4 - t281) - pkin(2) * t4 + pkin(7) * t2 + pkin(1) * (t2 * t239 - t240 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, 0, 0, 0, 0, 0, 0, t218 * t249 + t223 * t245, -t217 * t245 + t222 * t249, 0, -t161 * t249 + t162 * t245, 0, 0, 0, 0, 0, 0, t140 * t245 + t152 * t249, t144 * t245 + t156 * t249, t122 * t245 + t169 * t249, -t146 * t249 + t245 * t67, 0, 0, 0, 0, 0, 0, t45, t50, t30, t10 * t245 - t249 * t99, 0, 0, 0, 0, 0, 0, t45, t50, t30, t245 * t5 - t249 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t226, t232 - t233, t291, t226, t290, qJDD(3), -t161, -t162, 0, 0, t181 * t244 - t227 * t301, t152 * t244 - t156 * t248, t196 * t248 + t342, t197 * t244 + t248 * t255, t195 * t244 - t302, (-t244 * t271 + t301) * t227, pkin(3) * t152 + pkin(8) * t140 - t303, pkin(3) * t156 + pkin(8) * t144 + t309, pkin(3) * t169 + pkin(8) * t122 + t67, -pkin(3) * t146 + pkin(8) * t67, t64, t37, t79, t63, t80, t88, t244 * t65 + t248 * t47 + t316, t244 * t74 + t248 * t48 + t315, t13 * t248 + t15 * t244 + t317, -pkin(3) * t99 + pkin(8) * t10 - pkin(9) * t321 + t16 * t248, t64, t37, t79, t63, t80, t88, t18 * t248 + t24 * t244 + t316, t22 * t248 + t244 * t26 + t315, t244 * t7 + t248 * t6 + t317, -pkin(3) * t58 + pkin(8) * t5 + t1 * t248 + t244 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261, t191, t155, t261, -t151, -t206, -t100, -t101, 0, 0, t150, t148, t109, -t150, -t106, -t203, -t270, t286, t69, t330, t150, t148, t109, -t150, -t106, -t203, t343, t340, -t102 + t69, t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t148, t109, -t150, -t106, -t203, -t42, -t43, 0, 0, t150, t148, t109, -t150, -t106, -t203, t253, t258, -t102, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t335, t130, t58;];
tauJ_reg  = t14;

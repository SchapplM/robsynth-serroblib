% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:55:30
% EndTime: 2019-05-05 13:55:44
% DurationCPUTime: 6.12s
% Computational Cost: add. (28862->455), mult. (67544->682), div. (0->0), fcn. (49361->12), ass. (0->263)
t241 = sin(pkin(10));
t234 = t241 ^ 2;
t244 = cos(pkin(10));
t235 = t244 ^ 2;
t310 = qJD(1) ^ 2;
t222 = (t234 + t235) * t310;
t240 = sin(pkin(11));
t247 = sin(qJ(4));
t249 = cos(qJ(4));
t259 = t241 * t249 + t244 * t247;
t216 = t259 * qJD(1);
t243 = cos(pkin(11));
t200 = -t243 * qJD(4) + t216 * t240;
t202 = qJD(4) * t240 + t216 * t243;
t170 = t202 * t200;
t274 = t244 * qJDD(1);
t275 = t241 * qJDD(1);
t262 = t247 * t275 - t249 * t274;
t278 = t216 * qJD(4);
t190 = t262 + t278;
t317 = -t170 + t190;
t328 = t240 * t317;
t327 = t243 * t317;
t246 = sin(qJ(6));
t248 = cos(qJ(6));
t163 = t248 * t200 + t202 * t246;
t165 = -t200 * t246 + t202 * t248;
t130 = t165 * t163;
t184 = qJDD(6) + t190;
t319 = -t130 + t184;
t326 = t246 * t319;
t280 = qJD(1) * t244;
t286 = t241 * t247;
t214 = qJD(1) * t286 - t249 * t280;
t194 = t216 * t214;
t314 = qJDD(4) - t194;
t325 = t247 * t314;
t324 = t248 * t319;
t323 = t249 * t314;
t242 = sin(pkin(9));
t306 = sin(qJ(1));
t307 = cos(qJ(1));
t257 = t306 * g(1) - t307 * g(2);
t256 = qJDD(1) * pkin(1) + t257;
t258 = t307 * g(1) + t306 * g(2);
t220 = -t310 * pkin(1) - t258;
t245 = cos(pkin(9));
t285 = t245 * t220;
t254 = qJDD(1) * qJ(3) + t242 * t256 + t285;
t281 = -g(3) + qJDD(2);
t308 = 2 * qJD(3);
t167 = t244 * (-t310 * pkin(2) + t254) + t241 * t281 + t280 * t308;
t277 = t235 * t310;
t155 = -pkin(3) * t277 + pkin(7) * t274 + t167;
t255 = -t242 * t257 - t285;
t266 = t244 * t281;
t267 = pkin(1) * t242 + qJ(3);
t321 = t267 + pkin(7);
t252 = t266 + (-t321 * qJDD(1) + (-(2 * qJD(3)) + (t244 * pkin(3) + pkin(2)) * qJD(1)) * qJD(1) + t255) * t241;
t123 = t249 * t155 + t247 * t252;
t217 = t245 * t256;
t263 = -t242 * t220 + t217;
t181 = -qJDD(1) * pkin(2) - t310 * qJ(3) + qJDD(3) - t263;
t269 = pkin(1) * t245 + pkin(2);
t322 = -qJDD(1) * t269 + t222 * t267 + t181;
t213 = t259 * qJDD(1);
t279 = t214 * qJD(4);
t192 = t213 - t279;
t175 = qJDD(4) * t240 + t192 * t243;
t264 = -t243 * qJDD(4) + t192 * t240;
t115 = -t163 * qJD(6) + t248 * t175 - t246 * t264;
t209 = qJD(6) + t214;
t151 = t209 * t163;
t318 = -t151 + t115;
t179 = t214 * t200;
t140 = -t175 - t179;
t316 = t175 - t179;
t265 = t246 * t175 + t248 * t264;
t93 = (qJD(6) - t209) * t165 + t265;
t160 = t163 ^ 2;
t161 = t165 ^ 2;
t312 = t200 ^ 2;
t199 = t202 ^ 2;
t208 = t209 ^ 2;
t311 = t214 ^ 2;
t212 = t216 ^ 2;
t309 = qJD(4) ^ 2;
t183 = pkin(4) * t214 - qJ(5) * t216;
t100 = -t309 * pkin(4) + qJDD(4) * qJ(5) - t214 * t183 + t123;
t162 = -pkin(3) * t274 + t181 + (-t234 * t310 - t277) * pkin(7);
t107 = (-t192 + t279) * qJ(5) + (t190 + t278) * pkin(4) + t162;
t66 = 0.2e1 * qJD(5) * t202 + t240 * t100 - t243 * t107;
t48 = pkin(5) * t317 + pkin(8) * t140 - t66;
t174 = pkin(5) * t214 - pkin(8) * t202;
t67 = -0.2e1 * qJD(5) * t200 + t243 * t100 + t240 * t107;
t51 = -t312 * pkin(5) - t264 * pkin(8) - t214 * t174 + t67;
t24 = t246 * t51 - t248 * t48;
t25 = t246 * t48 + t248 * t51;
t12 = -t24 * t248 + t246 * t25;
t305 = t12 * t240;
t304 = t12 * t243;
t122 = t247 * t155 - t249 * t252;
t99 = -qJDD(4) * pkin(4) - t309 * qJ(5) + t216 * t183 + qJDD(5) + t122;
t303 = t240 * t99;
t72 = -t122 * t249 + t247 * t123;
t302 = t241 * t72;
t301 = t243 * t99;
t70 = t264 * pkin(5) - t312 * pkin(8) + t202 * t174 + t99;
t300 = t246 * t70;
t299 = t248 * t70;
t119 = t130 + t184;
t298 = t119 * t246;
t297 = t119 * t248;
t142 = t170 + t190;
t296 = t142 * t240;
t295 = t142 * t243;
t294 = t162 * t247;
t293 = t162 * t249;
t187 = qJDD(4) + t194;
t292 = t187 * t249;
t291 = t202 * t214;
t290 = t209 * t246;
t289 = t209 * t248;
t288 = t214 * t240;
t287 = t214 * t243;
t284 = t247 * t187;
t283 = t247 * t190;
t273 = t249 * t130;
t272 = t249 * t170;
t271 = t247 * t130;
t270 = t247 * t170;
t268 = -pkin(4) * t249 - pkin(3);
t13 = t24 * t246 + t248 * t25;
t39 = t240 * t66 + t243 * t67;
t73 = t122 * t247 + t249 * t123;
t166 = -t266 + ((-pkin(2) * qJD(1) + t308) * qJD(1) + t254) * t241;
t128 = t241 * t166 + t244 * t167;
t38 = t240 * t67 - t243 * t66;
t136 = t264 - t291;
t232 = t235 * qJDD(1);
t230 = t234 * qJDD(1);
t221 = t232 + t230;
t205 = -t212 - t309;
t204 = -t212 + t309;
t203 = t311 - t309;
t191 = t213 - 0.2e1 * t279;
t189 = t262 + 0.2e1 * t278;
t185 = -t311 - t309;
t182 = t249 * t190;
t177 = -t199 + t311;
t176 = -t311 + t312;
t171 = -t311 - t212;
t168 = -t199 + t312;
t159 = -t199 - t311;
t157 = -t247 * t205 - t292;
t156 = t205 * t249 - t284;
t154 = -t311 - t312;
t150 = t247 * t213 - t249 * t262;
t149 = -t213 * t249 - t247 * t262;
t148 = -t161 + t208;
t147 = t160 - t208;
t146 = t185 * t249 - t325;
t145 = t247 * t185 + t323;
t144 = -t199 - t312;
t135 = t264 + t291;
t134 = (-t200 * t243 + t202 * t240) * t214;
t133 = t175 * t243 - t202 * t288;
t132 = t200 * t287 + t240 * t264;
t131 = -t161 - t208;
t129 = t161 - t160;
t127 = -t156 * t241 + t157 * t244;
t126 = -t208 - t160;
t125 = t176 * t243 - t296;
t124 = -t177 * t240 + t327;
t117 = -t159 * t240 - t295;
t116 = t159 * t243 - t296;
t114 = -qJD(6) * t165 - t265;
t113 = (-t163 * t248 + t165 * t246) * t209;
t112 = (-t163 * t246 - t165 * t248) * t209;
t111 = -t149 * t241 + t150 * t244;
t110 = t154 * t243 - t328;
t109 = t154 * t240 + t327;
t108 = -t145 * t241 + t146 * t244;
t104 = -t160 - t161;
t103 = -t136 * t243 - t140 * t240;
t102 = -t135 * t243 - t240 * t316;
t101 = -t136 * t240 + t140 * t243;
t96 = t151 + t115;
t92 = (qJD(6) + t209) * t165 + t265;
t91 = t147 * t248 - t298;
t90 = -t148 * t246 + t324;
t89 = t147 * t246 + t297;
t88 = t148 * t248 + t326;
t87 = t115 * t248 - t165 * t290;
t86 = t115 * t246 + t165 * t289;
t85 = -t114 * t246 + t163 * t289;
t84 = t114 * t248 + t163 * t290;
t83 = t117 * t249 + t247 * t316;
t82 = t247 * t117 - t249 * t316;
t81 = -t131 * t246 - t297;
t80 = t131 * t248 - t298;
t79 = t110 * t249 + t247 * t135;
t78 = t247 * t110 - t135 * t249;
t77 = t103 * t249 + t247 * t144;
t76 = t247 * t103 - t144 * t249;
t75 = t126 * t248 - t326;
t74 = t126 * t246 + t324;
t71 = -t112 * t240 + t113 * t243;
t69 = -qJ(5) * t116 + t301;
t68 = -qJ(5) * t109 + t303;
t64 = t246 * t96 - t248 * t93;
t63 = -t246 * t318 - t248 * t92;
t62 = -t246 * t93 - t248 * t96;
t61 = -t246 * t92 + t248 * t318;
t60 = -t240 * t89 + t243 * t91;
t59 = -t240 * t88 + t243 * t90;
t58 = -t240 * t86 + t243 * t87;
t57 = -t240 * t84 + t243 * t85;
t56 = -t241 * t82 + t244 * t83;
t55 = -t240 * t80 + t243 * t81;
t54 = t240 * t81 + t243 * t80;
t53 = -t241 * t78 + t244 * t79;
t50 = -pkin(4) * t116 + t67;
t49 = -pkin(4) * t109 + t66;
t46 = -t240 * t74 + t243 * t75;
t45 = t240 * t75 + t243 * t74;
t44 = t244 * t73 - t302;
t43 = -pkin(8) * t80 + t299;
t42 = -pkin(8) * t74 + t300;
t41 = t247 * t318 + t249 * t55;
t40 = t247 * t55 - t249 * t318;
t37 = t247 * t92 + t249 * t46;
t36 = t247 * t46 - t249 * t92;
t35 = -pkin(5) * t318 + pkin(8) * t81 + t300;
t34 = -pkin(5) * t92 + pkin(8) * t75 - t299;
t33 = -t240 * t62 + t243 * t64;
t32 = -t240 * t61 + t243 * t63;
t31 = t240 * t64 + t243 * t62;
t30 = -qJ(5) * t101 - t38;
t29 = t247 * t99 + t249 * t39;
t28 = t247 * t39 - t249 * t99;
t27 = t247 * t104 + t249 * t33;
t26 = -t104 * t249 + t247 * t33;
t22 = -pkin(4) * t31 - pkin(5) * t62;
t21 = -t241 * t40 + t244 * t41;
t20 = -t241 * t36 + t244 * t37;
t19 = -pkin(4) * t54 - pkin(5) * t80 + t25;
t18 = -qJ(5) * t54 - t240 * t35 + t243 * t43;
t17 = -pkin(4) * t45 - pkin(5) * t74 + t24;
t15 = -qJ(5) * t45 - t240 * t34 + t243 * t42;
t14 = -t241 * t26 + t244 * t27;
t11 = -pkin(5) * t70 + pkin(8) * t13;
t10 = -pkin(8) * t62 - t12;
t9 = -pkin(5) * t104 + pkin(8) * t64 + t13;
t8 = t13 * t243 - t305;
t7 = t13 * t240 + t304;
t6 = t247 * t70 + t249 * t8;
t5 = t247 * t8 - t249 * t70;
t4 = -qJ(5) * t31 + t10 * t243 - t240 * t9;
t3 = -pkin(4) * t7 - pkin(5) * t12;
t2 = -pkin(8) * t304 - qJ(5) * t7 - t11 * t240;
t1 = -t241 * t5 + t244 * t6;
t16 = [0, 0, 0, 0, 0, qJDD(1), t257, t258, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t245 - t310 * t242) + t263, (-0.2e1 * qJDD(1) * t242 - t310 * t245) * pkin(1) + t255, 0, pkin(1) * (t242 ^ 2 * t256 + t245 * t217), t230, 0.2e1 * t241 * t274, 0, t232, 0, 0, -t322 * t244, t322 * t241, pkin(2) * t222 + qJ(3) * t221 + pkin(1) * (t221 * t242 + t222 * t245) + t128, -pkin(2) * t181 + qJ(3) * t128 + pkin(1) * (t128 * t242 - t181 * t245), t241 * (t192 * t249 - t247 * t278) + t244 * (t247 * t192 + t249 * t278), t241 * (-t189 * t249 - t247 * t191) + t244 * (-t247 * t189 + t191 * t249), t241 * (-t247 * t204 + t323) + t244 * (t204 * t249 + t325), t241 * (t249 * t279 + t283) + t244 * (t247 * t279 - t182), t241 * (t203 * t249 - t284) + t244 * (t247 * t203 + t292), (t241 * (-t214 * t249 + t216 * t247) + t244 * (-t214 * t247 - t216 * t249)) * qJD(4), t241 * (-pkin(7) * t145 + t294) + t244 * (-pkin(3) * t189 + pkin(7) * t146 - t293) - pkin(2) * t189 + qJ(3) * t108 + pkin(1) * (t108 * t242 - t189 * t245), t241 * (-pkin(7) * t156 + t293) + t244 * (-pkin(3) * t191 + pkin(7) * t157 + t294) - pkin(2) * t191 + qJ(3) * t127 + pkin(1) * (t127 * t242 - t191 * t245), t241 * (-pkin(7) * t149 - t72) + t244 * (-pkin(3) * t171 + pkin(7) * t150 + t73) - pkin(2) * t171 + qJ(3) * t111 + pkin(1) * (t111 * t242 - t171 * t245), -pkin(7) * t302 + t244 * (-pkin(3) * t162 + pkin(7) * t73) - pkin(2) * t162 + qJ(3) * t44 + pkin(1) * (-t162 * t245 + t242 * t44), t241 * (t133 * t249 + t270) + t244 * (t247 * t133 - t272), t241 * (t102 * t249 - t247 * t168) + t244 * (t247 * t102 + t168 * t249), t241 * (t124 * t249 - t247 * t140) + t244 * (t247 * t124 + t140 * t249), t241 * (t132 * t249 - t270) + t244 * (t247 * t132 + t272), t241 * (t125 * t249 - t247 * t136) + t244 * (t247 * t125 + t136 * t249), t241 * (t134 * t249 + t283) + t244 * (t134 * t247 - t182), t241 * (-pkin(7) * t78 - t247 * t49 + t249 * t68) + t244 * (-pkin(3) * t109 + pkin(7) * t79 + t247 * t68 + t249 * t49) - pkin(2) * t109 + qJ(3) * t53 + pkin(1) * (-t109 * t245 + t242 * t53), t241 * (-pkin(7) * t82 - t247 * t50 + t249 * t69) + t244 * (-pkin(3) * t116 + pkin(7) * t83 + t247 * t69 + t249 * t50) - pkin(2) * t116 + qJ(3) * t56 + pkin(1) * (-t116 * t245 + t242 * t56), t241 * (-pkin(7) * t76 + t249 * t30) + t244 * (pkin(7) * t77 + t247 * t30) + t267 * (-t241 * t76 + t244 * t77) + (pkin(4) * t286 + t244 * t268 - t269) * t101, (t241 * (pkin(4) * t247 - qJ(5) * t249) + t244 * (-qJ(5) * t247 + t268) - t269) * t38 + t321 * (-t241 * t28 + t244 * t29), t241 * (t249 * t58 + t271) + t244 * (t247 * t58 - t273), t241 * (t247 * t129 + t249 * t32) + t244 * (-t129 * t249 + t247 * t32), t241 * (t247 * t96 + t249 * t59) + t244 * (t247 * t59 - t249 * t96), t241 * (t249 * t57 - t271) + t244 * (t247 * t57 + t273), t241 * (-t247 * t93 + t249 * t60) + t244 * (t247 * t60 + t249 * t93), t241 * (t247 * t184 + t249 * t71) + t244 * (-t184 * t249 + t247 * t71), t241 * (-pkin(7) * t36 + t15 * t249 - t247 * t17) + t244 * (-pkin(3) * t45 + pkin(7) * t37 + t247 * t15 + t17 * t249) - pkin(2) * t45 + qJ(3) * t20 + pkin(1) * (t20 * t242 - t245 * t45), t241 * (-pkin(7) * t40 + t18 * t249 - t247 * t19) + t244 * (-pkin(3) * t54 + pkin(7) * t41 + t247 * t18 + t19 * t249) - pkin(2) * t54 + qJ(3) * t21 + pkin(1) * (t21 * t242 - t245 * t54), t241 * (-pkin(7) * t26 - t247 * t22 + t249 * t4) + t244 * (-pkin(3) * t31 + pkin(7) * t27 + t22 * t249 + t247 * t4) - pkin(2) * t31 + qJ(3) * t14 + pkin(1) * (t14 * t242 - t245 * t31), t241 * (-pkin(7) * t5 + t2 * t249 - t247 * t3) + t244 * (-pkin(3) * t7 + pkin(7) * t6 + t247 * t2 + t249 * t3) - pkin(2) * t7 + qJ(3) * t1 + pkin(1) * (t1 * t242 - t245 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t281, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166 * t244 + t167 * t241, 0, 0, 0, 0, 0, 0, t145 * t244 + t146 * t241, t156 * t244 + t157 * t241, t149 * t244 + t150 * t241, t241 * t73 + t244 * t72, 0, 0, 0, 0, 0, 0, t241 * t79 + t244 * t78, t241 * t83 + t244 * t82, t241 * t77 + t244 * t76, t241 * t29 + t244 * t28, 0, 0, 0, 0, 0, 0, t241 * t37 + t244 * t36, t241 * t41 + t244 * t40, t241 * t27 + t244 * t26, t241 * t6 + t244 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t274, t275, -t222, t181, 0, 0, 0, 0, 0, 0, t189, t191, t171, t162, 0, 0, 0, 0, 0, 0, t109, t116, t101, t38, 0, 0, 0, 0, 0, 0, t45, t54, t31, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, t212 - t311, t213, -t194, -t262, qJDD(4), -t122, -t123, 0, 0, t175 * t240 + t202 * t287, -t135 * t240 + t243 * t316, t177 * t243 + t328, t200 * t288 - t243 * t264, t176 * t240 + t295, (-t200 * t240 - t202 * t243) * t214, -pkin(4) * t135 + qJ(5) * t110 - t301, -pkin(4) * t316 + qJ(5) * t117 + t303, -pkin(4) * t144 + qJ(5) * t103 + t39, -pkin(4) * t99 + qJ(5) * t39, t240 * t87 + t243 * t86, t240 * t63 + t243 * t61, t240 * t90 + t243 * t88, t240 * t85 + t243 * t84, t240 * t91 + t243 * t89, t112 * t243 + t113 * t240, -pkin(4) * t92 + qJ(5) * t46 + t240 * t42 + t243 * t34, -pkin(4) * t318 + qJ(5) * t55 + t240 * t43 + t243 * t35, -pkin(4) * t104 + qJ(5) * t33 + t10 * t240 + t243 * t9, -pkin(4) * t70 - pkin(8) * t305 + qJ(5) * t8 + t11 * t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t316, t144, t99, 0, 0, 0, 0, 0, 0, t92, t318, t104, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t129, t96, -t130, -t93, t184, -t24, -t25, 0, 0;];
tauJ_reg  = t16;

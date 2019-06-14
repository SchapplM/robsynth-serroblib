% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 06:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRPP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:46:13
% EndTime: 2019-05-05 06:46:27
% DurationCPUTime: 5.37s
% Computational Cost: add. (7480->356), mult. (14519->409), div. (0->0), fcn. (10007->10), ass. (0->223)
t194 = sin(qJ(3));
t197 = cos(qJ(3));
t193 = sin(qJ(4));
t196 = cos(qJ(4));
t244 = qJD(2) * t194;
t154 = -t196 * qJD(3) + t193 * t244;
t239 = qJD(2) * qJD(3);
t234 = t197 * t239;
t238 = t194 * qJDD(2);
t159 = t234 + t238;
t219 = -qJDD(3) * t193 - t159 * t196;
t210 = qJD(4) * t154 + t219;
t177 = qJD(2) * t197 - qJD(4);
t251 = t177 * t154;
t298 = -t251 - t210;
t156 = qJD(3) * t193 + t196 * t244;
t230 = -t196 * qJDD(3) + t159 * t193;
t76 = (qJD(4) + t177) * t156 + t230;
t43 = -t193 * t298 + t196 * t76;
t152 = t156 ^ 2;
t286 = t154 ^ 2;
t91 = -t286 - t152;
t29 = -t194 * t91 + t197 * t43;
t358 = pkin(8) * t29;
t40 = t193 * t76 + t196 * t298;
t388 = -pkin(2) * t40 + t358;
t189 = sin(pkin(6));
t190 = cos(pkin(6));
t195 = sin(qJ(2));
t198 = cos(qJ(2));
t387 = t190 * (t194 * t43 + t197 * t91) + (t195 * t29 - t198 * t40) * t189;
t297 = t251 - t210;
t285 = t177 ^ 2;
t112 = -t285 - t152;
t122 = t156 * t154;
t180 = t194 * t239;
t237 = t197 * qJDD(2);
t160 = -t180 + t237;
t153 = -qJDD(4) + t160;
t303 = -t122 + t153;
t311 = t196 * t303;
t356 = t112 * t193 - t311;
t365 = t197 * t356;
t360 = -t194 * t297 + t365;
t317 = t193 * t303;
t60 = t112 * t196 + t317;
t376 = pkin(2) * t60;
t386 = pkin(8) * t360 + t376;
t366 = t194 * t356;
t369 = t198 * t60;
t384 = t190 * (t197 * t297 + t366) + (t195 * t360 + t369) * t189;
t374 = pkin(3) * t60;
t372 = pkin(9) * t60;
t295 = t286 - t285;
t328 = t194 * (t196 * t295 + t317);
t380 = t197 * t76 + t328;
t296 = -t285 - t286;
t302 = t122 + t153;
t312 = t196 * t302;
t332 = t193 * t296 - t312;
t348 = pkin(2) * t332;
t318 = t193 * t302;
t331 = t196 * t296 + t318;
t338 = t197 * t331;
t109 = qJD(4) * t156 + t230;
t136 = t156 * t177;
t74 = t109 - t136;
t364 = t194 * t74 + t338;
t379 = pkin(8) * t364 - t348;
t337 = t198 * t332;
t339 = t194 * t331;
t378 = t190 * (-t197 * t74 + t339) + (t195 * t364 - t337) * t189;
t375 = pkin(3) * t40;
t373 = pkin(9) * t40;
t371 = pkin(9) * t356;
t346 = pkin(9) * t331;
t368 = -pkin(3) * t74 + t346;
t330 = -pkin(3) * t91 - pkin(9) * t43;
t363 = pkin(3) * t297 + t371;
t130 = t152 - t285;
t361 = t194 * (-t130 * t193 + t312) + t197 * t298;
t347 = pkin(3) * t332;
t345 = pkin(9) * t332;
t355 = -t193 * t295 + t311;
t352 = t130 * t196 + t318;
t299 = t152 - t286;
t114 = t197 * t299;
t276 = t193 * t297;
t241 = qJD(4) - t177;
t75 = t241 * t156 + t230;
t351 = t194 * (t196 * t75 + t276) + t114;
t344 = qJ(5) * t76;
t343 = qJ(5) * t91;
t340 = qJ(5) * t297;
t266 = t196 * t297;
t335 = -t193 * t75 + t266;
t258 = sin(pkin(10));
t259 = cos(pkin(10));
t215 = t258 * g(1) - t259 * g(2);
t246 = -g(3) + qJDD(1);
t203 = -t189 * t215 + t190 * t246;
t126 = t197 * t203;
t199 = qJD(2) ^ 2;
t228 = -t197 * pkin(3) - t194 * pkin(9);
t164 = -t259 * g(1) - t258 * g(2);
t209 = t190 * t215;
t307 = t189 * t246 + t209;
t100 = t198 * t164 + t307 * t195;
t88 = -t199 * pkin(2) + qJDD(2) * pkin(8) + t100;
t231 = t199 * t228 + t88;
t283 = qJD(3) ^ 2;
t48 = -qJDD(3) * pkin(3) - t283 * pkin(9) + t231 * t194 - t126;
t324 = t109 * pkin(4) - t340 + t48;
t323 = qJ(5) * t303;
t322 = qJ(6) * t298;
t304 = -pkin(4) * t302 + qJ(5) * t296;
t243 = qJD(5) * t177;
t167 = -0.2e1 * t243;
t242 = qJD(6) * t154;
t301 = 0.2e1 * t242 + t167;
t168 = 0.2e1 * t243;
t300 = -0.2e1 * t242 + t168;
t116 = t197 * t122;
t249 = t177 * t196;
t236 = t154 * t249;
t208 = t194 * (t109 * t193 - t236) + t116;
t223 = pkin(5) * t177 - qJ(6) * t156;
t294 = t156 * t223 + qJDD(6);
t292 = -pkin(5) * t109 + t294;
t200 = t194 * t203;
t49 = -t283 * pkin(3) + qJDD(3) * pkin(9) + t231 * t197 + t200;
t221 = -t160 + t180;
t222 = t159 + t234;
t224 = t164 * t195 - t307 * t198;
t87 = -qJDD(2) * pkin(2) - t199 * pkin(8) + t224;
t52 = t221 * pkin(3) - t222 * pkin(9) + t87;
t27 = t193 * t49 - t196 * t52;
t217 = t153 * pkin(4) - qJ(5) * t285 + qJDD(5) + t27;
t206 = t153 * pkin(5) + t217 - t322;
t115 = pkin(4) * t154 - qJ(5) * t156;
t233 = -pkin(5) * t154 - t115;
t218 = (-0.2e1 * qJD(6) - t233) * t156;
t16 = t218 + t206;
t28 = t193 * t52 + t196 * t49;
t229 = -pkin(4) * t285 - t153 * qJ(5) - t154 * t115 + t28;
t214 = -pkin(5) * t286 + t109 * qJ(6) - t177 * t223 + t229;
t17 = t214 + t301;
t282 = pkin(4) + pkin(5);
t291 = qJ(5) * t17 - t282 * t16;
t290 = t282 * t298 + t344;
t250 = t177 * t193;
t127 = t156 * t250;
t260 = t194 * (-t196 * t210 + t127) - t116;
t289 = -t112 * t282 + t214 - t323;
t284 = 0.2e1 * t156;
t278 = t193 * t48;
t268 = t196 * t48;
t257 = qJ(5) * t193;
t256 = qJ(5) * t196;
t255 = t115 * t156;
t176 = t194 * t199 * t197;
t166 = qJDD(3) + t176;
t248 = t194 * t166;
t165 = -t176 + qJDD(3);
t247 = t197 * t165;
t15 = t193 * t27 + t196 * t28;
t66 = t194 * t88 - t126;
t67 = t197 * t88 + t200;
t33 = t194 * t66 + t197 * t67;
t23 = t167 + t229;
t24 = t217 + t255;
t227 = -pkin(4) * t24 + qJ(5) * t23;
t226 = -pkin(4) * t298 - t344;
t69 = -t156 * t249 - t193 * t210;
t225 = -t196 * t109 - t154 * t250;
t14 = t193 * t28 - t196 * t27;
t220 = -pkin(2) + t228;
t216 = (t154 * t193 + t156 * t196) * t177;
t139 = t197 * t153;
t212 = t194 * (-t127 + t236) + t139;
t207 = -pkin(4) * t112 + t229 - t323;
t204 = qJD(5) * t284 - t324;
t202 = -t24 + t304;
t201 = -t206 + t304;
t25 = (-pkin(4) * t177 - 0.2e1 * qJD(5)) * t156 + t324;
t22 = (-t74 + t136) * pkin(4) + t204;
t21 = pkin(4) * t136 + t204 + t340;
t186 = t197 ^ 2;
t185 = t194 ^ 2;
t184 = t186 * t199;
t182 = t185 * t199;
t173 = -t184 - t283;
t172 = -t182 - t283;
t163 = t182 + t184;
t162 = (t185 + t186) * qJDD(2);
t161 = -0.2e1 * t180 + t237;
t158 = 0.2e1 * t234 + t238;
t143 = qJD(6) * t284;
t125 = -t172 * t194 - t247;
t124 = t173 * t197 - t248;
t85 = t241 * t154 + t219;
t77 = -t109 - t136;
t46 = -qJ(5) * t74 - qJ(6) * t302;
t37 = -t194 * t85 - t365;
t36 = t194 * t75 + t338;
t32 = qJ(6) * t303 + t282 * t297;
t20 = t24 - t343;
t19 = -pkin(4) * t91 + t23;
t18 = qJ(6) * t286 + t25 - t292;
t13 = t21 + (-t112 - t286) * qJ(6) + t292;
t12 = (-t296 - t286) * qJ(6) + (-t109 - t74) * pkin(5) + t22 + t294;
t11 = t233 * t156 + t143 - t206 + t322 + t343;
t10 = t15 * t197 + t194 * t48;
t9 = -qJ(6) * t76 + t282 * t91 - t214 + t300;
t8 = t193 * t24 + t196 * t23;
t7 = t193 * t23 - t196 * t24;
t6 = -qJ(5) * t18 - qJ(6) * t16;
t5 = t16 * t193 + t17 * t196;
t4 = -t16 * t196 + t17 * t193;
t3 = t194 * t25 + t197 * t8;
t2 = -qJ(6) * t17 - t282 * t18;
t1 = t18 * t194 + t197 * t5;
t26 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t246, 0, 0, 0, 0, 0, 0, (qJDD(2) * t198 - t195 * t199) * t189, (-qJDD(2) * t195 - t198 * t199) * t189, 0, t190 ^ 2 * t246 + (t195 * t100 - t198 * t224 - t209) * t189, 0, 0, 0, 0, 0, 0, t190 * (t166 * t197 + t173 * t194) + (t124 * t195 + t161 * t198) * t189, t190 * (-t165 * t194 + t172 * t197) + (t125 * t195 - t158 * t198) * t189, (t162 * t195 + t163 * t198) * t189, t190 * (t194 * t67 - t197 * t66) + (t195 * t33 - t198 * t87) * t189, 0, 0, 0, 0, 0, 0, t190 * (-t197 * t75 + t339) + (t195 * t36 - t337) * t189, t190 * (t197 * t85 - t366) + (t195 * t37 - t369) * t189, -t387, t190 * (t15 * t194 - t197 * t48) + (t10 * t195 - t14 * t198) * t189, 0, 0, 0, 0, 0, 0, t378, -t387, t384, t190 * (t194 * t8 - t197 * t25) + (t195 * t3 - t198 * t7) * t189, 0, 0, 0, 0, 0, 0, t378, t384, t387, t190 * (-t18 * t197 + t194 * t5) + (t1 * t195 - t198 * t4) * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t224, -t100, 0, 0, t222 * t194, t158 * t197 + t161 * t194, t248 + t197 * (-t182 + t283), -t221 * t197, t194 * (t184 - t283) + t247, 0, pkin(2) * t161 + pkin(8) * t124 - t197 * t87, -pkin(2) * t158 + pkin(8) * t125 + t194 * t87, pkin(2) * t163 + pkin(8) * t162 + t33, -pkin(2) * t87 + pkin(8) * t33, t260, -t351, -t361, t208, t380, t212, t194 * (t278 - t345) + t197 * (t27 - t347) - t348 + pkin(8) * t36, t194 * (t268 - t372) + t197 * (t28 - t374) - t376 + pkin(8) * t37, -t14 * t194 - t220 * t40 - t358, pkin(8) * t10 + t14 * t220, t260, -t361, t351, t212, -t380, t208, t194 * (-t193 * t22 - t256 * t74 - t345) + t197 * (-t202 - t347) + t379, t194 * (-t19 * t193 + t196 * t20 + t373) + t197 * (-t226 + t375) - t388, t194 * (-pkin(4) * t276 + t196 * t21 + t372) + t197 * (t168 - t207 + t374) + t386, t194 * (-pkin(9) * t7 + (pkin(4) * t193 - t256) * t25) + t197 * (-pkin(3) * t7 - t227) - pkin(2) * t7 + pkin(8) * t3, t260, t194 * (t196 * t74 + t276) + t114, t361, t208, -t197 * t77 + t328, t139 + t194 * (t154 * t196 - t156 * t193) * t177, t194 * (-t193 * t12 + t196 * t46 - t345) + (pkin(5) * t302 - t201 + t218 - t347) * t197 + t379, t194 * (t13 * t196 - t193 * t32 + t372) + t197 * (-t289 + t300 + t374) + t386, t194 * (t11 * t196 - t193 * t9 - t373) + t197 * (-t290 - t375) + t388, t194 * (-pkin(9) * t4 - t193 * t2 + t196 * t6) + t197 * (-pkin(3) * t4 - t291) - pkin(2) * t4 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t182 - t184, t238, t176, t237, qJDD(3), -t66, -t67, 0, 0, t69, t335, -t352, t225, -t355, t216, -pkin(3) * t75 - t268 + t346, pkin(3) * t85 + t278 - t371, t15 + t330, -pkin(3) * t48 + pkin(9) * t15, t69, -t352, -t335, t216, t355, t225, t196 * t22 - t257 * t74 + t368, t19 * t196 + t193 * t20 + t330, pkin(4) * t266 + t193 * t21 + t363, pkin(9) * t8 + (-pkin(4) * t196 - pkin(3) - t257) * t25, t69, t193 * t74 - t266, t352, t225, -t355, t216, t12 * t196 + t193 * t46 + t368, t13 * t193 + t196 * t32 + t363, t11 * t193 + t196 * t9 - t330, -pkin(3) * t18 + pkin(9) * t5 + t193 * t6 + t196 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t299, t298, -t122, -t76, -t153, -t27, -t28, 0, 0, t122, t298, -t299, -t153, t76, -t122, t202, t226, t167 + t207, t227, t122, -t299, -t298, -t122, t77, -t153, -t255 + t143 + (-t302 - t122) * pkin(5) + t201, t289 + t301, t290, t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t298, t112, t24, 0, 0, 0, 0, 0, 0, t302, t112, -t298, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t297, t91, -t18;];
tauJ_reg  = t26;

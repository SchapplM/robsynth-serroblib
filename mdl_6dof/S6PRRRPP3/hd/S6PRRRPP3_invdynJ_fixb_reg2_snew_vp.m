% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRRPP3
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
% Datum: 2019-05-05 06:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRRPP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:53:56
% EndTime: 2019-05-05 06:54:08
% DurationCPUTime: 4.83s
% Computational Cost: add. (7484->362), mult. (14478->421), div. (0->0), fcn. (9963->10), ass. (0->226)
t189 = sin(pkin(6));
t190 = cos(pkin(6));
t195 = sin(qJ(2));
t197 = cos(qJ(3));
t193 = sin(qJ(4));
t196 = cos(qJ(4));
t194 = sin(qJ(3));
t244 = qJD(2) * t194;
t155 = qJD(3) * t193 + t196 * t244;
t239 = qJD(2) * qJD(3);
t232 = t197 * t239;
t238 = t194 * qJDD(2);
t159 = t232 + t238;
t227 = -qJDD(3) * t196 + t159 * t193;
t109 = qJD(4) * t155 + t227;
t177 = qJD(2) * t197 - qJD(4);
t138 = t155 * t177;
t295 = t109 - t138;
t198 = cos(qJ(2));
t180 = t194 * t239;
t237 = t197 * qJDD(2);
t160 = -t180 + t237;
t152 = -qJDD(4) + t160;
t153 = -qJD(3) * t196 + t193 * t244;
t251 = t155 * t153;
t212 = t152 + t251;
t263 = t196 * t212;
t175 = t177 ^ 2;
t285 = t153 ^ 2;
t293 = -t175 - t285;
t311 = t193 * t293 - t263;
t329 = t198 * t311;
t269 = t193 * t212;
t312 = t196 * t293 + t269;
t332 = t194 * t312;
t330 = t197 * t312;
t39 = t194 * t295 + t330;
t363 = t189 * (t195 * t39 - t329) + t190 * (-t197 * t295 + t332);
t284 = t155 ^ 2;
t291 = -t284 - t175;
t213 = t152 - t251;
t296 = t196 * t213;
t305 = t193 * t291 - t296;
t333 = t194 * t305;
t297 = t193 * t213;
t60 = t196 * t291 + t297;
t335 = t198 * t60;
t331 = t197 * t305;
t217 = qJDD(3) * t193 + t159 * t196;
t240 = qJD(4) - t177;
t83 = t153 * t240 - t217;
t38 = t194 * t83 + t331;
t362 = t189 * (t195 * t38 + t335) + t190 * (-t197 * t83 + t333);
t242 = qJD(5) * t177;
t342 = pkin(3) * t60;
t361 = 0.2e1 * t242 + t342;
t340 = pkin(2) * t311;
t360 = pkin(8) * t39 - t340;
t343 = pkin(2) * t60;
t359 = -pkin(8) * t38 - t343;
t210 = qJD(4) * t153 - t217;
t252 = t153 * t177;
t84 = t252 + t210;
t265 = t196 * t84;
t78 = (-qJD(4) - t177) * t155 - t227;
t327 = t193 * t78 + t265;
t92 = t284 + t285;
t318 = t194 * t92;
t272 = t193 * t84;
t334 = t196 * t78 - t272;
t345 = t197 * t334 - t318;
t358 = -pkin(2) * t327 + pkin(8) * t345;
t316 = t197 * t92;
t355 = t190 * (t194 * t334 + t316) + (t195 * t345 - t198 * t327) * t189;
t341 = pkin(9) * t60;
t353 = pkin(3) * t327;
t352 = pkin(9) * t327;
t130 = t285 - t175;
t298 = t194 * (t130 * t196 + t297);
t348 = -t197 * t78 + t298;
t322 = pkin(3) * t92;
t347 = pkin(9) * t334 + t322;
t131 = -t284 + t175;
t346 = t194 * (t131 * t193 + t263) - t197 * t84;
t339 = pkin(3) * t311;
t338 = pkin(9) * t305;
t337 = pkin(9) * t311;
t336 = pkin(9) * t312;
t313 = t131 * t196 - t269;
t292 = t252 - t210;
t273 = t193 * t292;
t117 = t285 - t284;
t314 = t197 * t117;
t75 = t155 * t240 + t227;
t323 = t194 * (t196 * t75 + t273) - t314;
t321 = (t295 - t138) * pkin(4);
t320 = qJ(5) * t92;
t315 = qJ(5) * t292;
t259 = sin(pkin(10));
t260 = cos(pkin(10));
t215 = g(1) * t259 - g(2) * t260;
t246 = -g(3) + qJDD(1);
t204 = -t189 * t215 + t190 * t246;
t127 = t197 * t204;
t199 = qJD(2) ^ 2;
t223 = -pkin(3) * t197 - pkin(9) * t194;
t164 = -g(1) * t260 - g(2) * t259;
t208 = t190 * t215;
t300 = t189 * t246 + t208;
t101 = t198 * t164 + t195 * t300;
t86 = -t199 * pkin(2) + qJDD(2) * pkin(8) + t101;
t228 = t199 * t223 + t86;
t283 = qJD(3) ^ 2;
t48 = -qJDD(3) * pkin(3) - pkin(9) * t283 + t194 * t228 - t127;
t207 = t109 * pkin(4) - t315 + t48;
t205 = 0.2e1 * qJD(5) * t155 - t207;
t266 = t196 * t292;
t310 = -t193 * t75 + t266;
t307 = t130 * t193 - t296;
t304 = qJ(5) * t293;
t294 = pkin(5) * t285 - 0.2e1 * qJD(6) * t153;
t278 = pkin(4) + qJ(6);
t73 = qJ(5) * t78;
t289 = t278 * t84 + t73;
t89 = qJ(5) * t213;
t288 = -t278 * t291 - t89;
t113 = pkin(4) * t153 - qJ(5) * t155;
t200 = t194 * t204;
t49 = -pkin(3) * t283 + qJDD(3) * pkin(9) + t197 * t228 + t200;
t219 = -t160 + t180;
t220 = t159 + t232;
t221 = t164 * t195 - t198 * t300;
t85 = -qJDD(2) * pkin(2) - pkin(8) * t199 + t221;
t52 = pkin(3) * t219 - pkin(9) * t220 + t85;
t27 = t193 * t49 - t196 * t52;
t24 = pkin(4) * t152 - qJ(5) * t175 + t113 * t155 + qJDD(5) + t27;
t206 = -t210 * pkin(5) + qJ(6) * t212 + t24;
t282 = 0.2e1 * qJD(6);
t16 = (-pkin(5) * t153 + t282) * t177 + t206;
t169 = -0.2e1 * t242;
t128 = pkin(5) * t155 + qJ(6) * t177;
t28 = t193 * t52 + t196 * t49;
t224 = -pkin(4) * t175 - qJ(5) * t152 - t113 * t153 + t28;
t209 = -pkin(5) * t109 - qJ(6) * t285 - t128 * t177 + qJDD(6) + t224;
t17 = t169 + t209;
t287 = qJ(5) * t17 - t16 * t278;
t250 = t177 * t193;
t129 = t155 * t250;
t235 = t197 * t251;
t225 = t194 * (-t196 * t210 + t129) - t235;
t286 = t212 * t278 - t304;
t281 = pkin(4) * t177;
t280 = pkin(4) * t193;
t279 = pkin(4) * t196;
t276 = t193 * t48;
t268 = t196 * t48;
t257 = qJ(5) * t196;
t256 = qJ(6) * t109;
t249 = t177 * t196;
t176 = t194 * t199 * t197;
t166 = qJDD(3) + t176;
t248 = t194 * t166;
t165 = -t176 + qJDD(3);
t247 = t197 * t165;
t236 = t153 * t249;
t233 = pkin(4) * t84 + t73;
t231 = qJ(5) * t193 + pkin(3);
t15 = t193 * t27 + t196 * t28;
t66 = t194 * t86 - t127;
t67 = t197 * t86 + t200;
t33 = t194 * t66 + t197 * t67;
t230 = -0.2e1 * qJD(5) - t281;
t226 = t194 * (t109 * t193 - t236) + t235;
t23 = t169 + t224;
t222 = -pkin(4) * t24 + qJ(5) * t23;
t71 = -t155 * t249 - t193 * t210;
t70 = -t109 * t196 - t153 * t250;
t14 = t193 * t28 - t196 * t27;
t218 = -pkin(2) + t223;
t216 = (t153 * t193 + t155 * t196) * t177;
t139 = t197 * t152;
t214 = t194 * (-t129 + t236) + t139;
t211 = -pkin(4) * t291 + t224 - t89;
t203 = pkin(4) * t212 + t24 - t304;
t202 = t205 + t294;
t201 = t177 * t282 + t206;
t186 = t197 ^ 2;
t185 = t194 ^ 2;
t184 = t186 * t199;
t182 = t185 * t199;
t174 = -t184 - t283;
t173 = -t182 - t283;
t163 = t182 + t184;
t162 = (t185 + t186) * qJDD(2);
t161 = -0.2e1 * t180 + t237;
t158 = 0.2e1 * t232 + t238;
t126 = -t173 * t194 - t247;
t125 = t174 * t197 - t248;
t76 = t109 + t138;
t47 = pkin(5) * t212 - qJ(5) * t295;
t45 = -t196 * t76 - t272;
t42 = -t193 * t76 + t265;
t35 = t194 * t75 + t330;
t34 = -t194 * t292 + t331;
t32 = -pkin(5) * t213 + t278 * t292;
t31 = t197 * t45 - t318;
t25 = t155 * t230 + t207;
t22 = -t205 + t321;
t21 = pkin(4) * t138 - qJ(5) * t83 + t205;
t20 = t24 + t320;
t19 = pkin(4) * t92 + t23;
t18 = t256 + (-t128 + t230) * t155 + t207 - t294;
t13 = (t128 + t281) * t155 + pkin(5) * t291 + t315 - t256 + t202;
t12 = t320 + (-t84 - t252) * pkin(5) + t201;
t11 = (-t109 - t295) * qJ(6) - t321 + pkin(5) * t293 + t128 * t155 + t202;
t10 = t15 * t197 + t194 * t48;
t9 = pkin(5) * t78 + t278 * t92 + t17;
t8 = t193 * t24 + t196 * t23;
t7 = t193 * t23 - t196 * t24;
t6 = pkin(5) * t16 - qJ(5) * t18;
t5 = t16 * t193 + t17 * t196;
t4 = -t16 * t196 + t17 * t193;
t3 = t194 * t25 + t197 * t8;
t2 = pkin(5) * t17 - t18 * t278;
t1 = t18 * t194 + t197 * t5;
t26 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t246, 0, 0, 0, 0, 0, 0, (qJDD(2) * t198 - t195 * t199) * t189, (-qJDD(2) * t195 - t198 * t199) * t189, 0, t190 ^ 2 * t246 + (t195 * t101 - t198 * t221 - t208) * t189, 0, 0, 0, 0, 0, 0, t190 * (t166 * t197 + t174 * t194) + (t125 * t195 + t161 * t198) * t189, t190 * (-t165 * t194 + t173 * t197) + (t126 * t195 - t158 * t198) * t189, (t162 * t195 + t163 * t198) * t189, t190 * (t194 * t67 - t197 * t66) + (t195 * t33 - t198 * t85) * t189, 0, 0, 0, 0, 0, 0, t190 * (-t197 * t75 + t332) + (t195 * t35 - t329) * t189, -t362, t190 * (t194 * t45 + t316) + (t195 * t31 - t198 * t42) * t189, t190 * (t15 * t194 - t197 * t48) + (t10 * t195 - t14 * t198) * t189, 0, 0, 0, 0, 0, 0, t355, -t363, t362, t190 * (t194 * t8 - t197 * t25) + (t195 * t3 - t198 * t7) * t189, 0, 0, 0, 0, 0, 0, t355, t190 * (t197 * t292 + t333) + (t195 * t34 + t335) * t189, t363, t190 * (-t18 * t197 + t194 * t5) + (t1 * t195 - t198 * t4) * t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t221, -t101, 0, 0, t220 * t194, t158 * t197 + t161 * t194, t248 + t197 * (-t182 + t283), -t219 * t197, t194 * (t184 - t283) + t247, 0, pkin(2) * t161 + pkin(8) * t125 - t197 * t85, -pkin(2) * t158 + pkin(8) * t126 + t194 * t85, pkin(2) * t163 + pkin(8) * t162 + t33, -pkin(2) * t85 + pkin(8) * t33, t225, -t323, -t346, t226, t348, t214, t194 * (t276 - t337) + t197 * (t27 - t339) - t340 + pkin(8) * t35, t194 * (t268 - t341) + t197 * (t28 - t342) + t359, pkin(8) * t31 - t14 * t194 + t218 * t42, pkin(8) * t10 + t14 * t218, t139 + t194 * (t153 * t196 - t155 * t193) * t177, t346, -t197 * t76 - t298, t225, t194 * (-t196 * t295 - t273) + t314, t226, t194 * (-t19 * t193 + t196 * t20 - t352) + t197 * (-t233 - t353) + t358, t194 * (-t193 * t22 + t257 * t295 + t337) + t197 * (-t203 + t339) - t360, t194 * (t196 * t21 + t280 * t83 + t341) + t197 * (-t211 + t361) - t359, t194 * (-pkin(9) * t7 + (-t257 + t280) * t25) + t197 * (-pkin(3) * t7 - t222) - pkin(2) * t7 + pkin(8) * t3, t214, -t348, -t346, t226, t323, t225, t194 * (t12 * t196 - t193 * t9 - t352) + t197 * (-t289 - t353) + t358, t194 * (t13 * t196 - t193 * t32 + t341) + t197 * (-t209 - t288 + t361) + t343 + pkin(8) * t34, t194 * (-t11 * t193 + t196 * t47 - t337) + t197 * (t16 + t286 - t339) + t360, t194 * (-pkin(9) * t4 - t193 * t2 + t196 * t6) + t197 * (-pkin(3) * t4 - t287) - pkin(2) * t4 + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t182 - t184, t238, t176, t237, qJDD(3), -t66, -t67, 0, 0, t71, t310, t313, t70, t307, t216, -pkin(3) * t75 - t268 + t336, pkin(3) * t83 + t276 - t338, pkin(9) * t45 + t15 + t322, -pkin(3) * t48 + pkin(9) * t15, t216, -t313, -t307, t71, -t193 * t295 + t266, t70, t19 * t196 + t193 * t20 + t347, t196 * t22 + t231 * t295 - t336, t338 + t193 * t21 - (pkin(3) + t279) * t83, pkin(9) * t8 + (-t231 - t279) * t25, t216, -t307, t313, t70, -t310, t71, t12 * t193 + t196 * t9 + t347, pkin(3) * t292 + t13 * t193 + t196 * t32 + t338, -pkin(3) * t295 + t11 * t196 + t193 * t47 + t336, -pkin(3) * t18 + pkin(9) * t5 + t193 * t6 + t196 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, -t117, -t84, -t251, t78, -t152, -t27, -t28, 0, 0, -t152, t84, t76, t251, -t117, -t251, t233, t203, t169 + t211, t222, -t152, -t78, -t84, -t251, t117, t251, t289, t17 + t288, pkin(5) * t252 - t201 - t286, t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t212, t291, t24, 0, 0, 0, 0, 0, 0, -t84, t291, t212, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t213, t293, t17;];
tauJ_reg  = t26;

% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:22:44
% EndTime: 2019-05-05 21:22:57
% DurationCPUTime: 4.61s
% Computational Cost: add. (7215->353), mult. (13952->396), div. (0->0), fcn. (8734->8), ass. (0->213)
t191 = sin(pkin(9));
t196 = sin(qJ(3));
t195 = sin(qJ(4));
t198 = cos(qJ(4));
t244 = qJD(1) * t196;
t154 = -t198 * qJD(3) + t195 * t244;
t199 = cos(qJ(3));
t240 = qJD(1) * qJD(3);
t234 = t199 * t240;
t239 = t196 * qJDD(1);
t161 = t234 + t239;
t217 = -qJDD(3) * t195 - t161 * t198;
t209 = qJD(4) * t154 + t217;
t177 = qJD(1) * t199 - qJD(4);
t251 = t177 * t154;
t295 = t251 - t209;
t156 = qJD(3) * t195 + t198 * t244;
t152 = t156 ^ 2;
t285 = t177 ^ 2;
t110 = -t285 - t152;
t120 = t156 * t154;
t181 = t196 * t240;
t238 = t199 * qJDD(1);
t162 = -t181 + t238;
t153 = -qJDD(4) + t162;
t301 = -t120 + t153;
t309 = t198 * t301;
t344 = t110 * t195 - t309;
t354 = t199 * t344;
t349 = -t196 * t295 + t354;
t192 = cos(pkin(9));
t313 = t195 * t301;
t56 = t110 * t198 + t313;
t356 = t192 * t56;
t360 = pkin(2) * t56;
t363 = t360 + pkin(7) * t349 + pkin(1) * (t191 * t349 + t356);
t359 = pkin(3) * t56;
t358 = pkin(8) * t56;
t357 = pkin(8) * t344;
t355 = t196 * t344;
t276 = t195 * t295;
t286 = t154 ^ 2;
t297 = t152 - t286;
t229 = -t198 * qJDD(3) + t161 * t195;
t107 = qJD(4) * t156 + t229;
t252 = t156 * t177;
t303 = t107 - t252;
t341 = t199 * t297 + t196 * (t198 * t303 + t276);
t89 = -t286 - t152;
t338 = pkin(3) * t89;
t307 = t251 + t209;
t275 = t195 * t307;
t73 = (qJD(4) + t177) * t156 + t229;
t42 = t198 * t73 + t275;
t353 = -pkin(8) * t42 - t338;
t352 = pkin(3) * t295 + t357;
t332 = t199 * t89;
t351 = t196 * t42 + t332;
t350 = t199 * t295 + t355;
t333 = t196 * t89;
t26 = t199 * t42 - t333;
t293 = t286 - t285;
t302 = t107 + t252;
t319 = t196 * (t198 * t293 + t313) + t199 * t302;
t127 = t152 - t285;
t300 = t120 + t153;
t310 = t198 * t300;
t320 = t196 * (t127 * t195 - t310) + t199 * t307;
t294 = -t285 - t286;
t325 = t195 * t294 - t310;
t314 = t195 * t300;
t324 = t198 * t294 + t314;
t340 = t196 * t303 + t199 * t324;
t339 = pkin(1) * (t191 * t340 - t192 * t325) + pkin(7) * t340 - pkin(2) * t325;
t336 = pkin(3) * t325;
t335 = pkin(8) * t325;
t347 = t307 * qJ(6);
t346 = t127 * t198 + t314;
t266 = t198 * t295;
t327 = -t195 * t303 + t266;
t343 = -t195 * t293 + t309;
t323 = -pkin(3) * t303 + pkin(8) * t324;
t322 = t196 * t324 - t199 * t303;
t334 = qJ(5) * t89;
t331 = qJ(5) * t295;
t265 = t198 * t307;
t39 = t195 * t73 - t265;
t246 = -g(3) + qJDD(2);
t180 = t199 * t246;
t201 = qJD(1) ^ 2;
t197 = sin(qJ(1));
t200 = cos(qJ(1));
t232 = t197 * g(1) - g(2) * t200;
t157 = qJDD(1) * pkin(1) + t232;
t222 = g(1) * t200 + g(2) * t197;
t158 = -pkin(1) * t201 - t222;
t245 = t191 * t157 + t192 * t158;
t101 = -pkin(2) * t201 + qJDD(1) * pkin(7) + t245;
t225 = -pkin(3) * t199 - pkin(8) * t196;
t228 = t201 * t225 + t101;
t283 = qJD(3) ^ 2;
t63 = -qJDD(3) * pkin(3) - t283 * pkin(8) + t196 * t228 - t180;
t318 = t107 * pkin(4) - t331 + t63;
t317 = qJ(5) * t301;
t304 = -pkin(4) * t300 + qJ(5) * t294;
t243 = qJD(5) * t177;
t168 = -0.2e1 * t243;
t242 = qJD(6) * t154;
t299 = 0.2e1 * t242 + t168;
t169 = 0.2e1 * t243;
t298 = -0.2e1 * t242 + t169;
t114 = t199 * t120;
t249 = t177 * t198;
t237 = t154 * t249;
t208 = t196 * (t107 * t195 - t237) + t114;
t221 = pkin(5) * t177 - qJ(6) * t156;
t292 = t156 * t221 + qJDD(6);
t290 = -pkin(5) * t107 + t292;
t230 = t192 * t157 - t191 * t158;
t100 = -qJDD(1) * pkin(2) - t201 * pkin(7) - t230;
t219 = -t162 + t181;
t220 = t161 + t234;
t55 = pkin(3) * t219 - pkin(8) * t220 + t100;
t231 = t196 * t246;
t64 = -t283 * pkin(3) + qJDD(3) * pkin(8) + t199 * t228 + t231;
t30 = t195 * t64 - t198 * t55;
t215 = t153 * pkin(4) - qJ(5) * t285 + qJDD(5) + t30;
t206 = t153 * pkin(5) + t215 + t347;
t113 = pkin(4) * t154 - qJ(5) * t156;
t233 = -pkin(5) * t154 - t113;
t216 = (-0.2e1 * qJD(6) - t233) * t156;
t14 = t216 + t206;
t31 = t195 * t55 + t198 * t64;
t226 = -pkin(4) * t285 - t153 * qJ(5) - t154 * t113 + t31;
t213 = -pkin(5) * t286 + t107 * qJ(6) - t177 * t221 + t226;
t15 = t213 + t299;
t282 = pkin(4) + pkin(5);
t289 = qJ(5) * t15 - t282 * t14;
t288 = qJ(5) * t73 - t282 * t307;
t250 = t177 * t195;
t124 = t156 * t250;
t260 = t196 * (-t198 * t209 + t124) - t114;
t287 = -t110 * t282 + t213 - t317;
t284 = 0.2e1 * t156;
t278 = t195 * t63;
t268 = t198 * t63;
t258 = qJ(5) * t195;
t257 = qJ(5) * t198;
t256 = t113 * t156;
t176 = t199 * t201 * t196;
t167 = qJDD(3) + t176;
t248 = t196 * t167;
t166 = -t176 + qJDD(3);
t247 = t199 * t166;
t235 = pkin(1) * t191 + pkin(7);
t17 = t195 * t30 + t198 * t31;
t85 = t101 * t196 - t180;
t86 = t199 * t101 + t231;
t45 = t196 * t85 + t199 * t86;
t227 = -t198 * t107 - t154 * t250;
t66 = -t156 * t249 - t195 * t209;
t23 = t168 + t226;
t24 = t215 + t256;
t224 = -pkin(4) * t24 + qJ(5) * t23;
t223 = pkin(4) * t307 - qJ(5) * t302;
t218 = t195 * t31 - t198 * t30;
t214 = (t154 * t195 + t156 * t198) * t177;
t212 = -pkin(1) * t192 - pkin(2) + t225;
t137 = t199 * t153;
t210 = t196 * (-t124 + t237) + t137;
t207 = -pkin(4) * t110 + t226 - t317;
t204 = qJD(5) * t284 - t318;
t203 = -t24 + t304;
t202 = -t206 + t304;
t25 = (-pkin(4) * t177 - 0.2e1 * qJD(5)) * t156 + t318;
t22 = (-t303 + t252) * pkin(4) + t204;
t21 = pkin(4) * t252 + t204 + t331;
t188 = t199 ^ 2;
t187 = t196 ^ 2;
t185 = t188 * t201;
t183 = t187 * t201;
t173 = -t185 - t283;
t172 = -t183 - t283;
t165 = t183 + t185;
t164 = (t187 + t188) * qJDD(1);
t163 = -0.2e1 * t181 + t238;
t160 = 0.2e1 * t234 + t239;
t143 = qJD(6) * t284;
t123 = -t172 * t196 - t247;
t122 = t173 * t199 - t248;
t82 = (qJD(4) - t177) * t154 + t217;
t46 = -qJ(5) * t303 - qJ(6) * t300;
t43 = -t198 * t302 - t275;
t40 = -t195 * t302 + t265;
t36 = -t196 * t82 - t354;
t32 = qJ(6) * t301 + t282 * t295;
t27 = t199 * t43 + t333;
t20 = t24 - t334;
t19 = -pkin(4) * t89 + t23;
t18 = qJ(6) * t286 + t25 - t290;
t13 = t21 + (-t110 - t286) * qJ(6) + t290;
t11 = t233 * t156 + t143 - t206 + t334 - t347;
t10 = (-t294 - t286) * qJ(6) + (-t107 - t303) * pkin(5) + t22 + t292;
t9 = t195 * t24 + t198 * t23;
t8 = t195 * t23 - t198 * t24;
t7 = -qJ(6) * t73 + t282 * t89 - t213 + t298;
t6 = -qJ(5) * t18 - qJ(6) * t14;
t5 = t196 * t25 + t199 * t9;
t4 = t14 * t195 + t15 * t198;
t3 = -t14 * t198 + t15 * t195;
t2 = -qJ(6) * t15 - t282 * t18;
t1 = t18 * t196 + t199 * t4;
t12 = [0, 0, 0, 0, 0, qJDD(1), t232, t222, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t192 - t191 * t201) + t230, pkin(1) * (-qJDD(1) * t191 - t192 * t201) - t245, 0, pkin(1) * (t191 * t245 + t192 * t230), t220 * t196, t160 * t199 + t163 * t196, t248 + t199 * (-t183 + t283), -t219 * t199, t196 * (t185 - t283) + t247, 0, -t199 * t100 + pkin(2) * t163 + pkin(7) * t122 + pkin(1) * (t122 * t191 + t163 * t192), t196 * t100 - pkin(2) * t160 + pkin(7) * t123 + pkin(1) * (t123 * t191 - t160 * t192), pkin(2) * t165 + pkin(7) * t164 + pkin(1) * (t164 * t191 + t165 * t192) + t45, -pkin(2) * t100 + pkin(7) * t45 + pkin(1) * (-t100 * t192 + t191 * t45), t260, -t341, t320, t208, t319, t210, t196 * (t278 - t335) + t199 * (t30 - t336) + t339, t196 * (t268 - t358) + t199 * (t31 - t359) - t360 + pkin(7) * t36 + pkin(1) * (t191 * t36 - t356), -t196 * t218 - t212 * t39 - t235 * t26, t235 * (t17 * t199 + t196 * t63) + t212 * t218, t260, t320, t341, t210, -t319, t208, t196 * (-t195 * t22 - t257 * t303 - t335) + t199 * (-t203 - t336) + t339, t196 * (-pkin(8) * t40 - t19 * t195 + t198 * t20) + t199 * (-pkin(3) * t40 - t223) - pkin(2) * t40 + pkin(7) * t27 + pkin(1) * (t191 * t27 - t192 * t40), t196 * (-pkin(4) * t276 + t198 * t21 + t358) + t199 * (t169 - t207 + t359) + t363, t196 * (-pkin(8) * t8 + (pkin(4) * t195 - t257) * t25) + t199 * (-pkin(3) * t8 - t224) - pkin(2) * t8 + pkin(7) * t5 + pkin(1) * (t191 * t5 - t192 * t8), t260, t341, -t320, t208, t319, t137 + t196 * (t154 * t198 - t156 * t195) * t177, t196 * (-t10 * t195 + t198 * t46 - t335) + (pkin(5) * t300 - t202 + t216 - t336) * t199 + t339, t196 * (t13 * t198 - t195 * t32 + t358) + t199 * (-t287 + t298 + t359) + t363, t196 * (-pkin(8) * t39 + t11 * t198 - t195 * t7) + t199 * (-pkin(3) * t39 - t288) - pkin(2) * t39 + pkin(7) * t26 + pkin(1) * (t191 * t26 - t192 * t39), t196 * (-pkin(8) * t3 - t195 * t2 + t198 * t6) + t199 * (-pkin(3) * t3 - t289) - pkin(2) * t3 + pkin(7) * t1 + pkin(1) * (t1 * t191 - t192 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, 0, 0, 0, 0, 0, 0, t167 * t199 + t173 * t196, -t166 * t196 + t172 * t199, 0, t196 * t86 - t199 * t85, 0, 0, 0, 0, 0, 0, t322, t199 * t82 - t355, -t351, t17 * t196 - t199 * t63, 0, 0, 0, 0, 0, 0, t322, t196 * t43 - t332, t350, t196 * t9 - t199 * t25, 0, 0, 0, 0, 0, 0, t322, t350, t351, -t18 * t199 + t196 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t183 - t185, t239, t176, t238, qJDD(3), -t85, -t86, 0, 0, t66, t327, -t346, t227, -t343, t214, -t268 + t323, pkin(3) * t82 + t278 - t357, t17 + t353, -pkin(3) * t63 + pkin(8) * t17, t66, -t346, -t327, t214, t343, t227, t198 * t22 - t258 * t303 + t323, pkin(8) * t43 + t19 * t198 + t195 * t20 - t338, pkin(4) * t266 + t195 * t21 + t352, pkin(8) * t9 + (-pkin(4) * t198 - pkin(3) - t258) * t25, t66, -t327, t346, t227, -t343, t214, t10 * t198 + t195 * t46 + t323, t13 * t195 + t198 * t32 + t352, t11 * t195 + t198 * t7 - t353, -pkin(3) * t18 + pkin(8) * t4 + t195 * t6 + t198 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t297, -t307, -t120, -t302, -t153, -t30, -t31, 0, 0, t120, -t307, -t297, -t153, t302, -t120, t203, t223, t168 + t207, t224, t120, -t297, t307, -t120, -t302, -t153, -t256 + t143 + (-t300 - t120) * pkin(5) + t202, t287 + t299, t288, t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t300, -t307, t110, t24, 0, 0, 0, 0, 0, 0, t300, t110, t307, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t303, t295, t89, -t18;];
tauJ_reg  = t12;

% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:39:23
% EndTime: 2019-05-05 17:39:42
% DurationCPUTime: 8.31s
% Computational Cost: add. (17706->409), mult. (37034->543), div. (0->0), fcn. (25022->10), ass. (0->257)
t240 = sin(pkin(9));
t242 = cos(pkin(9));
t245 = sin(qJ(3));
t248 = cos(qJ(3));
t239 = sin(pkin(10));
t241 = cos(pkin(10));
t290 = qJD(1) * t245;
t205 = -t241 * qJD(3) + t239 * t290;
t206 = t239 * qJD(3) + t241 * t290;
t244 = sin(qJ(5));
t247 = cos(qJ(5));
t177 = t247 * t205 + t244 * t206;
t287 = qJD(1) * qJD(3);
t280 = t248 * t287;
t286 = t245 * qJDD(1);
t212 = t280 + t286;
t185 = t239 * qJDD(3) + t241 * t212;
t261 = t241 * qJDD(3) - t239 * t212;
t253 = -t177 * qJD(5) + t247 * t185 + t244 * t261;
t288 = t248 * qJD(1);
t225 = -qJD(5) + t288;
t314 = t177 * t225;
t334 = t253 + t314;
t229 = t245 * t287;
t285 = t248 * qJDD(1);
t213 = -t229 + t285;
t207 = -qJDD(5) + t213;
t179 = -t244 * t205 + t247 * t206;
t313 = t179 * t177;
t255 = t207 - t313;
t301 = t244 * t255;
t176 = t179 ^ 2;
t322 = t225 ^ 2;
t336 = -t176 - t322;
t75 = -t247 * t336 - t301;
t296 = t247 * t255;
t82 = -t244 * t336 + t296;
t55 = t239 * t75 + t241 * t82;
t42 = t245 * t334 + t248 * t55;
t53 = t239 * t82 - t241 * t75;
t394 = -pkin(1) * (t240 * t42 - t242 * t53) + pkin(2) * t53 - pkin(7) * t42;
t390 = qJ(4) * t53;
t388 = -pkin(3) * t53 + pkin(4) * t75;
t387 = pkin(3) * t334 - qJ(4) * t55;
t386 = t245 * t55 - t248 * t334;
t163 = t179 * t225;
t276 = -t244 * t185 + t247 * t261;
t259 = t179 * qJD(5) - t276;
t102 = t163 + t259;
t324 = t177 ^ 2;
t154 = t324 - t322;
t93 = t244 * t154 - t296;
t97 = t247 * t154 + t301;
t385 = -t248 * t102 + t245 * (t239 * t93 - t241 * t97);
t335 = t176 - t324;
t337 = -t163 + t259;
t61 = -t244 * t337 + t247 * t334;
t303 = t244 * t334;
t63 = t247 * t337 + t303;
t384 = t245 * (t239 * t61 + t241 * t63) + t248 * t335;
t382 = pkin(8) * t75;
t381 = pkin(8) * t82;
t380 = t239 * t63 - t241 * t61;
t378 = t239 * t97 + t241 * t93;
t129 = t207 + t313;
t300 = t244 * t129;
t332 = -t322 - t324;
t340 = t247 * t332 + t300;
t125 = t247 * t129;
t341 = t244 * t332 - t125;
t352 = t239 * t340 + t241 * t341;
t353 = -t239 * t341 + t241 * t340;
t367 = t245 * t337 + t248 * t353;
t376 = pkin(1) * (t240 * t367 - t242 * t352) + pkin(7) * t367 - pkin(2) * t352;
t372 = qJ(4) * t352;
t370 = -pkin(3) * t352 - pkin(4) * t341;
t369 = -pkin(3) * t337 + qJ(4) * t353;
t368 = t245 * t353 - t248 * t337;
t155 = -t176 + t322;
t354 = t247 * t155 - t300;
t355 = -t244 * t155 - t125;
t366 = t239 * t355 + t241 * t354;
t333 = -t314 + t253;
t365 = t245 * (-t239 * t354 + t241 * t355) - t248 * t333;
t113 = -t324 - t176;
t364 = pkin(3) * t113;
t363 = pkin(4) * t113;
t361 = pkin(8) * t340;
t360 = pkin(8) * t341;
t358 = qJ(6) * t334;
t293 = -g(3) + qJDD(2);
t228 = t248 * t293;
t250 = qJD(1) ^ 2;
t246 = sin(qJ(1));
t249 = cos(qJ(1));
t278 = t246 * g(1) - t249 * g(2);
t208 = qJDD(1) * pkin(1) + t278;
t269 = t249 * g(1) + t246 * g(2);
t210 = -t250 * pkin(1) - t269;
t291 = t240 * t208 + t242 * t210;
t172 = -t250 * pkin(2) + qJDD(1) * pkin(7) + t291;
t265 = -t248 * pkin(3) - t245 * qJ(4);
t273 = t250 * t265 + t172;
t321 = qJD(3) ^ 2;
t139 = -qJDD(3) * pkin(3) - t321 * qJ(4) + t273 * t245 + qJDD(4) - t228;
t186 = -pkin(4) * t288 - t206 * pkin(8);
t323 = t205 ^ 2;
t79 = -t261 * pkin(4) - t323 * pkin(8) + t206 * t186 + t139;
t359 = pkin(5) * t259 - t358 + t79;
t357 = t245 * t113;
t356 = t248 * t113;
t311 = t206 * t205;
t257 = -t213 - t311;
t339 = t239 * t257;
t338 = t241 * t257;
t190 = t205 * t288;
t167 = -t185 + t190;
t191 = t206 * t288;
t165 = t261 - t191;
t140 = t177 * pkin(5) - t179 * qJ(6);
t275 = t242 * t208 - t240 * t210;
t171 = -qJDD(1) * pkin(2) - t250 * pkin(7) - t275;
t263 = t212 + t280;
t135 = -t263 * qJ(4) + (-t213 + t229) * pkin(3) + t171;
t277 = t245 * t293;
t141 = -t321 * pkin(3) + qJDD(3) * qJ(4) + t273 * t248 + t277;
t73 = 0.2e1 * qJD(4) * t206 - t241 * t135 + t239 * t141;
t69 = t257 * pkin(4) + t167 * pkin(8) - t73;
t74 = -0.2e1 * qJD(4) * t205 + t239 * t135 + t241 * t141;
t71 = -t323 * pkin(4) + t261 * pkin(8) + t186 * t288 + t74;
t36 = t244 * t69 + t247 * t71;
t274 = -t207 * qJ(6) - t177 * t140 + t36;
t329 = -pkin(5) * (t336 + t322) - qJ(6) * t255 + t274;
t256 = (t177 * t244 + t179 * t247) * t225;
t310 = t225 * t244;
t153 = t179 * t310;
t309 = t225 * t247;
t284 = t177 * t309;
t267 = -t153 + t284;
t328 = t239 * t267 + t241 * t256;
t260 = t244 * t259 - t284;
t268 = -t177 * t310 - t247 * t259;
t327 = t239 * t260 + t241 * t268;
t326 = t245 * (-t239 * t256 + t241 * t267) + t248 * t207;
t283 = t248 * t313;
t325 = t245 * (-t239 * t268 + t241 * t260) + t283;
t204 = t206 ^ 2;
t320 = pkin(5) * t247;
t35 = t244 * t71 - t247 * t69;
t16 = t244 * t36 - t247 * t35;
t319 = t239 * t16;
t318 = t241 * t16;
t317 = t244 * t79;
t316 = t247 * t79;
t315 = qJ(6) * t247;
t308 = t239 * t139;
t168 = t213 - t311;
t307 = t239 * t168;
t306 = t241 * t139;
t305 = t241 * t168;
t302 = t244 * t333;
t224 = t248 * t250 * t245;
t218 = qJDD(3) + t224;
t298 = t245 * t218;
t217 = -t224 + qJDD(3);
t294 = t248 * t217;
t289 = qJD(6) * t225;
t282 = t248 * t311;
t281 = pkin(1) * t240 + pkin(7);
t279 = -qJ(6) * t244 - pkin(4);
t46 = t239 * t73 + t241 * t74;
t17 = t244 * t35 + t247 * t36;
t150 = t245 * t172 - t228;
t151 = t248 * t172 + t277;
t112 = t245 * t150 + t248 * t151;
t89 = -t179 * t309 + t244 * t253;
t90 = t247 * t253 + t153;
t271 = t245 * (-t239 * t89 + t241 * t90) - t283;
t219 = -0.2e1 * t289;
t266 = t219 + t274;
t26 = -pkin(5) * t322 + t266;
t27 = t207 * pkin(5) - qJ(6) * t322 + t179 * t140 + qJDD(6) + t35;
t270 = -pkin(5) * t27 + qJ(6) * t26;
t264 = -pkin(5) * t333 - qJ(6) * t102;
t262 = t239 * t74 - t241 * t73;
t254 = -pkin(1) * t242 - pkin(2) + t265;
t252 = -pkin(5) * t129 + qJ(6) * t332 - t27;
t251 = 0.2e1 * qJD(6) * t179 - t359;
t236 = t248 ^ 2;
t235 = t245 ^ 2;
t232 = t236 * t250;
t231 = t235 * t250;
t222 = -t232 - t321;
t221 = -t231 - t321;
t216 = t231 + t232;
t215 = (t235 + t236) * qJDD(1);
t214 = -0.2e1 * t229 + t285;
t211 = 0.2e1 * t280 + t286;
t200 = t248 * t213;
t189 = -t204 - t232;
t188 = -t204 + t232;
t187 = -t232 + t323;
t183 = -t245 * t221 - t294;
t182 = t248 * t222 - t298;
t180 = -t232 - t323;
t166 = t185 + t190;
t164 = -t191 - t261;
t158 = -t204 - t323;
t146 = -t239 * t189 + t305;
t145 = t241 * t189 + t307;
t134 = t241 * t180 - t339;
t133 = t239 * t180 + t338;
t123 = t241 * t165 - t239 * t167;
t110 = t248 * t146 + t245 * t166;
t109 = t248 * t134 + t245 * t164;
t103 = (-qJD(5) - t225) * t179 + t276;
t99 = t247 * t333;
t66 = t247 * t103 + t302;
t64 = -t247 * t102 + t302;
t62 = t244 * t103 - t99;
t60 = -t244 * t102 - t99;
t58 = t239 * t90 + t241 * t89;
t52 = t316 + t382;
t51 = t317 - t360;
t44 = -pkin(4) * t334 + t317 + t381;
t40 = (-pkin(5) * t225 - 0.2e1 * qJD(6)) * t179 + t359;
t39 = -pkin(4) * t337 - t316 + t361;
t33 = -t239 * t62 + t241 * t66;
t32 = -t239 * t60 + t241 * t64;
t31 = t239 * t66 + t241 * t62;
t30 = t239 * t64 + t241 * t60;
t29 = t251 + (-t337 + t163) * pkin(5);
t28 = pkin(5) * t163 + t251 + t358;
t25 = t248 * t33 + t357;
t24 = t248 * t32 + t357;
t23 = -qJ(6) * t113 + t27;
t22 = (-t113 - t322) * pkin(5) + t266;
t21 = -t244 * t29 - t315 * t337 - t360;
t20 = -pkin(5) * t303 + t247 * t28 - t382;
t19 = t247 * t29 + t279 * t337 + t361;
t18 = -t381 + t244 * t28 + (pkin(4) + t320) * t334;
t15 = -pkin(4) * t79 + pkin(8) * t17;
t14 = -pkin(8) * t62 - t16;
t13 = t244 * t27 + t247 * t26;
t12 = t244 * t26 - t247 * t27;
t11 = pkin(8) * t66 + t17 - t363;
t10 = -pkin(8) * t60 - t244 * t22 + t247 * t23;
t9 = pkin(8) * t64 + t247 * t22 + t244 * t23 - t363;
t8 = t241 * t17 - t319;
t7 = t239 * t17 + t318;
t6 = -pkin(8) * t12 + (pkin(5) * t244 - t315) * t40;
t5 = t245 * t79 + t248 * t8;
t4 = -t239 * t12 + t241 * t13;
t3 = t241 * t12 + t239 * t13;
t2 = pkin(8) * t13 + (t279 - t320) * t40;
t1 = t245 * t40 + t248 * t4;
t34 = [0, 0, 0, 0, 0, qJDD(1), t278, t269, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t242 * qJDD(1) - t240 * t250) + t275, pkin(1) * (-t240 * qJDD(1) - t242 * t250) - t291, 0, pkin(1) * (t240 * t291 + t242 * t275), t263 * t245, t248 * t211 + t245 * t214, t298 + t248 * (-t231 + t321), -t245 * t280 + t200, t245 * (t232 - t321) + t294, 0, -t248 * t171 + pkin(2) * t214 + pkin(7) * t182 + pkin(1) * (t240 * t182 + t242 * t214), t245 * t171 - pkin(2) * t211 + pkin(7) * t183 + pkin(1) * (t240 * t183 - t242 * t211), pkin(2) * t216 + pkin(7) * t215 + pkin(1) * (t240 * t215 + t242 * t216) + t112, -pkin(2) * t171 + pkin(7) * t112 + pkin(1) * (t240 * t112 - t242 * t171), t245 * (t241 * t185 + t239 * t191) - t282, t245 * (-t241 * t164 - t239 * t166) + t248 * (-t204 + t323), t245 * (-t239 * t188 + t338) + t248 * t167, t245 * (-t241 * t190 - t239 * t261) + t282, t245 * (t241 * t187 + t307) - t248 * t165, t200 + t245 * (t205 * t241 - t206 * t239) * t288, t245 * (-qJ(4) * t133 + t308) + t248 * (-pkin(3) * t133 + t73) - pkin(2) * t133 + pkin(7) * t109 + pkin(1) * (t240 * t109 - t242 * t133), t245 * (-qJ(4) * t145 + t306) + t248 * (-pkin(3) * t145 + t74) - pkin(2) * t145 + pkin(7) * t110 + pkin(1) * (t240 * t110 - t242 * t145), -t245 * t262 + t281 * (t248 * t123 + t245 * t158) + t254 * (t239 * t165 + t241 * t167), t281 * (t245 * t139 + t248 * t46) + t254 * t262, t271, -t384, t365, t325, -t385, t326, t245 * (-t239 * t39 + t241 * t51 - t372) + t248 * (t35 + t370) + t376, t245 * (-t239 * t44 + t241 * t52 - t390) + t248 * (t36 + t388) - t394, t245 * (-qJ(4) * t31 - t239 * t11 + t241 * t14) + t248 * (-pkin(3) * t31 - pkin(4) * t62) - pkin(2) * t31 + pkin(7) * t25 + pkin(1) * (t240 * t25 - t242 * t31), t245 * (-pkin(8) * t318 - qJ(4) * t7 - t239 * t15) + t248 * (-pkin(3) * t7 - pkin(4) * t16) - pkin(2) * t7 + pkin(7) * t5 + pkin(1) * (t240 * t5 - t242 * t7), t271, t365, t384, t326, t385, t325, t245 * (-t239 * t19 + t241 * t21 - t372) + t248 * (-t252 + t370) + t376, t245 * (-qJ(4) * t30 + t241 * t10 - t239 * t9) + t248 * (-pkin(3) * t30 - pkin(4) * t60 - t264) - pkin(2) * t30 + pkin(7) * t24 + pkin(1) * (t240 * t24 - t242 * t30), t245 * (-t239 * t18 + t241 * t20 + t390) + t248 * (0.2e1 * t289 - t329 - t388) + t394, t245 * (-qJ(4) * t3 - t239 * t2 + t241 * t6) + t248 * (-pkin(3) * t3 - pkin(4) * t12 - t270) - pkin(2) * t3 + pkin(7) * t1 + pkin(1) * (t240 * t1 - t242 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, 0, 0, 0, 0, 0, 0, t248 * t218 + t245 * t222, -t245 * t217 + t248 * t221, 0, -t248 * t150 + t245 * t151, 0, 0, 0, 0, 0, 0, t245 * t134 - t248 * t164, t245 * t146 - t248 * t166, t245 * t123 - t248 * t158, -t248 * t139 + t245 * t46, 0, 0, 0, 0, 0, 0, t368, t386, t245 * t33 - t356, t245 * t8 - t248 * t79, 0, 0, 0, 0, 0, 0, t368, t245 * t32 - t356, -t386, t245 * t4 - t248 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224, t231 - t232, t286, t224, t285, qJDD(3), -t150, -t151, 0, 0, t239 * t185 - t241 * t191, -t239 * t164 + t241 * t166, t241 * t188 + t339, -t239 * t190 + t241 * t261, t239 * t187 - t305, (t205 * t239 + t206 * t241) * t288, -pkin(3) * t164 + qJ(4) * t134 - t306, -pkin(3) * t166 + qJ(4) * t146 + t308, -pkin(3) * t158 + qJ(4) * t123 + t46, -pkin(3) * t139 + qJ(4) * t46, t58, -t380, t366, t327, t378, t328, t239 * t51 + t241 * t39 + t369, t239 * t52 + t241 * t44 - t387, qJ(4) * t33 + t241 * t11 + t239 * t14 - t364, -pkin(3) * t79 - pkin(8) * t319 + qJ(4) * t8 + t241 * t15, t58, t366, t380, t328, -t378, t327, t241 * t19 + t239 * t21 + t369, qJ(4) * t32 + t239 * t10 + t241 * t9 - t364, t241 * t18 + t239 * t20 + t387, -pkin(3) * t40 + qJ(4) * t4 + t241 * t2 + t239 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t166, t158, t139, 0, 0, 0, 0, 0, 0, t337, t334, t113, t79, 0, 0, 0, 0, 0, 0, t337, t113, -t334, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t313, t335, t333, -t313, -t102, -t207, -t35, -t36, 0, 0, t313, t333, -t335, -t207, t102, -t313, t252, t264, t219 + t329, t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, t333, t336, t27;];
tauJ_reg  = t34;

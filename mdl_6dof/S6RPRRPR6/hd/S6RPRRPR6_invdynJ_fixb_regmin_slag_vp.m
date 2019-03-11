% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:17:46
% EndTime: 2019-03-09 05:18:01
% DurationCPUTime: 6.47s
% Computational Cost: add. (9399->472), mult. (22637->621), div. (0->0), fcn. (18095->16), ass. (0->250)
t244 = cos(pkin(10));
t253 = cos(qJ(3));
t324 = t253 * t244;
t222 = qJD(1) * t324;
t242 = sin(pkin(10));
t249 = sin(qJ(3));
t334 = t242 * t249;
t302 = qJD(1) * t334;
t191 = t222 - t302;
t178 = qJD(4) - t191;
t202 = t242 * t253 + t244 * t249;
t193 = t202 * qJD(1);
t248 = sin(qJ(4));
t252 = cos(qJ(4));
t312 = t252 * qJD(3);
t155 = t193 * t248 - t312;
t157 = qJD(3) * t248 + t193 * t252;
t241 = sin(pkin(11));
t243 = cos(pkin(11));
t106 = t155 * t241 - t157 * t243;
t247 = sin(qJ(6));
t251 = cos(qJ(6));
t313 = qJD(6) * t247;
t281 = -t155 * t243 - t157 * t241;
t326 = t251 * t281;
t309 = t244 * qJDD(1);
t310 = t242 * qJDD(1);
t304 = qJD(3) * t222 + t249 * t309 + t253 * t310;
t144 = -qJD(3) * t302 + t304;
t315 = qJD(4) * t248;
t98 = qJD(4) * t312 + t248 * qJDD(3) + t252 * t144 - t193 * t315;
t99 = qJD(4) * t157 - t252 * qJDD(3) + t248 * t144;
t49 = -t241 * t98 - t243 * t99;
t50 = -t241 * t99 + t243 * t98;
t12 = qJD(6) * t326 + t106 * t313 + t247 * t49 + t251 * t50;
t174 = qJD(6) + t178;
t58 = t106 * t247 + t326;
t359 = t174 * t58;
t399 = t12 - t359;
t384 = -t251 * t106 + t247 * t281;
t398 = t384 * t58;
t201 = t241 * t252 + t243 * t248;
t387 = t178 * t201;
t275 = t241 * t248 - t243 * t252;
t386 = t178 * t275;
t397 = t384 ^ 2 - t58 ^ 2;
t235 = qJ(4) + pkin(11) + qJ(6);
t224 = sin(t235);
t225 = cos(t235);
t254 = cos(qJ(1));
t240 = pkin(10) + qJ(3);
t232 = cos(t240);
t250 = sin(qJ(1));
t337 = t232 * t250;
t159 = t224 * t254 - t225 * t337;
t336 = t232 * t254;
t161 = t224 * t250 + t225 * t336;
t195 = t202 * qJD(3);
t283 = t249 * t310 - t253 * t309;
t145 = qJD(1) * t195 + t283;
t138 = qJDD(4) + t145;
t227 = pkin(2) * t244 + pkin(1);
t211 = -qJD(1) * t227 + qJD(2);
t114 = -t191 * pkin(3) - t193 * pkin(8) + t211;
t363 = pkin(7) + qJ(2);
t212 = t363 * t242;
t206 = qJD(1) * t212;
t213 = t363 * t244;
t207 = qJD(1) * t213;
t149 = -t249 * t206 + t253 * t207;
t140 = qJD(3) * pkin(8) + t149;
t81 = t114 * t248 + t140 * t252;
t210 = -qJDD(1) * t227 + qJDD(2);
t94 = t145 * pkin(3) - t144 * pkin(8) + t210;
t87 = t252 * t94;
t311 = qJD(1) * qJD(2);
t372 = t363 * qJDD(1) + t311;
t171 = t372 * t242;
t172 = t372 * t244;
t280 = -t249 * t171 + t253 * t172;
t379 = -t206 * t253 - t249 * t207;
t92 = qJDD(3) * pkin(8) + qJD(3) * t379 + t280;
t16 = t138 * pkin(4) - t98 * qJ(5) - qJD(4) * t81 - t157 * qJD(5) - t248 * t92 + t87;
t314 = qJD(4) * t252;
t265 = t114 * t314 - t140 * t315 + t248 * t94 + t252 * t92;
t18 = -qJ(5) * t99 - qJD(5) * t155 + t265;
t4 = t243 * t16 - t18 * t241;
t2 = pkin(5) * t138 - pkin(9) * t50 + t4;
t67 = -qJ(5) * t155 + t81;
t357 = t243 * t67;
t80 = t252 * t114 - t140 * t248;
t66 = -qJ(5) * t157 + t80;
t53 = pkin(4) * t178 + t66;
t30 = t241 * t53 + t357;
t383 = pkin(9) * t281;
t21 = t30 + t383;
t20 = t21 * t313;
t231 = sin(t240);
t365 = g(3) * t231;
t139 = -qJD(3) * pkin(3) - t379;
t100 = pkin(4) * t155 + qJD(5) + t139;
t54 = -pkin(5) * t281 + t100;
t396 = g(1) * t161 - g(2) * t159 - t247 * t2 + t225 * t365 - t54 * t58 + t20;
t13 = qJD(6) * t384 + t247 * t50 - t251 * t49;
t355 = t384 * t174;
t394 = -t13 + t355;
t362 = qJ(5) + pkin(8);
t296 = qJD(4) * t362;
t141 = pkin(3) * t193 - pkin(8) * t191;
t321 = t248 * t141 + t252 * t379;
t343 = t191 * t248;
t393 = qJ(5) * t343 + t252 * qJD(5) - t248 * t296 - t321;
t130 = t252 * t141;
t392 = -pkin(4) * t193 - t130 + (qJ(5) * t191 - t296) * t252 + (-qJD(5) + t379) * t248;
t289 = g(1) * t254 + g(2) * t250;
t260 = -g(3) * t232 + t231 * t289;
t317 = qJD(3) * t253;
t318 = qJD(3) * t249;
t272 = -t253 * t171 - t249 * t172 + t206 * t318 - t207 * t317;
t93 = -qJDD(3) * pkin(3) - t272;
t391 = -qJD(4) * pkin(8) * t178 + t260 - t93;
t390 = t315 - t343;
t158 = t224 * t337 + t225 * t254;
t160 = -t224 * t336 + t225 * t250;
t5 = t241 * t16 + t243 * t18;
t3 = pkin(9) * t49 + t5;
t303 = t251 * t2 - t247 * t3;
t389 = -g(1) * t160 + g(2) * t158 + t224 * t365 - t54 * t384 + t303;
t388 = pkin(9) * t106;
t147 = t201 * t251 - t247 * t275;
t360 = qJD(6) * t147 - t247 * t386 + t251 * t387;
t349 = qJDD(1) * pkin(1);
t378 = g(1) * t250 - g(2) * t254;
t273 = -qJDD(2) + t349 + t378;
t200 = -t324 + t334;
t194 = t200 * qJD(3);
t301 = t202 * t314;
t385 = -t194 * t248 + t301;
t352 = -t241 * t393 + t243 * t392;
t351 = t241 * t392 + t243 * t393;
t381 = t231 * t378;
t151 = t212 * t253 + t249 * t213;
t377 = pkin(4) * t390 - t149;
t325 = t252 * t254;
t329 = t248 * t250;
t179 = t232 * t329 + t325;
t327 = t250 * t252;
t328 = t248 * t254;
t181 = -t232 * t328 + t327;
t376 = -g(1) * t181 + g(2) * t179;
t375 = qJ(2) * qJDD(1);
t135 = qJDD(6) + t138;
t279 = -t201 * t247 - t251 * t275;
t361 = qJD(6) * t279 - t247 * t387 - t251 * t386;
t374 = -t147 * t135 - t361 * t174;
t373 = t289 * t232 + t365;
t371 = pkin(4) * t241;
t115 = -t200 * qJD(2) - qJD(3) * t151;
t142 = pkin(3) * t195 + pkin(8) * t194;
t131 = t252 * t142;
t143 = pkin(3) * t200 - pkin(8) * t202 - t227;
t152 = -t212 * t249 + t213 * t253;
t150 = t252 * t152;
t274 = qJ(5) * t194 - qJD(5) * t202;
t37 = t195 * pkin(4) - t248 * t115 + t131 + t274 * t252 + (-t150 + (qJ(5) * t202 - t143) * t248) * qJD(4);
t306 = t252 * t115 + t248 * t142 + t143 * t314;
t41 = -qJ(5) * t301 + (-qJD(4) * t152 + t274) * t248 + t306;
t11 = t241 * t37 + t243 * t41;
t60 = t241 * t67;
t34 = t243 * t66 - t60;
t133 = t252 * t143;
t340 = t202 * t252;
t71 = pkin(4) * t200 - qJ(5) * t340 - t152 * t248 + t133;
t320 = t248 * t143 + t150;
t341 = t202 * t248;
t82 = -qJ(5) * t341 + t320;
t43 = t241 * t71 + t243 * t82;
t358 = t193 * t58;
t29 = t243 * t53 - t60;
t19 = pkin(5) * t178 + t29 + t388;
t356 = t251 * t19;
t354 = t384 * t193;
t353 = t98 * t248;
t350 = pkin(5) * t387 + t377;
t347 = t155 * t178;
t346 = t155 * t193;
t345 = t157 * t178;
t344 = t157 * t193;
t330 = t248 * t138;
t127 = t252 * t138;
t214 = t362 * t248;
t215 = t362 * t252;
t154 = -t241 * t214 + t243 * t215;
t319 = t242 ^ 2 + t244 ^ 2;
t316 = qJD(4) * t202;
t229 = pkin(4) * t252 + pkin(3);
t300 = pkin(4) * t248 + t363;
t298 = qJD(6) * t19 + t3;
t10 = -t241 * t41 + t243 * t37;
t33 = -t241 * t66 - t357;
t42 = -t241 * t82 + t243 * t71;
t295 = t319 * qJD(1) ^ 2;
t294 = -qJD(4) * t114 - t92;
t153 = -t243 * t214 - t215 * t241;
t292 = t178 * t252;
t116 = qJD(2) * t202 - t212 * t318 + t213 * t317;
t291 = 0.2e1 * t319;
t290 = t279 * t135 - t174 * t360;
t287 = -t140 * t314 + t87;
t123 = -pkin(9) * t275 + t154;
t286 = pkin(5) * t193 - pkin(9) * t386 + qJD(6) * t123 - t352;
t122 = -pkin(9) * t201 + t153;
t285 = pkin(9) * t387 - qJD(6) * t122 - t351;
t284 = pkin(4) * t341 + t151;
t7 = t247 * t19 + t251 * t21;
t125 = t201 * t202;
t126 = t275 * t202;
t282 = -t251 * t125 + t126 * t247;
t89 = -t125 * t247 - t126 * t251;
t277 = t229 * t232 + t231 * t362;
t271 = -t178 * t390 + t127;
t226 = pkin(4) * t243 + pkin(5);
t270 = t226 * t247 + t251 * t371;
t269 = t226 * t251 - t247 * t371;
t268 = pkin(4) * t385 + t116;
t267 = t227 + t277;
t266 = -t194 * t252 - t202 * t315;
t264 = -pkin(8) * t138 + t139 * t178;
t263 = t273 + t349;
t257 = t291 * t311 - t289;
t45 = t99 * pkin(4) + qJDD(5) + t93;
t182 = t232 * t325 + t329;
t180 = -t232 * t327 + t328;
t165 = pkin(5) * t275 - t229;
t91 = pkin(4) * t157 - pkin(5) * t106;
t90 = pkin(5) * t125 + t284;
t85 = -t194 * t275 + t201 * t316;
t84 = t194 * t201 + t275 * t316;
t44 = -pkin(5) * t84 + t268;
t32 = -pkin(9) * t125 + t43;
t31 = pkin(5) * t200 + pkin(9) * t126 + t42;
t27 = qJD(6) * t89 - t247 * t85 - t251 * t84;
t26 = qJD(6) * t282 + t247 * t84 - t251 * t85;
t24 = t34 + t388;
t23 = t33 - t383;
t22 = -pkin(5) * t49 + t45;
t9 = pkin(9) * t84 + t11;
t8 = pkin(5) * t195 + pkin(9) * t85 + t10;
t6 = -t21 * t247 + t356;
t1 = [qJDD(1), t378, t289, t263 * t244, -t263 * t242, t291 * t375 + t257, pkin(1) * t273 + (t319 * t375 + t257) * qJ(2), t144 * t202 - t193 * t194, -t144 * t200 - t145 * t202 - t191 * t194 - t193 * t195, -qJD(3) * t194 + qJDD(3) * t202, -qJD(3) * t195 - qJDD(3) * t200, 0, -t116 * qJD(3) - t151 * qJDD(3) - t227 * t145 + t211 * t195 + t210 * t200 + t232 * t378, -t115 * qJD(3) - t152 * qJDD(3) - t227 * t144 - t211 * t194 + t210 * t202 - t381, t157 * t266 + t340 * t98 -(-t155 * t252 - t157 * t248) * t194 + (-t353 - t252 * t99 + (t155 * t248 - t157 * t252) * qJD(4)) * t202, t127 * t202 + t157 * t195 + t178 * t266 + t98 * t200, -t155 * t195 - t178 * t385 - t99 * t200 - t202 * t330, t138 * t200 + t178 * t195 (-t152 * t314 + t131) * t178 + t133 * t138 + t287 * t200 + t80 * t195 + t116 * t155 + t151 * t99 + t139 * t301 - g(1) * t180 - g(2) * t182 + ((-qJD(4) * t143 - t115) * t178 - t152 * t138 + t294 * t200 + t93 * t202 - t139 * t194) * t248 -(-t152 * t315 + t306) * t178 - t320 * t138 - t265 * t200 - t81 * t195 + t116 * t157 + t151 * t98 + t93 * t340 - g(1) * t179 - g(2) * t181 + t266 * t139, t10 * t106 + t11 * t281 - t5 * t125 + t4 * t126 + t29 * t85 + t30 * t84 - t42 * t50 + t43 * t49 + t381, t5 * t43 + t30 * t11 + t4 * t42 + t29 * t10 + t45 * t284 + t100 * t268 + (-g(1) * t300 - g(2) * t267) * t254 + (g(1) * t267 - g(2) * t300) * t250, t12 * t89 + t26 * t384, t12 * t282 - t13 * t89 + t26 * t58 - t27 * t384, t12 * t200 + t135 * t89 + t174 * t26 + t195 * t384, -t13 * t200 + t135 * t282 - t174 * t27 + t195 * t58, t135 * t200 + t174 * t195 (-t247 * t9 + t251 * t8) * t174 + (-t247 * t32 + t251 * t31) * t135 + t303 * t200 + t6 * t195 - t44 * t58 + t90 * t13 - t22 * t282 + t54 * t27 - g(1) * t159 - g(2) * t161 + ((-t247 * t31 - t251 * t32) * t174 - t7 * t200) * qJD(6), -g(1) * t158 - g(2) * t160 + t90 * t12 - t7 * t195 + t20 * t200 + t22 * t89 + t54 * t26 + t44 * t384 + (-(-qJD(6) * t32 + t8) * t174 - t31 * t135 - t2 * t200) * t247 + (-(qJD(6) * t31 + t9) * t174 - t32 * t135 - t298 * t200) * t251; 0, 0, 0, -t309, t310, -t295, -qJ(2) * t295 - t273, 0, 0, 0, 0, 0, 0.2e1 * t193 * qJD(3) + t283 (t191 - t302) * qJD(3) + t304, 0, 0, 0, 0, 0, t271 - t346, -t178 ^ 2 * t252 - t330 - t344, -t106 * t387 + t201 * t49 + t275 * t50 - t281 * t386, -t100 * t193 + t5 * t201 - t275 * t4 - t29 * t387 - t30 * t386 - t378, 0, 0, 0, 0, 0, t290 + t358, -t354 + t374; 0, 0, 0, 0, 0, 0, 0, -t193 * t191, -t191 ^ 2 + t193 ^ 2 (-t191 - t302) * qJD(3) + t304, -t283, qJDD(3), t149 * qJD(3) - t211 * t193 + t260 + t272, -t211 * t191 - t280 + t373, t157 * t292 + t353 (t98 - t347) * t252 + (-t99 - t345) * t248, t178 * t292 + t330 - t344, t271 + t346, -t178 * t193, -pkin(3) * t99 - t130 * t178 - t149 * t155 - t80 * t193 + (t178 * t379 + t264) * t248 + t391 * t252, -pkin(3) * t98 - t149 * t157 + t321 * t178 + t81 * t193 - t248 * t391 + t264 * t252, t106 * t352 - t153 * t50 + t154 * t49 - t4 * t201 - t275 * t5 + t281 * t351 + t29 * t386 - t30 * t387 - t373, t5 * t154 + t4 * t153 - t45 * t229 - g(3) * t277 + t351 * t30 + t352 * t29 + t377 * t100 + t289 * (t229 * t231 - t232 * t362) t12 * t147 + t361 * t384, t12 * t279 - t147 * t13 - t360 * t384 + t361 * t58, -t354 - t374, t290 - t358, -t174 * t193 (t122 * t251 - t123 * t247) * t135 + t165 * t13 - t22 * t279 - t6 * t193 - t350 * t58 + t360 * t54 + (t247 * t285 - t251 * t286) * t174 + t260 * t225 -(t122 * t247 + t123 * t251) * t135 + t165 * t12 + t22 * t147 + t7 * t193 + t350 * t384 + t361 * t54 + (t247 * t286 + t251 * t285) * t174 - t260 * t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157 * t155, -t155 ^ 2 + t157 ^ 2, t98 + t347, t345 - t99, t138, -t139 * t157 + t81 * t178 + (t294 + t365) * t248 + t287 + t376, g(1) * t182 - g(2) * t180 + t139 * t155 + t178 * t80 + t252 * t365 - t265 (t241 * t49 - t243 * t50) * pkin(4) + (-t34 + t29) * t281 + (-t30 - t33) * t106, -t29 * t33 - t30 * t34 + (-t100 * t157 + t5 * t241 + t4 * t243 + t248 * t365 + t376) * pkin(4), -t398, t397, t399, t394, t135, t269 * t135 - (t23 * t251 - t24 * t247) * t174 + t91 * t58 + (-t174 * t270 - t7) * qJD(6) + t389, -t270 * t135 - t251 * t3 + (t23 * t247 + t24 * t251) * t174 - t91 * t384 + (-t174 * t269 - t356) * qJD(6) + t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106 ^ 2 - t281 ^ 2, -t29 * t106 - t281 * t30 - t260 + t45, 0, 0, 0, 0, 0, t13 + t355, t12 + t359; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t398, t397, t399, t394, t135 (-qJD(6) + t174) * t7 + t389, t6 * t174 - t251 * t298 + t396;];
tau_reg  = t1;

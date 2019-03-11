% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:23
% EndTime: 2019-03-09 05:07:31
% DurationCPUTime: 4.51s
% Computational Cost: add. (4012->488), mult. (8139->616), div. (0->0), fcn. (5461->12), ass. (0->248)
t182 = sin(qJ(6));
t186 = cos(qJ(6));
t177 = qJ(1) + pkin(10);
t167 = sin(t177);
t168 = cos(t177);
t187 = cos(qJ(4));
t183 = sin(qJ(4));
t188 = cos(qJ(3));
t299 = t183 * t188;
t92 = -t167 * t187 + t168 * t299;
t296 = t187 * t188;
t93 = t167 * t183 + t168 * t296;
t227 = t182 * t93 - t186 * t92;
t276 = t187 * qJD(3);
t184 = sin(qJ(3));
t289 = qJD(1) * t184;
t120 = t183 * t289 - t276;
t287 = qJD(3) * t183;
t122 = t187 * t289 + t287;
t180 = sin(pkin(10));
t161 = pkin(1) * t180 + pkin(7);
t143 = t161 * qJD(1);
t97 = t188 * qJD(2) - t184 * t143;
t241 = qJD(3) * pkin(3) + t97;
t206 = qJ(5) * t122 + t241;
t336 = pkin(4) + pkin(5);
t24 = -t336 * t120 + t206;
t172 = t188 * qJDD(1);
t274 = qJD(1) * qJD(3);
t116 = t184 * t274 + qJDD(4) - t172;
t281 = qJD(4) * t187;
t283 = qJD(4) * t183;
t137 = t161 * qJDD(1);
t347 = qJD(3) * t97;
t48 = qJDD(3) * pkin(8) + t184 * qJDD(2) + t188 * t137 + t347;
t218 = pkin(3) * t188 + t184 * pkin(8) + pkin(2);
t181 = cos(pkin(10));
t333 = pkin(1) * t181;
t111 = -t218 - t333;
t239 = pkin(3) * t184 - pkin(8) * t188;
t131 = t239 * qJD(3);
t58 = qJD(1) * t131 + qJDD(1) * t111;
t98 = t184 * qJD(2) + t188 * t143;
t85 = qJD(3) * pkin(8) + t98;
t86 = t111 * qJD(1);
t253 = t183 * t48 - t187 * t58 + t85 * t281 + t86 * t283;
t233 = qJDD(5) + t253;
t258 = t188 * t274;
t272 = t184 * qJDD(1);
t282 = qJD(4) * t184;
t354 = qJD(1) * t282 - qJDD(3);
t56 = -qJD(4) * t276 + (-t258 - t272) * t187 + t354 * t183;
t3 = pkin(9) * t56 - t336 * t116 + t233;
t288 = qJD(1) * t188;
t57 = t183 * (qJD(3) * (qJD(4) + t288) + t272) + t354 * t187;
t109 = t116 * qJ(5);
t356 = qJD(4) - t288;
t147 = t356 * qJD(5);
t270 = t183 * t58 + t187 * t48 + t86 * t281;
t211 = -t283 * t85 + t270;
t9 = t109 + t147 + t211;
t5 = pkin(9) * t57 + t9;
t266 = -t182 * t5 + t186 * t3;
t90 = t167 * t299 + t168 * t187;
t91 = t167 * t296 - t168 * t183;
t348 = -t182 * t91 + t90 * t186;
t67 = t120 * t182 + t122 * t186;
t298 = t184 * t187;
t300 = t183 * t186;
t99 = t182 * t298 - t184 * t300;
t359 = -g(1) * t227 + g(2) * t348 - g(3) * t99 + t24 * t67 - t266;
t358 = -qJD(2) * qJD(3) - t137;
t326 = g(3) * t188;
t238 = g(1) * t168 + g(2) * t167;
t341 = t184 * t238;
t357 = t326 - t341;
t355 = qJD(6) - qJD(4);
t221 = -t186 * t120 + t122 * t182;
t353 = -t221 ^ 2 + t67 ^ 2;
t350 = t67 * t221;
t349 = t183 * qJD(5) + t98;
t36 = -t183 * t85 + t187 * t86;
t294 = qJD(5) - t36;
t259 = t184 * t281;
t285 = qJD(3) * t188;
t346 = t183 * t285 + t259;
t260 = t183 * t282;
t263 = t188 * t276;
t204 = -t260 + t263;
t275 = -qJD(6) + t356;
t345 = -t275 - qJD(6);
t278 = qJD(6) * t186;
t279 = qJD(6) * t182;
t13 = t120 * t278 - t122 * t279 + t182 * t57 - t186 * t56;
t343 = t221 * t275 - t13;
t330 = pkin(8) * t116;
t35 = pkin(4) * t120 - t206;
t342 = -t35 * t356 + t330;
t302 = t183 * qJ(5);
t340 = -t336 * t187 - t302;
t124 = t182 * t183 + t186 * t187;
t100 = t124 * t184;
t228 = t182 * t90 + t186 * t91;
t39 = t182 * t92 + t186 * t93;
t338 = -g(1) * t39 - g(2) * t228 - g(3) * t100 - t24 * t221;
t337 = t122 ^ 2;
t335 = pkin(8) - pkin(9);
t332 = pkin(4) * t116;
t331 = pkin(4) * t183;
t325 = -t120 * t263 - t57 * t298;
t110 = -qJDD(6) + t116;
t27 = qJD(6) * t100 + t182 * t204 - t186 * t346;
t324 = t99 * t110 + t27 * t275;
t207 = t124 * t188;
t323 = -qJD(1) * t207 - t124 * t355;
t265 = t183 * t288;
t303 = t182 * t187;
t322 = t182 * t281 + t183 * t278 - t187 * t279 - t288 * t303 + (t265 - t283) * t186;
t37 = t183 * t86 + t187 * t85;
t267 = t336 * t183;
t310 = qJ(5) * t187;
t212 = -t267 + t310;
t321 = t212 * t356 + t349;
t320 = pkin(8) * qJD(4);
t149 = t356 * qJ(5);
t30 = t149 + t37;
t319 = t356 * t30;
t318 = t356 * t37;
t26 = pkin(9) * t120 + t37;
t20 = t149 + t26;
t317 = t182 * t20;
t315 = t188 * t56;
t231 = -t310 + t331;
t314 = t231 * t356 - t349;
t313 = t111 * t281 + t183 * t131;
t128 = t239 * qJD(1);
t312 = t183 * t128 + t187 * t97;
t297 = t187 * t116;
t94 = t184 * t297;
t311 = t263 * t356 + t94;
t309 = t120 * qJ(5);
t308 = t120 * t356;
t307 = t122 * t120;
t306 = t122 * t356;
t305 = t356 * t187;
t304 = t161 * t183;
t301 = t183 * t184;
t295 = pkin(9) * t122 - t294;
t293 = qJDD(2) - g(3);
t127 = t161 * t296;
t292 = t183 * t111 + t127;
t291 = t238 * t298;
t178 = t184 ^ 2;
t290 = -t188 ^ 2 + t178;
t162 = -pkin(2) - t333;
t144 = qJD(1) * t162;
t286 = qJD(3) * t184;
t284 = qJD(4) * t120;
t280 = qJD(5) * t187;
t18 = -t336 * t356 - t295;
t271 = t18 * t278 + t182 * t3 + t186 * t5;
t269 = pkin(9) * t296;
t268 = qJ(5) * t286 + t313;
t44 = qJ(5) * t289 + t312;
t146 = t335 * t187;
t262 = t122 * t285;
t261 = t356 * t283;
t256 = -pkin(4) - t304;
t255 = -t182 * t56 - t186 * t57;
t252 = t122 * t286 + t315;
t82 = t183 * t97;
t251 = -t187 * t128 + t82;
t250 = -t56 + t284;
t126 = t161 * t299;
t249 = t111 * t187 - t126;
t248 = t275 ^ 2;
t246 = t188 * qJDD(2) - t143 * t285 + t184 * t358;
t245 = t122 * t259;
t244 = -g(1) * t90 + g(2) * t92;
t243 = g(1) * t91 - g(2) * t93;
t242 = t13 * t188 - t286 * t67;
t240 = -qJD(4) * t127 - t111 * t283 + t187 * t131;
t237 = g(1) * t167 - g(2) * t168;
t185 = sin(qJ(1));
t189 = cos(qJ(1));
t236 = g(1) * t185 - g(2) * t189;
t60 = -qJ(5) * t188 + t292;
t145 = t335 * t183;
t235 = pkin(9) * t265 - qJD(6) * t145 + t335 * t283 + t44;
t234 = (-t336 * t184 - t269) * qJD(1) + t251 + t355 * t146;
t232 = pkin(4) * t187 + t302;
t8 = t182 * t18 + t186 * t20;
t176 = t188 * pkin(4);
t40 = t188 * pkin(5) + t126 + t176 + (-pkin(9) * t184 - t111) * t187;
t50 = pkin(9) * t301 + t60;
t230 = -t182 * t50 + t186 * t40;
t229 = t182 * t40 + t186 * t50;
t29 = -pkin(4) * t356 + t294;
t226 = -t183 * t30 + t187 * t29;
t225 = t183 * t29 + t187 * t30;
t220 = -t300 + t303;
t28 = -t184 * t220 * t355 + qJD(3) * t207;
t224 = -t100 * t110 - t275 * t28;
t223 = qJ(5) * t186 - t182 * t336;
t222 = qJ(5) * t182 + t186 * t336;
t217 = pkin(3) + t232;
t216 = t320 * t356 + t326;
t214 = t20 * t279 - t271;
t14 = qJD(6) * t67 + t255;
t213 = t14 * t188 - t221 * t286;
t210 = t183 * t116 + t281 * t356;
t208 = qJDD(3) * pkin(3) + t246;
t198 = -qJ(5) * t56 + qJD(5) * t122 + t208;
t11 = pkin(4) * t57 - t198;
t209 = -t11 - t216;
t205 = -qJD(1) * t144 + t238;
t203 = -t241 * t356 - t330;
t202 = g(1) * t92 + g(2) * t90 + g(3) * t301 - t253;
t201 = 0.2e1 * qJD(3) * t144 - qJDD(3) * t161;
t200 = -g(3) * t184 - t188 * t238;
t191 = qJD(3) ^ 2;
t199 = -0.2e1 * qJDD(1) * t162 - t161 * t191 + t237;
t10 = t233 - t332;
t197 = qJD(4) * t226 + t10 * t183 + t9 * t187;
t196 = t122 * t35 + qJDD(5) - t202;
t195 = g(1) * t93 + g(2) * t91 + g(3) * t298 + t356 * t36 - t211;
t194 = -t116 * t301 + t120 * t286 - t188 * t57 - t346 * t356;
t192 = qJD(1) ^ 2;
t157 = qJ(5) * t298;
t136 = qJDD(3) * t188 - t184 * t191;
t135 = qJDD(3) * t184 + t188 * t191;
t112 = pkin(3) - t340;
t77 = -t157 + (t161 + t331) * t184;
t71 = pkin(4) * t122 + t309;
t68 = t157 + (-t161 - t267) * t184;
t61 = t176 - t249;
t46 = -pkin(4) * t289 + t251;
t43 = -t336 * t122 - t309;
t33 = (qJD(4) * t232 - t280) * t184 + (t161 + t231) * t285;
t31 = -t56 + t308;
t23 = t256 * t286 - t240;
t22 = (qJD(4) * t340 + t280) * t184 + (-t161 + t212) * t285;
t19 = -qJD(5) * t188 + (-t184 * t276 - t188 * t283) * t161 + t268;
t16 = (pkin(9) * qJD(4) - qJD(3) * t161) * t298 + (-qJD(5) + (pkin(9) * qJD(3) - qJD(4) * t161) * t183) * t188 + t268;
t15 = pkin(9) * t260 + (-t269 + (-pkin(5) + t256) * t184) * qJD(3) - t240;
t7 = t18 * t186 - t317;
t6 = -t336 * t57 + t198;
t1 = [qJDD(1), t236, g(1) * t189 + g(2) * t185 (t236 + (t180 ^ 2 + t181 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t178 + 0.2e1 * t184 * t258, 0.2e1 * t172 * t184 - 0.2e1 * t274 * t290, t135, t136, 0, t184 * t201 + t188 * t199, -t184 * t199 + t188 * t201, t122 * t204 - t298 * t56, -t245 + (-t262 + (t56 + t284) * t184) * t183 + t325, -t260 * t356 + t252 + t311 (-t287 * t356 + t57) * t188 + (-qJD(3) * t120 - t210) * t184, -t116 * t188 + t286 * t356, t240 * t356 + t249 * t116 + ((t120 * t161 - t183 * t241) * qJD(3) + t253) * t188 + (-t241 * t281 + t161 * t57 - t208 * t183 + (t304 * t356 + t36) * qJD(3)) * t184 + t243, -t313 * t356 - t292 * t116 + ((t161 * t356 - t85) * t283 + (t122 * t161 - t187 * t241) * qJD(3) + t270) * t188 + (t241 * t283 - t161 * t56 - t208 * t187 + (t161 * t305 - t37) * qJD(3)) * t184 + t244, -t116 * t61 + t120 * t33 - t356 * t23 + t57 * t77 + (t287 * t35 + t10) * t188 + (-qJD(3) * t29 + t11 * t183 + t281 * t35) * t184 + t243, -t120 * t19 + t122 * t23 - t56 * t61 - t57 * t60 + t226 * t285 + (-qJD(4) * t225 + t10 * t187 - t183 * t9 + t237) * t184, t116 * t60 - t122 * t33 + t356 * t19 + t56 * t77 + (-t276 * t35 - t9) * t188 + (qJD(3) * t30 - t11 * t187 + t283 * t35) * t184 - t244, t9 * t60 + t30 * t19 + t11 * t77 + t35 * t33 + t10 * t61 + t29 * t23 - g(1) * (-pkin(1) * t185 - pkin(4) * t91 - qJ(5) * t90) - g(2) * (pkin(1) * t189 + pkin(4) * t93 + qJ(5) * t92) + (-g(1) * pkin(7) - g(2) * t218) * t168 + (-g(2) * pkin(7) + g(1) * t218) * t167, t100 * t13 + t28 * t67, -t100 * t14 - t13 * t99 - t221 * t28 - t27 * t67, t224 + t242, -t213 + t324, -t110 * t188 + t275 * t286 -(t186 * t15 - t182 * t16) * t275 - t230 * t110 + t266 * t188 - t7 * t286 + t22 * t221 + t68 * t14 + t6 * t99 + t24 * t27 + g(1) * t228 - g(2) * t39 + (-t188 * t8 + t229 * t275) * qJD(6) (qJD(6) * t230 + t182 * t15 + t186 * t16) * t275 + t229 * t110 + t214 * t188 + t8 * t286 + t22 * t67 + t68 * t13 + t6 * t100 + t24 * t28 + g(1) * t348 + g(2) * t227; 0, 0, 0, t293, 0, 0, 0, 0, 0, t136, -t135, 0, 0, 0, 0, 0, t194, -t204 * t356 + t252 - t94, t194, t245 + (t184 * t250 + t262) * t183 + t325, -t315 + (-qJD(3) * t122 - t261) * t184 + t311, -g(3) + (qJD(3) * t225 - t11) * t188 + (qJD(3) * t35 + t197) * t184, 0, 0, 0, 0, 0, t213 + t324, -t224 + t242; 0, 0, 0, 0, -t184 * t192 * t188, t290 * t192, t272, t172, qJDD(3), qJD(3) * t98 + t184 * t205 + t246 - t326, t347 + (qJD(3) * t143 - t293) * t184 + (t205 + t358) * t188, t122 * t305 - t56 * t183 (-t56 - t308) * t187 + (-t306 - t57) * t183 (-t122 * t184 - t296 * t356) * qJD(1) + t210, -t261 + t297 + (t120 * t184 + t299 * t356) * qJD(1), -t356 * t289, -t36 * t289 - pkin(3) * t57 - t98 * t120 + t82 * t356 + (-t326 + t208 - (t128 + t320) * t356) * t187 + t203 * t183 + t291, pkin(3) * t56 + t312 * t356 + t37 * t289 - t98 * t122 + t203 * t187 + (t216 - t208 - t341) * t183, t314 * t120 - t183 * t342 + t209 * t187 - t217 * t57 + t29 * t289 + t356 * t46 + t291, t120 * t44 - t122 * t46 + (t9 + t356 * t29 + (qJD(4) * t122 - t57) * pkin(8)) * t187 + (pkin(8) * t250 + t10 - t319) * t183 + t200, -t30 * t289 - t217 * t56 - t356 * t44 - t314 * t122 + t342 * t187 + (t209 + t341) * t183, -t29 * t46 - t30 * t44 + t314 * t35 + (t197 + t200) * pkin(8) + (-t11 - t357) * t217, -t13 * t220 + t323 * t67, -t124 * t13 + t14 * t220 - t221 * t323 - t322 * t67, t110 * t220 - t275 * t323 + t67 * t289, t124 * t110 - t221 * t289 + t275 * t322, -t275 * t289 -(t145 * t186 - t146 * t182) * t110 + t112 * t14 + t321 * t221 + t322 * t24 - g(3) * t207 - (t182 * t235 - t186 * t234) * t275 + t7 * t289 + (t6 + t341) * t124 (t145 * t182 + t146 * t186) * t110 + t112 * t13 + t321 * t67 + t323 * t24 - (t182 * t234 + t186 * t235) * t275 - t8 * t289 + (-t6 + t357) * t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t307, -t120 ^ 2 + t337, t31, -t57 + t306, t116, t122 * t241 + t202 + t318, -t120 * t241 + t195, -t120 * t71 - t196 + t318 + 0.2e1 * t332, pkin(4) * t56 - qJ(5) * t57 + (t30 - t37) * t122 + (t29 - t294) * t120, -t120 * t35 + t122 * t71 + 0.2e1 * t109 + 0.2e1 * t147 - t195, t9 * qJ(5) - t10 * pkin(4) - t35 * t71 - t29 * t37 - g(1) * (-pkin(4) * t92 + qJ(5) * t93) - g(2) * (-pkin(4) * t90 + qJ(5) * t91) - g(3) * (-pkin(4) * t301 + t157) + t294 * t30, -t350, -t353, t343, t275 * t67 + t14, t110, t222 * t110 - t43 * t221 - (t182 * t295 - t186 * t26) * t275 + (t223 * t275 + t8) * qJD(6) + t359, t223 * t110 - t43 * t67 - (t182 * t26 + t186 * t295) * t275 + (-t222 * t275 - t317) * qJD(6) + t271 + t338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 + t307, t31, -t356 ^ 2 - t337, t196 - t319 - t332, 0, 0, 0, 0, 0, -t186 * t110 - t122 * t221 - t182 * t248, t182 * t110 - t122 * t67 - t186 * t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t350, t353, -t343, t345 * t67 - t255, -t110, t345 * t8 - t359, -t275 * t7 + t214 - t338;];
tau_reg  = t1;

% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:09:24
% EndTime: 2019-03-09 04:09:36
% DurationCPUTime: 5.32s
% Computational Cost: add. (5191->490), mult. (10717->660), div. (0->0), fcn. (7909->14), ass. (0->238)
t207 = sin(qJ(3));
t279 = t207 * qJD(1);
t185 = qJD(5) + t279;
t178 = qJD(6) + t185;
t205 = sin(qJ(6));
t209 = cos(qJ(6));
t203 = sin(pkin(10));
t211 = cos(qJ(3));
t291 = qJD(1) * t211;
t264 = t203 * t291;
t204 = cos(pkin(10));
t281 = t204 * qJD(3);
t154 = t264 - t281;
t263 = t204 * t291;
t282 = t203 * qJD(3);
t156 = t263 + t282;
t206 = sin(qJ(5));
t210 = cos(qJ(5));
t86 = t206 * t154 - t210 * t156;
t87 = t210 * t154 + t206 * t156;
t338 = t205 * t86 - t209 * t87;
t350 = t178 * t338;
t242 = pkin(3) * t211 + qJ(4) * t207;
t165 = t242 * qJD(1);
t213 = -pkin(1) - pkin(7);
t179 = t213 * qJD(1) + qJD(2);
t305 = t211 * t179;
t101 = t203 * t165 + t204 * t305;
t267 = t203 * t279;
t76 = pkin(8) * t267 + t101;
t349 = -t204 * qJD(4) + t76;
t334 = -t206 * t203 + t210 * t204;
t348 = qJD(5) * t334;
t237 = t205 * t87 + t209 * t86;
t347 = t237 * t338;
t346 = t86 * t185;
t345 = t87 * t185;
t344 = t178 * t237;
t337 = t334 * t207;
t299 = qJD(1) * t337 + t348;
t163 = t210 * t203 + t206 * t204;
t141 = t163 * qJD(1);
t333 = t163 * qJD(5);
t298 = t207 * t141 + t333;
t343 = t237 ^ 2 - t338 ^ 2;
t283 = qJD(6) * t209;
t284 = qJD(6) * t205;
t271 = t211 * qJDD(1);
t250 = -t204 * qJDD(3) + t203 * t271;
t277 = qJD(1) * qJD(3);
t261 = t207 * t277;
t109 = t203 * t261 - t250;
t295 = t203 * qJDD(3) + t204 * t271;
t110 = t204 * t261 - t295;
t286 = qJD(5) * t210;
t287 = qJD(5) * t206;
t31 = t206 * t109 - t210 * t110 - t154 * t286 - t156 * t287;
t32 = -qJD(5) * t86 - t210 * t109 - t206 * t110;
t6 = -t205 * t32 + t209 * t31 - t87 * t283 + t284 * t86;
t342 = t6 - t350;
t241 = t207 * pkin(3) - t211 * qJ(4);
t167 = qJ(2) + t241;
t146 = t167 * qJD(1);
t166 = t207 * t179;
t147 = qJD(3) * qJ(4) + t166;
t77 = t204 * t146 - t203 * t147;
t48 = pkin(4) * t279 - t156 * pkin(8) + t77;
t78 = t203 * t146 + t204 * t147;
t50 = -pkin(8) * t154 + t78;
t18 = t206 * t48 + t210 * t50;
t13 = -t87 * pkin(9) + t18;
t11 = t13 * t284;
t200 = pkin(10) + qJ(5);
t195 = qJ(6) + t200;
t187 = sin(t195);
t212 = cos(qJ(1));
t304 = t212 * t187;
t188 = cos(t195);
t208 = sin(qJ(1));
t310 = t208 * t188;
t112 = t207 * t310 + t304;
t303 = t212 * t188;
t311 = t208 * t187;
t114 = t207 * t303 - t311;
t328 = g(3) * t211;
t138 = -qJD(3) * pkin(3) + qJD(4) - t305;
t98 = t154 * pkin(4) + t138;
t42 = t87 * pkin(5) + t98;
t341 = g(1) * t112 - g(2) * t114 + t188 * t328 - t338 * t42 + t11;
t111 = -t207 * t311 + t303;
t113 = t207 * t304 + t310;
t260 = t211 * t277;
t272 = t207 * qJDD(1);
t226 = t260 + t272;
t161 = qJDD(5) + t226;
t134 = qJD(3) * t242 - t211 * qJD(4) + qJD(2);
t74 = qJD(1) * t134 + qJDD(1) * t167;
t177 = t213 * qJDD(1) + qJDD(2);
t97 = qJDD(3) * qJ(4) + t207 * t177 + (qJD(4) + t305) * qJD(3);
t33 = -t203 * t97 + t204 * t74;
t22 = pkin(4) * t226 + t110 * pkin(8) + t33;
t34 = t203 * t74 + t204 * t97;
t25 = pkin(8) * t109 + t34;
t256 = -t206 * t25 + t210 * t22;
t219 = -qJD(5) * t18 + t256;
t2 = t161 * pkin(5) - t31 * pkin(9) + t219;
t231 = t206 * t22 + t210 * t25 + t48 * t286 - t287 * t50;
t3 = -t32 * pkin(9) + t231;
t268 = t209 * t2 - t205 * t3;
t17 = -t206 * t50 + t210 * t48;
t12 = pkin(9) * t86 + t17;
t10 = t185 * pkin(5) + t12;
t321 = t209 * t13;
t5 = t205 * t10 + t321;
t340 = -g(1) * t111 - g(2) * t113 - qJD(6) * t5 + t187 * t328 + t42 * t237 + t268;
t218 = qJD(6) * t237 - t205 * t31 - t209 * t32;
t339 = t218 - t344;
t335 = g(1) * t208 - g(2) * t212;
t336 = t335 * t204;
t326 = pkin(8) + qJ(4);
t173 = t326 * t203;
t174 = t326 * t204;
t296 = -t206 * t173 + t210 * t174;
t329 = g(3) * t207;
t223 = t211 * t335 - t329;
t234 = t203 * qJD(4) + qJD(5) * t174;
t100 = t204 * t165 - t203 * t305;
t331 = pkin(8) * t204;
t270 = t207 * t331;
t59 = (pkin(4) * t211 + t270) * qJD(1) + t100;
t332 = t173 * t286 + t349 * t210 + (t234 + t59) * t206;
t92 = t205 * t163 - t209 * t334;
t325 = -qJD(6) * t92 - t205 * t298 + t209 * t299;
t93 = t209 * t163 + t205 * t334;
t324 = qJD(6) * t93 + t205 * t299 + t209 * t298;
t152 = t204 * t167;
t259 = -t203 * t213 + pkin(4);
t85 = t207 * t259 - t211 * t331 + t152;
t312 = t207 * t213;
t108 = t203 * t167 + t204 * t312;
t314 = t203 * t211;
t99 = -pkin(8) * t314 + t108;
t322 = t206 * t85 + t210 * t99;
t289 = qJD(3) * t211;
t320 = -qJD(1) * t334 - qJD(5) * t337 - t163 * t289;
t133 = t334 * t211;
t319 = -qJD(3) * t133 + t207 * t333 + t141;
t318 = pkin(1) * qJDD(1);
t290 = qJD(3) * t207;
t235 = -qJDD(3) * pkin(3) + t179 * t290 + qJDD(4);
t102 = -t211 * t177 + t235;
t317 = t102 * t203;
t316 = t102 * t204;
t315 = t102 * t211;
t192 = sin(t200);
t309 = t208 * t192;
t193 = cos(t200);
t308 = t208 * t193;
t302 = t212 * t192;
t301 = t212 * t193;
t215 = qJD(1) ^ 2;
t300 = t215 * qJ(2);
t288 = qJD(3) * t213;
t265 = t211 * t288;
t96 = t203 * t134 + t204 * t265;
t294 = t212 * pkin(1) + t208 * qJ(2);
t201 = t207 ^ 2;
t202 = t211 ^ 2;
t293 = t201 - t202;
t214 = qJD(3) ^ 2;
t292 = -t214 - t215;
t278 = -qJD(4) + t138;
t276 = qJDD(1) * qJ(2);
t275 = qJDD(1) * t203;
t274 = qJDD(1) * t204;
t273 = qJDD(3) * t207;
t269 = 0.2e1 * qJD(1) * qJD(2);
t189 = -t204 * pkin(4) - pkin(3);
t266 = t207 * t282;
t121 = -pkin(4) * t267 + t166;
t262 = pkin(5) * t298 - t121;
t257 = qJD(6) * t10 + t3;
t116 = t204 * t134;
t57 = t116 + (t211 * t259 + t270) * qJD(3);
t73 = pkin(8) * t266 + t96;
t254 = -t206 * t73 + t210 * t57;
t253 = -t206 * t99 + t210 * t85;
t251 = -t210 * t173 - t206 * t174;
t153 = pkin(4) * t314 - t211 * t213;
t56 = t210 * t59;
t68 = pkin(9) * t334 + t296;
t249 = pkin(5) * t291 + pkin(9) * t299 + t163 * qJD(4) + t296 * qJD(5) + qJD(6) * t68 - t206 * t76 + t56;
t67 = -t163 * pkin(9) + t251;
t248 = pkin(9) * t298 - qJD(6) * t67 + t332;
t247 = g(1) * t212 + g(2) * t208;
t245 = qJDD(2) - t335;
t130 = t163 * t207;
t244 = qJD(6) * t130 + t319;
t243 = qJD(6) * t337 - t320;
t240 = -t33 * t203 + t34 * t204;
t239 = t203 * t78 + t204 * t77;
t238 = -t77 * t203 + t78 * t204;
t131 = t163 * t211;
t61 = t209 * t131 + t205 * t133;
t62 = -t205 * t131 + t209 * t133;
t139 = -pkin(4) * t266 + t207 * t288;
t233 = -t177 + t335;
t232 = t335 * t203;
t230 = t206 * t57 + t210 * t73 + t85 * t286 - t287 * t99;
t227 = 0.2e1 * qJ(2) * t277 + qJDD(3) * t213;
t224 = t233 + t300;
t221 = -t207 * t335 - t328;
t49 = -t109 * pkin(4) + t102;
t220 = -t247 + t269 + 0.2e1 * t276;
t217 = -t213 * t214 + t220;
t197 = t212 * qJ(2);
t194 = qJDD(3) * t211;
t148 = qJDD(6) + t161;
t125 = t207 * t301 - t309;
t124 = t207 * t302 + t308;
t123 = t207 * t308 + t302;
t122 = -t207 * t309 + t301;
t120 = -pkin(5) * t334 + t189;
t107 = -t203 * t312 + t152;
t95 = -t203 * t265 + t116;
t94 = t131 * pkin(5) + t153;
t72 = -t206 * t207 * t281 - t210 * t266 + t211 * t348;
t70 = -qJD(3) * t337 - t211 * t333;
t47 = t72 * pkin(5) + t139;
t27 = -pkin(9) * t131 + t322;
t26 = t207 * pkin(5) - t133 * pkin(9) + t253;
t16 = qJD(6) * t62 + t205 * t70 + t209 * t72;
t15 = -qJD(6) * t61 - t205 * t72 + t209 * t70;
t14 = t32 * pkin(5) + t49;
t9 = -t72 * pkin(9) + t230;
t8 = pkin(5) * t289 - t70 * pkin(9) - qJD(5) * t322 + t254;
t4 = t209 * t10 - t205 * t13;
t1 = [qJDD(1), t335, t247, t245 - 0.2e1 * t318, t220 -(qJDD(2) - t318) * pkin(1) - g(1) * (-t208 * pkin(1) + t197) - g(2) * t294 + (t269 + t276) * qJ(2), t202 * qJDD(1) - 0.2e1 * t207 * t260, -0.2e1 * t207 * t271 + 0.2e1 * t277 * t293, -t214 * t207 + t194, -t214 * t211 - t273, 0, t207 * t217 + t211 * t227, -t207 * t227 + t211 * t217, t232 + (t317 + t213 * t109 + (qJD(1) * t107 + t77) * qJD(3)) * t211 + (t95 * qJD(1) + t107 * qJDD(1) + t33 - t247 * t204 + (-t138 * t203 + t154 * t213) * qJD(3)) * t207, t336 + (t316 + t213 * t110 + (-qJD(1) * t108 - t78) * qJD(3)) * t211 + (-t96 * qJD(1) - t108 * qJDD(1) - t34 + t247 * t203 + (-t138 * t204 + t156 * t213) * qJD(3)) * t207, t107 * t110 + t108 * t109 - t96 * t154 - t95 * t156 + t239 * t290 + (-t203 * t34 - t204 * t33 + t247) * t211, t34 * t108 + t78 * t96 + t33 * t107 + t77 * t95 - g(1) * (t212 * t241 + t197) - g(2) * (t212 * pkin(7) + t294) + (t138 * t290 - t315) * t213 + (-g(1) * t213 - g(2) * t241) * t208, t133 * t31 - t70 * t86, -t131 * t31 - t133 * t32 - t70 * t87 + t72 * t86, t133 * t161 + t70 * t185 + t31 * t207 - t289 * t86, -t131 * t161 - t72 * t185 - t32 * t207 - t289 * t87, t161 * t207 + t185 * t289, t254 * t185 + t253 * t161 + t256 * t207 + t17 * t289 + t139 * t87 + t153 * t32 + t49 * t131 + t98 * t72 - g(1) * t125 - g(2) * t123 + (-t18 * t207 - t185 * t322) * qJD(5), g(1) * t124 - g(2) * t122 + t49 * t133 - t139 * t86 + t153 * t31 - t322 * t161 - t18 * t289 - t230 * t185 - t231 * t207 + t98 * t70, -t15 * t237 + t6 * t62, t15 * t338 + t16 * t237 + t218 * t62 - t6 * t61, t62 * t148 + t15 * t178 + t6 * t207 - t237 * t289, -t61 * t148 - t16 * t178 + t207 * t218 + t289 * t338, t148 * t207 + t178 * t289 (-t205 * t9 + t209 * t8) * t178 + (-t205 * t27 + t209 * t26) * t148 + t268 * t207 + t4 * t289 - t47 * t338 - t94 * t218 + t14 * t61 + t42 * t16 - g(1) * t114 - g(2) * t112 + ((-t205 * t26 - t209 * t27) * t178 - t5 * t207) * qJD(6), -t5 * t289 + g(1) * t113 - g(2) * t111 + t11 * t207 + t14 * t62 + t42 * t15 - t47 * t237 + t94 * t6 + (-(-qJD(6) * t27 + t8) * t178 - t26 * t148 - t2 * t207) * t205 + (-(qJD(6) * t26 + t9) * t178 - t27 * t148 - t257 * t207) * t209; 0, 0, 0, qJDD(1), -t215, t245 - t300 - t318, 0, 0, 0, 0, 0, t207 * t292 + t194, t211 * t292 - t273, -t201 * t275 + t211 * t109 + (-t204 * t215 + (t154 - 0.2e1 * t264) * qJD(3)) * t207, -t201 * t274 + t211 * t110 + (t203 * t215 + (t156 - 0.2e1 * t263) * qJD(3)) * t207 (qJD(1) * t156 + t109 * t207 - t154 * t289) * t204 + (qJD(1) * t154 - t110 * t207 + t156 * t289) * t203, -t315 + t240 * t207 - t239 * qJD(1) + (t138 * t207 + t211 * t238) * qJD(3) - t335, 0, 0, 0, 0, 0, -t130 * t161 + t320 * t185 - t211 * t32 + t87 * t290, -t161 * t337 + t319 * t185 - t211 * t31 - t290 * t86, 0, 0, 0, 0, 0 (-t209 * t130 - t205 * t337) * t148 - t338 * t290 + t211 * t218 + (t205 * t244 - t209 * t243) * t178 -(-t205 * t130 + t209 * t337) * t148 - t237 * t290 - t211 * t6 + (t205 * t243 + t209 * t244) * t178; 0, 0, 0, 0, 0, 0, t211 * t215 * t207, -t293 * t215, t271, -t272, qJDD(3), -t211 * t224 + t329, t207 * t224 + t328, pkin(3) * t109 - t316 + (-t336 + (-qJ(4) * t282 - t77) * qJD(1)) * t211 + (-qJ(4) * t275 + g(3) * t204 - t154 * t179 + (t203 * t278 - t100) * qJD(1)) * t207, pkin(3) * t110 + t317 + (t232 + (-qJ(4) * t281 + t78) * qJD(1)) * t211 + (-qJ(4) * t274 - g(3) * t203 - t156 * t179 + (t204 * t278 + t101) * qJD(1)) * t207, t100 * t156 + t101 * t154 + (qJ(4) * t109 - qJD(4) * t154 - t279 * t77 + t34) * t204 + (-qJ(4) * t110 + qJD(4) * t156 - t279 * t78 - t33) * t203 + t221, -t138 * t166 - t77 * t100 - t78 * t101 + t238 * qJD(4) + (-t102 - t223) * pkin(3) + (t221 + t240) * qJ(4), t163 * t31 - t299 * t86, -t163 * t32 + t298 * t86 - t299 * t87 + t31 * t334, t163 * t161 + t185 * t299 + t291 * t86, t161 * t334 - t185 * t298 + t291 * t87, -t185 * t291, t251 * t161 + t189 * t32 - t49 * t334 - t17 * t291 - t121 * t87 + t298 * t98 + (-t56 - t234 * t210 + (qJD(5) * t173 + t349) * t206) * t185 - t223 * t193, t121 * t86 - t296 * t161 + t49 * t163 + t18 * t291 + t185 * t332 + t189 * t31 + t223 * t192 + t299 * t98, -t237 * t325 + t6 * t93, t218 * t93 + t237 * t324 + t325 * t338 - t6 * t92, t93 * t148 + t178 * t325 + t237 * t291, -t92 * t148 - t178 * t324 - t291 * t338, -t178 * t291 (-t205 * t68 + t209 * t67) * t148 - t120 * t218 + t14 * t92 - t4 * t291 + t324 * t42 - t262 * t338 + (t205 * t248 - t209 * t249) * t178 - t223 * t188 -(t205 * t67 + t209 * t68) * t148 + t120 * t6 + t14 * t93 + t5 * t291 + t325 * t42 - t262 * t237 + (t205 * t249 + t209 * t248) * t178 + t223 * t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t156 - t282) * t279 + t250 (-t154 - t281) * t279 + t295, -t154 ^ 2 - t156 ^ 2, t78 * t154 + t77 * t156 + t211 * t233 + t235 - t329, 0, 0, 0, 0, 0, t32 - t346, t31 - t345, 0, 0, 0, 0, 0, -t218 - t344, t6 + t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86 * t87, t86 ^ 2 - t87 ^ 2, t31 + t345, -t32 - t346, t161, -g(1) * t122 - g(2) * t124 + t18 * t185 + t192 * t328 + t86 * t98 + t219, g(1) * t123 - g(2) * t125 + t17 * t185 + t193 * t328 + t98 * t87 - t231, t347, t343, t342, t339, t148 -(-t205 * t12 - t321) * t178 + (t209 * t148 - t178 * t284 - t338 * t86) * pkin(5) + t340 (-t13 * t178 - t2) * t205 + (t12 * t178 - t257) * t209 + (-t205 * t148 - t178 * t283 - t237 * t86) * pkin(5) + t341; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t347, t343, t342, t339, t148, t5 * t178 + t340, t4 * t178 - t205 * t2 - t209 * t257 + t341;];
tau_reg  = t1;

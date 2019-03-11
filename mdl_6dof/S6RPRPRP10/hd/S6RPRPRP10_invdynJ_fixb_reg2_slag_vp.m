% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:40
% EndTime: 2019-03-09 03:32:47
% DurationCPUTime: 4.20s
% Computational Cost: add. (4168->492), mult. (7533->548), div. (0->0), fcn. (4128->6), ass. (0->259)
t170 = cos(qJ(5));
t269 = qJD(5) * t170;
t167 = sin(qJ(5));
t268 = t170 * qJD(3);
t168 = sin(qJ(3));
t274 = qJD(1) * t168;
t94 = t167 * t274 + t268;
t299 = qJD(5) * t94;
t171 = cos(qJ(3));
t265 = qJD(1) * qJD(3);
t245 = t171 * t265;
t262 = t168 * qJDD(1);
t352 = t245 + t262;
t39 = t167 * qJDD(3) - t170 * t352 + t299;
t35 = t39 * t167;
t248 = t170 * t274;
t273 = qJD(3) * t167;
t92 = -t248 + t273;
t353 = -t92 * t269 - t35;
t172 = cos(qJ(1));
t156 = g(2) * t172;
t169 = sin(qJ(1));
t337 = g(1) * t169 - t156;
t267 = t171 * qJD(1);
t281 = pkin(3) * t274 + qJD(1) * qJ(2);
t75 = -qJ(4) * t267 + t281;
t351 = -t75 * qJD(1) - t337;
t176 = qJD(1) ^ 2;
t289 = t176 * qJ(2);
t350 = t337 + t289;
t160 = qJDD(1) * qJ(2);
t128 = qJD(5) + t267;
t300 = qJD(3) * t94;
t146 = t171 * qJDD(1);
t246 = t168 * t265;
t336 = -t246 + t146;
t91 = -qJDD(5) - t336;
t72 = t167 * t91;
t207 = t300 - t72;
t238 = qJD(5) * t171 + qJD(1);
t272 = qJD(3) * t168;
t250 = t167 * t272;
t270 = qJD(5) * t167;
t38 = qJD(3) * t270 - qJD(5) * t248 - t170 * qJDD(3) - t167 * t352;
t312 = t168 * t38;
t349 = t128 * (t170 * t238 - t250) + t207 * t171 - t312;
t271 = qJD(3) * t171;
t310 = t170 * t38;
t348 = t168 * ((t167 * t94 + t170 * t92) * qJD(5) + t35 + t310) + (t167 * t92 - t170 * t94) * t271;
t204 = t94 * t128;
t205 = t92 * t128;
t347 = t167 * (t38 + t205) - t170 * (t39 + t204);
t174 = -pkin(1) - pkin(7);
t316 = pkin(4) - t174;
t162 = qJD(1) * qJD(2);
t208 = pkin(3) * t352 + qJ(4) * t246 + t160 + t162;
t241 = qJD(3) * pkin(8) - qJD(4);
t266 = qJ(4) * qJDD(1);
t22 = pkin(8) * t262 + (t241 * qJD(1) - t266) * t171 + t208;
t173 = -pkin(3) - pkin(8);
t116 = t174 * qJDD(1) + qJDD(2);
t237 = -t171 * t116 + qJDD(4);
t119 = t174 * qJD(1) + qJD(2);
t96 = t119 * t272;
t224 = t237 + t96;
t28 = pkin(4) * t336 + t173 * qJDD(3) + t224;
t103 = t171 * t119;
t334 = qJD(4) - t103;
t284 = pkin(4) * t267 + t334;
t51 = t173 * qJD(3) + t284;
t149 = t171 * qJ(4);
t318 = t168 * pkin(8);
t55 = (-t149 + t318) * qJD(1) + t281;
t244 = t167 * t22 - t170 * t28 + t55 * t269 + t51 * t270;
t320 = pkin(5) * t91;
t2 = qJDD(6) + t244 + t320;
t21 = t167 * t51 + t170 * t55;
t14 = t128 * qJ(6) + t21;
t313 = t14 * t128;
t346 = -t2 + t313;
t305 = t21 * t128;
t345 = -t244 + t305;
t344 = t128 * t269 - t72;
t73 = t170 * t91;
t343 = -t128 * t270 - t73;
t341 = t39 - t204;
t297 = t128 * t167;
t301 = qJD(3) * t92;
t339 = -t128 * t297 - t301 - t73;
t165 = t168 ^ 2;
t166 = t171 ^ 2;
t275 = t165 + t166;
t243 = t275 * t116;
t258 = 0.2e1 * t162;
t335 = t258 + 0.2e1 * t160;
t308 = t173 * t91;
t163 = qJD(3) * qJ(4);
t102 = t168 * t119;
t67 = -pkin(4) * t274 + t102;
t57 = t163 + t67;
t333 = t128 * t57 - t308;
t155 = g(3) * t168;
t261 = t167 * t28 + t170 * t22 + t51 * t269;
t291 = t172 * t167;
t79 = t169 * t170 + t171 * t291;
t290 = t172 * t170;
t294 = t169 * t171;
t81 = -t167 * t294 + t290;
t332 = -g(1) * t81 - g(2) * t79 - (qJD(5) * t55 + t155) * t167 + t261;
t319 = g(3) * t171;
t199 = qJD(5) * t128 * t173 + t319;
t331 = -t168 * t337 - t199;
t242 = -qJD(3) * pkin(3) + qJD(4);
t74 = -t103 + t242;
t76 = -t102 - t163;
t210 = t168 * t74 - t171 * t76;
t158 = qJDD(3) * qJ(4);
t161 = qJD(3) * qJD(4);
t99 = t168 * t116;
t43 = -t119 * t271 - t158 - t161 - t99;
t298 = qJDD(3) * pkin(3);
t44 = t224 - t298;
t183 = qJD(3) * t210 - t43 * t168 - t44 * t171;
t279 = t168 * pkin(3) - t149;
t105 = qJ(2) + t279;
t263 = qJDD(3) * t174;
t328 = (qJD(1) * t105 + t75) * qJD(3) + t263;
t107 = t316 * t171;
t232 = -t279 - t318;
t84 = qJ(2) - t232;
t314 = t167 * t107 + t170 * t84;
t247 = pkin(3) * t271 + qJ(4) * t272 + qJD(2);
t53 = t241 * t171 + t247;
t89 = t316 * t272;
t11 = -qJD(5) * t314 - t167 * t53 - t170 * t89;
t20 = -t167 * t55 + t170 * t51;
t216 = t167 * t20 - t170 * t21;
t3 = -t55 * t270 + t261;
t326 = -qJD(5) * t216 + t3 * t167 - t170 * t244;
t304 = t91 * qJ(6);
t1 = t128 * qJD(6) + t3 - t304;
t283 = qJD(6) - t20;
t13 = -t128 * pkin(5) + t283;
t218 = t13 * t167 + t14 * t170;
t325 = qJD(5) * t218 + t1 * t167 - t2 * t170;
t324 = t168 * (-t301 - t343) - t171 * (t128 * t268 - t39);
t323 = t94 ^ 2;
t322 = t128 ^ 2;
t317 = t94 * t92;
t98 = pkin(3) * t267 + qJ(4) * t274;
t66 = pkin(8) * t267 + t98;
t31 = t167 * t67 + t170 * t66;
t222 = pkin(5) * t170 + qJ(6) * t167;
t202 = -pkin(4) - t222;
t194 = t171 * t202;
t315 = -qJD(1) * t194 + qJD(5) * t222 - t170 * qJD(6) + t334;
t311 = t168 * t94;
t309 = t173 * t38;
t307 = t173 * t94;
t306 = t20 * t128;
t303 = pkin(1) * qJDD(1);
t302 = qJ(4) * t168;
t296 = t168 * t169;
t295 = t168 * t172;
t293 = t170 * t171;
t292 = t171 * t172;
t288 = t176 * t168;
t23 = t92 * pkin(5) - t94 * qJ(6) + t57;
t287 = t23 * qJD(3);
t286 = t57 * qJD(3);
t282 = pkin(3) * t294 + qJ(4) * t296;
t280 = (t258 + t160) * qJ(2);
t278 = t172 * pkin(1) + t169 * qJ(2);
t276 = t165 - t166;
t264 = qJDD(3) * t168;
t251 = t170 * t267;
t260 = t92 * t251 - t353;
t259 = t128 * t251 + t344;
t257 = t92 ^ 2 - t323;
t256 = t174 * t169;
t255 = t173 * t171;
t70 = t92 * t271;
t254 = pkin(8) * t294 + t282;
t253 = g(1) * t294 - g(2) * t292 - t155;
t252 = t172 * pkin(7) + t278;
t108 = t275 * qJDD(1);
t239 = qJDD(2) - t303;
t175 = qJD(3) ^ 2;
t110 = qJDD(3) * t171 - t175 * t168;
t236 = pkin(3) * t296 + t252;
t235 = t168 * t245;
t78 = t169 * t167 - t171 * t290;
t80 = t169 * t293 + t291;
t234 = -g(1) * t78 + g(2) * t80;
t233 = g(1) * t79 - g(2) * t81;
t231 = g(1) * t172 + g(2) * t169;
t225 = t173 * t353 - t253;
t223 = -pkin(3) * t171 - t302;
t221 = t167 * pkin(5) - t170 * qJ(6);
t219 = t13 * t170 - t14 * t167;
t217 = t167 * t21 + t170 * t20;
t30 = -t167 * t66 + t170 * t67;
t41 = t170 * t107 - t167 * t84;
t150 = t172 * qJ(2);
t201 = pkin(3) * t295 - qJ(4) * t292 + t150;
t104 = qJ(4) + t221;
t200 = t172 * pkin(4) + pkin(8) * t296 + t236;
t10 = t107 * t269 - t167 * t89 + t170 * t53 - t84 * t270;
t198 = t92 * t274 + t259;
t197 = pkin(8) * t295 + t201;
t196 = 0.2e1 * qJ(2) * t265 + t263;
t193 = -t174 * t175 - t231;
t192 = -t174 * t108 + t337;
t191 = t75 * t267 + t237 + t253;
t190 = t128 * t23 - t308;
t188 = (g(1) * t316 + g(2) * t149) * t169;
t187 = -g(1) * t296 - t199;
t185 = g(1) * t80 + g(2) * t78 - t170 * t155 - t244;
t184 = t167 * t204 - t260 + t310;
t29 = -pkin(4) * t352 - t43;
t182 = -t170 * t70 + (-t170 * t39 + t92 * t270) * t168;
t181 = t193 + t335;
t33 = (-qJD(1) * qJD(4) - t266) * t171 + t208;
t64 = -t171 * qJD(4) + t247;
t180 = -qJD(1) * t64 - qJDD(1) * t105 - t193 - t33;
t179 = t23 * t94 + qJDD(6) - t185;
t178 = t91 * t293 + t168 * t39 + t70 + (t238 * t167 + t168 * t268) * t128;
t177 = -t92 * t250 + (-qJD(1) * t94 + (t39 - t299) * t171) * t167 + (-t94 * t272 + qJD(1) * t92 + (qJD(5) * t92 - t38) * t171) * t170;
t140 = t168 * t174;
t126 = t174 * t271;
t124 = t171 * t288;
t112 = g(2) * t168 * t291;
t111 = t276 * t176;
t109 = t175 * t171 + t264;
t106 = -t168 * pkin(4) + t140;
t100 = t128 * t274;
t90 = -pkin(4) * t271 + t126;
t86 = t166 * qJDD(1) - 0.2e1 * t235;
t85 = t165 * qJDD(1) + 0.2e1 * t235;
t83 = t110 - t288;
t82 = t264 + (t175 + t176) * t171;
t54 = -0.2e1 * t168 * t146 + 0.2e1 * t276 * t265;
t52 = t168 * t202 + t140;
t48 = -t128 * t272 - t91 * t171;
t45 = pkin(5) * t94 + qJ(6) * t92;
t37 = -t171 * pkin(5) - t41;
t36 = t171 * qJ(6) + t314;
t27 = pkin(5) * t274 - t30;
t26 = -qJ(6) * t274 + t31;
t19 = t126 + (qJD(5) * t221 - qJD(6) * t167) * t168 + qJD(3) * t194;
t16 = -t38 + t205;
t15 = (-t171 * t297 + t311) * qJD(1) + t343;
t12 = -t297 * t94 - t310;
t9 = t269 * t311 + (t94 * t271 - t312) * t167;
t8 = pkin(5) * t272 - t11;
t7 = -qJ(6) * t272 + t171 * qJD(6) + t10;
t6 = (t128 * t273 - t38) * t171 + (-t300 + t344) * t168;
t5 = t39 * pkin(5) + t38 * qJ(6) - t94 * qJD(6) + t29;
t4 = [0, 0, 0, 0, 0, qJDD(1), t337, t231, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t337 - 0.2e1 * t303, -t231 + t335, -t239 * pkin(1) - g(1) * (-t169 * pkin(1) + t150) - g(2) * t278 + t280, t86, t54, t110, t85, -t109, 0, t168 * t181 + t171 * t196, -t168 * t196 + t171 * t181, t192 - t243, -g(1) * (t150 + t256) - g(2) * t252 + t174 * t243 + t280, 0, -t110, t109, t86, t54, t85, -t183 + t192, t180 * t168 - t171 * t328, t168 * t328 + t180 * t171, t33 * t105 + t75 * t64 - g(1) * (t256 + t201) - g(2) * (-t149 * t169 + t236) + t183 * t174, t9, -t348, t6, t182, -t324, t48, t106 * t39 + t11 * t128 - t41 * t91 + t90 * t92 + (-t268 * t57 - t244) * t171 + (-qJD(3) * t20 - t29 * t170 + t270 * t57) * t168 + t233, -t10 * t128 - t106 * t38 + t314 * t91 + t90 * t94 + (t273 * t57 - t3) * t171 + (qJD(3) * t21 + t29 * t167 + t269 * t57) * t168 + t234, -t10 * t92 - t11 * t94 + t41 * t38 - t314 * t39 - t216 * t271 + (-qJD(5) * t217 + t167 * t244 + t170 * t3 - t231) * t168, -g(1) * t197 - g(2) * t200 + t21 * t10 + t29 * t106 + t20 * t11 - t244 * t41 + t3 * t314 + t57 * t90 + t188, t9, t6, t348, t48, t324, t182, -t8 * t128 + t19 * t92 + t37 * t91 + t52 * t39 + (-t23 * t268 - t2) * t171 + (qJD(3) * t13 - t5 * t170 + t23 * t270) * t168 + t233, -t36 * t39 - t37 * t38 - t7 * t92 + t8 * t94 + t218 * t271 + (qJD(5) * t219 + t1 * t170 + t167 * t2 - t231) * t168, t7 * t128 - t19 * t94 - t36 * t91 + t52 * t38 + (-t23 * t273 + t1) * t171 + (-qJD(3) * t14 - t5 * t167 - t23 * t269) * t168 - t234, t1 * t36 + t14 * t7 + t5 * t52 + t23 * t19 + t2 * t37 + t13 * t8 - g(1) * (-t79 * pkin(5) - t78 * qJ(6) + t197) - g(2) * (t81 * pkin(5) + t80 * qJ(6) + t200) + t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t176, -t350 + t239, 0, 0, 0, 0, 0, 0, t83, -t82, -t108, t243 - t350, 0, 0, 0, 0, 0, 0, -t108, -t83, t82, t183 + t351, 0, 0, 0, 0, 0, 0, t178, t349, t177, t216 * qJD(1) + (qJD(3) * t217 + t29) * t168 + (t286 - t326) * t171 - t337, 0, 0, 0, 0, 0, 0, t178, t177, -t349, -t218 * qJD(1) + (-qJD(3) * t219 + t5) * t168 + (t287 - t325) * t171 - t337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, -t111, t146, -t124, -t262, qJDD(3) (t116 - t289) * t171 - t253, t168 * t350 + t319 - t99, 0, 0, qJDD(3), -t146, t262, t124, -t111, -t124, t223 * qJDD(1) + ((-t76 - t163) * t171 + (-t242 + t74) * t168) * qJD(1), t274 * t98 + t191 - 0.2e1 * t298, 0.2e1 * t158 + 0.2e1 * t161 + t99 + (qJD(1) * t98 - g(3)) * t171 + t351 * t168, -t44 * pkin(3) - g(1) * t282 + g(3) * t279 - t43 * qJ(4) - t76 * qJD(4) - t119 * t210 - t223 * t156 - t75 * t98, t12, t347, t15, t260, -t198, t100, t20 * t274 + qJ(4) * t39 - t30 * t128 + t112 + t284 * t92 + t333 * t170 + (t187 + t29) * t167, -t21 * t274 - qJ(4) * t38 + t31 * t128 + t284 * t94 - t333 * t167 + (t29 + t331) * t170, t30 * t94 + t31 * t92 + (t309 - t345) * t170 + (t20 * t267 - t3 + (t20 + t307) * qJD(5)) * t167 + t225, t29 * qJ(4) - t21 * t31 - t20 * t30 - g(1) * t254 - g(3) * t232 + t284 * t57 - (t255 - t302) * t156 + t326 * t173, t12, t15, -t347, t100, t198, t170 * t205 + t35, -t13 * t274 + t104 * t39 + t27 * t128 + t112 + t315 * t92 + t190 * t170 + (t187 + t5) * t167, t26 * t92 - t27 * t94 + (t309 - t346) * t170 + (-t13 * t267 - t1 + (-t13 + t307) * qJD(5)) * t167 + t225, t14 * t274 + t104 * t38 - t26 * t128 - t315 * t94 + t190 * t167 + (-t331 - t5) * t170, t5 * t104 - t14 * t26 - t13 * t27 - g(1) * (t221 * t296 + t254) - g(3) * (t171 * t221 + t232) + t315 * t23 + t325 * t173 - (-t104 * t168 + t255) * t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, qJDD(3) - t124, -t176 * t166 - t175, t76 * qJD(3) + t191 - t298 + t96, 0, 0, 0, 0, 0, 0, t339, -t322 * t170 - t207, t184, -t286 + t345 * t170 + (t3 - t306) * t167 + t253, 0, 0, 0, 0, 0, 0, t339, t184, t259 + t300, -t287 + t346 * t170 + (t128 * t13 + t1) * t167 + t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t317, -t257, t16, -t317, -t341, -t91, -t57 * t94 + t185 + t305, t57 * t92 + t306 - t332, 0, 0, t317, t16, t257, -t91, t341, -t317, -t45 * t92 - t179 + t305 - 0.2e1 * t320, pkin(5) * t38 - qJ(6) * t39 + (t14 - t21) * t94 + (t13 - t283) * t92, -0.2e1 * t304 - t23 * t92 + t45 * t94 + (0.2e1 * qJD(6) - t20) * t128 + t332, t1 * qJ(6) - t2 * pkin(5) - t23 * t45 - t13 * t21 - g(1) * (-pkin(5) * t80 + qJ(6) * t81) - g(2) * (-t78 * pkin(5) + t79 * qJ(6)) - t222 * t155 + t283 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91 + t317, t16, -t322 - t323, t179 - t313 + t320;];
tau_reg  = t4;

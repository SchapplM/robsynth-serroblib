% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:09:34
% EndTime: 2020-01-03 12:09:40
% DurationCPUTime: 3.32s
% Computational Cost: add. (6325->389), mult. (9824->486), div. (0->0), fcn. (6824->16), ass. (0->243)
t240 = cos(qJ(3));
t216 = t240 * qJD(4);
t236 = sin(qJ(3));
t234 = -qJ(4) - pkin(7);
t284 = qJD(3) * t234;
t152 = t236 * t284 + t216;
t153 = -t236 * qJD(4) + t240 * t284;
t232 = sin(pkin(9));
t233 = cos(pkin(9));
t163 = t232 * t240 + t233 * t236;
t241 = cos(qJ(2));
t303 = qJD(1) * t241;
t293 = pkin(1) * t303;
t323 = -t152 * t232 + t233 * t153 + t163 * t293;
t314 = t233 * t240;
t316 = t232 * t236;
t162 = -t314 + t316;
t322 = t233 * t152 + t232 * t153 + t162 * t293;
t237 = sin(qJ(2));
t302 = qJD(2) * t237;
t212 = pkin(1) * t302;
t345 = pkin(1) * t241;
t306 = -qJD(1) * t212 + qJDD(1) * t345;
t224 = qJDD(1) + qJDD(2);
t344 = pkin(2) * t224;
t139 = -t306 - t344;
t231 = qJ(1) + qJ(2);
t217 = sin(t231);
t218 = cos(t231);
t348 = g(2) * t218 + g(3) * t217;
t359 = t139 + t348;
t235 = sin(qJ(5));
t239 = cos(qJ(5));
t227 = qJD(1) + qJD(2);
t292 = t227 * t314;
t130 = t227 * t316 - t292;
t132 = t163 * t227;
t261 = t130 * t235 - t239 * t132;
t155 = t163 * qJD(3);
t312 = t240 * t224;
t313 = t236 * t224;
t263 = t232 * t313 - t233 * t312;
t79 = t155 * t227 + t263;
t300 = qJD(3) * t236;
t289 = t227 * t300;
t254 = t163 * t224 - t232 * t289;
t299 = qJD(3) * t240;
t288 = t227 * t299;
t80 = t233 * t288 + t254;
t248 = qJD(5) * t261 - t235 * t80 - t239 * t79;
t226 = qJD(3) + qJD(5);
t326 = t261 * t226;
t358 = t248 - t326;
t297 = qJD(5) * t239;
t298 = qJD(5) * t235;
t256 = -t130 * t297 - t132 * t298 - t235 * t79 + t239 * t80;
t74 = -t239 * t130 - t132 * t235;
t327 = t226 * t74;
t357 = t256 - t327;
t333 = t74 ^ 2;
t334 = t261 ^ 2;
t356 = -t333 + t334;
t332 = t74 * t261;
t156 = t162 * qJD(3);
t337 = pkin(8) * t156;
t355 = t337 + t323;
t147 = t155 * pkin(8);
t354 = -t147 + t322;
t229 = t236 ^ 2;
t230 = t240 ^ 2;
t304 = t229 + t230;
t349 = t241 * t304;
t353 = t227 * t349;
t296 = qJDD(1) * t237;
t301 = qJD(2) * t241;
t140 = t224 * pkin(7) + (qJD(1) * t301 + t296) * pkin(1);
t280 = t304 * t140;
t204 = g(3) * t218;
t347 = g(2) * t217 - t204;
t228 = qJ(3) + pkin(9);
t215 = qJ(5) + t228;
t197 = sin(t215);
t198 = cos(t215);
t257 = qJ(4) * t224 + qJD(4) * t227 + t140;
t329 = pkin(1) * qJD(1);
t295 = t237 * t329;
t173 = pkin(7) * t227 + t295;
t278 = qJ(4) * t227 + t173;
t260 = qJD(3) * t278;
t55 = qJDD(3) * pkin(3) - t236 * t257 - t240 * t260;
t58 = -t236 * t260 + t240 * t257;
t26 = -t232 * t58 + t233 * t55;
t16 = qJDD(3) * pkin(4) - pkin(8) * t80 + t26;
t27 = t232 * t55 + t233 * t58;
t17 = -pkin(8) * t79 + t27;
t338 = pkin(8) * t132;
t118 = t278 * t240;
t105 = t232 * t118;
t117 = t278 * t236;
t111 = qJD(3) * pkin(3) - t117;
t61 = t233 * t111 - t105;
t43 = qJD(3) * pkin(4) - t338 + t61;
t339 = pkin(8) * t130;
t315 = t233 * t118;
t62 = t232 * t111 + t315;
t45 = t62 - t339;
t4 = (qJD(5) * t43 + t17) * t239 + t235 * t16 - t45 * t298;
t340 = pkin(3) * t240;
t208 = pkin(2) + t340;
t129 = -t208 * t227 + qJD(4) - t293;
t81 = t130 * pkin(4) + t129;
t352 = g(1) * t197 + t347 * t198 - t74 * t81 - t4;
t318 = t197 * t218;
t319 = t197 * t217;
t19 = t235 * t43 + t239 * t45;
t5 = -qJD(5) * t19 + t239 * t16 - t235 * t17;
t351 = -g(1) * t198 + g(2) * t319 - g(3) * t318 + t261 * t81 + t5;
t343 = pkin(2) * t227;
t174 = -t293 - t343;
t350 = t173 * t349 + t174 * t237;
t346 = t132 ^ 2;
t342 = pkin(3) * t232;
t341 = pkin(3) * t236;
t336 = pkin(8) * t163;
t335 = g(1) * t240;
t186 = t234 * t236;
t219 = t240 * qJ(4);
t187 = pkin(7) * t240 + t219;
t109 = t233 * t186 - t187 * t232;
t84 = t109 - t336;
t110 = t232 * t186 + t233 * t187;
t159 = t162 * pkin(8);
t85 = -t159 + t110;
t38 = -t235 * t85 + t239 * t84;
t331 = qJD(5) * t38 + t355 * t235 + t354 * t239;
t39 = t235 * t84 + t239 * t85;
t330 = -qJD(5) * t39 - t354 * t235 + t355 * t239;
t207 = pkin(1) * t237 + pkin(7);
t311 = -qJ(4) - t207;
t275 = qJD(3) * t311;
t294 = pkin(1) * t301;
t100 = (-qJD(4) - t294) * t236 + t240 * t275;
t99 = t236 * t275 + t240 * t294 + t216;
t57 = t232 * t100 + t233 * t99;
t328 = pkin(1) * qJD(2);
t199 = pkin(3) * t233 + pkin(4);
t150 = t199 * t239 - t235 * t342;
t63 = t117 * t232 - t315;
t46 = t63 + t339;
t64 = -t233 * t117 - t105;
t47 = t64 - t338;
t325 = t150 * qJD(5) - t235 * t46 - t239 * t47;
t151 = t199 * t235 + t239 * t342;
t324 = -t151 * qJD(5) + t235 * t47 - t239 * t46;
t321 = t132 * t130;
t317 = t227 * t236;
t160 = t311 * t236;
t161 = t207 * t240 + t219;
t90 = t232 * t160 + t233 * t161;
t214 = cos(t228);
t287 = pkin(4) * t214 + t340;
t167 = pkin(2) + t287;
t225 = -pkin(8) + t234;
t310 = t217 * t167 + t218 * t225;
t309 = t217 * t208 + t218 * t234;
t308 = t218 * pkin(2) + t217 * pkin(7);
t305 = t229 - t230;
t211 = pkin(3) * t300;
t222 = t227 ^ 2;
t291 = t236 * t222 * t240;
t290 = t227 * t302;
t114 = pkin(4) * t155 + t211;
t56 = t233 * t100 - t232 * t99;
t279 = t304 * t224;
t89 = t233 * t160 - t161 * t232;
t277 = t218 * t167 - t217 * t225;
t276 = t218 * t208 - t217 * t234;
t86 = pkin(3) * t289 - t208 * t224 + qJDD(4) - t306;
t44 = t79 * pkin(4) + t86;
t48 = t235 * t155 + t239 * t156 + t162 * t297 + t163 * t298;
t94 = -t162 * t235 + t163 * t239;
t274 = g(2) * t318 + g(3) * t319 + t44 * t94 - t81 * t48;
t213 = sin(t228);
t273 = -t129 * t156 + t86 * t163 + t348 * t213;
t272 = t227 * t295;
t271 = -t347 + t280;
t270 = t174 * t299 + t359 * t236;
t269 = t306 - t348;
t268 = t236 * t288;
t267 = t114 - t295;
t238 = sin(qJ(1));
t242 = cos(qJ(1));
t264 = -g(2) * t242 - g(3) * t238;
t67 = t89 - t336;
t68 = -t159 + t90;
t34 = -t235 * t68 + t239 * t67;
t35 = t235 * t67 + t239 * t68;
t18 = -t235 * t45 + t239 * t43;
t49 = qJD(5) * t94 + t239 * t155 - t235 * t156;
t93 = t239 * t162 + t163 * t235;
t262 = t18 * t48 - t19 * t49 - t4 * t93 - t5 * t94 - t347;
t125 = pkin(4) * t162 - t208;
t259 = -t62 * t155 + t61 * t156 - t27 * t162 - t26 * t163 - t347;
t255 = -t174 * t227 - t140 + t347;
t253 = -t198 * t348 + t44 * t93 + t81 * t49;
t252 = t129 * t155 + t86 * t162 - t214 * t348;
t243 = qJD(3) ^ 2;
t251 = pkin(7) * t243 - t272 - t344;
t209 = -pkin(2) - t345;
t250 = pkin(1) * t290 + t207 * t243 + t209 * t224;
t249 = -pkin(7) * qJDD(3) + (t293 - t343) * qJD(3);
t247 = -qJDD(3) * t207 + (t209 * t227 - t294) * qJD(3);
t223 = qJDD(3) + qJDD(5);
t221 = t242 * pkin(1);
t220 = t238 * pkin(1);
t201 = t217 * pkin(2);
t184 = -t208 - t345;
t182 = qJDD(3) * t240 - t236 * t243;
t181 = qJDD(3) * t236 + t240 * t243;
t172 = t212 + t211;
t157 = t174 * t300;
t142 = t224 * t230 - 0.2e1 * t268;
t141 = t224 * t229 + 0.2e1 * t268;
t124 = t130 ^ 2;
t113 = t125 - t345;
t104 = -0.2e1 * qJD(3) * t227 * t305 + 0.2e1 * t236 * t312;
t103 = t114 + t212;
t98 = pkin(3) * t317 + pkin(4) * t132;
t92 = -qJD(3) * t156 + qJDD(3) * t163;
t91 = -qJD(3) * t155 - qJDD(3) * t162;
t41 = -t147 + t57;
t40 = t56 + t337;
t37 = -t132 * t156 + t163 * t80;
t36 = t130 * t155 + t162 * t79;
t29 = -t223 * t93 - t226 * t49;
t28 = t223 * t94 - t226 * t48;
t15 = t130 * t156 - t132 * t155 - t162 * t80 - t163 * t79;
t9 = -qJD(5) * t35 - t235 * t41 + t239 * t40;
t8 = qJD(5) * t34 + t235 * t40 + t239 * t41;
t7 = -t248 * t93 - t49 * t74;
t6 = t256 * t94 + t261 * t48;
t1 = t248 * t94 - t256 * t93 + t261 * t49 - t48 * t74;
t2 = [0, 0, 0, 0, 0, qJDD(1), t264, g(2) * t238 - g(3) * t242, 0, 0, 0, 0, 0, 0, 0, t224, (t224 * t241 - t290) * pkin(1) + t269, ((-qJDD(1) - t224) * t237 + (-qJD(1) - t227) * t301) * pkin(1) + t347, 0, (t264 + (t237 ^ 2 + t241 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t141, t104, t181, t142, t182, 0, t157 + t247 * t236 + (-t250 - t359) * t240, t236 * t250 + t240 * t247 + t270, t207 * t279 + t328 * t353 + t271, t139 * t209 - g(2) * (t221 + t308) - g(3) * (-pkin(7) * t218 + t201 + t220) + t207 * t280 + t350 * t328, t37, t15, t92, t36, t91, 0, t56 * qJD(3) + t89 * qJDD(3) + t172 * t130 + t184 * t79 + t252, -qJD(3) * t57 - qJDD(3) * t90 + t132 * t172 + t184 * t80 + t273, -t130 * t57 - t132 * t56 - t79 * t90 - t80 * t89 + t259, t27 * t90 + t62 * t57 + t26 * t89 + t61 * t56 + t86 * t184 + t129 * t172 - g(2) * (t221 + t276) - g(3) * (t220 + t309), t6, t1, t28, t7, t29, 0, -t103 * t74 - t113 * t248 + t34 * t223 + t9 * t226 + t253, -t103 * t261 + t113 * t256 - t223 * t35 - t226 * t8 + t274, t248 * t35 - t256 * t34 + t261 * t9 + t74 * t8 + t262, t4 * t35 + t19 * t8 + t5 * t34 + t18 * t9 + t44 * t113 + t81 * t103 - g(2) * (t221 + t277) - g(3) * (t220 + t310); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, t269 + t272, (-t296 + (-qJD(2) + t227) * t303) * pkin(1) + t347, 0, 0, t141, t104, t181, t142, t182, 0, t157 + t249 * t236 + (-t251 - t359) * t240, t236 * t251 + t240 * t249 + t270, pkin(7) * t279 - t329 * t353 + t271, -t139 * pkin(2) - g(2) * t308 - g(3) * t201 + (t280 + t204) * pkin(7) - t350 * t329, t37, t15, t92, t36, t91, 0, -t130 * t295 + qJDD(3) * t109 - t208 * t79 + (t130 * t341 + t323) * qJD(3) + t252, -t132 * t295 - qJDD(3) * t110 - t208 * t80 + (t132 * t341 - t322) * qJD(3) + t273, -t109 * t80 - t110 * t79 - t130 * t322 - t132 * t323 + t259, t27 * t110 + t26 * t109 - t86 * t208 - g(2) * t276 - g(3) * t309 + t322 * t62 + t323 * t61 + (t211 - t295) * t129, t6, t1, t28, t7, t29, 0, -t125 * t248 + t223 * t38 + t330 * t226 - t267 * t74 + t253, t125 * t256 - t39 * t223 - t331 * t226 - t261 * t267 + t274, t248 * t39 - t256 * t38 + t261 * t330 + t331 * t74 + t262, -g(2) * t277 - g(3) * t310 + t44 * t125 + t330 * t18 + t331 * t19 + t267 * t81 + t5 * t38 + t4 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t291, t305 * t222, t313, t291, t312, qJDD(3), t236 * t255 - t335, g(1) * t236 + t240 * t255, 0, 0, t321, -t124 + t346, (t130 + t292) * qJD(3) + t254, -t321, -t263, qJDD(3), -g(1) * t214 - t63 * qJD(3) - t129 * t132 + t347 * t213 + (qJDD(3) * t233 - t130 * t317) * pkin(3) + t26, g(1) * t213 + t64 * qJD(3) + t129 * t130 + t347 * t214 + (-qJDD(3) * t232 - t132 * t317) * pkin(3) - t27, (t62 + t63) * t132 + (-t61 + t64) * t130 + (-t232 * t79 - t233 * t80) * pkin(3), -t61 * t63 - t62 * t64 + (-t335 + t232 * t27 + t233 * t26 + (-t129 * t227 + t347) * t236) * pkin(3), t332, t356, t357, -t332, t358, t223, t150 * t223 + t226 * t324 + t74 * t98 + t351, -t151 * t223 - t325 * t226 + t261 * t98 + t352, -t150 * t256 + t151 * t248 + (t18 + t325) * t74 + (-t19 + t324) * t261, t4 * t151 + t5 * t150 - t81 * t98 - g(1) * t287 + t325 * t19 + t324 * t18 - t347 * (-pkin(4) * t213 - t341); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t132 * qJD(3) + t263, (-t130 + t292) * qJD(3) + t254, -t124 - t346, t62 * t130 + t61 * t132 + t348 + t86, 0, 0, 0, 0, 0, 0, -t248 - t326, t256 + t327, -t333 - t334, -t18 * t261 - t19 * t74 + t348 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t332, t356, t357, -t332, t358, t223, t19 * t226 + t351, t18 * t226 + t352, 0, 0;];
tau_reg = t2;

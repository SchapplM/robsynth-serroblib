% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:05
% EndTime: 2019-12-05 18:46:13
% DurationCPUTime: 3.88s
% Computational Cost: add. (6875->440), mult. (16704->529), div. (0->0), fcn. (12096->12), ass. (0->237)
t214 = sin(qJ(3));
t217 = cos(qJ(3));
t218 = cos(qJ(2));
t288 = t217 * t218;
t268 = qJD(1) * t288;
t215 = sin(qJ(2));
t283 = qJD(1) * t215;
t143 = -t214 * t283 + t268;
t289 = t214 * t218;
t154 = t217 * t215 + t289;
t144 = t154 * qJD(1);
t213 = sin(qJ(4));
t242 = t154 * qJD(2);
t230 = -qJD(3) * t154 - t242;
t228 = t230 * qJD(1);
t290 = t214 * t215;
t248 = -t288 + t290;
t241 = t248 * qJDD(1);
t224 = -t241 + t228;
t322 = cos(qJ(4));
t267 = qJD(4) * t322;
t280 = qJD(4) * t213;
t208 = qJD(2) + qJD(3);
t249 = t208 * t290;
t276 = t218 * qJDD(1);
t278 = qJD(1) * qJD(2);
t265 = t218 * t278;
t277 = t215 * qJDD(1);
t331 = -t265 - t277;
t258 = -qJD(3) * t268 - t214 * t276 + t331 * t217;
t83 = qJD(1) * t249 + t258;
t30 = -t143 * t267 + t144 * t280 - t213 * t224 + t322 * t83;
t201 = qJD(4) + t208;
t99 = -t322 * t143 + t213 * t144;
t302 = t99 * t201;
t338 = -t30 + t302;
t246 = t213 * t143 + t322 * t144;
t324 = t246 ^ 2;
t97 = t99 ^ 2;
t337 = -t97 + t324;
t212 = qJ(2) + qJ(3);
t202 = sin(t212);
t203 = cos(t212);
t216 = sin(qJ(1));
t219 = cos(qJ(1));
t253 = g(1) * t219 + g(2) * t216;
t336 = -g(3) * t203 + t253 * t202;
t206 = t218 * pkin(2);
t197 = t206 + pkin(1);
t172 = qJD(1) * t197;
t335 = t197 * qJDD(1);
t300 = t99 * qJ(5);
t334 = t99 * t246;
t298 = t246 * t201;
t31 = qJD(4) * t246 - t213 * t83 - t322 * t224;
t332 = -t31 + t298;
t122 = -t143 * pkin(3) - t172;
t205 = qJ(4) + t212;
t193 = sin(t205);
t194 = cos(t205);
t293 = t194 * t219;
t294 = t194 * t216;
t207 = qJDD(2) + qJDD(3);
t323 = pkin(7) + pkin(6);
t174 = t323 * t218;
t161 = qJD(1) * t174;
t149 = t217 * t161;
t173 = t323 * t215;
t159 = qJD(1) * t173;
t306 = qJD(2) * pkin(2);
t151 = -t159 + t306;
t108 = t214 * t151 + t149;
t120 = qJDD(2) * pkin(2) + t323 * t331;
t266 = t215 * t278;
t121 = t323 * (-t266 + t276);
t49 = -qJD(3) * t108 + t217 * t120 - t214 * t121;
t21 = t207 * pkin(3) + t83 * pkin(8) + t49;
t281 = qJD(3) * t217;
t282 = qJD(3) * t214;
t48 = t214 * t120 + t217 * t121 + t151 * t281 - t161 * t282;
t27 = t224 * pkin(8) + t48;
t145 = t214 * t161;
t107 = t217 * t151 - t145;
t137 = t144 * pkin(8);
t81 = t107 - t137;
t69 = t208 * pkin(3) + t81;
t314 = t143 * pkin(8);
t82 = t108 + t314;
t4 = t213 * t21 + t69 * t267 + t322 * t27 - t82 * t280;
t245 = g(1) * t293 + g(2) * t294 + g(3) * t193 - t4;
t237 = t122 * t99 + t245;
t264 = t99 * pkin(4) + qJD(5);
t62 = t122 + t264;
t310 = t62 * t246;
t196 = t217 * pkin(2) + pkin(3);
t292 = t213 * t214;
t105 = t196 * t267 + (-t214 * t280 + (t322 * t217 - t292) * qJD(3)) * pkin(2);
t269 = t322 * t214;
t140 = pkin(2) * t269 + t213 * t196;
t200 = qJDD(4) + t207;
t329 = -t105 * t201 - t140 * t200;
t94 = t246 * qJ(5);
t124 = -t214 * t173 + t217 * t174;
t259 = pkin(3) * t267;
t320 = pkin(3) * t213;
t328 = -t200 * t320 - t201 * t259;
t190 = t200 * pkin(4);
t327 = -t246 * qJD(5) + t190;
t128 = pkin(3) * t248 - t197;
t295 = t193 * t219;
t296 = t193 * t216;
t317 = g(3) * t194;
t326 = g(1) * t295 + g(2) * t296 - t317;
t304 = t30 * qJ(5);
t325 = t304 + t327;
t263 = t322 * t21 - t213 * t27;
t73 = t322 * t82;
t41 = t213 * t69 + t73;
t5 = -qJD(4) * t41 + t263;
t229 = t5 + t326;
t225 = -t122 * t246 + t229;
t321 = pkin(3) * t202;
t315 = g(3) * t218;
t313 = t144 * pkin(3);
t312 = t154 * pkin(8);
t311 = t215 * pkin(2);
t71 = t213 * t82;
t40 = t322 * t69 - t71;
t23 = t40 - t94;
t18 = t201 * pkin(4) + t23;
t308 = t18 - t23;
t307 = -t99 * t259 - t31 * t320;
t43 = t322 * t81 - t71;
t116 = t214 * t159 - t149;
t85 = t116 - t314;
t117 = -t217 * t159 - t145;
t86 = -t137 + t117;
t45 = t213 * t85 + t322 * t86;
t291 = t214 * t174;
t123 = -t217 * t173 - t291;
t95 = t123 - t312;
t96 = -pkin(8) * t248 + t124;
t56 = t213 * t95 + t322 * t96;
t303 = t31 * qJ(5);
t301 = pkin(6) * qJDD(1);
t297 = t144 * t143;
t287 = t99 * qJD(5);
t192 = pkin(3) * t203;
t286 = t192 + t206;
t210 = t215 ^ 2;
t211 = t218 ^ 2;
t285 = t210 - t211;
t284 = t210 + t211;
t209 = -pkin(8) - t323;
t106 = -t196 * t280 + (-t214 * t267 + (-t213 * t217 - t269) * qJD(3)) * pkin(2);
t275 = -t105 * t99 - t106 * t246 - t140 * t31;
t199 = t215 * t306;
t222 = qJD(1) ^ 2;
t274 = t215 * t222 * t218;
t270 = qJD(2) * t323;
t160 = t215 * t270;
t162 = t218 * t270;
t273 = -t217 * t160 - t214 * t162 - t173 * t281;
t188 = pkin(4) * t194;
t272 = t188 + t286;
t42 = -t213 * t81 - t73;
t44 = -t213 * t86 + t322 * t85;
t55 = -t213 * t96 + t322 * t95;
t261 = t208 * t215;
t189 = pkin(2) * t266;
t138 = t189 - t335;
t260 = -t138 + t335;
t257 = t215 * t265;
t139 = -pkin(2) * t292 + t322 * t196;
t256 = -g(1) * t296 + g(2) * t295;
t255 = g(1) * t294 - g(2) * t293;
t68 = pkin(4) * t246 + t313;
t254 = -pkin(4) * t193 - t321;
t252 = g(1) * t216 - g(2) * t219;
t24 = t41 - t300;
t251 = -t18 * t99 + t24 * t246;
t250 = t246 * t41 - t40 * t99;
t247 = -0.2e1 * pkin(1) * t278 - pkin(6) * qJDD(2);
t53 = -pkin(8) * t242 + (-t291 - t312) * qJD(3) + t273;
t118 = -qJD(2) * t288 - t218 * t281 + t249;
t67 = -t124 * qJD(3) + t214 * t160 - t217 * t162;
t54 = t118 * pkin(8) + t67;
t10 = t213 * t54 + t95 * t267 - t96 * t280 + t322 * t53;
t239 = t322 * t248;
t238 = t245 + t303;
t221 = qJD(2) ^ 2;
t236 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t221 + t252;
t235 = pkin(1) * t222 + t253 - t301;
t234 = g(3) * t202 + t172 * t143 + t253 * t203 - t48;
t233 = -0.2e1 * t172 * t208;
t115 = t322 * t154 - t213 * t248;
t11 = -t56 * qJD(4) - t213 * t53 + t322 * t54;
t231 = t62 * t99 + t238 + t287;
t104 = -pkin(3) * t230 + t199;
t227 = t172 * t144 + t336 + t49;
t226 = t229 + t304;
t61 = -pkin(3) * t228 + t128 * qJDD(1) + t189;
t12 = t31 * pkin(4) + qJDD(5) + t61;
t204 = -qJ(5) + t209;
t198 = pkin(2) * t283;
t195 = t322 * pkin(3) + pkin(4);
t158 = pkin(1) + t286;
t135 = pkin(4) + t139;
t130 = pkin(1) + t272;
t125 = t198 + t313;
t114 = t213 * t154 + t239;
t93 = t106 * t201;
t84 = -t143 ^ 2 + t144 ^ 2;
t74 = t114 * pkin(4) + t128;
t66 = -t174 * t282 + t273;
t63 = t198 + t68;
t60 = t144 * t208 + t224;
t59 = -qJD(1) * t214 * t261 - t143 * t208 - t258;
t47 = t115 * qJD(4) - t213 * t118 - t322 * t230;
t46 = qJD(4) * t239 + t322 * t118 + t154 * t280 - t213 * t230;
t38 = -t114 * qJ(5) + t56;
t37 = -t115 * qJ(5) + t55;
t36 = t47 * pkin(4) + t104;
t35 = -t94 + t45;
t34 = t44 + t300;
t33 = -t114 * t200 - t47 * t201;
t32 = t115 * t200 - t46 * t201;
t29 = -t94 + t43;
t28 = t42 + t300;
t9 = t31 * t114 + t99 * t47;
t8 = -t30 * t115 - t246 * t46;
t7 = t46 * qJ(5) - t115 * qJD(5) + t11;
t6 = -t47 * qJ(5) - t114 * qJD(5) + t10;
t3 = -t287 + t4 - t303;
t2 = t5 + t325;
t1 = t30 * t114 - t115 * t31 - t246 * t47 + t46 * t99;
t13 = [0, 0, 0, 0, 0, qJDD(1), t252, t253, 0, 0, t210 * qJDD(1) + 0.2e1 * t257, 0.2e1 * t215 * t276 - 0.2e1 * t285 * t278, qJDD(2) * t215 + t221 * t218, t211 * qJDD(1) - 0.2e1 * t257, qJDD(2) * t218 - t221 * t215, 0, t215 * t247 + t218 * t236, -t215 * t236 + t218 * t247, 0.2e1 * t284 * t301 - t253, -g(1) * (-t216 * pkin(1) + t219 * pkin(6)) - g(2) * (t219 * pkin(1) + t216 * pkin(6)) + (pkin(6) ^ 2 * t284 + pkin(1) ^ 2) * qJDD(1), -t144 * t118 - t83 * t154, -t118 * t143 + t144 * t230 + t224 * t154 + t248 * t83, -t118 * t208 + t154 * t207, t143 * t230 - t224 * t248, -t207 * t248 + t208 * t230, 0, t123 * t207 + t67 * t208 + t252 * t203 + (t214 * t233 + t217 * t260) * t218 + (-t143 * t306 - t214 * t260 + t217 * t233) * t215, t172 * t118 - t124 * t207 + t138 * t154 + t144 * t199 + t197 * t83 - t202 * t252 - t66 * t208, t107 * t118 + t108 * t230 + t123 * t83 + t224 * t124 + t66 * t143 - t67 * t144 - t49 * t154 - t248 * t48 - t253, t48 * t124 + t108 * t66 + t49 * t123 + t107 * t67 - t138 * t197 - t172 * t199 - g(1) * (-t216 * t197 + t219 * t323) - g(2) * (t219 * t197 + t216 * t323), t8, t1, t32, t9, t33, 0, t104 * t99 + t11 * t201 + t61 * t114 + t122 * t47 + t128 * t31 + t55 * t200 + t255, -t10 * t201 + t104 * t246 + t61 * t115 - t122 * t46 - t128 * t30 - t56 * t200 + t256, -t10 * t99 - t11 * t246 - t4 * t114 - t5 * t115 + t55 * t30 - t56 * t31 + t40 * t46 - t41 * t47 - t253, t4 * t56 + t41 * t10 + t5 * t55 + t40 * t11 + t61 * t128 + t122 * t104 - g(1) * (-t216 * t158 - t219 * t209) - g(2) * (t219 * t158 - t216 * t209), t8, t1, t32, t9, t33, 0, t12 * t114 + t37 * t200 + t7 * t201 + t74 * t31 + t36 * t99 + t62 * t47 + t255, t12 * t115 - t38 * t200 - t6 * t201 + t246 * t36 - t74 * t30 - t62 * t46 + t256, -t3 * t114 - t2 * t115 + t18 * t46 - t24 * t47 - t246 * t7 + t37 * t30 - t38 * t31 - t6 * t99 - t253, t3 * t38 + t24 * t6 + t2 * t37 + t18 * t7 + t12 * t74 + t62 * t36 - g(1) * (-t216 * t130 - t219 * t204) - g(2) * (t219 * t130 - t216 * t204); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t274, t285 * t222, t277, t274, t276, qJDD(2), t215 * t235 - t315, g(3) * t215 + t218 * t235, 0, 0, -t297, t84, t59, t297, t60, t207, -t116 * t208 + (t143 * t283 + t217 * t207 - t208 * t282) * pkin(2) + t227, t117 * t208 + (-t144 * t283 - t214 * t207 - t208 * t281) * pkin(2) + t234, (t108 + t116) * t144 + (t107 - t117) * t143 + ((qJD(3) * t143 + t83) * t217 + (qJD(3) * t144 - t241 + (-t208 * t289 - t217 * t261) * qJD(1)) * t214) * pkin(2), -t107 * t116 - t108 * t117 + (-t315 + t48 * t214 + t217 * t49 + (-t107 * t214 + t108 * t217) * qJD(3) + (qJD(1) * t172 + t253) * t215) * pkin(2), t334, t337, t338, -t334, t332, t200, -t125 * t99 + t139 * t200 - t44 * t201 + t225 + t93, -t125 * t246 + t45 * t201 + t237 + t329, t139 * t30 + t246 * t44 + t45 * t99 + t250 + t275, t4 * t140 + t5 * t139 - t122 * t125 - g(3) * t286 + (t105 - t45) * t41 + (t106 - t44) * t40 - t253 * (-t311 - t321), t334, t337, t338, -t334, t332, t200, t135 * t200 - t34 * t201 - t63 * t99 + t226 - t310 + t327 + t93, t35 * t201 - t246 * t63 + t231 + t329, t135 * t30 + t246 * t34 + t35 * t99 + t251 + t275, t3 * t140 + t2 * t135 - t62 * t63 - g(3) * t272 + (t105 - t35) * t24 + (t106 - t34) * t18 - t253 * (t254 - t311); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t297, t84, t59, t297, t60, t207, t108 * t208 + t227, t107 * t208 + t234, 0, 0, t334, t337, t338, -t334, t332, t200, -t42 * t201 + (-t144 * t99 + t322 * t200 - t201 * t280) * pkin(3) + t225, t43 * t201 - t246 * t313 + t237 + t328, t42 * t246 + t43 * t99 + (t246 * t280 + t322 * t30) * pkin(3) + t250 + t307, -t40 * t42 - t41 * t43 + (t322 * t5 - t122 * t144 + t213 * t4 + (-t213 * t40 + t322 * t41) * qJD(4) + t336) * pkin(3), t334, t337, t338, -t334, t332, t200, t195 * t200 - t28 * t201 - t310 - t68 * t99 + (-t73 + (-pkin(3) * t201 - t69) * t213) * qJD(4) + t263 + t325 + t326, t29 * t201 - t246 * t68 + t231 + t328, t195 * t30 + t29 * t99 + (pkin(3) * t280 + t28) * t246 + t251 + t307, t2 * t195 - t24 * t29 - t18 * t28 - t62 * t68 - g(3) * (t188 + t192) - t253 * t254 + (t3 * t213 + (-t18 * t213 + t322 * t24) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t334, t337, t338, -t334, t332, t200, t41 * t201 + t225, t40 * t201 + t237, 0, 0, t334, t337, t338, -t334, t332, t200, t24 * t201 + 0.2e1 * t190 + (-t264 - t62) * t246 + t226, -t324 * pkin(4) + t23 * t201 + (qJD(5) + t62) * t99 + t238, t30 * pkin(4) - t308 * t99, t308 * t24 + (t193 * t253 + t2 - t310 - t317) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 + t298, -t30 - t302, -t97 - t324, t18 * t246 + t24 * t99 + t12 - t252;];
tau_reg = t13;

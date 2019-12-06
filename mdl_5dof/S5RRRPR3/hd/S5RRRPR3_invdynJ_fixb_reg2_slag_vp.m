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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:43:03
% EndTime: 2019-12-05 18:43:10
% DurationCPUTime: 3.16s
% Computational Cost: add. (6325->390), mult. (9824->487), div. (0->0), fcn. (6824->16), ass. (0->240)
t232 = cos(qJ(3));
t210 = t232 * qJD(4);
t228 = sin(qJ(3));
t226 = -qJ(4) - pkin(7);
t276 = qJD(3) * t226;
t150 = t228 * t276 + t210;
t151 = -t228 * qJD(4) + t232 * t276;
t224 = sin(pkin(9));
t225 = cos(pkin(9));
t161 = t224 * t232 + t225 * t228;
t233 = cos(qJ(2));
t296 = qJD(1) * t233;
t286 = pkin(1) * t296;
t314 = -t150 * t224 + t225 * t151 + t161 * t286;
t305 = t225 * t232;
t307 = t224 * t228;
t160 = -t305 + t307;
t313 = t225 * t150 + t224 * t151 + t160 * t286;
t227 = sin(qJ(5));
t231 = cos(qJ(5));
t219 = qJD(1) + qJD(2);
t285 = t219 * t305;
t130 = t219 * t307 - t285;
t132 = t161 * t219;
t254 = t130 * t227 - t231 * t132;
t153 = t161 * qJD(3);
t216 = qJDD(1) + qJDD(2);
t303 = t232 * t216;
t304 = t228 * t216;
t256 = t224 * t304 - t225 * t303;
t79 = t153 * t219 + t256;
t293 = qJD(3) * t228;
t281 = t219 * t293;
t246 = t161 * t216 - t224 * t281;
t292 = qJD(3) * t232;
t280 = t219 * t292;
t80 = t225 * t280 + t246;
t240 = qJD(5) * t254 - t227 * t80 - t231 * t79;
t218 = qJD(3) + qJD(5);
t317 = t254 * t218;
t349 = t240 - t317;
t290 = qJD(5) * t231;
t291 = qJD(5) * t227;
t248 = -t130 * t290 - t132 * t291 - t227 * t79 + t231 * t80;
t74 = -t231 * t130 - t132 * t227;
t318 = t218 * t74;
t348 = t248 - t318;
t323 = t74 ^ 2;
t324 = t254 ^ 2;
t347 = -t323 + t324;
t322 = t74 * t254;
t154 = t160 * qJD(3);
t327 = pkin(8) * t154;
t346 = t327 + t314;
t147 = t153 * pkin(8);
t345 = -t147 + t313;
t221 = t228 ^ 2;
t222 = t232 ^ 2;
t297 = t221 + t222;
t340 = t297 * t233;
t344 = t219 * t340;
t229 = sin(qJ(2));
t289 = qJDD(1) * t229;
t294 = qJD(2) * t233;
t140 = t216 * pkin(7) + (qJD(1) * t294 + t289) * pkin(1);
t273 = t297 * t140;
t223 = qJ(1) + qJ(2);
t211 = sin(t223);
t212 = cos(t223);
t259 = g(2) * t212 + g(3) * t211;
t199 = g(2) * t211;
t339 = -g(3) * t212 + t199;
t220 = qJ(3) + pkin(9);
t209 = qJ(5) + t220;
t193 = sin(t209);
t194 = cos(t209);
t309 = t194 * t212;
t310 = t194 * t211;
t249 = qJ(4) * t216 + qJD(4) * t219 + t140;
t319 = pkin(1) * qJD(1);
t288 = t229 * t319;
t171 = pkin(7) * t219 + t288;
t270 = qJ(4) * t219 + t171;
t251 = qJD(3) * t270;
t55 = qJDD(3) * pkin(3) - t228 * t249 - t232 * t251;
t58 = -t228 * t251 + t232 * t249;
t26 = -t224 * t58 + t225 * t55;
t16 = qJDD(3) * pkin(4) - pkin(8) * t80 + t26;
t27 = t224 * t55 + t225 * t58;
t17 = -pkin(8) * t79 + t27;
t328 = pkin(8) * t132;
t118 = t270 * t232;
t105 = t224 * t118;
t117 = t270 * t228;
t111 = qJD(3) * pkin(3) - t117;
t61 = t225 * t111 - t105;
t43 = qJD(3) * pkin(4) - t328 + t61;
t329 = pkin(8) * t130;
t306 = t225 * t118;
t62 = t224 * t111 + t306;
t45 = t62 - t329;
t4 = (qJD(5) * t43 + t17) * t231 + t227 * t16 - t45 * t291;
t330 = pkin(3) * t232;
t202 = pkin(2) + t330;
t129 = -t202 * t219 + qJD(4) - t286;
t81 = t130 * pkin(4) + t129;
t343 = g(1) * t193 - g(2) * t310 + g(3) * t309 - t74 * t81 - t4;
t19 = t227 * t43 + t231 * t45;
t5 = -qJD(5) * t19 + t231 * t16 - t227 * t17;
t342 = -g(1) * t194 - t339 * t193 + t254 * t81 + t5;
t333 = pkin(2) * t219;
t172 = -t286 - t333;
t341 = t171 * t340 + t172 * t229;
t338 = t132 ^ 2;
t230 = sin(qJ(1));
t337 = pkin(1) * t230;
t336 = pkin(1) * t233;
t234 = cos(qJ(1));
t335 = pkin(1) * t234;
t334 = pkin(2) * t216;
t332 = pkin(3) * t224;
t331 = pkin(3) * t228;
t326 = pkin(8) * t161;
t325 = g(1) * t232;
t182 = t226 * t228;
t213 = t232 * qJ(4);
t183 = pkin(7) * t232 + t213;
t109 = t225 * t182 - t183 * t224;
t84 = t109 - t326;
t110 = t224 * t182 + t225 * t183;
t157 = t160 * pkin(8);
t85 = -t157 + t110;
t38 = -t227 * t85 + t231 * t84;
t321 = qJD(5) * t38 + t346 * t227 + t345 * t231;
t39 = t227 * t84 + t231 * t85;
t320 = -qJD(5) * t39 - t345 * t227 + t346 * t231;
t201 = pkin(1) * t229 + pkin(7);
t302 = -qJ(4) - t201;
t267 = qJD(3) * t302;
t287 = pkin(1) * t294;
t100 = (-qJD(4) - t287) * t228 + t232 * t267;
t99 = t228 * t267 + t232 * t287 + t210;
t57 = t224 * t100 + t225 * t99;
t195 = pkin(3) * t225 + pkin(4);
t148 = t195 * t231 - t227 * t332;
t63 = t117 * t224 - t306;
t46 = t63 + t329;
t64 = -t225 * t117 - t105;
t47 = t64 - t328;
t316 = t148 * qJD(5) - t227 * t46 - t231 * t47;
t149 = t195 * t227 + t231 * t332;
t315 = -t149 * qJD(5) + t227 * t47 - t231 * t46;
t312 = t132 * t130;
t308 = t219 * t228;
t295 = qJD(2) * t229;
t206 = pkin(1) * t295;
t299 = -qJD(1) * t206 + qJDD(1) * t336;
t139 = -t299 - t334;
t301 = t139 * t228 + t172 * t292;
t158 = t302 * t228;
t159 = t201 * t232 + t213;
t90 = t224 * t158 + t225 * t159;
t298 = t221 - t222;
t205 = pkin(3) * t293;
t214 = t219 ^ 2;
t284 = t228 * t214 * t232;
t283 = t172 * t293 + t259 * t232;
t282 = t219 * t295;
t208 = cos(t220);
t279 = pkin(4) * t208 + t330;
t114 = pkin(4) * t153 + t205;
t56 = t225 * t100 - t224 * t99;
t272 = t297 * t216;
t89 = t225 * t158 - t159 * t224;
t165 = pkin(2) + t279;
t217 = -pkin(8) + t226;
t269 = -t165 * t212 + t211 * t217;
t268 = -t202 * t212 + t211 * t226;
t86 = pkin(3) * t281 - t202 * t216 + qJDD(4) - t299;
t44 = t79 * pkin(4) + t86;
t94 = -t160 * t227 + t161 * t231;
t49 = qJD(5) * t94 + t231 * t153 - t227 * t154;
t93 = t231 * t160 + t161 * t227;
t266 = g(2) * t309 + g(3) * t310 + t44 * t93 + t81 * t49;
t265 = t129 * t153 + t86 * t160 + t259 * t208;
t264 = t219 * t288;
t263 = t339 + t273;
t262 = t299 + t259;
t261 = t228 * t280;
t260 = t114 - t288;
t257 = g(2) * t234 + g(3) * t230;
t67 = t89 - t326;
t68 = -t157 + t90;
t34 = -t227 * t68 + t231 * t67;
t35 = t227 * t67 + t231 * t68;
t18 = -t227 * t45 + t231 * t43;
t48 = t227 * t153 + t231 * t154 + t160 * t290 + t161 * t291;
t255 = t18 * t48 - t19 * t49 - t4 * t93 - t5 * t94 + t339;
t253 = -t165 * t211 - t212 * t217;
t252 = -t202 * t211 - t212 * t226;
t125 = pkin(4) * t160 - t202;
t250 = -t62 * t153 + t61 * t154 - t27 * t160 - t26 * t161 + t339;
t247 = -t172 * t219 - t140 - t339;
t245 = -t193 * t259 + t44 * t94 - t81 * t48;
t207 = sin(t220);
t244 = -t129 * t154 + t86 * t161 - t207 * t259;
t235 = qJD(3) ^ 2;
t243 = -pkin(7) * t235 + t264 + t334;
t203 = -pkin(2) - t336;
t242 = -pkin(1) * t282 - t201 * t235 - t203 * t216;
t241 = -pkin(7) * qJDD(3) + (t286 - t333) * qJD(3);
t239 = -qJDD(3) * t201 + (t203 * t219 - t287) * qJD(3);
t215 = qJDD(3) + qJDD(5);
t196 = t212 * pkin(7);
t180 = -t202 - t336;
t178 = qJDD(3) * t228 + t232 * t235;
t177 = qJDD(3) * t232 - t228 * t235;
t170 = t206 + t205;
t142 = t216 * t222 - 0.2e1 * t261;
t141 = t216 * t221 + 0.2e1 * t261;
t124 = t130 ^ 2;
t113 = t125 - t336;
t104 = -0.2e1 * qJD(3) * t219 * t298 + 0.2e1 * t228 * t303;
t103 = t114 + t206;
t98 = pkin(3) * t308 + pkin(4) * t132;
t92 = -qJD(3) * t154 + qJDD(3) * t161;
t91 = -qJD(3) * t153 - qJDD(3) * t160;
t41 = -t147 + t57;
t40 = t56 + t327;
t37 = -t132 * t154 + t161 * t80;
t36 = t130 * t153 + t160 * t79;
t29 = -t215 * t93 - t218 * t49;
t28 = t215 * t94 - t218 * t48;
t15 = t130 * t154 - t132 * t153 - t160 * t80 - t161 * t79;
t9 = -qJD(5) * t35 - t227 * t41 + t231 * t40;
t8 = qJD(5) * t34 + t227 * t40 + t231 * t41;
t7 = -t240 * t93 - t49 * t74;
t6 = t248 * t94 + t254 * t48;
t1 = t240 * t94 - t248 * t93 + t254 * t49 - t48 * t74;
t2 = [0, 0, 0, 0, 0, qJDD(1), t257, -g(2) * t230 + g(3) * t234, 0, 0, 0, 0, 0, 0, 0, t216, (t216 * t233 - t282) * pkin(1) + t262, ((-qJDD(1) - t216) * t229 + (-qJD(1) - t219) * t294) * pkin(1) - t339, 0, (t257 + (t229 ^ 2 + t233 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t141, t104, t178, t142, t177, 0, t239 * t228 + (-t139 + t242) * t232 + t283, t239 * t232 + (-t242 - t259) * t228 + t301, pkin(1) * qJD(2) * t344 + t201 * t272 + t263, t139 * t203 - g(2) * (-pkin(2) * t212 - pkin(7) * t211) - g(3) * (-pkin(2) * t211 + t196) + t201 * t273 + (t341 * qJD(2) + t257) * pkin(1), t37, t15, t92, t36, t91, 0, qJD(3) * t56 + qJDD(3) * t89 + t130 * t170 + t180 * t79 + t265, -t57 * qJD(3) - t90 * qJDD(3) + t170 * t132 + t180 * t80 + t244, -t130 * t57 - t132 * t56 - t79 * t90 - t80 * t89 + t250, t27 * t90 + t62 * t57 + t26 * t89 + t61 * t56 + t86 * t180 + t129 * t170 - g(2) * (t268 - t335) - g(3) * (t252 - t337), t6, t1, t28, t7, t29, 0, -t103 * t74 - t113 * t240 + t215 * t34 + t218 * t9 + t266, -t103 * t254 + t113 * t248 - t35 * t215 - t8 * t218 + t245, t240 * t35 - t248 * t34 + t254 * t9 + t74 * t8 + t255, t4 * t35 + t19 * t8 + t5 * t34 + t18 * t9 + t44 * t113 + t81 * t103 - g(2) * (t269 - t335) - g(3) * (t253 - t337); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, t262 + t264, (-t289 + (-qJD(2) + t219) * t296) * pkin(1) - t339, 0, 0, t141, t104, t178, t142, t177, 0, t241 * t228 + (-t139 + t243) * t232 + t283, t241 * t232 + (-t243 - t259) * t228 + t301, pkin(7) * t272 - t319 * t344 + t263, -g(3) * t196 + (-t139 + t259) * pkin(2) + (t273 + t199) * pkin(7) - t341 * t319, t37, t15, t92, t36, t91, 0, -t130 * t288 + qJDD(3) * t109 - t202 * t79 + (t130 * t331 + t314) * qJD(3) + t265, -t132 * t288 - qJDD(3) * t110 - t202 * t80 + (t132 * t331 - t313) * qJD(3) + t244, -t109 * t80 - t110 * t79 - t130 * t313 - t132 * t314 + t250, t27 * t110 + t26 * t109 - t86 * t202 - g(2) * t268 - g(3) * t252 + t313 * t62 + t314 * t61 + (t205 - t288) * t129, t6, t1, t28, t7, t29, 0, -t125 * t240 + t215 * t38 + t320 * t218 - t260 * t74 + t266, t125 * t248 - t215 * t39 - t321 * t218 - t254 * t260 + t245, t240 * t39 - t248 * t38 + t254 * t320 + t321 * t74 + t255, -g(2) * t269 - g(3) * t253 + t44 * t125 + t320 * t18 + t321 * t19 + t260 * t81 + t5 * t38 + t4 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t284, t298 * t214, t304, t284, t303, qJDD(3), t228 * t247 - t325, g(1) * t228 + t232 * t247, 0, 0, t312, -t124 + t338, (t130 + t285) * qJD(3) + t246, -t312, -t256, qJDD(3), -g(1) * t208 - qJD(3) * t63 - t129 * t132 - t339 * t207 + (qJDD(3) * t225 - t130 * t308) * pkin(3) + t26, g(1) * t207 + qJD(3) * t64 + t129 * t130 - t339 * t208 + (-qJDD(3) * t224 - t132 * t308) * pkin(3) - t27, (t62 + t63) * t132 + (-t61 + t64) * t130 + (-t224 * t79 - t225 * t80) * pkin(3), -t61 * t63 - t62 * t64 + (-t325 + t224 * t27 + t225 * t26 + (-t129 * t219 - t339) * t228) * pkin(3), t322, t347, t348, -t322, t349, t215, t148 * t215 + t218 * t315 + t74 * t98 + t342, -t149 * t215 - t218 * t316 + t254 * t98 + t343, -t148 * t248 + t149 * t240 + (t18 + t316) * t74 + (-t19 + t315) * t254, t4 * t149 + t5 * t148 - t81 * t98 - g(1) * t279 + t316 * t19 + t315 * t18 + t339 * (-pkin(4) * t207 - t331); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t132 * qJD(3) + t256, (-t130 + t285) * qJD(3) + t246, -t124 - t338, t130 * t62 + t132 * t61 - t259 + t86, 0, 0, 0, 0, 0, 0, -t240 - t317, t248 + t318, -t323 - t324, -t18 * t254 - t19 * t74 - t259 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, t347, t348, -t322, t349, t215, t19 * t218 + t342, t18 * t218 + t343, 0, 0;];
tau_reg = t2;

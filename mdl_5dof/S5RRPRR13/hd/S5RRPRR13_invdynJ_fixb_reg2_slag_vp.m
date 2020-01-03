% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR13
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR13_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:33:58
% EndTime: 2019-12-31 20:34:13
% DurationCPUTime: 7.57s
% Computational Cost: add. (9433->550), mult. (21776->739), div. (0->0), fcn. (15841->14), ass. (0->252)
t240 = sin(pkin(9));
t245 = sin(qJ(2));
t316 = qJD(1) * t245;
t300 = t240 * t316;
t241 = cos(pkin(9));
t308 = t241 * qJD(2);
t180 = t300 - t308;
t298 = t241 * t316;
t314 = qJD(2) * t240;
t182 = t298 + t314;
t244 = sin(qJ(4));
t247 = cos(qJ(4));
t107 = t180 * t244 - t182 * t247;
t108 = t247 * t180 + t182 * t244;
t243 = sin(qJ(5));
t355 = cos(qJ(5));
t54 = t355 * t107 + t243 * t108;
t376 = t54 ^ 2;
t58 = t243 * t107 - t355 * t108;
t375 = t58 ^ 2;
t248 = cos(qJ(2));
t315 = qJD(1) * t248;
t215 = -qJD(4) + t315;
t208 = -qJD(5) + t215;
t374 = t208 * t58;
t373 = t54 * t208;
t189 = t240 * t247 + t241 * t244;
t169 = t189 * qJD(4);
t266 = t189 * t248;
t321 = qJD(1) * t266 - t169;
t323 = t247 * t241;
t188 = t240 * t244 - t323;
t310 = qJD(4) * t247;
t311 = qJD(4) * t244;
t361 = -t240 * t311 + t241 * t310;
t366 = t188 * t315 + t361;
t304 = t245 * qJDD(1);
t221 = pkin(6) * t304;
t306 = qJD(1) * qJD(2);
t293 = t248 * t306;
t160 = -qJDD(2) * pkin(2) + pkin(6) * t293 + qJDD(3) + t221;
t235 = g(3) * t248;
t246 = sin(qJ(1));
t249 = cos(qJ(1));
t282 = g(1) * t249 + g(2) * t246;
t261 = t282 * t245 - t235;
t256 = t160 - t261;
t279 = pkin(2) * t245 - qJ(3) * t248;
t192 = t279 * qJD(1);
t170 = t240 * t192;
t327 = t241 * t245;
t328 = t240 * t248;
t267 = -pkin(6) * t327 - pkin(7) * t328;
t114 = qJD(1) * t267 + t170;
t342 = pkin(7) + qJ(3);
t201 = t342 * t240;
t202 = t342 * t241;
t129 = pkin(6) * t300 + t241 * t192;
t326 = t241 * t248;
t271 = pkin(3) * t245 - pkin(7) * t326;
t96 = qJD(1) * t271 + t129;
t337 = qJD(3) * t323 - t247 * t114 - t201 * t310 + (-t240 * qJD(3) - qJD(4) * t202 - t96) * t244;
t124 = -t244 * t201 + t247 * t202;
t336 = -t189 * qJD(3) - t124 * qJD(4) + t114 * t244 - t247 * t96;
t372 = t107 ^ 2;
t371 = t108 ^ 2;
t370 = -pkin(4) * t316 - t366 * pkin(8) + t336;
t369 = -t321 * pkin(8) - t337;
t343 = t58 * t54;
t368 = t107 * t215;
t367 = t108 * t215;
t365 = -t375 + t376;
t294 = qJD(5) * t355;
t309 = qJD(5) * t243;
t225 = t241 * qJDD(2);
t265 = t293 + t304;
t134 = t240 * t265 - t225;
t305 = qJDD(2) * t240;
t135 = t241 * t265 + t305;
t48 = t244 * t134 - t247 * t135 + t180 * t310 + t182 * t311;
t49 = t134 * t247 + t244 * t135 - t180 * t311 + t182 * t310;
t14 = -t107 * t309 + t108 * t294 + t243 * t49 + t355 * t48;
t364 = -t14 + t374;
t280 = pkin(2) * t248 + qJ(3) * t245;
t196 = -pkin(1) - t280;
t172 = t196 * qJD(1);
t223 = pkin(6) * t315;
t203 = qJD(2) * qJ(3) + t223;
t116 = t241 * t172 - t203 * t240;
t72 = -pkin(3) * t315 - pkin(7) * t182 + t116;
t117 = t240 * t172 + t241 * t203;
t75 = -pkin(7) * t180 + t117;
t34 = -t244 * t75 + t247 * t72;
t28 = pkin(8) * t107 + t34;
t26 = -pkin(4) * t215 + t28;
t35 = t244 * t72 + t247 * t75;
t29 = -pkin(8) * t108 + t35;
t229 = t248 * qJDD(1);
t292 = t245 * t306;
t264 = t292 - t229;
t187 = qJDD(4) + t264;
t165 = qJD(2) * t279 - qJD(3) * t245;
t106 = qJD(1) * t165 + qJDD(1) * t196;
t148 = -pkin(6) * t264 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t64 = t241 * t106 - t148 * t240;
t39 = pkin(3) * t264 - pkin(7) * t135 + t64;
t65 = t240 * t106 + t241 * t148;
t50 = -pkin(7) * t134 + t65;
t9 = -qJD(4) * t35 - t244 * t50 + t247 * t39;
t6 = pkin(4) * t187 + pkin(8) * t48 + t9;
t8 = t244 * t39 + t247 * t50 + t72 * t310 - t311 * t75;
t7 = -pkin(8) * t49 + t8;
t1 = t243 * t6 + t26 * t294 - t29 * t309 + t355 * t7;
t237 = pkin(9) + qJ(4);
t230 = qJ(5) + t237;
t218 = sin(t230);
t219 = cos(t230);
t324 = t246 * t248;
t137 = t218 * t249 - t219 * t324;
t322 = t248 * t249;
t139 = t218 * t246 + t219 * t322;
t344 = g(3) * t245;
t195 = -qJD(2) * pkin(2) + pkin(6) * t316 + qJD(3);
t128 = pkin(3) * t180 + t195;
t68 = pkin(4) * t108 + t128;
t363 = g(1) * t139 - g(2) * t137 + t219 * t344 - t58 * t68 - t1;
t179 = t241 * t196;
t115 = -pkin(7) * t327 + t179 + (-pkin(6) * t240 - pkin(3)) * t248;
t141 = pkin(6) * t326 + t240 * t196;
t329 = t240 * t245;
t122 = -pkin(7) * t329 + t141;
t63 = t244 * t115 + t247 * t122;
t238 = t245 ^ 2;
t239 = t248 ^ 2;
t318 = t238 - t239;
t285 = qJD(1) * t318;
t227 = sin(t237);
t228 = cos(t237);
t149 = t227 * t324 + t228 * t249;
t151 = -t227 * t322 + t228 * t246;
t360 = -g(1) * t151 + g(2) * t149 + t227 * t344;
t136 = t218 * t324 + t219 * t249;
t138 = -t218 * t322 + t219 * t246;
t301 = t355 * t29;
t11 = t243 * t26 + t301;
t2 = -qJD(5) * t11 - t243 * t7 + t355 * t6;
t359 = -g(1) * t138 + g(2) * t136 + t218 * t344 + t54 * t68 + t2;
t15 = -qJD(5) * t54 - t243 * t48 + t355 * t49;
t358 = -t15 + t373;
t357 = -0.2e1 * pkin(1);
t356 = pkin(4) * t49;
t354 = pkin(3) * t134;
t353 = pkin(3) * t240;
t352 = pkin(6) * t180;
t351 = pkin(6) * t182;
t349 = g(1) * t246;
t345 = g(2) * t249;
t123 = -t247 * t201 - t202 * t244;
t89 = -pkin(8) * t189 + t123;
t90 = -pkin(8) * t188 + t124;
t45 = -t243 * t90 + t355 * t89;
t341 = t45 * qJD(5) + t243 * t370 - t355 * t369;
t46 = t243 * t89 + t355 * t90;
t340 = -t46 * qJD(5) + t243 * t369 + t355 * t370;
t339 = t188 * t294 + t189 * t309 - t243 * t321 - t355 * t366;
t113 = -t243 * t188 + t355 * t189;
t338 = t113 * qJD(5) + t243 * t366 - t321 * t355;
t335 = t243 * t29;
t334 = pkin(6) * qJDD(1);
t333 = t107 * t108;
t332 = t134 * t241;
t331 = t135 * t240;
t251 = qJD(1) ^ 2;
t330 = t239 * t251;
t325 = t245 * t249;
t313 = qJD(2) * t245;
t302 = pkin(6) * t313;
t120 = t241 * t165 + t240 * t302;
t299 = t240 * t315;
t174 = pkin(3) * t299 + t223;
t312 = qJD(2) * t248;
t297 = t240 * t312;
t175 = pkin(3) * t297 + pkin(6) * t312;
t193 = pkin(3) * t329 + t245 * pkin(6);
t319 = t249 * pkin(1) + t246 * pkin(6);
t317 = t238 + t239;
t307 = qJD(3) - t195;
t220 = pkin(3) * t241 + pkin(2);
t291 = t240 * t304;
t290 = t240 * t229;
t289 = t241 * t304;
t288 = t241 * t229;
t287 = -pkin(4) * t321 - t174;
t62 = t247 * t115 - t122 * t244;
t284 = t248 * t292;
t216 = t245 * t349;
t283 = -g(2) * t325 + t216;
t281 = -t345 + t349;
t278 = t180 * t241 + t182 * t240;
t191 = pkin(4) * t228 + t220;
t236 = -pkin(8) - t342;
t277 = t191 * t248 - t236 * t245;
t275 = t220 * t248 + t245 * t342;
t273 = qJD(1) * (t180 + t308);
t272 = qJD(1) * (-t182 + t314);
t159 = t188 * t245;
t40 = -pkin(4) * t248 + pkin(8) * t159 + t62;
t158 = t189 * t245;
t47 = -pkin(8) * t158 + t63;
t22 = -t243 * t47 + t355 * t40;
t23 = t243 * t40 + t355 * t47;
t269 = -pkin(6) * qJDD(2) + t306 * t357;
t84 = -t243 * t158 - t355 * t159;
t85 = qJD(2) * t271 + t120;
t146 = t240 * t165;
t98 = qJD(2) * t267 + t146;
t24 = t115 * t310 - t122 * t311 + t244 * t85 + t247 * t98;
t263 = pkin(1) * t251 + t282;
t250 = qJD(2) ^ 2;
t262 = pkin(6) * t250 + qJDD(1) * t357 + t345;
t86 = t160 + t354;
t259 = -t248 * t282 - t344;
t25 = -t63 * qJD(4) - t244 * t98 + t247 * t85;
t255 = t256 + t354;
t233 = t249 * pkin(6);
t211 = t245 * t251 * t248;
t194 = pkin(4) * t227 + t353;
t177 = t239 * qJDD(1) - 0.2e1 * t284;
t173 = qJDD(5) + t187;
t152 = t227 * t246 + t228 * t322;
t150 = t227 * t249 - t228 * t324;
t145 = pkin(4) * t188 - t220;
t140 = -pkin(6) * t328 + t179;
t130 = -pkin(6) * t298 + t170;
t121 = -t241 * t302 + t146;
t119 = pkin(4) * t158 + t193;
t112 = t355 * t188 + t189 * t243;
t92 = qJD(2) * t266 + t245 * t361;
t91 = -t247 * t248 * t308 + t169 * t245 + t244 * t297;
t83 = t355 * t158 - t159 * t243;
t69 = pkin(4) * t92 + t175;
t32 = t86 + t356;
t31 = t84 * qJD(5) - t243 * t91 + t355 * t92;
t30 = t158 * t294 - t159 * t309 + t243 * t92 + t355 * t91;
t21 = -pkin(8) * t92 + t24;
t20 = pkin(4) * t313 + pkin(8) * t91 + t25;
t13 = t355 * t28 - t335;
t12 = -t243 * t28 - t301;
t10 = t355 * t26 - t335;
t4 = -t23 * qJD(5) + t355 * t20 - t243 * t21;
t3 = t22 * qJD(5) + t243 * t20 + t355 * t21;
t5 = [0, 0, 0, 0, 0, qJDD(1), t281, t282, 0, 0, t238 * qJDD(1) + 0.2e1 * t284, -0.2e1 * qJD(2) * t285 + 0.2e1 * t229 * t245, qJDD(2) * t245 + t248 * t250, t177, qJDD(2) * t248 - t245 * t250, 0, t269 * t245 + (-t262 + t349) * t248, t245 * t262 + t248 * t269 - t216, 0.2e1 * t317 * t334 - t282, -g(1) * (-pkin(1) * t246 + t233) - g(2) * t319 + (pkin(6) ^ 2 * t317 + pkin(1) ^ 2) * qJDD(1), (t245 * t135 + t182 * t312) * t241, (-t331 - t332) * t245 - t278 * t312, (-t135 - t289) * t248 + (t182 * t245 + t241 * t285) * qJD(2), (t245 * t134 + t180 * t312) * t240, (t134 + t291) * t248 + (-t180 * t245 - t240 * t285) * qJD(2), t177, -t282 * t240 + (pkin(6) * t134 + t160 * t240 + (qJD(1) * t140 + t116) * qJD(2)) * t245 + (-t120 * qJD(1) - t140 * qJDD(1) - t64 + t281 * t241 + (t195 * t240 + t352) * qJD(2)) * t248, -t282 * t241 + (pkin(6) * t135 + t160 * t241 + (-qJD(1) * t141 - t117) * qJD(2)) * t245 + (t121 * qJD(1) + t141 * qJDD(1) + t65 - t281 * t240 + (t195 * t241 + t351) * qJD(2)) * t248, -t120 * t182 - t121 * t180 - t134 * t141 - t135 * t140 + t216 + (-t116 * t241 - t117 * t240) * t312 + (-t240 * t65 - t241 * t64 - t345) * t245, t65 * t141 + t117 * t121 + t64 * t140 + t116 * t120 - g(1) * t233 - g(2) * (t249 * t280 + t319) - t196 * t349 + (t160 * t245 + t195 * t312) * pkin(6), t107 * t91 + t159 * t48, t107 * t92 + t108 * t91 + t158 * t48 + t159 * t49, -t107 * t313 - t159 * t187 + t215 * t91 + t248 * t48, t108 * t92 + t158 * t49, -t108 * t313 - t158 * t187 + t215 * t92 + t248 * t49, -t187 * t248 - t215 * t313, -g(1) * t150 - g(2) * t152 + t108 * t175 + t128 * t92 + t158 * t86 + t187 * t62 + t193 * t49 - t215 * t25 - t248 * t9 + t313 * t34, -g(1) * t149 - g(2) * t151 - t107 * t175 - t128 * t91 - t159 * t86 - t187 * t63 - t193 * t48 + t215 * t24 + t248 * t8 - t313 * t35, t107 * t25 - t108 * t24 - t158 * t8 + t159 * t9 + t34 * t91 - t35 * t92 + t48 * t62 - t49 * t63 + t283, t8 * t63 + t35 * t24 + t9 * t62 + t34 * t25 + t86 * t193 + t128 * t175 - g(1) * (t249 * t353 + t233) - g(2) * (t220 * t322 + t325 * t342 + t319) + (-g(1) * (-pkin(1) - t275) - g(2) * t353) * t246, -t14 * t84 + t30 * t54, t14 * t83 - t15 * t84 - t30 * t58 + t31 * t54, t14 * t248 + t173 * t84 + t208 * t30 - t313 * t54, t15 * t83 - t31 * t58, t15 * t248 - t173 * t83 + t208 * t31 + t313 * t58, -t173 * t248 - t208 * t313, -g(1) * t137 - g(2) * t139 + t10 * t313 + t119 * t15 + t173 * t22 - t2 * t248 - t208 * t4 + t31 * t68 + t32 * t83 - t58 * t69, -g(1) * t136 - g(2) * t138 + t1 * t248 - t11 * t313 - t119 * t14 - t173 * t23 + t208 * t3 - t30 * t68 + t32 * t84 - t54 * t69, -t1 * t83 + t10 * t30 - t11 * t31 + t14 * t22 - t15 * t23 - t2 * t84 + t3 * t58 + t4 * t54 + t283, t1 * t23 + t11 * t3 + t2 * t22 + t10 * t4 + t32 * t119 + t68 * t69 - g(1) * (t194 * t249 + t233) - g(2) * (t191 * t322 - t236 * t325 + t319) + (-g(1) * (-pkin(1) - t277) - g(2) * t194) * t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, t318 * t251, t304, t211, t229, qJDD(2), t245 * t263 - t221 - t235, t344 + (t263 - t334) * t248, 0, 0, -t182 * t241 * t315 + t331, -t134 * t240 + t135 * t241 + t278 * t315, t241 * t330 + t245 * t272 - t290, -t180 * t299 - t332, -t240 * t330 + t245 * t273 - t288, t211, qJ(3) * t290 - pkin(2) * t134 - t256 * t241 + ((-qJ(3) * t314 - t116) * t245 + (t240 * t307 + t129 - t352) * t248) * qJD(1), qJ(3) * t288 - pkin(2) * t135 + t256 * t240 + ((-qJ(3) * t308 + t117) * t245 + (t241 * t307 - t130 - t351) * t248) * qJD(1), t129 * t182 + t130 * t180 + (-qJ(3) * t134 - qJD(3) * t180 + t116 * t315 + t65) * t241 + (qJ(3) * t135 + qJD(3) * t182 + t117 * t315 - t64) * t240 + t259, -t195 * t223 - t116 * t129 - t117 * t130 + (-t116 * t240 + t117 * t241) * qJD(3) - t256 * pkin(2) + (-t64 * t240 + t65 * t241 + t259) * qJ(3), -t107 * t366 - t189 * t48, -t107 * t321 - t108 * t366 + t188 * t48 - t189 * t49, t107 * t316 + t187 * t189 - t215 * t366, -t108 * t321 + t188 * t49, t108 * t316 - t187 * t188 - t215 * t321, t215 * t316, -t108 * t174 + t123 * t187 - t128 * t321 + t188 * t86 - t215 * t336 - t220 * t49 + t228 * t261 - t316 * t34, t107 * t174 - t124 * t187 + t128 * t366 + t189 * t86 + t215 * t337 + t220 * t48 - t227 * t261 + t316 * t35, t107 * t336 - t108 * t337 + t123 * t48 - t124 * t49 - t188 * t8 - t189 * t9 + t321 * t35 - t34 * t366 + t259, -g(3) * t275 + t9 * t123 + t8 * t124 - t128 * t174 - t86 * t220 + t336 * t34 + t337 * t35 + t282 * (t220 * t245 - t248 * t342), -t113 * t14 + t339 * t54, t112 * t14 - t113 * t15 + t338 * t54 - t339 * t58, t113 * t173 + t208 * t339 + t316 * t54, t112 * t15 - t338 * t58, -t112 * t173 + t208 * t338 - t316 * t58, t208 * t316, -t10 * t316 + t112 * t32 + t145 * t15 + t173 * t45 - t208 * t340 + t219 * t261 - t287 * t58 + t338 * t68, t11 * t316 + t113 * t32 - t14 * t145 - t173 * t46 + t208 * t341 - t218 * t261 - t287 * t54 - t339 * t68, -t1 * t112 + t10 * t339 - t11 * t338 - t113 * t2 + t14 * t45 - t15 * t46 + t340 * t54 + t341 * t58 + t259, -g(3) * t277 + t1 * t46 + t10 * t340 + t11 * t341 + t32 * t145 + t2 * t45 + t287 * t68 + t282 * (t191 * t245 + t236 * t248); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248 * t272 - t225 + t291, t248 * t273 + t289 + t305, -t180 ^ 2 - t182 ^ 2, t116 * t182 + t117 * t180 + t256, 0, 0, 0, 0, 0, 0, t49 + t368, -t48 + t367, -t371 - t372, -t107 * t34 + t108 * t35 + t255, 0, 0, 0, 0, 0, 0, t15 + t373, -t14 - t374, -t375 - t376, -t10 * t54 - t11 * t58 + t255 + t356; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t333, -t371 + t372, -t48 - t367, t333, -t49 + t368, t187, t107 * t128 - t215 * t35 + t360 + t9, g(1) * t152 - g(2) * t150 + t108 * t128 - t215 * t34 + t228 * t344 - t8, 0, 0, t343, t365, t364, -t343, t358, t173, t12 * t208 + (-t107 * t58 + t355 * t173 + t208 * t309) * pkin(4) + t359, -t13 * t208 + (-t107 * t54 - t173 * t243 + t208 * t294) * pkin(4) + t363, t10 * t58 - t11 * t54 - t12 * t54 - t13 * t58 + (t355 * t14 - t15 * t243 + (-t243 * t54 + t355 * t58) * qJD(5)) * pkin(4), -t10 * t12 - t11 * t13 + (t1 * t243 + t2 * t355 + t68 * t107 + (-t10 * t243 + t11 * t355) * qJD(5) + t360) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t343, t365, t364, -t343, t358, t173, -t11 * t208 + t359, -t10 * t208 + t363, 0, 0;];
tau_reg = t5;

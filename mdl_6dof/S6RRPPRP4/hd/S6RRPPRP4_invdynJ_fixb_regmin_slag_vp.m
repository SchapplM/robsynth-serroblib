% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:45
% EndTime: 2019-03-09 08:39:59
% DurationCPUTime: 6.26s
% Computational Cost: add. (5223->584), mult. (11694->713), div. (0->0), fcn. (8089->8), ass. (0->261)
t216 = sin(pkin(9));
t217 = cos(pkin(9));
t223 = cos(qJ(2));
t307 = qJD(1) * qJD(2);
t290 = t223 * t307;
t220 = sin(qJ(2));
t305 = t220 * qJDD(1);
t246 = t290 + t305;
t108 = qJDD(2) * t216 + t217 * t246;
t219 = sin(qJ(5));
t222 = cos(qJ(5));
t199 = t217 * qJDD(2);
t236 = t216 * t246 - t199;
t317 = qJD(1) * t220;
t297 = t216 * t317;
t308 = t217 * qJD(2);
t147 = -t297 + t308;
t295 = t217 * t317;
t315 = qJD(2) * t216;
t148 = t295 + t315;
t73 = t147 * t219 - t148 * t222;
t233 = qJD(5) * t73 - t108 * t219 + t222 * t236;
t316 = qJD(1) * t223;
t184 = qJD(5) + t316;
t342 = t184 * t73;
t375 = t233 - t342;
t361 = t73 ^ 2;
t151 = t216 * t219 + t217 * t222;
t123 = t151 * t223;
t132 = t151 * qJD(5);
t324 = -qJD(1) * t123 - t132;
t294 = t217 * t316;
t296 = t216 * t316;
t310 = qJD(5) * t222;
t311 = qJD(5) * t219;
t323 = t216 * t310 - t217 * t311 - t219 * t294 + t222 * t296;
t221 = sin(qJ(1));
t224 = cos(qJ(1));
t273 = g(1) * t224 + g(2) * t221;
t252 = t273 * t220;
t194 = pkin(7) * t305;
t125 = -qJDD(2) * pkin(2) + pkin(7) * t290 + qJDD(3) + t194;
t212 = g(3) * t223;
t287 = -t125 - t212;
t230 = -t252 - t287;
t288 = -t216 * qJ(4) - pkin(2);
t374 = t217 * pkin(3) - t288;
t202 = t223 * qJDD(1);
t289 = t220 * t307;
t245 = -t289 + t202;
t150 = -qJDD(5) - t245;
t360 = t184 ^ 2;
t373 = -t148 * t73 - t150 * t219 + t222 * t360;
t372 = -2 * pkin(1);
t260 = t222 * t147 + t148 * t219;
t371 = t260 ^ 2;
t22 = -t222 * t108 + t147 * t310 + t148 * t311 - t219 * t236;
t343 = t184 * t260;
t370 = t22 - t343;
t140 = t148 * qJD(4);
t364 = pkin(3) * t236 - t108 * qJ(4) - t140;
t25 = t125 + t364;
t292 = -t25 - t212;
t336 = t216 * t222;
t152 = -t217 * t219 + t336;
t299 = -pkin(7) * t216 - pkin(3);
t333 = t217 * t223;
t239 = -pkin(8) * t333 + (-pkin(4) + t299) * t220;
t269 = pkin(2) * t220 - qJ(3) * t223;
t155 = t269 * qJD(1);
t338 = t155 * t217;
t49 = qJD(1) * t239 - t338;
t138 = t216 * t155;
t192 = qJ(4) * t317;
t334 = t217 * t220;
t335 = t216 * t223;
t249 = -pkin(7) * t334 + pkin(8) * t335;
t61 = qJD(1) * t249 + t138 + t192;
t349 = pkin(8) - qJ(3);
t161 = t349 * t216;
t162 = t349 * t217;
t90 = -t161 * t219 - t162 * t222;
t369 = qJD(3) * t152 - qJD(5) * t90 + t219 * t61 - t222 * t49;
t259 = -t161 * t222 + t162 * t219;
t368 = -qJD(3) * t151 - qJD(5) * t259 + t219 * t49 + t222 * t61;
t309 = t216 * qJD(4);
t197 = pkin(7) * t316;
t322 = qJ(4) * t294 - t197;
t359 = -pkin(3) - pkin(4);
t284 = -t296 * t359 + t309 - t322;
t367 = t151 * t220;
t248 = t220 * t152;
t145 = t148 ^ 2;
t366 = -t147 ^ 2 - t145;
t314 = qJD(2) * t220;
t365 = qJ(4) * t314 - qJD(4) * t223;
t363 = -pkin(5) * t233 + t22 * qJ(6) + qJD(6) * t73;
t203 = t220 * qJ(3);
t208 = t223 * pkin(2);
t319 = t208 + t203;
t158 = -pkin(1) - t319;
t177 = pkin(7) * t335;
t207 = t223 * pkin(3);
t69 = pkin(4) * t223 + t177 + t207 + (-pkin(8) * t220 - t158) * t217;
t337 = t216 * t220;
t110 = pkin(7) * t333 + t216 * t158;
t92 = -qJ(4) * t223 + t110;
t79 = pkin(8) * t337 + t92;
t262 = t219 * t69 + t222 * t79;
t130 = qJD(2) * t269 - qJD(3) * t220;
t339 = t130 * t217;
t43 = qJD(2) * t239 - t339;
t113 = t216 * t130;
t44 = qJD(2) * t249 + t113 + t365;
t362 = -qJD(5) * t262 - t219 * t44 + t222 * t43;
t226 = qJD(1) ^ 2;
t357 = pkin(1) * t226;
t356 = pkin(5) * t150;
t355 = pkin(7) * t148;
t354 = g(1) * t221;
t300 = -pkin(1) - t208;
t255 = t300 - t203;
t142 = t255 * qJD(1);
t165 = qJD(2) * qJ(3) + t197;
t80 = t142 * t217 - t216 * t165;
t59 = pkin(3) * t316 + qJD(4) - t80;
t36 = pkin(4) * t316 - pkin(8) * t148 + t59;
t81 = t216 * t142 + t217 * t165;
t68 = -qJ(4) * t316 + t81;
t40 = -pkin(8) * t147 + t68;
t12 = t219 * t36 + t222 * t40;
t9 = qJ(6) * t184 + t12;
t351 = t184 * t9;
t350 = t73 * t260;
t348 = t323 * pkin(5) - t324 * qJ(6) - qJD(6) * t152 + t284;
t347 = qJ(6) * t317 - t368;
t346 = -pkin(5) * t317 - t369;
t218 = qJD(2) * pkin(2);
t344 = t12 * t184;
t114 = pkin(7) * t245 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t71 = qJD(1) * t130 + qJDD(1) * t255;
t35 = t217 * t114 + t216 * t71;
t340 = qJ(6) * t150;
t332 = t220 * t221;
t331 = t220 * t224;
t330 = t221 * t217;
t329 = t221 * t223;
t328 = t223 * t224;
t327 = t224 * t216;
t11 = -t219 * t40 + t222 * t36;
t325 = qJD(6) - t11;
t293 = t223 * t308;
t321 = -qJ(4) * t293 - qJD(4) * t334;
t320 = g(1) * t332 - g(2) * t331;
t213 = t220 ^ 2;
t214 = t223 ^ 2;
t318 = t213 - t214;
t313 = qJD(2) * t223;
t306 = qJDD(1) * qJ(4);
t304 = pkin(7) * t314;
t303 = qJ(4) * t289 + t35;
t302 = t222 * t335;
t301 = g(1) * t328 + g(2) * t329 + g(3) * t220;
t298 = pkin(3) * t216 + pkin(7);
t157 = pkin(7) * t317 + qJD(3) - t218;
t291 = qJ(3) * t202;
t34 = -t216 * t114 + t217 * t71;
t253 = pkin(3) * t202 + qJDD(4) - t34;
t28 = -pkin(3) * t289 + t253;
t15 = pkin(4) * t245 - pkin(8) * t108 + t28;
t266 = t216 * t305 - t199;
t18 = t266 * pkin(8) + (-t306 + (pkin(8) * t315 - qJD(4)) * qJD(1)) * t223 + t303;
t286 = -t222 * t15 + t219 * t18 + t40 * t310 + t36 * t311;
t109 = t158 * t217 - t177;
t139 = t217 * pkin(4) + t374;
t283 = g(1) * t217 * t331 + g(2) * t220 * t330 + qJD(3) * t296 + t216 * t291;
t282 = pkin(3) * t333 + qJ(4) * t335 + t319;
t281 = t224 * pkin(1) + pkin(2) * t328 + t221 * pkin(7) + qJ(3) * t331;
t280 = t217 * t291;
t279 = t216 * t359 - pkin(7);
t134 = t216 * t329 + t217 * t224;
t135 = t217 * t329 - t327;
t64 = t134 * t222 - t135 * t219;
t136 = t223 * t327 - t330;
t137 = t216 * t221 + t217 * t328;
t66 = -t136 * t222 + t137 * t219;
t278 = g(1) * t64 + g(2) * t66;
t63 = t134 * t219 + t135 * t222;
t67 = t136 * t219 + t137 * t222;
t277 = g(1) * t63 - g(2) * t67;
t276 = t220 * t299;
t275 = -g(1) * t134 + g(2) * t136;
t274 = g(1) * t135 - g(2) * t137;
t272 = -g(2) * t224 + t354;
t209 = t224 * pkin(7);
t271 = -t135 * pkin(3) - qJ(4) * t134 + t209;
t270 = -t301 + (-qJ(3) * t236 + qJD(3) * t147) * t217;
t120 = t219 * t334 - t220 * t336;
t268 = -pkin(5) * t120 + qJ(6) * t367;
t263 = -t219 * t79 + t222 * t69;
t54 = -t147 * pkin(3) - t148 * qJ(4) + t157;
t258 = qJ(3) * t108 + qJD(3) * t148;
t225 = qJD(2) ^ 2;
t256 = qJDD(2) * t223 - t220 * t225;
t98 = -pkin(7) * t295 + t138;
t87 = -t217 * t304 + t113;
t251 = t219 * t15 + t222 * t18 + t36 * t310 - t311 * t40;
t250 = t219 * t43 + t222 * t44 + t69 * t310 - t311 * t79;
t174 = qJ(4) * t334;
t91 = t220 * t279 + t174;
t39 = pkin(4) * t147 - t54;
t243 = t137 * pkin(3) + qJ(4) * t136 + t281;
t242 = t255 * t354;
t241 = -g(1) * t136 - g(2) * t134 - g(3) * t337;
t240 = t22 + t343;
t104 = t221 * t248;
t106 = t224 * t248;
t122 = t219 * t333 - t302;
t238 = -g(1) * t106 - g(2) * t104 - g(3) * t122 - t90 * t150;
t105 = t221 * t367;
t107 = t224 * t367;
t237 = g(1) * t107 + g(2) * t105 - g(3) * t123 - t150 * t259;
t58 = t279 * t313 - t321;
t234 = g(1) * t66 - g(2) * t64 + g(3) * t120 - t286;
t232 = (-t148 + t315) * t316 + t266;
t231 = -t148 * t260 - t150 * t222 - t219 * t360;
t229 = -g(1) * t67 - g(2) * t63 - g(3) * t367 + t251;
t228 = t233 + t342;
t10 = pkin(5) * t260 + qJ(6) * t73 + t39;
t227 = -t10 * t73 + qJDD(6) - t234;
t19 = -pkin(4) * t236 - t25;
t182 = qJ(3) * t328;
t178 = qJ(3) * t329;
t126 = t147 * t316;
t115 = t220 * t298 - t174;
t103 = pkin(3) * t296 - t322;
t97 = pkin(7) * t297 + t338;
t96 = -t109 + t207;
t86 = t216 * t304 + t339;
t85 = qJD(1) * t276 - t338;
t84 = t192 + t98;
t83 = t298 * t313 + t321;
t72 = qJD(2) * t276 - t339;
t57 = t87 + t365;
t55 = -t126 + t108;
t53 = qJD(2) * t123 + qJD(5) * t248;
t52 = -qJD(2) * t302 + t132 * t220 + t219 * t293;
t48 = pkin(5) * t151 - qJ(6) * t152 + t139;
t31 = -t268 + t91;
t30 = -pkin(5) * t73 + qJ(6) * t260;
t27 = -pkin(5) * t223 - t263;
t26 = qJ(6) * t223 + t262;
t24 = (-qJD(1) * qJD(4) - t306) * t223 + t303;
t7 = -pkin(5) * t184 + t325;
t6 = pkin(5) * t52 - qJ(6) * t53 - qJD(6) * t367 + t58;
t5 = pkin(5) * t314 - t362;
t4 = -qJ(6) * t314 + qJD(6) * t223 + t250;
t3 = t19 + t363;
t2 = qJDD(6) + t286 + t356;
t1 = qJD(6) * t184 + t251 - t340;
t8 = [qJDD(1), t272, t273, qJDD(1) * t213 + 0.2e1 * t223 * t289, 0.2e1 * t202 * t220 - 0.2e1 * t307 * t318, qJDD(2) * t220 + t223 * t225, t256, 0 (-pkin(7) * qJDD(2) + t307 * t372) * t220 + (0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t225 + t272) * t223, -pkin(7) * t256 + t246 * t372 - t320 (t125 * t216 + (qJD(1) * t109 + t80) * qJD(2) + t266 * pkin(7)) * t220 + (-t86 * qJD(1) - t109 * qJDD(1) - t34 + (t157 * t216 + (-t147 + t297) * pkin(7)) * qJD(2)) * t223 + t274 (pkin(7) * t108 + t125 * t217 + (-qJD(1) * t110 - t81) * qJD(2)) * t220 + (qJD(1) * t87 + qJDD(1) * t110 + t35 + (t157 * t217 + t355) * qJD(2)) * t223 + t275, -t109 * t108 + t110 * t199 + t87 * t147 - t86 * t148 + (-t220 * t34 - t313 * t80) * t217 + (-t110 * t246 - t35 * t220 - t313 * t81) * t216 + t320, t35 * t110 + t81 * t87 + t34 * t109 + t80 * t86 - g(1) * t209 - g(2) * t281 - t242 + (t125 * t220 + t157 * t313) * pkin(7), -t115 * t199 - t83 * t147 + ((qJDD(1) * t115 + t25) * t216 + (-qJD(1) * t96 - t59) * qJD(2)) * t220 + (t54 * t315 + t96 * qJDD(1) + t28 + (t115 * t315 + t72) * qJD(1)) * t223 + t274, t96 * t108 + t57 * t147 + t72 * t148 + t92 * t199 + (t220 * t28 + t313 * t59) * t217 + (-t24 * t220 - t246 * t92 - t313 * t68) * t216 + t320, -t108 * t115 - t148 * t83 + (-t25 * t217 + (qJD(1) * t92 + t68) * qJD(2)) * t220 + (-qJD(1) * t57 - qJDD(1) * t92 - t308 * t54 - t24) * t223 - t275, -g(1) * t271 - g(2) * t243 + t25 * t115 + t24 * t92 + t28 * t96 + t54 * t83 + t68 * t57 + t59 * t72 - t242, -t22 * t367 - t53 * t73, t120 * t22 + t233 * t367 - t260 * t53 + t52 * t73, -t150 * t367 + t184 * t53 - t22 * t223 + t314 * t73, t120 * t150 - t184 * t52 + t223 * t233 + t260 * t314, -t150 * t223 - t184 * t314, -t11 * t314 + t19 * t120 - t263 * t150 + t184 * t362 - t286 * t223 - t233 * t91 + t58 * t260 + t39 * t52 + t277, t12 * t314 + t150 * t262 - t184 * t250 + t19 * t367 - t91 * t22 - t223 * t251 + t39 * t53 - t58 * t73 + t278, t10 * t52 + t120 * t3 + t150 * t27 - t184 * t5 - t2 * t223 - t233 * t31 + t260 * t6 + t314 * t7 + t277, -t1 * t120 + t2 * t367 - t22 * t27 + t233 * t26 - t260 * t4 - t5 * t73 - t52 * t9 + t53 * t7 - t320, t1 * t223 - t10 * t53 - t150 * t26 + t184 * t4 + t22 * t31 - t3 * t367 - t314 * t9 + t6 * t73 - t278, t1 * t26 + t9 * t4 + t3 * t31 + t10 * t6 + t2 * t27 + t7 * t5 - g(1) * (-pkin(4) * t135 - pkin(5) * t63 + qJ(6) * t64 + t271) - g(2) * (pkin(4) * t137 + pkin(5) * t67 - pkin(8) * t331 + qJ(6) * t66 + t243) - (t220 * t349 + t300) * t354; 0, 0, 0, -t220 * t226 * t223, t318 * t226, t305, t202, qJDD(2), -t194 - t212 + (t273 + t357) * t220 (-pkin(7) * qJDD(1) + t357) * t223 + t301, -pkin(2) * t266 + t287 * t217 + ((-qJ(3) * t315 - t80) * t220 + (pkin(7) * t147 + t97 + (-t157 - t218) * t216) * t223) * qJD(1) + t283, t280 - pkin(2) * t108 + t230 * t216 + ((-qJ(3) * t308 + t81) * t220 + (-t355 - t98 + (qJD(3) - t157) * t217) * t223) * qJD(1), -t147 * t98 + t148 * t97 + (t316 * t80 + t35) * t217 + (t316 * t81 + t258 - t34) * t216 + t270, -t125 * pkin(2) - t81 * t98 - t80 * t97 - t157 * t197 - g(1) * (-pkin(2) * t331 + t182) - g(2) * (-pkin(2) * t332 + t178) - g(3) * t319 + (-t80 * t216 + t81 * t217) * qJD(3) + (-t216 * t34 + t217 * t35) * qJ(3), -t374 * t266 + t292 * t217 + (t103 + t309) * t147 + (t59 * t220 - t85 * t223 + (-t223 * t54 + (-t223 * t374 - t203) * qJD(2)) * t216) * qJD(1) + t283, -t147 * t84 - t148 * t85 + (-t316 * t59 + t24) * t217 + (t316 * t68 + t258 + t28) * t216 + t270, -t280 + t103 * t148 + t108 * t374 + (t140 + t292 + t252) * t216 + (-t220 * t68 + t223 * t84 + (qJ(3) * t314 + (-qJD(3) + t54) * t223) * t217) * qJD(1), -t68 * t84 - t54 * t103 - t59 * t85 - g(1) * t182 - g(2) * t178 - g(3) * t282 + (t24 * qJ(3) + t68 * qJD(3)) * t217 + (t28 * qJ(3) + t59 * qJD(3) - t54 * qJD(4)) * t216 + (-t25 + t252) * t374, -t152 * t22 - t324 * t73, t151 * t22 + t152 * t233 - t260 * t324 + t323 * t73, -t150 * t152 + t184 * t324 - t317 * t73, t150 * t151 - t184 * t323 - t260 * t317, t184 * t317, t11 * t317 - t139 * t233 + t19 * t151 + t184 * t369 + t284 * t260 + t323 * t39 + t237, -t12 * t317 - t139 * t22 + t19 * t152 + t184 * t368 - t284 * t73 + t324 * t39 - t238, t10 * t323 + t151 * t3 - t184 * t346 - t233 * t48 + t260 * t348 - t317 * t7 + t237, -t1 * t151 + t152 * t2 + t22 * t259 + t233 * t90 - t260 * t347 - t323 * t9 + t324 * t7 - t346 * t73 + t301, -t10 * t324 - t152 * t3 + t184 * t347 + t22 * t48 + t317 * t9 + t348 * t73 + t238, t1 * t90 + t3 * t48 - t2 * t259 - g(1) * (-pkin(5) * t107 - pkin(8) * t328 + qJ(6) * t106 + t182) - g(2) * (-pkin(5) * t105 - pkin(8) * t329 + qJ(6) * t104 + t178) - g(3) * (pkin(4) * t333 + pkin(5) * t123 + qJ(6) * t122 + t282) + t347 * t9 + t346 * t7 + t348 * t10 + (g(3) * pkin(8) + t273 * (-t217 * t359 - t288)) * t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, t55, t366, -t147 * t81 + t148 * t80 + t230, t232, t366, -t55, -t147 * t68 - t148 * t59 + t230 + t364, 0, 0, 0, 0, 0, t228, t240, t228, t371 + t361, -t240 (t216 * t290 - t199) * pkin(4) - t9 * t260 - t7 * t73 + (qJDD(1) * t216 * pkin(4) - t273) * t220 - t363 - t292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147 * t148 + t245, t126 + t108, -t214 * t226 - t145, t148 * t54 + (-pkin(3) * t314 + t223 * t68) * qJD(1) + t241 + t253, 0, 0, 0, 0, 0, t231, -t373, t231, t219 * t375 + t370 * t222, t373, -t10 * t148 + (-t2 + t351) * t222 + (t184 * t7 + t1) * t219 + t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t350, t361 - t371, -t370, t375, -t150, t39 * t73 + t234 + t344, t11 * t184 + t260 * t39 - t229, -t260 * t30 - t227 + t344 - 0.2e1 * t356, pkin(5) * t22 + qJ(6) * t233 - (-t12 + t9) * t73 + (t7 - t325) * t260, -0.2e1 * t340 - t10 * t260 - t30 * t73 + (0.2e1 * qJD(6) - t11) * t184 + t229, t1 * qJ(6) - t2 * pkin(5) - t10 * t30 - t7 * t12 - g(1) * (-pkin(5) * t66 + qJ(6) * t67) - g(2) * (pkin(5) * t64 + qJ(6) * t63) - g(3) * t268 + t325 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150 - t350, -t370, -t360 - t361, t227 - t351 + t356;];
tau_reg  = t8;

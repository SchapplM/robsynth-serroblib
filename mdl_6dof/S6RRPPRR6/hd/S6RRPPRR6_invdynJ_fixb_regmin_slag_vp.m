% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:24
% EndTime: 2019-03-09 09:15:33
% DurationCPUTime: 5.19s
% Computational Cost: add. (5175->470), mult. (11526->593), div. (0->0), fcn. (8647->12), ass. (0->252)
t208 = qJDD(2) - qJDD(5);
t209 = qJD(2) - qJD(5);
t219 = sin(qJ(6));
t223 = cos(qJ(6));
t216 = sin(pkin(10));
t217 = cos(pkin(10));
t225 = cos(qJ(2));
t302 = qJD(1) * t225;
t221 = sin(qJ(2));
t303 = qJD(1) * t221;
t122 = t216 * t303 + t217 * t302;
t285 = t216 * t302;
t125 = t217 * t303 - t285;
t220 = sin(qJ(5));
t224 = cos(qJ(5));
t298 = qJD(5) * t224;
t299 = qJD(5) * t220;
t141 = t216 * t221 + t217 * t225;
t294 = qJD(1) * qJD(2);
t284 = t221 * t294;
t236 = qJDD(1) * t141 - t217 * t284;
t283 = t225 * t294;
t80 = t216 * t283 + t236;
t131 = t141 * qJD(2);
t196 = t221 * qJDD(1);
t293 = t225 * qJDD(1);
t256 = t217 * t196 - t216 * t293;
t81 = qJD(1) * t131 + t256;
t243 = t122 * t298 + t125 * t299 + t220 * t80 - t224 * t81;
t296 = qJD(6) * t223;
t297 = qJD(6) * t219;
t68 = t122 * t220 - t125 * t224;
t10 = -t219 * t208 - t209 * t296 - t223 * t243 + t297 * t68;
t252 = t209 * t219 + t223 * t68;
t11 = -t252 * qJD(6) + t223 * t208 - t219 * t243;
t341 = t224 * t122 + t125 * t220;
t350 = qJD(6) + t341;
t327 = t252 * t350;
t47 = t223 * t209 - t219 * t68;
t328 = t47 * t350;
t369 = (t10 - t328) * t223 + (-t11 + t327) * t219;
t26 = -qJD(5) * t68 + t220 * t81 + t224 * t80;
t24 = qJDD(6) + t26;
t21 = t223 * t24;
t275 = t350 ^ 2;
t366 = -t219 * t275 + t21;
t323 = t10 * t219;
t362 = t223 * t350;
t365 = t252 * t362 - t323;
t322 = t219 * t24;
t355 = t252 * t68;
t364 = t350 * t362 + t322 - t355;
t357 = t47 * t68;
t361 = t357 - t366;
t331 = pkin(8) * t125;
t191 = pkin(7) * t303;
t148 = qJ(4) * t303 - t191;
t334 = pkin(2) + pkin(3);
t286 = t334 * qJD(2);
t106 = qJD(3) - t286 - t148;
t192 = pkin(7) * t302;
t150 = -qJ(4) * t302 + t192;
t212 = qJD(2) * qJ(3);
t136 = t150 + t212;
t57 = t217 * t106 - t136 * t216;
t42 = -qJD(2) * pkin(4) - t331 + t57;
t332 = pkin(8) * t122;
t58 = t216 * t106 + t217 * t136;
t43 = t58 - t332;
t15 = t220 * t42 + t224 * t43;
t13 = -pkin(9) * t209 + t15;
t137 = -qJD(1) * pkin(1) - pkin(2) * t302 - qJ(3) * t303;
t100 = pkin(3) * t302 + qJD(4) - t137;
t67 = pkin(4) * t122 + t100;
t22 = pkin(5) * t341 + pkin(9) * t68 + t67;
t5 = -t13 * t219 + t22 * t223;
t360 = t5 * t68;
t6 = t13 * t223 + t219 * t22;
t359 = t6 * t68;
t358 = t68 * pkin(5);
t356 = t350 * t68;
t354 = t68 * t209;
t353 = t68 * t341;
t352 = t209 * t341;
t351 = t283 + t196;
t349 = t341 ^ 2 - t68 ^ 2;
t347 = -t26 + t354;
t292 = pkin(10) + qJ(5);
t195 = sin(t292);
t272 = cos(t292);
t238 = t221 * t195 + t225 * t272;
t222 = sin(qJ(1));
t254 = t221 * t272;
t339 = -t225 * t195 + t254;
t96 = t339 * t222;
t226 = cos(qJ(1));
t312 = t225 * t226;
t98 = t195 * t312 - t226 * t254;
t242 = g(1) * t98 - g(2) * t96 + g(3) * t238;
t187 = pkin(7) * t196;
t282 = pkin(7) * t283 + qJDD(3) + t187;
t295 = t221 * qJD(4);
t63 = -t351 * qJ(4) - qJD(1) * t295 - t334 * qJDD(2) + t282;
t301 = qJD(2) * t221;
t326 = pkin(7) - qJ(4);
t118 = -t225 * qJD(4) - t326 * t301;
t188 = pkin(7) * t293;
t210 = qJDD(2) * qJ(3);
t211 = qJD(2) * qJD(3);
t289 = t188 + t210 + t211;
t66 = -qJ(4) * t293 + qJD(1) * t118 + t289;
t277 = t216 * t66 - t217 * t63;
t19 = -qJDD(2) * pkin(4) - pkin(8) * t81 - t277;
t34 = t216 * t63 + t217 * t66;
t20 = -pkin(8) * t80 + t34;
t251 = -t224 * t19 + t220 * t20 + t43 * t298 + t42 * t299;
t346 = t67 * t68 + t242 - t251;
t244 = t220 * t19 + t224 * t20 + t42 * t298 - t43 * t299;
t97 = t238 * t222;
t261 = g(2) * t97 + g(3) * t339;
t99 = t238 * t226;
t345 = g(1) * t99 + t341 * t67 - t244 + t261;
t343 = t243 + t352;
t306 = t225 * pkin(2) + t221 * qJ(3);
t342 = -pkin(1) - t306;
t151 = -qJ(3) * t216 - t217 * t334;
t147 = -pkin(4) + t151;
t152 = qJ(3) * t217 - t216 * t334;
t308 = t220 * t147 + t224 * t152;
t207 = g(1) * t226;
t340 = g(2) * t222 + t207;
t2 = pkin(5) * t208 + t251;
t240 = -t2 + t242;
t185 = qJ(3) * t302;
t115 = -t334 * t303 + t185;
t79 = -t125 * pkin(4) + t115;
t85 = -pkin(9) + t308;
t338 = (-pkin(9) * t341 + qJD(6) * t85 + t358 + t79) * t350 + t240;
t337 = (t350 * pkin(9) - t358) * t350 - t240;
t321 = pkin(7) * qJDD(2);
t336 = (qJD(1) * t342 + t137) * qJD(2) - t321;
t14 = -t220 * t43 + t224 * t42;
t12 = pkin(5) * t209 - t14;
t279 = -pkin(9) * t208 + qJD(6) * t22 + t244;
t248 = t216 * t225 - t217 * t221;
t82 = t224 * t141 - t220 * t248;
t83 = -t141 * t220 - t224 * t248;
t200 = t225 * pkin(3);
t288 = t200 + t306;
t139 = pkin(1) + t288;
t90 = pkin(4) * t141 + t139;
t30 = pkin(5) * t82 - pkin(9) * t83 + t90;
t158 = t326 * t221;
t159 = t326 * t225;
t91 = t217 * t158 - t159 * t216;
t61 = pkin(8) * t248 + t91;
t92 = t216 * t158 + t217 * t159;
t62 = -pkin(8) * t141 + t92;
t32 = t220 * t61 + t224 * t62;
t300 = qJD(2) * t225;
t130 = t216 * t300 - t217 * t301;
t37 = -qJD(5) * t82 - t220 * t130 + t224 * t131;
t31 = t220 * t62 - t224 * t61;
t119 = qJD(2) * t159 - t295;
t59 = -t118 * t216 + t217 * t119;
t44 = -pkin(8) * t131 + t59;
t60 = t217 * t118 + t216 * t119;
t45 = -pkin(8) * t130 + t60;
t7 = -t31 * qJD(5) + t220 * t44 + t224 * t45;
t335 = t12 * t37 + t2 * t83 - t32 * t24 - (qJD(6) * t30 + t7) * t350 - t279 * t82 + t207;
t333 = g(1) * t97;
t206 = g(1) * t222;
t330 = g(2) * t226;
t329 = t12 * t83;
t140 = t216 * t220 - t217 * t224;
t250 = t147 * t224 - t152 * t220;
t86 = -t148 * t216 + t217 * t150;
t50 = t86 - t332;
t87 = t217 * t148 + t216 * t150;
t51 = t87 + t331;
t325 = -qJD(3) * t140 + qJD(5) * t250 - t220 * t50 - t224 * t51;
t143 = t216 * t224 + t217 * t220;
t324 = t143 * qJD(3) + t308 * qJD(5) - t220 * t51 + t224 * t50;
t320 = qJD(6) * t13;
t215 = qJDD(1) * pkin(1);
t319 = qJDD(2) * pkin(2);
t317 = t221 * t222;
t316 = t221 * t226;
t229 = qJD(1) ^ 2;
t315 = t221 * t229;
t314 = t222 * t225;
t310 = t209 * t140;
t309 = t209 * t143;
t197 = t221 * qJD(3);
t307 = qJ(3) * t300 + t197;
t213 = t221 ^ 2;
t214 = t225 ^ 2;
t304 = t213 - t214;
t291 = t83 * t297;
t290 = t225 * t315;
t287 = -g(1) * t316 - g(2) * t317 + g(3) * t225;
t280 = t206 - t330;
t269 = -qJD(2) * pkin(2) + qJD(3);
t268 = qJD(3) * t216 + t86;
t267 = qJD(3) * t217 - t87;
t265 = t226 * pkin(1) + pkin(2) * t312 + t222 * pkin(7) + qJ(3) * t316;
t264 = -t187 - t287;
t263 = t221 * t286;
t262 = t30 * t24 + t333;
t228 = qJD(2) ^ 2;
t260 = pkin(7) * t228 + t330;
t258 = t24 * t83 + t350 * t37;
t257 = -t320 - t330;
t255 = t216 * t57 - t217 * t58;
t253 = pkin(2) * t293 + t351 * qJ(3) + qJD(1) * t197 + t215;
t155 = t191 + t269;
t157 = t192 + t212;
t249 = t155 * t225 - t157 * t221;
t247 = qJD(6) * t143 + t303;
t114 = t282 - t319;
t245 = -0.2e1 * pkin(1) * t294 - t321;
t95 = -t263 + t307;
t241 = -t260 + 0.2e1 * t215;
t239 = t261 - t279;
t237 = -pkin(9) * t24 + (t12 + t14) * t350;
t65 = t130 * pkin(4) + t95;
t234 = -t85 * t24 + (-t12 - t325) * t350;
t120 = pkin(2) * t301 - t307;
t78 = pkin(2) * t284 - t253;
t233 = -qJD(1) * t120 - qJDD(1) * t342 - t260 - t78;
t105 = -pkin(7) * t284 + t289;
t231 = qJD(2) * t249 + t105 * t225 + t114 * t221;
t46 = pkin(3) * t293 - qJD(1) * t263 + qJDD(4) + t253;
t35 = t80 * pkin(4) + t46;
t202 = t226 * pkin(7);
t178 = g(1) * t314;
t172 = qJ(3) * t312;
t170 = qJ(3) * t314;
t149 = pkin(2) * t303 - t185;
t113 = t141 * t226;
t112 = t248 * t226;
t111 = t141 * t222;
t110 = t248 * t222;
t89 = -t219 * t222 + t223 * t99;
t88 = -t219 * t99 - t222 * t223;
t84 = pkin(5) - t250;
t38 = qJD(5) * t83 + t224 * t130 + t220 * t131;
t9 = t38 * pkin(5) - t37 * pkin(9) + t65;
t8 = t32 * qJD(5) + t220 * t45 - t224 * t44;
t4 = t26 * pkin(5) + pkin(9) * t243 + t35;
t3 = t223 * t4;
t1 = [qJDD(1), t280, t340, qJDD(1) * t213 + 0.2e1 * t221 * t283, 0.2e1 * t221 * t293 - 0.2e1 * t304 * t294, qJDD(2) * t221 + t225 * t228, qJDD(2) * t225 - t221 * t228, 0, t221 * t245 + t225 * t241 + t178, t245 * t225 + (-t241 - t206) * t221, t336 * t221 + t233 * t225 + t178 (t213 + t214) * qJDD(1) * pkin(7) + t231 - t340, -t336 * t225 + (t233 + t206) * t221, t231 * pkin(7) - g(1) * t202 - g(2) * t265 + t137 * t120 + (-t206 + t78) * t342, g(1) * t111 - g(2) * t113 - qJD(2) * t59 - qJDD(2) * t91 + t100 * t130 + t122 * t95 + t139 * t80 + t141 * t46, -g(1) * t110 + g(2) * t112 + qJD(2) * t60 + qJDD(2) * t92 + t100 * t131 + t125 * t95 + t139 * t81 - t248 * t46, -t122 * t60 - t125 * t59 - t130 * t58 - t131 * t57 - t141 * t34 - t248 * t277 - t80 * t92 - t81 * t91 + t340, t34 * t92 + t58 * t60 - t277 * t91 + t57 * t59 + t46 * t139 + t100 * t95 - g(1) * (-t226 * qJ(4) + t202) - g(2) * (pkin(3) * t312 + t265) + (-g(1) * (t342 - t200) + g(2) * qJ(4)) * t222, -t243 * t83 - t37 * t68, t243 * t82 - t26 * t83 - t341 * t37 + t38 * t68, -t208 * t83 - t209 * t37, t208 * t82 + t209 * t38, 0, -g(2) * t99 + t208 * t31 + t209 * t8 + t26 * t90 + t341 * t65 + t35 * t82 + t38 * t67 + t333, g(1) * t96 + g(2) * t98 + t208 * t32 + t209 * t7 - t243 * t90 + t35 * t83 + t37 * t67 - t65 * t68, t252 * t291 + (t10 * t83 - t252 * t37) * t223 (t219 * t252 - t223 * t47) * t37 + (-t323 - t11 * t223 + (t219 * t47 + t223 * t252) * qJD(6)) * t83, t10 * t82 + t223 * t258 - t252 * t38 - t291 * t350, -t296 * t350 * t83 - t11 * t82 - t219 * t258 - t47 * t38, t24 * t82 + t350 * t38, -g(2) * t89 + t31 * t11 + t3 * t82 + t5 * t38 + t8 * t47 + (t9 * t350 + (-t13 * t82 - t32 * t350 + t329) * qJD(6) + t262) * t223 + t335 * t219, -g(2) * t88 + t31 * t10 - t6 * t38 - t8 * t252 + (-(-qJD(6) * t32 + t9) * t350 - (t4 - t320) * t82 - qJD(6) * t329 - t262) * t219 + t335 * t223; 0, 0, 0, -t290, t304 * t229, t196, t293, qJDD(2), pkin(1) * t315 + t264, g(3) * t221 - t188 + (pkin(1) * t229 + t340) * t225, 0.2e1 * t319 - qJDD(3) + (-t137 * t221 + t149 * t225) * qJD(1) + t264 (-pkin(2) * t221 + qJ(3) * t225) * qJDD(1) + ((t157 - t212) * t221 + (-t155 + t269) * t225) * qJD(1), t188 + 0.2e1 * t210 + 0.2e1 * t211 + (qJD(1) * t149 - g(3)) * t221 + (qJD(1) * t137 - t340) * t225, t105 * qJ(3) + t157 * qJD(3) - t114 * pkin(2) - t137 * t149 - g(1) * (-pkin(2) * t316 + t172) - g(2) * (-pkin(2) * t317 + t170) - g(3) * t306 - t249 * qJD(1) * pkin(7), -g(1) * t112 - g(2) * t110 - g(3) * t141 + qJD(2) * t268 - t151 * qJDD(2) + t100 * t125 - t115 * t122 + t277, -g(1) * t113 - g(2) * t111 + g(3) * t248 + qJD(2) * t267 + t152 * qJDD(2) - t100 * t122 - t115 * t125 + t34, -t151 * t81 - t152 * t80 + (t268 - t58) * t125 + (-t267 + t57) * t122, t340 * t221 * t334 - g(1) * t172 - g(2) * t170 - g(3) * t288 - t255 * qJD(3) - t100 * t115 - t151 * t277 + t34 * t152 - t57 * t86 - t58 * t87, t353, t349, t343, -t347, t208, -t250 * t208 + t324 * t209 - t341 * t79 - t346, t308 * t208 + t325 * t209 + t68 * t79 - t345, t365, -t369, -t364, t361, -t356, t84 * t11 + t234 * t219 - t338 * t223 + t324 * t47 - t360, t84 * t10 + t338 * t219 + t234 * t223 - t252 * t324 + t359; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t290, t196, -t213 * t229 - t228, -qJD(2) * t157 + t137 * t303 + t114 + t287, -qJDD(2) * t217 - t122 * t303 - t216 * t228, qJDD(2) * t216 - t125 * t303 - t217 * t228, -t216 * t80 - t217 * t81 + (t122 * t217 - t125 * t216) * qJD(2), qJD(2) * t255 - t100 * t303 + t34 * t216 - t217 * t277 + t287, 0, 0, 0, 0, 0, t140 * t208 - t209 * t309 - t303 * t341, t143 * t208 + t209 * t310 + t303 * t68, 0, 0, 0, 0, 0, -t143 * t322 + t140 * t11 - t309 * t47 + (-t219 * t310 - t223 * t247) * t350, -t143 * t21 + t140 * t10 + t309 * t252 + (t219 * t247 - t223 * t310) * t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t125 + t285) * qJD(2) + t236 (qJD(1) * t141 + t122) * qJD(2) + t256, -t122 ^ 2 - t125 ^ 2, t58 * t122 + t57 * t125 + t280 + t46, 0, 0, 0, 0, 0, t26 + t354, -t243 + t352, 0, 0, 0, 0, 0, t357 + t366, -t223 * t275 - t322 - t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t353, -t349, -t343, t347, -t208, -t15 * t209 + t346, -t14 * t209 + t345, -t365, t369, t364, -t361, t356, -pkin(5) * t11 - t15 * t47 + t237 * t219 - t337 * t223 + t360, -pkin(5) * t10 + t15 * t252 + t337 * t219 + t237 * t223 - t359; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t252 * t47, t252 ^ 2 - t47 ^ 2, t10 + t328, -t11 - t327, t24, -g(1) * t88 + t12 * t252 + t219 * t239 + t223 * t257 + t350 * t6 + t3, g(1) * t89 + t12 * t47 + t5 * t350 + (-t257 - t4) * t219 + t239 * t223;];
tau_reg  = t1;

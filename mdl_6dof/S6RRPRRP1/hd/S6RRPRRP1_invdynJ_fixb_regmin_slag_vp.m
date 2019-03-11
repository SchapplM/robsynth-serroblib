% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:41:19
% EndTime: 2019-03-09 11:41:31
% DurationCPUTime: 5.17s
% Computational Cost: add. (9221->442), mult. (22004->556), div. (0->0), fcn. (16951->14), ass. (0->245)
t215 = sin(pkin(10));
t216 = cos(pkin(10));
t221 = sin(qJ(2));
t225 = cos(qJ(2));
t166 = -t215 * t221 + t216 * t225;
t156 = t166 * qJD(1);
t167 = t215 * t225 + t216 * t221;
t158 = t167 * qJD(1);
t220 = sin(qJ(4));
t224 = cos(qJ(4));
t110 = t224 * t156 - t158 * t220;
t223 = cos(qJ(5));
t288 = qJD(5) * t223;
t380 = -t110 * t223 + t288;
t212 = qJ(2) + pkin(10);
t206 = qJ(4) + t212;
t199 = sin(t206);
t222 = sin(qJ(1));
t226 = cos(qJ(1));
t253 = g(1) * t226 + g(2) * t222;
t379 = t253 * t199;
t219 = sin(qJ(5));
t191 = g(3) * t199;
t200 = cos(t206);
t276 = t253 * t200 + t191;
t157 = t167 * qJD(2);
t119 = -qJD(1) * t157 + t166 * qJDD(1);
t285 = qJD(1) * qJD(2);
t271 = t225 * t285;
t272 = t221 * t285;
t120 = t167 * qJDD(1) - t215 * t272 + t216 * t271;
t290 = qJD(4) * t224;
t291 = qJD(4) * t220;
t235 = t220 * t119 + t224 * t120 + t156 * t290 - t158 * t291;
t262 = -t224 * t119 + t220 * t120;
t247 = t156 * t220 + t224 * t158;
t351 = qJD(4) * t247;
t47 = t262 + t351;
t209 = t225 * pkin(2);
t203 = t209 + pkin(1);
t234 = pkin(2) * t272 - qJDD(1) * t203 + qJDD(3);
t90 = -t119 * pkin(3) + t234;
t15 = t47 * pkin(4) - pkin(9) * t235 + t90;
t289 = qJD(5) * t219;
t283 = qJD(2) + qJD(4);
t218 = -qJ(3) - pkin(7);
t187 = t218 * t225;
t176 = qJD(1) * t187;
t161 = t215 * t176;
t185 = t218 * t221;
t175 = qJD(1) * t185;
t320 = qJD(2) * pkin(2);
t165 = t175 + t320;
t117 = t216 * t165 + t161;
t339 = pkin(8) * t158;
t87 = qJD(2) * pkin(3) + t117 - t339;
t301 = t216 * t176;
t118 = t215 * t165 - t301;
t340 = pkin(8) * t156;
t91 = t118 + t340;
t51 = t220 * t87 + t224 * t91;
t42 = pkin(9) * t283 + t51;
t178 = -qJD(1) * t203 + qJD(3);
t126 = -pkin(3) * t156 + t178;
t56 = -pkin(4) * t110 - pkin(9) * t247 + t126;
t282 = qJDD(2) + qJDD(4);
t266 = qJD(2) * t218;
t152 = -qJD(3) * t221 + t225 * t266;
t116 = qJDD(2) * pkin(2) + t152 * qJD(1) + qJDD(1) * t185;
t151 = qJD(3) * t225 + t221 * t266;
t123 = t151 * qJD(1) - qJDD(1) * t187;
t73 = t216 * t116 - t123 * t215;
t49 = qJDD(2) * pkin(3) - pkin(8) * t120 + t73;
t74 = t215 * t116 + t216 * t123;
t54 = pkin(8) * t119 + t74;
t345 = -t224 * (qJD(4) * t87 + t54) - t220 * t49 + t91 * t291;
t9 = pkin(9) * t282 - t345;
t242 = t219 * t15 + t223 * t9 + t56 * t288 - t42 * t289;
t95 = t219 * t283 + t223 * t247;
t32 = qJD(5) * t95 + t219 * t235 - t223 * t282;
t256 = t223 * t283;
t93 = t219 * t247 - t256;
t3 = -qJ(6) * t32 - qJD(6) * t93 + t242;
t14 = t223 * t15;
t23 = t219 * t56 + t223 * t42;
t31 = -qJD(5) * t256 - t219 * t282 - t223 * t235 + t247 * t289;
t46 = qJDD(5) + t47;
t1 = pkin(5) * t46 + qJ(6) * t31 - t23 * qJD(5) - qJD(6) * t95 - t219 * t9 + t14;
t18 = -qJ(6) * t93 + t23;
t367 = qJD(5) - t110;
t349 = t18 * t367 + t1;
t378 = -t219 * t349 + t3 * t223 - t276;
t201 = pkin(2) * t216 + pkin(3);
t343 = pkin(2) * t215;
t260 = t201 * t224 - t220 * t343;
t124 = -t175 * t215 + t301;
t97 = t124 - t340;
t125 = t216 * t175 + t161;
t98 = t125 - t339;
t359 = -t260 * qJD(4) + t220 * t97 + t224 * t98;
t357 = t110 * t283;
t376 = t235 - t357;
t29 = t31 * t219;
t375 = t380 * t95 - t29;
t318 = t247 * t95;
t39 = t219 * t46;
t374 = t380 * t367 - t318 + t39;
t50 = -t220 * t91 + t224 * t87;
t41 = -pkin(4) * t283 - t50;
t373 = t41 * t110;
t371 = t247 * t110;
t341 = pkin(5) * t223;
t202 = pkin(4) + t341;
t217 = -qJ(6) - pkin(9);
t261 = -t199 * t217 + t200 * t202;
t369 = t261 + pkin(3) * cos(t212) + t209;
t366 = -t110 ^ 2 + t247 ^ 2;
t72 = pkin(4) * t247 - pkin(9) * t110;
t365 = -t126 * t110 + t276 + t345;
t208 = t223 * qJ(6);
t364 = -pkin(5) * t247 + t110 * t208;
t243 = -t223 * t31 - t95 * t289;
t314 = t223 * t93;
t316 = t219 * t95;
t250 = t314 + t316;
t328 = -t219 * t32 - t93 * t288;
t363 = t110 * t250 + t243 + t328;
t319 = t247 * t93;
t294 = t220 * t201 + t224 * t343;
t312 = t294 * qJD(4) - t220 * t98 + t224 * t97;
t358 = t367 * t247;
t355 = g(1) * t222 - g(2) * t226;
t356 = t355 * t199;
t287 = t247 * qJD(2);
t354 = t287 - t262;
t330 = t221 * pkin(2);
t131 = pkin(3) * t158 + qJD(1) * t330;
t61 = t131 + t72;
t353 = t219 * t61 + t223 * t359;
t296 = t223 * t226;
t300 = t219 * t222;
t145 = t200 * t300 + t296;
t298 = t222 * t223;
t299 = t219 * t226;
t147 = -t200 * t299 + t298;
t352 = -g(1) * t147 + g(2) * t145;
t333 = g(3) * t200;
t350 = t333 - t379;
t22 = -t219 * t42 + t223 * t56;
t37 = t41 * t289;
t348 = -t22 * t247 + t223 * t379 + t37;
t264 = t220 * t54 - t224 * t49 + t91 * t290 + t87 * t291;
t10 = -pkin(4) * t282 + t264;
t332 = g(3) * t219;
t347 = t10 * t219 + t200 * t332 + t23 * t247 + t41 * t288;
t346 = -t126 * t247 - t264 - t350;
t344 = t95 ^ 2;
t342 = pkin(5) * t219;
t331 = g(3) * t225;
t17 = -qJ(6) * t95 + t22;
t11 = pkin(5) * t367 + t17;
t329 = -t17 + t11;
t326 = t219 * t72 + t223 * t50;
t127 = t216 * t185 + t187 * t215;
t102 = -pkin(8) * t167 + t127;
t128 = t215 * t185 - t216 * t187;
t103 = pkin(8) * t166 + t128;
t68 = t102 * t220 + t103 * t224;
t65 = t223 * t68;
t122 = t166 * t220 + t167 * t224;
t134 = -pkin(3) * t166 - t203;
t246 = t224 * t166 - t167 * t220;
t66 = -pkin(4) * t246 - pkin(9) * t122 + t134;
t324 = t219 * t66 + t65;
t207 = t223 * qJD(6);
t150 = pkin(9) + t294;
t295 = -qJ(6) - t150;
t259 = qJD(5) * t295;
t306 = t110 * t219;
t323 = qJ(6) * t306 + t219 * t259 + t207 - t353;
t58 = t223 * t61;
t322 = t223 * t259 - t58 + (-qJD(6) + t359) * t219 + t364;
t317 = t11 * t223;
t40 = t223 * t46;
t160 = t166 * qJD(2);
t75 = qJD(4) * t246 - t157 * t220 + t160 * t224;
t315 = t223 * t75;
t313 = t223 * t95;
t265 = qJD(5) * t217;
t311 = t207 - t326 + (qJ(6) * t110 + t265) * t219;
t70 = t223 * t72;
t310 = t223 * t265 - t70 + (-qJD(6) + t50) * t219 + t364;
t308 = t367 * t219;
t304 = t122 * t219;
t101 = t216 * t151 + t215 * t152;
t213 = t221 ^ 2;
t292 = -t225 ^ 2 + t213;
t284 = t225 * qJDD(1);
t248 = t102 * t224 - t103 * t220;
t100 = -t151 * t215 + t216 * t152;
t78 = -pkin(8) * t160 + t100;
t79 = -pkin(8) * t157 + t101;
t26 = qJD(4) * t248 + t220 * t78 + t224 * t79;
t205 = t221 * t320;
t132 = pkin(3) * t157 + t205;
t76 = qJD(4) * t122 + t224 * t157 + t160 * t220;
t35 = pkin(4) * t76 - pkin(9) * t75 + t132;
t280 = t219 * t35 + t223 * t26 + t66 * t288;
t279 = qJD(5) * pkin(9) * t367;
t274 = t122 * t288;
t273 = -t10 - t333;
t270 = pkin(8) - t218 + t342;
t268 = -qJD(5) * t56 - t9;
t258 = t367 * t223;
t149 = -pkin(4) - t260;
t255 = -t42 * t288 + t14;
t254 = -pkin(9) * t46 - t373;
t251 = -t150 * t46 - t373;
t249 = -qJ(6) * t75 - qJD(6) * t122;
t245 = t199 * t202 + t200 * t217;
t244 = t40 + (-t289 + t306) * t367;
t240 = -0.2e1 * pkin(1) * t285 - pkin(7) * qJDD(2);
t239 = t219 * t75 + t274;
t238 = -pkin(1) - t369;
t227 = qJD(2) ^ 2;
t232 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t227 + t355;
t228 = qJD(1) ^ 2;
t231 = pkin(1) * t228 - pkin(7) * qJDD(1) + t253;
t6 = t32 * pkin(5) + qJDD(6) + t10;
t27 = qJD(4) * t68 + t220 * t79 - t224 * t78;
t186 = pkin(9) * t223 + t208;
t184 = t217 * t219;
t148 = t200 * t296 + t300;
t146 = -t200 * t298 + t299;
t130 = t150 * t223 + t208;
t129 = t295 * t219;
t92 = t93 ^ 2;
t64 = t223 * t66;
t36 = t93 * pkin(5) + qJD(6) + t41;
t34 = t223 * t35;
t24 = -qJ(6) * t304 + t324;
t20 = -pkin(5) * t246 - t122 * t208 - t219 * t68 + t64;
t5 = -qJ(6) * t274 + (-qJD(5) * t68 + t249) * t219 + t280;
t4 = pkin(5) * t76 - t219 * t26 + t34 + t249 * t223 + (-t65 + (qJ(6) * t122 - t66) * t219) * qJD(5);
t2 = [qJDD(1), t355, t253, qJDD(1) * t213 + 0.2e1 * t221 * t271, 0.2e1 * t221 * t284 - 0.2e1 * t292 * t285, qJDD(2) * t221 + t225 * t227, qJDD(2) * t225 - t221 * t227, 0, t221 * t240 + t225 * t232, -t221 * t232 + t225 * t240, -t100 * t158 + t101 * t156 - t117 * t160 - t118 * t157 + t119 * t128 - t120 * t127 + t166 * t74 - t167 * t73 - t253, t74 * t128 + t118 * t101 + t73 * t127 + t117 * t100 - t234 * t203 + t178 * t205 - g(1) * (-t203 * t222 - t218 * t226) - g(2) * (t203 * t226 - t218 * t222) t122 * t235 + t247 * t75, t110 * t75 - t122 * t47 + t235 * t246 - t247 * t76, t122 * t282 + t283 * t75, t246 * t282 - t283 * t76, 0, -t110 * t132 + t126 * t76 + t134 * t47 + t200 * t355 - t246 * t90 + t248 * t282 - t27 * t283, t90 * t122 + t126 * t75 + t132 * t247 + t134 * t235 - t26 * t283 - t282 * t68 - t356, t122 * t243 + t75 * t313, -t250 * t75 + (t29 - t223 * t32 + (t219 * t93 - t313) * qJD(5)) * t122, t122 * t40 + t246 * t31 + t76 * t95 + (-t122 * t289 + t315) * t367, -t239 * t367 + t246 * t32 - t46 * t304 - t76 * t93, -t246 * t46 + t367 * t76 (-t288 * t68 + t34) * t367 + t64 * t46 - t255 * t246 + t22 * t76 + t27 * t93 - t248 * t32 + t41 * t274 - g(1) * t146 - g(2) * t148 + ((-qJD(5) * t66 - t26) * t367 - t68 * t46 - t268 * t246 + t10 * t122 + t41 * t75) * t219 -(-t289 * t68 + t280) * t367 - t324 * t46 + t242 * t246 - t23 * t76 + t27 * t95 + t248 * t31 + t41 * t315 - g(1) * t145 - g(2) * t147 + (t10 * t223 - t37) * t122, t20 * t31 - t24 * t32 - t4 * t95 - t5 * t93 + (-t18 * t219 - t317) * t75 + t356 + (-t1 * t223 - t219 * t3 + (t11 * t219 - t18 * t223) * qJD(5)) * t122, t3 * t24 + t18 * t5 + t1 * t20 + t11 * t4 + t6 * (pkin(5) * t304 - t248) + t36 * (pkin(5) * t239 + t27) + (-g(1) * t270 + g(2) * t238) * t226 + (-g(1) * t238 - g(2) * t270) * t222; 0, 0, 0, -t221 * t228 * t225, t292 * t228, t221 * qJDD(1), t284, qJDD(2), t221 * t231 - t331, g(3) * t221 + t225 * t231 (t118 + t124) * t158 + (t117 - t125) * t156 + (t119 * t215 - t120 * t216) * pkin(2), -t117 * t124 - t118 * t125 + (-t331 + t215 * t74 + t216 * t73 + (-qJD(1) * t178 + t253) * t221) * pkin(2), -t371, t366, t376, t354, t282, t110 * t131 + t260 * t282 - t312 * t283 + t346, -t131 * t247 - t282 * t294 + t359 * t283 + t365, t375, t363, t374, t244 + t319, -t358, t149 * t32 + t312 * t93 + t273 * t223 + t251 * t219 + (-t150 * t288 + t359 * t219 - t58) * t367 + t348, -t149 * t31 + t312 * t95 + t251 * t223 - t219 * t379 + (t150 * t289 + t353) * t367 + t347, -t11 * t258 + t129 * t31 - t130 * t32 - t322 * t95 - t323 * t93 + t378, t3 * t130 + t1 * t129 + t6 * (t149 - t341) - g(3) * t369 + (pkin(5) * t308 + t312) * t36 + t323 * t18 + t322 * t11 + t253 * (pkin(3) * sin(t212) + t330 + t245); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156 ^ 2 - t158 ^ 2, t117 * t158 - t118 * t156 + t234 - t355, 0, 0, 0, 0, 0, t262 + t287 + 0.2e1 * t351, t235 + t357, 0, 0, 0, 0, 0, t244 - t319, -t258 * t367 - t318 - t39 (t314 - t316) * t110 - t243 + t328, -t247 * t36 + t349 * t223 + (-t11 * t367 + t3) * t219 - t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t371, t366, t376, t354, t282, t283 * t51 + t346, t283 * t50 + t365, t375, t363, t374, -t308 * t367 + t319 + t40, -t358, -pkin(4) * t32 - t70 * t367 - t51 * t93 + (t367 * t50 + t254) * t219 + (t273 - t279) * t223 + t348, pkin(4) * t31 + t326 * t367 - t51 * t95 + t254 * t223 + (-t379 + t279) * t219 + t347, t184 * t31 - t186 * t32 - t310 * t95 - t311 * t93 - t367 * t317 + t378, t3 * t186 + t1 * t184 - t6 * t202 - g(3) * t261 + (t342 * t367 - t51) * t36 + t311 * t18 + t310 * t11 + t253 * t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * t93, -t92 + t344, t367 * t93 - t31, t367 * t95 - t32, t46, t367 * t23 - t41 * t95 + (t268 + t191) * t219 + t255 + t352, g(1) * t148 - g(2) * t146 + t223 * t191 + t22 * t367 + t41 * t93 - t242, pkin(5) * t31 - t329 * t93, t329 * t18 + (t199 * t332 - t36 * t95 + t1 + t352) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92 - t344, t11 * t95 + t18 * t93 + t350 + t6;];
tau_reg  = t2;

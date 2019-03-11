% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:19
% EndTime: 2019-03-09 06:05:32
% DurationCPUTime: 5.14s
% Computational Cost: add. (6926->507), mult. (14194->645), div. (0->0), fcn. (9583->14), ass. (0->250)
t199 = sin(qJ(4));
t203 = cos(qJ(3));
t300 = qJD(1) * t203;
t276 = t199 * t300;
t347 = pkin(9) + pkin(8);
t277 = qJD(4) * t347;
t196 = sin(pkin(10));
t178 = pkin(1) * t196 + pkin(7);
t156 = t178 * qJD(1);
t200 = sin(qJ(3));
t109 = qJD(2) * t203 - t200 * t156;
t252 = pkin(3) * t200 - pkin(8) * t203;
t141 = t252 * qJD(1);
t202 = cos(qJ(4));
t322 = t202 * t109 + t199 * t141;
t371 = pkin(9) * t276 - t199 * t277 - t322;
t124 = t202 * t141;
t309 = t202 * t203;
t245 = pkin(4) * t200 - pkin(9) * t309;
t370 = t245 * qJD(1) - t199 * t109 + t202 * t277 + t124;
t188 = t203 * qJDD(1);
t288 = qJD(1) * qJD(3);
t132 = t200 * t288 + qJDD(4) - t188;
t126 = qJDD(5) + t132;
t195 = qJ(4) + qJ(5);
t190 = sin(t195);
t192 = qJ(1) + pkin(10);
t186 = sin(t192);
t187 = cos(t192);
t251 = g(1) * t187 + g(2) * t186;
t237 = t251 * t200;
t339 = g(3) * t203;
t159 = t347 * t199;
t160 = t347 * t202;
t198 = sin(qJ(5));
t346 = cos(qJ(5));
t89 = -t198 * t159 + t346 * t160;
t369 = t89 * t126 + (t237 - t339) * t190;
t295 = qJD(3) * t202;
t301 = qJD(1) * t200;
t136 = -t199 * t301 + t295;
t297 = qJD(3) * t199;
t137 = t202 * t301 + t297;
t235 = t198 * t136 + t346 * t137;
t99 = -qJD(3) * pkin(3) - t109;
t73 = -pkin(4) * t136 + t99;
t79 = -t346 * t136 + t137 * t198;
t25 = pkin(5) * t79 - qJ(6) * t235 + t73;
t368 = t25 * t79;
t367 = t73 * t79;
t338 = t235 * t79;
t314 = t198 * t199;
t233 = t346 * t202 - t314;
t353 = qJD(4) + qJD(5);
t268 = t346 * qJD(5);
t355 = t346 * qJD(4) + t268;
t324 = t202 * t355 - t233 * t300 - t314 * t353;
t139 = t198 * t202 + t346 * t199;
t85 = t353 * t139;
t323 = -t139 * t300 + t85;
t151 = t178 * qJDD(1);
t366 = qJD(2) * qJD(3) + t151;
t260 = qJD(4) + t300;
t283 = t200 * qJDD(1);
t221 = t260 * qJD(3) + t283;
t287 = qJD(1) * qJD(4);
t266 = t200 * t287;
t246 = -qJDD(3) + t266;
t230 = t246 * t199;
t211 = t221 * t202 - t230;
t244 = qJD(3) * qJD(4) + t283;
t267 = t203 * t288;
t219 = t244 + t267;
t259 = t219 * t199 + t202 * t266;
t228 = t202 * qJDD(3) - t259;
t289 = qJD(5) * t198;
t21 = -t136 * t268 + t137 * t289 - t198 * t228 - t346 * t211;
t296 = qJD(3) * t200;
t264 = t203 * t21 + t235 * t296;
t112 = t233 * t200;
t176 = -qJD(4) + t300;
t166 = -qJD(5) + t176;
t294 = qJD(3) * t203;
t258 = t346 * t294;
t274 = t199 * t294;
t47 = t198 * t274 + t85 * t200 - t202 * t258;
t335 = t112 * t126 + t47 * t166;
t364 = t264 - t335;
t290 = qJD(4) * t202;
t363 = t200 * t290 + t274;
t362 = t109 * qJD(3);
t348 = t235 ^ 2;
t361 = -t79 ^ 2 + t348;
t13 = -t166 * t79 - t21;
t40 = pkin(5) * t235 + qJ(6) * t79;
t234 = -t346 * t159 - t198 * t160;
t359 = t234 * qJD(5) - t370 * t198 + t371 * t346;
t358 = t89 * qJD(5) + t371 * t198 + t370 * t346;
t110 = t200 * qJD(2) + t203 * t156;
t100 = qJD(3) * pkin(8) + t110;
t197 = cos(pkin(10));
t180 = -pkin(1) * t197 - pkin(2);
t127 = -pkin(3) * t203 - pkin(8) * t200 + t180;
t101 = t127 * qJD(1);
t55 = t202 * t100 + t199 * t101;
t46 = pkin(9) * t136 + t55;
t22 = t235 * qJD(5) + t198 * t211 - t346 * t228;
t263 = -t203 * t22 + t79 * t296;
t292 = qJD(4) * t199;
t257 = -t110 + (-t276 + t292) * pkin(4);
t118 = t126 * qJ(6);
t148 = t166 * qJD(6);
t357 = t118 - t148;
t111 = t139 * t200;
t291 = qJD(4) * t200;
t273 = t199 * t291;
t313 = t199 * t200;
t48 = t199 * t258 - t198 * t273 - t289 * t313 + (t198 * t294 + t200 * t355) * t202;
t334 = -t111 * t126 + t48 * t166;
t356 = t263 + t334;
t312 = t199 * t203;
t104 = t186 * t312 + t187 * t202;
t106 = t186 * t202 - t187 * t312;
t354 = -g(1) * t106 + g(2) * t104;
t119 = t126 * pkin(5);
t352 = t119 - qJDD(6);
t61 = qJDD(3) * pkin(8) + t200 * qJDD(2) + t203 * t151 + t362;
t144 = t252 * qJD(3);
t74 = qJD(1) * t144 + t127 * qJDD(1);
t282 = t101 * t290 + t199 * t74 + t202 * t61;
t231 = -t100 * t292 + t282;
t11 = t228 * pkin(9) + t231;
t54 = -t100 * t199 + t202 * t101;
t45 = -pkin(9) * t137 + t54;
t38 = -pkin(4) * t176 + t45;
t284 = t199 * qJDD(3);
t69 = t202 * t74;
t8 = -t199 * t61 + t69 - (t284 + (t267 + t283) * t202) * pkin(9) + t132 * pkin(4) - t46 * qJD(4);
t270 = t198 * t11 + t46 * t268 + t38 * t289 - t346 * t8;
t318 = t190 * t200;
t191 = cos(t195);
t317 = t190 * t203;
t94 = t186 * t317 + t187 * t191;
t96 = -t186 * t191 + t187 * t317;
t224 = g(1) * t96 + g(2) * t94 + g(3) * t318 - t270;
t212 = t235 * t25 - t224 - t352;
t351 = -t73 * t235 + t224;
t350 = -t166 * t235 - t22;
t114 = t202 * t127;
t311 = t200 * t202;
t319 = t178 * t199;
t63 = -pkin(9) * t311 + t114 + (-pkin(4) - t319) * t203;
t140 = t178 * t309;
t303 = t199 * t127 + t140;
t72 = -pkin(9) * t313 + t303;
t238 = t198 * t63 + t346 * t72;
t304 = t202 * t144 + t296 * t319;
t32 = t245 * qJD(3) + (-t140 + (pkin(9) * t200 - t127) * t199) * qJD(4) + t304;
t305 = t127 * t290 + t199 * t144;
t35 = (-t200 * t295 - t203 * t292) * t178 - t363 * pkin(9) + t305;
t349 = -t238 * qJD(5) - t198 * t35 + t346 * t32;
t340 = g(3) * t200;
t337 = -t112 * t22 + t47 * t79;
t336 = t323 * pkin(5) - t324 * qJ(6) - qJD(6) * t139 + t257;
t333 = -qJ(6) * t301 + t359;
t332 = pkin(5) * t301 + t358;
t280 = t346 * t46;
t17 = t198 * t38 + t280;
t330 = t166 * t17;
t329 = t198 * t46;
t19 = t346 * t45 - t329;
t325 = pkin(4) * t268 + qJD(6) - t19;
t321 = qJDD(3) * pkin(3);
t320 = t137 * t176;
t316 = t191 * t200;
t315 = t191 * t203;
t310 = t200 * t347;
t308 = t203 * t176;
t16 = t346 * t38 - t329;
t307 = qJD(6) - t16;
t306 = qJDD(2) - g(3);
t177 = pkin(4) * t313;
t117 = t200 * t178 + t177;
t193 = t200 ^ 2;
t302 = -t203 ^ 2 + t193;
t157 = qJD(1) * t180;
t298 = qJD(3) * t136;
t293 = qJD(4) * t136;
t279 = t156 * t294 + t366 * t200;
t86 = t363 * pkin(4) + t178 * t294;
t185 = pkin(4) * t202 + pkin(3);
t278 = pkin(4) * t199 + pkin(7);
t275 = t176 * t297;
t272 = t176 * t291;
t265 = t346 * t11 + t198 * t8 + t38 * t268 - t46 * t289;
t262 = -qJD(4) * t101 - t61;
t261 = t176 * t178 + t100;
t18 = t198 * t45 + t280;
t256 = pkin(4) * t289 - t18;
t255 = -pkin(5) * t318 + qJ(6) * t316;
t254 = -g(1) * t94 + g(2) * t96;
t95 = t186 * t315 - t187 * t190;
t97 = t186 * t190 + t187 * t315;
t253 = g(1) * t95 - g(2) * t97;
t250 = g(1) * t186 - g(2) * t187;
t201 = sin(qJ(1));
t204 = cos(qJ(1));
t249 = g(1) * t201 - g(2) * t204;
t248 = -t111 * t21 + t235 * t48;
t247 = -t176 + t260;
t243 = t185 * t203 + pkin(2) + t310;
t242 = pkin(5) * t191 + qJ(6) * t190 + t185;
t240 = -t198 * t72 + t346 * t63;
t236 = t198 * t32 + t63 * t268 - t72 * t289 + t346 * t35;
t232 = -t199 * t132 + t176 * t290;
t229 = t246 * t203;
t227 = -qJD(1) * t157 + t251;
t226 = -g(3) * t315 + t234 * t126 + t251 * t316;
t225 = -pkin(8) * t132 - t176 * t99;
t62 = -t203 * qJDD(2) + t279 - t321;
t223 = g(1) * t97 + g(2) * t95 + g(3) * t316 - t265;
t222 = 0.2e1 * t157 * qJD(3) - t178 * qJDD(3);
t220 = -qJD(4) * pkin(8) * t176 - t321 + t339 + t62;
t218 = -g(1) * (-t96 * pkin(5) + t97 * qJ(6)) - g(2) * (-t94 * pkin(5) + t95 * qJ(6));
t206 = qJD(3) ^ 2;
t214 = -0.2e1 * qJDD(1) * t180 - t178 * t206 + t250;
t213 = -t16 * t166 + t223;
t34 = -t228 * pkin(4) + t62;
t207 = qJD(1) ^ 2;
t184 = -t346 * pkin(4) - pkin(5);
t179 = pkin(4) * t198 + qJ(6);
t147 = qJDD(3) * t203 - t200 * t206;
t146 = qJDD(3) * t200 + t203 * t206;
t115 = t137 * t296;
t107 = t186 * t199 + t187 * t309;
t105 = -t186 * t309 + t187 * t199;
t75 = -pkin(5) * t233 - qJ(6) * t139 - t185;
t51 = pkin(5) * t111 - qJ(6) * t112 + t117;
t31 = pkin(4) * t137 + t40;
t28 = t203 * pkin(5) - t240;
t27 = -qJ(6) * t203 + t238;
t15 = -t166 * qJ(6) + t17;
t14 = t166 * pkin(5) + t307;
t12 = pkin(5) * t48 + qJ(6) * t47 - qJD(6) * t112 + t86;
t5 = -pkin(5) * t296 - t349;
t4 = qJ(6) * t296 - qJD(6) * t203 + t236;
t3 = t22 * pkin(5) + t21 * qJ(6) - qJD(6) * t235 + t34;
t2 = t270 - t352;
t1 = t265 + t357;
t6 = [qJDD(1), t249, g(1) * t204 + g(2) * t201 (t249 + (t196 ^ 2 + t197 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t193 + 0.2e1 * t200 * t267, 0.2e1 * t200 * t188 - 0.2e1 * t302 * t288, t146, t147, 0, t200 * t222 + t203 * t214, -t200 * t214 + t203 * t222, -t137 * t273 + (t137 * t294 - t200 * t230 + t219 * t311) * t202 (t136 * t202 - t137 * t199) * t294 + ((-t137 * qJD(4) + t228) * t202 + (-t293 - t211) * t199) * t200, t115 + (t229 + t272) * t199 + ((t132 - t188) * t200 + (-t176 - t260) * t294) * t202 (-t228 + t275) * t203 + (t232 + t298) * t200, -t132 * t203 - t176 * t296 -(-t127 * t292 + t304) * t176 + t114 * t132 - g(1) * t105 - g(2) * t107 + (-t178 * t298 - t69 + t261 * t290 + (qJD(3) * t99 - t132 * t178 - t262) * t199) * t203 + (t54 * qJD(3) - t178 * t228 + t62 * t199 + t290 * t99) * t200, t305 * t176 - t303 * t132 - g(1) * t104 - g(2) * t106 + (-t261 * t292 + (t137 * t178 + t202 * t99) * qJD(3) + t282) * t203 + (-t99 * t292 + t62 * t202 + (t284 + (qJDD(1) * t202 - t199 * t287) * t200) * t178 + (t178 * t202 * t247 - t55) * qJD(3)) * t200, -t112 * t21 - t235 * t47, -t248 + t337, t264 + t335, t334 - t263, -t126 * t203 - t166 * t296, t34 * t111 + t117 * t22 + t240 * t126 + t16 * t296 - t166 * t349 + t270 * t203 + t73 * t48 + t86 * t79 + t253, t34 * t112 - t117 * t21 - t126 * t238 + t166 * t236 - t17 * t296 + t203 * t265 + t235 * t86 - t73 * t47 + t254, t111 * t3 + t12 * t79 - t126 * t28 - t14 * t296 + t166 * t5 + t2 * t203 + t22 * t51 + t25 * t48 + t253, -t1 * t111 + t2 * t112 - t14 * t47 - t15 * t48 + t200 * t250 - t28 * t21 - t27 * t22 + t235 * t5 - t4 * t79, -t1 * t203 - t112 * t3 - t12 * t235 + t126 * t27 + t15 * t296 - t166 * t4 + t21 * t51 + t25 * t47 - t254, t1 * t27 + t15 * t4 + t3 * t51 + t25 * t12 + t2 * t28 + t14 * t5 - g(1) * (-t201 * pkin(1) - t95 * pkin(5) - t94 * qJ(6)) - g(2) * (t204 * pkin(1) + t97 * pkin(5) + t96 * qJ(6)) + (-g(1) * t278 - g(2) * t243) * t187 + (g(1) * t243 - g(2) * t278) * t186; 0, 0, 0, t306, 0, 0, 0, 0, 0, t147, -t146, 0, 0, 0, 0, 0 (t228 + t275) * t203 + (t232 - t298) * t200, t115 + (t229 - t272) * t199 + ((-t132 - t188) * t200 - t247 * t294) * t202, 0, 0, 0, 0, 0, t356, t364, t356, t248 + t337, -t364, t1 * t112 + t111 * t2 + t14 * t48 - t15 * t47 - t203 * t3 + t25 * t296 - g(3); 0, 0, 0, 0, -t203 * t207 * t200, t302 * t207, t283, t188, qJDD(3), t110 * qJD(3) + t227 * t200 + t306 * t203 - t279, t362 + (qJD(3) * t156 - t306) * t200 + (t227 - t366) * t203, -t246 * t199 ^ 2 + (t199 * t221 - t320) * t202 (-t259 + t320) * t199 + (t293 + 0.2e1 * t284 + t244 * t202 + (-t273 + (-t136 + t295) * t203) * qJD(1)) * t202 (-t137 * t200 + t202 * t308) * qJD(1) - t232, t176 * t292 + t202 * t132 + (-t136 * t200 - t199 * t308) * qJD(1), t176 * t301, -pkin(3) * t259 + t124 * t176 - t54 * t301 + t110 * t136 + (-t109 * t176 + t225) * t199 + (t237 - t220) * t202, -t322 * t176 + t55 * t301 - t110 * t137 + (-pkin(3) * t219 + t225) * t202 + ((pkin(3) * t287 - t251) * t200 + t220) * t199, -t21 * t139 + t235 * t324, -t139 * t22 - t21 * t233 - t235 * t323 - t324 * t79, t139 * t126 - t166 * t324 - t235 * t301, t126 * t233 + t166 * t323 + t301 * t79, t166 * t301, -t16 * t301 + t166 * t358 - t185 * t22 - t233 * t34 + t257 * t79 + t323 * t73 + t226, t34 * t139 + t359 * t166 + t17 * t301 + t185 * t21 + t257 * t235 + t324 * t73 - t369, t14 * t301 + t166 * t332 + t75 * t22 - t233 * t3 + t25 * t323 + t336 * t79 + t226, t1 * t233 + t2 * t139 + t14 * t324 - t15 * t323 - t203 * t251 + t21 * t234 - t89 * t22 + t235 * t332 - t333 * t79 - t340, -t3 * t139 - t15 * t301 - t333 * t166 + t75 * t21 - t336 * t235 - t324 * t25 + t369, -g(3) * t310 + t1 * t89 + t14 * t332 + t15 * t333 - t2 * t234 - t242 * t339 + t25 * t336 + t3 * t75 + t251 * (t242 * t200 - t203 * t347); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137 * t136, -t136 ^ 2 + t137 ^ 2, t136 * t176 + t211, t228 - t320, t132, -t100 * t290 - t99 * t137 - t55 * t176 + t69 + (t262 + t340) * t199 + t354, g(1) * t107 - g(2) * t105 + g(3) * t311 - t136 * t99 - t176 * t54 - t231, t338, t361, t13, t350, t126, -t18 * t166 + (t126 * t346 - t137 * t79 + t166 * t289) * pkin(4) + t351, -t19 * t166 + t367 + (-t198 * t126 - t137 * t235 + t166 * t268) * pkin(4) + t223, -t184 * t126 + t166 * t256 - t31 * t79 - t212, -t179 * t22 - t184 * t21 + (t15 + t256) * t235 + (t14 - t325) * t79, t179 * t126 - t166 * t325 + t235 * t31 - t223 + t357 - t368, t1 * t179 + t2 * t184 - t25 * t31 - t14 * t18 - g(3) * (-t177 + t255) + t325 * t15 + (t14 * t289 + t354) * pkin(4) + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338, t361, t13, t350, t126, -t330 + t351, t213 + t367, -t40 * t79 + t119 - t212 - t330, pkin(5) * t21 - t22 * qJ(6) + (t15 - t17) * t235 + (t14 - t307) * t79, t235 * t40 + 0.2e1 * t118 - 0.2e1 * t148 - t213 - t368, -t2 * pkin(5) - g(3) * t255 + t1 * qJ(6) - t14 * t17 + t15 * t307 - t25 * t40 + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 + t338, t13, -t166 ^ 2 - t348, t15 * t166 + t212;];
tau_reg  = t6;

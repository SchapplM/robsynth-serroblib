% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:10
% EndTime: 2019-03-09 03:03:17
% DurationCPUTime: 4.21s
% Computational Cost: add. (7324->524), mult. (15875->621), div. (0->0), fcn. (11130->14), ass. (0->280)
t357 = 2 * qJD(3);
t192 = sin(pkin(9));
t168 = pkin(1) * t192 + pkin(7);
t154 = t168 * qJDD(1);
t250 = -(qJD(2) * qJD(3)) - t154;
t191 = sin(pkin(10));
t197 = sin(qJ(3));
t200 = cos(qJ(3));
t316 = cos(pkin(10));
t146 = t191 * t200 + t197 * t316;
t136 = t146 * qJD(3);
t137 = t146 * qJD(1);
t258 = t316 * t200;
t224 = -t191 * t197 + t258;
t163 = qJD(1) * t258;
t280 = qJD(1) * qJD(3);
t263 = t197 * t280;
t214 = qJDD(1) * t146 - t191 * t263;
t351 = qJD(3) * t163 + t214;
t356 = t137 * t136 - t224 * t351;
t196 = sin(qJ(5));
t199 = cos(qJ(5));
t116 = qJD(3) * t196 + t137 * t199;
t277 = qJD(3) * qJD(5);
t282 = qJD(5) * t196;
t223 = t196 * qJDD(3) - t137 * t282 + (t351 + t277) * t199;
t326 = qJ(6) * t223;
t276 = t197 * qJDD(1);
t238 = -qJDD(1) * t258 + t191 * t276;
t102 = qJD(1) * t136 + t238;
t101 = qJDD(5) + t102;
t341 = pkin(5) * t101;
t182 = t200 * qJDD(2);
t279 = qJD(1) * qJD(4);
t156 = t168 * qJD(1);
t251 = qJ(4) * qJD(1) + t156;
t353 = t200 * t251;
t66 = qJDD(3) * pkin(3) + t182 - qJD(3) * t353 + (-qJ(4) * qJDD(1) + t250 - t279) * t197;
t275 = t200 * qJDD(1);
t267 = -t197 * qJDD(2) + t200 * t250;
t283 = qJD(3) * t197;
t95 = -t156 * t283 - t267;
t70 = t200 * t279 + (-t263 + t275) * qJ(4) + t95;
t34 = t191 * t66 + t316 * t70;
t32 = qJDD(3) * pkin(8) + t34;
t184 = t200 * qJD(2);
t119 = -t197 * t251 + t184;
t111 = qJD(3) * pkin(3) + t119;
t284 = qJD(2) * t197;
t120 = t284 + t353;
t259 = t316 * t120;
t69 = t191 * t111 + t259;
t63 = qJD(3) * pkin(8) + t69;
t193 = cos(pkin(9));
t170 = -t193 * pkin(1) - pkin(2);
t185 = t200 * pkin(3);
t352 = t170 - t185;
t133 = qJD(1) * t352 + qJD(4);
t285 = qJD(1) * t197;
t134 = t191 * t285 - t163;
t79 = pkin(4) * t134 - pkin(8) * t137 + t133;
t37 = t196 * t79 + t199 * t63;
t315 = pkin(1) * qJDD(1);
t269 = t193 * t315;
t274 = pkin(3) * t263 + qJDD(4);
t46 = -qJDD(1) * pkin(2) - pkin(3) * t275 + t102 * pkin(4) - pkin(8) * t351 - t269 + t274;
t4 = -qJD(5) * t37 - t196 * t32 + t199 * t46;
t1 = -qJD(6) * t116 - t326 + t341 + t4;
t131 = qJD(5) + t134;
t114 = -t199 * qJD(3) + t137 * t196;
t22 = -qJ(6) * t114 + t37;
t324 = t131 * t22;
t355 = t1 + t324;
t253 = t131 * t196;
t354 = t116 * t253;
t187 = qJ(3) + pkin(10);
t175 = sin(t187);
t188 = qJ(1) + pkin(9);
t176 = sin(t188);
t178 = cos(t188);
t244 = g(1) * t178 + g(2) * t176;
t227 = t244 * t175;
t339 = g(1) * t176;
t261 = g(2) * t178 - t339;
t350 = qJD(5) * t116;
t291 = t178 * t199;
t177 = cos(t187);
t293 = t177 * t196;
t124 = t176 * t293 + t291;
t292 = t178 * t196;
t294 = t176 * t199;
t126 = -t177 * t292 + t294;
t334 = g(3) * t175;
t349 = -g(1) * t126 + g(2) * t124 + t196 * t334;
t139 = t224 * qJD(3);
t300 = t139 * t199;
t225 = t146 * t282 - t300;
t94 = t199 * t101;
t348 = -t131 * t225 + t146 * t94;
t216 = -t244 * t177 - t334;
t347 = t116 ^ 2;
t346 = t137 ^ 2;
t198 = sin(qJ(1));
t345 = pkin(1) * t198;
t344 = pkin(3) * t191;
t343 = pkin(3) * t197;
t342 = pkin(4) * t177;
t337 = g(1) * t198;
t333 = g(3) * t177;
t332 = g(3) * t200;
t331 = t114 * pkin(5);
t330 = t199 * pkin(5);
t36 = -t196 * t63 + t199 * t79;
t21 = -qJ(6) * t116 + t36;
t20 = pkin(5) * t131 + t21;
t329 = -t21 + t20;
t298 = t146 * t199;
t206 = -t199 * qJDD(3) + t196 * t351;
t58 = t206 + t350;
t328 = -t114 * t300 - t58 * t298;
t33 = -t191 * t70 + t316 * t66;
t108 = t191 * t120;
t72 = t119 * t316 - t108;
t89 = pkin(3) * t285 + pkin(4) * t137 + pkin(8) * t134;
t42 = t196 * t89 + t199 * t72;
t327 = t116 * t136 - t223 * t224;
t290 = qJ(4) + t168;
t142 = t290 * t200;
t256 = t290 * t197;
t99 = t142 * t316 - t191 * t256;
t92 = t199 * t99;
t93 = -pkin(4) * t224 - pkin(8) * t146 + t352;
t52 = t196 * t93 + t92;
t325 = qJ(6) * t58;
t323 = t131 * t36;
t322 = t131 * t37;
t321 = t196 * t223;
t320 = -t146 * t102 - t139 * t134;
t281 = qJD(5) * t199;
t319 = -t114 * t281 - t196 * t58;
t167 = pkin(8) + t344;
t289 = qJ(6) + t167;
t254 = qJD(5) * t289;
t304 = t134 * t196;
t318 = -qJ(6) * t304 + qJD(6) * t199 - t196 * t254 - t42;
t41 = -t196 * t72 + t199 * t89;
t317 = -pkin(5) * t137 - qJD(6) * t196 - t41 + (-qJ(6) * t134 - t254) * t199;
t314 = t101 * t196;
t313 = t114 * t131;
t312 = t114 * t134;
t311 = t114 * t137;
t310 = t114 * t196;
t309 = t116 * t114;
t308 = t116 * t131;
t307 = t116 * t137;
t306 = t116 * t199;
t305 = t131 * t137;
t303 = t137 * t134;
t301 = t139 * t196;
t299 = t146 * t196;
t297 = t156 * t197;
t296 = t156 * t200;
t295 = t175 * t178;
t288 = (g(1) * t291 + g(2) * t294) * t175;
t174 = t185 + pkin(2);
t201 = cos(qJ(1));
t186 = t201 * pkin(1);
t287 = t178 * t174 + t186;
t189 = t197 ^ 2;
t190 = t200 ^ 2;
t286 = t189 - t190;
t157 = qJD(1) * t170;
t155 = qJDD(1) * t170;
t255 = qJD(3) * t290;
t121 = qJD(4) * t200 - t197 * t255;
t218 = -qJD(4) * t197 - t200 * t255;
t78 = t121 * t316 + t191 * t218;
t272 = pkin(3) * t283;
t90 = pkin(4) * t136 - pkin(8) * t139 + t272;
t273 = t196 * t90 + t199 * t78 + t93 * t281;
t270 = t116 * t301;
t203 = qJD(1) ^ 2;
t268 = t197 * t203 * t200;
t266 = t316 * pkin(3);
t265 = t146 * t281;
t31 = -qJDD(3) * pkin(4) - t33;
t15 = pkin(5) * t58 + qJDD(6) + t31;
t264 = -t15 - t333;
t195 = -qJ(4) - pkin(7);
t262 = pkin(5) * t196 - t195;
t260 = -t196 * t78 + t199 * t90;
t51 = -t196 * t99 + t199 * t93;
t3 = t196 * t46 + t199 * t32 + t79 * t281 - t63 * t282;
t257 = qJD(6) + t331;
t71 = t119 * t191 + t259;
t77 = t121 * t191 - t316 * t218;
t98 = t142 * t191 + t316 * t256;
t252 = t131 * t199;
t249 = t200 * t263;
t248 = -g(2) * t295 + t175 * t339;
t169 = -t266 - pkin(4);
t247 = -pkin(8) * t175 - t342;
t246 = -g(1) * t124 - g(2) * t126;
t125 = -t177 * t294 + t292;
t127 = t176 * t196 + t177 * t291;
t245 = -g(1) * t125 - g(2) * t127;
t242 = -g(2) * t201 + t337;
t241 = qJD(5) * t131 * t167 + t31;
t2 = -qJD(6) * t114 + t3 - t325;
t240 = -t131 * t20 + t2;
t239 = -t178 * t195 - t345;
t237 = -t196 * t22 - t199 * t20;
t236 = t196 * t20 - t199 * t22;
t235 = -t196 * t37 - t199 * t36;
t234 = t196 * t36 - t199 * t37;
t68 = t111 * t316 - t108;
t233 = -t114 * t136 + t224 * t58;
t173 = pkin(4) + t330;
t194 = -qJ(6) - pkin(8);
t232 = t173 * t177 - t175 * t194;
t231 = -qJ(6) * t139 - qJD(6) * t146;
t129 = t284 + t296;
t229 = t94 + (-t282 - t304) * t131;
t226 = t265 + t301;
t62 = -qJD(3) * pkin(4) - t68;
t222 = -t101 * t167 + t131 * t62;
t220 = -qJD(1) * t157 + t244;
t217 = -qJDD(3) * t168 + t157 * t357;
t215 = g(1) * t127 - g(2) * t125 + t199 * t334 - t3;
t118 = qJDD(1) * t352 + t274;
t202 = qJD(3) ^ 2;
t213 = -t168 * t202 - 0.2e1 * t155 - t261;
t212 = qJD(5) * t235 - t4 * t196 + t3 * t199;
t211 = -t101 * t299 - t131 * t226;
t128 = t184 - t297;
t96 = -qJD(3) * t129 - t154 * t197 + t182;
t210 = -t96 * t197 + t95 * t200 + (-t128 * t200 - t129 * t197) * qJD(3);
t209 = t4 + t349;
t205 = t137 * t281 + t196 * t277 + t206;
t161 = g(3) * t293;
t153 = qJDD(3) * t200 - t197 * t202;
t152 = qJDD(3) * t197 + t200 * t202;
t150 = t169 - t330;
t141 = t289 * t199;
t140 = t289 * t196;
t132 = t134 ^ 2;
t113 = t114 ^ 2;
t104 = qJD(3) * t139 + qJDD(3) * t146;
t103 = -qJD(3) * t136 + qJDD(3) * t224;
t76 = pkin(5) * t299 + t98;
t55 = -t113 + t347;
t54 = -pkin(5) * t304 + t71;
t50 = -t101 * t224 + t131 * t136;
t48 = t257 + t62;
t47 = pkin(5) * t226 + t77;
t43 = -qJ(6) * t299 + t52;
t40 = -t205 + t308;
t39 = t223 + t313;
t38 = -pkin(5) * t224 - qJ(6) * t298 + t51;
t29 = -t131 ^ 2 * t199 - t307 - t314;
t28 = t131 * t252 - t307 + t314;
t27 = t229 + t311;
t26 = t229 - t311;
t24 = t114 * t253 - t199 * t58;
t23 = t116 * t252 + t321;
t19 = t114 * t226 + t299 * t58;
t18 = -t116 * t225 + t223 * t298;
t17 = -qJD(5) * t52 + t260;
t16 = -t282 * t99 + t273;
t14 = -qJ(6) * t265 + (-qJD(5) * t99 + t231) * t196 + t273;
t13 = t211 - t233;
t12 = t211 + t233;
t11 = t327 - t348;
t10 = t327 + t348;
t9 = pkin(5) * t136 + t231 * t199 + (-t92 + (qJ(6) * t146 - t93) * t196) * qJD(5) + t260;
t8 = (t223 - t312) * t199 - t354 + t319;
t7 = (-t223 - t312) * t199 + t354 + t319;
t6 = -t270 + (-t321 + (-t306 + t310) * qJD(5)) * t146 + t328;
t5 = t270 + (t321 + (t306 + t310) * qJD(5)) * t146 + t328;
t25 = [0, 0, 0, 0, 0, qJDD(1), t242, g(1) * t201 + g(2) * t198, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t261 + 0.2e1 * t269, -0.2e1 * t192 * t315 + t244, 0 (t242 + (t192 ^ 2 + t193 ^ 2) * t315) * pkin(1), qJDD(1) * t189 + 0.2e1 * t249, 0.2e1 * t197 * t275 - 0.2e1 * t280 * t286, t152, qJDD(1) * t190 - 0.2e1 * t249, t153, 0, t197 * t217 + t200 * t213, -t197 * t213 + t200 * t217 (t189 + t190) * t154 + t210 - t244, t155 * t170 - g(1) * (-pkin(2) * t176 + pkin(7) * t178 - t345) - g(2) * (pkin(2) * t178 + pkin(7) * t176 + t186) + t210 * t168, t137 * t139 + t146 * t351, t320 - t356, t104, -t102 * t224 + t134 * t136, t103, 0, -qJDD(3) * t98 + t102 * t352 - t118 * t224 + t133 * t136 - t261 * t177 + (t134 * t343 - t77) * qJD(3), -t78 * qJD(3) - t99 * qJDD(3) + t118 * t146 + t133 * t139 + t137 * t272 + t351 * t352 - t248, -t99 * t102 - t78 * t134 - t69 * t136 + t77 * t137 - t68 * t139 - t33 * t146 + t224 * t34 + t351 * t98 - t244, t34 * t99 + t69 * t78 - t33 * t98 - t68 * t77 + t118 * t352 + t133 * t272 - g(1) * (-t174 * t176 + t239) - g(2) * (-t176 * t195 + t287) t18, t6, t10, t19, t12, t50, t62 * t301 + t101 * t51 + t114 * t77 + t131 * t17 + t136 * t36 - t224 * t4 + t58 * t98 + (t196 * t31 + t281 * t62) * t146 + t245, t62 * t300 - t101 * t52 + t116 * t77 - t131 * t16 - t136 * t37 + t224 * t3 + t223 * t98 + (t199 * t31 - t282 * t62) * t146 + t246, -t114 * t16 - t116 * t17 - t51 * t223 - t52 * t58 + t235 * t139 + (qJD(5) * t234 - t196 * t3 - t199 * t4) * t146 + t248, t3 * t52 + t37 * t16 + t4 * t51 + t36 * t17 + t31 * t98 + t62 * t77 - g(1) * t239 - g(2) * (pkin(8) * t295 + t178 * t342 + t287) + (-g(1) * (-t174 + t247) + g(2) * t195) * t176, t18, t6, t10, t19, t12, t50, t48 * t301 - t1 * t224 + t101 * t38 + t114 * t47 + t131 * t9 + t136 * t20 + t58 * t76 + (t15 * t196 + t281 * t48) * t146 + t245, t48 * t300 - t101 * t43 + t116 * t47 - t131 * t14 - t136 * t22 + t224 * t2 + t223 * t76 + (t15 * t199 - t282 * t48) * t146 + t246, -t114 * t14 - t116 * t9 - t38 * t223 - t43 * t58 + t237 * t139 + (qJD(5) * t236 - t1 * t199 - t196 * t2) * t146 + t248, t2 * t43 + t22 * t14 + t1 * t38 + t20 * t9 + t15 * t76 + t48 * t47 + pkin(1) * t337 - g(2) * t287 + (-g(1) * t262 - g(2) * t232) * t178 + (-g(1) * (-t174 - t232) - g(2) * t262) * t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t153, -t152, 0, t197 * t95 + t200 * t96 - g(3) + (-t128 * t197 + t129 * t200) * qJD(3), 0, 0, 0, 0, 0, 0, t103, -t104, t320 + t356, -t136 * t68 + t139 * t69 + t146 * t34 + t224 * t33 - g(3), 0, 0, 0, 0, 0, 0, t13, t11, t5, t136 * t62 - t139 * t234 + t146 * t212 - t224 * t31 - g(3), 0, 0, 0, 0, 0, 0, t13, t11, t5, t136 * t48 - t224 * t15 - g(3) - t236 * t139 + (qJD(5) * t237 - t1 * t196 + t199 * t2) * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t268, t286 * t203, t276, t268, t275, qJDD(3), -t332 + t182 + (t129 - t296) * qJD(3) + (t220 + t250) * t197, g(3) * t197 + (t128 + t297) * qJD(3) + t220 * t200 + t267, 0, 0, t303, -t132 + t346 (t163 + t134) * qJD(3) + t214, -t303, -t238, qJDD(3), -t333 + t71 * qJD(3) - t133 * t137 + t227 + (qJDD(3) * t316 - t134 * t285) * pkin(3) + t33, qJD(3) * t72 + t133 * t134 + (-qJDD(3) * t191 - t137 * t285) * pkin(3) - t34 - t216, -t102 * t344 - t351 * t266 - (-t69 + t71) * t137 + (t72 - t68) * t134, t68 * t71 - t69 * t72 + (t316 * t33 - t332 + t191 * t34 + (-qJD(1) * t133 + t244) * t197) * pkin(3), t23, t8, t28, t24, t27, -t305, -t114 * t71 - t131 * t41 - t137 * t36 + t169 * t58 + (-t241 - t333) * t199 + t222 * t196 + t288, -t116 * t71 + t131 * t42 + t137 * t37 + t169 * t223 + t161 + t222 * t199 + (-t227 + t241) * t196, t114 * t42 + t116 * t41 + (-t134 * t36 - t167 * t58 + t3 + (t116 * t167 - t36) * qJD(5)) * t199 + (-t134 * t37 + t167 * t223 - t4 + (t114 * t167 - t37) * qJD(5)) * t196 + t216, t31 * t169 - t37 * t42 - t36 * t41 - t62 * t71 - g(3) * (t185 - t247) + t212 * t167 + t244 * (pkin(4) * t175 - pkin(8) * t177 + t343) t23, t8, t28, t24, t27, -t305, -t101 * t140 - t114 * t54 - t137 * t20 + t150 * t58 + t264 * t199 + t317 * t131 + (t134 * t48 + (t48 + t331) * qJD(5)) * t196 + t288, -t101 * t141 - t116 * t54 + t137 * t22 + t150 * t223 + t161 + t48 * t252 - t318 * t131 + (pkin(5) * t350 + t15 - t227) * t196, -t318 * t114 - t317 * t116 + t140 * t223 - t141 * t58 - t196 * t355 + t240 * t199 + t216, t2 * t141 - t1 * t140 + t15 * t150 - g(3) * (t185 + t232) + (pkin(5) * t282 - t54) * t48 + t318 * t22 + t317 * t20 + t244 * (t173 * t175 + t177 * t194 + t343); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137 * t357 + t238 (t163 - t134) * qJD(3) + t214, -t132 - t346, t134 * t69 + t137 * t68 + t118 + t261, 0, 0, 0, 0, 0, 0, t26, t29, t7, -t137 * t62 + (t4 + t322) * t199 + (t3 - t323) * t196 + t261, 0, 0, 0, 0, 0, 0, t26, t29, t7, -t137 * t48 + t240 * t196 + t199 * t355 + t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t309, t55, t39, -t309, t40, t101, -t116 * t62 + t209 + t322, t114 * t62 + t215 + t323, 0, 0, t309, t55, t39, -t309, t40, t101, 0.2e1 * t341 - t326 + t324 + (-t257 - t48) * t116 + t209, -pkin(5) * t347 + t325 + t131 * t21 + (qJD(6) + t48) * t114 + t215, -pkin(5) * t223 - t114 * t329, t329 * t22 + (-t48 * t116 + t1 + t349) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205 + t308, t223 - t313, -t113 - t347, t114 * t22 + t116 * t20 - t227 - t264;];
tau_reg  = t25;

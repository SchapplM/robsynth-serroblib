% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:14:01
% EndTime: 2019-03-09 02:14:07
% DurationCPUTime: 4.02s
% Computational Cost: add. (6760->470), mult. (13563->552), div. (0->0), fcn. (9573->10), ass. (0->249)
t181 = sin(qJ(5));
t184 = cos(qJ(5));
t176 = sin(pkin(9));
t180 = -pkin(1) - qJ(3);
t144 = t180 * qJD(1) + qJD(2);
t232 = -pkin(7) * qJD(1) + t144;
t104 = t232 * t176;
t177 = cos(pkin(9));
t105 = t232 * t177;
t182 = sin(qJ(4));
t185 = cos(qJ(4));
t69 = t104 * t185 + t105 * t182;
t66 = qJD(4) * pkin(8) + t69;
t124 = t176 * t185 + t177 * t182;
t115 = t124 * qJD(1);
t269 = t177 * t185;
t240 = qJD(1) * t269;
t257 = qJD(1) * t176;
t117 = -t182 * t257 + t240;
t158 = qJD(1) * qJ(2) + qJD(3);
t131 = pkin(3) * t257 + t158;
t67 = pkin(4) * t115 - pkin(8) * t117 + t131;
t38 = t181 * t67 + t184 * t66;
t251 = t184 * qJD(4);
t94 = t117 * t181 - t251;
t28 = -qJ(6) * t94 + t38;
t325 = qJD(5) + t115;
t330 = t325 * t28;
t37 = -t181 * t66 + t184 * t67;
t329 = t325 * t37;
t328 = t325 * t38;
t172 = pkin(9) + qJ(4);
t159 = sin(t172);
t160 = cos(t172);
t304 = g(3) * t160;
t183 = sin(qJ(1));
t186 = cos(qJ(1));
t320 = g(1) * t183 - g(2) * t186;
t199 = t159 * t320 + t304;
t125 = -t182 * t176 + t269;
t229 = t184 * t325;
t248 = qJD(1) * qJD(4);
t238 = t182 * t248;
t195 = t124 * qJDD(1) - t176 * t238;
t237 = t185 * t248;
t84 = t177 * t237 + t195;
t81 = qJDD(5) + t84;
t291 = t181 * t81;
t326 = -t229 * t325 - t291;
t170 = t176 ^ 2;
t171 = t177 ^ 2;
t258 = t170 + t171;
t323 = t144 * t258;
t301 = -pkin(7) + t180;
t127 = t301 * t176;
t128 = t301 * t177;
t322 = -t127 * t182 + t185 * t128;
t164 = t176 * pkin(3);
t179 = -pkin(7) - qJ(3);
t321 = t186 * t164 + t183 * t179;
t173 = qJDD(1) * qJ(2);
t174 = qJD(1) * qJD(2);
t319 = t173 + t174;
t246 = t177 * qJDD(1);
t247 = t176 * qJDD(1);
t197 = -t176 * t237 - t177 * t238 - t182 * t247 + t185 * t246;
t96 = qJD(4) * t181 + t117 * t184;
t281 = qJD(5) * t96;
t55 = -t184 * qJDD(4) + t181 * t197 + t281;
t134 = qJDD(3) + t319;
t222 = g(1) * t186 + g(2) * t183;
t196 = -t222 + t134;
t264 = t184 * t186;
t268 = t181 * t183;
t108 = -t159 * t268 + t264;
t265 = t183 * t184;
t267 = t181 * t186;
t110 = t159 * t267 + t265;
t318 = -g(1) * t108 - g(2) * t110 + t181 * t304;
t253 = qJD(5) * t181;
t54 = -qJD(5) * t251 - t181 * qJDD(4) + t117 * t253 - t184 * t197;
t297 = qJ(6) * t54;
t310 = pkin(5) * t81;
t314 = -qJD(1) * qJD(3) + qJDD(1) * t180;
t126 = qJDD(2) + t314;
t230 = -pkin(7) * qJDD(1) + t126;
t101 = t230 * t176;
t102 = t230 * t177;
t254 = qJD(4) * t185;
t242 = -t185 * t101 - t182 * t102 - t105 * t254;
t255 = qJD(4) * t182;
t39 = -t104 * t255 - t242;
t35 = qJDD(4) * pkin(8) + t39;
t153 = pkin(3) * t247;
t122 = t153 + t134;
t47 = t84 * pkin(4) - pkin(8) * t197 + t122;
t6 = -qJD(5) * t38 - t181 * t35 + t184 * t47;
t1 = -qJD(6) * t96 + t297 + t310 + t6;
t317 = t1 + t330;
t316 = t6 + t328;
t119 = t124 * qJD(4);
t315 = -t117 * t119 + t125 * t197;
t305 = g(3) * t159;
t194 = -t320 * t160 + t305;
t313 = t96 ^ 2;
t312 = t117 ^ 2;
t311 = 0.2e1 * t174;
t309 = pkin(5) * t94;
t308 = pkin(5) * t181;
t303 = g(3) * t184;
t302 = t96 * t94;
t300 = qJ(6) + pkin(8);
t27 = -qJ(6) * t96 + t37;
t22 = pkin(5) * t325 + t27;
t299 = -t27 + t22;
t252 = qJD(5) * t184;
t298 = -t181 * t55 - t94 * t252;
t103 = t182 * t104;
t68 = t105 * t185 - t103;
t83 = pkin(4) * t117 + pkin(8) * t115;
t46 = t181 * t83 + t184 * t68;
t151 = qJ(2) + t164;
t82 = pkin(4) * t124 - pkin(8) * t125 + t151;
t87 = t127 * t185 + t128 * t182;
t85 = t184 * t87;
t52 = t181 * t82 + t85;
t296 = qJ(6) * t55;
t295 = t325 * t94;
t294 = t117 * t94;
t293 = t117 * t96;
t292 = t181 * t54;
t290 = t181 * t94;
t289 = t181 * t96;
t288 = t184 * t55;
t74 = t184 * t81;
t287 = t184 * t94;
t286 = t184 * t96;
t285 = t96 * t325;
t233 = qJD(5) * t300;
t279 = t115 * t181;
t284 = -qJ(6) * t279 + qJD(6) * t184 - t181 * t233 - t46;
t45 = -t181 * t68 + t184 * t83;
t283 = -pkin(5) * t117 - qJD(6) * t181 - t45 + (-qJ(6) * t115 - t233) * t184;
t282 = pkin(1) * qJDD(1);
t280 = t325 * t117;
t278 = t117 * t115;
t276 = t119 * t181;
t275 = t119 * t184;
t274 = t125 * t181;
t273 = t125 * t184;
t271 = t160 * t183;
t270 = t160 * t186;
t263 = -t119 * qJD(4) + t125 * qJDD(4);
t243 = g(2) * t270;
t262 = t159 * t303 + t184 * t243;
t261 = g(1) * t270 + g(2) * t271;
t260 = t186 * pkin(1) + t183 * qJ(2);
t256 = qJD(4) * t115;
t61 = -t124 * qJD(3) + t322 * qJD(4);
t120 = -t176 * t255 + t177 * t254;
t79 = pkin(4) * t120 + pkin(8) * t119 + qJD(2);
t245 = t181 * t79 + t184 * t61 + t82 * t252;
t244 = g(1) * t271;
t241 = t183 * t164 + t260;
t239 = t125 * t252;
t166 = t186 * qJ(2);
t236 = -pkin(1) * t183 + t166;
t235 = -qJD(6) - t309;
t234 = -t181 * t61 + t184 * t79;
t51 = -t181 * t87 + t184 * t82;
t5 = t181 * t47 + t184 * t35 + t67 * t252 - t66 * t253;
t40 = -t182 * t101 + t185 * t102 - t104 * t254 - t105 * t255;
t231 = t258 * t126;
t228 = qJDD(2) - t282;
t227 = qJD(5) * t124 + qJD(1);
t36 = -qJDD(4) * pkin(4) - t40;
t226 = -pkin(8) * qJD(5) * t325 - t36;
t225 = pkin(4) * t159 - pkin(8) * t160;
t224 = g(1) * t110 - g(2) * t108;
t109 = t159 * t265 + t267;
t111 = t159 * t264 - t268;
t223 = -g(1) * t111 - g(2) * t109;
t3 = -qJD(6) * t94 - t296 + t5;
t220 = -t22 * t325 + t3;
t219 = t5 - t329;
t218 = t181 * t28 + t184 * t22;
t217 = t181 * t22 - t184 * t28;
t216 = t181 * t38 + t184 * t37;
t215 = t181 * t37 - t184 * t38;
t214 = -t287 + t289;
t213 = t287 + t289;
t212 = t286 + t290;
t211 = t120 * t115 + t124 * t84;
t157 = pkin(5) * t184 + pkin(4);
t209 = t157 * t159 - t160 * t300;
t208 = qJ(6) * t119 - qJD(6) * t125;
t207 = t236 + t321;
t206 = t74 + (-t253 - t279) * t325;
t205 = -qJD(4) * t120 - qJDD(4) * t124;
t204 = -t243 - t305;
t203 = -t179 * t186 + t241;
t202 = -t184 * t54 - t96 * t253;
t65 = -qJD(4) * pkin(4) - t68;
t201 = -pkin(8) * t81 + t325 * t65;
t200 = t239 - t276;
t193 = g(1) * t109 - g(2) * t111 + t160 * t303 - t5;
t15 = pkin(5) * t55 + qJDD(6) + t36;
t191 = t196 + t319;
t190 = -t216 * qJD(5) - t6 * t181 + t5 * t184;
t189 = t68 * t119 - t69 * t120 - t39 * t124 - t40 * t125 + t320;
t62 = qJD(3) * t125 + qJD(4) * t87;
t188 = t6 + t318;
t187 = qJD(1) ^ 2;
t133 = t300 * t184;
t132 = t300 * t181;
t130 = t181 * t244;
t114 = t115 ^ 2;
t93 = t94 ^ 2;
t63 = pkin(5) * t274 - t322;
t56 = -pkin(5) * t279 + t69;
t50 = -t235 + t65;
t49 = -t93 + t313;
t48 = t120 * t325 + t124 * t81;
t44 = pkin(5) * t200 + t62;
t41 = -qJ(6) * t274 + t52;
t33 = pkin(5) * t124 - qJ(6) * t273 + t51;
t31 = t285 - t55;
t30 = -t54 + t295;
t26 = -t293 + t326;
t25 = -t293 - t326;
t24 = t206 + t294;
t23 = t206 - t294;
t21 = t290 * t325 - t288;
t20 = t96 * t229 - t292;
t19 = -t52 * qJD(5) + t234;
t18 = -t87 * t253 + t245;
t17 = -t298 * t125 - t94 * t276;
t16 = t125 * t202 - t96 * t275;
t14 = -qJ(6) * t239 + (-qJD(5) * t87 + t208) * t181 + t245;
t13 = pkin(5) * t120 + t208 * t184 + (-t85 + (qJ(6) * t125 - t82) * t181) * qJD(5) + t234;
t12 = -t120 * t94 - t124 * t55 - t200 * t325 - t81 * t274;
t11 = t81 * t273 + t120 * t96 - t124 * t54 + (-t125 * t253 - t275) * t325;
t10 = -t124 * t291 + t119 * t94 - t125 * t55 + (-t120 * t181 - t227 * t184) * t325;
t9 = -t124 * t74 + t119 * t96 + t125 * t54 + (-t120 * t184 + t227 * t181) * t325;
t8 = -t213 * t115 + t202 + t298;
t7 = t214 * t115 - t202 + t298;
t4 = t213 * t119 + (t292 - t288 + (-t286 + t290) * qJD(5)) * t125;
t2 = t214 * t120 + t212 * qJD(1) + (t212 * qJD(5) - t288 - t292) * t124;
t29 = [0, 0, 0, 0, 0, qJDD(1), t320, t222, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t320 - 0.2e1 * t282, 0.2e1 * t173 + t311 - t222, -t228 * pkin(1) - g(1) * t236 - g(2) * t260 + (t173 + t311) * qJ(2), t171 * qJDD(1), -0.2e1 * t176 * t246, 0, t170 * qJDD(1), 0, 0, t191 * t176, t191 * t177, t320 + t258 * (-t126 - t314) t134 * qJ(2) + t158 * qJD(2) - g(1) * (t180 * t183 + t166) - g(2) * (qJ(3) * t186 + t260) + t180 * t231 - qJD(3) * t323, t315, t119 * t115 - t117 * t120 - t124 * t197 - t125 * t84, t263, t211, t205, 0, qJD(2) * t115 - qJD(4) * t62 + qJDD(4) * t322 + t120 * t131 + t122 * t124 + t151 * t84 - t159 * t222, qJD(2) * t117 - t61 * qJD(4) - t87 * qJDD(4) - t131 * t119 + t122 * t125 + t151 * t197 - t261, -t61 * t115 + t62 * t117 - t197 * t322 - t87 * t84 + t189, -g(1) * t207 - g(2) * t203 + t131 * qJD(2) + t122 * t151 + t322 * t40 + t39 * t87 + t69 * t61 - t68 * t62, t16, t4, t11, t17, t12, t48, -t65 * t276 + t325 * t19 + t120 * t37 + t124 * t6 + t51 * t81 - t55 * t322 + t62 * t94 + (t181 * t36 + t252 * t65) * t125 + t223, -t65 * t275 - t325 * t18 - t120 * t38 - t124 * t5 - t52 * t81 + t54 * t322 + t62 * t96 + (t184 * t36 - t253 * t65) * t125 + t224, -t18 * t94 - t19 * t96 + t51 * t54 - t52 * t55 + t216 * t119 + (qJD(5) * t215 - t181 * t5 - t184 * t6) * t125 + t261, t5 * t52 + t38 * t18 + t6 * t51 + t37 * t19 - t36 * t322 + t65 * t62 - g(1) * (t186 * t225 + t207) - g(2) * (t183 * t225 + t203) t16, t4, t11, t17, t12, t48, -t50 * t276 + t1 * t124 + t325 * t13 + t120 * t22 + t33 * t81 + t44 * t94 + t55 * t63 + (t15 * t181 + t252 * t50) * t125 + t223, -t50 * t275 - t325 * t14 - t120 * t28 - t124 * t3 - t41 * t81 + t44 * t96 - t54 * t63 + (t15 * t184 - t253 * t50) * t125 + t224, -t13 * t96 - t14 * t94 + t33 * t54 - t41 * t55 + t218 * t119 + (qJD(5) * t217 - t1 * t184 - t181 * t3) * t125 + t261, t3 * t41 + t28 * t14 + t1 * t33 + t22 * t13 + t15 * t63 + t50 * t44 - g(1) * (t166 + t321) - g(2) * t241 + (-g(1) * t209 - g(2) * (-t179 + t308)) * t186 + (-g(1) * (-pkin(1) - t308) - g(2) * t209) * t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t187, -qJ(2) * t187 + t228 - t320, 0, 0, 0, 0, 0, 0, -t187 * t176, -t187 * t177, -t258 * qJDD(1), -qJD(1) * t158 + t231 - t320, 0, 0, 0, 0, 0, 0, -qJD(1) * t115 + t263, -qJD(1) * t117 + t205, -t211 - t315, -qJD(1) * t131 - t189, 0, 0, 0, 0, 0, 0, t10, t9, t2, -qJD(1) * t216 + t119 * t65 - t120 * t215 + t124 * t190 - t125 * t36 - t320, 0, 0, 0, 0, 0, 0, t10, t9, t2, t119 * t50 - t125 * t15 - t217 * t120 - t218 * qJD(1) + (-qJD(5) * t218 - t1 * t181 + t184 * t3) * t124 - t320; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t247, t246, -t258 * t187, qJD(1) * t323 + t196, 0, 0, 0, 0, 0, 0 (t117 + t240) * qJD(4) + t195, t197 - t256, -t114 - t312, t115 * t69 + t117 * t68 + t153 + t196, 0, 0, 0, 0, 0, 0, t23, t26, t7, -t117 * t65 + t219 * t181 + t316 * t184 - t222, 0, 0, 0, 0, 0, 0, t23, t26, t7, -t117 * t50 + t220 * t181 + t317 * t184 - t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, -t114 + t312, t197 + t256, -t278 (t117 - t240) * qJD(4) - t195, qJDD(4), qJD(4) * t69 - t117 * t131 + t194 + t40, t115 * t131 + (t68 + t103) * qJD(4) + t199 + t242, 0, 0, t20, t8, t25, t21, t24, -t280, -pkin(4) * t55 - t325 * t45 - t117 * t37 - t69 * t94 + (t226 - t244) * t184 + t201 * t181 + t262, pkin(4) * t54 + t325 * t46 + t117 * t38 - t69 * t96 + t130 + t201 * t184 + (t204 - t226) * t181, t45 * t96 + t46 * t94 + ((-t55 + t281) * pkin(8) + t219) * t184 + ((qJD(5) * t94 - t54) * pkin(8) - t316) * t181 - t199, -t37 * t45 - t38 * t46 - t65 * t69 + (-t36 + t194) * pkin(4) + (t190 - t199) * pkin(8), t20, t8, t25, t21, t24, -t280, -t117 * t22 - t132 * t81 - t157 * t55 - t56 * t94 + (-t15 - t244) * t184 + t283 * t325 + (t115 * t50 + (t50 + t309) * qJD(5)) * t181 + t262, t117 * t28 - t133 * t81 + t157 * t54 - t56 * t96 + t130 + t50 * t229 - t284 * t325 + (pkin(5) * t281 + t15 + t204) * t181, -t132 * t54 - t133 * t55 - t317 * t181 + t220 * t184 - t283 * t96 - t284 * t94 - t199, t3 * t133 - t1 * t132 - t15 * t157 + g(3) * t209 + (pkin(5) * t253 - t56) * t50 + t284 * t28 + t283 * t22 - t320 * (t157 * t160 + t159 * t300); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t49, t30, -t302, t31, t81, -t65 * t96 + t188 + t328, t65 * t94 + t193 + t329, 0, 0, t302, t49, t30, -t302, t31, t81, 0.2e1 * t310 + t297 + t330 + (t235 - t50) * t96 + t188, -pkin(5) * t313 + t296 + t325 * t27 + (qJD(6) + t50) * t94 + t193, t54 * pkin(5) - t299 * t94, t299 * t28 + (-t50 * t96 + t1 + t318) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 + t285, -t54 - t295, -t93 - t313, t22 * t96 + t28 * t94 + t15 - t194;];
tau_reg  = t29;

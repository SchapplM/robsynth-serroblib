% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR2
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
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:31
% EndTime: 2019-12-05 18:28:42
% DurationCPUTime: 5.73s
% Computational Cost: add. (33736->448), mult. (81180->644), div. (0->0), fcn. (59383->10), ass. (0->273)
t315 = -2 * qJD(3);
t239 = sin(pkin(9));
t240 = cos(pkin(9));
t247 = cos(qJ(2));
t243 = sin(qJ(2));
t277 = qJD(1) * t243;
t211 = -qJD(1) * t240 * t247 + t239 * t277;
t213 = (t239 * t247 + t240 * t243) * qJD(1);
t196 = t213 * t211;
t306 = qJDD(2) - t196;
t314 = t239 * t306;
t313 = t240 * t306;
t241 = sin(qJ(5));
t242 = sin(qJ(4));
t246 = cos(qJ(4));
t187 = t211 * t246 + t213 * t242;
t189 = -t211 * t242 + t213 * t246;
t245 = cos(qJ(5));
t156 = t187 * t245 + t189 * t241;
t158 = -t187 * t241 + t189 * t245;
t120 = t158 * t156;
t235 = qJDD(2) + qJDD(4);
t227 = qJDD(5) + t235;
t308 = -t120 + t227;
t312 = t241 * t308;
t164 = t189 * t187;
t307 = -t164 + t235;
t311 = t242 * t307;
t310 = t245 * t308;
t309 = t246 * t307;
t229 = t243 * qJDD(1);
t269 = qJD(1) * qJD(2);
t267 = t247 * t269;
t219 = t229 + t267;
t230 = t247 * qJDD(1);
t268 = t243 * t269;
t220 = t230 - t268;
t197 = -t219 * t239 + t220 * t240;
t198 = t219 * t240 + t220 * t239;
t256 = t242 * t197 + t246 * t198;
t144 = -qJD(4) * t187 + t256;
t236 = qJD(2) + qJD(4);
t183 = t236 * t187;
t130 = t144 + t183;
t208 = qJD(2) * t211;
t179 = t198 + t208;
t249 = qJD(1) ^ 2;
t283 = t243 * t249;
t244 = sin(qJ(1));
t305 = cos(qJ(1));
t255 = g(1) * t305 + g(2) * t244;
t298 = qJDD(1) * pkin(6);
t216 = -pkin(1) * t249 - t255 + t298;
t285 = t243 * t216;
t172 = qJDD(2) * pkin(2) - t219 * qJ(3) - t285 + (pkin(2) * t283 + qJ(3) * t269 - g(3)) * t247;
t201 = -t243 * g(3) + t247 * t216;
t238 = t247 ^ 2;
t232 = t238 * t249;
t252 = qJD(2) * pkin(2) - qJ(3) * t277;
t173 = -pkin(2) * t232 + t220 * qJ(3) - qJD(2) * t252 + t201;
t259 = t172 * t240 - t239 * t173 + t213 * t315;
t136 = t172 * t239 + t173 * t240 + t211 * t315;
t154 = t156 ^ 2;
t155 = t158 ^ 2;
t185 = t187 ^ 2;
t186 = t189 ^ 2;
t209 = t211 ^ 2;
t210 = t213 ^ 2;
t228 = qJD(5) + t236;
t226 = t228 ^ 2;
t234 = t236 ^ 2;
t109 = pkin(3) * t306 - pkin(7) * t179 + t259;
t258 = qJD(2) * pkin(3) - pkin(7) * t213;
t110 = -t209 * pkin(3) + t197 * pkin(7) - qJD(2) * t258 + t136;
t67 = -t109 * t246 + t242 * t110;
t68 = t242 * t109 + t246 * t110;
t39 = t242 * t68 - t246 * t67;
t304 = t239 * t39;
t303 = t240 * t39;
t266 = t244 * g(1) - g(2) * t305;
t253 = qJDD(1) * pkin(1) + t266;
t174 = t220 * pkin(2) - qJDD(3) - t252 * t277 + (qJ(3) * t238 + pkin(6)) * t249 + t253;
t134 = t197 * pkin(3) + t209 * pkin(7) - t213 * t258 + t174;
t263 = -t197 * t246 + t242 * t198;
t143 = -qJD(4) * t189 - t263;
t260 = pkin(4) * t236 - pkin(8) * t189;
t82 = t143 * pkin(4) + t185 * pkin(8) - t189 * t260 + t134;
t302 = t241 * t82;
t54 = pkin(4) * t307 - pkin(8) * t130 - t67;
t56 = -t185 * pkin(4) + t143 * pkin(8) - t236 * t260 + t68;
t30 = t241 * t56 - t245 * t54;
t31 = t241 * t54 + t245 * t56;
t17 = t241 * t31 - t245 * t30;
t301 = t242 * t17;
t300 = t245 * t82;
t299 = t246 * t17;
t297 = t228 * t241;
t296 = t228 * t245;
t295 = t236 * t242;
t294 = t236 * t246;
t293 = t239 * t174;
t193 = qJDD(2) + t196;
t292 = t239 * t193;
t291 = t240 * t174;
t290 = t240 * t193;
t115 = t120 + t227;
t289 = t241 * t115;
t288 = t242 * t134;
t161 = t164 + t235;
t287 = t242 * t161;
t101 = t136 * t239 + t240 * t259;
t286 = t243 * t101;
t225 = t247 * t283;
t284 = t243 * (qJDD(2) + t225);
t282 = t245 * t115;
t281 = t246 * t134;
t280 = t246 * t161;
t279 = t247 * (qJDD(2) - t225);
t276 = qJD(2) * t213;
t275 = qJD(2) * t239;
t274 = qJD(2) * t240;
t271 = qJD(4) + t236;
t270 = qJD(5) + t228;
t18 = t241 * t30 + t245 * t31;
t40 = t242 * t67 + t246 * t68;
t102 = t136 * t240 - t239 * t259;
t265 = -t143 * t245 + t241 * t144;
t200 = t247 * g(3) + t285;
t262 = t243 * t200 + t201 * t247;
t113 = -t226 - t154;
t84 = t113 * t241 + t310;
t261 = pkin(4) * t84 - t30;
t257 = t241 * t143 + t245 * t144;
t147 = -t155 - t226;
t98 = t147 * t245 - t289;
t254 = pkin(4) * t98 - t31;
t177 = t197 + t276;
t251 = (-qJD(5) + t228) * t158 - t265;
t250 = (-qJD(4) + t236) * t189 - t263;
t90 = -qJD(5) * t156 + t257;
t248 = qJD(2) ^ 2;
t237 = t243 ^ 2;
t231 = t237 * t249;
t221 = t230 - 0.2e1 * t268;
t218 = t229 + 0.2e1 * t267;
t215 = pkin(6) * t249 + t253;
t204 = -t210 - t248;
t203 = -t210 + t248;
t202 = t209 - t248;
t191 = -t248 - t209;
t182 = -t186 + t234;
t181 = t185 - t234;
t180 = -t186 - t234;
t178 = t198 - t208;
t176 = -t197 + t276;
t175 = -t209 - t210;
t168 = -t204 * t239 - t290;
t167 = t204 * t240 - t292;
t166 = t191 * t240 - t314;
t165 = t191 * t239 + t313;
t163 = t186 - t185;
t159 = -t234 - t185;
t152 = t228 * t156;
t151 = -t155 + t226;
t150 = t154 - t226;
t149 = (-t187 * t246 + t189 * t242) * t236;
t148 = (-t187 * t242 - t189 * t246) * t236;
t146 = t177 * t240 + t179 * t239;
t145 = t177 * t239 - t179 * t240;
t141 = -t185 - t186;
t140 = t181 * t246 - t287;
t139 = -t182 * t242 + t309;
t138 = t181 * t242 + t280;
t137 = t182 * t246 + t311;
t133 = -t180 * t242 - t280;
t132 = t180 * t246 - t287;
t129 = t144 - t183;
t128 = -t187 * t271 + t256;
t125 = t189 * t271 + t263;
t124 = t144 * t246 - t189 * t295;
t123 = t144 * t242 + t189 * t294;
t122 = -t143 * t242 + t187 * t294;
t121 = t143 * t246 + t187 * t295;
t119 = t155 - t154;
t118 = t159 * t246 - t311;
t117 = t159 * t242 + t309;
t112 = (-t156 * t245 + t158 * t241) * t228;
t111 = (-t156 * t241 - t158 * t245) * t228;
t107 = -t154 - t155;
t106 = t150 * t245 - t289;
t105 = -t151 * t241 + t310;
t104 = t150 * t241 + t282;
t103 = t151 * t245 + t312;
t100 = -pkin(7) * t132 - t281;
t99 = -t147 * t241 - t282;
t96 = -t132 * t239 + t133 * t240;
t95 = t132 * t240 + t133 * t239;
t94 = t130 * t242 + t246 * t250;
t93 = -t125 * t246 - t129 * t242;
t92 = -t130 * t246 + t242 * t250;
t91 = -t125 * t242 + t129 * t246;
t89 = -qJD(5) * t158 - t265;
t88 = -pkin(7) * t117 - t288;
t87 = -t117 * t239 + t118 * t240;
t86 = t117 * t240 + t118 * t239;
t85 = t113 * t245 - t312;
t81 = -t111 * t242 + t112 * t246;
t80 = t111 * t246 + t112 * t242;
t79 = -t156 * t270 + t257;
t78 = t152 + t90;
t77 = -t152 + t90;
t74 = t158 * t270 + t265;
t73 = -t158 * t297 + t245 * t90;
t72 = t158 * t296 + t241 * t90;
t71 = t156 * t296 - t241 * t89;
t70 = t156 * t297 + t245 * t89;
t69 = -pkin(3) * t128 + pkin(7) * t133 - t288;
t65 = -pkin(3) * t125 + pkin(7) * t118 + t281;
t64 = -t104 * t242 + t106 * t246;
t63 = -t103 * t242 + t105 * t246;
t62 = t104 * t246 + t106 * t242;
t61 = t103 * t246 + t105 * t242;
t60 = -t242 * t98 + t246 * t99;
t59 = t242 * t99 + t246 * t98;
t58 = -t239 * t92 + t240 * t94;
t57 = t239 * t94 + t240 * t92;
t55 = -pkin(8) * t98 - t300;
t52 = -t242 * t84 + t246 * t85;
t51 = t242 * t85 + t246 * t84;
t50 = -pkin(8) * t84 - t302;
t49 = t241 * t78 + t245 * t251;
t48 = -t241 * t77 - t245 * t74;
t47 = t241 * t251 - t245 * t78;
t46 = -t241 * t74 + t245 * t77;
t45 = pkin(4) * t47;
t44 = -t242 * t72 + t246 * t73;
t43 = -t242 * t70 + t246 * t71;
t42 = t242 * t73 + t246 * t72;
t41 = t242 * t71 + t246 * t70;
t38 = pkin(3) * t134 + pkin(7) * t40;
t37 = -pkin(4) * t79 + pkin(8) * t99 - t302;
t36 = -t239 * t59 + t240 * t60;
t35 = t239 * t60 + t240 * t59;
t34 = -pkin(4) * t74 + pkin(8) * t85 + t300;
t33 = -pkin(7) * t92 - t39;
t32 = -pkin(3) * t141 + pkin(7) * t94 + t40;
t28 = -t239 * t51 + t240 * t52;
t27 = t239 * t52 + t240 * t51;
t26 = -t242 * t47 + t246 * t49;
t25 = -t242 * t46 + t246 * t48;
t24 = t242 * t49 + t246 * t47;
t23 = t242 * t48 + t246 * t46;
t22 = t240 * t40 - t304;
t21 = t239 * t40 + t303;
t20 = -pkin(7) * t59 - t242 * t37 + t246 * t55;
t19 = -pkin(7) * t51 - t242 * t34 + t246 * t50;
t16 = pkin(4) * t17;
t15 = -pkin(3) * t79 + pkin(7) * t60 + t242 * t55 + t246 * t37;
t14 = pkin(4) * t82 + pkin(8) * t18;
t13 = -pkin(3) * t74 + pkin(7) * t52 + t242 * t50 + t246 * t34;
t12 = -t239 * t24 + t240 * t26;
t11 = t239 * t26 + t24 * t240;
t10 = -pkin(8) * t47 - t17;
t9 = -pkin(4) * t107 + pkin(8) * t49 + t18;
t8 = t18 * t246 - t301;
t7 = t18 * t242 + t299;
t6 = -pkin(7) * t24 + t10 * t246 - t242 * t9;
t5 = -pkin(3) * t107 + pkin(7) * t26 + t10 * t242 + t246 * t9;
t4 = -t239 * t7 + t240 * t8;
t3 = t239 * t8 + t240 * t7;
t2 = -pkin(7) * t7 - pkin(8) * t299 - t14 * t242;
t1 = pkin(3) * t82 + pkin(7) * t8 - pkin(8) * t301 + t14 * t246;
t29 = [0, 0, 0, 0, 0, qJDD(1), t266, t255, 0, 0, (t219 + t267) * t243, t218 * t247 + t221 * t243, t284 + t247 * (-t231 + t248), (t220 - t268) * t247, t243 * (t232 - t248) + t279, 0, t247 * t215 + pkin(1) * t221 + pkin(6) * (t247 * (-t232 - t248) - t284), -t243 * t215 - pkin(1) * t218 + pkin(6) * (-t279 - t243 * (-t231 - t248)), pkin(1) * (t231 + t232) + (t237 + t238) * t298 + t262, pkin(1) * t215 + pkin(6) * t262, t243 * (t198 * t240 - t213 * t275) + t247 * (t198 * t239 + t213 * t274), t243 * (-t176 * t240 - t178 * t239) + t247 * (-t176 * t239 + t178 * t240), t243 * (-t203 * t239 + t313) + t247 * (t203 * t240 + t314), t243 * (-t197 * t239 + t211 * t274) + t247 * (t197 * t240 + t211 * t275), t243 * (t202 * t240 - t292) + t247 * (t202 * t239 + t290), (t243 * (-t211 * t240 + t213 * t239) + t247 * (-t211 * t239 - t213 * t240)) * qJD(2), t243 * (-qJ(3) * t165 - t293) + t247 * (-pkin(2) * t176 + qJ(3) * t166 + t291) - pkin(1) * t176 + pkin(6) * (-t165 * t243 + t166 * t247), t243 * (-qJ(3) * t167 - t291) + t247 * (-pkin(2) * t178 + qJ(3) * t168 - t293) - pkin(1) * t178 + pkin(6) * (-t167 * t243 + t168 * t247), t243 * (-qJ(3) * t145 - t101) + t247 * (-pkin(2) * t175 + qJ(3) * t146 + t102) - pkin(1) * t175 + pkin(6) * (-t145 * t243 + t146 * t247), -qJ(3) * t286 + t247 * (pkin(2) * t174 + qJ(3) * t102) + pkin(1) * t174 + pkin(6) * (t102 * t247 - t286), t243 * (-t123 * t239 + t124 * t240) + t247 * (t123 * t240 + t124 * t239), t243 * (-t239 * t91 + t240 * t93) + t247 * (t239 * t93 + t240 * t91), t243 * (-t137 * t239 + t139 * t240) + t247 * (t137 * t240 + t139 * t239), t243 * (-t121 * t239 + t122 * t240) + t247 * (t121 * t240 + t122 * t239), t243 * (-t138 * t239 + t140 * t240) + t247 * (t138 * t240 + t140 * t239), t243 * (-t148 * t239 + t149 * t240) + t247 * (t148 * t240 + t149 * t239), t243 * (-qJ(3) * t86 - t239 * t65 + t240 * t88) + t247 * (-pkin(2) * t125 + qJ(3) * t87 + t239 * t88 + t240 * t65) - pkin(1) * t125 + pkin(6) * (-t243 * t86 + t247 * t87), t243 * (-qJ(3) * t95 + t100 * t240 - t239 * t69) + t247 * (-pkin(2) * t128 + qJ(3) * t96 + t100 * t239 + t240 * t69) - pkin(1) * t128 + pkin(6) * (-t243 * t95 + t247 * t96), t243 * (-qJ(3) * t57 - t239 * t32 + t240 * t33) + t247 * (-pkin(2) * t141 + qJ(3) * t58 + t239 * t33 + t240 * t32) - pkin(1) * t141 + pkin(6) * (-t243 * t57 + t247 * t58), t243 * (-pkin(7) * t303 - qJ(3) * t21 - t239 * t38) + t247 * (pkin(2) * t134 - pkin(7) * t304 + qJ(3) * t22 + t240 * t38) + pkin(1) * t134 + pkin(6) * (-t21 * t243 + t22 * t247), t243 * (-t239 * t42 + t240 * t44) + t247 * (t239 * t44 + t240 * t42), t243 * (-t23 * t239 + t240 * t25) + t247 * (t23 * t240 + t239 * t25), t243 * (-t239 * t61 + t240 * t63) + t247 * (t239 * t63 + t240 * t61), t243 * (-t239 * t41 + t240 * t43) + t247 * (t239 * t43 + t240 * t41), t243 * (-t239 * t62 + t240 * t64) + t247 * (t239 * t64 + t240 * t62), t243 * (-t239 * t80 + t240 * t81) + t247 * (t239 * t81 + t240 * t80), t243 * (-qJ(3) * t27 - t13 * t239 + t19 * t240) + t247 * (-pkin(2) * t74 + qJ(3) * t28 + t13 * t240 + t19 * t239) - pkin(1) * t74 + pkin(6) * (-t243 * t27 + t247 * t28), t243 * (-qJ(3) * t35 - t15 * t239 + t20 * t240) + t247 * (-pkin(2) * t79 + qJ(3) * t36 + t15 * t240 + t20 * t239) - pkin(1) * t79 + pkin(6) * (-t243 * t35 + t247 * t36), t243 * (-qJ(3) * t11 - t239 * t5 + t240 * t6) + t247 * (-pkin(2) * t107 + qJ(3) * t12 + t239 * t6 + t240 * t5) - pkin(1) * t107 + pkin(6) * (-t11 * t243 + t12 * t247), t243 * (-qJ(3) * t3 - t1 * t239 + t2 * t240) + t247 * (pkin(2) * t82 + qJ(3) * t4 + t1 * t240 + t2 * t239) + pkin(1) * t82 + pkin(6) * (-t243 * t3 + t247 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225, t231 - t232, t229, t225, t230, qJDD(2), -t200, -t201, 0, 0, t196, t210 - t209, t179, -t196, t177, qJDD(2), pkin(2) * t165 + t259, pkin(2) * t167 - t136, pkin(2) * t145, pkin(2) * t101, t164, t163, t130, -t164, t250, t235, pkin(2) * t86 + pkin(3) * t117 - t67, pkin(2) * t95 + pkin(3) * t132 - t68, pkin(2) * t57 + pkin(3) * t92, pkin(2) * t21 + pkin(3) * t39, t120, t119, t78, -t120, t251, t227, pkin(2) * t27 + pkin(3) * t51 + t261, pkin(2) * t35 + pkin(3) * t59 + t254, pkin(2) * t11 + pkin(3) * t24 + t45, pkin(2) * t3 + pkin(3) * t7 + t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t178, t175, -t174, 0, 0, 0, 0, 0, 0, t125, t128, t141, -t134, 0, 0, 0, 0, 0, 0, t74, t79, t107, -t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, t163, t130, -t164, t250, t235, -t67, -t68, 0, 0, t120, t119, t78, -t120, t251, t227, t261, t254, t45, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t119, t78, -t120, t251, t227, -t30, -t31, 0, 0;];
tauJ_reg = t29;

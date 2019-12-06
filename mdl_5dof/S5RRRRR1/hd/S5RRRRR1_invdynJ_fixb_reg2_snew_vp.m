% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:51:02
% EndTime: 2019-12-05 18:51:17
% DurationCPUTime: 4.62s
% Computational Cost: add. (16648->352), mult. (36934->511), div. (0->0), fcn. (28909->10), ass. (0->235)
t221 = cos(qJ(3));
t250 = qJDD(2) + qJDD(3);
t222 = cos(qJ(2));
t216 = sin(qJ(3));
t217 = sin(qJ(2));
t262 = t216 * t217;
t186 = (-t222 * t221 + t262) * qJD(1);
t187 = (-t222 * t216 - t217 * t221) * qJD(1);
t271 = t186 * t187;
t285 = t250 + t271;
t292 = t221 * t285;
t215 = sin(qJ(4));
t220 = cos(qJ(4));
t218 = sin(qJ(1));
t223 = cos(qJ(1));
t233 = t223 * g(1) + t218 * g(2);
t225 = qJD(1) ^ 2;
t277 = t225 * pkin(1);
t194 = -t233 - t277;
t179 = -t222 * g(3) + t217 * t194;
t202 = t222 * t225 * t217;
t252 = qJDD(2) + t202;
t174 = t252 * pkin(2) - t179;
t180 = t217 * g(3) + t222 * t194;
t224 = qJD(2) ^ 2;
t213 = t222 ^ 2;
t267 = t213 * t225;
t176 = (-t224 - t267) * pkin(2) + t180;
t144 = t216 * t174 + t221 * t176;
t184 = t186 ^ 2;
t211 = qJD(2) + qJD(3);
t210 = t211 ^ 2;
t280 = -t184 - t210;
t130 = t280 * pkin(3) + t144;
t143 = -t221 * t174 + t216 * t176;
t227 = pkin(3) * t285 - t143;
t91 = t215 * t130 - t220 * t227;
t92 = t220 * t130 + t215 * t227;
t48 = t215 * t92 - t220 * t91;
t291 = t221 * t48;
t214 = sin(qJ(5));
t168 = t215 * t186 + t220 * t187;
t255 = qJD(1) * qJD(2);
t243 = t222 * t255;
t254 = t217 * qJDD(1);
t195 = -t243 - t254;
t244 = t217 * t255;
t253 = t222 * qJDD(1);
t198 = t244 - t253;
t235 = t216 * t195 - t221 * t198;
t155 = -t187 * qJD(3) - t235;
t156 = t186 * qJD(3) + t221 * t195 + t216 * t198;
t236 = -t220 * t155 + t215 * t156;
t110 = -t168 * qJD(4) - t236;
t109 = qJDD(5) - t110;
t207 = qJD(4) + t211;
t219 = cos(qJ(5));
t150 = t214 * t168 - t219 * t207;
t152 = t219 * t168 + t214 * t207;
t124 = t152 * t150;
t284 = t109 - t124;
t290 = t214 * t284;
t166 = -t220 * t186 + t215 * t187;
t137 = t168 * t166;
t206 = qJDD(4) + t250;
t283 = -t137 + t206;
t289 = t215 * t283;
t288 = t219 * t284;
t287 = t220 * t283;
t286 = -t198 - t244;
t242 = t218 * g(1) - t223 * g(2);
t193 = qJDD(1) * pkin(1) + t242;
t173 = t286 * pkin(2) + t193;
t268 = t211 * t187;
t122 = t173 + (-t155 + t268) * pkin(3);
t111 = -t166 * qJD(4) + t215 * t155 + t220 * t156;
t160 = t207 * t166;
t282 = -t160 + t111;
t42 = -t282 * pkin(6) + (t207 * t168 - t110) * pkin(4) + t122;
t135 = t166 * pkin(4) - t168 * pkin(6);
t279 = t207 ^ 2;
t70 = -t279 * pkin(4) + t206 * pkin(6) - t166 * t135 + t92;
t19 = t214 * t70 - t219 * t42;
t20 = t214 * t42 + t219 * t70;
t14 = t214 * t19 + t219 * t20;
t181 = t211 * t186;
t281 = t181 + t156;
t163 = qJD(5) + t166;
t237 = t214 * t111 - t219 * t206;
t77 = (qJD(5) - t163) * t152 + t237;
t148 = t150 ^ 2;
t149 = t152 ^ 2;
t162 = t163 ^ 2;
t164 = t166 ^ 2;
t165 = t168 ^ 2;
t185 = t187 ^ 2;
t278 = pkin(4) * t215;
t69 = -t206 * pkin(4) - t279 * pkin(6) + t168 * t135 + t91;
t276 = -pkin(4) * t69 + pkin(6) * t14;
t66 = t214 * t69;
t86 = t109 + t124;
t275 = t214 * t86;
t67 = t219 * t69;
t274 = t219 * t86;
t273 = t163 * t214;
t272 = t163 * t219;
t270 = t207 * t215;
t269 = t207 * t220;
t266 = t215 * t122;
t133 = t137 + t206;
t265 = t215 * t133;
t170 = -t271 + t250;
t264 = t216 * t170;
t263 = t216 * t186;
t261 = t220 * t122;
t260 = t220 * t133;
t259 = t221 * t173;
t258 = t221 * t187;
t256 = qJD(5) + t163;
t9 = t215 * t14 - t220 * t69;
t251 = pkin(3) * t9 + t276;
t119 = -t149 - t162;
t55 = -t214 * t119 - t274;
t231 = -t219 * t111 - t214 * t206;
t82 = t256 * t150 + t231;
t249 = pkin(4) * t82 + pkin(6) * t55 + t66;
t112 = -t162 - t148;
t52 = t219 * t112 - t290;
t78 = -t256 * t152 - t237;
t248 = pkin(4) * t78 + pkin(6) * t52 - t67;
t247 = t215 * t124;
t246 = t220 * t124;
t245 = -pkin(4) * t220 - pkin(3);
t241 = t215 * t91 + t220 * t92;
t30 = t215 * t55 + t220 * t82;
t240 = pkin(3) * t30 + t249;
t28 = t215 * t52 + t220 * t78;
t239 = pkin(3) * t28 + t248;
t108 = t148 + t149;
t129 = t163 * t150;
t89 = -t150 * qJD(5) - t231;
t81 = t129 + t89;
t39 = t214 * t81 - t219 * t77;
t238 = pkin(4) * t108 + pkin(6) * t39 + t14;
t22 = t220 * t108 + t215 * t39;
t234 = pkin(3) * t22 + t238;
t131 = -t279 - t164;
t106 = t215 * t131 + t287;
t232 = pkin(3) * t106 - t91;
t13 = -t219 * t19 + t214 * t20;
t230 = t221 * t143 - t216 * t144;
t228 = (-qJD(4) + t207) * t168 - t236;
t139 = -(qJD(3) - t211) * t187 - t235;
t153 = -t165 - t279;
t114 = t220 * t153 - t265;
t226 = pkin(3) * t114 - t92;
t212 = t217 ^ 2;
t197 = -0.2e1 * t244 + t253;
t196 = -0.2e1 * t243 - t254;
t178 = -t185 + t210;
t177 = t184 - t210;
t172 = t185 - t184;
t159 = -t165 + t279;
t158 = t164 - t279;
t157 = -t184 - t185;
t142 = -t181 + t156;
t138 = -(-qJD(3) - t211) * t187 + t235;
t136 = t165 - t164;
t128 = -t149 + t162;
t127 = t148 - t162;
t126 = (-t166 * t220 + t168 * t215) * t207;
t125 = (-t166 * t215 - t168 * t220) * t207;
t123 = t149 - t148;
t120 = -t164 - t165;
t118 = t220 * t158 - t265;
t117 = -t215 * t159 + t287;
t116 = t215 * t158 + t260;
t115 = t220 * t159 + t289;
t104 = (-t150 * t219 + t152 * t214) * t163;
t103 = (-t150 * t214 - t152 * t219) * t163;
t101 = t160 + t111;
t97 = (qJD(4) + t207) * t168 + t236;
t96 = t220 * t111 - t168 * t270;
t95 = t215 * t111 + t168 * t269;
t94 = -t215 * t110 + t166 * t269;
t93 = t220 * t110 + t166 * t270;
t88 = -t152 * qJD(5) - t237;
t84 = -pkin(3) * t97 - t261;
t83 = -pkin(3) * t282 + t266;
t80 = -t129 + t89;
t74 = -t152 * t273 + t219 * t89;
t73 = t152 * t272 + t214 * t89;
t72 = t150 * t272 - t214 * t88;
t71 = t150 * t273 + t219 * t88;
t65 = t220 * t104 + t215 * t109;
t64 = t215 * t104 - t220 * t109;
t63 = t219 * t127 - t275;
t62 = -t214 * t128 + t288;
t61 = t214 * t127 + t274;
t60 = t219 * t128 + t290;
t59 = -t215 * t282 - t220 * t97;
t58 = -t220 * t101 + t215 * t228;
t57 = -t215 * t97 + t220 * t282;
t56 = pkin(3) * t58;
t54 = t219 * t119 - t275;
t51 = t214 * t112 + t288;
t47 = pkin(3) * t48;
t46 = t220 * t74 + t247;
t45 = t220 * t72 - t247;
t44 = t215 * t74 - t246;
t43 = t215 * t72 + t246;
t40 = -pkin(3) * t120 + t241;
t38 = -t214 * t80 + t219 * t78;
t37 = -t214 * t77 - t219 * t81;
t36 = t214 * t78 + t219 * t80;
t34 = -t215 * t77 + t220 * t63;
t33 = t215 * t81 + t220 * t62;
t32 = t215 * t63 + t220 * t77;
t31 = t215 * t62 - t220 * t81;
t26 = t215 * t123 + t220 * t38;
t25 = -t220 * t123 + t215 * t38;
t24 = -pkin(6) * t54 + t67;
t23 = -pkin(6) * t51 + t66;
t16 = -pkin(4) * t54 + t20;
t15 = -pkin(4) * t51 + t19;
t11 = -t215 * t16 + t220 * t24;
t10 = -t215 * t15 + t220 * t23;
t7 = -pkin(6) * t37 - t13;
t6 = -pkin(3) * t54 + t220 * t16 + t215 * t24;
t5 = -pkin(3) * t51 + t220 * t15 + t215 * t23;
t4 = t220 * t7 + t37 * t278;
t3 = (-pkin(6) * t220 + t278) * t13;
t2 = t215 * t7 + t245 * t37;
t1 = (-pkin(6) * t215 + t245) * t13;
t8 = [0, 0, 0, 0, 0, qJDD(1), t242, t233, 0, 0, (-t195 + t243) * t217, -t222 * t196 + t217 * t197, -t217 * t252 - t222 * (-t212 * t225 + t224), t286 * t222, -t217 * (-t224 + t267) - t222 * (qJDD(2) - t202), 0, pkin(1) * t197 + t222 * t193, pkin(1) * t196 - t217 * t193, -t217 * t179 - t222 * t180 + (-t212 - t213) * t277, pkin(1) * t193, -t217 * (t221 * t156 - t216 * t268) - t222 * (t216 * t156 + t211 * t258), -t217 * (-t221 * t138 - t216 * t281) - t222 * (-t216 * t138 + t221 * t281), -t217 * (-t216 * t178 + t292) - t222 * (t221 * t178 + t216 * t285), -t217 * (-t216 * t155 - t181 * t221) - t222 * (t221 * t155 - t211 * t263), -t217 * (t221 * t177 - t264) - t222 * (t221 * t170 + t216 * t177), (-t217 * (t186 * t221 + t187 * t216) - t222 * (-t258 + t263)) * t211, -t173 * t262 - t222 * (-pkin(2) * t138 - t259) + pkin(1) * t138, -t217 * t259 - t222 * (-pkin(2) * t281 + t216 * t173) + pkin(1) * t281, -t217 * t230 - t222 * (-pkin(2) * t157 + t216 * t143 + t221 * t144) + pkin(1) * t157, (t222 * pkin(2) + pkin(1)) * t173, -t217 * (-t216 * t95 + t221 * t96) - t222 * (t216 * t96 + t221 * t95), -t217 * (-t216 * t57 + t221 * t59) - t222 * (t216 * t59 + t221 * t57), -t217 * (-t216 * t115 + t221 * t117) - t222 * (t221 * t115 + t216 * t117), -t217 * (-t216 * t93 + t221 * t94) - t222 * (t216 * t94 + t221 * t93), -t217 * (-t216 * t116 + t221 * t118) - t222 * (t221 * t116 + t216 * t118), -t217 * (-t216 * t125 + t221 * t126) - t222 * (t221 * t125 + t216 * t126), -t217 * (-t216 * t84 + t221 * t266) - t222 * (-pkin(2) * t97 + t216 * t266 + t221 * t84) + pkin(1) * t97, -t217 * (-t216 * t83 + t221 * t261) - t222 * (-pkin(2) * t282 + t216 * t261 + t221 * t83) + pkin(1) * t282, -t217 * (-t216 * t40 - t291) - t222 * (-pkin(2) * t120 - t216 * t48 + t221 * t40) + pkin(1) * t120, (-pkin(3) * t262 - t222 * (-pkin(3) * t221 - pkin(2)) + pkin(1)) * t122, -t217 * (-t216 * t44 + t221 * t46) - t222 * (t216 * t46 + t221 * t44), -t217 * (-t216 * t25 + t221 * t26) - t222 * (t216 * t26 + t221 * t25), -t217 * (-t216 * t31 + t221 * t33) - t222 * (t216 * t33 + t221 * t31), -t217 * (-t216 * t43 + t221 * t45) - t222 * (t216 * t45 + t221 * t43), -t217 * (-t216 * t32 + t221 * t34) - t222 * (t216 * t34 + t221 * t32), -t217 * (-t216 * t64 + t221 * t65) - t222 * (t216 * t65 + t221 * t64), -t217 * (t221 * t10 - t216 * t5) - t222 * (-pkin(2) * t51 + t216 * t10 + t221 * t5) + pkin(1) * t51, -t217 * (t221 * t11 - t216 * t6) - t222 * (-pkin(2) * t54 + t216 * t11 + t221 * t6) + pkin(1) * t54, -t217 * (-t216 * t2 + t221 * t4) - t222 * (-pkin(2) * t37 + t221 * t2 + t216 * t4) + pkin(1) * t37, -t217 * (-t216 * t1 + t221 * t3) - t222 * (-pkin(2) * t13 + t221 * t1 + t216 * t3) + pkin(1) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t202, (t212 - t213) * t225, -t254, t202, -t253, qJDD(2), -t179, -t180, 0, 0, -t271, t172, t142, t271, t139, t250, pkin(2) * (t216 * t280 + t292) - t143, pkin(2) * (-t264 + t221 * (-t185 - t210)) - t144, pkin(2) * (t139 * t216 - t221 * t142), -pkin(2) * t230, t137, t136, t101, -t137, t228, t206, pkin(2) * (t216 * (t220 * t131 - t289) + t221 * t106) + t232, pkin(2) * (t216 * (-t215 * t153 - t260) + t221 * t114) + t226, pkin(2) * (t216 * (t215 * t101 + t220 * t228) + t221 * t58) + t56, pkin(2) * (t216 * t241 + t291) + t47, t73, t36, t60, t71, t61, t103, pkin(2) * (t216 * (-t215 * t78 + t220 * t52) + t221 * t28) + t239, pkin(2) * (t216 * (-t215 * t82 + t220 * t55) + t221 * t30) + t240, pkin(2) * (t216 * (-t215 * t108 + t220 * t39) + t221 * t22) + t234, pkin(2) * (t216 * (t220 * t14 + t215 * t69) + t221 * t9) + t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t271, t172, t142, t271, t139, t250, -t143, -t144, 0, 0, t137, t136, t101, -t137, t228, t206, t232, t226, t56, t47, t73, t36, t60, t71, t61, t103, t239, t240, t234, t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, t136, t101, -t137, t228, t206, -t91, -t92, 0, 0, t73, t36, t60, t71, t61, t103, t248, t249, t238, t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t123, t81, -t124, -t77, t109, -t19, -t20, 0, 0;];
tauJ_reg = t8;

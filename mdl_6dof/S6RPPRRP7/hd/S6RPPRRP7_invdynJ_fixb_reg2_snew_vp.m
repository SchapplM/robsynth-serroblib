% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRP7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:04:27
% EndTime: 2019-05-05 15:04:38
% DurationCPUTime: 3.62s
% Computational Cost: add. (12037->334), mult. (27594->436), div. (0->0), fcn. (19633->8), ass. (0->226)
t194 = sin(qJ(4));
t190 = sin(pkin(9));
t191 = cos(pkin(9));
t197 = cos(qJ(4));
t215 = t190 * t197 + t191 * t194;
t172 = t215 * qJD(1);
t248 = t190 * t194;
t174 = (t191 * t197 - t248) * qJD(1);
t249 = t174 * t172;
t288 = qJDD(4) - t249;
t290 = t194 * t288;
t289 = t197 * t288;
t193 = sin(qJ(5));
t196 = cos(qJ(5));
t157 = -qJD(4) * t196 + t174 * t193;
t159 = qJD(4) * t193 + t174 * t196;
t128 = t159 * t157;
t233 = t174 * qJD(4);
t277 = t215 * qJDD(1);
t147 = -t277 - t233;
t138 = qJDD(5) - t147;
t278 = -t128 + t138;
t287 = pkin(5) * t278;
t261 = pkin(1) + qJ(3);
t229 = t191 * qJDD(1);
t230 = t190 * qJDD(1);
t171 = -t194 * t230 + t197 * t229;
t234 = t172 * qJD(4);
t149 = t171 - t234;
t118 = -qJD(5) * t157 + qJDD(4) * t193 + t149 * t196;
t167 = qJD(5) + t172;
t134 = t167 * t157;
t97 = t118 + t134;
t286 = qJ(6) * t97;
t246 = t193 * t278;
t242 = t196 * t278;
t199 = qJD(1) ^ 2;
t195 = sin(qJ(1));
t198 = cos(qJ(1));
t216 = t195 * g(1) - t198 * g(2);
t212 = qJDD(2) - t216;
t207 = -t199 * qJ(2) + t212;
t228 = -0.2e1 * qJD(3) * qJD(1);
t285 = -qJDD(1) * t261 + t207 + t228;
t156 = t159 ^ 2;
t166 = t167 ^ 2;
t122 = -t156 - t166;
t155 = t157 ^ 2;
t222 = -qJDD(4) * t196 + t193 * t149;
t117 = -qJD(5) * t159 - t222;
t130 = pkin(5) * t167 - qJ(6) * t159;
t141 = pkin(4) * t172 - pkin(8) * t174;
t281 = -pkin(7) - t261;
t203 = (t228 + (-pkin(3) * t190 - qJ(2)) * t199 + t281 * qJDD(1) + t212) * t191;
t262 = t190 * g(3);
t202 = t203 + t262;
t140 = -t191 * g(3) + t190 * t285;
t186 = t190 ^ 2;
t129 = -pkin(3) * t186 * t199 - pkin(7) * t230 + t140;
t241 = t197 * t129;
t274 = qJD(4) ^ 2;
t68 = -pkin(4) * t274 + qJDD(4) * pkin(8) - t172 * t141 + t194 * t202 + t241;
t231 = qJD(2) * qJD(1);
t185 = 0.2e1 * t231;
t188 = qJDD(1) * qJ(2);
t217 = t198 * g(1) + t195 * g(2);
t213 = -t188 + t217;
t210 = -qJDD(3) + t213;
t187 = t191 ^ 2;
t237 = t186 + t187;
t283 = pkin(3) * t230 - (pkin(7) * t237 + t261) * t199;
t81 = t185 + (-t149 + t234) * pkin(8) + (-t147 + t233) * pkin(4) - t210 + t283;
t46 = t193 * t81 + t196 * t68;
t211 = qJ(6) * t117 - 0.2e1 * qJD(6) * t157 - t167 * t130 + t46;
t284 = -t211 + (t122 + t155) * pkin(5);
t280 = t237 * t199;
t279 = t118 - t134;
t94 = (qJD(5) - t167) * t159 + t222;
t169 = t172 ^ 2;
t170 = t174 ^ 2;
t119 = -t166 - t155;
t72 = t119 * t193 + t242;
t273 = pkin(4) * t72;
t105 = t128 + t138;
t247 = t193 * t105;
t77 = t122 * t196 - t247;
t272 = pkin(4) * t77;
t235 = qJD(6) * t159;
t152 = -0.2e1 * t235;
t45 = t193 * t68 - t196 * t81;
t208 = -t286 - t45 + t287;
t29 = t152 + t208;
t271 = pkin(5) * t29;
t270 = pkin(5) * t97;
t113 = -t155 - t156;
t62 = t193 * t97 - t196 * t94;
t48 = -t113 * t197 + t194 * t62;
t269 = pkin(7) * t48;
t73 = t119 * t196 - t246;
t93 = (qJD(5) + t167) * t159 + t222;
t51 = t194 * t73 - t197 * t93;
t268 = pkin(7) * t51;
t243 = t196 * t105;
t78 = -t122 * t193 - t243;
t55 = t194 * t78 - t197 * t279;
t267 = pkin(7) * t55;
t60 = -t193 * t94 - t196 * t97;
t266 = pkin(8) * t60;
t265 = pkin(8) * t72;
t264 = pkin(8) * t77;
t263 = pkin(4) * t194;
t260 = -pkin(4) * t93 + pkin(8) * t73;
t259 = -pkin(4) * t279 + pkin(8) * t78;
t100 = t129 * t194 - t197 * t202;
t101 = g(3) * t248 + t194 * t203 + t241;
t63 = -t100 * t197 + t101 * t194;
t258 = t191 * t63;
t257 = t193 * t29;
t67 = -qJDD(4) * pkin(4) - pkin(8) * t274 + t141 * t174 + t100;
t256 = t193 * t67;
t255 = t196 * t29;
t254 = t196 * t67;
t253 = -pkin(4) * t113 + pkin(8) * t62;
t252 = qJDD(1) * pkin(1);
t251 = t167 * t193;
t250 = t167 * t196;
t205 = t210 - 0.2e1 * t231;
t136 = t205 - t283;
t245 = t194 * t136;
t144 = qJDD(4) + t249;
t244 = t194 * t144;
t240 = t197 * t136;
t239 = t197 * t144;
t227 = t194 * t128;
t226 = t197 * t128;
t225 = -pkin(4) * t197 - pkin(3);
t52 = t194 * t93 + t197 * t73;
t224 = -pkin(3) * t72 + pkin(7) * t52;
t56 = t194 * t279 + t197 * t78;
t223 = -pkin(3) * t77 + pkin(7) * t56;
t21 = t193 * t45 + t196 * t46;
t64 = t100 * t194 + t101 * t197;
t161 = t199 * t261 + t205;
t221 = -t161 + t188;
t49 = t113 * t194 + t197 * t62;
t23 = t190 * t49 + t191 * t48;
t220 = qJ(2) * t60 - t23 * t261;
t26 = t190 * t52 + t191 * t51;
t219 = qJ(2) * t72 - t26 * t261;
t28 = t190 * t56 + t191 * t55;
t218 = qJ(2) * t77 - t261 * t28;
t3 = t191 * (t194 * t21 - t197 * t67) + t190 * (t194 * t67 + t197 * t21);
t20 = t193 * t46 - t196 * t45;
t107 = t191 * (t191 * t285 + t262) + t190 * t140;
t206 = t208 + t287;
t41 = -pkin(5) * t117 - qJ(6) * t155 + t130 * t159 + qJDD(6) + t67;
t177 = t237 * qJDD(1);
t176 = t190 * t280;
t175 = t191 * t280;
t168 = -t207 + t252;
t164 = -t170 - t274;
t163 = -t170 + t274;
t162 = t169 - t274;
t153 = 0.2e1 * t235;
t148 = t171 - 0.2e1 * t234;
t146 = t277 + 0.2e1 * t233;
t142 = -t274 - t169;
t132 = -t156 + t166;
t131 = t155 - t166;
t127 = -t169 - t170;
t125 = t156 - t155;
t121 = -t164 * t194 - t239;
t120 = t164 * t197 - t244;
t112 = t171 * t194 - t197 * t277;
t111 = -t171 * t197 - t194 * t277;
t109 = t142 * t197 - t290;
t108 = t142 * t194 + t289;
t103 = (-t157 * t196 + t159 * t193) * t167;
t102 = (-t157 * t193 - t159 * t196) * t167;
t90 = t118 * t196 - t159 * t251;
t89 = t118 * t193 + t159 * t250;
t88 = -t117 * t193 + t157 * t250;
t87 = t117 * t196 + t157 * t251;
t86 = t120 * t191 + t121 * t190;
t85 = t131 * t196 - t247;
t84 = -t132 * t193 + t242;
t83 = t131 * t193 + t243;
t82 = t132 * t196 + t246;
t74 = t111 * t191 + t112 * t190;
t71 = t108 * t191 + t109 * t190;
t65 = -pkin(5) * t279 - qJ(6) * t105;
t61 = -t193 * t279 - t196 * t93;
t59 = -t193 * t93 + t196 * t279;
t53 = t191 * (t103 * t197 + t138 * t194) - t190 * (t103 * t194 - t138 * t197);
t47 = pkin(7) * t49;
t44 = t254 - t264;
t42 = t256 - t265;
t40 = -pkin(4) * t60 + t270;
t39 = t191 * (t197 * t90 + t227) - t190 * (t194 * t90 - t226);
t38 = t191 * (t197 * t88 - t227) - t190 * (t194 * t88 + t226);
t37 = -qJ(6) * t122 + t41;
t36 = t190 * t64 + t258;
t35 = t46 - t272;
t34 = t45 - t273;
t33 = -pkin(5) * t155 + t211;
t32 = -pkin(5) * t93 + qJ(6) * t119 - t41;
t31 = t191 * (-t194 * t94 + t197 * t85) - t190 * (t194 * t85 + t197 * t94);
t30 = t191 * (t194 * t97 + t197 * t84) - t190 * (t194 * t84 - t197 * t97);
t24 = t191 * (t125 * t194 + t197 * t61) - t190 * (-t125 * t197 + t194 * t61);
t19 = t153 - t208 + t286;
t18 = -qJ(6) * t94 + (-t113 - t155) * pkin(5) + t211;
t17 = -t272 - t284;
t16 = -t193 * t65 + t196 * t37 - t264;
t15 = -qJ(6) * t242 - t193 * t32 - t265;
t14 = t153 - t206 - t273;
t11 = -pkin(5) * t41 + qJ(6) * t33;
t10 = -t20 - t266;
t9 = t196 * t33 - t257;
t8 = t193 * t33 + t255;
t7 = t194 * t41 + t197 * t9;
t6 = t194 * t9 - t197 * t41;
t5 = -t18 * t193 + t19 * t196 - t266;
t4 = -pkin(4) * t8 - t271;
t2 = -pkin(8) * t8 - qJ(6) * t255 - t11 * t193;
t1 = t190 * t7 + t191 * t6;
t12 = [0, 0, 0, 0, 0, qJDD(1), t216, t217, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t212 - 0.2e1 * t252, t185 + 0.2e1 * t188 - t217, pkin(1) * t168 + qJ(2) * (-t199 * pkin(1) + t185 - t213), t187 * qJDD(1), -0.2e1 * t190 * t229, 0, t186 * qJDD(1), 0, 0, t176 * t261 + t190 * t221, t175 * t261 + t191 * t221, -qJ(2) * t280 + t177 * t261 - t107, -qJ(2) * t161 - t107 * t261, t191 * (t149 * t197 - t194 * t233) - t190 * (t149 * t194 + t197 * t233), t191 * (-t146 * t197 - t148 * t194) - t190 * (-t146 * t194 + t148 * t197), t191 * (-t163 * t194 + t289) - t190 * (t163 * t197 + t290), t191 * (-t147 * t194 + t197 * t234) - t190 * (t147 * t197 + t194 * t234), t191 * (t162 * t197 - t244) - t190 * (t162 * t194 + t239), (t191 * (-t172 * t197 + t174 * t194) - t190 * (-t172 * t194 - t174 * t197)) * qJD(4), t191 * (-pkin(7) * t108 - t245) - t190 * (-pkin(3) * t146 + pkin(7) * t109 + t240) + qJ(2) * t146 - t261 * t71, t191 * (-pkin(7) * t120 - t240) - t190 * (-pkin(3) * t148 + pkin(7) * t121 - t245) + qJ(2) * t148 - t261 * t86, t191 * (-pkin(7) * t111 - t63) - t190 * (-pkin(3) * t127 + pkin(7) * t112 + t64) + qJ(2) * t127 - t261 * t74, -pkin(7) * t258 - t190 * (pkin(3) * t136 + pkin(7) * t64) - qJ(2) * t136 - t261 * t36, t39, t24, t30, t38, t31, t53, t191 * (-t194 * t34 + t197 * t42 - t268) - t190 * (t194 * t42 + t197 * t34 + t224) + t219, t191 * (-t194 * t35 + t197 * t44 - t267) - t190 * (t194 * t44 + t197 * t35 + t223) + t218, t191 * (t10 * t197 + t263 * t60 - t269) - t190 * (t194 * t10 + t225 * t60 + t47) + t220, (t191 * (-pkin(8) * t197 + t263) - t190 * (-pkin(8) * t194 + t225) + qJ(2)) * t20 + t281 * t3, t39, t24, t30, t38, t31, t53, t191 * (-t14 * t194 + t15 * t197 - t268) - t190 * (t14 * t197 + t15 * t194 + t224) + t219, t191 * (t16 * t197 - t17 * t194 - t267) - t190 * (t16 * t194 + t17 * t197 + t223) + t218, t191 * (-t194 * t40 + t197 * t5 - t269) - t190 * (-pkin(3) * t60 + t194 * t5 + t197 * t40 + t47) + t220, t191 * (-pkin(7) * t6 - t194 * t4 + t197 * t2) - t190 * (-pkin(3) * t8 + pkin(7) * t7 + t194 * t2 + t197 * t4) + qJ(2) * t8 - t261 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t199, -t168, 0, 0, 0, 0, 0, 0, -t176, -t175, -t177, t107, 0, 0, 0, 0, 0, 0, t71, t86, t74, t36, 0, 0, 0, 0, 0, 0, t26, t28, t23, t3, 0, 0, 0, 0, 0, 0, t26, t28, t23, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, t229, -t280, -t161, 0, 0, 0, 0, 0, 0, t146, t148, t127, -t136, 0, 0, 0, 0, 0, 0, t72, t77, t60, t20, 0, 0, 0, 0, 0, 0, t72, t77, t60, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, t170 - t169, t171, -t249, -t277, qJDD(4), -t100, -t101, 0, 0, t89, t59, t82, t87, t83, t102, -t254 + t260, t256 + t259, t21 + t253, -pkin(4) * t67 + pkin(8) * t21, t89, t59, t82, t87, t83, t102, -qJ(6) * t246 + t196 * t32 + t260, t193 * t37 + t196 * t65 + t259, t18 * t196 + t19 * t193 + t253, -pkin(4) * t41 + pkin(8) * t9 - qJ(6) * t257 + t11 * t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t125, t97, -t128, -t94, t138, -t45, -t46, 0, 0, t128, t125, t97, -t128, -t94, t138, t152 + t206, t284, -t270, t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t279, t113, t41;];
tauJ_reg  = t12;

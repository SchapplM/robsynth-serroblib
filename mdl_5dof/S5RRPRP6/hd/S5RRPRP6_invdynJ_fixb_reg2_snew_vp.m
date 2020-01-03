% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:29
% EndTime: 2019-12-31 19:58:37
% DurationCPUTime: 2.83s
% Computational Cost: add. (10575->348), mult. (24399->458), div. (0->0), fcn. (16843->8), ass. (0->221)
t193 = sin(pkin(8));
t194 = cos(pkin(8));
t200 = cos(qJ(2));
t243 = qJD(1) * t200;
t197 = sin(qJ(2));
t244 = qJD(1) * t197;
t169 = t193 * t243 + t194 * t244;
t196 = sin(qJ(4));
t199 = cos(qJ(4));
t152 = -qJD(2) * t199 + t169 * t196;
t154 = qJD(2) * t196 + t169 * t199;
t124 = t154 * t152;
t183 = t197 * qJDD(1);
t234 = qJD(1) * qJD(2);
t224 = t200 * t234;
t174 = t183 + t224;
t185 = t200 * qJDD(1);
t225 = t197 * t234;
t175 = t185 - t225;
t218 = t174 * t193 - t175 * t194;
t141 = qJDD(4) + t218;
t282 = -t124 + t141;
t289 = pkin(4) * t282;
t144 = t174 * t194 + t175 * t193;
t112 = -qJD(4) * t152 + qJDD(2) * t196 + t144 * t199;
t167 = t193 * t244 - t194 * t243;
t164 = qJD(4) + t167;
t133 = t164 * t152;
t94 = t112 + t133;
t288 = qJ(5) * t94;
t142 = t169 * t167;
t281 = qJDD(2) - t142;
t287 = t193 * t281;
t286 = t194 * t281;
t257 = t282 * t196;
t256 = t282 * t199;
t151 = t154 ^ 2;
t163 = t164 ^ 2;
t115 = -t151 - t163;
t150 = t152 ^ 2;
t220 = -qJDD(2) * t199 + t196 * t144;
t111 = -qJD(4) * t154 - t220;
t129 = pkin(4) * t164 - qJ(5) * t154;
t135 = pkin(3) * t167 - pkin(7) * t169;
t201 = qJD(2) ^ 2;
t202 = qJD(1) ^ 2;
t198 = sin(qJ(1));
t273 = cos(qJ(1));
t213 = g(1) * t273 + g(2) * t198;
t260 = qJDD(1) * pkin(6);
t206 = -pkin(1) * t202 - t213 + t260;
t156 = -t197 * g(3) + t200 * t206;
t191 = t200 ^ 2;
t188 = t191 * t202;
t211 = qJD(2) * pkin(2) - qJ(3) * t244;
t119 = -pkin(2) * t188 + t175 * qJ(3) - qJD(2) * t211 + t156;
t204 = t197 * t206;
t248 = t197 * t202;
t203 = -t204 - t174 * qJ(3) + qJDD(2) * pkin(2) + (pkin(2) * t248 + qJ(3) * t234 - g(3)) * t200;
t79 = -0.2e1 * qJD(3) * t167 + t119 * t194 + t193 * t203;
t65 = -pkin(3) * t201 + qJDD(2) * pkin(7) - t135 * t167 + t79;
t222 = g(1) * t198 - g(2) * t273;
t212 = qJDD(1) * pkin(1) + t222;
t122 = pkin(2) * t175 + (qJ(3) * t191 + pkin(6)) * t202 - t211 * t244 - qJDD(3) + t212;
t242 = qJD(2) * t169;
t125 = t218 + t242;
t162 = qJD(2) * t167;
t127 = t144 - t162;
t72 = pkin(3) * t125 - t127 * pkin(7) - t122;
t36 = t196 * t72 + t199 * t65;
t210 = qJ(5) * t111 - 0.2e1 * qJD(5) * t152 - t164 * t129 + t36;
t285 = -t210 + (t115 + t150) * pkin(4);
t283 = t112 - t133;
t91 = (qJD(4) - t164) * t154 + t220;
t165 = t167 ^ 2;
t166 = t169 ^ 2;
t107 = -t163 - t150;
t68 = t107 * t196 + t256;
t280 = pkin(3) * t68;
t101 = t124 + t141;
t259 = t101 * t196;
t75 = t115 * t199 - t259;
t279 = pkin(3) * t75;
t236 = qJD(5) * t154;
t146 = -0.2e1 * t236;
t35 = t196 * t65 - t199 * t72;
t209 = -t288 - t35 + t289;
t21 = t146 + t209;
t278 = pkin(4) * t21;
t277 = pkin(4) * t94;
t59 = -t196 * t91 - t199 * t94;
t276 = pkin(7) * t59;
t275 = pkin(7) * t68;
t274 = pkin(7) * t75;
t272 = pkin(3) * t193;
t106 = -t150 - t151;
t61 = t196 * t94 - t199 * t91;
t43 = -t106 * t194 + t193 * t61;
t44 = t106 * t193 + t194 * t61;
t271 = pkin(6) * (-t197 * t43 + t200 * t44) - pkin(1) * t59;
t69 = t107 * t199 - t257;
t90 = (qJD(4) + t164) * t154 + t220;
t49 = t193 * t69 - t194 * t90;
t50 = t193 * t90 + t194 * t69;
t270 = pkin(6) * (-t197 * t49 + t200 * t50) - pkin(1) * t68;
t258 = t101 * t199;
t76 = -t115 * t196 - t258;
t54 = t193 * t76 - t194 * t283;
t55 = t193 * t283 + t194 * t76;
t269 = pkin(6) * (-t197 * t54 + t200 * t55) - pkin(1) * t75;
t268 = qJ(3) * t43;
t267 = qJ(3) * t49;
t266 = qJ(3) * t54;
t265 = t196 * t21;
t221 = t193 * t119 - t194 * t203;
t214 = -qJDD(2) * pkin(3) - pkin(7) * t201 + t221;
t217 = (0.2e1 * qJD(3) + t135) * t169;
t64 = t217 + t214;
t264 = t196 * t64;
t238 = qJD(3) * t169;
t78 = t221 + 0.2e1 * t238;
t45 = t193 * t79 - t194 * t78;
t263 = t197 * t45;
t262 = t199 * t21;
t261 = t199 * t64;
t255 = t122 * t193;
t254 = t122 * t194;
t138 = qJDD(2) + t142;
t253 = t138 * t193;
t252 = t138 * t194;
t251 = t164 * t196;
t250 = t164 * t199;
t182 = t200 * t248;
t249 = t197 * (qJDD(2) + t182);
t247 = t200 * (qJDD(2) - t182);
t241 = qJD(2) * t193;
t240 = qJD(2) * t194;
t233 = pkin(2) * t49 - pkin(3) * t90 + pkin(7) * t69;
t232 = pkin(2) * t54 - pkin(3) * t283 + pkin(7) * t76;
t231 = pkin(2) * t43 - pkin(3) * t106 + pkin(7) * t61;
t230 = t193 * t124;
t229 = t194 * t124;
t228 = -pkin(3) * t194 - pkin(2);
t227 = -pkin(2) * t68 + qJ(3) * t50;
t226 = -pkin(2) * t75 + qJ(3) * t55;
t46 = t193 * t78 + t194 * t79;
t18 = t196 * t35 + t199 * t36;
t155 = t200 * g(3) + t204;
t219 = t197 * t155 + t156 * t200;
t17 = t196 * t36 - t199 * t35;
t126 = -t218 + t242;
t208 = t209 + t289;
t205 = -pkin(4) * t111 - qJ(5) * t150 + t129 * t154 + qJDD(5) + t214;
t32 = t217 + t205;
t190 = t197 ^ 2;
t186 = t190 * t202;
t176 = t185 - 0.2e1 * t225;
t173 = t183 + 0.2e1 * t224;
t171 = pkin(6) * t202 + t212;
t161 = -0.2e1 * t238;
t159 = -t166 - t201;
t158 = -t166 + t201;
t157 = t165 - t201;
t147 = 0.2e1 * t236;
t136 = -t201 - t165;
t131 = -t151 + t163;
t130 = t150 - t163;
t128 = t144 + t162;
t123 = -t165 - t166;
t120 = t151 - t150;
t114 = -t159 * t193 - t252;
t113 = t159 * t194 - t253;
t104 = t136 * t194 - t287;
t103 = t136 * t193 + t286;
t99 = (-t152 * t199 + t154 * t196) * t164;
t98 = (-t152 * t196 - t154 * t199) * t164;
t97 = t126 * t194 + t128 * t193;
t96 = t126 * t193 - t128 * t194;
t87 = t112 * t199 - t154 * t251;
t86 = t112 * t196 + t154 * t250;
t85 = -t111 * t196 + t152 * t250;
t84 = t111 * t199 + t152 * t251;
t83 = t130 * t199 - t259;
t82 = -t131 * t196 + t256;
t81 = t130 * t196 + t258;
t80 = t131 * t199 + t257;
t63 = -pkin(4) * t283 - qJ(5) * t101;
t60 = -t196 * t283 - t199 * t90;
t58 = -t196 * t90 + t199 * t283;
t51 = t197 * (t141 * t193 + t194 * t99) + t200 * (-t141 * t194 + t193 * t99);
t41 = qJ(3) * t44;
t40 = -pkin(3) * t59 + t277;
t39 = t261 - t274;
t38 = t197 * (t194 * t87 + t230) + t200 * (t193 * t87 - t229);
t37 = t197 * (t194 * t85 - t230) + t200 * (t193 * t85 + t229);
t33 = t264 - t275;
t31 = -qJ(5) * t115 + t32;
t30 = t36 - t279;
t29 = t197 * (-t193 * t91 + t194 * t83) + t200 * (t193 * t83 + t194 * t91);
t28 = t197 * (t193 * t94 + t194 * t82) + t200 * (t193 * t82 - t194 * t94);
t27 = t35 - t280;
t26 = -pkin(4) * t150 + t210;
t25 = -pkin(4) * t90 + qJ(5) * t107 - t135 * t169 + t161 - t205;
t22 = t197 * (t120 * t193 + t194 * t60) + t200 * (-t120 * t194 + t193 * t60);
t19 = t147 - t209 + t288;
t16 = -qJ(5) * t91 + (-t106 - t150) * pkin(4) + t210;
t15 = -t279 - t285;
t14 = -t196 * t63 + t199 * t31 - t274;
t13 = -qJ(5) * t256 - t196 * t25 - t275;
t12 = t147 - t208 - t280;
t11 = -pkin(4) * t32 + qJ(5) * t26;
t9 = t18 * t193 - t194 * t64;
t8 = -t17 - t276;
t7 = t199 * t26 - t265;
t6 = t196 * t26 + t262;
t5 = t193 * t32 + t194 * t7;
t4 = t193 * t7 - t194 * t32;
t3 = -t16 * t196 + t19 * t199 - t276;
t2 = -pkin(3) * t6 - t278;
t1 = -pkin(7) * t6 - qJ(5) * t262 - t11 * t196;
t10 = [0, 0, 0, 0, 0, qJDD(1), t222, t213, 0, 0, (t174 + t224) * t197, t173 * t200 + t176 * t197, t249 + t200 * (-t186 + t201), (t175 - t225) * t200, t197 * (t188 - t201) + t247, 0, t200 * t171 + pkin(1) * t176 + pkin(6) * (t200 * (-t188 - t201) - t249), -t197 * t171 - pkin(1) * t173 + pkin(6) * (-t247 - t197 * (-t186 - t201)), pkin(1) * (t186 + t188) + (t190 + t191) * t260 + t219, pkin(1) * t171 + pkin(6) * t219, t197 * (t144 * t194 - t169 * t241) + t200 * (t144 * t193 + t169 * t240), t197 * (-t125 * t194 - t127 * t193) + t200 * (-t125 * t193 + t127 * t194), t197 * (-t158 * t193 + t286) + t200 * (t158 * t194 + t287), t197 * (t167 * t240 + t193 * t218) + t200 * (t167 * t241 - t194 * t218), t197 * (t157 * t194 - t253) + t200 * (t157 * t193 + t252), (t197 * (-t167 * t194 + t169 * t193) + t200 * (-t167 * t193 - t169 * t194)) * qJD(2), t197 * (-qJ(3) * t103 - t255) + t200 * (-pkin(2) * t125 + qJ(3) * t104 + t254) - pkin(1) * t125 + pkin(6) * (-t103 * t197 + t104 * t200), t197 * (-qJ(3) * t113 - t254) + t200 * (-pkin(2) * t127 + qJ(3) * t114 - t255) - pkin(1) * t127 + pkin(6) * (-t113 * t197 + t114 * t200), t197 * (-qJ(3) * t96 - t45) + t200 * (-pkin(2) * t123 + qJ(3) * t97 + t46) - pkin(1) * t123 + pkin(6) * (-t197 * t96 + t200 * t97), -qJ(3) * t263 + t200 * (pkin(2) * t122 + qJ(3) * t46) + pkin(1) * t122 + pkin(6) * (t200 * t46 - t263), t38, t22, t28, t37, t29, t51, t197 * (-t193 * t27 + t194 * t33 - t267) + t200 * (t193 * t33 + t194 * t27 + t227) + t270, t197 * (-t193 * t30 + t194 * t39 - t266) + t200 * (t193 * t39 + t194 * t30 + t226) + t269, t197 * (t194 * t8 + t272 * t59 - t268) + t200 * (t193 * t8 + t228 * t59 + t41) + t271, (t197 * (-pkin(7) * t194 + t272) + t200 * (-pkin(7) * t193 + t228) - pkin(1)) * t17 + (pkin(6) + qJ(3)) * ((t18 * t194 + t193 * t64) * t200 - t197 * t9), t38, t22, t28, t37, t29, t51, t197 * (-t12 * t193 + t13 * t194 - t267) + t200 * (t12 * t194 + t13 * t193 + t227) + t270, t197 * (t14 * t194 - t15 * t193 - t266) + t200 * (t14 * t193 + t15 * t194 + t226) + t269, t197 * (-t193 * t40 + t194 * t3 - t268) + t200 * (-pkin(2) * t59 + t193 * t3 + t194 * t40 + t41) + t271, t197 * (-qJ(3) * t4 + t1 * t194 - t193 * t2) + t200 * (-pkin(2) * t6 + qJ(3) * t5 + t1 * t193 + t194 * t2) - pkin(1) * t6 + pkin(6) * (-t197 * t4 + t200 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t186 - t188, t183, t182, t185, qJDD(2), -t155, -t156, 0, 0, t142, t166 - t165, t128, -t142, t126, qJDD(2), pkin(2) * t103 + t161 - t221, pkin(2) * t113 - t79, pkin(2) * t96, pkin(2) * t45, t86, t58, t80, t84, t81, t98, t233 - t261, t232 + t264, t18 + t231, pkin(2) * t9 - pkin(3) * t64 + pkin(7) * t18, t86, t58, t80, t84, t81, t98, -qJ(5) * t257 + t199 * t25 + t233, t196 * t31 + t199 * t63 + t232, t16 * t199 + t19 * t196 + t231, pkin(2) * t4 - pkin(3) * t32 + pkin(7) * t7 - qJ(5) * t265 + t11 * t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t127, t123, -t122, 0, 0, 0, 0, 0, 0, t68, t75, t59, t17, 0, 0, 0, 0, 0, 0, t68, t75, t59, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t120, t94, -t124, -t91, t141, -t35, -t36, 0, 0, t124, t120, t94, -t124, -t91, t141, t146 + t208, t285, -t277, t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t283, t106, t32;];
tauJ_reg = t10;

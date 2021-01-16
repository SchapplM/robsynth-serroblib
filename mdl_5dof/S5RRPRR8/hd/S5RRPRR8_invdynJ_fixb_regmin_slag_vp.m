% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR8
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
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:35:11
% EndTime: 2021-01-15 21:35:22
% DurationCPUTime: 2.98s
% Computational Cost: add. (4039->336), mult. (9756->446), div. (0->0), fcn. (7588->14), ass. (0->188)
t171 = cos(qJ(5));
t219 = qJD(5) * t171;
t164 = sin(pkin(9));
t165 = cos(pkin(9));
t169 = sin(qJ(2));
t173 = cos(qJ(2));
t122 = t164 * t173 + t165 * t169;
t113 = t122 * qJD(1);
t168 = sin(qJ(4));
t230 = t165 * t173;
t210 = qJD(1) * t230;
t222 = qJD(1) * t169;
t111 = -t164 * t222 + t210;
t172 = cos(qJ(4));
t96 = t172 * t111;
t68 = -t168 * t113 + t96;
t273 = t171 * t68;
t281 = t219 - t273;
t161 = qJ(2) + pkin(9);
t158 = qJ(4) + t161;
t150 = sin(t158);
t170 = sin(qJ(1));
t174 = cos(qJ(1));
t195 = g(1) * t174 + g(2) * t170;
t280 = t195 * t150;
t160 = qJD(2) + qJD(4);
t233 = t68 * t160;
t221 = qJD(4) * t168;
t112 = t122 * qJD(2);
t216 = t173 * qJDD(1);
t217 = t169 * qJDD(1);
t192 = t164 * t217 - t165 * t216;
t77 = qJD(1) * t112 + t192;
t218 = qJD(1) * qJD(2);
t209 = t169 * t218;
t182 = qJDD(1) * t122 - t164 * t209;
t208 = t173 * t218;
t78 = t165 * t208 + t182;
t25 = qJD(4) * t96 - t113 * t221 - t168 * t77 + t172 * t78;
t279 = t25 - t233;
t159 = qJDD(2) + qJDD(4);
t167 = sin(qJ(5));
t191 = t168 * t111 + t172 * t113;
t220 = qJD(5) * t167;
t13 = t167 * t159 + t160 * t219 + t171 * t25 - t191 * t220;
t54 = t167 * t160 + t171 * t191;
t14 = qJD(5) * t54 - t171 * t159 + t167 * t25;
t52 = -t171 * t160 + t167 * t191;
t278 = t13 * t171 - t167 * t14 - t281 * t52;
t11 = t13 * t167;
t277 = t281 * t54 + t11;
t26 = qJD(4) * t191 + t168 * t78 + t172 * t77;
t24 = qJDD(5) + t26;
t18 = t167 * t24;
t244 = t54 * t191;
t271 = qJD(5) - t68;
t57 = t271 * t219;
t276 = -t271 * t273 + t18 - t244 + t57;
t247 = t113 * pkin(7);
t166 = -qJ(3) - pkin(6);
t137 = t166 * t173;
t129 = qJD(1) * t137;
t116 = t164 * t129;
t136 = t166 * t169;
t128 = qJD(1) * t136;
t238 = qJD(2) * pkin(2);
t120 = t128 + t238;
t75 = t165 * t120 + t116;
t47 = qJD(2) * pkin(3) - t247 + t75;
t248 = t111 * pkin(7);
t231 = t165 * t129;
t76 = t164 * t120 - t231;
t51 = t76 + t248;
t28 = -t168 * t51 + t172 * t47;
t20 = -t160 * pkin(4) - t28;
t275 = t20 * t68;
t151 = cos(t158);
t250 = g(3) * t151;
t200 = qJD(2) * t166;
t109 = -t169 * qJD(3) + t173 * t200;
t74 = qJDD(2) * pkin(2) + qJD(1) * t109 + qJDD(1) * t136;
t108 = t173 * qJD(3) + t169 * t200;
t81 = qJD(1) * t108 - qJDD(1) * t137;
t40 = -t164 * t81 + t165 * t74;
t27 = qJDD(2) * pkin(3) - t78 * pkin(7) + t40;
t29 = t168 * t47 + t172 * t51;
t41 = t164 * t74 + t165 * t81;
t30 = -t77 * pkin(7) + t41;
t257 = qJD(4) * t29 + t168 * t30 - t172 * t27;
t3 = -t159 * pkin(4) + t257;
t274 = t3 + t250;
t272 = t191 * t68;
t234 = t191 * t160;
t269 = -t26 + t234;
t267 = t191 ^ 2 - t68 ^ 2;
t143 = g(3) * t150;
t256 = (qJD(4) * t47 + t30) * t172 + t168 * t27 - t51 * t221;
t153 = t173 * pkin(2) + pkin(1);
t130 = -qJD(1) * t153 + qJD(3);
t84 = -t111 * pkin(3) + t130;
t266 = t195 * t151 - t84 * t68 + t143 - t256;
t264 = pkin(4) * t191;
t243 = t191 * t52;
t152 = t165 * pkin(2) + pkin(3);
t254 = pkin(2) * t164;
t224 = t168 * t152 + t172 * t254;
t107 = pkin(8) + t224;
t87 = pkin(2) * t222 + t113 * pkin(3);
t263 = (-t68 * pkin(8) + qJD(5) * t107 + t264 + t87) * t271;
t262 = (t271 * pkin(8) + t264) * t271;
t261 = t271 * t191;
t21 = t160 * pkin(8) + t29;
t31 = -pkin(4) * t68 - pkin(8) * t191 + t84;
t6 = -t167 * t21 + t171 * t31;
t260 = t171 * t280 - t6 * t191 + t20 * t220;
t7 = t167 * t31 + t171 * t21;
t259 = t274 * t167 + t7 * t191 + t20 * t219;
t258 = -t191 * t84 - t250 - t257 + t280;
t206 = t159 * pkin(8) + qJD(5) * t31 + t256;
t121 = t164 * t169 - t230;
t79 = t172 * t121 + t168 * t122;
t80 = -t168 * t121 + t172 * t122;
t90 = t121 * pkin(3) - t153;
t35 = t79 * pkin(4) - t80 * pkin(8) + t90;
t85 = t165 * t136 + t164 * t137;
t60 = -t122 * pkin(7) + t85;
t86 = t164 * t136 - t165 * t137;
t61 = -t121 * pkin(7) + t86;
t37 = t168 * t60 + t172 * t61;
t115 = t121 * qJD(2);
t42 = -qJD(4) * t79 - t168 * t112 - t172 * t115;
t36 = t168 * t61 - t172 * t60;
t58 = -t164 * t108 + t165 * t109;
t44 = t115 * pkin(7) + t58;
t59 = t165 * t108 + t164 * t109;
t45 = -t112 * pkin(7) + t59;
t8 = -qJD(4) * t36 + t168 * t44 + t172 * t45;
t255 = t20 * t42 - (qJD(5) * t35 + t8) * t271 - t206 * t79 - t37 * t24 + t3 * t80;
t253 = pkin(2) * t169;
t249 = g(3) * t173;
t246 = t20 * t80;
t245 = t35 * t24;
t188 = t172 * t152 - t168 * t254;
t82 = -t164 * t128 + t231;
t55 = t82 - t248;
t83 = t165 * t128 + t116;
t56 = t83 - t247;
t240 = -t188 * qJD(4) + t168 * t55 + t172 * t56;
t239 = t224 * qJD(4) - t168 * t56 + t172 * t55;
t237 = t167 * t54;
t229 = t170 * t167;
t228 = t170 * t171;
t227 = t174 * t167;
t226 = t174 * t171;
t162 = t169 ^ 2;
t223 = -t173 ^ 2 + t162;
t155 = t169 * t238;
t212 = t80 * t220;
t88 = t112 * pkin(3) + t155;
t105 = pkin(2) * t209 - qJDD(1) * t153 + qJDD(3);
t50 = t77 * pkin(3) + t105;
t5 = t26 * pkin(4) - t25 * pkin(8) + t50;
t205 = qJD(5) * t21 - t5;
t197 = t167 * t271;
t194 = g(1) * t170 - g(2) * t174;
t193 = t24 * t80 + t271 * t42;
t19 = t171 * t24;
t190 = t19 - (-t167 * t68 + t220) * t271;
t189 = -t206 + t143;
t186 = -0.2e1 * pkin(1) * t218 - pkin(6) * qJDD(2);
t185 = -pkin(8) * t24 + t271 * t28 - t275;
t183 = -t107 * t24 + t240 * t271 - t275;
t175 = qJD(2) ^ 2;
t181 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t175 + t194;
t176 = qJD(1) ^ 2;
t180 = pkin(1) * t176 - pkin(6) * qJDD(1) + t195;
t157 = cos(t161);
t156 = sin(t161);
t106 = -pkin(4) - t188;
t104 = t151 * t226 + t229;
t103 = -t151 * t227 + t228;
t102 = -t151 * t228 + t227;
t101 = t151 * t229 + t226;
t43 = qJD(4) * t80 + t172 * t112 - t168 * t115;
t15 = t43 * pkin(4) - t42 * pkin(8) + t88;
t9 = qJD(4) * t37 + t168 * t45 - t172 * t44;
t4 = t171 * t5;
t1 = [qJDD(1), t194, t195, t162 * qJDD(1) + 0.2e1 * t169 * t208, 0.2e1 * t169 * t216 - 0.2e1 * t218 * t223, qJDD(2) * t169 + t175 * t173, qJDD(2) * t173 - t175 * t169, 0, t169 * t186 + t173 * t181, -t169 * t181 + t173 * t186, t85 * qJDD(2) + t105 * t121 + t130 * t112 - t153 * t77 + t194 * t157 + (-t111 * t253 + t58) * qJD(2), -t86 * qJDD(2) + t105 * t122 - t130 * t115 - t153 * t78 - t194 * t156 + (t113 * t253 - t59) * qJD(2), t59 * t111 - t76 * t112 - t58 * t113 + t75 * t115 - t41 * t121 - t40 * t122 - t86 * t77 - t85 * t78 - t195, t41 * t86 + t76 * t59 + t40 * t85 + t75 * t58 - t105 * t153 + t130 * t155 - g(1) * (-t170 * t153 - t166 * t174) - g(2) * (t174 * t153 - t170 * t166), t191 * t42 + t25 * t80, -t191 * t43 - t25 * t79 - t80 * t26 + t42 * t68, t80 * t159 + t42 * t160, -t79 * t159 - t43 * t160, 0, t151 * t194 - t36 * t159 - t9 * t160 + t90 * t26 + t84 * t43 + t50 * t79 - t68 * t88, -t150 * t194 - t37 * t159 - t8 * t160 + t191 * t88 + t90 * t25 + t84 * t42 + t50 * t80, -t54 * t212 + (t13 * t80 + t42 * t54) * t171, (-t171 * t52 - t237) * t42 + (-t11 - t14 * t171 + (t167 * t52 - t171 * t54) * qJD(5)) * t80, t13 * t79 + t171 * t193 - t212 * t271 + t54 * t43, -t14 * t79 - t167 * t193 - t52 * t43 - t57 * t80, t24 * t79 + t271 * t43, -g(1) * t102 - g(2) * t104 + t36 * t14 + t4 * t79 + t6 * t43 + t9 * t52 + (t15 * t271 + t245 + (-t21 * t79 - t271 * t37 + t246) * qJD(5)) * t171 + t255 * t167, -g(1) * t101 - g(2) * t103 + t36 * t13 - t7 * t43 + t9 * t54 + (-(-qJD(5) * t37 + t15) * t271 - t245 + t205 * t79 - qJD(5) * t246) * t167 + t255 * t171; 0, 0, 0, -t169 * t176 * t173, t223 * t176, t217, t216, qJDD(2), t169 * t180 - t249, g(3) * t169 + t173 * t180, -g(3) * t157 - t82 * qJD(2) - t130 * t113 + t195 * t156 + (qJDD(2) * t165 + t111 * t222) * pkin(2) + t40, g(3) * t156 + t83 * qJD(2) - t130 * t111 + t195 * t157 + (-qJDD(2) * t164 - t113 * t222) * pkin(2) - t41, (t76 + t82) * t113 + (t75 - t83) * t111 + (-t164 * t77 - t165 * t78) * pkin(2), -t75 * t82 - t76 * t83 + (-t249 + t164 * t41 + t165 * t40 + (-qJD(1) * t130 + t195) * t169) * pkin(2), -t272, t267, t279, t269, t159, t159 * t188 - t160 * t239 + t68 * t87 + t258, -t159 * t224 + t160 * t240 - t191 * t87 + t266, t277, -t237 * t271 + t278, t276, t190 + t243, -t261, t106 * t14 + t239 * t52 + (-t274 - t263) * t171 + t183 * t167 + t260, t106 * t13 + t239 * t54 + t183 * t171 + (-t280 + t263) * t167 + t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t113 * qJD(2) + t192, (t111 + t210) * qJD(2) + t182, -t111 ^ 2 - t113 ^ 2, -t76 * t111 + t75 * t113 + t105 - t194, 0, 0, 0, 0, 0, t26 + t234, t25 + t233, 0, 0, 0, 0, 0, t190 - t243, -t171 * t271 ^ 2 - t18 - t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t272, t267, t279, t269, t159, t29 * t160 + t258, t28 * t160 + t266, t277, -t197 * t54 + t278, t276, -t197 * t271 + t19 + t243, -t261, -pkin(4) * t14 - t29 * t52 + t185 * t167 + (-t274 - t262) * t171 + t260, -pkin(4) * t13 - t29 * t54 + t185 * t171 + (-t280 + t262) * t167 + t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t52, -t52 ^ 2 + t54 ^ 2, t271 * t52 + t13, t271 * t54 - t14, t24, -g(1) * t103 + g(2) * t101 + t167 * t189 - t20 * t54 - t21 * t219 + t271 * t7 + t4, g(1) * t104 - g(2) * t102 + t167 * t205 + t171 * t189 + t20 * t52 + t271 * t6;];
tau_reg = t1;

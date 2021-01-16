% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:06:19
% EndTime: 2021-01-16 01:06:34
% DurationCPUTime: 4.00s
% Computational Cost: add. (3253->409), mult. (7213->584), div. (0->0), fcn. (5939->22), ass. (0->225)
t170 = cos(qJ(4));
t270 = cos(pkin(12));
t218 = t270 * t170;
t132 = qJD(2) * t218;
t158 = sin(pkin(12));
t167 = sin(qJ(4));
t245 = qJD(2) * t167;
t98 = t158 * t245 - t132;
t97 = qJD(6) + t98;
t305 = t97 - qJD(6);
t304 = qJ(5) + pkin(8);
t168 = sin(qJ(2));
t171 = cos(qJ(2));
t161 = sin(pkin(6));
t246 = qJD(2) * t161;
t222 = qJD(1) * t246;
t240 = qJDD(1) * t161;
t303 = t168 * t240 + t171 * t222;
t155 = qJ(2) + pkin(11);
t150 = cos(t155);
t160 = sin(pkin(10));
t125 = t160 * t150;
t163 = cos(pkin(10));
t126 = t163 * t150;
t148 = sin(t155);
t164 = cos(pkin(6));
t255 = t163 * t164;
t263 = t160 * t164;
t292 = -g(2) * (t148 * t255 + t125) + g(1) * (t148 * t263 - t126);
t154 = qJ(4) + pkin(12);
t149 = cos(t154);
t288 = g(3) * t161;
t237 = t150 * t288;
t110 = t158 * t170 + t270 * t167;
t144 = pkin(4) * t170 + pkin(3);
t162 = cos(pkin(11));
t290 = pkin(2) * t162;
t120 = -t144 - t290;
t188 = -t158 * t167 + t218;
t48 = -pkin(5) * t188 - pkin(9) * t110 + t120;
t100 = t110 * qJD(4);
t239 = t167 * qJDD(2);
t203 = -qJDD(2) * t218 + t158 * t239;
t56 = qJD(2) * t100 + t203;
t55 = qJDD(6) + t56;
t302 = ((g(1) * t163 + g(2) * t160) * t148 - t237) * t149 + t48 * t55;
t133 = t164 * qJDD(1) + qJDD(3);
t116 = t170 * t133;
t135 = qJD(1) * t164 + qJD(3);
t159 = sin(pkin(11));
t130 = t171 * t240;
t88 = qJDD(2) * pkin(2) - t168 * t222 + t130;
t44 = t159 * t88 + t303 * t162;
t40 = qJDD(2) * pkin(8) + t44;
t182 = qJ(5) * qJDD(2) + qJD(2) * qJD(5) + qJD(4) * t135 + t40;
t247 = qJD(1) * t161;
t223 = t171 * t247;
t113 = qJD(2) * pkin(2) + t223;
t224 = t168 * t247;
t118 = t162 * t224;
t69 = t159 * t113 + t118;
t212 = t304 * qJD(2) + t69;
t197 = t212 * qJD(4);
t12 = qJDD(4) * pkin(4) - t167 * t182 - t170 * t197 + t116;
t13 = (t133 - t197) * t167 + t182 * t170;
t3 = t270 * t12 - t158 * t13;
t1 = -qJDD(4) * pkin(5) - t3;
t101 = t110 * qJD(2);
t127 = t164 * t149;
t139 = pkin(4) * t158 + pkin(9);
t147 = sin(t154);
t266 = t148 * t161;
t205 = g(1) * t160 - g(2) * t163;
t297 = t205 * t161;
t301 = t292 * t147 + t149 * t297 + (pkin(4) * t245 + pkin(5) * t101 + pkin(9) * t98 + qJD(6) * t139) * t97 + g(3) * (-t147 * t266 + t127) + t1;
t169 = cos(qJ(6));
t166 = sin(qJ(6));
t243 = qJD(6) * t166;
t194 = -t169 * t55 + t97 * t243;
t279 = qJD(4) * pkin(4);
t80 = t159 * t223 + t118;
t298 = t167 * t279 - t80;
t138 = t164 * t170;
t94 = (t159 * t171 + t162 * t168) * t161;
t217 = -t167 * t94 + t138;
t241 = qJD(2) * qJD(4);
t221 = t167 * t241;
t179 = qJDD(2) * t110 - t158 * t221;
t57 = qJD(4) * t132 + t179;
t73 = qJD(4) * t166 + t101 * t169;
t25 = qJD(6) * t73 - t169 * qJDD(4) + t166 * t57;
t103 = t188 * qJD(4);
t42 = t135 * t167 + t170 * t212;
t276 = t158 * t42;
t41 = t170 * t135 - t167 * t212;
t36 = t41 + t279;
t18 = t270 * t36 - t276;
t16 = -qJD(4) * pkin(5) - t18;
t4 = t158 * t12 + t270 * t13;
t2 = qJDD(4) * pkin(9) + t4;
t117 = t159 * t224;
t68 = t113 * t162 - t117;
t54 = -qJD(2) * t144 + qJD(5) - t68;
t23 = pkin(5) * t98 - pkin(9) * t101 + t54;
t140 = pkin(2) * t159 + pkin(8);
t250 = qJ(5) + t140;
t208 = qJD(4) * t250;
t185 = -qJD(5) * t167 - t170 * t208;
t79 = qJD(5) * t170 - t167 * t208;
t82 = t162 * t223 - t117;
t281 = t158 * t185 - t188 * t82 + t270 * t79;
t108 = t250 * t170;
t214 = t250 * t167;
t52 = t270 * t108 - t158 * t214;
t295 = (qJD(6) * t23 + t2) * t188 + t1 * t110 + t16 * t103 + (t164 * t205 - t288) * t148 + (-qJD(6) * t48 - t281) * t97 - t52 * t55;
t269 = t103 * t169;
t294 = -t110 * t194 + t97 * t269;
t146 = pkin(6) - t155;
t291 = sin(t146);
t289 = pkin(4) * t167;
t287 = g(3) * t164;
t242 = t169 * qJD(4);
t71 = t101 * t166 - t242;
t285 = t71 * t97;
t284 = t73 * t97;
t24 = qJD(6) * t242 + t166 * qJDD(4) - t101 * t243 + t169 * t57;
t283 = t73 * t100 - t188 * t24;
t34 = t270 * t42;
t19 = t158 * t36 + t34;
t282 = -t110 * t82 + t158 * t79 - t270 * t185;
t280 = pkin(5) * t100 - pkin(9) * t103 + t298;
t278 = t101 * t71;
t277 = t101 * t73;
t275 = t16 * t110;
t274 = t166 * t24;
t273 = t166 * t55;
t272 = t166 * t97;
t215 = t169 * t97;
t268 = t147 * t164;
t267 = t148 * t160;
t265 = t149 * t150;
t264 = t160 * t161;
t262 = t160 * t168;
t261 = t161 * t163;
t260 = t161 * t166;
t259 = t161 * t169;
t258 = t161 * t170;
t257 = t161 * t171;
t256 = t163 * t148;
t254 = t164 * t166;
t253 = t164 * t168;
t252 = t164 * t169;
t251 = t164 * t171;
t249 = qJDD(1) - g(3);
t156 = t167 ^ 2;
t248 = -t170 ^ 2 + t156;
t244 = qJD(4) * t170;
t238 = t170 * qJDD(2);
t236 = cos(t146) / 0.2e1;
t229 = t166 * t265;
t228 = t169 * t265;
t227 = t149 * t254;
t226 = t149 * t252;
t225 = t163 * t251;
t213 = pkin(6) + t155;
t195 = sin(t213) / 0.2e1;
t104 = t195 - t291 / 0.2e1;
t210 = t104 * t163 + t125;
t209 = -t104 * t160 + t126;
t43 = -t303 * t159 + t162 * t88;
t204 = cos(t213);
t202 = -t100 * t71 + t188 * t25;
t17 = qJD(4) * pkin(9) + t19;
t6 = t166 * t23 + t169 * t17;
t201 = t166 * t17 - t169 * t23;
t59 = t164 * t167 + t170 * t94;
t27 = t158 * t217 + t270 * t59;
t198 = t159 * t168 - t162 * t171;
t93 = t198 * t161;
t200 = t166 * t93 + t169 * t27;
t199 = -t166 * t27 + t169 * t93;
t196 = -t98 * t272 - t194;
t189 = -t160 * t251 - t163 * t168;
t105 = t291 / 0.2e1 + t195;
t184 = t204 / 0.2e1 + t236;
t62 = -t163 * t184 + t267;
t65 = t160 * t184 + t256;
t187 = g(1) * t65 + g(2) * t62 - g(3) * t105;
t186 = -g(1) * t264 + g(2) * t261 - t287;
t21 = t270 * t41 - t276;
t183 = -t139 * t55 + (t16 + t21) * t97;
t142 = -pkin(3) - t290;
t60 = -qJD(2) * pkin(3) - t68;
t181 = -qJDD(4) * t140 + (qJD(2) * t142 + t60 + t82) * qJD(4);
t83 = t198 * t246;
t180 = -qJD(4) * t59 + t167 * t83;
t178 = g(3) * t266 - t60 * qJD(2) - t292 - t40;
t177 = -g(1) * t189 - g(3) * t257;
t176 = (-qJD(6) * t215 - t273) * t110 - t103 * t272;
t28 = pkin(4) * t221 - qJDD(2) * t144 + qJDD(5) - t43;
t172 = qJD(4) ^ 2;
t175 = -g(1) * (t150 * t263 + t256) + g(2) * (t150 * t255 - t267) - qJD(2) * t80 + t140 * t172 + t237 - t43 + (-pkin(3) + t142) * qJDD(2);
t173 = qJD(2) ^ 2;
t141 = -t270 * pkin(4) - pkin(5);
t123 = pkin(2) * t225;
t122 = qJDD(4) * t170 - t167 * t172;
t121 = t167 * qJDD(4) + t170 * t172;
t106 = t236 - t204 / 0.2e1;
t87 = t160 * t226 - t163 * t166;
t86 = t160 * t166 + t163 * t226;
t85 = t160 * t227 + t163 * t169;
t84 = t160 * t169 - t163 * t227;
t81 = qJD(2) * t94;
t78 = t149 * t266 + t268;
t76 = t147 * t259 + t150 * t254;
t75 = -t147 * t260 + t150 * t252;
t51 = t108 * t158 + t270 * t214;
t30 = t217 * qJD(4) - t170 * t83;
t26 = t158 * t59 - t270 * t217;
t20 = t158 * t41 + t34;
t15 = t158 * t180 + t270 * t30;
t14 = t158 * t30 - t270 * t180;
t10 = pkin(5) * t56 - pkin(9) * t57 + t28;
t7 = t169 * t10;
t5 = [t249, 0, (qJDD(2) * t171 - t168 * t173) * t161, (-qJDD(2) * t168 - t171 * t173) * t161, t133 * t164 - t43 * t93 + t44 * t94 - t68 * t81 - t69 * t83 - g(3), 0, 0, 0, 0, 0, t217 * qJDD(4) + (-t81 * qJD(2) - t93 * qJDD(2)) * t170 + (-t94 * t244 + (t93 * qJD(2) - qJD(4) * t164 + t83) * t167) * qJD(4), t93 * t239 - qJD(4) * t30 - qJDD(4) * t59 + (t167 * t81 + t244 * t93) * qJD(2), -qJD(4) * t14 - qJDD(4) * t26 + t56 * t93 + t81 * t98, -qJD(4) * t15 - qJDD(4) * t27 + t101 * t81 + t57 * t93, t101 * t14 - t15 * t98 + t26 * t57 - t27 * t56, -t14 * t18 + t15 * t19 - t26 * t3 + t27 * t4 + t28 * t93 + t54 * t81 - g(3), 0, 0, 0, 0, 0, (-qJD(6) * t200 - t15 * t166 + t169 * t81) * t97 + t199 * t55 + t14 * t71 + t26 * t25, -(qJD(6) * t199 + t15 * t169 + t166 * t81) * t97 - t200 * t55 + t14 * t73 + t26 * t24; 0, qJDD(2), t130 - g(2) * (t225 - t262) + t177, -g(1) * (t160 * t253 - t163 * t171) - g(2) * (-t160 * t171 - t163 * t253) - t249 * t168 * t161, -g(2) * t123 + t68 * t80 - t69 * t82 + (g(2) * t262 + t44 * t159 + t43 * t162 + t177) * pkin(2), qJDD(2) * t156 + 0.2e1 * t170 * t221, 0.2e1 * t167 * t238 - 0.2e1 * t241 * t248, t121, t122, 0, t167 * t181 - t170 * t175, t167 * t175 + t170 * t181, -qJDD(4) * t51 + t100 * t54 - t188 * t28 + t120 * t56 - t80 * t98 + t187 * t149 + (t98 * t289 - t282) * qJD(4), -qJDD(4) * t52 - t101 * t80 + t103 * t54 + t110 * t28 + t120 * t57 - t187 * t147 + (t101 * t289 - t281) * qJD(4), -g(1) * t209 - g(2) * t210 - g(3) * t106 - t100 * t19 + t282 * t101 - t103 * t18 - t110 * t3 + t188 * t4 - t281 * t98 + t51 * t57 - t52 * t56, t4 * t52 - t3 * t51 + t28 * t120 - g(1) * (pkin(2) * t189 - t144 * t65 + t209 * t304) - g(2) * (-pkin(2) * t262 - t144 * t62 + t210 * t304 + t123) - g(3) * (pkin(2) * t257 + t105 * t144 + t106 * t304) + t298 * t54 + t281 * t19 - t282 * t18, t73 * t269 + (t169 * t24 - t243 * t73) * t110, (-t166 * t73 - t169 * t71) * t103 + (-t274 - t169 * t25 + (t166 * t71 - t169 * t73) * qJD(6)) * t110, t283 + t294, t176 + t202, t100 * t97 - t188 * t55, -t201 * t100 - t7 * t188 + t51 * t25 + t282 * t71 + (g(1) * t87 - g(2) * t86) * t150 + (t280 * t97 + (t17 * t188 - t52 * t97 + t275) * qJD(6) + t302) * t169 + t295 * t166, -t6 * t100 + t51 * t24 + t282 * t73 + (-g(1) * t85 - g(2) * t84) * t150 + ((-qJD(6) * t17 + t10) * t188 - qJD(6) * t275 + (qJD(6) * t52 - t280) * t97 - t302) * t166 + t295 * t169; 0, 0, 0, 0, t186 + t133, 0, 0, 0, 0, 0, t122, -t121, -qJD(4) * t100 + qJDD(4) * t188, -qJD(4) * t103 - qJDD(4) * t110, t100 * t101 - t103 * t98 - t110 * t56 - t188 * t57, -t100 * t18 + t103 * t19 + t110 * t4 + t188 * t3 + t186, 0, 0, 0, 0, 0, t176 - t202, t283 - t294; 0, 0, 0, 0, 0, -t167 * t173 * t170, t248 * t173, t239, t238, qJDD(4), -g(3) * t138 + t167 * t178 - t170 * t297 + t116, (-t133 + t297 + t287) * t167 + t178 * t170, t20 * qJD(4) - t54 * t101 - g(1) * (-t147 * t209 + t149 * t264) - g(2) * (-t147 * t210 - t149 * t261) - g(3) * (-t106 * t147 + t127) + (t270 * qJDD(4) - t98 * t245) * pkin(4) + t3, t21 * qJD(4) + t54 * t98 - g(1) * (-t147 * t264 - t149 * t209) - g(2) * (t147 * t261 - t149 * t210) - g(3) * (-t106 * t149 - t268) + (-t158 * qJDD(4) - t101 * t245) * pkin(4) - t4, (-t18 + t21) * t98 + (t19 - t20) * t101 + (-t158 * t56 - t270 * t57) * pkin(4), t18 * t20 - t19 * t21 + (t4 * t158 + t3 * t270 - t54 * t245 - g(1) * (t160 * t258 - t167 * t209) - g(2) * (-t163 * t258 - t167 * t210) - g(3) * (-t106 * t167 + t138)) * pkin(4), t215 * t73 + t274, (t24 - t285) * t169 + (-t25 - t284) * t166, t215 * t97 + t273 - t277, t196 + t278, -t97 * t101, t201 * t101 + t141 * t25 + t183 * t166 - t301 * t169 - t20 * t71, t6 * t101 + t141 * t24 + t301 * t166 + t183 * t169 - t20 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t101 * qJD(4) + t203, (t132 - t98) * qJD(4) + t179, -t101 ^ 2 - t98 ^ 2, t101 * t18 + t19 * t98 - t187 + t28, 0, 0, 0, 0, 0, t196 - t278, -t97 ^ 2 * t169 - t273 - t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, t24 + t285, -t25 + t284, t55, -t166 * t2 + t7 - t16 * t73 - g(1) * (t148 * t85 + t160 * t75 - t163 * t229) - g(2) * (t148 * t84 - t160 * t229 - t163 * t75) - g(3) * (-t150 * t259 - t166 * t78) + t305 * t6, -t169 * t2 - t166 * t10 + t16 * t71 - g(1) * (t148 * t87 - t160 * t76 - t163 * t228) - g(2) * (-t148 * t86 - t160 * t228 + t163 * t76) - g(3) * (t150 * t260 - t169 * t78) - t305 * t201;];
tau_reg = t5;

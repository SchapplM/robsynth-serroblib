% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:55
% EndTime: 2019-03-09 02:04:01
% DurationCPUTime: 3.35s
% Computational Cost: add. (4094->445), mult. (7326->520), div. (0->0), fcn. (4305->10), ass. (0->233)
t143 = cos(qJ(4));
t130 = t143 * pkin(8);
t140 = sin(qJ(4));
t311 = pkin(4) * t140 - t130;
t247 = qJD(1) * t140;
t108 = qJD(5) + t247;
t139 = sin(qJ(5));
t204 = qJD(5) * t140 + qJD(1);
t142 = cos(qJ(5));
t238 = t142 * qJD(4);
t216 = t143 * t238;
t239 = qJD(5) * t143;
t214 = t139 * t239;
t165 = t140 * t238 + t214;
t229 = t143 * qJDD(1);
t41 = qJD(1) * t165 - qJD(5) * t238 - t139 * qJDD(4) - t142 * t229;
t267 = t143 * t41;
t237 = qJD(1) * qJD(4);
t211 = t143 * t237;
t230 = t140 * qJDD(1);
t84 = qJDD(5) + t211 + t230;
t270 = t142 * t84;
t244 = qJD(4) * t139;
t246 = qJD(1) * t143;
t89 = t142 * t246 + t244;
t310 = t108 * (t139 * t204 - t216) + (qJD(4) * t89 - t270) * t140 + t267;
t180 = t108 * t89;
t87 = t139 * t246 - t238;
t181 = t108 * t87;
t243 = qJD(4) * t140;
t218 = t139 * t243;
t42 = -qJD(1) * t218 + qJD(5) * t89 - t142 * qJDD(4) + t139 * t229;
t309 = (t42 + t180) * t139 + t142 * (t41 + t181);
t307 = t42 - t180;
t73 = t87 * t243;
t305 = t143 * t42 - t73;
t263 = pkin(1) * qJDD(1);
t304 = t142 * t204;
t132 = qJ(1) + pkin(9);
t123 = cos(t132);
t116 = g(2) * t123;
t122 = sin(t132);
t117 = g(1) * t122;
t303 = t117 - t116;
t137 = sin(pkin(9));
t113 = pkin(1) * t137 + qJ(3);
t78 = t113 + t311;
t138 = cos(pkin(9));
t118 = -pkin(1) * t138 - pkin(2);
t302 = t118 * qJDD(1);
t245 = qJD(2) * t143;
t107 = -pkin(7) + t118;
t82 = t107 * qJD(1) + qJD(3);
t57 = t140 * t82 + t245;
t54 = qJD(4) * pkin(8) + t57;
t58 = t78 * qJD(1);
t20 = -t139 * t54 + t142 * t58;
t253 = qJD(6) - t20;
t13 = -pkin(5) * t108 + t253;
t21 = t139 * t58 + t142 * t54;
t14 = qJ(6) * t108 + t21;
t189 = t13 * t139 + t14 * t142;
t236 = qJD(2) * qJD(4);
t79 = t107 * qJDD(1) + qJDD(3);
t26 = -t140 * qJDD(2) - t82 * t243 + (-t236 + t79) * t143;
t24 = -qJDD(4) * pkin(4) - t26;
t5 = pkin(5) * t42 + qJ(6) * t41 - qJD(6) * t89 + t24;
t301 = qJD(4) * t189 - t5;
t254 = t142 * t143;
t66 = t84 * t254;
t300 = -t108 * t165 + t66;
t56 = -t140 * qJD(2) + t143 * t82;
t53 = -qJD(4) * pkin(4) - t56;
t19 = pkin(5) * t87 - qJ(6) * t89 + t53;
t287 = pkin(8) * t84;
t299 = t108 * t19 - t287;
t240 = qJD(5) * t142;
t170 = t108 * t240 + t139 * t84;
t81 = t108 * t218;
t298 = -t143 * t170 + t81;
t187 = t139 * t20 - t142 * t21;
t297 = qJD(4) * t187 + t24;
t225 = pkin(8) * qJD(5) * t108;
t284 = g(3) * t140;
t294 = t143 * t303 + t225 - t284;
t268 = t142 * t87;
t273 = t139 * t89;
t213 = t142 * t239;
t256 = t139 * t143;
t280 = t89 * t213 - t41 * t256;
t293 = -(t268 + t273) * t243 + t280;
t241 = qJD(5) * t139;
t257 = t139 * t140;
t292 = qJD(1) * (t108 * t257 - t143 * t87) + t108 * t241 - t270;
t242 = qJD(4) * t143;
t226 = -t143 * qJDD(2) - t140 * t79 - t82 * t242;
t25 = -t140 * t236 - t226;
t155 = -(t140 * t56 - t143 * t57) * qJD(4) + t25 * t140 + t26 * t143;
t96 = qJD(1) * t113;
t291 = 0.2e1 * qJD(4) * t96 + qJDD(4) * t107;
t255 = t140 * t142;
t277 = t107 * t255 + t139 * t78;
t199 = pkin(4) * t143 + pkin(8) * t140;
t85 = qJD(4) * t199 + qJD(3);
t290 = -qJD(5) * t277 + t142 * t85;
t289 = t89 ^ 2;
t288 = pkin(5) * t84;
t286 = pkin(8) * t89;
t283 = g(3) * t143;
t282 = t89 * t87;
t37 = t140 * t42;
t281 = t87 * t242 + t37;
t91 = t199 * qJD(1);
t35 = t139 * t91 + t142 * t56;
t38 = t140 * t41;
t72 = t89 * t242;
t279 = t72 - t38;
t192 = pkin(5) * t139 - qJ(6) * t142;
t278 = qJD(5) * t192 - qJD(6) * t139 - t245 - (-qJD(1) * t192 + t82) * t140;
t276 = qJ(6) * t84;
t275 = t108 * t14;
t274 = t108 * t21;
t272 = t142 * t42;
t271 = t142 * t78;
t265 = t143 * t89;
t258 = t123 * t143;
t227 = g(2) * t258;
t264 = g(3) * t255 + t142 * t227;
t262 = qJD(1) * t96;
t261 = t122 * t140;
t260 = t122 * t143;
t259 = t123 * t140;
t252 = g(1) * t258 + g(2) * t260;
t119 = t137 * t263;
t250 = qJDD(1) * qJ(3) + t119;
t134 = t140 ^ 2;
t135 = t143 ^ 2;
t249 = t134 - t135;
t145 = qJD(4) ^ 2;
t146 = qJD(1) ^ 2;
t248 = -t145 - t146;
t235 = qJD(3) * qJD(1);
t234 = qJDD(1) * t113;
t232 = qJDD(4) * t140;
t231 = qJDD(4) * t143;
t228 = t107 * t216 + t139 * t85 + t78 * t240;
t102 = g(2) * t259;
t224 = t87 ^ 2 - t289;
t223 = t89 * t243;
t222 = t108 * t255;
t221 = t143 * t146 * t140;
t144 = cos(qJ(1));
t219 = t144 * pkin(1) + t123 * pkin(2) + t122 * qJ(3);
t217 = t139 * t242;
t215 = t107 * t241;
t212 = t108 * t246;
t141 = sin(qJ(1));
t210 = -pkin(1) * t141 + t123 * qJ(3);
t209 = t107 * t139 - pkin(5);
t208 = -t102 + t283;
t23 = qJDD(4) * pkin(8) + t25;
t40 = t85 * qJD(1) + t311 * qJDD(1) + t250;
t4 = -t139 * t23 + t142 * t40 - t54 * t240 - t58 * t241;
t207 = qJD(5) * t87 - t41;
t205 = (-t134 - t135) * qJDD(1);
t203 = t123 * pkin(7) + t219;
t202 = t140 * t211;
t61 = t122 * t257 - t123 * t142;
t63 = t122 * t142 + t123 * t257;
t201 = g(1) * t63 + g(2) * t61;
t62 = t122 * t255 + t123 * t139;
t64 = -t122 * t139 + t123 * t255;
t200 = -g(1) * t64 - g(2) * t62;
t198 = g(1) * t123 + g(2) * t122;
t196 = g(1) * t141 - g(2) * t144;
t195 = -t303 - t262;
t194 = t207 * pkin(8);
t193 = pkin(5) * t142 + qJ(6) * t139;
t191 = -t262 - t117;
t190 = t13 * t142 - t139 * t14;
t188 = t139 * t21 + t142 * t20;
t34 = -t139 * t56 + t142 * t91;
t83 = t235 + t250;
t183 = t96 * qJD(3) + t83 * t113;
t30 = t42 * t254;
t179 = t87 * t214 - t30;
t177 = pkin(4) + t193;
t175 = -g(1) * (pkin(4) * t260 + pkin(8) * t261) - g(3) * t130;
t174 = -t107 + t192;
t3 = t139 * t40 + t142 * t23 + t58 * t240 - t54 * t241;
t168 = qJDD(3) + t302;
t167 = -g(1) * t260 - t225;
t166 = (-pkin(2) - pkin(7)) * t122 + t210;
t164 = -pkin(8) * t272 - g(1) * t261 - t208;
t163 = t108 * t53 - t287;
t162 = pkin(4) * t261 - pkin(8) * t260 + t203;
t161 = g(1) * t61 - g(2) * t63 + g(3) * t256 + t4;
t160 = -t107 * t145 + t234 + t235 + t83;
t159 = t139 * t181 - t272;
t158 = pkin(4) * t259 - pkin(8) * t258 + t166;
t1 = qJD(6) * t108 + t276 + t3;
t2 = qJDD(6) - t4 - t288;
t157 = qJD(5) * t190 + t1 * t142 + t2 * t139;
t156 = -qJD(5) * t188 - t4 * t139 + t3 * t142;
t154 = t19 * t89 + qJDD(6) - t161;
t153 = t139 * t305 + t87 * t213;
t152 = -g(1) * t62 + g(2) * t64 - g(3) * t254 + t3;
t151 = qJD(4) * t19 + t157;
t150 = qJD(4) * t53 + t156;
t149 = -t84 * t257 + (-t217 - t304) * t108 - t305;
t148 = -t42 * t255 - t87 * t216 + t89 * t304 + (qJD(1) * t87 + t207 * t140 + t72) * t139;
t136 = qJDD(2) - g(3);
t95 = -t140 * t145 + t231;
t94 = -t143 * t145 - t232;
t55 = t174 * t143;
t49 = t108 * t242 + t140 * t84;
t48 = pkin(5) * t89 + qJ(6) * t87;
t45 = -t107 * t257 + t271;
t44 = t209 * t140 - t271;
t43 = qJ(6) * t140 + t277;
t28 = -pkin(5) * t246 - t34;
t27 = qJ(6) * t246 + t35;
t22 = (qJD(5) * t193 - qJD(6) * t142) * t143 - t174 * t243;
t16 = t181 - t41;
t15 = (t222 - t265) * qJD(1) + t170;
t12 = -t107 * t217 + t290;
t11 = -t140 * t215 + t228;
t10 = t209 * t242 - t290;
t9 = qJ(6) * t242 + (qJD(6) - t215) * t140 + t228;
t8 = -t139 * t41 + t142 * t180;
t7 = -t89 * t214 + (-t223 - t267) * t142;
t6 = t279 + t300;
t17 = [0, 0, 0, 0, 0, qJDD(1), t196, g(1) * t144 + g(2) * t141, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t138 * t263 + t303, -0.2e1 * t119 + t198, 0 (t196 + (t137 ^ 2 + t138 ^ 2) * t263) * pkin(1), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - t303 + 0.2e1 * t302, -t198 + t234 + 0.2e1 * t235 + t250, t168 * t118 - g(1) * (-pkin(2) * t122 + t210) - g(2) * t219 + t183, qJDD(1) * t135 - 0.2e1 * t202, -0.2e1 * t140 * t229 + 0.2e1 * t237 * t249, t95, qJDD(1) * t134 + 0.2e1 * t202, t94, 0, t291 * t143 + (t160 - t198) * t140, -t140 * t291 + t160 * t143 - t252, t107 * t205 - t155 + t303, -g(1) * t166 - g(2) * t203 + t107 * t155 + t183, t7, t179 - t293, t6, t153, -t37 + t81 + (-qJD(4) * t87 - t170) * t143, t49, t108 * t12 + t45 * t84 + (t4 + (t107 * t87 - t139 * t53) * qJD(4)) * t140 + (qJD(4) * t20 - t107 * t42 + t139 * t24 + t240 * t53) * t143 + t200, -t108 * t11 - t277 * t84 + (-t3 + (t107 * t89 - t142 * t53) * qJD(4)) * t140 + (-qJD(4) * t21 + t107 * t41 + t142 * t24 - t241 * t53) * t143 + t201, -t11 * t87 - t12 * t89 + t41 * t45 - t42 * t277 + t188 * t243 + (qJD(5) * t187 - t139 * t3 - t142 * t4) * t143 + t252, t3 * t277 + t21 * t11 + t4 * t45 + t20 * t12 - g(1) * t158 - g(2) * t162 + (-t143 * t24 + t243 * t53) * t107, t7, t6 (-t241 * t87 + t272) * t143 + t293, t49, t281 - t298, t153, -t10 * t108 + t22 * t87 + t42 * t55 - t44 * t84 + (-t19 * t244 - t2) * t140 + (-qJD(4) * t13 + t139 * t5 + t19 * t240) * t143 + t200, t10 * t89 - t41 * t44 - t42 * t43 - t87 * t9 - t190 * t243 + (-qJD(5) * t189 - t1 * t139 + t142 * t2) * t143 + t252, t108 * t9 - t22 * t89 + t41 * t55 + t43 * t84 + (t19 * t238 + t1) * t140 + (qJD(4) * t14 - t142 * t5 + t19 * t241) * t143 - t201, t1 * t43 + t14 * t9 + t5 * t55 + t19 * t22 + t2 * t44 + t13 * t10 - g(1) * (pkin(5) * t64 + qJ(6) * t63 + t158) - g(2) * (pkin(5) * t62 + qJ(6) * t61 + t162); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, 0, 0, 0, 0, 0, t94, -t95, 0, -t140 * t26 + t143 * t25 - g(3) + (-t140 * t57 - t143 * t56) * qJD(4), 0, 0, 0, 0, 0, 0, t281 + t298, t279 - t300, -t30 + (t239 * t89 + t73) * t142 + (t143 * t207 - t223) * t139, t140 * t297 + t150 * t143 - g(3), 0, 0, 0, 0, 0, 0, -t84 * t256 + (-t213 + t218) * t108 + t281 (t268 - t273) * t243 + t179 + t280, -t108 * t214 + t38 + t66 + (-t222 - t265) * qJD(4), -t140 * t301 + t151 * t143 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t146, t168 + t195, 0, 0, 0, 0, 0, 0, t140 * t248 + t231, t143 * t248 - t232, t205, t155 + t195, 0, 0, 0, 0, 0, 0, t149, t310, t148, -t188 * qJD(1) + t150 * t140 - t143 * t297 - t303, 0, 0, 0, 0, 0, 0, t149, t148, -t310, t190 * qJD(1) + t151 * t140 + t143 * t301 - t303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, -t249 * t146, t229, -t221, -t230, qJDD(4), t284 + qJD(4) * t57 + (t191 + t116) * t143 + t26, qJD(4) * t56 + (-t191 + t236) * t140 + t208 + t226, 0, 0, t8, -t309, t15, t159, -t292, -t212, -t20 * t246 - pkin(4) * t42 - t108 * t34 - t57 * t87 + (t167 - t24) * t142 + t163 * t139 + t264, t21 * t246 + pkin(4) * t41 + t108 * t35 - t57 * t89 + t163 * t142 + (t24 + t294) * t139, t34 * t89 + t35 * t87 + (-t20 * t247 + t3 + (-t20 + t286) * qJD(5)) * t142 + (t194 - t4 - t274) * t139 + t164, -t21 * t35 - t20 * t34 - t53 * t57 + (t227 - t24 + t284) * pkin(4) + (t156 + t102) * pkin(8) + t175, t8, t15, t309, -t212, t292, t159, t13 * t246 + t108 * t28 - t42 * t177 + t278 * t87 + (t167 - t5) * t142 + t299 * t139 + t264, t27 * t87 - t28 * t89 + (t13 * t247 + t1 + (t13 + t286) * qJD(5)) * t142 + (t194 + t2 - t275) * t139 + t164, -t14 * t246 - t108 * t27 - t41 * t177 - t278 * t89 - t299 * t142 + (-t294 - t5) * t139, -t14 * t27 - t13 * t28 + t278 * t19 - t117 * t193 * t143 + (t157 + t102) * pkin(8) + t175 + (t116 * t143 + t284 - t5) * t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t282, -t224, t16, -t282, -t307, t84, -t53 * t89 + t161 + t274, t108 * t20 + t53 * t87 - t152, 0, 0, t282, t16, t224, t84, t307, -t282, -t48 * t87 - t154 + t274 + 0.2e1 * t288, pkin(5) * t41 - qJ(6) * t42 + (t14 - t21) * t89 + (t13 - t253) * t87, 0.2e1 * t276 - t19 * t87 + t48 * t89 + (0.2e1 * qJD(6) - t20) * t108 + t152, t1 * qJ(6) - t2 * pkin(5) - t19 * t48 - t13 * t21 - g(1) * (-pkin(5) * t61 + qJ(6) * t62) - g(2) * (pkin(5) * t63 - qJ(6) * t64) + t192 * t283 + t253 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84 + t282, t16, -t108 ^ 2 - t289, t154 - t275 - t288;];
tau_reg  = t17;

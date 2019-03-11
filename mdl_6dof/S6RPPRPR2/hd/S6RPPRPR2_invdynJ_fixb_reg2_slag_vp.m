% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:42
% EndTime: 2019-03-09 01:42:48
% DurationCPUTime: 3.59s
% Computational Cost: add. (5241->415), mult. (11072->490), div. (0->0), fcn. (8011->14), ass. (0->217)
t157 = sin(pkin(10));
t159 = cos(pkin(10));
t163 = sin(qJ(4));
t282 = cos(qJ(4));
t116 = t282 * t157 + t163 * t159;
t291 = t116 * qJD(1);
t298 = qJD(6) + t291;
t221 = t282 * t159;
t208 = qJD(1) * t221;
t241 = t163 * t157;
t220 = qJD(1) * t241;
t105 = -t208 + t220;
t162 = sin(qJ(6));
t165 = cos(qJ(6));
t82 = qJD(4) * t162 - t165 * t105;
t214 = t298 * t82;
t233 = qJD(6) * t165;
t234 = qJD(6) * t162;
t110 = t116 * qJD(4);
t217 = qJDD(1) * t282;
t228 = t157 * qJDD(1);
t198 = -t159 * t217 + t163 * t228;
t62 = qJD(1) * t110 + t198;
t29 = qJD(4) * t234 - t165 * qJDD(4) - t105 * t233 - t162 * t62;
t302 = t29 - t214;
t158 = sin(pkin(9));
t136 = pkin(1) * t158 + qJ(3);
t127 = t136 * qJD(1);
t146 = t159 * qJD(2);
t268 = pkin(7) * qJD(1);
t80 = t146 + (-t127 - t268) * t157;
t94 = t157 * qJD(2) + t159 * t127;
t81 = t159 * t268 + t94;
t35 = t163 * t81 - t282 * t80;
t240 = -qJD(5) - t35;
t239 = pkin(5) * t291 - t240;
t284 = pkin(4) + pkin(8);
t22 = -t284 * qJD(4) + t239;
t141 = t159 * pkin(3) + pkin(2);
t160 = cos(pkin(9));
t275 = t160 * pkin(1);
t126 = -t141 - t275;
t101 = qJD(1) * t126 + qJD(3);
t175 = -qJ(5) * t291 + t101;
t31 = t284 * t105 + t175;
t8 = -t162 * t31 + t165 * t22;
t306 = t8 * t298;
t142 = -pkin(2) - t275;
t230 = qJDD(1) * t142;
t121 = qJDD(3) + t230;
t156 = qJ(1) + pkin(9);
t148 = sin(t156);
t150 = cos(t156);
t218 = -g(1) * t148 + g(2) * t150;
t305 = -t121 - t218;
t227 = t159 * qJDD(1);
t223 = qJD(4) * t208 + t157 * t217 + t163 * t227;
t61 = qJD(4) * t220 - t223;
t261 = t61 * qJ(5);
t98 = qJDD(1) * t126 + qJDD(3);
t170 = -qJD(5) * t291 + t261 + t98;
t10 = t284 * t62 + t170;
t219 = qJD(4) * t282;
t235 = qJD(4) * t163;
t113 = qJD(1) * qJD(3) + t136 * qJDD(1);
t144 = t159 * qJDD(2);
t77 = t144 + (-pkin(7) * qJDD(1) - t113) * t157;
t86 = t157 * qJDD(2) + t159 * t113;
t78 = pkin(7) * t227 + t86;
t215 = t163 * t78 + t81 * t219 + t80 * t235 - t282 * t77;
t200 = qJDD(5) + t215;
t6 = -t61 * pkin(5) - t284 * qJDD(4) + t200;
t9 = t162 * t22 + t165 * t31;
t2 = -qJD(6) * t9 - t162 * t10 + t165 * t6;
t290 = t298 * t9 + t2;
t84 = qJD(4) * t165 + t105 * t162;
t30 = qJD(6) * t84 + t162 * qJDD(4) - t165 * t62;
t304 = t298 * t84 - t30;
t295 = t165 * t298;
t303 = t84 * t295;
t213 = t162 * t298;
t60 = -qJDD(6) + t61;
t53 = t165 * t60;
t187 = -t213 * t298 - t53;
t155 = pkin(10) + qJ(4);
t147 = sin(t155);
t149 = cos(t155);
t206 = g(1) * t150 + g(2) * t148;
t178 = -g(3) * t147 - t206 * t149;
t216 = -t163 * t77 - t80 * t219 + t81 * t235 - t282 * t78;
t301 = t178 - t216;
t115 = -t221 + t241;
t192 = -t110 * t291 + t115 * t61;
t109 = t157 * t235 - t159 * t219;
t294 = t109 * t105 - t116 * t62;
t300 = -t192 + t294;
t299 = t192 + t294;
t297 = 0.2e1 * t291 * qJD(4) + t198;
t258 = pkin(1) * qJDD(1);
t237 = t149 * pkin(4) + t147 * qJ(5);
t104 = t291 ^ 2;
t285 = t105 ^ 2;
t293 = -t285 - t104;
t292 = -t285 + t104;
t276 = t105 * pkin(5);
t36 = t163 * t80 + t282 * t81;
t34 = -qJD(4) * qJ(5) - t36;
t23 = -t34 - t276;
t289 = t23 * t298 + t284 * t60;
t251 = t110 * t162;
t265 = t162 * t60;
t288 = -t115 * (t233 * t298 - t265) - t298 * t251;
t250 = t110 * t165;
t264 = t165 * t29;
t287 = -t115 * (t84 * t234 + t264) + t84 * t250;
t272 = pkin(7) + t136;
t111 = t272 * t157;
t112 = t272 * t159;
t37 = (qJD(3) * t157 + qJD(4) * t112) * t163 - qJD(3) * t221 + t111 * t219;
t55 = -t163 * t111 + t282 * t112;
t286 = t37 * qJD(4) - t55 * qJDD(4) + t147 * t218;
t283 = t62 * pkin(4);
t281 = pkin(8) * t149;
t138 = g(3) * t149;
t164 = sin(qJ(1));
t274 = t164 * pkin(1);
t273 = t84 * t82;
t271 = -t84 * t109 - t29 * t116;
t267 = t105 * t82;
t28 = t165 * t30;
t262 = t30 * t162;
t260 = t84 * t105;
t257 = qJD(4) * t36;
t255 = qJDD(4) * pkin(4);
t254 = t105 * qJ(5);
t253 = t105 * t291;
t249 = t147 * t148;
t248 = t147 * t150;
t247 = t148 * t149;
t246 = t148 * t162;
t245 = t148 * t165;
t244 = t149 * t150;
t243 = t150 * t162;
t242 = t150 * t165;
t166 = cos(qJ(1));
t152 = t166 * pkin(1);
t238 = t150 * t141 + t152;
t153 = t157 ^ 2;
t154 = t159 ^ 2;
t236 = t153 + t154;
t229 = qJDD(4) * qJ(5);
t226 = -t82 * t251 + (-t233 * t82 - t262) * t115;
t222 = -g(1) * t248 - g(2) * t249 + t138;
t54 = t282 * t111 + t112 * t163;
t209 = pkin(4) * t244 + qJ(5) * t248 + t238;
t207 = g(1) * t247 - g(2) * t244;
t204 = g(1) * t164 - g(2) * t166;
t203 = t162 * t9 + t165 * t8;
t202 = -t162 * t8 + t165 * t9;
t1 = qJD(6) * t8 + t165 * t10 + t162 * t6;
t201 = t1 - t306;
t161 = -pkin(7) - qJ(3);
t199 = -t150 * t161 - t274;
t197 = -t109 * t82 + t116 * t30;
t85 = -t113 * t157 + t144;
t196 = -t85 * t157 + t86 * t159;
t195 = t157 * (-t127 * t157 + t146) - t159 * t94;
t181 = -t116 * qJ(5) + t126;
t39 = t284 * t115 + t181;
t41 = pkin(5) * t116 + t54;
t18 = t162 * t41 + t165 * t39;
t17 = -t162 * t39 + t165 * t41;
t194 = t105 * t110 + t115 * t62;
t193 = -t109 * t291 - t116 * t61;
t191 = t109 * qJ(5) - t116 * qJD(5);
t66 = qJD(4) * t109 - qJDD(4) * t116;
t65 = qJD(4) * t110 + qJDD(4) * t115;
t188 = -t141 - t237;
t184 = t298 * t250 + (-t234 * t298 - t53) * t115;
t183 = -t215 - t222;
t182 = -t295 * t298 + t265;
t179 = -t230 + t305;
t13 = -qJD(4) * qJD(5) + t216 - t229;
t38 = qJD(3) * t116 + qJD(4) * t55;
t177 = -qJD(4) * t38 - qJDD(4) * t54 + t207;
t45 = t105 * pkin(4) + t175;
t176 = t291 * t45 + qJDD(5) - t183;
t174 = t202 * qJD(6) + t1 * t162 + t2 * t165;
t173 = t218 + t98;
t7 = -pkin(5) * t62 - t13;
t172 = qJD(6) * t284 * t298 + t178 + t7;
t171 = t105 * t37 + t291 * t38 - t54 * t61 - t55 * t62 - t206;
t119 = qJ(5) * t244;
t117 = qJ(5) * t247;
t97 = qJD(4) * t105;
t90 = -t147 * t246 + t242;
t89 = t147 * t245 + t243;
t88 = t147 * t243 + t245;
t87 = t147 * t242 - t246;
t59 = pkin(4) * t291 + t254;
t52 = pkin(4) * t115 + t181;
t46 = pkin(4) * t110 + t191;
t43 = t61 - t97;
t42 = -t115 * pkin(5) + t55;
t40 = t284 * t291 + t254;
t33 = -qJD(4) * pkin(4) - t240;
t32 = t284 * t110 + t191;
t27 = -t109 * pkin(5) + t38;
t26 = -t110 * pkin(5) - t37;
t25 = t36 - t276;
t19 = t170 + t283;
t14 = t200 - t255;
t12 = t162 * t25 + t165 * t40;
t11 = -t162 * t40 + t165 * t25;
t4 = -t18 * qJD(6) - t162 * t32 + t165 * t27;
t3 = t17 * qJD(6) + t162 * t27 + t165 * t32;
t5 = [0, 0, 0, 0, 0, qJDD(1), t204, g(1) * t166 + g(2) * t164, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t160 * t258 - t218, -0.2e1 * t158 * t258 + t206, 0 (t204 + (t158 ^ 2 + t160 ^ 2) * t258) * pkin(1), t153 * qJDD(1), 0.2e1 * t157 * t227, 0, t154 * qJDD(1), 0, 0, t179 * t159, -t179 * t157, t113 * t236 + t196 - t206, t121 * t142 - g(1) * (-pkin(2) * t148 + qJ(3) * t150 - t274) - g(2) * (pkin(2) * t150 + qJ(3) * t148 + t152) + t196 * t136 - t195 * qJD(3), t193, t299, -t66, t194, -t65, 0, t101 * t110 + t115 * t98 + t126 * t62 + t177, -t101 * t109 + t98 * t116 - t126 * t61 + t286, -t109 * t35 - t110 * t36 + t115 * t216 + t116 * t215 + t171, -t216 * t55 - t36 * t37 + t215 * t54 + t35 * t38 + t98 * t126 - g(1) * (-t148 * t141 + t199) - g(2) * (-t148 * t161 + t238) 0, t66, t65, t193, t299, t194, -t109 * t33 + t110 * t34 + t115 * t13 + t116 * t14 + t171, -t105 * t46 - t110 * t45 - t115 * t19 - t52 * t62 - t177, t45 * t109 - t19 * t116 - t291 * t46 + t52 * t61 - t286, t19 * t52 + t45 * t46 - t13 * t55 + t34 * t37 + t14 * t54 + t33 * t38 - g(1) * t199 - g(2) * t209 + (-g(1) * t188 + g(2) * t161) * t148, t84 * t251 + (-t162 * t29 + t233 * t84) * t115, t226 + t287, t271 - t288, -t82 * t250 + (t234 * t82 - t28) * t115, t184 - t197, -t109 * t298 - t116 * t60, -t23 * t250 - g(1) * t90 - g(2) * t88 - t8 * t109 + t2 * t116 - t17 * t60 + t26 * t82 + t42 * t30 + t4 * t298 + (-t7 * t165 + t23 * t234) * t115, t23 * t251 + g(1) * t89 - g(2) * t87 - t1 * t116 + t9 * t109 + t18 * t60 + t26 * t84 - t42 * t29 - t3 * t298 + (t7 * t162 + t23 * t233) * t115, t17 * t29 - t18 * t30 - t3 * t82 - t4 * t84 + t202 * t110 + (-qJD(6) * t203 + t1 * t165 - t2 * t162) * t115 + t207, t1 * t18 + t9 * t3 + t2 * t17 + t8 * t4 + t7 * t42 + t23 * t26 - g(1) * (t150 * pkin(5) + t199) - g(2) * (pkin(8) * t244 + t209) + (-g(1) * (t188 - t281) - g(2) * (pkin(5) - t161)) * t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t157 * t86 + t159 * t85 - g(3), 0, 0, 0, 0, 0, 0, -t65, t66, t300, -t109 * t36 + t110 * t35 + t115 * t215 - t116 * t216 - g(3), 0, 0, 0, 0, 0, 0, t300, t65, -t66, t109 * t34 + t110 * t33 + t115 * t14 - t116 * t13 - g(3), 0, 0, 0, 0, 0, 0, t184 + t197, t271 + t288, t226 - t287, -t23 * t109 + t110 * t203 + t115 * t174 + t7 * t116 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227, t228, -t236 * qJD(1) ^ 2, qJD(1) * t195 - t305, 0, 0, 0, 0, 0, 0, t297 (-t105 - t220) * qJD(4) + t223, t293, t36 * t105 - t291 * t35 + t173, 0, 0, 0, 0, 0, 0, t293, -t297, t61 + t97, t283 + t261 - t34 * t105 + (-qJD(5) - t33) * t291 + t173, 0, 0, 0, 0, 0, 0, t182 + t267, t260 - t187, -t302 * t162 - t28 + t303, t23 * t105 - t290 * t162 + t201 * t165 + t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, t292 (t105 - t220) * qJD(4) + t223, -t253, -t198, qJDD(4), -t101 * t291 + t183 + t257, -t35 * qJD(4) + t101 * t105 - t301, 0, 0, qJDD(4), t43, t198, t253, t292, -t253, pkin(4) * t61 - qJ(5) * t62 + (-t34 - t36) * t291 + (t33 + t240) * t105, t105 * t59 + t176 - 0.2e1 * t255 - t257, 0.2e1 * t229 - t45 * t105 + t59 * t291 + (0.2e1 * qJD(5) + t35) * qJD(4) + t301, -t13 * qJ(5) - t14 * pkin(4) - t45 * t59 - t33 * t36 - g(1) * (-pkin(4) * t248 + t119) - g(2) * (-pkin(4) * t249 + t117) - g(3) * t237 + t240 * t34, -t213 * t84 - t264, -t28 - t303 + (t29 + t214) * t162, t187 + t260, t295 * t82 + t262, t182 - t267, t298 * t105, qJ(5) * t30 + t8 * t105 - t11 * t298 + t172 * t162 + t165 * t289 + t239 * t82, -qJ(5) * t29 - t9 * t105 + t12 * t298 - t162 * t289 + t172 * t165 + t239 * t84, t11 * t84 + t12 * t82 + (-t291 * t9 - t284 * t29 - t2 + (t284 * t82 - t9) * qJD(6)) * t165 + (t291 * t8 + t284 * t30 - t1 + (-t284 * t84 + t8) * qJD(6)) * t162 - t222, t7 * qJ(5) - t9 * t12 - t8 * t11 - g(1) * t119 - g(2) * t117 - g(3) * (t237 + t281) + t239 * t23 + (t206 * t147 - t174) * t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, qJDD(4) - t253, -qJD(4) ^ 2 - t104, qJD(4) * t34 + t176 - t255, 0, 0, 0, 0, 0, 0, -qJD(4) * t82 + t187, -qJD(4) * t84 + t182, t304 * t162 + t302 * t165, -t23 * qJD(4) + t201 * t162 + t165 * t290 + t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, -t82 ^ 2 + t84 ^ 2, -t302, -t273, t304, -t60, -g(1) * t87 - g(2) * t89 + t165 * t138 - t23 * t84 + t290, g(1) * t88 - g(2) * t90 + t23 * t82 + t306 + (-qJD(6) * t22 - t10) * t165 + (qJD(6) * t31 - t138 - t6) * t162, 0, 0;];
tau_reg  = t5;

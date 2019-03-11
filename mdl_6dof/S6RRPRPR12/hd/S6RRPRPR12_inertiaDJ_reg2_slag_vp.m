% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPR12_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:21:48
% EndTime: 2019-03-09 11:22:04
% DurationCPUTime: 6.12s
% Computational Cost: add. (7923->428), mult. (19612->763), div. (0->0), fcn. (18747->10), ass. (0->239)
t271 = pkin(2) + pkin(9);
t108 = sin(pkin(6));
t111 = sin(qJ(2));
t114 = cos(qJ(2));
t252 = cos(pkin(6));
t218 = pkin(1) * t252;
t196 = t114 * t218;
t150 = -t252 * pkin(2) - t196;
t133 = -t252 * pkin(9) + t150;
t270 = pkin(3) + pkin(8);
t225 = t108 * t270;
t129 = t111 * t225 + t133;
t238 = qJD(3) * t111;
t249 = qJ(3) * t114;
t280 = -qJD(4) * t129 - t108 * (-t238 + (t271 * t111 - t249) * qJD(2));
t107 = sin(pkin(11));
t251 = cos(pkin(11));
t113 = cos(qJ(4));
t204 = t251 * t113;
t110 = sin(qJ(4));
t237 = qJD(4) * t110;
t75 = -qJD(4) * t204 + t107 * t237;
t82 = t107 * t113 + t251 * t110;
t76 = t82 * qJD(4);
t279 = (t107 * t75 + t251 * t76) * pkin(4);
t278 = t111 * t237;
t250 = qJ(3) * t111;
t180 = -pkin(2) * t114 - t250;
t166 = -pkin(1) + t180;
t154 = pkin(9) * t114 - t166;
t142 = qJD(4) * t154;
t215 = t114 * t237;
t240 = qJD(2) * t111;
t157 = t113 * t240 + t215;
t127 = qJ(5) * t157 + t110 * t142;
t239 = qJD(2) * t114;
t117 = t251 * t127 + (t113 * t142 - t110 * (t271 * t240 - t238) - t270 * t278 + (t110 * qJ(3) + t113 * t270 + pkin(4)) * t239) * t107;
t197 = t111 * t218;
t143 = t114 * t225 + t197;
t138 = qJD(2) * t143;
t230 = -t110 * t138 + t280 * t113;
t245 = t108 * t114;
t93 = t113 * t245;
t99 = t252 * t110;
t263 = t99 + t93;
t205 = t252 * t113;
t97 = qJD(4) * t205;
t145 = -t97 * qJ(5) - t263 * qJD(5) - t230;
t214 = t108 * t240;
t59 = t263 * qJD(4) - t110 * t214;
t219 = t110 * t245;
t77 = t205 - t219;
t170 = t59 * qJ(5) - t77 * qJD(5);
t203 = qJD(2) * t252;
t193 = t111 * t203;
t120 = t107 * (pkin(1) * t113 * t193 - t133 * t237 + t170) + t251 * t145;
t206 = t252 * qJ(3);
t60 = t206 + t143;
t42 = t263 * pkin(4) + t60;
t45 = t107 * t77 + t251 * t263;
t46 = -t107 * t263 + t251 * t77;
t126 = t45 * pkin(5) - t46 * pkin(10) + t42;
t277 = -(pkin(10) * t239 + t117) * t108 - t120 - qJD(6) * t126;
t109 = sin(qJ(6));
t105 = t109 ^ 2;
t112 = cos(qJ(6));
t106 = t112 ^ 2;
t242 = t105 - t106;
t200 = qJD(6) * t242;
t246 = t108 * t111;
t229 = pkin(8) * t246;
t152 = -t196 + t229;
t73 = t152 * qJD(2);
t146 = t157 * t108;
t139 = t97 - t146;
t208 = t251 * t59;
t130 = t107 * t139 + t208;
t83 = -t107 * t110 + t204;
t276 = -t130 * t83 - t46 * t76;
t151 = t108 * (-t271 * t114 - pkin(1) - t250);
t141 = qJD(4) * t151;
t18 = t110 * t141 + t230;
t19 = (t138 - t141) * t113 + t280 * t110;
t52 = t113 * t129;
t32 = -t110 * t151 + t52;
t33 = t110 * t129 + t113 * t151;
t134 = -qJD(4) * (t32 * t110 - t33 * t113) - t18 * t110 + t19 * t113;
t243 = qJ(5) + t271;
t153 = t113 * qJD(5) - t243 * t237;
t202 = t243 * t113;
t72 = -qJD(4) * t202 - t110 * qJD(5);
t131 = -t107 * t153 + t251 * t72;
t102 = pkin(4) * t110 + qJ(3);
t155 = pkin(5) * t82 - pkin(10) * t83 + t102;
t144 = t112 * t155;
t236 = qJD(4) * t113;
t94 = pkin(4) * t236 + qJD(3);
t149 = -pkin(5) * t75 + pkin(10) * t76 + t94;
t234 = qJD(6) * t109;
t84 = t243 * t110;
t57 = -t107 * t202 - t251 * t84;
t15 = -qJD(6) * t144 - t109 * t149 - t112 * t131 + t234 * t57;
t29 = t109 * t155 + t112 * t57;
t16 = -qJD(6) * t29 - t109 * t131 + t112 * t149;
t28 = -t109 * t57 + t144;
t174 = t109 * t28 - t112 * t29;
t275 = qJD(6) * t174 + t109 * t15 - t112 * t16;
t147 = t107 * t157;
t161 = -t107 * t97 - t208;
t192 = t114 * t203;
t201 = t252 * qJD(3);
t31 = -t107 * t59 + t251 * t139;
t119 = pkin(1) * t192 + t201 + t97 * pkin(4) - t161 * pkin(10) + t31 * pkin(5) + (-pkin(4) * t157 - pkin(10) * t147 - t270 * t240) * t108;
t25 = -t77 * qJ(5) + t52 + (t111 * pkin(4) + t110 * t154) * t108;
t27 = -t263 * qJ(5) + t33;
t14 = t107 * t25 + t251 * t27;
t12 = pkin(10) * t246 + t14;
t1 = -t109 * t119 + t277 * t112 + t12 * t234;
t7 = -t109 * t12 + t112 * t126;
t8 = t109 * t126 + t112 * t12;
t188 = t109 * t7 - t112 * t8;
t233 = qJD(6) * t112;
t2 = t277 * t109 + t112 * t119 - t12 * t233;
t274 = qJD(6) * t188 + t1 * t109 - t112 * t2;
t81 = t83 ^ 2;
t273 = 0.2e1 * t108;
t272 = 0.2e1 * qJD(3);
t125 = t107 * t127;
t95 = t108 * t239;
t199 = pkin(4) * t95;
t118 = t170 + t19 + t199;
t136 = -t107 * t145 + t251 * t118;
t4 = (-pkin(5) * t239 + t125) * t108 - t136;
t269 = t112 * t4;
t41 = t107 * t72 + t251 * t153;
t56 = -t107 * t84 + t251 * t202;
t267 = t56 * t41;
t266 = t82 * t75;
t265 = t83 * t76;
t38 = t109 * t246 + t112 * t46;
t21 = qJD(6) * t38 - t109 * t130 - t112 * t95;
t220 = t112 * t246;
t37 = t109 * t46 - t220;
t264 = -t109 * t21 - t37 * t233;
t217 = t251 * pkin(4);
t101 = -t217 - pkin(5);
t262 = t101 * t76;
t20 = -qJD(6) * t220 - t109 * t95 + t112 * t130 + t234 * t46;
t261 = t109 * t20;
t260 = t109 * t37;
t259 = t109 * t38;
t79 = pkin(8) * t245 + t197;
t74 = t79 * qJD(2);
t258 = t111 * t74;
t257 = t112 * t21;
t256 = t112 * t37;
t255 = t112 * t38;
t254 = t112 * t76;
t253 = t113 * t59;
t100 = pkin(4) * t107 + pkin(10);
t248 = t100 * t109;
t247 = t100 * t112;
t244 = t108 * t271;
t241 = t105 + t106;
t235 = qJD(4) * t271;
t232 = 0.2e1 * t45 * t31;
t231 = -0.2e1 * t266;
t228 = 0.2e1 * qJD(6) * t101;
t227 = t109 * t254;
t226 = t82 ^ 2 + t81;
t224 = t82 * t234;
t223 = t83 * t234;
t222 = t82 * t233;
t221 = t83 * t233;
t104 = t108 ^ 2;
t216 = t104 * t239;
t213 = t109 * t233;
t212 = t110 * t239;
t211 = t110 * t236;
t210 = t74 * t252;
t209 = t263 * t113;
t207 = t241 * t75;
t198 = t81 * t213;
t195 = t111 * t216;
t13 = -t107 * t27 + t251 * t25;
t11 = -pkin(5) * t246 - t13;
t194 = -t11 * t76 + t4 * t83;
t191 = -t1 * t112 - t109 * t2;
t189 = t109 * t8 + t112 * t7;
t187 = -t20 * t83 - t38 * t76;
t186 = t21 * t83 - t37 * t76;
t185 = t31 * t82 - t45 * t75;
t184 = t83 * t31 - t76 * t45;
t183 = t41 * t83 - t56 * t76;
t182 = t75 * t83 + t76 * t82;
t181 = t265 + t266;
t179 = t100 * t75 - t262;
t178 = t100 * t82 - t101 * t83;
t176 = -t257 - t261;
t175 = t109 * t29 + t112 * t28;
t173 = t255 + t260;
t169 = t108 * t193;
t168 = t108 * t192;
t167 = t113 * t214 - t97;
t165 = -t112 * t20 - t234 * t38;
t164 = t109 * t31 + t233 * t45;
t48 = -t109 * t75 + t222;
t163 = t109 * t76 - t221;
t162 = -t223 - t254;
t159 = 0.2e1 * t181;
t158 = t113 * t239 - t278;
t156 = (t111 * t75 - t239 * t82) * t108;
t124 = t127 * t108;
t5 = -t107 * t124 + t136;
t6 = t108 * t117 + t120;
t148 = t13 * t76 + t14 * t75 - t5 * t83 - t6 * t82;
t140 = t111 * t146;
t137 = -qJD(6) * t189 + t191;
t135 = -qJD(6) * t175 - t109 * t16 - t112 * t15;
t122 = t167 * t110 + t253 + (-t209 + (t77 + t219) * t110) * qJD(4);
t121 = t131 * t82 - t57 * t75 - t183;
t103 = qJ(3) * t272;
t86 = -0.2e1 * t195;
t85 = 0.2e1 * t195;
t71 = 0.2e1 * (-t111 ^ 2 + t114 ^ 2) * t104 * qJD(2);
t68 = t150 + t229;
t67 = t166 * t108;
t66 = -t206 - t79;
t65 = t158 * t108;
t64 = (-t111 * t236 - t212) * t108;
t62 = -t201 + t73;
t61 = (-t238 + (pkin(2) * t111 - t249) * qJD(2)) * t108;
t54 = t201 + (-t270 * t246 + t196) * qJD(2);
t47 = t112 * t75 + t224;
t43 = (-t111 * t76 + t239 * t83) * t108;
t35 = t200 * t83 + t227;
t34 = t201 - (t108 * t215 - t97) * pkin(4) + (t196 + (-pkin(4) * t113 - t270) * t246) * qJD(2);
t23 = t112 * t31 - t234 * t45;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t71, 0.2e1 * t168, t86, -0.2e1 * t169, 0, -0.2e1 * pkin(1) * t104 * t240 - 0.2e1 * t210, -0.2e1 * pkin(1) * t216 + 0.2e1 * t73 * t252 (t258 - t114 * t73 + (-t111 * t79 + t114 * t152) * qJD(2)) * t273, 0.2e1 * t152 * t74 - 0.2e1 * t73 * t79, 0, -0.2e1 * t168, 0.2e1 * t169, t85, t71, t86 (t258 - t114 * t62 + (t111 * t66 + t114 * t68) * qJD(2)) * t273, 0.2e1 * t210 + 0.2e1 * (t114 * t61 - t240 * t67) * t108, -0.2e1 * t62 * t252 + 0.2e1 * (-t111 * t61 - t239 * t67) * t108, 0.2e1 * t61 * t67 + 0.2e1 * t62 * t66 + 0.2e1 * t68 * t74, -0.2e1 * t77 * t59, -0.2e1 * t77 * t139 + 0.2e1 * t59 * t263 (-t111 * t59 + t239 * t77) * t273, 0.2e1 * t263 * t139 (-t97 * t111 - t263 * t239 + t140) * t273, t85, 0.2e1 * t54 * t263 + 0.2e1 * t60 * t97 + 0.2e1 * (t19 * t111 - t157 * t60 + t239 * t32) * t108, 0.2e1 * t54 * t77 - 0.2e1 * t59 * t60 + 0.2e1 * (t111 * t18 - t239 * t33) * t108, -0.2e1 * t33 * t139 + 0.2e1 * t18 * t263 - 0.2e1 * t19 * t77 + 0.2e1 * t32 * t59, -0.2e1 * t18 * t33 + 0.2e1 * t19 * t32 + 0.2e1 * t54 * t60, -0.2e1 * t46 * t130, 0.2e1 * t130 * t45 - 0.2e1 * t46 * t31 (t107 * t140 + t111 * t161 + t239 * t46) * t273, t232 (-t111 * t31 - t239 * t45) * t273, t85, 0.2e1 * t31 * t42 + 0.2e1 * t34 * t45 + 0.2e1 * (t111 * t5 + t13 * t239) * t108, 0.2e1 * t34 * t46 + 0.2e1 * t42 * t161 + 0.2e1 * (-t6 * t111 - t14 * t239 + t147 * t42) * t108, 0.2e1 * t13 * t130 - 0.2e1 * t14 * t31 - 0.2e1 * t6 * t45 - 0.2e1 * t5 * t46, 0.2e1 * t13 * t5 + 0.2e1 * t14 * t6 + 0.2e1 * t34 * t42, -0.2e1 * t38 * t20, 0.2e1 * t20 * t37 - 0.2e1 * t21 * t38, -0.2e1 * t20 * t45 + 0.2e1 * t31 * t38, 0.2e1 * t37 * t21, -0.2e1 * t21 * t45 - 0.2e1 * t31 * t37, t232, 0.2e1 * t11 * t21 + 0.2e1 * t2 * t45 + 0.2e1 * t31 * t7 + 0.2e1 * t37 * t4, 0.2e1 * t1 * t45 - 0.2e1 * t11 * t20 - 0.2e1 * t31 * t8 + 0.2e1 * t38 * t4, 0.2e1 * t1 * t37 - 0.2e1 * t2 * t38 + 0.2e1 * t20 * t7 - 0.2e1 * t21 * t8, -0.2e1 * t1 * t8 + 0.2e1 * t11 * t4 + 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, -t214, 0, -t74, t73, 0, 0, 0, -t95, t214, 0, 0, 0 (qJD(2) * t180 + qJD(3) * t114) * t108, t74, 0.2e1 * t201 - t73, -pkin(2) * t74 - qJ(3) * t62 - qJD(3) * t66, -t237 * t77 - t253, t113 * t167 + t59 * t110 + (-t77 * t113 + (t99 + 0.2e1 * t93) * t110) * qJD(4), t65, qJD(4) * t209 + t110 * t139, t64, 0, qJD(3) * t263 + qJ(3) * t97 + t54 * t110 + t60 * t236 + (-qJ(3) * t157 - t158 * t271) * t108, t212 * t244 - qJ(3) * t59 + qJD(3) * t77 + t54 * t113 + (t111 * t113 * t244 - t110 * t60) * qJD(4), -t122 * t271 - t134, t54 * qJ(3) + t60 * qJD(3) - t134 * t271, t276, t130 * t82 + t46 * t75 - t184, t43, t185, t156, 0, t102 * t31 + t34 * t82 - t42 * t75 + t45 * t94 + (-t111 * t41 - t239 * t56) * t108, t94 * t46 + t102 * t161 + t34 * t83 - t42 * t76 + (t102 * t147 - t111 * t131 - t239 * t57) * t108, -t130 * t56 - t131 * t45 - t57 * t31 + t41 * t46 + t148, t34 * t102 - t13 * t41 + t131 * t14 + t42 * t94 - t5 * t56 + t6 * t57, t112 * t187 - t223 * t38 (t256 + t259) * t76 + (t261 - t257 + (-t255 + t260) * qJD(6)) * t83, t112 * t184 - t20 * t82 - t223 * t45 - t38 * t75, t109 * t186 + t221 * t37, -t109 * t184 - t21 * t82 - t221 * t45 + t37 * t75, t185, t109 * t194 + t11 * t221 + t16 * t45 + t2 * t82 + t21 * t56 + t28 * t31 + t37 * t41 - t7 * t75, t1 * t82 + t11 * t162 + t15 * t45 - t20 * t56 + t83 * t269 - t29 * t31 + t38 * t41 + t75 * t8, t15 * t37 - t16 * t38 + t189 * t76 + t20 * t28 - t21 * t29 + t274 * t83, -t1 * t29 + t11 * t41 - t15 * t8 + t16 * t7 + t2 * t28 + t4 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, t103, -0.2e1 * t211, 0.2e1 * (t110 ^ 2 - t113 ^ 2) * qJD(4), 0, 0.2e1 * t211, 0, 0, 0.2e1 * qJ(3) * t236 + 0.2e1 * t110 * qJD(3), -0.2e1 * qJ(3) * t237 + 0.2e1 * qJD(3) * t113, 0, t103, -0.2e1 * t265, 0.2e1 * t182, 0, t231, 0, 0, -0.2e1 * t102 * t75 + 0.2e1 * t82 * t94, -0.2e1 * t102 * t76 + 0.2e1 * t83 * t94, -0.2e1 * t121, 0.2e1 * t102 * t94 + 0.2e1 * t131 * t57 + 0.2e1 * t267, -0.2e1 * t106 * t265 - 0.2e1 * t198, 0.2e1 * t200 * t81 + 0.4e1 * t227 * t83, -0.2e1 * t112 * t182 - 0.2e1 * t223 * t82, -0.2e1 * t105 * t265 + 0.2e1 * t198, 0.2e1 * t109 * t182 - 0.2e1 * t221 * t82, t231, 0.2e1 * t109 * t183 + 0.2e1 * t16 * t82 + 0.2e1 * t221 * t56 - 0.2e1 * t28 * t75, 0.2e1 * t112 * t183 + 0.2e1 * t15 * t82 - 0.2e1 * t223 * t56 + 0.2e1 * t29 * t75, 0.2e1 * t175 * t76 + 0.2e1 * t275 * t83, -0.2e1 * t15 * t29 + 0.2e1 * t16 * t28 + 0.2e1 * t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, 0, t74, 0, 0, 0, 0, 0, 0, t65, t64, t122, t134, 0, 0, 0, 0, 0, 0, t43, t156, -t185 - t276, -t148, 0, 0, 0, 0, 0, 0, -t109 * t185 - t222 * t45 - t186, -t112 * t185 + t224 * t45 - t187 (t256 - t259) * t75 + (qJD(6) * t173 + t176) * t82, t137 * t82 + t188 * t75 - t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t121, 0, 0, 0, 0, 0, 0, t109 * t159 - t226 * t233, t112 * t159 + t226 * t234, 0, t135 * t82 + t174 * t75 - t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t181, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t207 * t82 - 0.2e1 * t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, -t139, t95, t19, t18, 0, 0, 0, 0, -t130, 0, -t31, t95 (t217 * t239 - t125) * t108 + t136, -t251 * (t124 + t145) + (-t199 - t118) * t107 (-t107 * t31 + t251 * t130) * pkin(4) (t107 * t6 + t251 * t5) * pkin(4), t233 * t38 - t261, t165 + t264, t164, t234 * t37 - t257, t23, 0, -t31 * t248 + t101 * t21 - t269 + (t109 * t11 - t45 * t247) * qJD(6), -t31 * t247 - t101 * t20 + t109 * t4 + (t11 * t112 + t248 * t45) * qJD(6), t176 * t100 + (t100 * t173 - t189) * qJD(6) + t191, t100 * t137 + t101 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, 0, -t236, 0, t110 * t235, t113 * t235, 0, 0, 0, 0, -t76, 0, t75, 0, -t41, -t131, t279 (t131 * t107 - t41 * t251) * pkin(4), -t35, -0.4e1 * t83 * t213 + t242 * t76, t48, t35, -t47, 0, -t112 * t41 + t179 * t109 + (t109 * t56 - t112 * t178) * qJD(6), t109 * t41 + t179 * t112 + (t109 * t178 + t112 * t56) * qJD(6), t135, t100 * t135 + t101 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t237, -t236, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t75, 0, -t279, 0, 0, 0, 0, 0, 0, t162, t163, -t207, -t100 * t207 + t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t213, -0.2e1 * t200, 0, -0.2e1 * t213, 0, 0, t109 * t228, t112 * t228, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t130, 0, t34, 0, 0, 0, 0, 0, 0, t23, -t164, -t165 + t264, -t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t76, 0, t94, 0, 0, 0, 0, 0, 0, -t47, -t48, t241 * t76, -t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, -t21, t31, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, 0, t163, -t75, t16, t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233, 0, -t234, 0, -t100 * t233, t100 * t234, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, -t233, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;

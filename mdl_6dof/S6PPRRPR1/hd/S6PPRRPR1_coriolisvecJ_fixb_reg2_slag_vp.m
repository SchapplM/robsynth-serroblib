% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:38
% EndTime: 2019-03-08 18:47:46
% DurationCPUTime: 4.80s
% Computational Cost: add. (6056->401), mult. (16427->620), div. (0->0), fcn. (14303->14), ass. (0->210)
t163 = sin(qJ(4));
t166 = cos(qJ(4));
t193 = pkin(4) * t163 - qJ(5) * t166;
t112 = qJD(4) * t193 - qJD(5) * t163;
t157 = sin(pkin(6));
t155 = sin(pkin(12));
t164 = sin(qJ(3));
t167 = cos(qJ(3));
t159 = cos(pkin(12));
t160 = cos(pkin(7));
t234 = t159 * t160;
t183 = t155 * t167 + t164 * t234;
t175 = t183 * t157;
t161 = cos(pkin(6));
t144 = qJD(1) * t161 + qJD(2);
t156 = sin(pkin(7));
t239 = t156 * t164;
t215 = t144 * t239;
t80 = qJD(1) * t175 + t215;
t274 = t112 - t80;
t154 = sin(pkin(13));
t220 = qJD(4) * t163;
t216 = pkin(9) * t220;
t140 = t154 * t216;
t158 = cos(pkin(13));
t241 = t154 * t166;
t227 = qJD(1) * t157;
t210 = t159 * t227;
t240 = t155 * t164;
t78 = t167 * (t144 * t156 + t160 * t210) - t227 * t240;
t255 = -t158 * t274 - t241 * t78 - t140;
t235 = t158 * t166;
t273 = t154 * t274 - t235 * t78;
t221 = qJD(4) * t158;
t225 = qJD(3) * t163;
t124 = -t154 * t225 + t221;
t206 = t158 * t225;
t222 = qJD(4) * t154;
t125 = t206 + t222;
t162 = sin(qJ(6));
t165 = cos(qJ(6));
t84 = -t165 * t124 + t125 * t162;
t272 = t84 ^ 2;
t87 = t124 * t162 + t125 * t165;
t271 = t87 ^ 2;
t187 = pkin(5) * t163 - pkin(10) * t235;
t180 = t187 * qJD(4);
t270 = -t180 + t255;
t236 = t158 * t163;
t269 = -(-pkin(9) * t236 - pkin(10) * t241) * qJD(4) - t273;
t223 = qJD(3) * t166;
t146 = -qJD(6) + t223;
t268 = t146 * t84;
t238 = t156 * t167;
t265 = (t167 * t234 - t240) * t157;
t267 = t161 * t238 + t265;
t266 = qJD(1) * t265 + t144 * t238;
t218 = qJD(6) * t165;
t242 = t154 * t162;
t264 = -qJD(6) * t242 + t158 * t218;
t217 = qJD(3) * qJD(4);
t203 = t166 * t217;
t195 = t154 * t203;
t219 = qJD(4) * t166;
t233 = t165 * t158;
t196 = t219 * t233;
t58 = t162 * (qJD(6) * t125 + t195) - qJD(3) * t196 - t124 * t218;
t101 = t144 * t160 - t156 * t210;
t72 = t266 * qJD(3);
t246 = t101 * t219 + t166 * t72;
t77 = qJD(3) * pkin(9) + t80;
t68 = t163 * t77;
t30 = (qJD(5) - t68) * qJD(4) + t246;
t57 = (t80 + t112) * qJD(3);
t14 = -t154 * t30 + t158 * t57;
t10 = qJD(3) * t180 + t14;
t15 = t154 * t57 + t158 * t30;
t11 = -pkin(10) * t195 + t15;
t245 = t101 * t163;
t50 = t166 * t77 + t245;
t48 = qJD(4) * qJ(5) + t50;
t133 = -pkin(4) * t166 - qJ(5) * t163 - pkin(3);
t64 = qJD(3) * t133 - t78;
t26 = -t154 * t48 + t158 * t64;
t22 = -pkin(5) * t223 - pkin(10) * t125 + t26;
t27 = t154 * t64 + t158 * t48;
t23 = pkin(10) * t124 + t27;
t191 = t162 * t23 - t165 * t22;
t1 = -qJD(6) * t191 + t10 * t162 + t11 * t165;
t127 = -t233 + t242;
t129 = t193 * qJD(3);
t49 = t101 * t166 - t68;
t39 = t158 * t129 - t154 * t49;
t33 = qJD(3) * t187 + t39;
t209 = t154 * t223;
t40 = t154 * t129 + t158 * t49;
t34 = -pkin(10) * t209 + t40;
t258 = pkin(10) + qJ(5);
t136 = t258 * t154;
t137 = t258 * t158;
t98 = -t136 * t165 - t137 * t162;
t263 = -qJD(5) * t127 + qJD(6) * t98 - t162 * t33 - t165 * t34;
t128 = t154 * t165 + t158 * t162;
t99 = -t136 * t162 + t137 * t165;
t262 = -qJD(5) * t128 - qJD(6) * t99 + t162 * t34 - t165 * t33;
t113 = -t156 * t157 * t159 + t160 * t161;
t90 = t161 * t239 + t175;
t190 = t113 * t166 - t163 * t90;
t32 = qJD(4) * t50 + t163 * t72;
t261 = t32 * t190;
t73 = qJD(3) * t80;
t260 = t73 * t267;
t259 = t87 * t84;
t123 = t158 * t133;
t88 = -pkin(10) * t236 + t123 + (-pkin(9) * t154 - pkin(5)) * t166;
t103 = pkin(9) * t235 + t154 * t133;
t97 = -pkin(10) * t154 * t163 + t103;
t52 = t162 * t88 + t165 * t97;
t257 = qJD(6) * t52 - t269 * t162 + t270 * t165;
t51 = -t162 * t97 + t165 * t88;
t256 = -qJD(6) * t51 + t270 * t162 + t269 * t165;
t199 = t158 * t216;
t254 = -t199 + t273;
t253 = qJD(3) * pkin(3);
t117 = -t166 * t160 + t163 * t239;
t252 = t117 * t32;
t251 = t154 * t32;
t250 = t158 * t32;
t249 = t163 * t32;
t248 = t163 * t78;
t247 = t167 * t73;
t244 = t124 * t158;
t153 = t166 ^ 2;
t169 = qJD(3) ^ 2;
t243 = t153 * t169;
t237 = t156 * t169;
t116 = t128 * qJD(6);
t181 = t128 * t166;
t232 = qJD(3) * t181 - t116;
t231 = -t127 * t223 - t264;
t152 = t163 ^ 2;
t229 = t152 - 0.2e1 * t153;
t228 = t152 - t153;
t226 = qJD(3) * t154;
t224 = qJD(3) * t164;
t212 = t164 * t237;
t211 = pkin(5) * t154 + pkin(9);
t208 = t156 * t224;
t207 = qJD(3) * t238;
t202 = t163 * t217;
t201 = -qJD(4) * pkin(4) + qJD(5);
t200 = -t124 + t221;
t198 = t163 * t207;
t197 = t166 * t207;
t194 = t166 * t202;
t6 = t162 * t22 + t165 * t23;
t63 = t113 * t163 + t166 * t90;
t37 = -t154 * t63 - t158 * t267;
t38 = -t154 * t267 + t158 * t63;
t12 = -t162 * t38 + t165 * t37;
t13 = t162 * t37 + t165 * t38;
t118 = t160 * t163 + t166 * t239;
t93 = -t118 * t154 - t158 * t238;
t94 = t118 * t158 - t154 * t238;
t55 = -t162 * t94 + t165 * t93;
t56 = t162 * t93 + t165 * t94;
t189 = qJD(3) * t200;
t188 = qJD(3) * (-t125 + t222);
t168 = qJD(4) ^ 2;
t186 = pkin(9) * t168;
t76 = -t78 - t253;
t185 = qJD(4) * (t76 + t78 - t253);
t182 = t166 * t188;
t46 = t201 - t49;
t174 = qJD(4) * t181;
t173 = -qJ(5) * t220 + (t201 - t46) * t166;
t2 = -qJD(6) * t6 + t165 * t10 - t11 * t162;
t28 = pkin(5) * t195 + t32;
t31 = -t220 * t77 + t246;
t170 = t249 + t166 * t31 + (-t163 * t50 - t166 * t49) * qJD(4);
t59 = qJD(3) * t174 + qJD(6) * t87;
t151 = t158 ^ 2;
t150 = t154 ^ 2;
t149 = -pkin(5) * t158 - pkin(4);
t143 = t166 * t169 * t163;
t139 = -0.2e1 * t194;
t130 = t211 * t163;
t121 = t211 * t219;
t110 = t127 * t163;
t109 = t128 * t163;
t102 = -pkin(9) * t241 + t123;
t96 = qJD(4) * t118 + t198;
t95 = -qJD(4) * t117 + t197;
t82 = t90 * qJD(3);
t81 = t267 * qJD(3);
t75 = t163 * t264 + t174;
t74 = t116 * t163 + t219 * t242 - t196;
t71 = t154 * t208 + t158 * t95;
t70 = -t154 * t95 + t158 * t208;
t44 = t245 + (pkin(5) * t226 + t77) * t166;
t41 = -pkin(5) * t124 + t46;
t36 = qJD(4) * t190 + t166 * t81;
t35 = qJD(4) * t63 + t163 * t81;
t25 = t154 * t82 + t158 * t36;
t24 = -t154 * t36 + t158 * t82;
t19 = -qJD(6) * t56 - t162 * t71 + t165 * t70;
t18 = qJD(6) * t55 + t162 * t70 + t165 * t71;
t4 = -qJD(6) * t13 - t162 * t25 + t165 * t24;
t3 = qJD(6) * t12 + t162 * t24 + t165 * t25;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82 * qJD(3), -t81 * qJD(3), 0, t72 * t90 - t78 * t82 + t80 * t81 - t260, 0, 0, 0, 0, 0, 0, -qJD(4) * t35 + (-t166 * t82 - t220 * t267) * qJD(3), -qJD(4) * t36 + (t163 * t82 - t219 * t267) * qJD(3) (t163 * t35 + t166 * t36 + (-t163 * t63 - t166 * t190) * qJD(4)) * qJD(3), t31 * t63 - t35 * t49 + t36 * t50 + t76 * t82 - t260 - t261, 0, 0, 0, 0, 0, 0, -t124 * t35 + (-t166 * t24 + (t163 * t37 - t190 * t241) * qJD(4)) * qJD(3), t125 * t35 + (t166 * t25 + (-t163 * t38 - t190 * t235) * qJD(4)) * qJD(3), t124 * t25 - t125 * t24 + (-t154 * t38 - t158 * t37) * t203, t14 * t37 + t15 * t38 + t24 * t26 + t25 * t27 + t35 * t46 - t261, 0, 0, 0, 0, 0, 0, t12 * t202 - t146 * t4 - t190 * t59 + t35 * t84, -t13 * t202 + t146 * t3 + t190 * t58 + t35 * t87, t12 * t58 - t13 * t59 - t3 * t84 - t4 * t87, t1 * t13 + t12 * t2 - t190 * t28 - t191 * t4 + t3 * t6 + t35 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t212, -t167 * t237, 0 (t164 * t72 - t247 + (-t164 * t78 + t167 * t80) * qJD(3)) * t156, 0, 0, 0, 0, 0, 0, -t166 * t212 + (-t96 - t198) * qJD(4), t163 * t212 + (-t95 - t197) * qJD(4) (t163 * t96 + t166 * t95 + (t117 * t166 - t118 * t163) * qJD(4)) * qJD(3), t252 + t118 * t31 - t49 * t96 + t50 * t95 + (t224 * t76 - t247) * t156, 0, 0, 0, 0, 0, 0, -t124 * t96 + (-t166 * t70 + (t117 * t241 + t163 * t93) * qJD(4)) * qJD(3), t125 * t96 + (t166 * t71 + (t117 * t235 - t163 * t94) * qJD(4)) * qJD(3), t124 * t71 - t125 * t70 + (-t154 * t94 - t158 * t93) * t203, t14 * t93 + t15 * t94 + t26 * t70 + t27 * t71 + t46 * t96 + t252, 0, 0, 0, 0, 0, 0, t117 * t59 - t146 * t19 + t202 * t55 + t84 * t96, -t117 * t58 + t146 * t18 - t202 * t56 + t87 * t96, -t18 * t84 - t19 * t87 + t55 * t58 - t56 * t59, t1 * t56 + t117 * t28 + t18 * t6 - t19 * t191 + t2 * t55 + t41 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t183 * t227 - t215 + t80) * qJD(3) (-t266 + t78) * qJD(3), 0, 0, 0.2e1 * t194, -0.2e1 * t228 * t217, t168 * t166, t139, -t168 * t163, 0, t163 * t185 - t166 * t186, t163 * t186 + t166 * t185 (-t152 - t153) * t78 * qJD(3) + t170, -pkin(3) * t73 - t76 * t80 + (t163 * t49 - t166 * t50) * t78 + t170 * pkin(9) (t125 * t158 + t151 * t225) * t219 (t244 + (-t125 - 0.2e1 * t206) * t154) * t219 (qJD(3) * t158 * t229 + t125 * t163) * qJD(4) (-t124 * t154 + t150 * t225) * t219 (t124 * t163 - t226 * t229) * qJD(4), t139 (t124 * t78 + t251 + (qJD(3) * t102 + t26) * qJD(4)) * t163 + (-t14 + (-pkin(9) * t124 + t154 * t46) * qJD(4) + (t140 + t255) * qJD(3)) * t166 (-t125 * t78 + t250 + (-qJD(3) * t103 - t27) * qJD(4)) * t163 + (t15 + (pkin(9) * t125 + t158 * t46) * qJD(4) + (t199 + t254) * qJD(3)) * t166 (-t14 * t158 - t15 * t154) * t163 + t255 * t125 + t254 * t124 + (-t154 * t27 - t158 * t26 + (-t102 * t158 - t103 * t154) * qJD(3)) * t219, -t46 * t248 + t102 * t14 + t103 * t15 + t254 * t27 - t255 * t26 + (t219 * t46 + t249) * pkin(9), t110 * t58 - t74 * t87, t109 * t58 + t110 * t59 + t74 * t84 - t75 * t87, t146 * t74 + t166 * t58 + (-qJD(3) * t110 + t87) * t220, t109 * t59 + t75 * t84, t146 * t75 + t166 * t59 + (-qJD(3) * t109 - t84) * t220 (-t146 - t223) * t220, t109 * t28 + t121 * t84 + t130 * t59 - t166 * t2 + t41 * t75 + t257 * t146 + (-t78 * t84 + (qJD(3) * t51 - t191) * qJD(4)) * t163, t1 * t166 - t110 * t28 + t121 * t87 - t130 * t58 - t41 * t74 - t256 * t146 + (-t78 * t87 + (-qJD(3) * t52 - t6) * qJD(4)) * t163, -t1 * t109 + t110 * t2 - t191 * t74 + t256 * t84 + t257 * t87 + t51 * t58 - t52 * t59 - t6 * t75, t1 * t52 + t130 * t28 + t2 * t51 - t256 * t6 + t257 * t191 + (t121 - t248) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t228 * t169, 0, t143, 0, 0 (-qJD(3) * t76 - t72) * t163, -t76 * t223 + (t49 + t68) * qJD(4) - t246, 0, 0, t158 * t182 (-t244 + t125 * t154 + (-t150 + t151) * qJD(4)) * t223, t158 * t243 + t163 * t188, -t200 * t209, -t154 * t243 + t163 * t189, t143, t124 * t50 - t250 + (t154 * t173 - t163 * t26 + t166 * t39) * qJD(3), -t125 * t50 + t251 + (t158 * t173 + t163 * t27 - t166 * t40) * qJD(3), -t124 * t40 + t125 * t39 + (qJD(5) * t124 + t223 * t26 + t15) * t158 + (qJD(5) * t125 + t223 * t27 - t14) * t154, -pkin(4) * t32 - t26 * t39 - t27 * t40 - t46 * t50 + (-t154 * t26 + t158 * t27) * qJD(5) + (-t14 * t154 + t15 * t158) * qJ(5), -t128 * t58 - t231 * t87, t127 * t58 - t128 * t59 + t231 * t84 + t232 * t87, t231 * t146 + (qJD(4) * t128 - t87) * t225, t127 * t59 - t232 * t84, -t232 * t146 + (-qJD(4) * t127 + t84) * t225, t146 * t225, t127 * t28 + t149 * t59 - t44 * t84 - t232 * t41 - t262 * t146 + (qJD(4) * t98 + t191) * t225, t128 * t28 - t149 * t58 - t44 * t87 - t231 * t41 + t263 * t146 + (-qJD(4) * t99 + t6) * t225, -t1 * t127 - t128 * t2 - t191 * t231 + t232 * t6 - t262 * t87 - t263 * t84 + t58 * t98 - t59 * t99, t1 * t99 + t149 * t28 - t191 * t262 + t2 * t98 + t263 * t6 - t41 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, t166 * t189, -t124 ^ 2 - t125 ^ 2, -t124 * t27 + t125 * t26 + t32, 0, 0, 0, 0, 0, 0, -t146 * t87 + t59, -t58 + t268, -t271 - t272, -t191 * t87 + t6 * t84 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, t271 - t272, -t58 - t268, -t259, -t128 * t203 + (-qJD(6) - t146) * t87, t202, -t146 * t6 - t41 * t87 + t2, t146 * t191 + t41 * t84 - t1, 0, 0;];
tauc_reg  = t5;

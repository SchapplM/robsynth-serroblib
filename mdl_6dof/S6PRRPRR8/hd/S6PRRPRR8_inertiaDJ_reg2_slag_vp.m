% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_inertiaDJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:38
% EndTime: 2019-03-08 22:40:52
% DurationCPUTime: 5.35s
% Computational Cost: add. (4478->409), mult. (12667->730), div. (0->0), fcn. (12361->12), ass. (0->222)
t100 = cos(qJ(5));
t101 = cos(qJ(3));
t102 = cos(qJ(2));
t218 = cos(pkin(7));
t171 = t102 * t218;
t93 = sin(pkin(6));
t157 = t93 * t171;
t97 = sin(qJ(3));
t98 = sin(qJ(2));
t198 = t93 * t97 * t98;
t92 = sin(pkin(7));
t221 = t101 * t92;
t94 = cos(pkin(6));
t36 = -t101 * t157 - t221 * t94 + t198;
t57 = -t93 * t102 * t92 + t218 * t94;
t96 = sin(qJ(5));
t143 = t36 * t100 - t57 * t96;
t256 = qJD(5) * t143;
t215 = qJD(3) * t97;
t189 = t92 * t215;
t168 = qJD(5) * t218;
t184 = qJD(5) * t221;
t228 = -t100 * t184 - t96 * t168;
t129 = t189 * t96 + t228;
t255 = t129 * t100;
t160 = t218 * t97 * pkin(2);
t246 = pkin(4) + pkin(9);
t193 = t246 * t92;
t123 = t101 * t193 + t160;
t247 = pkin(3) + pkin(10);
t156 = -t247 * t101 - pkin(2);
t223 = qJ(4) * t97;
t125 = t92 * (t156 - t223);
t252 = qJD(3) * t123 - qJD(5) * t125;
t172 = t101 * t218;
t158 = pkin(2) * t172;
t124 = -pkin(3) * t218 - t158;
t115 = -pkin(10) * t218 + t124;
t109 = t193 * t97 + t115;
t214 = qJD(4) * t97;
t217 = qJ(4) * t101;
t253 = -qJD(5) * t109 - t92 * (-t214 + (t247 * t97 - t217) * qJD(3));
t10 = t253 * t100 - t252 * t96;
t169 = t218 * qJ(4);
t43 = t169 + t123;
t59 = t100 * t221 + t218 * t96;
t60 = t100 * t218 - t221 * t96;
t108 = t59 * pkin(5) - t60 * pkin(11) + t43;
t206 = qJD(3) * t101;
t84 = t92 * t206;
t254 = -pkin(11) * t84 - qJD(6) * t108 + t10;
t95 = sin(qJ(6));
t88 = t95 ^ 2;
t99 = cos(qJ(6));
t90 = t99 ^ 2;
t227 = t88 - t90;
t173 = qJD(6) * t227;
t230 = t92 * t97;
t200 = pkin(9) * t230;
t126 = -t158 + t200;
t55 = t126 * qJD(3);
t216 = qJD(2) * t93;
t188 = t98 * t216;
t165 = t92 * t188;
t164 = t97 * t188;
t178 = t102 * t216;
t27 = -t218 * t164 - qJD(3) * t198 + (t178 + (t92 * t94 + t157) * qJD(3)) * t101;
t127 = t101 * t98 + t171 * t97;
t37 = t127 * t93 + t94 * t230;
t17 = t37 * t27;
t26 = t94 * t189 + (t127 * qJD(3) + (t102 * t97 + t172 * t98) * qJD(2)) * t93;
t251 = 0.2e1 * t57 * t165 + 0.2e1 * t26 * t36 + 0.2e1 * t17;
t250 = t97 ^ 2;
t249 = 0.2e1 * t92;
t248 = 0.2e1 * qJD(4);
t245 = pkin(5) * t96;
t244 = pkin(11) * t96;
t243 = pkin(11) * t100;
t29 = t100 * t57 + t36 * t96;
t12 = qJD(5) * t29 - t26 * t100 + t165 * t96;
t242 = t12 * t143;
t197 = t99 * t230;
t208 = qJD(6) * t95;
t18 = -qJD(6) * t197 - t129 * t99 + t208 * t60 - t84 * t95;
t241 = t18 * t95;
t240 = t18 * t99;
t41 = t95 * t230 + t99 * t60;
t19 = qJD(6) * t41 + t129 * t95 - t84 * t99;
t239 = t19 * t95;
t238 = t19 * t99;
t66 = pkin(9) * t221 + t160;
t56 = t66 * qJD(3);
t237 = t36 * t56;
t39 = -t184 * t96 + (t168 - t189) * t100;
t236 = t39 * t96;
t40 = t60 * t95 - t197;
t235 = t40 * t95;
t234 = t40 * t99;
t233 = t41 * t95;
t232 = t41 * t99;
t231 = t56 * t97;
t229 = t96 * t97;
t24 = t100 * t125 + t96 * t109;
t226 = t88 + t90;
t89 = t96 ^ 2;
t91 = t100 ^ 2;
t225 = t89 - t91;
t224 = t89 + t91;
t222 = t100 * t12;
t220 = t95 * t247;
t219 = t99 * t247;
t213 = qJD(5) * t40;
t212 = qJD(5) * t41;
t211 = qJD(5) * t95;
t210 = qJD(5) * t96;
t209 = qJD(5) * t99;
t207 = qJD(6) * t99;
t205 = qJD(5) * t100;
t204 = qJD(5) * t247;
t203 = qJD(6) * t100;
t202 = qJD(6) * t247;
t201 = 0.2e1 * t59 * t39;
t199 = -0.2e1 * pkin(5) * qJD(6);
t196 = t247 * t230;
t195 = t96 * t220;
t194 = t96 * t219;
t192 = t59 * t211;
t191 = t59 * t209;
t190 = t95 * t210;
t187 = t95 * t207;
t186 = t96 * t209;
t87 = t92 ^ 2;
t185 = t87 * t206;
t183 = t96 * t204;
t182 = t99 * t203;
t181 = t95 * t202;
t180 = t95 * t203;
t179 = t96 * t205;
t177 = t100 * t204;
t175 = t100 * t226;
t174 = t56 * t218;
t170 = t225 * qJD(5);
t167 = t218 * qJD(4);
t83 = 0.2e1 * t179;
t166 = t87 * t188;
t163 = t95 * t186;
t162 = t91 * t187;
t161 = t97 * t185;
t159 = t247 * t84;
t21 = pkin(11) * t230 + t24;
t8 = t108 * t99 - t95 * t21;
t9 = t108 * t95 + t99 * t21;
t155 = t8 * t99 + t9 * t95;
t154 = t8 * t95 - t9 * t99;
t153 = qJD(3) * t92 * t218;
t152 = -t243 + t245;
t151 = pkin(5) * t100 + t244;
t150 = t99 * pkin(5) - t220;
t15 = -t29 * t95 + t37 * t99;
t16 = t29 * t99 + t37 * t95;
t149 = t15 * t99 + t16 * t95;
t148 = t15 * t95 - t16 * t99;
t147 = t233 + t234;
t139 = qJ(4) + t152;
t128 = t99 * t139;
t47 = t128 + t195;
t48 = t139 * t95 - t194;
t146 = t47 * t99 + t48 * t95;
t145 = t47 * t95 - t48 * t99;
t144 = -pkin(3) * t101 - t223;
t142 = t27 * qJ(4) + t37 * qJD(4);
t141 = t97 * t153;
t140 = t101 * t153;
t110 = t100 * t115;
t130 = t96 * qJ(4) + t100 * t246;
t136 = t96 * t156;
t20 = -t110 + (t136 + (-pkin(5) - t130) * t97) * t92;
t11 = t252 * t100 + t253 * t96;
t7 = -pkin(5) * t84 - t11;
t138 = t20 * t207 + t7 * t95;
t137 = t20 * t208 - t7 * t99;
t135 = t12 * t95 - t143 * t207;
t134 = -t12 * t99 - t143 * t208;
t133 = t207 * t59 + t39 * t95;
t132 = t208 * t59 - t99 * t39;
t131 = t205 * t59 + t236;
t105 = t167 - t228 * pkin(11) + t39 * pkin(5) + (t158 + (-t244 - t246) * t230) * qJD(3);
t2 = -t105 * t95 + t208 * t21 + t254 * t99;
t3 = t105 * t99 - t207 * t21 + t254 * t95;
t120 = -qJD(6) * t155 - t2 * t99 - t3 * t95;
t13 = t100 * t165 + t26 * t96 + t256;
t4 = -qJD(6) * t16 - t13 * t95 + t27 * t99;
t5 = qJD(6) * t15 + t13 * t99 + t27 * t95;
t119 = -qJD(6) * t149 - t4 * t95 + t5 * t99;
t118 = -t241 - t238 + (t232 + t235) * qJD(6);
t30 = -t96 * t181 - t95 * (qJD(5) * t151 + qJD(4)) - qJD(6) * t128 + t99 * t177;
t31 = t99 * qJD(4) - t48 * qJD(6) + (t100 * t220 + t151 * t99) * qJD(5);
t117 = -qJD(6) * t146 - t30 * t99 - t31 * t95;
t23 = t110 + (t130 * t97 - t136) * t92;
t116 = -t10 * t96 + t11 * t100 + (t24 * t100 - t23 * t96) * qJD(5);
t1 = -t222 + t13 * t96 + (t100 * t29 - t143 * t96) * qJD(5);
t113 = (t101 * t27 + t26 * t97 + (t101 * t36 - t37 * t97) * qJD(3)) * t92;
t112 = t101 * t166 - t189 * t57 + t218 * t26;
t111 = -t164 * t87 + t218 * t27 - t57 * t84;
t86 = qJ(4) * t248;
t71 = -0.2e1 * t161;
t70 = 0.2e1 * t161;
t64 = -t182 + t190;
t63 = t205 * t95 + t207 * t96;
t62 = -t180 - t186;
t61 = -t205 * t99 + t208 * t96;
t54 = 0.2e1 * (t101 ^ 2 - t250) * t87 * qJD(3);
t53 = t124 + t200;
t52 = (-pkin(2) + t144) * t92;
t51 = -t169 - t66;
t50 = (-t205 * t97 - t206 * t96) * t92;
t49 = (t100 * t206 - t210 * t97) * t92;
t46 = -t167 + t55;
t45 = t100 * t173 + t163;
t44 = (-t214 + (pkin(3) * t97 - t217) * qJD(3)) * t92;
t35 = t167 + (-t246 * t230 + t158) * qJD(3);
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t13 * t29 + 0.2e1 * t17 - 0.2e1 * t242, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t15 * t4 + 0.2e1 * t16 * t5 - 0.2e1 * t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188, -t178, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t111, t113, -pkin(2) * t166 + t126 * t26 + t27 * t66 - t37 * t55 + t237, 0, 0, 0, 0, 0, 0, t113, t112, t111, t165 * t52 + t26 * t53 - t27 * t51 - t37 * t46 + t44 * t57 + t237, 0, 0, 0, 0, 0, 0, t27 * t59 + t37 * t39 + (-t12 * t97 + t143 * t206) * t92, t27 * t60 + t37 * t228 + (-t13 * t97 + (-t101 * t29 + t37 * t229) * qJD(3)) * t92, t12 * t60 - t129 * t143 - t13 * t59 - t29 * t39, -t10 * t29 + t11 * t143 - t12 * t23 + t13 * t24 + t27 * t43 + t35 * t37, 0, 0, 0, 0, 0, 0, t12 * t40 - t143 * t19 + t15 * t39 + t4 * t59, t12 * t41 + t143 * t18 - t16 * t39 - t5 * t59, t15 * t18 - t16 * t19 - t4 * t41 - t40 * t5, t12 * t20 - t143 * t7 + t15 * t3 - t16 * t2 + t4 * t8 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t54, 0.2e1 * t140, t71, -0.2e1 * t141, 0, -0.2e1 * pkin(2) * t215 * t87 - 0.2e1 * t174, -0.2e1 * pkin(2) * t185 + 0.2e1 * t218 * t55 (-t101 * t55 + t231 + (t101 * t126 - t66 * t97) * qJD(3)) * t249, 0.2e1 * t126 * t56 - 0.2e1 * t55 * t66, 0, -0.2e1 * t140, 0.2e1 * t141, t70, t54, t71 (-t101 * t46 + t231 + (t101 * t53 + t51 * t97) * qJD(3)) * t249, 0.2e1 * t174 + 0.2e1 * (t101 * t44 - t215 * t52) * t92, -0.2e1 * t46 * t218 + 0.2e1 * (-t206 * t52 - t44 * t97) * t92, 0.2e1 * t44 * t52 + 0.2e1 * t46 * t51 + 0.2e1 * t53 * t56, 0.2e1 * t60 * t129, -0.2e1 * t129 * t59 - 0.2e1 * t60 * t39 (t228 * t97 + (t92 * t96 * t250 + t101 * t60) * qJD(3)) * t249, t201 (-t206 * t59 - t39 * t97) * t249, t70, 0.2e1 * t35 * t59 + 0.2e1 * t39 * t43 + 0.2e1 * (t11 * t97 + t206 * t23) * t92, 0.2e1 * t35 * t60 + 0.2e1 * t43 * t228 + 0.2e1 * (t10 * t97 + (-t101 * t24 + t43 * t229) * qJD(3)) * t92, 0.2e1 * t10 * t59 - 0.2e1 * t11 * t60 - 0.2e1 * t129 * t23 - 0.2e1 * t24 * t39, -0.2e1 * t10 * t24 + 0.2e1 * t11 * t23 + 0.2e1 * t35 * t43, -0.2e1 * t41 * t18, 0.2e1 * t18 * t40 - 0.2e1 * t19 * t41, -0.2e1 * t18 * t59 + 0.2e1 * t39 * t41, 0.2e1 * t40 * t19, -0.2e1 * t19 * t59 - 0.2e1 * t39 * t40, t201, 0.2e1 * t19 * t20 + 0.2e1 * t3 * t59 + 0.2e1 * t39 * t8 + 0.2e1 * t40 * t7, -0.2e1 * t18 * t20 + 0.2e1 * t2 * t59 - 0.2e1 * t39 * t9 + 0.2e1 * t41 * t7, 0.2e1 * t18 * t8 - 0.2e1 * t19 * t9 + 0.2e1 * t2 * t40 - 0.2e1 * t3 * t41, -0.2e1 * t2 * t9 + 0.2e1 * t20 * t7 + 0.2e1 * t3 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t27, -pkin(3) * t26 + t142, 0, 0, 0, 0, 0, 0, t205 * t37 + t27 * t96, t100 * t27 - t210 * t37, -t1, -t1 * t247 + t142, 0, 0, 0, 0, 0, 0 (t143 * t211 + t4) * t96 + (qJD(5) * t15 + t135) * t100 (t143 * t209 - t5) * t96 + (-qJD(5) * t16 - t134) * t100, t149 * t210 + (qJD(6) * t148 - t4 * t99 - t5 * t95) * t100, t15 * t31 - t16 * t30 + t4 * t47 + t5 * t48 - (-t143 * t210 - t222) * t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, -t189, 0, -t56, t55, 0, 0, 0, -t84, t189, 0, 0, 0 (qJD(3) * t144 + qJD(4) * t101) * t92, t56, 0.2e1 * t167 - t55, -pkin(3) * t56 - qJ(4) * t46 - qJD(4) * t51, -t210 * t60 + t255, -t100 * t39 - t129 * t96 + (-t100 * t60 + t59 * t96) * qJD(5), t49, t131, t50, 0, -t100 * t159 + qJ(4) * t39 + qJD(4) * t59 + t35 * t96 + (t100 * t43 + t196 * t96) * qJD(5), t96 * t159 + qJD(4) * t60 + qJ(4) * t129 + t35 * t100 + (t100 * t196 - t43 * t96) * qJD(5) (t247 * t39 + t10) * t96 + (t129 * t247 - t11) * t100 + ((-t247 * t60 + t23) * t96 + (t247 * t59 - t24) * t100) * qJD(5), t35 * qJ(4) + t43 * qJD(4) - t116 * t247, -t41 * t186 + (-t208 * t41 - t240) * t100, t147 * t210 + (t241 - t238 + (-t232 + t235) * qJD(6)) * t100 (-t18 - t191) * t96 + (-t132 + t212) * t100, -t40 * t190 + (t207 * t40 + t239) * t100 (-t19 + t192) * t96 + (-t133 - t213) * t100, t131, t31 * t59 + t47 * t39 + (t3 + (-t20 * t95 - t247 * t40) * qJD(5)) * t96 + (qJD(5) * t8 + t19 * t247 + t138) * t100, t30 * t59 - t48 * t39 + (t2 + (-t20 * t99 - t247 * t41) * qJD(5)) * t96 + (-qJD(5) * t9 - t18 * t247 - t137) * t100, t18 * t47 - t19 * t48 + t30 * t40 - t31 * t41 + t155 * t210 + (qJD(6) * t154 + t2 * t95 - t3 * t99) * t100, -t2 * t48 + t3 * t47 - t9 * t30 + t8 * t31 - (-t100 * t7 + t20 * t210) * t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, t86, -0.2e1 * t179, 0.2e1 * t170, 0, t83, 0, 0, 0.2e1 * qJ(4) * t205 + 0.2e1 * qJD(4) * t96, -0.2e1 * qJ(4) * t210 + 0.2e1 * qJD(4) * t100, 0, t86, -0.2e1 * t179 * t90 - 0.2e1 * t162, 0.4e1 * t100 * t163 + 0.2e1 * t173 * t91, -0.2e1 * t180 * t96 - 0.2e1 * t209 * t225, -0.2e1 * t179 * t88 + 0.2e1 * t162, 0.2e1 * t170 * t95 - 0.2e1 * t182 * t96, t83, 0.2e1 * t91 * t99 * t202 + 0.2e1 * t31 * t96 + 0.2e1 * (t47 - 0.2e1 * t195) * t205, -0.2e1 * t91 * t181 + 0.2e1 * t30 * t96 + 0.2e1 * (-t48 - 0.2e1 * t194) * t205, 0.2e1 * t146 * t210 + 0.2e1 * (qJD(6) * t145 + t30 * t95 - t31 * t99) * t100, -0.2e1 * t179 * t247 ^ 2 - 0.2e1 * t48 * t30 + 0.2e1 * t47 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-qJD(5) * t148 - t12) * t100 + (t119 - t256) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, 0, t56, 0, 0, 0, 0, 0, 0, t49, t50, -t236 - t255 + (-t100 * t59 + t60 * t96) * qJD(5), t116, 0, 0, 0, 0, 0, 0 (-t19 - t192) * t100 + (-t133 + t213) * t96 (t18 - t191) * t100 + (t132 + t212) * t96 (t233 - t234) * t205 + t118 * t96 (-qJD(5) * t154 - t7) * t100 + (qJD(5) * t20 + t120) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224 * t207, t224 * t208, 0, -t145 * t205 + (t117 + 0.2e1 * t177) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t226) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, 0, 0, 0, 0, 0, 0, 0, 0, t134, t135, t119, -pkin(5) * t12 + pkin(11) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, 0, -t39, t84, t11, t10, 0, 0, t207 * t41 - t241, -qJD(6) * t147 - t239 - t240, t133, t208 * t40 - t238, -t132, 0, -pkin(5) * t19 - pkin(11) * t133 + t137, pkin(5) * t18 + pkin(11) * t132 + t138, pkin(11) * t118 + t120, -pkin(5) * t7 + pkin(11) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, 0, -t205, 0, t183, t177, 0, 0, -t45, -0.4e1 * t180 * t99 + t210 * t227, t63, t45, -t61, 0 (-t100 * t150 - t99 * t244) * qJD(6) + (t152 * t95 + t194) * qJD(5) (t95 * t244 + (t95 * pkin(5) + t219) * t100) * qJD(6) + (t150 * t96 - t99 * t243) * qJD(5), t117, pkin(5) * t183 + pkin(11) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210, -t205, 0, 0, 0, 0, 0, 0, 0, 0, t62, t64, qJD(5) * t175 (pkin(11) * t175 - t245) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t187, -0.2e1 * t173, 0, -0.2e1 * t187, 0, 0, t95 * t199, t99 * t199, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, -t19, t39, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t64, t205, t31, t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t61, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207, 0, -t208, 0, -pkin(11) * t207, pkin(11) * t208, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;

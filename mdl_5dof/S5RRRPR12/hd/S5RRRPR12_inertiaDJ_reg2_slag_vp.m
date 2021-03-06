% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR12_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:35
% EndTime: 2019-12-31 21:40:49
% DurationCPUTime: 4.72s
% Computational Cost: add. (5292->412), mult. (14537->782), div. (0->0), fcn. (13765->10), ass. (0->200)
t101 = sin(qJ(2));
t211 = qJD(2) * t101;
t98 = sin(pkin(5));
t189 = t98 * t211;
t102 = cos(qJ(3));
t218 = cos(pkin(5));
t171 = t218 * qJD(3);
t103 = cos(qJ(2));
t210 = qJD(2) * t103;
t188 = t98 * t210;
t100 = sin(qJ(3));
t223 = t101 * t98;
t193 = qJD(3) * t223;
t76 = t100 * t193;
t241 = -(t171 + t188) * t102 + t76;
t97 = sin(pkin(10));
t99 = cos(pkin(10));
t119 = t97 * t189 - t241 * t99;
t246 = t119 * t97;
t116 = t119 * t100;
t208 = qJD(3) * t102;
t174 = t218 * t100;
t214 = t101 * t102;
t139 = t98 * t214 + t174;
t220 = t103 * t98;
t48 = t139 * t99 - t97 * t220;
t245 = t48 * t208 + t116;
t169 = t99 * t189;
t45 = t241 * t97;
t147 = t45 + t169;
t228 = t139 * t97 + t99 * t220;
t244 = t147 * t100 - t228 * t208;
t232 = sin(qJ(5));
t178 = qJD(5) * t232;
t233 = cos(qJ(5));
t179 = qJD(5) * t233;
t64 = t97 * t178 - t99 * t179;
t207 = qJD(3) * t103;
t243 = t100 * t207 + t102 * t211;
t209 = qJD(3) * t100;
t132 = t139 * qJD(3);
t49 = t100 * t188 + t132;
t173 = t218 * t102;
t66 = t100 * t223 - t173;
t29 = -t102 * t49 + t66 * t209;
t72 = t232 * t97 - t233 * t99;
t161 = t218 * t103 * pkin(1);
t203 = pkin(7) * t223;
t135 = -t161 + t203;
t131 = t135 * qJD(2);
t96 = t102 ^ 2;
t242 = (t100 ^ 2 - t96) * qJD(3);
t94 = t99 ^ 2;
t240 = -0.2e1 * t64;
t198 = t233 * t97;
t73 = t232 * t99 + t198;
t65 = t73 * qJD(5);
t239 = 0.2e1 * t65;
t238 = pkin(3) * t99;
t237 = pkin(8) * t97;
t236 = pkin(8) * t98;
t235 = pkin(8) * t99;
t234 = t66 * pkin(3);
t231 = pkin(3) * t100;
t153 = pkin(2) * t101 - pkin(8) * t103;
t175 = t101 * t218;
t162 = pkin(1) * t175;
t133 = t218 * pkin(8) + t162;
t202 = pkin(7) * t220;
t126 = t133 + t202;
t123 = qJD(3) * t126;
t184 = -pkin(8) * t101 - pkin(1);
t146 = -pkin(2) * t103 + t184;
t142 = t146 * t98;
t136 = qJD(3) * t142;
t199 = t102 * t123 + (-t131 + t136) * t100;
t216 = qJD(2) * t98;
t18 = (-t101 * pkin(3) - t102 * t153) * t216 + t199;
t230 = t18 * t97;
t229 = qJ(4) + pkin(9);
t130 = t102 * t133;
t112 = t130 + (t100 * t184 + (-t100 * pkin(2) + t102 * pkin(7) - qJ(4)) * t103) * t98;
t134 = -t218 * pkin(2) - t161;
t59 = t134 + t203;
t114 = -t139 * qJ(4) + t234 + t59;
t17 = t99 * t112 + t97 * t114;
t217 = qJ(4) * t100;
t145 = -pkin(3) * t102 - pkin(2) - t217;
t87 = t102 * t235;
t57 = t97 * t145 + t87;
t226 = qJ(4) * t49;
t225 = t100 * t97;
t224 = t100 * t99;
t221 = t102 * t97;
t219 = t18 * t100;
t215 = qJD(4) * t99;
t213 = t102 * t103;
t212 = t103 * qJ(4);
t206 = qJD(4) * t100;
t205 = qJD(4) * t102;
t36 = 0.2e1 * t66 * t49;
t204 = -0.2e1 * pkin(2) * qJD(3);
t201 = pkin(8) * t221;
t200 = pkin(8) * t224;
t91 = pkin(8) * t208;
t90 = -pkin(4) * t99 - pkin(3);
t93 = t98 ^ 2;
t195 = t93 * t210;
t191 = t97 * t208;
t190 = t97 * t206;
t187 = qJD(4) * t220;
t186 = t99 * t208;
t185 = t99 * t206;
t182 = t100 * t208;
t180 = qJD(4) * t232;
t177 = t229 * t102;
t176 = t233 * qJD(4);
t172 = qJD(2) * t218;
t170 = 0.2e1 * t182;
t168 = t97 * t186;
t167 = t101 * t195;
t92 = t97 ^ 2;
t166 = 0.2e1 * (t92 + t94) * qJD(4);
t165 = -0.2e1 * t242;
t113 = -t100 * t123 + t102 * t136;
t111 = t113 - t187;
t125 = t49 * pkin(3) - t139 * qJD(4);
t115 = -(t102 * t171 - t76) * qJ(4) + t125;
t124 = -pkin(7) * t214 + t101 * qJ(4) + t100 * t153;
t144 = pkin(7) * t103 - t102 * t212;
t156 = t103 * t173;
t7 = -t97 * t111 + t99 * t115 + ((-t97 * t156 + t99 * t175) * pkin(1) + (-t97 * t124 + t99 * t144) * t98) * qJD(2);
t8 = t99 * t111 + t97 * t115 + ((t99 * t156 + t97 * t175) * pkin(1) + (t99 * t124 + t97 * t144) * t98) * qJD(2);
t160 = -t7 * t97 + t8 * t99;
t159 = t229 * t232;
t158 = t233 * t228;
t157 = t232 * t228;
t154 = t101 * t172;
t151 = -qJ(4) * t102 + t231;
t46 = -t185 + (pkin(8) * t225 + t99 * t151) * qJD(3);
t47 = -t190 + (t97 * t151 - t200) * qJD(3);
t152 = -t46 * t97 + t47 * t99;
t140 = t153 * t216;
t19 = -t100 * t140 + t102 * t131 - t113;
t20 = t102 * t140 - t199;
t150 = -t20 * t100 - t19 * t102;
t149 = pkin(4) + t237 + t238;
t148 = t229 * t198;
t143 = t147 * t99;
t138 = t100 * t211 - t102 * t207;
t68 = t162 + t202;
t34 = -t100 * t126 + t102 * t142;
t28 = pkin(3) * t220 - t34;
t127 = pkin(7) * t213 + t100 * t146;
t122 = (-t229 * t100 - pkin(2)) * t99 - t149 * t102;
t121 = t233 * t122;
t120 = t232 * t122;
t118 = t190 - (-t200 + (-t177 + t231) * t97) * qJD(3);
t117 = -t185 + (t149 * t100 - t99 * t177) * qJD(3);
t110 = pkin(1) * t154 + pkin(7) * t188 + qJ(4) * t241 + t125;
t109 = qJ(4) * t189 - t187 - t19;
t108 = t66 * pkin(4) - t48 * pkin(9) - t112 * t97 + t114 * t99;
t107 = t233 * t108;
t106 = t232 * t108;
t105 = t147 * pkin(9) + t99 * t109 + t97 * t110;
t104 = t49 * pkin(4) - pkin(9) * t119 - t109 * t97 + t110 * t99;
t83 = -0.2e1 * t182;
t79 = t229 * t99;
t75 = -0.2e1 * t167;
t74 = (pkin(4) * t97 + pkin(8)) * t100;
t69 = pkin(4) * t191 + t91;
t63 = qJD(2) * t68;
t62 = t72 * t100;
t61 = t73 * t100;
t56 = t99 * t145 - t201;
t52 = -t97 * t159 + t233 * t79;
t51 = -t232 * t79 - t148;
t50 = -pkin(9) * t225 + t57;
t40 = -t64 * t100 + t73 * t208;
t39 = t100 * t65 - t233 * t186 + t232 * t191;
t38 = -t79 * t179 - t99 * t180 + (qJD(5) * t159 - t176) * t97;
t37 = qJD(5) * t148 - t99 * t176 + t79 * t178 + t97 * t180;
t35 = t127 * t98 + t130;
t25 = t233 * t48 - t157;
t24 = t232 * t48 + t158;
t23 = t233 * t50 + t120;
t22 = -t232 * t50 + t121;
t21 = t228 * pkin(4) + t28;
t16 = -t97 * t130 + t99 * (-qJ(4) * t174 + t134 + t234) + (-t97 * (t127 - t212) + t99 * (pkin(7) * t101 - qJ(4) * t214)) * t98;
t15 = -t45 * pkin(4) + (pkin(8) * t213 + (-pkin(2) * t102 + t90) * t101) * t216 + t199;
t14 = -t228 * pkin(9) + t17;
t13 = -qJD(5) * t120 + t233 * t117 + t232 * t118 - t50 * t179;
t12 = -qJD(5) * t121 - t232 * t117 + t233 * t118 + t50 * t178;
t11 = -qJD(5) * t157 + t232 * t119 - t233 * t147 + t48 * t179;
t10 = qJD(5) * t158 - t233 * t119 - t232 * t147 + t48 * t178;
t4 = t233 * t14 + t106;
t3 = -t232 * t14 + t107;
t2 = -qJD(5) * t106 + t233 * t104 - t232 * t105 - t14 * t179;
t1 = -qJD(5) * t107 - t232 * t104 - t233 * t105 + t14 * t178;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t167, 0.2e1 * (-t101 ^ 2 + t103 ^ 2) * t93 * qJD(2), 0.2e1 * t172 * t220, t75, -0.2e1 * t98 * t154, 0, -0.2e1 * t93 * pkin(1) * t211 - 0.2e1 * t63 * t218, -0.2e1 * pkin(1) * t195 + 0.2e1 * t218 * t131, -0.2e1 * t131 * t220 + 0.2e1 * t135 * t188 - 0.2e1 * t68 * t189 + 0.2e1 * t63 * t223, -0.2e1 * t131 * t68 + 0.2e1 * t135 * t63, -0.2e1 * t139 * t241, -0.2e1 * t139 * t49 + 0.2e1 * t241 * t66, 0.2e1 * t139 * t189 + 0.2e1 * t220 * t241, t36, 0.2e1 * (t103 * t49 - t66 * t211) * t98, t75, 0.2e1 * t59 * t49 + 0.2e1 * t63 * t66 + 0.2e1 * (-t103 * t20 + t34 * t211) * t98, 0.2e1 * t139 * t63 - 0.2e1 * t189 * t35 - 0.2e1 * t19 * t220 - 0.2e1 * t241 * t59, -0.2e1 * t139 * t20 + 0.2e1 * t19 * t66 + 0.2e1 * t241 * t34 - 0.2e1 * t35 * t49, -0.2e1 * t19 * t35 + 0.2e1 * t20 * t34 + 0.2e1 * t59 * t63, 0.2e1 * t48 * t119, -0.2e1 * t119 * t228 + 0.2e1 * t48 * t147, 0.2e1 * t119 * t66 + 0.2e1 * t48 * t49, -0.2e1 * t228 * t147, 0.2e1 * t147 * t66 - 0.2e1 * t228 * t49, t36, -0.2e1 * t28 * t147 + 0.2e1 * t16 * t49 + 0.2e1 * t18 * t228 + 0.2e1 * t7 * t66, 0.2e1 * t119 * t28 - 0.2e1 * t17 * t49 + 0.2e1 * t18 * t48 - 0.2e1 * t8 * t66, -0.2e1 * t119 * t16 + 0.2e1 * t147 * t17 - 0.2e1 * t228 * t8 - 0.2e1 * t7 * t48, 0.2e1 * t16 * t7 + 0.2e1 * t17 * t8 + 0.2e1 * t18 * t28, -0.2e1 * t25 * t10, 0.2e1 * t10 * t24 - 0.2e1 * t11 * t25, -0.2e1 * t10 * t66 + 0.2e1 * t25 * t49, 0.2e1 * t24 * t11, -0.2e1 * t11 * t66 - 0.2e1 * t24 * t49, t36, 0.2e1 * t11 * t21 + 0.2e1 * t15 * t24 + 0.2e1 * t2 * t66 + 0.2e1 * t3 * t49, 0.2e1 * t1 * t66 - 0.2e1 * t10 * t21 + 0.2e1 * t15 * t25 - 0.2e1 * t4 * t49, 0.2e1 * t1 * t24 + 0.2e1 * t10 * t3 - 0.2e1 * t11 * t4 - 0.2e1 * t2 * t25, -0.2e1 * t1 * t4 + 0.2e1 * t15 * t21 + 0.2e1 * t2 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, 0, -t189, 0, -t63, t131, 0, 0, t96 * t193 + (-t76 + (0.2e1 * t171 + t188) * t102) * t100, -t102 * t241 - t66 * t208 + (-t132 - t49) * t100, t138 * t98, t29, t243 * t98, 0, -pkin(2) * t49 - t63 * t102 - t138 * t236 + t59 * t209, pkin(2) * t241 + t63 * t100 + t59 * t208 - t243 * t236, t139 * t91 - t34 * t208 - t35 * t209 + t150 + (-t100 * t241 + t29) * pkin(8), -pkin(2) * t63 + ((-t35 * t100 - t34 * t102) * qJD(3) + t150) * pkin(8), t245 * t99, -t116 * t97 - t191 * t48 + t244 * t99, -t102 * t119 + t186 * t66 + t209 * t48 + t49 * t224, -t244 * t97, -t49 * t225 - t147 * t102 + (-t228 * t100 - t66 * t221) * qJD(3), t29, -t7 * t102 + t46 * t66 + t56 * t49 + (-pkin(8) * t147 + t230) * t100 + (t16 * t100 + (pkin(8) * t228 + t28 * t97) * t102) * qJD(3), t245 * pkin(8) + t8 * t102 - t17 * t209 + t186 * t28 + t219 * t99 - t47 * t66 - t57 * t49, -t119 * t56 + t147 * t57 - t16 * t186 - t17 * t191 - t224 * t7 - t225 * t8 - t228 * t47 - t46 * t48, t16 * t46 + t17 * t47 + t56 * t7 + t57 * t8 + (t208 * t28 + t219) * pkin(8), t10 * t62 - t25 * t39, t10 * t61 + t11 * t62 + t24 * t39 - t25 * t40, t10 * t102 + t209 * t25 - t39 * t66 - t49 * t62, t11 * t61 + t24 * t40, t102 * t11 - t209 * t24 - t40 * t66 - t49 * t61, t29, -t102 * t2 + t11 * t74 + t13 * t66 + t15 * t61 + t209 * t3 + t21 * t40 + t22 * t49 + t24 * t69, -t1 * t102 - t10 * t74 + t12 * t66 - t15 * t62 - t209 * t4 - t21 * t39 - t23 * t49 + t25 * t69, t1 * t61 + t10 * t22 - t11 * t23 + t12 * t24 - t13 * t25 + t2 * t62 + t3 * t39 - t4 * t40, -t1 * t23 - t12 * t4 + t13 * t3 + t15 * t74 + t2 * t22 + t21 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170, t165, 0, t83, 0, 0, t100 * t204, t102 * t204, 0, 0, t94 * t170, -0.4e1 * t100 * t168, 0.2e1 * t99 * t242, t92 * t170, t97 * t165, t83, -0.2e1 * t102 * t46 + 0.2e1 * (t56 + 0.2e1 * t201) * t209, 0.2e1 * t102 * t47 + 0.2e1 * (-t57 + 0.2e1 * t87) * t209, 0.2e1 * (-t46 * t99 - t47 * t97) * t100 + 0.2e1 * (-t56 * t99 - t57 * t97) * t208, 0.2e1 * pkin(8) ^ 2 * t182 + 0.2e1 * t46 * t56 + 0.2e1 * t47 * t57, 0.2e1 * t62 * t39, 0.2e1 * t39 * t61 + 0.2e1 * t40 * t62, 0.2e1 * t102 * t39 - 0.2e1 * t209 * t62, 0.2e1 * t61 * t40, 0.2e1 * t102 * t40 - 0.2e1 * t209 * t61, t83, -0.2e1 * t102 * t13 + 0.2e1 * t209 * t22 + 0.2e1 * t40 * t74 + 0.2e1 * t61 * t69, -0.2e1 * t102 * t12 - 0.2e1 * t209 * t23 - 0.2e1 * t39 * t74 - 0.2e1 * t62 * t69, 0.2e1 * t12 * t61 + 0.2e1 * t13 * t62 + 0.2e1 * t22 * t39 - 0.2e1 * t23 * t40, -0.2e1 * t12 * t23 + 0.2e1 * t13 * t22 + 0.2e1 * t69 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, 0, -t49, t189, t20, t19, 0, 0, t246, -t241 * t94 + (t45 + 0.2e1 * t169) * t97, t97 * t49, t143, t99 * t49, 0, pkin(3) * t147 - t18 * t99 + (-qJD(4) * t66 - t226) * t97, -pkin(3) * t119 - t215 * t66 - t226 * t99 + t230, t97 * qJD(4) * t48 - t215 * t228 + t160 + (t143 + t246) * qJ(4), -pkin(3) * t18 + (-t16 * t97 + t17 * t99) * qJD(4) + t160 * qJ(4), -t10 * t73 - t25 * t64, t10 * t72 - t11 * t73 + t24 * t64 - t25 * t65, t49 * t73 - t64 * t66, t11 * t72 + t24 * t65, -t49 * t72 - t65 * t66, 0, t11 * t90 + t15 * t72 + t21 * t65 + t38 * t66 + t49 * t51, -t10 * t90 + t15 * t73 - t21 * t64 + t37 * t66 - t49 * t52, t1 * t72 + t10 * t51 - t11 * t52 - t2 * t73 + t24 * t37 - t25 * t38 + t3 * t64 - t4 * t65, -t1 * t52 + t15 * t90 + t2 * t51 + t3 * t38 - t37 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, 0, -t209, 0, -t91, pkin(8) * t209, 0, 0, t168, (-t92 + t94) * t208, t97 * t209, -t168, t99 * t209, 0, t97 * t205 + (-t97 * t217 + (-pkin(3) * t97 - t235) * t102) * qJD(3), t99 * t205 + (-t99 * t217 + (t237 - t238) * t102) * qJD(3), t152, -pkin(3) * t91 + (-t56 * t97 + t57 * t99) * qJD(4) + t152 * qJ(4), -t39 * t73 + t62 * t64, t39 * t72 - t40 * t73 + t61 * t64 + t62 * t65, t102 * t64 + t209 * t73, t40 * t72 + t61 * t65, t102 * t65 - t209 * t72, 0, -t102 * t38 + t209 * t51 + t40 * t90 + t65 * t74 + t69 * t72, -t102 * t37 - t209 * t52 - t39 * t90 - t64 * t74 + t69 * t73, t12 * t72 - t13 * t73 + t22 * t64 - t23 * t65 + t37 * t61 + t38 * t62 + t39 * t51 - t40 * t52, -t12 * t52 + t13 * t51 + t22 * t38 - t23 * t37 + t69 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, qJ(4) * t166, t73 * t240, 0.2e1 * t64 * t72 - 0.2e1 * t65 * t73, 0, t72 * t239, 0, 0, t90 * t239, t90 * t240, 0.2e1 * t37 * t72 - 0.2e1 * t38 * t73 + 0.2e1 * t51 * t64 - 0.2e1 * t52 * t65, -0.2e1 * t37 * t52 + 0.2e1 * t38 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t119, 0, t18, 0, 0, 0, 0, 0, 0, t11, -t10, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, t186, 0, t91, 0, 0, 0, 0, 0, 0, t40, -t39, 0, t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t64, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, -t11, t49, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, 0, -t40, t209, t13, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, 0, -t65, 0, t38, t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;

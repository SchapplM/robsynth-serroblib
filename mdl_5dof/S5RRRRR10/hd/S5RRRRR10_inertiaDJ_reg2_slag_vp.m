% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR10_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:35:58
% EndTime: 2019-12-31 22:36:13
% DurationCPUTime: 5.14s
% Computational Cost: add. (7856->387), mult. (20574->716), div. (0->0), fcn. (19934->10), ass. (0->191)
t125 = cos(qJ(2));
t217 = cos(pkin(5));
t189 = pkin(1) * t217;
t118 = sin(pkin(5));
t122 = sin(qJ(2));
t215 = t118 * t122;
t143 = pkin(7) * t215 - t125 * t189;
t138 = t143 * qJD(2);
t153 = (-pkin(2) * t125 - pkin(8) * t122 - pkin(1)) * t118;
t241 = -qJD(3) * t153 + t138;
t208 = qJD(2) * t125;
t181 = t118 * t208;
t240 = t217 * qJD(3) + t181;
t119 = sin(qJ(5));
t116 = t119 ^ 2;
t123 = cos(qJ(5));
t117 = t123 ^ 2;
t211 = t116 - t117;
t238 = t211 * qJD(5);
t121 = sin(qJ(3));
t124 = cos(qJ(3));
t176 = t217 * t124;
t147 = t121 * t215 - t176;
t233 = cos(qJ(4));
t137 = t233 * t147;
t120 = sin(qJ(4));
t182 = qJD(3) * t215;
t198 = t240 * t121 + t124 * t182;
t204 = qJD(4) * t120;
t64 = -t121 * t182 + t240 * t124;
t80 = t217 * t121 + t124 * t215;
t140 = t120 * t198 + t80 * t204 - t233 * t64;
t131 = qJD(4) * t137 + t140;
t55 = -t120 * t147 + t233 * t80;
t31 = qJD(4) * t55 + t120 * t64 + t233 * t198;
t214 = t118 * t125;
t83 = pkin(7) * t214 + t122 * t189;
t77 = qJD(2) * t83;
t51 = t198 * pkin(3) + t77;
t127 = t31 * pkin(4) + pkin(10) * t131 + t51;
t203 = qJD(5) * t119;
t209 = qJD(2) * t122;
t104 = t118 * t209;
t54 = t120 * t80 + t137;
t74 = -t217 * pkin(2) + t143;
t57 = pkin(3) * t147 + t74;
t130 = t54 * pkin(4) - t55 * pkin(10) + t57;
t174 = pkin(3) * t104;
t210 = qJD(2) * t118;
t152 = (pkin(2) * t122 - pkin(8) * t125) * t210;
t206 = qJD(3) * t124;
t75 = t217 * pkin(8) + t83;
t33 = t241 * t121 + t124 * t152 - t75 * t206;
t128 = -t64 * pkin(9) + t174 + t33;
t207 = qJD(3) * t121;
t32 = -t121 * t152 + t241 * t124 + t75 * t207;
t136 = t198 * pkin(9) + t32;
t177 = qJD(4) * t233;
t199 = pkin(3) * t214;
t52 = -t121 * t75 + t124 * t153;
t38 = -t80 * pkin(9) - t199 + t52;
t53 = t121 * t153 + t124 * t75;
t42 = -pkin(9) * t147 + t53;
t9 = -t120 * t128 + t233 * t136 - t38 * t177 + t42 * t204;
t239 = -pkin(10) * t104 - qJD(5) * t130 + t9;
t39 = t233 * t42;
t229 = t120 * t38 + t39;
t25 = -pkin(10) * t214 + t229;
t2 = -t119 * t127 + t239 * t123 + t25 * t203;
t1 = t2 * t123;
t11 = -t119 * t25 + t123 * t130;
t12 = t119 * t130 + t123 * t25;
t162 = t11 * t123 + t119 * t12;
t114 = qJD(5) * t123;
t3 = -t25 * t114 + t239 * t119 + t123 * t127;
t135 = -t162 * qJD(5) - t3 * t119 - t1;
t237 = qJD(3) + qJD(4);
t133 = -t120 * t136 - t233 * t128;
t10 = -qJD(4) * t229 - t133;
t236 = 0.2e1 * t118;
t235 = -pkin(9) - pkin(8);
t26 = -t120 * t42 + t233 * t38;
t24 = pkin(4) * t214 - t26;
t23 = t24 * t114;
t8 = -pkin(4) * t104 - t10;
t234 = t8 * t119 + t23;
t232 = pkin(8) * t118;
t184 = t233 * t124;
t169 = qJD(3) * t184;
t213 = t120 * t121;
t99 = t235 * t124;
t66 = t235 * t213 - t233 * t99;
t188 = qJD(3) * t235;
t89 = t121 * t188;
t46 = t66 * qJD(4) + t120 * t89 - t235 * t169;
t185 = t233 * t121;
t65 = -t120 * t99 - t235 * t185;
t231 = t46 * t65;
t62 = -t124 * t177 + t237 * t213 - t169;
t212 = t120 * t124;
t87 = t185 + t212;
t230 = t87 * t62;
t61 = t65 * t114;
t228 = t46 * t119 + t61;
t227 = pkin(3) * qJD(4);
t47 = t119 * t55 + t123 * t214;
t18 = qJD(5) * t47 - t119 * t104 + t123 * t131;
t226 = t119 * t18;
t225 = t119 * t47;
t224 = t120 * t65;
t223 = t121 * t64;
t190 = t119 * t214;
t19 = -qJD(5) * t190 - t123 * t104 + t55 * t114 - t119 * t131;
t222 = t123 * t19;
t48 = t123 * t55 - t190;
t221 = t123 * t48;
t220 = t123 * t62;
t219 = t77 * t121;
t112 = -t233 * pkin(3) - pkin(4);
t194 = pkin(3) * t204;
t218 = t112 * t114 + t119 * t194;
t205 = qJD(3) * t125;
t202 = 0.2e1 * t54 * t31;
t63 = t237 * t87;
t86 = -t184 + t213;
t201 = 0.2e1 * t86 * t63;
t200 = -0.2e1 * pkin(2) * qJD(3);
t197 = pkin(3) * t207;
t196 = pkin(4) * t203;
t195 = pkin(4) * t114;
t193 = t119 * t220;
t192 = t87 * t203;
t191 = t87 * t114;
t22 = t24 * t203;
t60 = t65 * t203;
t113 = -pkin(3) * t124 - pkin(2);
t187 = t119 * t233;
t186 = t123 * t233;
t115 = t118 ^ 2;
t183 = t115 * t208;
t180 = t119 * t114;
t179 = t121 * t206;
t178 = -t8 * t123 + t22;
t173 = pkin(3) * t177;
t84 = t87 ^ 2;
t172 = t84 * t180;
t171 = t122 * t183;
t170 = t198 * t124;
t168 = -t24 * t62 + t8 * t87;
t167 = t217 * t210;
t166 = t87 * t31 - t62 * t54;
t165 = t31 * t86 + t54 * t63;
t164 = t46 * t87 - t62 * t65;
t163 = t62 * t86 - t63 * t87;
t111 = t120 * pkin(3) + pkin(10);
t161 = t111 * t86 - t112 * t87;
t151 = -pkin(4) * t86 + pkin(10) * t87 - t113;
t142 = t123 * t151;
t40 = -t119 * t66 - t142;
t41 = -t119 * t151 + t123 * t66;
t160 = t119 * t41 + t123 * t40;
t159 = t119 * t48 + t123 * t47;
t158 = t112 * t203 - t123 * t194;
t21 = t54 * t114 + t119 * t31;
t156 = -t123 * t31 + t54 * t203;
t155 = t192 + t220;
t154 = -t123 * t63 + t86 * t203;
t150 = (t116 + t117) * t233;
t149 = t121 * t205 + t124 * t209;
t148 = t121 * t209 - t124 * t205;
t146 = (t120 * t87 - t233 * t86) * qJD(4);
t145 = pkin(4) * t63 + pkin(10) * t62 + t197;
t141 = t147 * t121;
t45 = t65 * qJD(4) - t188 * t212 - t233 * t89;
t16 = qJD(5) * t142 - t119 * t145 + t123 * t45 + t66 * t203;
t17 = -qJD(5) * t41 + t119 * t45 + t123 * t145;
t4 = -t160 * qJD(5) - t119 * t17 - t123 * t16;
t134 = -t33 * t121 - t32 * t124 + (-t53 * t121 - t52 * t124) * qJD(3);
t132 = pkin(3) * t146 - t111 * t63 - t112 * t62;
t103 = -0.2e1 * t180;
t102 = 0.2e1 * t180;
t90 = -0.2e1 * t171;
t85 = -0.2e1 * t238;
t79 = t150 * t227;
t50 = t86 * t114 + t119 * t63;
t37 = t87 * t238 + t193;
t34 = -0.4e1 * t87 * t180 + t211 * t62;
t15 = t47 * t203 - t222;
t14 = t48 * t114 - t226;
t5 = -t159 * qJD(5) - t119 * t19 - t123 * t18;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t171, 0.2e1 * (-t122 ^ 2 + t125 ^ 2) * t115 * qJD(2), 0.2e1 * t125 * t167, t90, -0.2e1 * t122 * t167, 0, -0.2e1 * t115 * pkin(1) * t209 - 0.2e1 * t77 * t217, -0.2e1 * pkin(1) * t183 + 0.2e1 * t217 * t138, -0.2e1 * t104 * t83 - 0.2e1 * t138 * t214 + 0.2e1 * t143 * t181 + 0.2e1 * t77 * t215, -0.2e1 * t138 * t83 + 0.2e1 * t143 * t77, 0.2e1 * t80 * t64, -0.2e1 * t147 * t64 - 0.2e1 * t198 * t80, (-t125 * t64 + t209 * t80) * t236, 0.2e1 * t147 * t198, -0.2e1 * t104 * t147 + 0.2e1 * t198 * t214, t90, -0.2e1 * t77 * t176 + 0.2e1 * t74 * t198 + 0.2e1 * (-t33 * t125 + (qJD(2) * t52 + t219) * t122) * t118, 0.2e1 * t74 * t64 + 0.2e1 * t77 * t80 + 0.2e1 * (-t125 * t32 - t209 * t53) * t118, 0.2e1 * t147 * t32 - 0.2e1 * t198 * t53 - 0.2e1 * t33 * t80 - 0.2e1 * t52 * t64, -0.2e1 * t32 * t53 + 0.2e1 * t33 * t52 + 0.2e1 * t74 * t77, -0.2e1 * t55 * t131, 0.2e1 * t131 * t54 - 0.2e1 * t55 * t31, 0.2e1 * t104 * t55 + 0.2e1 * t131 * t214, t202, (t125 * t31 - t209 * t54) * t236, t90, 0.2e1 * t57 * t31 + 0.2e1 * t51 * t54 + 0.2e1 * (-t10 * t125 + t209 * t26) * t118, -0.2e1 * t104 * t229 - 0.2e1 * t131 * t57 - 0.2e1 * t9 * t214 + 0.2e1 * t51 * t55, -0.2e1 * t10 * t55 + 0.2e1 * t131 * t26 - 0.2e1 * t229 * t31 + 0.2e1 * t9 * t54, 0.2e1 * t10 * t26 - 0.2e1 * t229 * t9 + 0.2e1 * t51 * t57, -0.2e1 * t48 * t18, 0.2e1 * t18 * t47 - 0.2e1 * t19 * t48, -0.2e1 * t18 * t54 + 0.2e1 * t31 * t48, 0.2e1 * t47 * t19, -0.2e1 * t19 * t54 - 0.2e1 * t31 * t47, t202, 0.2e1 * t11 * t31 + 0.2e1 * t19 * t24 + 0.2e1 * t3 * t54 + 0.2e1 * t47 * t8, -0.2e1 * t12 * t31 - 0.2e1 * t18 * t24 + 0.2e1 * t2 * t54 + 0.2e1 * t48 * t8, 0.2e1 * t11 * t18 - 0.2e1 * t12 * t19 + 0.2e1 * t2 * t47 - 0.2e1 * t3 * t48, 0.2e1 * t11 * t3 - 0.2e1 * t12 * t2 + 0.2e1 * t24 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, 0, -t104, 0, -t77, t138, 0, 0, t206 * t80 + t223, -t121 * t198 + t64 * t124 + (-t80 * t121 - t124 * t147) * qJD(3), t148 * t118, qJD(3) * t141 - t170, t149 * t118, 0, -pkin(2) * t198 - t77 * t124 - t148 * t232 + t207 * t74, -pkin(2) * t64 - t149 * t232 + t206 * t74 + t219, (-t170 + t223 + (t124 * t80 + t141) * qJD(3)) * pkin(8) + t134, -pkin(2) * t77 + pkin(8) * t134, -t131 * t87 - t55 * t62, t131 * t86 - t55 * t63 - t166, (t125 * t62 + t209 * t87) * t118, t165, (t125 * t63 - t209 * t86) * t118, 0, t54 * t197 + t113 * t31 + t51 * t86 + t57 * t63 + (t125 * t46 - t209 * t65) * t118, -t104 * t66 - t113 * t131 + t197 * t55 - t45 * t214 + t51 * t87 - t57 * t62, -t10 * t87 - t131 * t65 - t229 * t63 + t26 * t62 - t66 * t31 + t45 * t54 + t46 * t55 + t9 * t86, -t10 * t65 + t113 * t51 + t197 * t57 - t229 * t45 - t26 * t46 - t66 * t9, -t48 * t192 + (-t18 * t87 - t48 * t62) * t123, t159 * t62 + (t226 - t222 + (-t221 + t225) * qJD(5)) * t87, t123 * t166 - t18 * t86 - t192 * t54 + t48 * t63, t47 * t191 + (t19 * t87 - t47 * t62) * t119, -t119 * t166 - t19 * t86 - t191 * t54 - t47 * t63, t165, t11 * t63 + t119 * t168 + t17 * t54 + t19 * t65 + t23 * t87 + t3 * t86 + t31 * t40 + t46 * t47, -t12 * t63 + t123 * t168 + t16 * t54 - t18 * t65 + t2 * t86 - t22 * t87 - t31 * t41 + t46 * t48, t16 * t47 - t17 * t48 + t18 * t40 - t19 * t41 + t162 * t62 + (t119 * t2 - t123 * t3 + (t11 * t119 - t12 * t123) * qJD(5)) * t87, t11 * t17 - t12 * t16 - t2 * t41 + t24 * t46 + t3 * t40 + t65 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t179, 0.2e1 * (-t121 ^ 2 + t124 ^ 2) * qJD(3), 0, -0.2e1 * t179, 0, 0, t121 * t200, t124 * t200, 0, 0, -0.2e1 * t230, 0.2e1 * t163, 0, t201, 0, 0, 0.2e1 * t113 * t63 + 0.2e1 * t197 * t86, -0.2e1 * t113 * t62 + 0.2e1 * t197 * t87, 0.2e1 * t45 * t86 - 0.2e1 * t63 * t66 + 0.2e1 * t164, 0.2e1 * t113 * t197 - 0.2e1 * t45 * t66 + 0.2e1 * t231, -0.2e1 * t117 * t230 - 0.2e1 * t172, 0.4e1 * t193 * t87 + 0.2e1 * t84 * t238, -0.2e1 * t123 * t163 - 0.2e1 * t192 * t86, -0.2e1 * t116 * t230 + 0.2e1 * t172, 0.2e1 * t119 * t163 - 0.2e1 * t191 * t86, t201, 0.2e1 * t119 * t164 + 0.2e1 * t17 * t86 + 0.2e1 * t40 * t63 + 0.2e1 * t61 * t87, 0.2e1 * t123 * t164 + 0.2e1 * t16 * t86 - 0.2e1 * t41 * t63 - 0.2e1 * t60 * t87, 0.2e1 * t160 * t62 + 0.2e1 * (t119 * t16 - t123 * t17 + (t119 * t40 - t123 * t41) * qJD(5)) * t87, -0.2e1 * t16 * t41 + 0.2e1 * t17 * t40 + 0.2e1 * t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, -t198, t104, t33, t32, 0, 0, 0, 0, -t131, 0, -t31, t104, t233 * t174 + (-t39 + (-t38 + t199) * t120) * qJD(4) - t133, -t120 * t174 + t173 * t214 + t9, (-t120 * t31 + t233 * t140 + (t120 * t55 + (-t54 + t137) * t233) * qJD(4)) * pkin(3), (t233 * t10 - t120 * t9 + (-t120 * t26 + t229 * t233) * qJD(4)) * pkin(3), t14, t5, t21, t15, -t156, 0, t112 * t19 - t21 * t111 + (t120 * t47 - t187 * t54) * t227 + t178, -t112 * t18 + t156 * t111 + (t120 * t48 - t186 * t54) * t227 + t234, -t1 + (-t47 * t173 - t111 * t19 + (t111 * t48 - t11) * qJD(5)) * t123 + (t48 * t173 - t111 * t18 - t3 + (t111 * t47 - t12) * qJD(5)) * t119, t8 * t112 + (-t11 * t187 + t12 * t186 + t120 * t24) * t227 + t135 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, 0, -t207, 0, -pkin(8) * t206, pkin(8) * t207, 0, 0, 0, 0, -t62, 0, -t63, 0, -t46, t45, (-t120 * t63 + t233 * t62 + t146) * pkin(3), (-t233 * t46 - t120 * t45 + (t233 * t66 + t224) * qJD(4)) * pkin(3), -t37, t34, t50, t37, -t154, 0, t60 + (-qJD(5) * t161 - t46) * t123 + t132 * t119, t123 * t132 + t161 * t203 + t228, t4, t46 * t112 + (t186 * t41 - t187 * t40 + t224) * t227 + t4 * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t194, -0.2e1 * t173, 0, 0, t102, t85, 0, t103, 0, 0, 0.2e1 * t158, 0.2e1 * t218, 0.2e1 * t79, 0.2e1 * (t111 * t150 + t112 * t120) * t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, 0, -t31, t104, t10, t9, 0, 0, t14, t5, t21, t15, -t156, 0, -pkin(4) * t19 - pkin(10) * t21 + t178, pkin(4) * t18 + pkin(10) * t156 + t234, (-t226 - t222 + (t221 + t225) * qJD(5)) * pkin(10) + t135, -pkin(4) * t8 + pkin(10) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, 0, -t63, 0, -t46, t45, 0, 0, -t37, t34, t50, t37, -t154, 0, t60 + (pkin(4) * t62 - pkin(10) * t63) * t119 + (-t46 + (-pkin(4) * t87 - pkin(10) * t86) * qJD(5)) * t123, pkin(4) * t155 + pkin(10) * t154 + t228, t4, -pkin(4) * t46 + pkin(10) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t194, -t173, 0, 0, t102, t85, 0, t103, 0, 0, t158 - t196, -t195 + t218, t79, (-pkin(4) * t120 + pkin(10) * t150) * t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t85, 0, t103, 0, 0, -0.2e1 * t196, -0.2e1 * t195, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, -t19, t31, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, 0, t119 * t62 - t191, t63, t17, t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, -t203, 0, -t111 * t114 - t119 * t173, t111 * t203 - t123 * t173, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, -t203, 0, -pkin(10) * t114, pkin(10) * t203, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
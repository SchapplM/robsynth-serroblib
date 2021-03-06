% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPPR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:36:01
% EndTime: 2019-03-09 15:36:14
% DurationCPUTime: 5.08s
% Computational Cost: add. (5435->376), mult. (13077->656), div. (0->0), fcn. (11532->8), ass. (0->170)
t124 = sin(qJ(3));
t126 = cos(qJ(3));
t214 = sin(pkin(10));
t215 = cos(pkin(10));
t96 = t214 * t124 - t215 * t126;
t145 = -t215 * t124 - t214 * t126;
t228 = t145 * qJD(3);
t92 = t96 * qJD(3);
t237 = 0.2e1 * t145 * t228 - 0.2e1 * t92 * t96;
t236 = -0.2e1 * t228 * t96;
t125 = sin(qJ(2));
t127 = cos(qJ(2));
t202 = t127 * qJD(2);
t55 = t125 * t92 + t145 * t202;
t175 = qJD(2) * t214;
t176 = qJD(2) * t215;
t210 = t126 * t127;
t212 = t124 * t127;
t56 = -t125 * t228 + t175 * t212 - t176 * t210;
t86 = t145 * t125;
t87 = t96 * t125;
t235 = t145 * t55 + t228 * t87 - t56 * t96 + t92 * t86;
t234 = t228 * t86 - t55 * t96;
t116 = t125 * qJD(2);
t233 = t96 * t116 + t127 * t228;
t185 = t126 * t202;
t207 = qJD(3) * t124;
t230 = -t125 * t207 + t185;
t205 = qJD(3) * t127;
t190 = t126 * t205;
t229 = t124 * t116 - t190;
t120 = t125 ^ 2;
t173 = (-t127 ^ 2 + t120) * qJD(2);
t227 = -0.2e1 * qJD(3);
t226 = 2 * qJD(5);
t225 = pkin(4) + pkin(5);
t224 = cos(qJ(6));
t223 = pkin(2) * t127;
t222 = t124 * pkin(7);
t221 = t125 * pkin(3);
t220 = t126 * pkin(2);
t217 = -qJ(4) - pkin(8);
t193 = -pkin(3) - t222;
t211 = t125 * t126;
t166 = -t125 * pkin(8) - t223;
t160 = -pkin(1) + t166;
t99 = t126 * t160;
t65 = -qJ(4) * t211 + t193 * t127 + t99;
t213 = t124 * t125;
t111 = pkin(7) * t210;
t78 = t124 * t160 + t111;
t70 = -qJ(4) * t213 + t78;
t34 = t214 * t65 + t215 * t70;
t152 = qJD(3) * t160;
t165 = pkin(2) * t125 - pkin(8) * t127;
t155 = t165 * qJD(2);
t216 = -t124 * t155 - t126 * t152;
t72 = t96 * t217;
t103 = pkin(3) * t213 + t125 * pkin(7);
t119 = t124 ^ 2;
t121 = t126 ^ 2;
t209 = t119 - t121;
t206 = qJD(3) * t126;
t123 = sin(qJ(6));
t204 = qJD(6) * t123;
t203 = t125 * qJD(4);
t201 = t127 * qJD(5);
t200 = 0.2e1 * t86 * t55;
t199 = -0.2e1 * pkin(1) * qJD(2);
t198 = pkin(2) * t227;
t197 = pkin(7) * t212;
t196 = pkin(3) * t207;
t195 = pkin(7) * t202;
t115 = -pkin(3) * t126 - pkin(2);
t191 = t124 * t205;
t188 = t124 * t206;
t187 = t125 * t202;
t186 = t126 * t203;
t181 = qJD(3) * t217;
t148 = qJD(4) * t126 + t124 * t181;
t149 = -qJD(4) * t124 + t126 * t181;
t57 = t214 * t148 - t215 * t149;
t58 = t215 * t148 + t214 * t149;
t71 = t145 * t217;
t184 = t57 * t71 + t72 * t58;
t183 = qJD(6) * t224;
t30 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t211 + (-t203 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t127) * t124 - t216;
t182 = t215 * t30;
t174 = qJD(3) * t209;
t172 = 0.2e1 * t187;
t171 = t125 * t185;
t170 = t120 * t188;
t26 = -qJ(5) * t127 + t34;
t114 = -t215 * pkin(3) - pkin(4);
t169 = -qJ(5) * t87 - t103;
t164 = -t55 * t87 - t56 * t86;
t132 = -t186 + (-t111 + (-t217 * t125 + pkin(1) + t223) * t124) * qJD(3) + (t217 * t210 + (-t193 + t220) * t125) * qJD(2);
t163 = -t215 * t132 + t214 * t30;
t77 = t99 - t197;
t162 = -t124 * t78 - t126 * t77;
t161 = -pkin(5) + t114;
t159 = -qJ(5) * t145 - t115;
t13 = t214 * t132 + t182;
t130 = qJ(5) * t116 + t13;
t10 = t130 - t201;
t129 = -t55 * pkin(9) + t10;
t147 = t214 * t70 - t215 * t65;
t136 = t87 * pkin(9) + t225 * t127 + t147;
t133 = t224 * t136;
t139 = t56 * pkin(9) - t225 * t116 + t163;
t18 = -t86 * pkin(9) + t26;
t1 = -qJD(6) * t133 - t123 * t139 - t224 * t129 + t18 * t204;
t51 = -t123 * t86 - t224 * t87;
t64 = t123 * t96 - t145 * t224;
t158 = -t86 * t116 + t127 * t55;
t157 = -t71 * t116 + t57 * t127;
t156 = t72 * t116 - t127 * t58;
t40 = -pkin(4) * t228 + t92 * qJ(5) + qJD(5) * t145 + t196;
t154 = t72 * t55 - t71 * t56 - t57 * t87 + t58 * t86;
t151 = t126 * t116 + t191;
t150 = t124 * t202 + t125 * t206;
t146 = t224 * t161;
t144 = -0.2e1 * t145 * t57 + 0.2e1 * t228 * t72 - 0.2e1 * t58 * t96 - 0.2e1 * t71 * t92;
t143 = pkin(9) * t145 + t71;
t140 = t123 * t143;
t112 = t214 * pkin(3) + qJ(5);
t80 = t224 * t112 + t123 * t161;
t44 = t151 * pkin(7) + t216;
t45 = -t78 * qJD(3) + (pkin(7) * t213 + t126 * t165) * qJD(2);
t138 = t162 * qJD(3) - t124 * t45 - t126 * t44;
t137 = t224 * t143;
t75 = t150 * pkin(3) + t195;
t135 = t92 * pkin(9) + t57;
t134 = t123 * t136;
t16 = -t55 * pkin(4) + t56 * qJ(5) + t87 * qJD(5) + t75;
t8 = t224 * t18 + t134;
t131 = -pkin(9) * t228 + t58;
t128 = -t123 * t129 + t224 * t139;
t108 = -0.2e1 * t187;
t79 = -t123 * t112 + t146;
t73 = -t124 * t185 + t125 * t174;
t69 = t123 * qJD(5) + t80 * qJD(6);
t68 = -t224 * qJD(5) - qJD(6) * t146 + t112 * t204;
t67 = 0.2e1 * t145 * t92;
t63 = -t123 * t145 - t224 * t96;
t61 = pkin(4) * t96 - t159;
t60 = -t116 * t145 + t127 * t92;
t59 = pkin(9) * t96 + t72;
t50 = -t123 * t87 + t224 * t86;
t43 = -t225 * t96 + t159;
t42 = -pkin(4) * t86 - t169;
t38 = 0.2e1 * t87 * t56;
t29 = -0.2e1 * t87 * t116 + 0.2e1 * t127 * t56;
t28 = t127 * pkin(4) + t147;
t27 = t225 * t86 + t169;
t25 = -pkin(5) * t228 + t40;
t24 = t64 * qJD(6) - t123 * t92 + t224 * t228;
t23 = t123 * t228 - t145 * t204 - t96 * t183 + t224 * t92;
t21 = t224 * t59 + t140;
t20 = -t123 * t59 + t137;
t19 = t145 * t56 + t87 * t92;
t15 = t51 * qJD(6) - t123 * t56 + t224 * t55;
t14 = t123 * t55 + t86 * t183 - t87 * t204 + t224 * t56;
t11 = -pkin(4) * t116 + t163;
t9 = -t55 * pkin(5) + t16;
t7 = -t123 * t18 + t133;
t4 = qJD(6) * t140 + t123 * t131 - t224 * t135 + t59 * t183;
t3 = -qJD(6) * t137 - t123 * t135 - t224 * t131 + t59 * t204;
t2 = -t8 * qJD(6) + t128;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, -0.2e1 * t173, 0, t108, 0, 0, t125 * t199, t127 * t199, 0, 0, 0.2e1 * t121 * t187 - 0.2e1 * t170, 0.2e1 * t120 * t174 - 0.4e1 * t124 * t171, 0.2e1 * t125 * t191 + 0.2e1 * t126 * t173, 0.2e1 * t119 * t187 + 0.2e1 * t170, -0.2e1 * t124 * t173 + 0.2e1 * t125 * t190, t108, 0.2e1 * t77 * t116 - 0.2e1 * t127 * t45 + 0.2e1 * (t120 * t206 + t124 * t172) * pkin(7), -0.2e1 * t78 * t116 - 0.2e1 * t127 * t44 + 0.2e1 * (-t120 * t207 + 0.2e1 * t171) * pkin(7), 0.2e1 * t162 * t202 + 0.2e1 * (t124 * t44 - t126 * t45 + (t124 * t77 - t126 * t78) * qJD(3)) * t125, 0.2e1 * pkin(7) ^ 2 * t187 - 0.2e1 * t78 * t44 + 0.2e1 * t77 * t45, t38, 0.2e1 * t164, t29, t200, -0.2e1 * t158, t108, -0.2e1 * t103 * t55 - 0.2e1 * t116 * t147 + 0.2e1 * t127 * t163 - 0.2e1 * t75 * t86, -0.2e1 * t103 * t56 - 0.2e1 * t116 * t34 + 0.2e1 * t127 * t13 - 0.2e1 * t75 * t87, 0.2e1 * t13 * t86 - 0.2e1 * t147 * t56 - 0.2e1 * t163 * t87 + 0.2e1 * t34 * t55, 0.2e1 * t103 * t75 + 0.2e1 * t13 * t34 + 0.2e1 * t147 * t163, t38, t29, -0.2e1 * t164, t108, 0.2e1 * t158, t200, 0.2e1 * t11 * t127 - 0.2e1 * t116 * t28 - 0.2e1 * t16 * t86 - 0.2e1 * t42 * t55, 0.2e1 * t10 * t86 - 0.2e1 * t11 * t87 + 0.2e1 * t26 * t55 - 0.2e1 * t28 * t56, -0.2e1 * t10 * t127 + 0.2e1 * t116 * t26 + 0.2e1 * t16 * t87 + 0.2e1 * t42 * t56, 0.2e1 * t10 * t26 + 0.2e1 * t11 * t28 + 0.2e1 * t16 * t42, -0.2e1 * t51 * t14, 0.2e1 * t14 * t50 - 0.2e1 * t15 * t51, -0.2e1 * t116 * t51 - 0.2e1 * t127 * t14, 0.2e1 * t50 * t15, 0.2e1 * t116 * t50 - 0.2e1 * t127 * t15, t108, -0.2e1 * t116 * t7 + 0.2e1 * t127 * t2 + 0.2e1 * t27 * t15 - 0.2e1 * t9 * t50, 0.2e1 * t1 * t127 + 0.2e1 * t116 * t8 - 0.2e1 * t27 * t14 - 0.2e1 * t9 * t51, 0.2e1 * t1 * t50 + 0.2e1 * t14 * t7 - 0.2e1 * t15 * t8 - 0.2e1 * t2 * t51, -0.2e1 * t1 * t8 + 0.2e1 * t2 * t7 - 0.2e1 * t27 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, 0, -t116, 0, -t195, pkin(7) * t116, 0, 0, -t73, -0.4e1 * t125 * t188 - t202 * t209, t229, t73, t151, 0 (pkin(8) * t210 + (-t220 + t222) * t125) * qJD(3) + (t124 * t166 - t111) * qJD(2) (pkin(7) * t211 + t124 * t165) * qJD(3) + (t126 * t166 + t197) * qJD(2), t138, -pkin(2) * t195 + pkin(8) * t138, t19, -t235, t60, t234, -t233, 0, -t103 * t228 - t115 * t55 - t196 * t86 + t75 * t96 + t157, -t103 * t92 - t115 * t56 - t145 * t75 - t196 * t87 - t156, -t13 * t96 - t145 * t163 - t147 * t92 + t228 * t34 + t154, t103 * t196 + t115 * t75 + t13 * t72 + t147 * t57 + t163 * t71 + t34 * t58, t19, t60, t235, 0, t233, t234, t16 * t96 - t228 * t42 - t40 * t86 - t61 * t55 + t157, -t10 * t96 - t11 * t145 + t228 * t26 - t28 * t92 + t154, t145 * t16 + t40 * t87 + t42 * t92 + t61 * t56 + t156, t10 * t72 + t11 * t71 + t16 * t61 + t26 * t58 + t28 * t57 + t40 * t42, -t14 * t64 - t23 * t51, t14 * t63 - t15 * t64 + t23 * t50 - t24 * t51, -t116 * t64 - t127 * t23, t15 * t63 + t24 * t50, t116 * t63 - t127 * t24, 0, -t116 * t20 - t127 * t4 + t43 * t15 + t27 * t24 - t25 * t50 - t9 * t63, t116 * t21 + t127 * t3 - t43 * t14 - t27 * t23 - t25 * t51 - t9 * t64, t1 * t63 + t14 * t20 - t15 * t21 - t2 * t64 + t23 * t7 - t24 * t8 + t3 * t50 + t4 * t51, -t1 * t21 + t2 * t20 - t25 * t27 - t3 * t8 - t4 * t7 - t43 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t188, t209 * t227, 0, -0.2e1 * t188, 0, 0, t124 * t198, t126 * t198, 0, 0, t67, -t237, 0, t236, 0, 0, -0.2e1 * t115 * t228 + 0.2e1 * t196 * t96, -0.2e1 * t115 * t92 - 0.2e1 * t145 * t196, t144, 0.2e1 * t115 * t196 + 0.2e1 * t184, t67, 0, t237, 0, 0, t236, -0.2e1 * t228 * t61 + 0.2e1 * t40 * t96, t144, 0.2e1 * t145 * t40 + 0.2e1 * t61 * t92, 0.2e1 * t40 * t61 + 0.2e1 * t184, -0.2e1 * t64 * t23, 0.2e1 * t23 * t63 - 0.2e1 * t24 * t64, 0, 0.2e1 * t63 * t24, 0, 0, 0.2e1 * t24 * t43 - 0.2e1 * t25 * t63, -0.2e1 * t23 * t43 - 0.2e1 * t25 * t64, 0.2e1 * t20 * t23 - 0.2e1 * t21 * t24 + 0.2e1 * t3 * t63 + 0.2e1 * t4 * t64, -0.2e1 * t20 * t4 - 0.2e1 * t21 * t3 - 0.2e1 * t25 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, 0, -t150, t116, t45, t44, 0, 0, 0, 0, -t56, 0, t55, t116, t176 * t221 - t163, -t182 - t214 * (pkin(7) * t229 - qJ(4) * t230 - t124 * t152 + t126 * t155 - t186) - 0.2e1 * t175 * t221 (t214 * t55 + t215 * t56) * pkin(3) (t13 * t214 - t163 * t215) * pkin(3), 0, -t56, 0, t116, -t55, 0 (pkin(4) - t114) * t116 - t163, qJD(5) * t86 + t112 * t55 - t114 * t56, t112 * t116 + t130 - 0.2e1 * t201, qJD(5) * t26 + t10 * t112 + t11 * t114, 0, 0, t14, 0, t15, t116, qJD(6) * t134 - t116 * t79 - t69 * t127 + t18 * t183 - t128, t116 * t80 + t127 * t68 - t1, t14 * t79 - t15 * t80 + t50 * t68 + t51 * t69, -t1 * t80 + t2 * t79 - t68 * t8 - t69 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206, 0, -t207, 0, -pkin(8) * t206, pkin(8) * t207, 0, 0, 0, 0, -t92, 0, t228, 0, -t57, -t58 (t214 * t228 + t215 * t92) * pkin(3) (t214 * t58 - t215 * t57) * pkin(3), 0, -t92, 0, 0, -t228, 0, -t57, -qJD(5) * t96 + t112 * t228 - t114 * t92, t58, qJD(5) * t72 + t112 * t58 + t114 * t57, 0, 0, t23, 0, t24, 0, t4, -t3, t23 * t79 - t24 * t80 + t63 * t68 + t64 * t69, -t20 * t69 - t21 * t68 - t3 * t80 - t4 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, t112 * t226, 0, 0, 0, 0, 0, 0, 0.2e1 * t69, -0.2e1 * t68, 0, -0.2e1 * t68 * t80 - 0.2e1 * t69 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t56, 0, t75, 0, 0, 0, 0, 0, 0, -t55, 0, t56, t16, 0, 0, 0, 0, 0, 0, -t15, t14, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t228, -t92, 0, t196, 0, 0, 0, 0, 0, 0, -t228, 0, t92, t40, 0, 0, 0, 0, 0, 0, -t24, t23, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t56, 0, t11, 0, 0, 0, 0, 0, 0, -t224 * t116 - t127 * t204, t116 * t123 - t127 * t183, t224 * t14 - t123 * t15 + (t123 * t51 - t224 * t50) * qJD(6), t2 * t224 - t1 * t123 + (-t123 * t7 + t224 * t8) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, t224 * t23 - t123 * t24 + (t123 * t64 - t224 * t63) * qJD(6), -t4 * t224 - t3 * t123 + (-t123 * t20 + t224 * t21) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, t183, 0, -t69 * t224 - t68 * t123 + (-t123 * t79 + t224 * t80) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, -t15, -t116, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, -t24, 0, -t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t68, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, -t183, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;

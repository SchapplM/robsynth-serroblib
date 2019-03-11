% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:40
% EndTime: 2019-03-09 11:55:51
% DurationCPUTime: 4.12s
% Computational Cost: add. (6949->333), mult. (15227->547), div. (0->0), fcn. (14754->8), ass. (0->163)
t115 = sin(qJ(4));
t117 = cos(qJ(4));
t190 = qJD(4) * t117;
t116 = sin(qJ(2));
t118 = cos(qJ(2));
t197 = sin(pkin(10));
t198 = cos(pkin(10));
t93 = t198 * t116 + t197 * t118;
t178 = t93 * t190;
t92 = t197 * t116 - t198 * t118;
t88 = t92 * qJD(2);
t221 = -t115 * t88 + t178;
t114 = sin(qJ(5));
t212 = cos(qJ(5));
t177 = t212 * t115;
t191 = qJD(4) * t115;
t179 = t93 * t191;
t196 = t114 * t115;
t180 = t93 * t196;
t171 = t212 * qJD(5);
t217 = t212 * qJD(4) + t171;
t22 = -t88 * t177 - t114 * t179 - qJD(5) * t180 + (-t114 * t88 + t217 * t93) * t117;
t96 = t114 * t117 + t177;
t57 = t96 * t93;
t216 = qJD(4) + qJD(5);
t67 = -t217 * t117 + t216 * t196;
t204 = -t96 * t22 + t67 * t57;
t68 = t216 * t96;
t176 = t212 * t117;
t95 = -t176 + t196;
t21 = t68 * t93 - t95 * t88;
t58 = t93 * t176 - t180;
t218 = -t21 * t95 + t58 * t68;
t220 = t218 - t204;
t219 = t218 + t204;
t201 = t67 * t95 - t96 * t68;
t112 = t115 ^ 2;
t113 = t117 ^ 2;
t194 = t112 - t113;
t166 = qJD(4) * t194;
t206 = -qJ(3) - pkin(7);
t170 = qJD(2) * t206;
t134 = -t116 * qJD(3) + t118 * t170;
t86 = t118 * qJD(3) + t116 * t170;
t124 = t197 * t134 + t198 * t86;
t193 = qJD(2) * t116;
t183 = pkin(2) * t193;
t87 = t93 * qJD(2);
t147 = t87 * pkin(3) + t183;
t133 = t88 * pkin(8) + t147;
t110 = -t118 * pkin(2) - pkin(1);
t149 = -t92 * pkin(3) - t110;
t136 = -t93 * pkin(8) - t149;
t61 = t117 * t136;
t98 = t206 * t116;
t99 = t206 * t118;
t70 = t197 * t98 - t198 * t99;
t17 = -qJD(4) * t61 - t115 * t133 - t117 * t124 + t70 * t191;
t123 = t115 * t124;
t66 = t117 * t70;
t42 = t115 * t136 + t66;
t18 = -t42 * qJD(4) + t117 * t133 - t123;
t41 = -t115 * t70 + t61;
t215 = t17 * t115 - t18 * t117 + (t115 * t41 - t117 * t42) * qJD(4);
t214 = -pkin(9) - pkin(8);
t120 = -t123 + t87 * pkin(4) + (-t214 * t88 + t147) * t117 + (-t66 + (-t214 * t93 + t149) * t115) * qJD(4);
t125 = -pkin(9) * t221 - t17;
t122 = t114 * t125 - t212 * t120;
t213 = t92 * pkin(4);
t28 = -t117 * t93 * pkin(9) + t213 + t41;
t199 = t115 * t93;
t32 = -pkin(9) * t199 + t42;
t31 = t212 * t32;
t205 = t114 * t28 + t31;
t4 = -qJD(5) * t205 - t122;
t119 = 2 * qJD(6);
t84 = t87 * pkin(5);
t55 = -t198 * t134 + t197 * t86;
t69 = -t197 * t99 - t198 * t98;
t210 = t69 * t55;
t209 = t93 * t88;
t208 = t95 * t68;
t195 = t115 * t117;
t192 = qJD(2) * t118;
t189 = qJD(5) * t114;
t188 = 0.2e1 * t57 * t22;
t65 = 0.2e1 * t92 * t87;
t187 = 0.2e1 * t208;
t186 = -0.2e1 * pkin(1) * qJD(2);
t107 = -t198 * pkin(2) - pkin(3);
t185 = 0.2e1 * qJD(4) * t107;
t184 = t212 * pkin(4);
t182 = pkin(4) * t191;
t181 = pkin(4) * t189;
t175 = t115 * t190;
t174 = t116 * t192;
t164 = t197 * pkin(2) + pkin(8);
t150 = pkin(9) + t164;
t132 = t150 * t212;
t126 = qJD(4) * t132;
t129 = t115 * t132;
t138 = t114 * t150;
t131 = qJD(4) * t138;
t90 = t150 * t117;
t35 = qJD(5) * t129 + t115 * t126 + t117 * t131 + t90 * t189;
t135 = t115 * t138;
t36 = -qJD(5) * t135 - t115 * t131 + t117 * t126 + t90 * t171;
t63 = t114 * t90 + t129;
t64 = t212 * t90 - t135;
t173 = -t64 * t35 + t63 * t36;
t167 = -0.4e1 * t93 * t195;
t91 = t93 ^ 2;
t165 = t91 * t175;
t163 = t21 * t57 - t58 * t22;
t161 = t22 * t95 + t57 * t68;
t160 = t92 * t22 + t87 * t57;
t159 = -t35 * t92 + t64 * t87;
t158 = -t36 * t92 - t63 * t87;
t157 = t55 * t93 - t69 * t88;
t156 = t67 * t92 - t96 * t87;
t155 = t92 * t68 + t87 * t95;
t154 = t93 * t87 - t88 * t92;
t152 = -t42 * t115 - t41 * t117;
t52 = t96 * t67;
t148 = -0.2e1 * t52 + 0.2e1 * t208;
t48 = pkin(4) * t199 + t69;
t146 = t115 * t164;
t145 = t117 * t164;
t13 = -t114 * t32 + t212 * t28;
t143 = t115 * t87 + t92 * t190;
t3 = -t114 * t120 - t212 * t125 - t28 * t171 + t32 * t189;
t141 = qJD(4) * t164;
t140 = -t63 * t21 - t64 * t22 + t35 * t57 + t36 * t58;
t139 = -t35 * t96 + t36 * t95 + t63 * t68 - t64 * t67;
t34 = t221 * pkin(4) + t55;
t97 = -t117 * pkin(4) + t107;
t137 = -t68 * pkin(5) - t67 * qJ(6) + t96 * qJD(6);
t82 = t87 * qJ(6);
t89 = t92 * qJD(6);
t1 = -t3 + t82 + t89;
t130 = 0.2e1 * t35 * t95 + 0.2e1 * t36 * t96 - 0.2e1 * t63 * t67 - 0.2e1 * t64 * t68;
t128 = -t107 * t88 - t164 * t87;
t127 = -t107 * t93 + t164 * t92;
t121 = (-t31 + (-t28 - t213) * t114) * qJD(5) - t122;
t111 = pkin(4) * t171;
t109 = -t184 - pkin(5);
t106 = t114 * pkin(4) + qJ(6);
t104 = -0.2e1 * t181;
t100 = t111 + qJD(6);
t62 = t117 * t87 - t92 * t191;
t59 = t95 * pkin(5) - t96 * qJ(6) + t97;
t49 = -0.2e1 * t52;
t45 = t93 * t166 + t88 * t195;
t33 = -t137 + t182;
t23 = t57 * pkin(5) - t58 * qJ(6) + t48;
t15 = -0.2e1 * t58 * t21;
t12 = -0.2e1 * t21 * t92 + 0.2e1 * t58 * t87;
t10 = -t21 * t96 - t58 * t67;
t9 = -t92 * pkin(5) - t13;
t8 = t92 * qJ(6) + t205;
t5 = t22 * pkin(5) + t21 * qJ(6) - t58 * qJD(6) + t34;
t2 = -t4 - t84;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t174, 0.2e1 * (-t116 ^ 2 + t118 ^ 2) * qJD(2), 0, -0.2e1 * t174, 0, 0, t116 * t186, t118 * t186, 0, 0, -0.2e1 * t209, -0.2e1 * t154, 0, t65, 0, 0, 0.2e1 * t110 * t87 + 0.2e1 * t92 * t183, -0.2e1 * t110 * t88 + 0.2e1 * t93 * t183, -0.2e1 * t124 * t92 - 0.2e1 * t70 * t87 + 0.2e1 * t157, 0.2e1 * t110 * t183 + 0.2e1 * t124 * t70 + 0.2e1 * t210, -0.2e1 * t113 * t209 - 0.2e1 * t165, 0.2e1 * t166 * t91 - t167 * t88, 0.2e1 * t117 * t154 - 0.2e1 * t179 * t92, -0.2e1 * t112 * t209 + 0.2e1 * t165, -0.2e1 * t115 * t154 - 0.2e1 * t178 * t92, t65, 0.2e1 * t115 * t157 + 0.2e1 * t178 * t69 + 0.2e1 * t18 * t92 + 0.2e1 * t41 * t87, 0.2e1 * t117 * t157 + 0.2e1 * t17 * t92 - 0.2e1 * t179 * t69 - 0.2e1 * t42 * t87, -0.2e1 * t152 * t88 + 0.2e1 * t215 * t93, -0.2e1 * t42 * t17 + 0.2e1 * t41 * t18 + 0.2e1 * t210, t15, 0.2e1 * t163, t12, t188, -0.2e1 * t160, t65, 0.2e1 * t13 * t87 + 0.2e1 * t48 * t22 + 0.2e1 * t34 * t57 + 0.2e1 * t4 * t92, -0.2e1 * t205 * t87 - 0.2e1 * t48 * t21 + 0.2e1 * t3 * t92 + 0.2e1 * t34 * t58, 0.2e1 * t13 * t21 - 0.2e1 * t205 * t22 + 0.2e1 * t3 * t57 - 0.2e1 * t4 * t58, 0.2e1 * t13 * t4 - 0.2e1 * t205 * t3 + 0.2e1 * t48 * t34, t15, t12, -0.2e1 * t163, t65, 0.2e1 * t160, t188, -0.2e1 * t2 * t92 + 0.2e1 * t23 * t22 + 0.2e1 * t5 * t57 - 0.2e1 * t9 * t87, -0.2e1 * t1 * t57 + 0.2e1 * t2 * t58 - 0.2e1 * t9 * t21 - 0.2e1 * t8 * t22, 0.2e1 * t1 * t92 + 0.2e1 * t23 * t21 - 0.2e1 * t5 * t58 + 0.2e1 * t8 * t87, 0.2e1 * t8 * t1 + 0.2e1 * t9 * t2 + 0.2e1 * t23 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, 0, -t193, 0, -pkin(7) * t192, pkin(7) * t193, 0, 0, 0, 0, -t88, 0, -t87, 0, -t55, -t124 (-t197 * t87 + t198 * t88) * pkin(2) (t124 * t197 - t55 * t198) * pkin(2), -t45, qJD(4) * t167 + t194 * t88, t143, t45, t62, 0, -t55 * t117 + t128 * t115 + (t69 * t115 - t117 * t127) * qJD(4), t55 * t115 + t128 * t117 + (t115 * t127 + t69 * t117) * qJD(4), qJD(4) * t152 - t18 * t115 - t17 * t117, -t17 * t145 - t18 * t146 + t55 * t107 + (-t145 * t41 - t146 * t42) * qJD(4), t10, -t220, -t156, t161, -t155, 0, t57 * t182 + t97 * t22 + t34 * t95 + t48 * t68 + t158, t58 * t182 - t97 * t21 + t34 * t96 - t48 * t67 - t159, t13 * t67 - t205 * t68 + t3 * t95 - t4 * t96 + t140, -t13 * t36 + t48 * t182 - t205 * t35 - t3 * t64 + t34 * t97 - t4 * t63, t10, -t156, t220, 0, t155, t161, t59 * t22 + t23 * t68 + t33 * t57 + t5 * t95 + t158, -t1 * t95 + t2 * t96 - t9 * t67 - t8 * t68 + t140, t59 * t21 + t23 * t67 - t33 * t58 - t5 * t96 + t159, t1 * t64 + t2 * t63 + t23 * t33 - t8 * t35 + t9 * t36 + t5 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t175, -0.2e1 * t166, 0, -0.2e1 * t175, 0, 0, t115 * t185, t117 * t185, 0, 0, t49, 0.2e1 * t201, 0, t187, 0, 0, 0.2e1 * t95 * t182 + 0.2e1 * t97 * t68, 0.2e1 * t96 * t182 - 0.2e1 * t97 * t67, t130, 0.2e1 * t97 * t182 + 0.2e1 * t173, t49, 0, -0.2e1 * t201, 0, 0, t187, 0.2e1 * t33 * t95 + 0.2e1 * t59 * t68, t130, -0.2e1 * t33 * t96 + 0.2e1 * t59 * t67, 0.2e1 * t59 * t33 + 0.2e1 * t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t88, 0, t183, 0, 0, 0, 0, 0, 0, t62, -t143 -(-t112 - t113) * t88, -t215, 0, 0, 0, 0, 0, 0, -t155, t156, t219, -t13 * t68 - t205 * t67 - t3 * t96 - t4 * t95, 0, 0, 0, 0, 0, 0, -t155, t219, -t156, t1 * t96 + t2 * t95 - t8 * t67 + t9 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117 * t88 - t179, 0, -t221, t87, t18, t17, 0, 0, 0, 0, -t21, 0, -t22, t87, t87 * t184 + t121 (-t114 * t87 - t171 * t92) * pkin(4) + t3 (t212 * t21 - t114 * t22 + (t114 * t58 - t212 * t57) * qJD(5)) * pkin(4) (t212 * t4 - t114 * t3 + (-t114 * t13 + t205 * t212) * qJD(5)) * pkin(4), 0, -t21, 0, t87, t22, 0, -t109 * t87 + t121 + t84, -t100 * t57 - t106 * t22 - t109 * t21 + t58 * t181, t100 * t92 + t106 * t87 + t1, t1 * t106 + t8 * t100 + t2 * t109 + t181 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, 0, -t191, 0, -t117 * t141, t115 * t141, 0, 0, 0, 0, -t67, 0, -t68, 0, -t36, t35 (t212 * t67 - t114 * t68 + (t114 * t96 - t212 * t95) * qJD(5)) * pkin(4) (-t212 * t36 - t114 * t35 + (t114 * t63 + t212 * t64) * qJD(5)) * pkin(4), 0, -t67, 0, 0, t68, 0, -t36, -t100 * t95 - t106 * t68 - t109 * t67 + t96 * t181, -t35, t64 * t100 - t35 * t106 + t36 * t109 + t181 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, -t190, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t67, 0 (-t212 * t68 - t114 * t67 + (t114 * t95 + t212 * t96) * qJD(5)) * pkin(4), 0, 0, 0, 0, 0, 0, -t68, 0, -t67, t96 * t100 - t67 * t106 + t68 * t109 + t181 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, -0.2e1 * t111, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, 0.2e1 * t100, 0.2e1 * t106 * t100 + 0.2e1 * t109 * t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, -t22, t87, t4, t3, 0, 0, 0, -t21, 0, t87, t22, 0, t4 + 0.2e1 * t84, pkin(5) * t21 - t22 * qJ(6) - t57 * qJD(6), -t3 + 0.2e1 * t82 + 0.2e1 * t89, -t2 * pkin(5) + t1 * qJ(6) + t8 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, 0, -t68, 0, -t36, t35, 0, 0, 0, -t67, 0, 0, t68, 0, -t36, pkin(5) * t67 - t68 * qJ(6) - t95 * qJD(6), -t35, -t36 * pkin(5) - t35 * qJ(6) + t64 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t67, 0, 0, 0, 0, 0, 0, 0, 0, -t68, 0, -t67, t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, -t111, 0, 0, 0, 0, 0, 0, 0, 0, -t181, 0, t119 + t111, -pkin(5) * t181 + t100 * qJ(6) + t106 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, qJ(6) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t21, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;

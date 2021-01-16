% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:35:53
% EndTime: 2021-01-16 00:36:05
% DurationCPUTime: 2.78s
% Computational Cost: add. (2841->300), mult. (7903->561), div. (0->0), fcn. (7189->8), ass. (0->153)
t85 = sin(qJ(2));
t88 = cos(qJ(2));
t116 = pkin(2) * t85 - pkin(8) * t88;
t82 = sin(pkin(5));
t155 = qJD(2) * t82;
t101 = t116 * t155;
t87 = cos(qJ(3));
t149 = qJD(3) * t87;
t84 = sin(qJ(3));
t151 = qJD(3) * t84;
t156 = cos(pkin(5));
t123 = pkin(1) * t156;
t167 = t82 * t88;
t97 = pkin(7) * t167 + t123 * t85;
t48 = pkin(8) * t156 + t97;
t49 = (-pkin(2) * t88 - pkin(8) * t85 - pkin(1)) * t82;
t168 = t82 * t85;
t96 = -pkin(7) * t168 + t123 * t88;
t95 = qJD(2) * t96;
t12 = -t101 * t84 - t149 * t49 + t151 * t48 - t87 * t95;
t154 = qJD(2) * t85;
t130 = t82 * t154;
t47 = -pkin(2) * t156 - t96;
t53 = -t156 * t87 + t84 * t168;
t54 = t156 * t84 + t87 * t168;
t93 = t53 * pkin(3) - t54 * pkin(9) + t47;
t187 = -pkin(9) * t130 - qJD(4) * t93 + t12;
t186 = -0.4e1 * t84;
t153 = qJD(2) * t88;
t129 = t82 * t153;
t33 = -qJD(3) * t53 + t129 * t87;
t185 = -qJD(4) * t167 + t33;
t184 = t84 * qJD(5) + (pkin(8) * qJD(4) + qJ(5) * qJD(3)) * t87;
t183 = 0.2e1 * t82;
t182 = 0.2e1 * qJD(4);
t181 = pkin(8) * t82;
t83 = sin(qJ(4));
t180 = pkin(8) * t83;
t32 = qJD(3) * t54 + t129 * t84;
t179 = t32 * pkin(4);
t86 = cos(qJ(4));
t34 = t86 * t167 + t54 * t83;
t178 = t34 * pkin(4);
t141 = -t48 * t149 - t49 * t151 - t84 * t95;
t11 = (-t85 * pkin(3) - t116 * t87) * t155 - t141;
t76 = qJD(4) * t86;
t19 = -t130 * t86 + t185 * t83 + t54 * t76;
t5 = t19 * pkin(4) + t11;
t177 = t5 * t83;
t176 = t5 * t86;
t175 = t86 * pkin(4);
t147 = qJD(4) * t83;
t20 = t130 * t83 - t147 * t54 + t185 * t86;
t174 = t20 * t83;
t35 = -t83 * t167 + t54 * t86;
t173 = t35 * t83;
t139 = pkin(8) * t149;
t56 = t149 * t83 + t76 * t84;
t42 = pkin(4) * t56 + t139;
t172 = t42 * t83;
t171 = t42 * t86;
t163 = qJ(5) + pkin(9);
t66 = t163 * t83;
t170 = t66 * t84;
t67 = t163 * t86;
t169 = t67 * t84;
t166 = t83 * t87;
t165 = t84 * t86;
t164 = t86 * t87;
t162 = t87 * t48 + t84 * t49;
t27 = -pkin(9) * t167 + t162;
t9 = t86 * t27 + t83 * t93;
t115 = -t87 * pkin(3) - t84 * pkin(9);
t109 = -pkin(2) + t115;
t102 = qJD(4) * t109;
t114 = pkin(3) * t84 - pkin(9) * t87;
t103 = t114 * qJD(3);
t161 = -t86 * t102 - t83 * t103;
t128 = t83 * t151;
t160 = pkin(8) * t128 + t86 * t103;
t72 = pkin(8) * t164;
t45 = t83 * t109 + t72;
t78 = t83 ^ 2;
t80 = t86 ^ 2;
t159 = t78 - t80;
t79 = t84 ^ 2;
t158 = -t87 ^ 2 + t79;
t157 = qJ(5) * t84;
t152 = qJD(3) * t83;
t150 = qJD(3) * t86;
t148 = qJD(3) * t88;
t146 = qJD(4) * t87;
t144 = pkin(8) * t166;
t143 = -0.2e1 * pkin(2) * qJD(3);
t142 = -0.2e1 * pkin(3) * qJD(4);
t140 = pkin(4) * t151;
t138 = pkin(4) * t147;
t77 = t82 ^ 2;
t137 = t77 * t153;
t135 = t84 * t147;
t134 = t83 * t146;
t133 = t86 * t146;
t122 = -t84 * t48 + t87 * t49;
t26 = pkin(3) * t167 - t122;
t14 = t26 + t178;
t132 = t14 * t147;
t62 = (pkin(4) * t83 + pkin(8)) * t84;
t131 = t62 * t147;
t127 = t83 * t76;
t126 = t84 * t149;
t125 = t86 * t149;
t74 = -pkin(3) - t175;
t124 = -t74 + t175;
t8 = -t83 * t27 + t86 * t93;
t121 = t159 * qJD(4);
t120 = t158 * qJD(3);
t119 = 0.2e1 * t126;
t118 = t85 * t137;
t117 = t83 * t125;
t6 = t53 * pkin(4) - t35 * qJ(5) + t8;
t7 = -t34 * qJ(5) + t9;
t113 = -t6 * t86 - t7 * t83;
t112 = t156 * t155;
t111 = pkin(4) * t78 + t74 * t86;
t110 = -t34 * t86 - t173;
t107 = t11 * t83 + t26 * t76;
t106 = -t11 * t86 + t147 * t26;
t105 = t83 * t32 + t53 * t76;
t104 = t147 * t53 - t86 * t32;
t50 = t97 * qJD(2);
t94 = t32 * pkin(3) - t33 * pkin(9) + t50;
t3 = t147 * t27 + t187 * t86 - t83 * t94;
t100 = t148 * t84 + t154 * t87;
t99 = -t148 * t87 + t154 * t84;
t55 = t125 - t135;
t98 = t150 * t84 + t134;
t91 = qJ(5) * t135 - t83 * t102 - t184 * t86 + t160;
t4 = t187 * t83 - t27 * t76 + t86 * t94;
t89 = -t20 * qJ(5) - t35 * qJD(5) + t4;
t61 = t86 * t109;
t52 = -t83 * qJD(5) - t163 * t76;
t51 = -t86 * qJD(5) + t163 * t147;
t44 = t61 - t144;
t36 = -t157 * t83 + t45;
t30 = -t86 * t157 + t61 + (-pkin(4) - t180) * t87;
t29 = -qJD(4) * t45 + t160;
t28 = t98 * pkin(8) + t161;
t21 = (pkin(8) * qJD(3) + qJ(5) * qJD(4)) * t165 + t184 * t83 + t161;
t17 = t91 + t140;
t13 = t101 * t87 + t141;
t2 = t19 * qJ(5) + t34 * qJD(5) + t3;
t1 = t89 + t179;
t10 = [0, 0, 0, 0.2e1 * t118, 0.2e1 * (-t85 ^ 2 + t88 ^ 2) * t77 * qJD(2), 0.2e1 * t88 * t112, -0.2e1 * t85 * t112, 0, -0.2e1 * pkin(1) * t154 * t77 - 0.2e1 * t156 * t50, -0.2e1 * pkin(1) * t137 - 0.2e1 * t156 * t95, 0.2e1 * t54 * t33, -0.2e1 * t54 * t32 - 0.2e1 * t33 * t53, (t154 * t54 - t33 * t88) * t183, (-t154 * t53 + t32 * t88) * t183, -0.2e1 * t118, 0.2e1 * t47 * t32 + 0.2e1 * t50 * t53 + 0.2e1 * (t122 * t154 - t13 * t88) * t82, 0.2e1 * t47 * t33 + 0.2e1 * t50 * t54 + 0.2e1 * (-t12 * t88 - t162 * t154) * t82, 0.2e1 * t35 * t20, -0.2e1 * t35 * t19 - 0.2e1 * t20 * t34, 0.2e1 * t20 * t53 + 0.2e1 * t35 * t32, -0.2e1 * t19 * t53 - 0.2e1 * t34 * t32, 0.2e1 * t53 * t32, 0.2e1 * t11 * t34 + 0.2e1 * t26 * t19 + 0.2e1 * t8 * t32 + 0.2e1 * t4 * t53, 0.2e1 * t11 * t35 + 0.2e1 * t26 * t20 + 0.2e1 * t3 * t53 - 0.2e1 * t9 * t32, 0.2e1 * t1 * t53 + 0.2e1 * t14 * t19 + 0.2e1 * t6 * t32 + 0.2e1 * t5 * t34, 0.2e1 * t14 * t20 + 0.2e1 * t2 * t53 - 0.2e1 * t7 * t32 + 0.2e1 * t5 * t35, -0.2e1 * t1 * t35 - 0.2e1 * t7 * t19 + 0.2e1 * t2 * t34 - 0.2e1 * t6 * t20, 0.2e1 * t6 * t1 + 0.2e1 * t14 * t5 - 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, t129, -t130, 0, -t50, -t95, t149 * t54 + t33 * t84, -t84 * t32 + t33 * t87 + (-t53 * t87 - t54 * t84) * qJD(3), t99 * t82, t100 * t82, 0, -pkin(2) * t32 + t151 * t47 - t99 * t181 - t50 * t87, -pkin(2) * t33 - t100 * t181 + t149 * t47 + t50 * t84, t20 * t165 + t35 * t55, t110 * t149 + (-t19 * t86 - t174 + (t34 * t83 - t35 * t86) * qJD(4)) * t84, (t150 * t53 - t20) * t87 + (qJD(3) * t35 - t104) * t84, (-t152 * t53 + t19) * t87 + (-qJD(3) * t34 - t105) * t84, t151 * t53 - t32 * t87, t29 * t53 + t44 * t32 + (-t4 + (pkin(8) * t34 + t26 * t83) * qJD(3)) * t87 + (pkin(8) * t19 + qJD(3) * t8 + t107) * t84, t28 * t53 - t45 * t32 + (-t3 + (pkin(8) * t35 + t26 * t86) * qJD(3)) * t87 + (pkin(8) * t20 - qJD(3) * t9 - t106) * t84, t17 * t53 + t62 * t19 + t30 * t32 + t42 * t34 + (t14 * t152 - t1) * t87 + (qJD(3) * t6 + t14 * t76 + t177) * t84, t62 * t20 + t21 * t53 - t36 * t32 + t42 * t35 + (t14 * t150 - t2) * t87 + (-qJD(3) * t7 - t132 + t176) * t84, -t17 * t35 - t36 * t19 - t30 * t20 + t21 * t34 + t113 * t149 + (-t1 * t86 + t2 * t83 + (t6 * t83 - t7 * t86) * qJD(4)) * t84, t1 * t30 + t14 * t42 + t6 * t17 - t2 * t36 - t7 * t21 + t5 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -0.2e1 * t120, 0, 0, 0, t84 * t143, t87 * t143, 0.2e1 * t126 * t80 - 0.2e1 * t127 * t79, t159 * t79 * t182 + t117 * t186, 0.2e1 * t134 * t84 + 0.2e1 * t150 * t158, -0.2e1 * t120 * t83 + 0.2e1 * t133 * t84, -0.2e1 * t126, 0.2e1 * t44 * t151 - 0.2e1 * t29 * t87 + 0.2e1 * (t119 * t83 + t76 * t79) * pkin(8), -0.2e1 * t45 * t151 - 0.2e1 * t28 * t87 + 0.2e1 * (t119 * t86 - t147 * t79) * pkin(8), 0.2e1 * (t152 * t62 - t17) * t87 + 0.2e1 * (qJD(3) * t30 + t62 * t76 + t172) * t84, 0.2e1 * (t150 * t62 - t21) * t87 + 0.2e1 * (-qJD(3) * t36 - t131 + t171) * t84, 0.2e1 * (-t30 * t86 - t36 * t83) * t149 + 0.2e1 * (-t17 * t86 + t21 * t83 + (t30 * t83 - t36 * t86) * qJD(4)) * t84, 0.2e1 * t30 * t17 - 0.2e1 * t36 * t21 + 0.2e1 * t62 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, t130, t13, t12, t35 * t76 + t174, qJD(4) * t110 - t83 * t19 + t20 * t86, t105, -t104, 0, -pkin(3) * t19 - pkin(9) * t105 + t106, -pkin(3) * t20 + pkin(9) * t104 + t107, t74 * t19 - t66 * t32 - t176 + t52 * t53 + (t14 + t178) * t147, t74 * t20 - t67 * t32 + t177 + t51 * t53 + (pkin(4) * t173 + t14 * t86) * qJD(4), qJD(4) * t113 - t1 * t83 - t67 * t19 - t2 * t86 + t66 * t20 + t51 * t34 - t52 * t35, pkin(4) * t132 - t1 * t66 - t2 * t67 + t5 * t74 - t7 * t51 + t6 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, -t151, 0, -t139, pkin(8) * t151, -t121 * t84 + t117, t127 * t186 - t149 * t159, t128 - t133, t98, 0, (pkin(9) * t164 + (-pkin(3) * t86 + t180) * t84) * qJD(4) + (t115 * t83 - t72) * qJD(3), (pkin(8) * t165 + t114 * t83) * qJD(4) + (t115 * t86 + t144) * qJD(3), -t171 - t52 * t87 + (t74 * t166 - t170) * qJD(3) + (t111 * t84 + t62 * t83) * qJD(4), t172 - t51 * t87 + (t74 * t164 - t169) * qJD(3) + (t124 * t83 * t84 + t62 * t86) * qJD(4), (t66 * t149 - t52 * t84 - t21 + (-t30 - t169) * qJD(4)) * t86 + (-t67 * t149 + t51 * t84 - t17 + (-t36 - t170) * qJD(4)) * t83, pkin(4) * t131 - t17 * t66 - t21 * t67 + t30 * t52 - t36 * t51 + t42 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t127, -0.2e1 * t121, 0, 0, 0, t83 * t142, t86 * t142, -0.2e1 * t124 * t147, t111 * t182, -0.2e1 * t51 * t86 - 0.2e1 * t52 * t83 + 0.2e1 * (t66 * t86 - t67 * t83) * qJD(4), 0.2e1 * t138 * t74 - 0.2e1 * t67 * t51 - 0.2e1 * t66 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t32, t4, t3, t89 + 0.2e1 * t179, t2, -t20 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t56, t151, t29, t28, t91 + 0.2e1 * t140, t21, -t55 * pkin(4), t17 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t147, 0, -pkin(9) * t76, pkin(9) * t147, t52, t51, -pkin(4) * t76, t52 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t55, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t76, 0, t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;

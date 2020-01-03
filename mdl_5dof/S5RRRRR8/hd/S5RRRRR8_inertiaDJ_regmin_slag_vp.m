% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:47
% EndTime: 2019-12-31 22:25:52
% DurationCPUTime: 1.30s
% Computational Cost: add. (2007->194), mult. (4803->332), div. (0->0), fcn. (4519->8), ass. (0->150)
t103 = sin(qJ(4));
t106 = cos(qJ(5));
t102 = sin(qJ(5));
t107 = cos(qJ(4));
t156 = t102 * t107;
t77 = t106 * t103 + t156;
t104 = sin(qJ(3));
t105 = sin(qJ(2));
t108 = cos(qJ(3));
t109 = cos(qJ(2));
t78 = t104 * t109 + t108 * t105;
t43 = t77 * t78;
t101 = t107 ^ 2;
t154 = t103 ^ 2 - t101;
t128 = t154 * qJD(4);
t182 = qJD(2) + qJD(3);
t181 = qJD(4) + qJD(5);
t180 = pkin(6) + pkin(7);
t179 = -pkin(8) - pkin(9);
t76 = t104 * t105 - t108 * t109;
t178 = pkin(4) * t76;
t55 = t182 * t78;
t177 = t55 * pkin(4);
t94 = t104 * pkin(2) + pkin(8);
t176 = -pkin(9) - t94;
t175 = pkin(4) * t106;
t98 = qJD(4) * t107;
t137 = t78 * t98;
t54 = t182 * t76;
t161 = t103 * t54;
t113 = t137 - t161;
t151 = qJD(4) * t103;
t149 = t105 * qJD(2);
t146 = pkin(2) * t149;
t28 = pkin(3) * t55 + pkin(8) * t54 + t146;
t134 = qJD(2) * t180;
t124 = t109 * t134;
t125 = t105 * t134;
t152 = qJD(3) * t108;
t153 = qJD(3) * t104;
t87 = t180 * t105;
t89 = t180 * t109;
t33 = t104 * t124 + t108 * t125 + t87 * t152 + t89 * t153;
t97 = -t109 * pkin(2) - pkin(1);
t49 = t76 * pkin(3) - t78 * pkin(8) + t97;
t58 = -t104 * t87 + t108 * t89;
t8 = -t103 * t28 + t107 * t33 + t58 * t151 - t49 * t98;
t6 = -t113 * pkin(9) - t8;
t174 = t106 * t6;
t173 = t108 * pkin(2);
t172 = t78 * t54;
t34 = t58 * qJD(3) - t104 * t125 + t108 * t124;
t14 = t113 * pkin(4) + t34;
t160 = t103 * t78;
t57 = t104 * t89 + t108 * t87;
t39 = pkin(4) * t160 + t57;
t53 = t181 * t77;
t75 = t102 * t103 - t106 * t107;
t171 = t14 * t75 + t39 * t53;
t52 = t181 * t75;
t170 = t14 * t77 - t39 * t52;
t51 = t57 * t98;
t169 = t34 * t103 + t51;
t56 = t107 * t58;
t168 = t103 * t49 + t56;
t141 = pkin(4) * t151;
t143 = pkin(2) * t153;
t82 = t141 + t143;
t96 = -t107 * pkin(4) - pkin(3);
t85 = t96 - t173;
t167 = t85 * t53 + t82 * t75;
t166 = -t85 * t52 + t82 * t77;
t165 = t75 * t141 + t96 * t53;
t164 = t77 * t141 - t96 * t52;
t95 = -pkin(3) - t173;
t163 = t103 * t143 + t95 * t98;
t22 = -pkin(9) * t160 + t168;
t162 = t102 * t22;
t159 = t106 * t22;
t158 = t107 * t54;
t157 = t107 * t78;
t155 = t103 * t107;
t150 = qJD(5) * t102;
t148 = t109 * qJD(2);
t147 = -0.2e1 * pkin(1) * qJD(2);
t145 = pkin(3) * t151;
t144 = pkin(3) * t98;
t142 = pkin(2) * t152;
t140 = pkin(4) * t150;
t139 = qJD(5) * t175;
t138 = t78 * t151;
t50 = t57 * t151;
t46 = t107 * t49;
t15 = -pkin(9) * t157 - t103 * t58 + t178 + t46;
t136 = -t15 - t178;
t130 = t103 * t33 + t107 * t28;
t5 = pkin(9) * t158 + t177 + (-t56 + (pkin(9) * t78 - t49) * t103) * qJD(4) + t130;
t135 = -t102 * t6 + t106 * t5;
t133 = qJD(4) * t179;
t132 = t103 * t98;
t131 = qJD(4) * t176;
t129 = -0.4e1 * t78 * t155;
t127 = t103 * t142;
t126 = t107 * t142;
t123 = t34 * t78 - t54 * t57;
t122 = -t54 * t76 + t55 * t78;
t121 = t76 * t94 - t78 * t95;
t120 = -t106 * t15 + t162;
t119 = t102 * t15 + t159;
t70 = t176 * t103;
t99 = t107 * pkin(9);
t71 = t107 * t94 + t99;
t118 = t102 * t71 - t106 * t70;
t117 = t102 * t70 + t106 * t71;
t86 = t179 * t103;
t88 = t107 * pkin(8) + t99;
t116 = t102 * t88 - t106 * t86;
t115 = t102 * t86 + t106 * t88;
t114 = -t107 * t143 + t95 * t151;
t112 = t138 + t158;
t111 = -t107 * t55 + t76 * t151;
t110 = -t54 * t95 - t55 * t94 + (t104 * t78 - t108 * t76) * qJD(3) * pkin(2);
t91 = 0.2e1 * t132;
t81 = t107 * t133;
t80 = t103 * t133;
t74 = -0.2e1 * t128;
t73 = t78 ^ 2;
t62 = t107 * t131 - t127;
t61 = t103 * t131 + t126;
t44 = t75 * t78;
t38 = -0.2e1 * t77 * t52;
t37 = 0.2e1 * t76 * t55;
t36 = t103 * t55 + t76 * t98;
t32 = -t115 * qJD(5) - t102 * t80 + t106 * t81;
t31 = t116 * qJD(5) - t102 * t81 - t106 * t80;
t27 = -t78 * t128 - t54 * t155;
t21 = -t53 * t76 - t55 * t75;
t20 = -t52 * t76 + t55 * t77;
t19 = qJD(4) * t129 + t154 * t54;
t18 = 0.2e1 * t52 * t75 - 0.2e1 * t53 * t77;
t17 = -t117 * qJD(5) - t102 * t61 + t106 * t62;
t16 = t118 * qJD(5) - t102 * t62 - t106 * t61;
t11 = -t54 * t156 - t102 * t138 - t150 * t160 + (t181 * t157 - t161) * t106;
t10 = -t181 * t43 + t75 * t54;
t9 = -t168 * qJD(4) + t130;
t7 = t10 * t77 + t44 * t52;
t3 = -t10 * t75 - t11 * t77 + t43 * t52 + t44 * t53;
t2 = -t119 * qJD(5) + t135;
t1 = t120 * qJD(5) - t102 * t5 - t174;
t4 = [0, 0, 0, 0.2e1 * t105 * t148, 0.2e1 * (-t105 ^ 2 + t109 ^ 2) * qJD(2), 0, 0, 0, t105 * t147, t109 * t147, -0.2e1 * t172, -0.2e1 * t122, 0, 0, 0, 0.2e1 * t76 * t146 + 0.2e1 * t55 * t97, 0.2e1 * t78 * t146 - 0.2e1 * t54 * t97, -0.2e1 * t101 * t172 - 0.2e1 * t73 * t132, 0.2e1 * t73 * t128 - t54 * t129, 0.2e1 * t122 * t107 - 0.2e1 * t76 * t138, -0.2e1 * t122 * t103 - 0.2e1 * t76 * t137, t37, 0.2e1 * t78 * t51 + 0.2e1 * t46 * t55 + 0.2e1 * t9 * t76 + 0.2e1 * (-t55 * t58 + t123) * t103, 0.2e1 * t123 * t107 - 0.2e1 * t168 * t55 - 0.2e1 * t78 * t50 + 0.2e1 * t8 * t76, -0.2e1 * t44 * t10, -0.2e1 * t10 * t43 + 0.2e1 * t11 * t44, 0.2e1 * t10 * t76 - 0.2e1 * t44 * t55, -0.2e1 * t11 * t76 - 0.2e1 * t43 * t55, t37, 0.2e1 * t39 * t11 - 0.2e1 * t120 * t55 + 0.2e1 * t14 * t43 + 0.2e1 * t2 * t76, 0.2e1 * t1 * t76 + 0.2e1 * t39 * t10 - 0.2e1 * t119 * t55 - 0.2e1 * t14 * t44; 0, 0, 0, 0, 0, t148, -t149, 0, -pkin(6) * t148, pkin(6) * t149, 0, 0, -t54, -t55, 0, -t34, t33, t27, t19, t36, -t111, 0, t50 + (-t121 * qJD(4) - t34) * t107 + t110 * t103, t110 * t107 + t121 * t151 + t169, t7, t3, t20, t21, 0, t85 * t11 - t118 * t55 + t17 * t76 + t82 * t43 + t171, t85 * t10 - t117 * t55 + t16 * t76 - t82 * t44 + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t143, -0.2e1 * t142, t91, t74, 0, 0, 0, 0.2e1 * t114, 0.2e1 * t163, t38, t18, 0, 0, 0, 0.2e1 * t167, 0.2e1 * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t55, 0, -t34, t33, t27, t19, t36, -t111, 0, t50 + (pkin(3) * t54 - pkin(8) * t55) * t103 + (-t34 + (-pkin(3) * t78 - pkin(8) * t76) * qJD(4)) * t107, t112 * pkin(3) + t111 * pkin(8) + t169, t7, t3, t20, t21, 0, t96 * t11 - t116 * t55 + t141 * t43 + t32 * t76 + t171, t96 * t10 - t115 * t55 - t141 * t44 + t31 * t76 + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, -t142, t91, t74, 0, 0, 0, t114 - t145, -t144 + t163, t38, t18, 0, 0, 0, t165 + t167, t164 + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t74, 0, 0, 0, -0.2e1 * t145, -0.2e1 * t144, t38, t18, 0, 0, 0, 0.2e1 * t165, 0.2e1 * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, -t113, t55, t9, t8, 0, 0, t10, -t11, t55, t55 * t175 + (t102 * t136 - t159) * qJD(5) + t135, -t174 + (-t5 - t177) * t102 + (t106 * t136 + t162) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, -t151, 0, -t94 * t98 - t127, t94 * t151 - t126, 0, 0, -t52, -t53, 0, t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, -t151, 0, -pkin(8) * t98, pkin(8) * t151, 0, 0, -t52, -t53, 0, t32, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t140, -0.2e1 * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t55, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t53, 0, t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t53, 0, t32, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, -t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;

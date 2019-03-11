% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:20:28
% EndTime: 2019-03-08 23:20:33
% DurationCPUTime: 1.88s
% Computational Cost: add. (2399->240), mult. (6450->455), div. (0->0), fcn. (6114->12), ass. (0->145)
t118 = sin(qJ(3));
t117 = sin(qJ(4));
t122 = cos(qJ(3));
t156 = t122 * qJD(3);
t151 = t117 * t156;
t121 = cos(qJ(4));
t159 = qJD(4) * t121;
t178 = t118 * t159 + t151;
t166 = t121 * t122;
t101 = pkin(8) * t166;
t138 = -pkin(3) * t122 - pkin(9) * t118;
t91 = -pkin(2) + t138;
t172 = t117 * t91 + t101;
t110 = t121 ^ 2;
t165 = t117 ^ 2 - t110;
t141 = t165 * qJD(4);
t112 = sin(pkin(12));
t177 = pkin(4) * t112;
t176 = pkin(8) * t117;
t175 = -qJ(5) - pkin(9);
t114 = cos(pkin(12));
t157 = qJD(5) * t121;
t106 = t118 * qJD(3);
t152 = t117 * t106;
t137 = pkin(3) * t118 - pkin(9) * t122;
t88 = t137 * qJD(3);
t173 = pkin(8) * t152 + t121 * t88;
t20 = -t118 * t157 + (pkin(4) * t118 - qJ(5) * t166) * qJD(3) + (-t101 + (qJ(5) * t118 - t91) * t117) * qJD(4) + t173;
t167 = t118 * t121;
t174 = -t117 * t88 - t91 * t159;
t26 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t167 + (-qJD(5) * t118 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t122) * t117 - t174;
t12 = t112 * t20 + t114 * t26;
t84 = t121 * t91;
t49 = -qJ(5) * t167 + t84 + (-pkin(4) - t176) * t122;
t168 = t117 * t118;
t56 = -qJ(5) * t168 + t172;
t28 = t112 * t49 + t114 * t56;
t142 = qJD(4) * t175;
t71 = t117 * t142 + t157;
t72 = -qJD(5) * t117 + t121 * t142;
t42 = t112 * t72 + t114 * t71;
t92 = t175 * t117;
t93 = t175 * t121;
t58 = t112 * t92 - t114 * t93;
t115 = cos(pkin(6));
t113 = sin(pkin(6));
t119 = sin(qJ(2));
t170 = t113 * t119;
t76 = -t115 * t122 + t118 * t170;
t171 = qJD(3) * t76;
t123 = cos(qJ(2));
t169 = t113 * t123;
t89 = pkin(4) * t168 + t118 * pkin(8);
t109 = t118 ^ 2;
t164 = -t122 ^ 2 + t109;
t163 = qJD(2) * t119;
t162 = qJD(3) * t121;
t161 = qJD(4) * t117;
t160 = qJD(4) * t118;
t158 = qJD(4) * t122;
t155 = -0.2e1 * pkin(2) * qJD(3);
t154 = -0.2e1 * pkin(3) * qJD(4);
t104 = pkin(8) * t156;
t60 = t178 * pkin(4) + t104;
t105 = pkin(4) * t161;
t153 = t76 * t161;
t103 = -pkin(4) * t121 - pkin(3);
t150 = t121 * t156;
t149 = t117 * t158;
t147 = t121 * t158;
t146 = t113 * t163;
t145 = qJD(2) * t169;
t144 = t117 * t159;
t143 = t118 * t156;
t10 = -t112 * t26 + t114 * t20;
t27 = -t112 * t56 + t114 * t49;
t41 = -t112 * t71 + t114 * t72;
t57 = t112 * t93 + t114 * t92;
t140 = t164 * qJD(3);
t139 = t121 * t143;
t116 = sin(qJ(6));
t120 = cos(qJ(6));
t70 = -t112 * t168 + t114 * t167;
t16 = -pkin(5) * t122 - pkin(10) * t70 + t27;
t82 = t112 * t121 + t114 * t117;
t69 = t82 * t118;
t17 = -pkin(10) * t69 + t28;
t136 = t116 * t17 - t120 * t16;
t135 = t116 * t16 + t120 * t17;
t77 = t115 * t118 + t122 * t170;
t125 = t117 * t169 - t77 * t121;
t54 = -t77 * t117 - t121 * t169;
t30 = t112 * t125 + t114 * t54;
t31 = t112 * t54 - t114 * t125;
t134 = t116 * t31 - t120 * t30;
t133 = t116 * t30 + t120 * t31;
t43 = -pkin(10) * t82 + t57;
t130 = t112 * t117 - t114 * t121;
t44 = -pkin(10) * t130 + t58;
t132 = t116 * t44 - t120 * t43;
t131 = t116 * t43 + t120 * t44;
t37 = t116 * t70 + t120 * t69;
t38 = -t116 * t69 + t120 * t70;
t47 = t116 * t82 + t120 * t130;
t48 = -t116 * t130 + t120 * t82;
t102 = pkin(4) * t114 + pkin(5);
t129 = t102 * t116 + t120 * t177;
t128 = -t102 * t120 + t116 * t177;
t53 = qJD(3) * t77 + t118 * t145;
t127 = t117 * t53 + t76 * t159;
t126 = -t121 * t53 + t153;
t124 = t121 * t106 + t149;
t40 = t112 * t151 - t114 * t150 + t82 * t160;
t8 = pkin(5) * t106 + pkin(10) * t40 + t10;
t39 = t130 * t160 - t82 * t156;
t9 = pkin(10) * t39 + t12;
t2 = -t135 * qJD(6) - t116 * t9 + t120 * t8;
t1 = t136 * qJD(6) - t116 * t8 - t120 * t9;
t98 = -0.2e1 * t143;
t75 = -t112 * t161 + t114 * t159;
t74 = t82 * qJD(4);
t67 = t129 * qJD(6);
t66 = t128 * qJD(6);
t62 = pkin(5) * t130 + t103;
t59 = pkin(5) * t74 + t105;
t52 = t122 * t145 - t171;
t50 = pkin(5) * t69 + t89;
t35 = -t172 * qJD(4) + t173;
t34 = pkin(8) * t124 + t174;
t33 = -pkin(10) * t74 + t42;
t32 = -pkin(10) * t75 + t41;
t29 = -pkin(5) * t39 + t60;
t25 = t54 * qJD(4) + t117 * t146 + t52 * t121;
t24 = qJD(4) * t125 - t52 * t117 + t121 * t146;
t22 = t48 * qJD(6) + t116 * t75 + t120 * t74;
t21 = -t47 * qJD(6) - t116 * t74 + t120 * t75;
t15 = t38 * qJD(6) - t116 * t40 - t120 * t39;
t14 = -t37 * qJD(6) + t116 * t39 - t120 * t40;
t13 = t112 * t24 + t114 * t25;
t11 = -t112 * t25 + t114 * t24;
t6 = -t131 * qJD(6) - t116 * t33 + t120 * t32;
t5 = t132 * qJD(6) - t116 * t32 - t120 * t33;
t4 = -t133 * qJD(6) + t11 * t120 - t116 * t13;
t3 = t134 * qJD(6) - t11 * t116 - t120 * t13;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t11 * t30 + 0.2e1 * t13 * t31 + 0.2e1 * t53 * t76, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t146, -t145, 0, 0, 0, 0, 0 (-t123 * t106 - t122 * t163) * t113 (t118 * t163 - t123 * t156) * t113, 0, 0, 0, 0, 0 (t117 * t171 - t24) * t122 + (qJD(3) * t54 + t127) * t118 (t76 * t162 + t25) * t122 + (qJD(3) * t125 - t126) * t118, -t11 * t70 - t13 * t69 + t30 * t40 + t31 * t39, t10 * t30 + t11 * t27 + t12 * t31 + t13 * t28 + t53 * t89 + t60 * t76, 0, 0, 0, 0, 0, -t134 * t106 - t4 * t122 + t76 * t15 + t53 * t37, -t133 * t106 - t3 * t122 + t76 * t14 + t53 * t38; 0, 0, 0, 0, 0.2e1 * t143, -0.2e1 * t140, 0, 0, 0, t118 * t155, t122 * t155, -0.2e1 * t109 * t144 + 0.2e1 * t110 * t143, 0.2e1 * t109 * t141 - 0.4e1 * t117 * t139, 0.2e1 * t118 * t149 + 0.2e1 * t164 * t162, -0.2e1 * t117 * t140 + 0.2e1 * t118 * t147, t98, 0.2e1 * t84 * t106 - 0.2e1 * t35 * t122 + 0.2e1 * (t109 * t159 + t117 * t143) * pkin(8), -0.2e1 * t34 * t122 - 0.2e1 * t172 * t106 + 0.2e1 * (-t109 * t161 + 0.2e1 * t139) * pkin(8), -0.2e1 * t10 * t70 - 0.2e1 * t12 * t69 + 0.2e1 * t27 * t40 + 0.2e1 * t28 * t39, 0.2e1 * t10 * t27 + 0.2e1 * t12 * t28 + 0.2e1 * t60 * t89, 0.2e1 * t38 * t14, -0.2e1 * t14 * t37 - 0.2e1 * t15 * t38, 0.2e1 * t38 * t106 - 0.2e1 * t122 * t14, -0.2e1 * t37 * t106 + 0.2e1 * t122 * t15, t98, -0.2e1 * t136 * t106 - 0.2e1 * t2 * t122 + 0.2e1 * t50 * t15 + 0.2e1 * t29 * t37, -0.2e1 * t1 * t122 - 0.2e1 * t135 * t106 + 0.2e1 * t50 * t14 + 0.2e1 * t29 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t52, 0, 0, 0, 0, 0, t126, t127, -t11 * t82 - t13 * t130 - t30 * t75 - t31 * t74, pkin(4) * t153 + t103 * t53 + t11 * t57 + t13 * t58 + t30 * t41 + t31 * t42, 0, 0, 0, 0, 0, t22 * t76 + t47 * t53, t21 * t76 + t48 * t53; 0, 0, 0, 0, 0, 0, t156, -t106, 0, -t104, pkin(8) * t106, t117 * t150 - t118 * t141, -0.4e1 * t118 * t144 - t165 * t156, -t147 + t152, t124, 0 (pkin(9) * t166 + (-pkin(3) * t121 + t176) * t118) * qJD(4) + (t138 * t117 - t101) * qJD(3) (pkin(8) * t167 + t117 * t137) * qJD(4) + (t121 * t138 + t122 * t176) * qJD(3), -t10 * t82 - t12 * t130 - t27 * t75 - t28 * t74 + t39 * t58 + t40 * t57 - t41 * t70 - t42 * t69, t10 * t57 + t103 * t60 + t105 * t89 + t12 * t58 + t27 * t41 + t28 * t42, t14 * t48 + t21 * t38, -t14 * t47 - t15 * t48 - t21 * t37 - t22 * t38, t48 * t106 - t122 * t21, -t47 * t106 + t122 * t22, 0, -t132 * t106 - t6 * t122 + t62 * t15 + t50 * t22 + t29 * t47 + t59 * t37, -t131 * t106 - t5 * t122 + t62 * t14 + t50 * t21 + t29 * t48 + t59 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t144, -0.2e1 * t141, 0, 0, 0, t117 * t154, t121 * t154, -0.2e1 * t130 * t42 - 0.2e1 * t41 * t82 - 0.2e1 * t57 * t75 - 0.2e1 * t58 * t74, 0.2e1 * t103 * t105 + 0.2e1 * t41 * t57 + 0.2e1 * t42 * t58, 0.2e1 * t48 * t21, -0.2e1 * t21 * t47 - 0.2e1 * t22 * t48, 0, 0, 0, 0.2e1 * t22 * t62 + 0.2e1 * t47 * t59, 0.2e1 * t21 * t62 + 0.2e1 * t48 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t25, 0 (t11 * t114 + t112 * t13) * pkin(4), 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117 * t160 + t150, -t178, t106, t35, t34 (t112 * t39 + t114 * t40) * pkin(4) (t10 * t114 + t112 * t12) * pkin(4), 0, 0, t14, -t15, t106, -t128 * t106 + t67 * t122 + t2, -t129 * t106 - t66 * t122 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, -t161, 0, -pkin(9) * t159, pkin(9) * t161 (-t112 * t74 - t114 * t75) * pkin(4) (t112 * t42 + t114 * t41) * pkin(4), 0, 0, t21, -t22, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t67, 0.2e1 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, t15, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, 0, 0, 0, 0, 0, t22, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, t106, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t22, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;

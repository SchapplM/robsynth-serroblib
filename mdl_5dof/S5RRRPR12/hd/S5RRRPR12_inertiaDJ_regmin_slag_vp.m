% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR12_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:35
% EndTime: 2019-12-31 21:40:41
% DurationCPUTime: 1.65s
% Computational Cost: add. (2698->277), mult. (7687->541), div. (0->0), fcn. (7245->10), ass. (0->128)
t113 = cos(pkin(10));
t106 = -t113 * pkin(4) - pkin(3);
t165 = 0.2e1 * t106;
t112 = sin(pkin(5));
t164 = 0.2e1 * t112;
t117 = sin(qJ(2));
t163 = pkin(1) * t117;
t162 = pkin(8) * t112;
t161 = pkin(9) + qJ(4);
t111 = sin(pkin(10));
t120 = cos(qJ(2));
t150 = qJD(2) * t117;
t116 = sin(qJ(3));
t119 = cos(qJ(3));
t147 = qJD(3) * t119;
t148 = qJD(3) * t116;
t114 = cos(pkin(5));
t155 = t112 * t120;
t139 = pkin(7) * t155;
t73 = t139 + (pkin(8) + t163) * t114;
t74 = (-pkin(2) * t120 - pkin(8) * t117 - pkin(1)) * t112;
t77 = (pkin(2) * t117 - pkin(8) * t120) * t112 * qJD(2);
t135 = t112 * t150;
t149 = qJD(2) * t120;
t78 = -t114 * pkin(1) * t149 + pkin(7) * t135;
t25 = -t116 * t77 + t119 * t78 - t74 * t147 + t73 * t148;
t20 = (qJ(4) * t150 - qJD(4) * t120) * t112 - t25;
t134 = t112 * t149;
t156 = t112 * t117;
t84 = t114 * t116 + t119 * t156;
t54 = qJD(3) * t84 + t116 * t134;
t83 = -t114 * t119 + t116 * t156;
t55 = -qJD(3) * t83 + t119 * t134;
t79 = (t114 * t163 + t139) * qJD(2);
t24 = t54 * pkin(3) - t55 * qJ(4) - t84 * qJD(4) + t79;
t8 = t111 * t24 + t113 * t20;
t72 = pkin(7) * t156 + (-pkin(1) * t120 - pkin(2)) * t114;
t35 = t83 * pkin(3) - t84 * qJ(4) + t72;
t160 = t116 * t74 + t119 * t73;
t36 = -qJ(4) * t155 + t160;
t18 = t111 * t35 + t113 * t36;
t138 = pkin(8) * t148;
t80 = -t116 * qJD(4) + (pkin(3) * t116 - qJ(4) * t119) * qJD(3);
t50 = t111 * t138 + t113 * t80;
t26 = t116 * t78 + t119 * t77 - t73 * t147 - t74 * t148;
t23 = -pkin(3) * t135 - t26;
t159 = t23 * t111;
t158 = t23 * t113;
t153 = t113 * t119;
t103 = pkin(8) * t153;
t128 = -t119 * pkin(3) - t116 * qJ(4);
t93 = -pkin(2) + t128;
t66 = t111 * t93 + t103;
t157 = t111 * t119;
t154 = t113 * t116;
t115 = sin(qJ(5));
t152 = t115 * t111;
t118 = cos(qJ(5));
t151 = t118 * t113;
t146 = qJD(3) * t120;
t145 = qJD(4) * t119;
t144 = qJD(5) * t115;
t143 = qJD(5) * t116;
t142 = qJD(5) * t118;
t141 = -0.2e1 * pkin(2) * qJD(3);
t140 = pkin(8) * t157;
t107 = pkin(8) * t147;
t109 = t112 ^ 2;
t137 = t109 * t149;
t136 = t111 * t147;
t133 = t116 * t147;
t7 = -t111 * t20 + t113 * t24;
t17 = -t111 * t36 + t113 * t35;
t132 = -t116 * t73 + t119 * t74;
t42 = t55 * t111 - t113 * t135;
t131 = t42 * pkin(9) - t8;
t130 = t117 * t137;
t129 = 0.2e1 * (t111 ^ 2 + t113 ^ 2) * qJD(4);
t37 = pkin(3) * t155 - t132;
t53 = -t111 * t155 + t84 * t113;
t10 = t83 * pkin(4) - t53 * pkin(9) + t17;
t52 = t84 * t111 + t113 * t155;
t15 = -t52 * pkin(9) + t18;
t4 = t115 * t10 + t118 * t15;
t70 = t111 * t80;
t51 = -t113 * t138 + t70;
t127 = -t50 * t111 + t51 * t113;
t88 = t113 * t93;
t48 = -pkin(9) * t154 + t88 + (-pkin(8) * t111 - pkin(4)) * t119;
t56 = -t111 * t116 * pkin(9) + t66;
t29 = t115 * t48 + t118 * t56;
t30 = t115 * t53 + t118 * t52;
t31 = -t115 * t52 + t118 * t53;
t95 = t161 * t111;
t96 = t161 * t113;
t58 = -t115 * t95 + t118 * t96;
t126 = -qJ(4) * t54 - qJD(4) * t83;
t90 = t118 * t111 + t115 * t113;
t89 = -t151 + t152;
t125 = pkin(4) * t116 - pkin(9) * t153;
t124 = -pkin(8) * t154 - pkin(9) * t157;
t123 = t116 * t146 + t119 * t150;
t122 = t116 * t150 - t119 * t146;
t43 = t111 * t135 + t55 * t113;
t121 = t54 * pkin(4) - t43 * pkin(9) + t7;
t91 = (pkin(4) * t111 + pkin(8)) * t116;
t85 = pkin(4) * t136 + t107;
t82 = t90 * qJD(5);
t81 = t89 * qJD(5);
t76 = t89 * t116;
t75 = t90 * t116;
t65 = t88 - t140;
t57 = -t115 * t96 - t118 * t95;
t45 = t142 * t154 - t143 * t152 + t90 * t147;
t44 = -t90 * t143 - t89 * t147;
t41 = -qJD(4) * t90 - qJD(5) * t58;
t40 = -qJD(4) * t151 + t95 * t142 + (qJD(4) * t111 + qJD(5) * t96) * t115;
t28 = -t115 * t56 + t118 * t48;
t27 = t52 * pkin(4) + t37;
t16 = t42 * pkin(4) + t23;
t14 = -t115 * t70 + t118 * t50 - t29 * qJD(5) + (-t115 * t124 + t118 * t125) * qJD(3);
t13 = t56 * t144 - t115 * (qJD(3) * t125 + t50) - t118 * (qJD(3) * t124 + t70) - t48 * t142;
t12 = qJD(5) * t31 + t115 * t43 + t118 * t42;
t11 = -qJD(5) * t30 - t115 * t42 + t118 * t43;
t3 = t118 * t10 - t115 * t15;
t2 = -qJD(5) * t4 + t115 * t131 + t118 * t121;
t1 = -t10 * t142 - t115 * t121 + t118 * t131 + t15 * t144;
t5 = [0, 0, 0, 0.2e1 * t130, 0.2e1 * (-t117 ^ 2 + t120 ^ 2) * t109 * qJD(2), 0.2e1 * t114 * t134, -0.2e1 * t114 * t135, 0, -0.2e1 * t109 * pkin(1) * t150 - 0.2e1 * t79 * t114, -0.2e1 * pkin(1) * t137 + 0.2e1 * t78 * t114, 0.2e1 * t84 * t55, -0.2e1 * t84 * t54 - 0.2e1 * t55 * t83, (-t120 * t55 + t84 * t150) * t164, (t120 * t54 - t83 * t150) * t164, -0.2e1 * t130, 0.2e1 * t72 * t54 + 0.2e1 * t79 * t83 + 0.2e1 * (-t26 * t120 + t132 * t150) * t112, 0.2e1 * t72 * t55 + 0.2e1 * t79 * t84 + 0.2e1 * (-t25 * t120 - t160 * t150) * t112, 0.2e1 * t17 * t54 + 0.2e1 * t23 * t52 + 0.2e1 * t37 * t42 + 0.2e1 * t7 * t83, -0.2e1 * t18 * t54 + 0.2e1 * t23 * t53 + 0.2e1 * t37 * t43 - 0.2e1 * t8 * t83, -0.2e1 * t17 * t43 - 0.2e1 * t18 * t42 - 0.2e1 * t8 * t52 - 0.2e1 * t7 * t53, 0.2e1 * t17 * t7 + 0.2e1 * t18 * t8 + 0.2e1 * t37 * t23, 0.2e1 * t31 * t11, -0.2e1 * t11 * t30 - 0.2e1 * t31 * t12, 0.2e1 * t11 * t83 + 0.2e1 * t31 * t54, -0.2e1 * t12 * t83 - 0.2e1 * t30 * t54, 0.2e1 * t83 * t54, 0.2e1 * t27 * t12 + 0.2e1 * t16 * t30 + 0.2e1 * t2 * t83 + 0.2e1 * t3 * t54, 0.2e1 * t1 * t83 + 0.2e1 * t27 * t11 + 0.2e1 * t16 * t31 - 0.2e1 * t4 * t54; 0, 0, 0, 0, 0, t134, -t135, 0, -t79, t78, t55 * t116 + t84 * t147, -t116 * t54 + t55 * t119 + (-t116 * t84 - t119 * t83) * qJD(3), t122 * t112, t123 * t112, 0, -pkin(2) * t54 - t79 * t119 - t122 * t162 + t72 * t148, -pkin(2) * t55 + t79 * t116 - t123 * t162 + t72 * t147, -t7 * t119 + t50 * t83 + t65 * t54 + (pkin(8) * t42 + t159) * t116 + (t116 * t17 + (pkin(8) * t52 + t111 * t37) * t119) * qJD(3), t8 * t119 - t51 * t83 - t66 * t54 + (pkin(8) * t43 + t158) * t116 + (-t116 * t18 + (pkin(8) * t53 + t113 * t37) * t119) * qJD(3), -t66 * t42 - t65 * t43 - t50 * t53 - t51 * t52 + (-t111 * t8 - t113 * t7) * t116 + (-t111 * t18 - t113 * t17) * t147, t17 * t50 + t18 * t51 + t7 * t65 + t8 * t66 + (t116 * t23 + t37 * t147) * pkin(8), -t11 * t76 + t31 * t44, -t11 * t75 + t76 * t12 - t44 * t30 - t31 * t45, -t11 * t119 + t148 * t31 + t44 * t83 - t76 * t54, t12 * t119 - t148 * t30 - t45 * t83 - t75 * t54, -t54 * t119 + t148 * t83, -t2 * t119 + t91 * t12 + t14 * t83 + t148 * t3 + t16 * t75 + t27 * t45 + t28 * t54 + t85 * t30, -t1 * t119 + t91 * t11 + t13 * t83 - t148 * t4 - t16 * t76 + t27 * t44 - t29 * t54 + t85 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t133, 0.2e1 * (-t116 ^ 2 + t119 ^ 2) * qJD(3), 0, 0, 0, t116 * t141, t119 * t141, -0.2e1 * t50 * t119 + 0.2e1 * (t65 + 0.2e1 * t140) * t148, 0.2e1 * t51 * t119 + 0.2e1 * (-t66 + 0.2e1 * t103) * t148, 0.2e1 * (-t111 * t51 - t113 * t50) * t116 + 0.2e1 * (-t111 * t66 - t113 * t65) * t147, 0.2e1 * pkin(8) ^ 2 * t133 + 0.2e1 * t65 * t50 + 0.2e1 * t66 * t51, -0.2e1 * t76 * t44, -0.2e1 * t44 * t75 + 0.2e1 * t76 * t45, -0.2e1 * t44 * t119 - 0.2e1 * t148 * t76, 0.2e1 * t45 * t119 - 0.2e1 * t148 * t75, -0.2e1 * t133, -0.2e1 * t14 * t119 + 0.2e1 * t148 * t28 + 0.2e1 * t91 * t45 + 0.2e1 * t85 * t75, -0.2e1 * t13 * t119 - 0.2e1 * t148 * t29 + 0.2e1 * t91 * t44 - 0.2e1 * t85 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t54, t135, t26, t25, -pkin(3) * t42 + t111 * t126 - t158, -pkin(3) * t43 + t113 * t126 + t159, (-qJ(4) * t42 - qJD(4) * t52 + t8) * t113 + (qJ(4) * t43 + qJD(4) * t53 - t7) * t111, -t23 * pkin(3) + (-t111 * t17 + t113 * t18) * qJD(4) + (-t7 * t111 + t8 * t113) * qJ(4), t11 * t90 - t31 * t81, -t11 * t89 - t90 * t12 + t81 * t30 - t31 * t82, t90 * t54 - t81 * t83, -t89 * t54 - t82 * t83, 0, t106 * t12 + t16 * t89 + t27 * t82 + t41 * t83 + t57 * t54, t106 * t11 + t16 * t90 - t27 * t81 + t40 * t83 - t58 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, -t148, 0, -t107, t138, t111 * t145 + (t111 * t128 - t103) * qJD(3), t113 * t145 + (t113 * t128 + t140) * qJD(3), t127, -pkin(3) * t107 + (-t111 * t65 + t113 * t66) * qJD(4) + t127 * qJ(4), t44 * t90 + t76 * t81, -t44 * t89 - t90 * t45 + t81 * t75 + t76 * t82, t81 * t119 + t148 * t90, t82 * t119 - t148 * t89, 0, t106 * t45 - t41 * t119 + t148 * t57 + t91 * t82 + t85 * t89, t106 * t44 - t40 * t119 - t148 * t58 - t91 * t81 + t85 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, qJ(4) * t129, -0.2e1 * t90 * t81, 0.2e1 * t81 * t89 - 0.2e1 * t90 * t82, 0, 0, 0, t82 * t165, -t81 * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t43, 0, t23, 0, 0, 0, 0, 0, t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t113 * t147, 0, t107, 0, 0, 0, 0, 0, t45, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t12, t54, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t45, t148, t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t82, 0, t41, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;

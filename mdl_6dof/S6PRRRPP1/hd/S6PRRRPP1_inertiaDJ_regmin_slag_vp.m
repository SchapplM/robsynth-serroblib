% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRPP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:46:54
% EndTime: 2019-03-08 22:46:58
% DurationCPUTime: 1.63s
% Computational Cost: add. (2218->236), mult. (5901->422), div. (0->0), fcn. (5297->10), ass. (0->129)
t101 = cos(qJ(3));
t151 = qJD(3) * t101;
t97 = sin(qJ(4));
t138 = t97 * t151;
t100 = cos(qJ(4));
t149 = qJD(4) * t100;
t98 = sin(qJ(3));
t176 = t98 * t149 + t138;
t119 = -pkin(3) * t101 - pkin(9) * t98;
t74 = -pkin(2) + t119;
t154 = t100 * t101;
t82 = pkin(8) * t154;
t165 = t97 * t74 + t82;
t102 = cos(qJ(2));
t96 = sin(pkin(6));
t160 = t102 * t96;
t134 = qJD(2) * t160;
t159 = cos(pkin(6));
t99 = sin(qJ(2));
t170 = t96 * t99;
t59 = -t159 * t101 + t98 * t170;
t104 = -t59 * qJD(3) + t101 * t134;
t175 = -qJD(4) * t160 + t104;
t110 = t101 * t170 + t159 * t98;
t157 = qJD(2) * t99;
t139 = t96 * t157;
t174 = qJD(4) * t110 - t139;
t158 = cos(pkin(11));
t123 = t158 * t100;
t155 = qJD(4) * t97;
t95 = sin(pkin(11));
t58 = qJD(4) * t123 - t95 * t155;
t93 = t100 ^ 2;
t164 = t97 ^ 2 - t93;
t125 = t164 * qJD(4);
t173 = 2 * qJD(6);
t172 = pkin(8) * t97;
t171 = t95 * t97;
t169 = t97 * t98;
t168 = -qJ(5) - pkin(9);
t147 = qJD(5) * t100;
t162 = qJ(5) * t98;
t156 = qJD(3) * t98;
t142 = t97 * t156;
t120 = pkin(3) * t98 - pkin(9) * t101;
t68 = t120 * qJD(3);
t166 = pkin(8) * t142 + t100 * t68;
t18 = -t98 * t147 + (pkin(4) * t98 - qJ(5) * t154) * qJD(3) + (-t82 + (-t74 + t162) * t97) * qJD(4) + t166;
t161 = t100 * t98;
t167 = -t74 * t149 - t97 * t68;
t23 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t161 + (-qJD(5) * t98 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t101) * t97 - t167;
t9 = t158 * t23 + t95 * t18;
t66 = t100 * t74;
t45 = -qJ(5) * t161 + t66 + (-pkin(4) - t172) * t101;
t49 = -t97 * t162 + t165;
t27 = t158 * t49 + t95 * t45;
t69 = pkin(4) * t169 + t98 * pkin(8);
t92 = t98 ^ 2;
t163 = -t101 ^ 2 + t92;
t152 = qJD(3) * t100;
t150 = qJD(3) * t102;
t148 = qJD(4) * t101;
t146 = qJD(6) * t101;
t145 = -0.2e1 * pkin(2) * qJD(3);
t144 = -0.2e1 * pkin(3) * qJD(4);
t143 = qJ(6) * t156 + t9;
t88 = pkin(8) * t151;
t52 = t176 * pkin(4) + t88;
t89 = pkin(4) * t155;
t141 = t59 * t155;
t136 = t97 * t148;
t135 = t59 * t149;
t132 = t97 * t149;
t131 = t98 * t151;
t87 = -pkin(4) * t100 - pkin(3);
t130 = t100 * t151;
t129 = t100 * t148;
t7 = t158 * t18 - t23 * t95;
t126 = qJD(4) * t168;
t108 = -qJD(5) * t97 + t100 * t126;
t56 = t97 * t126 + t147;
t40 = -t158 * t108 + t56 * t95;
t41 = t95 * t108 + t158 * t56;
t127 = t158 * t97;
t75 = t168 * t100;
t50 = -t168 * t127 - t75 * t95;
t51 = -t158 * t75 + t168 * t171;
t128 = t40 * t50 + t51 * t41;
t124 = t163 * qJD(3);
t122 = t98 * t130;
t118 = t158 * t151;
t103 = t174 * t100 + t175 * t97;
t20 = t175 * t100 - t174 * t97;
t10 = -t95 * t103 + t158 * t20;
t105 = -t100 * t160 - t110 * t97;
t48 = t110 * t100 - t97 * t160;
t28 = -t158 * t105 + t48 * t95;
t29 = t95 * t105 + t158 * t48;
t47 = t110 * qJD(3) + t98 * t134;
t8 = t158 * t103 + t20 * t95;
t116 = 0.2e1 * t29 * t10 + 0.2e1 * t28 * t8 + 0.2e1 * t59 * t47;
t115 = -t100 * t47 + t141;
t38 = -t97 * t118 - t95 * t130 - t58 * t98;
t64 = t95 * t100 + t127;
t57 = t64 * qJD(4);
t39 = -t100 * t118 + t95 * t138 + t98 * t57;
t54 = t64 * t98;
t55 = t98 * t123 - t95 * t169;
t114 = -t10 * t54 - t28 * t39 + t29 * t38 + t55 * t8;
t113 = t10 * t51 + t28 * t40 + t29 * t41 + t50 * t8;
t63 = -t123 + t171;
t112 = -t10 * t63 + t28 * t58 - t29 * t57 + t64 * t8;
t26 = t158 * t45 - t95 * t49;
t111 = t51 * t38 - t39 * t50 + t40 * t55 - t41 * t54;
t109 = t98 * t152 + t136;
t107 = 0.2e1 * t40 * t64 - 0.2e1 * t41 * t63 + 0.2e1 * t50 * t58 - 0.2e1 * t51 * t57;
t85 = -t158 * pkin(4) - pkin(5);
t83 = pkin(4) * t95 + qJ(6);
t42 = pkin(5) * t63 - qJ(6) * t64 + t87;
t35 = -t165 * qJD(4) + t166;
t34 = t109 * pkin(8) + t167;
t33 = pkin(5) * t54 - qJ(6) * t55 + t69;
t31 = pkin(5) * t57 - qJ(6) * t58 - qJD(6) * t64 + t89;
t22 = t101 * pkin(5) - t26;
t21 = -qJ(6) * t101 + t27;
t11 = -pkin(5) * t38 + qJ(6) * t39 - qJD(6) * t55 + t52;
t6 = -pkin(5) * t156 - t7;
t5 = t143 - t146;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, 0, 0, 0, t116; 0, 0, -t139, -t134, 0, 0, 0, 0, 0 (-t101 * t157 - t98 * t150) * t96 (-t101 * t150 + t98 * t157) * t96, 0, 0, 0, 0, 0, t103 * t101 + t105 * t156 + t98 * t135 + t59 * t138 + t47 * t169 (t59 * t152 + t20) * t101 + (-qJD(3) * t48 - t115) * t98, t114, t10 * t27 - t26 * t8 - t28 * t7 + t29 * t9 + t47 * t69 + t52 * t59, t101 * t8 - t28 * t156 - t38 * t59 + t47 * t54, t114, -t10 * t101 + t29 * t156 + t39 * t59 - t47 * t55, t10 * t21 + t11 * t59 + t22 * t8 + t28 * t6 + t29 * t5 + t33 * t47; 0, 0, 0, 0, 0.2e1 * t131, -0.2e1 * t124, 0, 0, 0, t98 * t145, t101 * t145, 0.2e1 * t93 * t131 - 0.2e1 * t92 * t132, -0.4e1 * t97 * t122 + 0.2e1 * t92 * t125, 0.2e1 * t98 * t136 + 0.2e1 * t163 * t152, -0.2e1 * t97 * t124 + 0.2e1 * t98 * t129, -0.2e1 * t131, 0.2e1 * t66 * t156 - 0.2e1 * t35 * t101 + 0.2e1 * (t97 * t131 + t92 * t149) * pkin(8), -0.2e1 * t34 * t101 - 0.2e1 * t165 * t156 + 0.2e1 * (-t92 * t155 + 0.2e1 * t122) * pkin(8), 0.2e1 * t26 * t39 + 0.2e1 * t27 * t38 - 0.2e1 * t54 * t9 - 0.2e1 * t55 * t7, 0.2e1 * t26 * t7 + 0.2e1 * t27 * t9 + 0.2e1 * t52 * t69, 0.2e1 * t101 * t6 + 0.2e1 * t11 * t54 - 0.2e1 * t22 * t156 - 0.2e1 * t33 * t38, 0.2e1 * t21 * t38 - 0.2e1 * t22 * t39 - 0.2e1 * t5 * t54 + 0.2e1 * t55 * t6, -0.2e1 * t101 * t5 - 0.2e1 * t11 * t55 + 0.2e1 * t21 * t156 + 0.2e1 * t33 * t39, 0.2e1 * t11 * t33 + 0.2e1 * t21 * t5 + 0.2e1 * t22 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t104, 0, 0, 0, 0, 0, t115, t47 * t97 + t135, t112, pkin(4) * t141 + t47 * t87 + t113, t47 * t63 + t57 * t59, t112, -t47 * t64 - t58 * t59, t31 * t59 + t42 * t47 + t113; 0, 0, 0, 0, 0, 0, t151, -t156, 0, -t88, pkin(8) * t156, -t98 * t125 + t97 * t130, -0.4e1 * t98 * t132 - t164 * t151, -t129 + t142, t109, 0 (pkin(8) * t169 - t120 * t100) * qJD(4) + (-pkin(9) * t169 + (-pkin(3) * t97 - pkin(8) * t100) * t101) * qJD(3) (pkin(8) * t161 + t120 * t97) * qJD(4) + (t119 * t100 + t101 * t172) * qJD(3), -t26 * t58 - t27 * t57 - t63 * t9 - t64 * t7 + t111, -t26 * t40 + t27 * t41 - t50 * t7 + t51 * t9 + t52 * t87 + t69 * t89, t101 * t40 + t11 * t63 - t50 * t156 + t31 * t54 + t33 * t57 - t38 * t42, -t21 * t57 + t22 * t58 - t5 * t63 + t6 * t64 + t111, -t101 * t41 - t11 * t64 + t51 * t156 - t31 * t55 - t33 * t58 + t39 * t42, t11 * t42 + t21 * t41 + t22 * t40 + t31 * t33 + t5 * t51 + t50 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t132, -0.2e1 * t125, 0, 0, 0, t97 * t144, t100 * t144, t107, 0.2e1 * t87 * t89 + 0.2e1 * t128, 0.2e1 * t31 * t63 + 0.2e1 * t42 * t57, t107, -0.2e1 * t31 * t64 - 0.2e1 * t42 * t58, 0.2e1 * t31 * t42 + 0.2e1 * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t20, 0 (t10 * t95 - t158 * t8) * pkin(4), -t8, 0, t10, qJD(6) * t29 + t10 * t83 + t8 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98 * t155 + t130, -t176, t156, t35, t34 (t158 * t39 + t38 * t95) * pkin(4) (t158 * t7 + t9 * t95) * pkin(4) (pkin(5) - t85) * t156 + t7, -qJD(6) * t54 + t38 * t83 - t39 * t85, t83 * t156 + t143 - 0.2e1 * t146, qJD(6) * t21 + t5 * t83 + t6 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, -t155, 0, -pkin(9) * t149, pkin(9) * t155 (-t158 * t58 - t57 * t95) * pkin(4) (-t158 * t40 + t41 * t95) * pkin(4), -t40, -qJD(6) * t63 - t57 * t83 + t58 * t85, t41, qJD(6) * t51 + t40 * t85 + t41 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t83 * t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t38, 0, t39, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t57, 0, -t58, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t39, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t1;

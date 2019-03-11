% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:40:46
% EndTime: 2019-03-09 04:40:51
% DurationCPUTime: 1.45s
% Computational Cost: add. (3490->205), mult. (7920->362), div. (0->0), fcn. (7698->8), ass. (0->111)
t106 = sin(qJ(4));
t108 = cos(qJ(4));
t139 = qJD(4) * t108;
t104 = sin(pkin(9));
t109 = cos(qJ(3));
t105 = cos(pkin(9));
t107 = sin(qJ(3));
t146 = t107 * t105;
t83 = t109 * t104 + t146;
t134 = t83 * t139;
t81 = t107 * t104 - t109 * t105;
t73 = t81 * qJD(3);
t162 = -t106 * t73 + t134;
t103 = sin(pkin(10));
t150 = cos(pkin(10));
t127 = t150 * t108;
t140 = qJD(4) * t106;
t72 = qJD(4) * t127 - t103 * t140;
t149 = t103 * t106;
t161 = t127 - t149;
t128 = t150 * t106;
t82 = t103 * t108 + t128;
t102 = t108 ^ 2;
t143 = t106 ^ 2 - t102;
t125 = t143 * qJD(4);
t96 = -t105 * pkin(2) - pkin(1);
t160 = 0.2e1 * t96;
t159 = 2 * qJD(6);
t118 = qJ(5) * t73 - qJD(5) * t83;
t141 = qJD(3) * t109;
t142 = qJD(2) * t109;
t157 = pkin(7) + qJ(2);
t85 = t157 * t104;
t86 = t157 * t105;
t34 = t85 * t141 - t105 * t142 + (qJD(2) * t104 + qJD(3) * t86) * t107;
t74 = t83 * qJD(3);
t49 = t74 * pkin(3) + t73 * pkin(8);
t51 = t81 * pkin(3) - t83 * pkin(8) + t96;
t136 = t106 * t49 - t108 * t34 + t51 * t139;
t151 = t107 * t85;
t57 = t109 * t86 - t151;
t10 = -qJ(5) * t134 + (-qJD(4) * t57 + t118) * t106 + t136;
t130 = t106 * t34 + t108 * t49;
t154 = qJ(5) * t83;
t52 = t108 * t57;
t8 = t74 * pkin(4) + t118 * t108 + (-t52 + (-t51 + t154) * t106) * qJD(4) + t130;
t4 = t150 * t10 + t103 * t8;
t158 = t83 * t73;
t156 = -qJ(5) - pkin(8);
t47 = t108 * t51;
t20 = t81 * pkin(4) - t106 * t57 - t108 * t154 + t47;
t152 = t106 * t83;
t155 = t106 * t51 + t52;
t23 = -qJ(5) * t152 + t155;
t14 = t103 * t20 + t150 * t23;
t147 = t106 * t108;
t145 = t81 * qJD(6);
t144 = t82 * qJD(6);
t138 = t74 * qJ(6) + t4;
t137 = -0.2e1 * pkin(3) * qJD(4);
t98 = pkin(4) * t140;
t135 = t83 * t140;
t97 = -t108 * pkin(4) - pkin(3);
t132 = t106 * t139;
t3 = -t103 * t10 + t150 * t8;
t129 = qJD(4) * t156;
t111 = -t106 * qJD(5) + t108 * t129;
t70 = t108 * qJD(5) + t106 * t129;
t39 = t103 * t70 - t150 * t111;
t40 = t103 * t111 + t150 * t70;
t87 = t156 * t108;
t58 = -t103 * t87 - t156 * t128;
t59 = t156 * t149 - t150 * t87;
t131 = t58 * t39 + t59 * t40;
t56 = t107 * t86 + t109 * t85;
t126 = -0.4e1 * t83 * t147;
t35 = qJD(2) * t146 - qJD(3) * t151 + t104 * t142 + t86 * t141;
t124 = 0.2e1 * (t104 ^ 2 + t105 ^ 2) * qJD(2);
t36 = pkin(4) * t152 + t56;
t123 = pkin(3) * t73 - pkin(8) * t74;
t122 = pkin(3) * t83 + pkin(8) * t81;
t120 = t35 * t83 - t56 * t73;
t119 = -t73 * t81 + t83 * t74;
t26 = t162 * pkin(4) + t35;
t71 = t82 * qJD(4);
t117 = -0.2e1 * t161 * t71 + 0.2e1 * t82 * t72;
t116 = t106 * t74 + t81 * t139;
t13 = -t103 * t23 + t150 * t20;
t24 = -t72 * t83 + t82 * t73;
t25 = t161 * t73 + t71 * t83;
t42 = t82 * t83;
t43 = t161 * t83;
t115 = t59 * t24 - t58 * t25 + t39 * t43 - t40 * t42;
t114 = t161 * t25 + t82 * t24 - t72 * t42 + t71 * t43;
t113 = -t161 * t39 + t82 * t40 + t71 * t58 + t72 * t59;
t110 = 0.2e1 * t161 * t40 + 0.2e1 * t39 * t82 + 0.2e1 * t58 * t72 - 0.2e1 * t59 * t71;
t95 = -t150 * pkin(4) - pkin(5);
t92 = t103 * pkin(4) + qJ(6);
t78 = t83 ^ 2;
t50 = -pkin(5) * t161 - t82 * qJ(6) + t97;
t48 = t108 * t74 - t81 * t140;
t29 = t71 * pkin(5) - t72 * qJ(6) - t144 + t98;
t17 = t42 * pkin(5) - t43 * qJ(6) + t36;
t16 = -t155 * qJD(4) + t130;
t15 = t57 * t140 - t136;
t12 = -t81 * pkin(5) - t13;
t11 = t81 * qJ(6) + t14;
t5 = -t24 * pkin(5) + t25 * qJ(6) - t43 * qJD(6) + t26;
t2 = -t74 * pkin(5) - t3;
t1 = t138 + t145;
t6 = [0, 0, 0, 0, 0, t124, qJ(2) * t124, -0.2e1 * t158, -0.2e1 * t119, 0, 0, 0, t74 * t160, -t73 * t160, -0.2e1 * t102 * t158 - 0.2e1 * t78 * t132, 0.2e1 * t78 * t125 - t73 * t126, 0.2e1 * t108 * t119 - 0.2e1 * t81 * t135, -0.2e1 * t106 * t119 - 0.2e1 * t81 * t134, 0.2e1 * t81 * t74, 0.2e1 * t56 * t134 + 0.2e1 * t16 * t81 + 0.2e1 * t47 * t74 + 0.2e1 * (-t57 * t74 + t120) * t106, 0.2e1 * t120 * t108 - 0.2e1 * t56 * t135 + 0.2e1 * t15 * t81 - 0.2e1 * t155 * t74, 0.2e1 * t13 * t25 + 0.2e1 * t14 * t24 - 0.2e1 * t3 * t43 - 0.2e1 * t4 * t42, 0.2e1 * t13 * t3 + 0.2e1 * t14 * t4 + 0.2e1 * t36 * t26, -0.2e1 * t12 * t74 - 0.2e1 * t17 * t24 - 0.2e1 * t2 * t81 + 0.2e1 * t5 * t42, -0.2e1 * t1 * t42 + 0.2e1 * t11 * t24 - 0.2e1 * t12 * t25 + 0.2e1 * t2 * t43, 0.2e1 * t1 * t81 + 0.2e1 * t11 * t74 + 0.2e1 * t17 * t25 - 0.2e1 * t5 * t43, 0.2e1 * t11 * t1 + 0.2e1 * t12 * t2 + 0.2e1 * t17 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t73, 0, 0, 0, 0, 0, t48, -t116, t114, -t13 * t71 + t14 * t72 + t161 * t3 + t4 * t82, t161 * t74 - t71 * t81, t114, t72 * t81 + t82 * t74, t1 * t82 + t11 * t72 + t12 * t71 - t161 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, 0, 0, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t74, 0, -t35, t34, -t125 * t83 - t73 * t147, qJD(4) * t126 + t143 * t73, t116, t48, 0, -t35 * t108 + t123 * t106 + (t106 * t56 - t108 * t122) * qJD(4), t35 * t106 + t123 * t108 + (t106 * t122 + t108 * t56) * qJD(4), -t13 * t72 - t14 * t71 + t161 * t4 - t3 * t82 + t115, -t13 * t39 + t14 * t40 + t26 * t97 - t3 * t58 + t36 * t98 + t4 * t59, -t161 * t5 + t17 * t71 - t50 * t24 + t29 * t42 - t39 * t81 - t58 * t74, t1 * t161 - t11 * t71 + t12 * t72 + t2 * t82 + t115, -t17 * t72 + t50 * t25 - t29 * t43 + t40 * t81 - t5 * t82 + t59 * t74, t1 * t59 + t11 * t40 + t12 * t39 + t17 * t29 + t2 * t58 + t5 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, 0, 0, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t132, -0.2e1 * t125, 0, 0, 0, t106 * t137, t108 * t137, t110, 0.2e1 * t97 * t98 + 0.2e1 * t131, -0.2e1 * t161 * t29 + 0.2e1 * t50 * t71, t110, -0.2e1 * t29 * t82 - 0.2e1 * t50 * t72, 0.2e1 * t50 * t29 + 0.2e1 * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108 * t73 - t135, -t162, t74, t16, t15 (t103 * t24 + t150 * t25) * pkin(4) (t103 * t4 + t150 * t3) * pkin(4) (pkin(5) - t95) * t74 + t3, -qJD(6) * t42 + t92 * t24 - t95 * t25, t92 * t74 + t138 + 0.2e1 * t145, t11 * qJD(6) + t1 * t92 + t2 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, -t139, 0 (t103 * t72 - t150 * t71) * pkin(4), -t71, 0, t72, t71 * t95 + t72 * t92 + t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t140, 0, -pkin(8) * t139, pkin(8) * t140 (-t103 * t71 - t150 * t72) * pkin(4) (t103 * t40 - t150 * t39) * pkin(4), -t39, qJD(6) * t161 - t92 * t71 + t95 * t72, t40, t59 * qJD(6) + t39 * t95 + t40 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t92 * t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t24, 0, t25, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t71, 0, -t72, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t25, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t6;

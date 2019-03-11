% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:25:01
% EndTime: 2019-03-09 05:25:05
% DurationCPUTime: 1.80s
% Computational Cost: add. (2219->224), mult. (5228->413), div. (0->0), fcn. (4710->8), ass. (0->132)
t112 = cos(qJ(4));
t113 = cos(qJ(3));
t147 = qJD(4) * t113;
t137 = t112 * t147;
t109 = sin(qJ(4));
t110 = sin(qJ(3));
t145 = t110 * qJD(3);
t141 = t109 * t145;
t161 = t137 - t141;
t103 = t110 ^ 2;
t105 = t113 ^ 2;
t130 = (t103 - t105) * qJD(3);
t104 = t112 ^ 2;
t152 = t109 ^ 2 - t104;
t131 = t152 * qJD(4);
t160 = 2 * qJD(2);
t106 = sin(pkin(10));
t159 = pkin(4) * t106;
t158 = -qJ(5) - pkin(8);
t107 = cos(pkin(10));
t126 = pkin(3) * t110 - pkin(8) * t113;
t86 = qJ(2) + t126;
t77 = t109 * t86;
t114 = -pkin(1) - pkin(7);
t92 = t112 * t110 * t114;
t157 = t92 + t77;
t127 = pkin(3) * t113 + pkin(8) * t110;
t79 = t127 * qJD(3) + qJD(2);
t115 = -t157 * qJD(4) + t112 * t79;
t155 = t109 * t114;
t133 = pkin(4) - t155;
t140 = t112 * t145;
t144 = t112 * qJD(5);
t149 = qJD(4) * t109;
t18 = qJ(5) * t140 + (qJ(5) * t149 + t133 * qJD(3) - t144) * t113 + t115;
t101 = t113 * qJD(3);
t134 = t114 * t101;
t148 = qJD(4) * t112;
t142 = t109 * t79 + t112 * t134 + t86 * t148;
t146 = qJD(4) * t114;
t20 = -qJ(5) * t137 + (-qJD(5) * t113 + (qJ(5) * qJD(3) - t146) * t110) * t109 + t142;
t9 = t106 * t18 + t107 * t20;
t154 = t112 * t113;
t78 = t112 * t86;
t43 = -qJ(5) * t154 + t133 * t110 + t78;
t156 = t109 * t113;
t48 = -qJ(5) * t156 + t157;
t24 = t106 * t43 + t107 * t48;
t132 = qJD(4) * t158;
t65 = t109 * t132 + t144;
t66 = -t109 * qJD(5) + t112 * t132;
t38 = t106 * t66 + t107 * t65;
t89 = t158 * t109;
t90 = t158 * t112;
t50 = t106 * t89 - t107 * t90;
t153 = t113 * t114;
t150 = t103 + t105;
t143 = -0.2e1 * pkin(3) * qJD(4);
t100 = pkin(4) * t149;
t99 = -pkin(4) * t112 - pkin(3);
t139 = t109 * t147;
t138 = t109 * t146;
t136 = t109 * t148;
t135 = t110 * t101;
t8 = -t106 * t20 + t107 * t18;
t23 = -t106 * t48 + t107 * t43;
t37 = -t106 * t65 + t107 * t66;
t49 = t106 * t90 + t107 * t89;
t82 = pkin(4) * t156 - t153;
t129 = t109 * t134;
t128 = t109 * t140;
t108 = sin(qJ(6));
t111 = cos(qJ(6));
t64 = -t106 * t156 + t107 * t154;
t14 = pkin(5) * t110 - pkin(9) * t64 + t23;
t81 = t106 * t112 + t107 * t109;
t62 = t81 * t113;
t15 = -pkin(9) * t62 + t24;
t125 = t108 * t15 - t111 * t14;
t124 = t108 * t14 + t111 * t15;
t39 = -pkin(9) * t81 + t49;
t119 = t106 * t109 - t107 * t112;
t40 = -pkin(9) * t119 + t50;
t123 = t108 * t40 - t111 * t39;
t122 = t108 * t39 + t111 * t40;
t61 = t81 * t110;
t63 = t119 * t110;
t121 = -t108 * t63 + t111 * t61;
t120 = -t108 * t61 - t111 * t63;
t31 = t108 * t64 + t111 * t62;
t32 = -t108 * t62 + t111 * t64;
t44 = t108 * t81 + t111 * t119;
t45 = -t108 * t119 + t111 * t81;
t98 = pkin(4) * t107 + pkin(5);
t118 = t108 * t98 + t111 * t159;
t117 = t108 * t159 - t111 * t98;
t96 = t114 * t145;
t51 = t161 * pkin(4) + t96;
t116 = qJD(4) * t119;
t70 = t81 * qJD(4);
t36 = -t106 * t141 + t107 * t140 + t81 * t147;
t6 = pkin(5) * t101 + pkin(9) * t36 + t8;
t34 = t113 * t116 + t81 * t145;
t7 = pkin(9) * t34 + t9;
t2 = -t124 * qJD(6) - t108 * t7 + t111 * t6;
t1 = t125 * qJD(6) - t108 * t6 - t111 * t7;
t94 = 0.2e1 * t135;
t74 = t109 * t101 + t110 * t148;
t73 = -t139 - t140;
t72 = -t112 * t101 + t110 * t149;
t71 = -t106 * t149 + t107 * t148;
t59 = t118 * qJD(6);
t58 = t117 * qJD(6);
t54 = pkin(5) * t119 + t99;
t52 = pkin(5) * t70 + t100;
t47 = t62 * pkin(5) + t82;
t35 = -t119 * t101 - t110 * t70;
t33 = -qJD(3) * t62 + t110 * t116;
t29 = t115 - t129;
t28 = t110 * t138 - t142;
t27 = -pkin(9) * t70 + t38;
t26 = -pkin(9) * t71 + t37;
t25 = -pkin(5) * t34 + t51;
t22 = qJD(6) * t45 + t108 * t71 + t111 * t70;
t21 = -qJD(6) * t44 - t108 * t70 + t111 * t71;
t13 = qJD(6) * t32 - t108 * t36 - t111 * t34;
t12 = -qJD(6) * t120 - t108 * t35 + t111 * t33;
t11 = -qJD(6) * t31 + t108 * t34 - t111 * t36;
t10 = t121 * qJD(6) - t108 * t33 - t111 * t35;
t4 = -t122 * qJD(6) - t108 * t27 + t111 * t26;
t3 = t123 * qJD(6) - t108 * t26 - t111 * t27;
t5 = [0, 0, 0, 0, t160, qJ(2) * t160, -0.2e1 * t135, 0.2e1 * t130, 0, 0, 0, 0.2e1 * qJ(2) * t101 + 0.2e1 * qJD(2) * t110, -0.2e1 * qJ(2) * t145 + 0.2e1 * qJD(2) * t113, -0.2e1 * t104 * t135 - 0.2e1 * t105 * t136, 0.2e1 * t105 * t131 + 0.4e1 * t113 * t128, -0.2e1 * t110 * t139 - 0.2e1 * t112 * t130, 0.2e1 * t109 * t130 - 0.2e1 * t110 * t137, t94, -0.2e1 * t105 * t112 * t146 + 0.2e1 * t78 * t101 + 0.2e1 * (t29 + t129) * t110, 0.2e1 * t105 * t138 + 0.2e1 * t28 * t110 + 0.2e1 * (-t77 + t92) * t101, 0.2e1 * t23 * t36 + 0.2e1 * t24 * t34 - 0.2e1 * t62 * t9 - 0.2e1 * t64 * t8, 0.2e1 * t23 * t8 + 0.2e1 * t24 * t9 + 0.2e1 * t51 * t82, 0.2e1 * t32 * t11, -0.2e1 * t11 * t31 - 0.2e1 * t13 * t32, 0.2e1 * t101 * t32 + 0.2e1 * t11 * t110, -0.2e1 * t101 * t31 - 0.2e1 * t110 * t13, t94, -0.2e1 * t101 * t125 + 0.2e1 * t2 * t110 + 0.2e1 * t47 * t13 + 0.2e1 * t25 * t31, 0.2e1 * t1 * t110 - 0.2e1 * t101 * t124 + 0.2e1 * t47 * t11 + 0.2e1 * t25 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150 * t148, t150 * t149, -t33 * t64 - t34 * t63 - t35 * t62 - t36 * t61, -t113 * t51 + t145 * t82 + t23 * t33 + t24 * t35 - t61 * t8 - t63 * t9, 0, 0, 0, 0, 0, t12 * t110 - t113 * t13 + (t110 * t31 - t113 * t121) * qJD(3), t10 * t110 - t113 * t11 + (t110 * t32 - t113 * t120) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t33 * t61 - 0.2e1 * t35 * t63 - 0.2e1 * t135, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t101, 0, -t96, -t134, -t113 * t131 - t128, -0.4e1 * t113 * t136 + t152 * t145, t74, -t72, 0 (-t109 * t153 - t127 * t112) * qJD(4) + (t126 * t109 - t92) * qJD(3) (t127 * t109 - t112 * t153) * qJD(4) + (-pkin(8) * t154 + (pkin(3) * t112 + t155) * t110) * qJD(3), -t119 * t9 - t23 * t71 - t24 * t70 + t34 * t50 + t36 * t49 - t37 * t64 - t38 * t62 - t8 * t81, t100 * t82 + t23 * t37 + t24 * t38 + t49 * t8 + t50 * t9 + t51 * t99, t11 * t45 + t21 * t32, -t11 * t44 - t13 * t45 - t21 * t31 - t22 * t32, t101 * t45 + t110 * t21, -t101 * t44 - t110 * t22, 0, -t101 * t123 + t4 * t110 + t54 * t13 + t47 * t22 + t25 * t44 + t52 * t31, -t101 * t122 + t54 * t11 + t3 * t110 + t47 * t21 + t25 * t45 + t52 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t101, 0, 0, 0, 0, 0, t73, -t161, -t119 * t35 - t33 * t81 + t61 * t71 + t63 * t70, -pkin(4) * t139 + t145 * t99 + t33 * t49 + t35 * t50 - t37 * t61 - t38 * t63, 0, 0, 0, 0, 0, -t113 * t22 + t145 * t44, -t113 * t21 + t145 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t136, -0.2e1 * t131, 0, 0, 0, t109 * t143, t112 * t143, -0.2e1 * t119 * t38 - 0.2e1 * t37 * t81 - 0.2e1 * t49 * t71 - 0.2e1 * t50 * t70, 0.2e1 * t100 * t99 + 0.2e1 * t37 * t49 + 0.2e1 * t38 * t50, 0.2e1 * t45 * t21, -0.2e1 * t21 * t44 - 0.2e1 * t22 * t45, 0, 0, 0, 0.2e1 * t22 * t54 + 0.2e1 * t44 * t52, 0.2e1 * t21 * t54 + 0.2e1 * t45 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t161, t101, t29, t28 (t106 * t34 + t107 * t36) * pkin(4) (t106 * t9 + t107 * t8) * pkin(4), 0, 0, t11, -t13, t101, -t101 * t117 - t59 * t110 + t2, -t101 * t118 + t58 * t110 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t72, 0 (t106 * t35 + t107 * t33) * pkin(4), 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, -t149, 0, -pkin(8) * t148, pkin(8) * t149 (-t106 * t70 - t107 * t71) * pkin(4) (t106 * t38 + t107 * t37) * pkin(4), 0, 0, t21, -t22, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t59, 0.2e1 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, t13, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, 0, 0, 0, 0, 0, t22, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t13, t101, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t22, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;

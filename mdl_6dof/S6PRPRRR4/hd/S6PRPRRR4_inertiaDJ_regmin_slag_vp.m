% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:38:56
% EndTime: 2019-03-08 20:38:59
% DurationCPUTime: 1.18s
% Computational Cost: add. (1657->186), mult. (4288->341), div. (0->0), fcn. (4458->12), ass. (0->129)
t77 = sin(pkin(12));
t79 = cos(pkin(12));
t83 = sin(qJ(4));
t87 = cos(qJ(4));
t56 = t87 * t77 + t83 * t79;
t81 = sin(qJ(6));
t82 = sin(qJ(5));
t85 = cos(qJ(6));
t86 = cos(qJ(5));
t58 = t81 * t86 + t85 * t82;
t33 = t58 * t56;
t134 = t87 * t79;
t55 = t83 * t77 - t134;
t70 = -t79 * pkin(3) - pkin(2);
t39 = t55 * pkin(4) - t56 * pkin(9) + t70;
t133 = pkin(8) + qJ(3);
t62 = t133 * t77;
t63 = t133 * t79;
t45 = -t83 * t62 + t87 * t63;
t40 = t86 * t45;
t132 = t82 * t39 + t40;
t78 = sin(pkin(6));
t84 = sin(qJ(2));
t144 = t78 * t84;
t80 = cos(pkin(6));
t48 = -t77 * t144 + t80 * t79;
t49 = t79 * t144 + t80 * t77;
t97 = t87 * t48 - t83 * t49;
t89 = t97 * qJD(4);
t127 = qJD(5) * t82;
t50 = t55 * qJD(4);
t137 = t86 * t50;
t92 = -t56 * t127 - t137;
t76 = t86 ^ 2;
t130 = t82 ^ 2 - t76;
t107 = t130 * qJD(5);
t123 = qJD(5) + qJD(6);
t152 = 0.2e1 * t70;
t151 = pkin(9) + pkin(10);
t51 = t56 * qJD(4);
t150 = t51 * pkin(5);
t149 = t55 * pkin(5);
t126 = qJD(5) * t86;
t135 = t87 * t62;
t25 = qJD(4) * t135 - qJD(3) * t134 + (qJD(3) * t77 + qJD(4) * t63) * t83;
t38 = t51 * pkin(4) + t50 * pkin(9);
t8 = -t39 * t126 + t45 * t127 + t86 * t25 - t82 * t38;
t140 = t82 * t50;
t93 = t56 * t126 - t140;
t7 = -t93 * pkin(10) - t8;
t148 = t85 * t7;
t147 = t56 * t50;
t146 = t56 * t82;
t145 = t56 * t86;
t88 = cos(qJ(2));
t143 = t78 * t88;
t15 = -pkin(10) * t146 + t132;
t142 = t81 * t15;
t141 = t81 * t82;
t139 = t82 * t51;
t138 = t85 * t15;
t136 = t86 * t51;
t131 = t77 ^ 2 + t79 ^ 2;
t129 = qJD(2) * t78;
t128 = qJD(2) * t84;
t125 = qJD(6) * t81;
t124 = qJD(6) * t85;
t122 = -0.2e1 * pkin(4) * qJD(5);
t121 = pkin(5) * t127;
t120 = pkin(5) * t125;
t119 = pkin(5) * t124;
t117 = t78 * t128;
t116 = t88 * t129;
t115 = t82 * t126;
t109 = t82 * t25 + t86 * t38;
t6 = pkin(10) * t137 + t150 + (-t40 + (pkin(10) * t56 - t39) * t82) * qJD(5) + t109;
t114 = t85 * t6 - t81 * t7;
t108 = t86 * t39 - t82 * t45;
t14 = -pkin(10) * t145 + t108 + t149;
t113 = -t14 - t149;
t112 = qJD(5) * t151;
t111 = t131 * t88;
t110 = -0.4e1 * t82 * t145;
t44 = t83 * t63 + t135;
t106 = 0.2e1 * t131 * qJD(3);
t105 = pkin(4) * t50 - pkin(9) * t51;
t104 = pkin(4) * t56 + pkin(9) * t55;
t103 = t85 * t14 - t142;
t102 = t81 * t14 + t138;
t29 = t83 * t48 + t87 * t49;
t23 = -t82 * t143 + t86 * t29;
t94 = t86 * t143 + t82 * t29;
t101 = -t81 * t23 - t85 * t94;
t100 = t85 * t23 - t81 * t94;
t42 = t123 * t141 - t86 * t124 - t85 * t126;
t99 = t42 * t55 - t58 * t51;
t98 = -t48 * t77 + t49 * t79;
t64 = t151 * t82;
t65 = t151 * t86;
t96 = -t85 * t64 - t81 * t65;
t95 = -t81 * t64 + t85 * t65;
t57 = -t85 * t86 + t141;
t91 = t55 * t126 + t139;
t90 = t55 * t88;
t26 = t56 * qJD(3) + t45 * qJD(4);
t72 = -t86 * pkin(5) - pkin(4);
t60 = t86 * t112;
t59 = t82 * t112;
t53 = t56 ^ 2;
t43 = t123 * t58;
t41 = 0.2e1 * t55 * t51;
t37 = -t55 * t127 + t136;
t34 = t57 * t56;
t27 = pkin(5) * t146 + t44;
t21 = -t95 * qJD(6) + t81 * t59 - t85 * t60;
t20 = -t96 * qJD(6) + t85 * t59 + t81 * t60;
t19 = t29 * qJD(4) + t56 * t116;
t17 = -t43 * t55 - t57 * t51;
t16 = t93 * pkin(5) + t26;
t13 = -t125 * t146 + (t123 * t145 - t140) * t85 + t92 * t81;
t12 = -t123 * t33 + t57 * t50;
t11 = -t29 * t126 - t82 * t89 + (t88 * t127 + (t82 * t90 + t84 * t86) * qJD(2)) * t78;
t10 = -t86 * (-t90 * t129 + t89) - t82 * t117 + t94 * qJD(5);
t9 = -t132 * qJD(5) + t109;
t4 = -t100 * qJD(6) + t81 * t10 + t85 * t11;
t3 = -t101 * qJD(6) + t85 * t10 - t81 * t11;
t2 = -t102 * qJD(6) + t114;
t1 = -t103 * qJD(6) - t81 * t6 - t148;
t5 = [0, 0, 0, 0, 0, 0, 0, 0.2e1 * (t98 - t144) * t116, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t117, -t116, -t79 * t117, t77 * t117, t111 * t129, t98 * qJD(3) + (-pkin(2) * t84 + qJ(3) * t111) * t129, 0, 0, 0, 0, 0 (t55 * t128 - t51 * t88) * t78 (t56 * t128 + t50 * t88) * t78, 0, 0, 0, 0, 0, t11 * t55 + t19 * t146 - t51 * t94 - t93 * t97, t10 * t55 + t19 * t145 - t23 * t51 - t92 * t97, 0, 0, 0, 0, 0, t101 * t51 - t13 * t97 + t19 * t33 + t4 * t55, -t100 * t51 - t12 * t97 - t19 * t34 + t3 * t55; 0, 0, 0, 0, 0, 0, t106, qJ(3) * t106, -0.2e1 * t147, 0.2e1 * t50 * t55 - 0.2e1 * t56 * t51, 0, 0, 0, t51 * t152, -t50 * t152, -0.2e1 * t53 * t115 - 0.2e1 * t76 * t147, 0.2e1 * t53 * t107 - t50 * t110, 0.2e1 * t56 * t136 + 0.2e1 * t92 * t55, -0.2e1 * t56 * t139 - 0.2e1 * t93 * t55, t41, 0.2e1 * t108 * t51 + 0.2e1 * t26 * t146 + 0.2e1 * t93 * t44 + 0.2e1 * t9 * t55, -0.2e1 * t132 * t51 + 0.2e1 * t26 * t145 + 0.2e1 * t92 * t44 + 0.2e1 * t8 * t55, -0.2e1 * t34 * t12, -0.2e1 * t12 * t33 + 0.2e1 * t34 * t13, 0.2e1 * t12 * t55 - 0.2e1 * t34 * t51, -0.2e1 * t13 * t55 - 0.2e1 * t33 * t51, t41, 0.2e1 * t103 * t51 + 0.2e1 * t27 * t13 + 0.2e1 * t16 * t33 + 0.2e1 * t2 * t55, 0.2e1 * t1 * t55 - 0.2e1 * t102 * t51 + 0.2e1 * t27 * t12 - 0.2e1 * t16 * t34; 0, 0, 0, 0, 0, 0, 0, t117, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t50, 0, 0, 0, 0, 0, t37, -t91, 0, 0, 0, 0, 0, t17, t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t55 * t116 - t89, 0, 0, 0, 0, 0, -t127 * t97 - t19 * t86, -t126 * t97 + t19 * t82, 0, 0, 0, 0, 0, t19 * t57 - t43 * t97, t19 * t58 + t42 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t51, 0, -t26, t25, -t56 * t107 - t82 * t137, qJD(5) * t110 + t130 * t50, t91, t37, 0, -t26 * t86 + t105 * t82 + (-t104 * t86 + t44 * t82) * qJD(5), t26 * t82 + t105 * t86 + (t104 * t82 + t44 * t86) * qJD(5), t12 * t58 + t34 * t42, -t12 * t57 - t58 * t13 + t42 * t33 + t34 * t43, -t99, t17, 0, t33 * t121 + t72 * t13 + t16 * t57 + t21 * t55 + t27 * t43 + t96 * t51, t72 * t12 - t121 * t34 + t16 * t58 + t20 * t55 - t27 * t42 - t51 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t115, -0.2e1 * t107, 0, 0, 0, t82 * t122, t86 * t122, -0.2e1 * t58 * t42, 0.2e1 * t42 * t57 - 0.2e1 * t58 * t43, 0, 0, 0, 0.2e1 * t57 * t121 + 0.2e1 * t72 * t43, 0.2e1 * t121 * t58 - 0.2e1 * t72 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t10, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, -t93, t51, t9, t8, 0, 0, t12, -t13, t51, t85 * t150 + (t113 * t81 - t138) * qJD(6) + t114, -t148 + (-t6 - t150) * t81 + (t113 * t85 + t142) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, -t126, 0, 0, 0, 0, 0, -t43, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, -t127, 0, -pkin(9) * t126, pkin(9) * t127, 0, 0, -t42, -t43, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t120, -0.2e1 * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t13, t51, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t43, 0, t21, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;

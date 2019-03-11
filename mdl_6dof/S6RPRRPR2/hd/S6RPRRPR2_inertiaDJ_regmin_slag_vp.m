% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:02:44
% EndTime: 2019-03-09 05:02:48
% DurationCPUTime: 1.70s
% Computational Cost: add. (2238->211), mult. (5398->389), div. (0->0), fcn. (4879->10), ass. (0->128)
t111 = sin(qJ(4));
t112 = sin(qJ(3));
t114 = cos(qJ(3));
t146 = t114 * qJD(3);
t113 = cos(qJ(4));
t149 = qJD(4) * t113;
t168 = t111 * t146 + t112 * t149;
t167 = -0.4e1 * t112;
t110 = sin(qJ(6));
t166 = cos(qJ(6));
t107 = sin(pkin(11));
t108 = cos(pkin(11));
t81 = t107 * t113 + t108 * t111;
t62 = t81 * t112;
t156 = t112 * t113;
t158 = t111 * t112;
t63 = -t107 * t158 + t108 * t156;
t35 = t110 * t63 + t166 * t62;
t104 = t112 ^ 2;
t131 = (-t114 ^ 2 + t104) * qJD(3);
t165 = pkin(4) * t107;
t164 = -qJ(5) - pkin(8);
t147 = qJD(5) * t113;
t155 = t113 * t114;
t102 = t112 * qJD(3);
t136 = t111 * t102;
t97 = sin(pkin(10)) * pkin(1) + pkin(7);
t79 = t97 * t136;
t128 = pkin(3) * t112 - pkin(8) * t114;
t86 = t128 * qJD(3);
t162 = t113 * t86 + t79;
t129 = -pkin(3) * t114 - t112 * pkin(8);
t99 = -cos(pkin(10)) * pkin(1) - pkin(2);
t78 = t129 + t99;
t85 = t97 * t155;
t17 = -t112 * t147 + (pkin(4) * t112 - qJ(5) * t155) * qJD(3) + (-t85 + (qJ(5) * t112 - t78) * t111) * qJD(4) + t162;
t163 = -t111 * t86 - t78 * t149;
t21 = (-qJ(5) * qJD(4) - qJD(3) * t97) * t156 + (-qJD(5) * t112 + (-qJ(5) * qJD(3) - qJD(4) * t97) * t114) * t111 - t163;
t9 = t107 * t17 + t108 * t21;
t159 = t111 * t97;
t65 = t113 * t78;
t38 = -qJ(5) * t156 + t65 + (-pkin(4) - t159) * t114;
t64 = t111 * t78;
t161 = t85 + t64;
t44 = -qJ(5) * t158 + t161;
t20 = t107 * t38 + t108 * t44;
t133 = qJD(4) * t164;
t66 = t111 * t133 + t147;
t67 = -qJD(5) * t111 + t113 * t133;
t40 = t107 * t67 + t108 * t66;
t88 = t164 * t111;
t89 = t164 * t113;
t50 = t107 * t88 - t108 * t89;
t68 = pkin(4) * t158 + t112 * t97;
t157 = t111 * t114;
t105 = t113 ^ 2;
t154 = t111 ^ 2 - t105;
t152 = qJD(4) * t104;
t151 = qJD(4) * t111;
t150 = qJD(4) * t112;
t148 = qJD(4) * t114;
t145 = -0.2e1 * pkin(3) * qJD(4);
t144 = 0.2e1 * qJD(3) * t99;
t87 = t97 * t146;
t48 = t168 * pkin(4) + t87;
t101 = pkin(4) * t151;
t142 = t97 * t152;
t100 = -pkin(4) * t113 - pkin(3);
t140 = t113 * t146;
t139 = t111 * t148;
t137 = t113 * t148;
t135 = t111 * t149;
t134 = t112 * t146;
t8 = -t107 * t21 + t108 * t17;
t19 = -t107 * t44 + t108 * t38;
t39 = -t107 * t66 + t108 * t67;
t49 = t107 * t89 + t108 * t88;
t132 = t154 * qJD(4);
t130 = t111 * t140;
t127 = t107 * t111 - t108 * t113;
t13 = -pkin(5) * t114 - t63 * pkin(9) + t19;
t14 = -pkin(9) * t62 + t20;
t126 = t110 * t14 - t166 * t13;
t125 = t110 * t13 + t166 * t14;
t41 = -pkin(9) * t81 + t49;
t42 = -pkin(9) * t127 + t50;
t124 = t110 * t42 - t166 * t41;
t123 = t110 * t41 + t166 * t42;
t36 = -t110 * t62 + t166 * t63;
t122 = -t110 * t81 - t127 * t166;
t47 = -t110 * t127 + t166 * t81;
t70 = t81 * qJD(4);
t71 = -t107 * t151 + t108 * t149;
t23 = t47 * qJD(6) + t110 * t71 + t166 * t70;
t121 = -t102 * t122 - t114 * t23;
t120 = t81 * t114;
t119 = t127 * t112;
t98 = pkin(4) * t108 + pkin(5);
t118 = t110 * t165 - t166 * t98;
t117 = t110 * t98 + t166 * t165;
t72 = t111 * t150 - t140;
t73 = t113 * t102 + t139;
t116 = qJD(4) * t119;
t37 = t127 * t146 + t81 * t150;
t4 = pkin(5) * t102 + pkin(9) * t37 + t8;
t115 = -qJD(3) * t120 + t116;
t5 = pkin(9) * t115 + t9;
t2 = -qJD(6) * t125 - t110 * t5 + t166 * t4;
t1 = t126 * qJD(6) - t110 * t4 - t166 * t5;
t94 = -0.2e1 * t134;
t75 = t136 - t137;
t59 = t117 * qJD(6);
t58 = t118 * qJD(6);
t54 = pkin(5) * t127 + t100;
t51 = pkin(5) * t70 + t101;
t45 = pkin(5) * t62 + t68;
t28 = -pkin(9) * t70 + t40;
t27 = -pkin(9) * t71 + t39;
t26 = -t161 * qJD(4) + t162;
t25 = t73 * t97 + t163;
t24 = -pkin(5) * t115 + t48;
t22 = t122 * qJD(6) - t110 * t70 + t166 * t71;
t12 = t47 * t102 - t114 * t22;
t11 = t36 * qJD(6) - t110 * t37 - t166 * t115;
t10 = t35 * qJD(6) - t110 * t115 + t166 * t37;
t7 = -t123 * qJD(6) - t110 * t28 + t166 * t27;
t6 = t124 * qJD(6) - t110 * t27 - t166 * t28;
t3 = [0, 0, 0, 0, 0.2e1 * t134, -0.2e1 * t131, 0, 0, 0, t112 * t144, t114 * t144, -0.2e1 * t104 * t135 + 0.2e1 * t105 * t134, t130 * t167 + 0.2e1 * t154 * t152, 0.2e1 * t112 * t139 + 0.2e1 * t113 * t131, -0.2e1 * t111 * t131 + 0.2e1 * t112 * t137, t94, 0.2e1 * t113 * t142 + 0.2e1 * t65 * t102 + 0.2e1 * (-t26 + t79) * t114, -0.2e1 * t111 * t142 - 0.2e1 * t25 * t114 + 0.2e1 * (-t64 + t85) * t102, 0.2e1 * t115 * t20 + 0.2e1 * t19 * t37 - 0.2e1 * t9 * t62 - 0.2e1 * t8 * t63, 0.2e1 * t19 * t8 + 0.2e1 * t20 * t9 + 0.2e1 * t48 * t68, -0.2e1 * t36 * t10, 0.2e1 * t10 * t35 - 0.2e1 * t11 * t36, 0.2e1 * t10 * t114 + 0.2e1 * t36 * t102, -0.2e1 * t35 * t102 + 0.2e1 * t11 * t114, t94, -0.2e1 * t126 * t102 + 0.2e1 * t45 * t11 - 0.2e1 * t2 * t114 + 0.2e1 * t24 * t35, -0.2e1 * t1 * t114 - 0.2e1 * t45 * t10 - 0.2e1 * t125 * t102 + 0.2e1 * t24 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t114 - t20 * t37 - t8 * t62 + t9 * t63 + t19 * t116 + (t68 * t112 - t120 * t19) * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t63 * t37 - 0.2e1 * t62 * t116 + 0.2e1 * (t62 * t81 - t112) * t146, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t146, -t102, 0, -t87, t97 * t102, -t112 * t132 + t130, t135 * t167 - t154 * t146, t75, t73, 0 (pkin(8) * t155 + (-pkin(3) * t113 + t159) * t112) * qJD(4) + (t111 * t129 - t85) * qJD(3) (t111 * t128 + t97 * t156) * qJD(4) + (t113 * t129 + t97 * t157) * qJD(3), t115 * t50 - t127 * t9 - t19 * t71 - t20 * t70 + t49 * t37 - t39 * t63 - t40 * t62 - t8 * t81, t100 * t48 + t68 * t101 + t19 * t39 + t20 * t40 + t49 * t8 + t50 * t9, -t10 * t47 + t22 * t36, -t10 * t122 - t11 * t47 - t22 * t35 - t23 * t36, t12, -t121, 0, -t124 * t102 + t54 * t11 - t7 * t114 - t122 * t24 + t45 * t23 + t51 * t35, -t54 * t10 - t123 * t102 - t6 * t114 + t45 * t22 + t24 * t47 + t51 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t146, 0, 0, 0, 0, 0, -t73, t75, -t115 * t81 + t127 * t37 + t62 * t71 - t63 * t70, -t37 * t50 - t62 * t39 + t63 * t40 + (-pkin(4) * t157 + t119 * t49) * qJD(4) + (t112 * t100 - t120 * t49) * qJD(3), 0, 0, 0, 0, 0, t121, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t135, -0.2e1 * t132, 0, 0, 0, t111 * t145, t113 * t145, -0.2e1 * t127 * t40 - 0.2e1 * t39 * t81 - 0.2e1 * t49 * t71 - 0.2e1 * t50 * t70, 0.2e1 * t100 * t101 + 0.2e1 * t39 * t49 + 0.2e1 * t40 * t50, 0.2e1 * t47 * t22, 0.2e1 * t122 * t22 - 0.2e1 * t23 * t47, 0, 0, 0, -0.2e1 * t122 * t51 + 0.2e1 * t23 * t54, 0.2e1 * t22 * t54 + 0.2e1 * t47 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t168, t102, t26, t25 (t108 * t37 + (t107 * t72 - t108 * t168) * t107) * pkin(4) (t107 * t9 + t108 * t8) * pkin(4), 0, 0, -t10, -t11, t102, -t118 * t102 + t59 * t114 + t2, -t117 * t102 - t58 * t114 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t168, t72, 0 (-t168 * t108 ^ 2 + (t108 * t72 - t37) * t107) * pkin(4), 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, -t151, 0, -pkin(8) * t149, pkin(8) * t151 (-t107 * t70 - t108 * t71) * pkin(4) (t107 * t40 + t108 * t39) * pkin(4), 0, 0, t22, -t23, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t59, 0.2e1 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, t11, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, 0, 0, 0, 0, t23, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, t102, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t23, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;

% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:07
% EndTime: 2019-03-09 10:19:13
% DurationCPUTime: 1.85s
% Computational Cost: add. (3962->219), mult. (9015->415), div. (0->0), fcn. (8959->10), ass. (0->134)
t118 = sin(pkin(10));
t120 = cos(pkin(10));
t123 = sin(qJ(2));
t125 = cos(qJ(2));
t102 = t118 * t125 + t120 * t123;
t124 = cos(qJ(4));
t155 = qJD(4) * t124;
t149 = t102 * t155;
t122 = sin(qJ(4));
t158 = qJD(2) * t125;
t159 = qJD(2) * t123;
t95 = -t118 * t159 + t120 * t158;
t171 = t122 * t95;
t178 = t149 + t171;
t117 = sin(pkin(11));
t119 = cos(pkin(11));
t101 = t117 * t124 + t119 * t122;
t121 = sin(qJ(6));
t175 = cos(qJ(6));
t177 = -t117 * t122 + t119 * t124;
t62 = t121 * t101 - t175 * t177;
t176 = 0.2e1 * qJD(4);
t174 = t117 * pkin(4);
t173 = -qJ(3) - pkin(7);
t139 = -qJ(5) * t95 - qJD(5) * t102;
t145 = qJD(2) * t173;
t90 = t125 * qJD(3) + t123 * t145;
t91 = -t123 * qJD(3) + t125 * t145;
t51 = t118 * t91 + t120 * t90;
t114 = pkin(2) * t159;
t93 = t102 * qJD(2);
t52 = t93 * pkin(3) - t95 * pkin(8) + t114;
t147 = -t122 * t51 + t124 * t52;
t100 = t118 * t123 - t120 * t125;
t150 = -t125 * pkin(2) - pkin(1);
t61 = t100 * pkin(3) - t102 * pkin(8) + t150;
t105 = t173 * t123;
t106 = t173 * t125;
t67 = t118 * t105 - t120 * t106;
t65 = t124 * t67;
t16 = t93 * pkin(4) + t139 * t124 + (-t65 + (qJ(5) * t102 - t61) * t122) * qJD(4) + t147;
t153 = t122 * t52 + t124 * t51 + t61 * t155;
t18 = -qJ(5) * t149 + (-qJD(4) * t67 + t139) * t122 + t153;
t7 = t117 * t16 + t119 * t18;
t146 = -t122 * t67 + t124 * t61;
t166 = t102 * t124;
t27 = t100 * pkin(4) - qJ(5) * t166 + t146;
t167 = t102 * t122;
t172 = t122 * t61 + t65;
t30 = -qJ(5) * t167 + t172;
t20 = t117 * t27 + t119 * t30;
t110 = t118 * pkin(2) + pkin(8);
t161 = qJ(5) + t110;
t143 = qJD(4) * t161;
t74 = t124 * qJD(5) - t122 * t143;
t75 = -t122 * qJD(5) - t124 * t143;
t42 = t117 * t75 + t119 * t74;
t96 = t161 * t122;
t97 = t161 * t124;
t60 = -t117 * t96 + t119 * t97;
t170 = t124 * t95;
t50 = t118 * t90 - t120 * t91;
t169 = t50 * t122;
t168 = t50 * t124;
t162 = t122 * t124;
t116 = t124 ^ 2;
t160 = t122 ^ 2 - t116;
t157 = qJD(4) * t102;
t156 = qJD(4) * t122;
t154 = -0.2e1 * pkin(1) * qJD(2);
t112 = -t120 * pkin(2) - pkin(3);
t152 = t112 * t176;
t113 = pkin(4) * t156;
t148 = t122 * t155;
t6 = -t117 * t18 + t119 * t16;
t19 = -t117 * t30 + t119 * t27;
t41 = -t117 * t74 + t119 * t75;
t59 = -t117 * t97 - t119 * t96;
t144 = -0.4e1 * t102 * t162;
t66 = -t120 * t105 - t118 * t106;
t142 = t160 * qJD(4);
t48 = pkin(4) * t167 + t66;
t128 = t101 * qJD(4);
t94 = -t117 * t156 + t119 * t155;
t37 = qJD(6) * t62 + t121 * t128 - t175 * t94;
t63 = t175 * t101 + t121 * t177;
t141 = t37 * t100 - t63 * t93;
t140 = -t110 * t93 + t112 * t95;
t36 = t178 * pkin(4) + t50;
t138 = t100 * t110 - t102 * t112;
t104 = -t124 * pkin(4) + t112;
t55 = t177 * t102;
t12 = t100 * pkin(5) - t55 * pkin(9) + t19;
t54 = t101 * t102;
t15 = -t54 * pkin(9) + t20;
t136 = -t175 * t12 + t121 * t15;
t135 = t121 * t12 + t175 * t15;
t43 = -t101 * pkin(9) + t59;
t44 = pkin(9) * t177 + t60;
t134 = t121 * t44 - t175 * t43;
t133 = t121 * t43 + t175 * t44;
t132 = -t121 * t55 - t175 * t54;
t34 = -t121 * t54 + t175 * t55;
t131 = t100 * t155 + t122 * t93;
t129 = -t102 * t156 + t170;
t111 = t119 * pkin(4) + pkin(5);
t127 = -t175 * t111 + t121 * t174;
t126 = t121 * t111 + t175 * t174;
t32 = t101 * t157 - t177 * t95;
t4 = t93 * pkin(5) + t32 * pkin(9) + t6;
t31 = -t101 * t95 - t157 * t177;
t5 = t31 * pkin(9) + t7;
t2 = -t135 * qJD(6) - t121 * t5 + t175 * t4;
t1 = t136 * qJD(6) - t121 * t4 - t175 * t5;
t98 = t102 ^ 2;
t79 = t126 * qJD(6);
t78 = t127 * qJD(6);
t72 = pkin(5) * t128 + t113;
t71 = -pkin(5) * t177 + t104;
t64 = 0.2e1 * t100 * t93;
t58 = -t100 * t156 + t124 * t93;
t40 = -pkin(9) * t128 + t42;
t39 = -t94 * pkin(9) + t41;
t38 = t63 * qJD(6) + t121 * t94 + t175 * t128;
t35 = t54 * pkin(5) + t48;
t24 = -t38 * t100 - t62 * t93;
t23 = -t172 * qJD(4) + t147;
t22 = t67 * t156 - t153;
t21 = -t31 * pkin(5) + t36;
t11 = -t133 * qJD(6) - t121 * t40 + t175 * t39;
t10 = t134 * qJD(6) - t121 * t39 - t175 * t40;
t9 = t34 * qJD(6) - t121 * t32 - t175 * t31;
t8 = t132 * qJD(6) + t121 * t31 - t175 * t32;
t3 = [0, 0, 0, 0.2e1 * t123 * t158, 0.2e1 * (-t123 ^ 2 + t125 ^ 2) * qJD(2), 0, 0, 0, t123 * t154, t125 * t154, -0.2e1 * t51 * t100 + 0.2e1 * t50 * t102 + 0.2e1 * t66 * t95 - 0.2e1 * t67 * t93, 0.2e1 * t150 * t114 + 0.2e1 * t66 * t50 + 0.2e1 * t67 * t51, 0.2e1 * t116 * t102 * t95 - 0.2e1 * t98 * t148, t160 * t98 * t176 + t95 * t144, 0.2e1 * t129 * t100 + 0.2e1 * t93 * t166, -0.2e1 * t100 * t178 - 0.2e1 * t93 * t167, t64, 0.2e1 * t23 * t100 + 0.2e1 * t146 * t93 + 0.2e1 * t66 * t171 + 0.2e1 * (t66 * t155 + t169) * t102, 0.2e1 * t22 * t100 - 0.2e1 * t172 * t93 + 0.2e1 * t66 * t170 + 0.2e1 * (-t66 * t156 + t168) * t102, 0.2e1 * t19 * t32 + 0.2e1 * t20 * t31 - 0.2e1 * t7 * t54 - 0.2e1 * t6 * t55, 0.2e1 * t19 * t6 + 0.2e1 * t20 * t7 + 0.2e1 * t48 * t36, 0.2e1 * t34 * t8, 0.2e1 * t132 * t8 - 0.2e1 * t34 * t9, 0.2e1 * t8 * t100 + 0.2e1 * t34 * t93, -0.2e1 * t9 * t100 + 0.2e1 * t132 * t93, t64, 0.2e1 * t2 * t100 - 0.2e1 * t132 * t21 - 0.2e1 * t136 * t93 + 0.2e1 * t35 * t9, 0.2e1 * t1 * t100 - 0.2e1 * t135 * t93 + 0.2e1 * t21 * t34 + 0.2e1 * t35 * t8; 0, 0, 0, 0, 0, t158, -t159, 0, -pkin(7) * t158, pkin(7) * t159 (-t118 * t93 - t120 * t95) * pkin(2) (t118 * t51 - t120 * t50) * pkin(2), -t102 * t142 + t95 * t162, qJD(4) * t144 - t160 * t95, t131, t58, 0, -t168 + t140 * t122 + (t122 * t66 - t138 * t124) * qJD(4), t169 + t140 * t124 + (t138 * t122 + t124 * t66) * qJD(4), -t6 * t101 - t128 * t20 + t177 * t7 - t19 * t94 + t60 * t31 + t59 * t32 - t41 * t55 - t42 * t54, t36 * t104 + t113 * t48 + t19 * t41 + t20 * t42 + t6 * t59 + t7 * t60, -t34 * t37 + t8 * t63, -t132 * t37 - t34 * t38 - t8 * t62 - t63 * t9, -t141, t24, 0, t11 * t100 - t132 * t72 - t134 * t93 + t21 * t62 + t35 * t38 + t71 * t9, t10 * t100 - t133 * t93 + t21 * t63 + t72 * t34 - t35 * t37 + t71 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t148, -0.2e1 * t142, 0, 0, 0, t122 * t152, t124 * t152, -0.2e1 * t41 * t101 - 0.2e1 * t128 * t60 + 0.2e1 * t177 * t42 - 0.2e1 * t59 * t94, 0.2e1 * t104 * t113 + 0.2e1 * t59 * t41 + 0.2e1 * t60 * t42, -0.2e1 * t63 * t37, 0.2e1 * t37 * t62 - 0.2e1 * t63 * t38, 0, 0, 0, 0.2e1 * t71 * t38 + 0.2e1 * t72 * t62, -0.2e1 * t71 * t37 + 0.2e1 * t72 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, 0, 0, 0, 0, t58, -t131, t101 * t31 + t128 * t55 + t177 * t32 - t94 * t54, t7 * t101 - t128 * t19 + t177 * t6 + t20 * t94, 0, 0, 0, 0, 0, t24, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t101 - t128 * t59 + t177 * t41 + t60 * t94, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t101 * t94 - 0.2e1 * t128 * t177, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, -t178, t93, t23, t22 (t117 * t31 + t119 * t32) * pkin(4) (t117 * t7 + t119 * t6) * pkin(4), 0, 0, t8, -t9, t93, -t79 * t100 - t127 * t93 + t2, t78 * t100 - t126 * t93 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t156, 0, -t110 * t155, t110 * t156 (-t117 * t128 - t119 * t94) * pkin(4) (t117 * t42 + t119 * t41) * pkin(4), 0, 0, -t37, -t38, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t155, 0 (t94 * t117 - t119 * t128) * pkin(4), 0, 0, 0, 0, 0, -t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t79, 0.2e1 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, 0, 0, 0, 0, t38, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, t93, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t38, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;

% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRPR11_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:03
% EndTime: 2019-03-09 11:15:09
% DurationCPUTime: 1.77s
% Computational Cost: add. (2396->204), mult. (5357->376), div. (0->0), fcn. (4689->8), ass. (0->133)
t125 = cos(qJ(2));
t110 = qJD(2) * t125;
t122 = sin(qJ(2));
t120 = sin(qJ(6));
t123 = cos(qJ(6));
t118 = sin(pkin(10));
t119 = cos(pkin(10));
t121 = sin(qJ(4));
t124 = cos(qJ(4));
t136 = t118 * t121 - t119 * t124;
t137 = t118 * t124 + t119 * t121;
t186 = -t120 * t137 - t123 * t136;
t169 = qJD(4) * t124;
t156 = t119 * t169;
t170 = qJD(4) * t121;
t70 = t118 * t170 - t156;
t71 = -t118 * t169 - t119 * t170;
t23 = qJD(6) * t186 + t120 * t71 - t123 * t70;
t45 = -t120 * t136 + t123 * t137;
t188 = -t110 * t45 - t122 * t23;
t182 = pkin(3) + pkin(7);
t187 = (t118 * t70 - t119 * t71) * pkin(4);
t168 = qJD(4) * t125;
t158 = t121 * t168;
t171 = qJD(2) * t124;
t159 = t122 * t171;
t132 = t158 + t159;
t117 = t125 ^ 2;
t152 = qJD(2) * (t122 ^ 2 - t117);
t114 = t121 ^ 2;
t174 = -t124 ^ 2 + t114;
t151 = t174 * qJD(4);
t127 = -qJD(6) * t45 + t120 * t70 + t123 * t71;
t185 = t110 * t186 + t122 * t127;
t126 = -pkin(2) - pkin(8);
t177 = t122 * qJ(3);
t138 = -t125 * t126 + t177;
t165 = t125 * qJD(3);
t95 = t182 * t125;
t184 = qJD(2) * t138 - qJD(4) * t95 - t165;
t183 = 0.2e1 * qJD(3);
t181 = pkin(4) * t118;
t80 = -pkin(1) - t138;
t153 = qJ(5) * t125 - t80;
t172 = qJD(2) * t122;
t150 = pkin(2) * t172 - t122 * qJD(3);
t178 = qJ(3) * t125;
t52 = (pkin(8) * t122 - t178) * qJD(2) + t150;
t107 = pkin(7) * t110;
t87 = pkin(3) * t110 + t107;
t76 = t124 * t87;
t94 = t182 * t122;
t16 = pkin(4) * t110 + t76 + t153 * t169 + (-qJ(5) * t172 - qJD(4) * t94 + qJD(5) * t125 - t52) * t121;
t166 = t124 * qJD(5);
t25 = -t121 * t87 - t124 * t52 - t169 * t94 + t170 * t80;
t18 = qJ(5) * t132 - t125 * t166 - t25;
t9 = t118 * t16 + t119 * t18;
t85 = t124 * t94;
t39 = t122 * pkin(4) + t121 * t153 + t85;
t176 = t124 * t125;
t180 = t121 * t94 + t124 * t80;
t43 = -qJ(5) * t176 + t180;
t20 = t118 * t39 + t119 * t43;
t175 = qJ(5) - t126;
t62 = t170 * t175 - t166;
t89 = t175 * t124;
t63 = -qJD(4) * t89 - t121 * qJD(5);
t34 = t118 * t62 + t119 * t63;
t88 = t175 * t121;
t48 = -t118 * t89 - t119 * t88;
t106 = pkin(4) * t121 + qJ(3);
t167 = qJD(4) * t126;
t97 = pkin(4) * t169 + qJD(3);
t164 = -0.2e1 * pkin(1) * qJD(2);
t163 = pkin(7) * t172;
t72 = pkin(4) * t176 + t95;
t160 = t121 * t172;
t157 = t124 * t168;
t155 = t122 * t110;
t154 = t121 * t169;
t8 = -t118 * t18 + t119 * t16;
t19 = -t118 * t43 + t119 * t39;
t33 = -t118 * t63 + t119 * t62;
t47 = t118 * t88 - t119 * t89;
t149 = t121 * t159;
t148 = -t136 * t71 - t137 * t70;
t147 = -pkin(2) * t125 - t177;
t65 = t137 * t125;
t12 = pkin(5) * t122 + pkin(9) * t65 + t19;
t64 = t136 * t125;
t13 = pkin(9) * t64 + t20;
t145 = t12 * t123 - t120 * t13;
t144 = t12 * t120 + t123 * t13;
t35 = pkin(9) * t136 + t47;
t36 = -pkin(9) * t137 + t48;
t143 = t120 * t36 - t123 * t35;
t142 = t120 * t35 + t123 * t36;
t141 = t120 * t65 + t123 * t64;
t38 = t120 * t64 - t123 * t65;
t105 = pkin(4) * t119 + pkin(5);
t134 = t105 * t120 + t123 * t181;
t133 = -t105 * t123 + t120 * t181;
t86 = t182 * t172;
t131 = -t86 + (-t122 * t126 - t178) * qJD(4);
t130 = -t136 * t8 + t137 * t9 + t19 * t71 - t20 * t70;
t129 = -t136 * t33 + t137 * t34 + t47 * t71 - t48 * t70;
t128 = qJD(2) * t147 + t165;
t41 = -t118 * t132 - t119 * t160 + t125 * t156;
t4 = pkin(5) * t110 + pkin(9) * t41 + t8;
t40 = qJD(4) * t65 - t136 * t172;
t5 = pkin(9) * t40 + t9;
t2 = -qJD(6) * t144 - t120 * t5 + t123 * t4;
t49 = -pkin(4) * t158 + (-pkin(4) * t124 - t182) * t172;
t1 = -qJD(6) * t145 - t120 * t4 - t123 * t5;
t96 = 0.2e1 * t155;
t90 = -pkin(1) + t147;
t74 = -t110 * t121 - t122 * t169;
t73 = t110 * t124 - t122 * t170;
t67 = -qJ(3) * t110 + t150;
t59 = t134 * qJD(6);
t58 = t133 * qJD(6);
t57 = pkin(5) * t137 + t106;
t51 = -pkin(5) * t70 + t97;
t44 = -pkin(5) * t64 + t72;
t29 = pkin(9) * t70 + t34;
t28 = -pkin(9) * t71 + t33;
t27 = -pkin(5) * t40 + t49;
t26 = -qJD(4) * t180 - t121 * t52 + t76;
t11 = qJD(6) * t38 - t120 * t41 - t123 * t40;
t10 = qJD(6) * t141 + t120 * t40 - t123 * t41;
t7 = -qJD(6) * t142 - t120 * t29 + t123 * t28;
t6 = qJD(6) * t143 - t120 * t28 - t123 * t29;
t3 = [0, 0, 0, t96, -0.2e1 * t152, 0, 0, 0, t122 * t164, t125 * t164, 0, 0.2e1 * t125 * t67 - 0.2e1 * t172 * t90, -0.2e1 * t110 * t90 - 0.2e1 * t122 * t67, 0.2e1 * t90 * t67, -0.2e1 * t114 * t155 + 0.2e1 * t117 * t154, -0.2e1 * t117 * t151 - 0.4e1 * t125 * t149, 0.2e1 * t121 * t152 - 0.2e1 * t122 * t157, 0.2e1 * t122 * t158 + 0.2e1 * t124 * t152, t96, 0.2e1 * (-t171 * t95 + t26) * t122 + 0.2e1 * ((-t121 * t80 + t85) * qJD(2) - t86 * t124 - t95 * t170) * t125, 0.2e1 * (qJD(2) * t121 * t95 + t25) * t122 + 0.2e1 * (-qJD(2) * t180 + t86 * t121 - t169 * t95) * t125, 0.2e1 * t19 * t41 + 0.2e1 * t20 * t40 + 0.2e1 * t64 * t9 + 0.2e1 * t65 * t8, 0.2e1 * t19 * t8 + 0.2e1 * t20 * t9 + 0.2e1 * t49 * t72, 0.2e1 * t38 * t10, 0.2e1 * t10 * t141 - 0.2e1 * t11 * t38, 0.2e1 * t10 * t122 + 0.2e1 * t110 * t38, -0.2e1 * t11 * t122 + 0.2e1 * t110 * t141, t96, 0.2e1 * t11 * t44 + 0.2e1 * t110 * t145 + 0.2e1 * t122 * t2 - 0.2e1 * t141 * t27, 0.2e1 * t1 * t122 + 0.2e1 * t10 * t44 - 0.2e1 * t110 * t144 + 0.2e1 * t27 * t38; 0, 0, 0, 0, 0, t110, -t172, 0, -t107, t163, t128, t107, -t163, t128 * pkin(7), t125 * t151 + t149, 0.4e1 * t125 * t154 - t172 * t174, t73, t74, 0, t121 * t131 - t124 * t184, t121 * t184 + t124 * t131, t33 * t65 + t34 * t64 + t40 * t48 + t41 * t47 - t130, t106 * t49 + t19 * t33 + t20 * t34 + t47 * t8 + t48 * t9 + t72 * t97, t10 * t186 + t127 * t38, -t10 * t45 - t11 * t186 + t127 * t141 - t23 * t38, t185, t188, 0, t11 * t57 - t110 * t143 + t122 * t7 - t141 * t51 + t23 * t44 + t27 * t45, t10 * t57 - t110 * t142 + t122 * t6 + t127 * t44 + t186 * t27 + t38 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, qJ(3) * t183, -0.2e1 * t154, 0.2e1 * t151, 0, 0, 0, 0.2e1 * qJ(3) * t169 + 0.2e1 * qJD(3) * t121, -0.2e1 * qJ(3) * t170 + 0.2e1 * qJD(3) * t124, -0.2e1 * t129, 0.2e1 * t106 * t97 + 0.2e1 * t33 * t47 + 0.2e1 * t34 * t48, 0.2e1 * t186 * t127, -0.2e1 * t127 * t45 - 0.2e1 * t186 * t23, 0, 0, 0, 0.2e1 * t23 * t57 + 0.2e1 * t45 * t51, 0.2e1 * t127 * t57 + 0.2e1 * t186 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, 0, 0, t107, 0, 0, 0, 0, 0, t73, t74, -t136 * t41 + t137 * t40 - t64 * t70 + t65 * t71, t130, 0, 0, 0, 0, 0, t185, t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t148, t129, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t148, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157 + t160, t132, t110, t26, t25 (t118 * t40 + t119 * t41) * pkin(4) (t118 * t9 + t119 * t8) * pkin(4), 0, 0, t10, -t11, t110, -t110 * t133 - t59 * t122 + t2, -t110 * t134 + t58 * t122 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, -t169, 0, -t121 * t167, -t124 * t167, t187 (t118 * t34 + t119 * t33) * pkin(4), 0, 0, t127, -t23, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, -t169, 0, -t187, 0, 0, 0, 0, 0, t127, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t59, 0.2e1 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, t23, t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t110, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, -t23, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;

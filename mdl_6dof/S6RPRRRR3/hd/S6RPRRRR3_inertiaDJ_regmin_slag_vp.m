% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:30
% EndTime: 2019-03-09 07:02:37
% DurationCPUTime: 2.24s
% Computational Cost: add. (2907->243), mult. (6979->417), div. (0->0), fcn. (6382->10), ass. (0->148)
t116 = sin(qJ(3));
t115 = sin(qJ(4));
t120 = cos(qJ(3));
t169 = qJD(3) * t120;
t146 = t115 * t169;
t119 = cos(qJ(4));
t167 = qJD(4) * t119;
t191 = t116 * t167 + t146;
t190 = -0.4e1 * t116;
t110 = t119 ^ 2;
t172 = t115 ^ 2 - t110;
t134 = t172 * qJD(4);
t161 = qJD(4) + qJD(5);
t102 = sin(pkin(11)) * pkin(1) + pkin(7);
t166 = qJD(4) * t120;
t149 = t115 * t166;
t170 = qJD(3) * t119;
t67 = t116 * t170 + t149;
t103 = -cos(pkin(11)) * pkin(1) - pkin(2);
t131 = -t120 * pkin(3) - t116 * pkin(8);
t75 = t103 + t131;
t130 = pkin(3) * t116 - pkin(8) * t120;
t86 = t130 * qJD(3);
t29 = t67 * t102 - t115 * t86 - t75 * t167;
t189 = pkin(8) + pkin(9);
t118 = cos(qJ(5));
t188 = t118 * pkin(4);
t187 = t120 * pkin(5);
t114 = sin(qJ(5));
t175 = t116 * t119;
t179 = t102 * t115;
t63 = t119 * t75;
t41 = -pkin(9) * t175 + t63 + (-pkin(4) - t179) * t120;
t176 = t115 * t116;
t62 = t115 * t75;
t173 = t119 * t120;
t83 = t102 * t173;
t182 = t83 + t62;
t46 = -pkin(9) * t176 + t182;
t42 = t118 * t46;
t186 = t114 * t41 + t42;
t107 = qJD(3) * t116;
t143 = t102 * t107;
t184 = t115 * t143 + t119 * t86;
t94 = t189 * t115;
t95 = t189 * t119;
t183 = -t114 * t94 + t118 * t95;
t174 = t118 * t115;
t181 = t114 * t175 + t116 * t174;
t117 = cos(qJ(6));
t15 = -t181 * pkin(10) + t186;
t180 = t117 * t15;
t64 = pkin(4) * t176 + t116 * t102;
t113 = sin(qJ(6));
t178 = t113 * t114;
t177 = t114 * t117;
t109 = t116 ^ 2;
t171 = -t120 ^ 2 + t109;
t168 = qJD(4) * t115;
t165 = qJD(5) * t114;
t164 = qJD(5) * t118;
t163 = qJD(6) * t113;
t162 = qJD(6) * t117;
t160 = -0.2e1 * pkin(3) * qJD(4);
t23 = (pkin(4) * t116 - pkin(9) * t173) * qJD(3) + (-t83 + (pkin(9) * t116 - t75) * t115) * qJD(4) + t184;
t25 = -pkin(9) * t191 - t29;
t159 = t114 * t23 + t118 * t25 + t41 * t164;
t87 = t102 * t169;
t55 = pkin(4) * t191 + t87;
t158 = 0.2e1 * qJD(3) * t103;
t157 = pkin(5) * t107;
t156 = pkin(4) * t168;
t155 = pkin(4) * t165;
t154 = pkin(4) * t164;
t153 = pkin(5) * t163;
t152 = pkin(5) * t162;
t106 = -t119 * pkin(4) - pkin(3);
t142 = t119 * t169;
t79 = t114 * t119 + t174;
t54 = t161 * t79;
t31 = t114 * t146 + t54 * t116 - t118 * t142;
t138 = -t114 * t25 + t118 * t23;
t9 = -t186 * qJD(5) + t138;
t4 = t31 * pkin(10) + t157 + t9;
t129 = t161 * t176;
t135 = t114 * t142 + t118 * t191 + t164 * t175;
t5 = -t135 * pkin(10) + (pkin(10) * t129 - qJD(5) * t46) * t114 + t159;
t151 = -t113 * t5 + t117 * t4;
t150 = qJD(4) * t189;
t147 = t119 * t166;
t145 = t115 * t167;
t144 = t116 * t169;
t38 = t118 * t41;
t78 = t114 * t115 - t118 * t119;
t61 = t78 * t116;
t14 = t61 * pkin(10) - t114 * t46 - t187 + t38;
t141 = -t14 + t187;
t13 = t15 * t163;
t140 = -t113 * t4 + t13;
t139 = t117 * t181;
t137 = -t114 * t95 - t118 * t94;
t105 = pkin(5) + t188;
t136 = qJD(6) * (-pkin(5) - t105);
t133 = t171 * qJD(3);
t132 = t115 * t142;
t128 = t113 * t14 + t180;
t43 = -t79 * pkin(10) + t137;
t44 = -t78 * pkin(10) + t183;
t127 = t113 * t44 - t117 * t43;
t126 = t113 * t43 + t117 * t44;
t51 = t113 * t79 + t117 * t78;
t52 = -t113 * t78 + t117 * t79;
t53 = t161 * t78;
t17 = t52 * qJD(6) - t113 * t53 + t117 * t54;
t125 = t51 * t107 - t120 * t17;
t124 = t78 * t107 - t120 * t54;
t8 = t46 * t165 - t159;
t84 = t115 * t150;
t85 = t119 * t150;
t32 = t114 * t85 + t118 * t84 + t94 * t164 + t95 * t165;
t40 = -t113 * t181 - t117 * t61;
t123 = -t114 * t129 + t135;
t2 = -t128 * qJD(6) + t151;
t33 = -t183 * qJD(5) + t114 * t84 - t118 * t85;
t122 = (t114 * t163 + (-t117 * t118 + t178) * qJD(5)) * pkin(4);
t121 = (-t114 * t162 + (-t113 * t118 - t177) * qJD(5)) * pkin(4);
t100 = -0.2e1 * t144;
t69 = t115 * t107 - t147;
t66 = t116 * t168 - t142;
t59 = t78 * pkin(5) + t106;
t49 = -t105 * t163 + t121;
t48 = -t105 * t162 + t122;
t47 = t181 * pkin(5) + t64;
t45 = t54 * pkin(5) + t156;
t39 = -t113 * t61 + t139;
t34 = t79 * t107 + t53 * t120;
t30 = -t182 * qJD(4) + t184;
t20 = t123 * pkin(5) + t55;
t19 = t53 * pkin(10) + t33;
t18 = -t54 * pkin(10) - t32;
t16 = -t51 * qJD(6) - t113 * t54 - t117 * t53;
t12 = t52 * t107 - t16 * t120;
t11 = t40 * qJD(6) - t113 * t31 + t117 * t123;
t10 = qJD(6) * t139 + t113 * t123 + t117 * t31 - t61 * t163;
t7 = -t126 * qJD(6) - t113 * t18 + t117 * t19;
t6 = t127 * qJD(6) - t113 * t19 - t117 * t18;
t1 = (-qJD(6) * t14 - t5) * t117 + t140;
t3 = [0, 0, 0, 0, 0.2e1 * t144, -0.2e1 * t133, 0, 0, 0, t116 * t158, t120 * t158, -0.2e1 * t109 * t145 + 0.2e1 * t110 * t144, 0.2e1 * t109 * t134 + t132 * t190, 0.2e1 * t116 * t149 + 0.2e1 * t171 * t170, -0.2e1 * t115 * t133 + 0.2e1 * t116 * t147, t100, 0.2e1 * t63 * t107 - 0.2e1 * t30 * t120 + 0.2e1 * (t109 * t167 + t115 * t144) * t102, -0.2e1 * t109 * t102 * t168 - 0.2e1 * t29 * t120 + 0.2e1 * (-t62 + t83) * t107, 0.2e1 * t61 * t31, 0.2e1 * t61 * t123 + 0.2e1 * t31 * t181, -0.2e1 * t61 * t107 + 0.2e1 * t31 * t120, -0.2e1 * t181 * t107 + 0.2e1 * t123 * t120, t100, -0.2e1 * t9 * t120 + 0.2e1 * t55 * t181 + 0.2e1 * t64 * t135 + 0.2e1 * (t38 * qJD(3) + (-t64 * t161 * t115 - t46 * qJD(3)) * t114) * t116, -0.2e1 * t186 * t107 - 0.2e1 * t8 * t120 - 0.2e1 * t64 * t31 - 0.2e1 * t55 * t61, -0.2e1 * t40 * t10, 0.2e1 * t10 * t39 - 0.2e1 * t40 * t11, 0.2e1 * t10 * t120 + 0.2e1 * t40 * t107, -0.2e1 * t39 * t107 + 0.2e1 * t120 * t11, t100, -0.2e1 * t2 * t120 + 0.2e1 * (-t113 * t15 + t117 * t14) * t107 + 0.2e1 * t20 * t39 + 0.2e1 * t47 * t11, -0.2e1 * t1 * t120 - 0.2e1 * t47 * t10 - 0.2e1 * t107 * t128 + 0.2e1 * t20 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t169, -t107, 0, -t87, t143, -t116 * t134 + t132, t145 * t190 - t172 * t169, t69, t67, 0 (pkin(8) * t173 + (-pkin(3) * t119 + t179) * t116) * qJD(4) + (t131 * t115 - t83) * qJD(3) (t102 * t175 + t130 * t115) * qJD(4) + (t131 * t119 + t120 * t179) * qJD(3), -t31 * t79 + t61 * t53, -t79 * t123 + t53 * t181 + t31 * t78 + t61 * t54, t34, -t124, 0, -t33 * t120 + t137 * t107 + t106 * t135 + t55 * t78 + t64 * t54 + (-t106 * t161 * t116 * t114 + qJD(4) * pkin(4) * t181) * t115, -t106 * t31 - t183 * t107 - t32 * t120 - t61 * t156 - t64 * t53 + t55 * t79, -t10 * t52 + t40 * t16, t10 * t51 - t52 * t11 - t16 * t39 - t40 * t17, t12, -t125, 0, -t127 * t107 + t59 * t11 - t7 * t120 + t47 * t17 + t20 * t51 + t45 * t39, -t59 * t10 - t107 * t126 - t6 * t120 + t47 * t16 + t20 * t52 + t45 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t169, 0, 0, 0, 0, 0, -t67, t69, 0, 0, 0, 0, 0, t124, t34, 0, 0, 0, 0, 0, t125, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t145, -0.2e1 * t134, 0, 0, 0, t115 * t160, t119 * t160, -0.2e1 * t79 * t53, 0.2e1 * t53 * t78 - 0.2e1 * t79 * t54, 0, 0, 0, 0.2e1 * t106 * t54 + 0.2e1 * t78 * t156, -0.2e1 * t106 * t53 + 0.2e1 * t79 * t156, 0.2e1 * t52 * t16, -0.2e1 * t16 * t51 - 0.2e1 * t52 * t17, 0, 0, 0, 0.2e1 * t59 * t17 + 0.2e1 * t45 * t51, 0.2e1 * t59 * t16 + 0.2e1 * t45 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t191, t107, t30, t29, 0, 0, -t31, -t123, t107, t107 * t188 + (-t42 + (t120 * pkin(4) - t41) * t114) * qJD(5) + t138 (-t114 * t107 + t120 * t164) * pkin(4) + t8, 0, 0, -t10, -t11, t107, -t49 * t120 + (-pkin(4) * t178 + t117 * t105) * t107 + t2, -t48 * t120 - (pkin(4) * t177 + t113 * t105) * t107 - t117 * t5 - t14 * t162 + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, t66, 0, 0, 0, 0, 0, -t123, t31, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, -t168, 0, -pkin(8) * t167, pkin(8) * t168, 0, 0, -t53, -t54, 0, t33, t32, 0, 0, t16, -t17, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t155, -0.2e1 * t154, 0, 0, 0, 0, 0, 0.2e1 * t49, 0.2e1 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t123, t107, t9, t8, 0, 0, -t10, -t11, t107, t117 * t157 + (t113 * t141 - t180) * qJD(6) + t151, t13 + (-t4 - t157) * t113 + (qJD(6) * t141 - t5) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, t31, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t54, 0, t33, t32, 0, 0, t16, -t17, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, -t154, 0, 0, 0, 0, 0, t113 * t136 + t121, t117 * t136 + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t153, -0.2e1 * t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, t107, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, 0, t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, -t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
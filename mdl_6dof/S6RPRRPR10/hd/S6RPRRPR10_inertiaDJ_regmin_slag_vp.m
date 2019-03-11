% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:37:10
% EndTime: 2019-03-09 05:37:15
% DurationCPUTime: 1.82s
% Computational Cost: add. (1213->226), mult. (2748->390), div. (0->0), fcn. (2168->6), ass. (0->140)
t90 = sin(qJ(4));
t155 = qJD(4) * t90;
t89 = sin(qJ(6));
t169 = t89 * t90;
t92 = cos(qJ(6));
t93 = cos(qJ(4));
t51 = t92 * t93 + t169;
t83 = qJD(4) * t93;
t22 = t51 * qJD(6) - t89 * t155 - t92 * t83;
t91 = sin(qJ(3));
t183 = t22 * t91;
t159 = qJ(5) * t93;
t175 = pkin(4) + pkin(5);
t103 = -t175 * t90 + t159;
t96 = -pkin(1) - pkin(7);
t182 = t96 + t103;
t152 = qJD(6) * t89;
t181 = t93 * t152 + t92 * t155 - t89 * t83;
t85 = t90 ^ 2;
t87 = t93 ^ 2;
t163 = t85 - t87;
t123 = t163 * qJD(4);
t158 = t90 * qJ(5);
t180 = -t175 * t93 - t158;
t112 = pkin(4) * t90 - t159;
t149 = t90 * qJD(5);
t33 = t112 * qJD(4) - t149;
t94 = cos(qJ(3));
t165 = t94 * t33;
t113 = t93 * pkin(4) + t158;
t54 = -pkin(3) - t113;
t170 = t54 * t91;
t171 = t94 * pkin(8);
t104 = t112 - t96;
t28 = t104 * t94;
t179 = (t170 + t171) * qJD(3) - qJD(4) * t28 - t165;
t147 = t93 * qJD(5);
t178 = t113 * qJD(4) - t147;
t172 = pkin(9) * t94;
t116 = t91 * pkin(3) - t171;
t53 = qJ(2) + t116;
t168 = t91 * t96;
t69 = t90 * t168;
t20 = t69 + (-t53 - t172) * t93 - t175 * t91;
t49 = t90 * t53;
t70 = t93 * t168;
t25 = t91 * qJ(5) + t49 + t70;
t21 = t90 * t172 + t25;
t111 = t89 * t20 + t92 * t21;
t82 = t94 * qJD(3);
t126 = t96 * t82;
t153 = qJD(4) * t96;
t137 = t93 * t153;
t173 = pkin(8) * t91;
t117 = pkin(3) * t94 + t173;
t50 = t117 * qJD(3) + qJD(2);
t121 = t90 * t126 + t91 * t137 + t53 * t155 - t93 * t50;
t154 = qJD(4) * t94;
t140 = t90 * t154;
t5 = pkin(9) * t140 + (pkin(9) * t91 * t93 - t175 * t94) * qJD(3) + t121;
t148 = t91 * qJD(3);
t133 = t90 * t148;
t138 = t93 * t154;
t41 = -t133 + t138;
t139 = t90 * t153;
t16 = -t93 * t126 + t91 * t139 - t90 * t50 - t53 * t83;
t145 = qJ(5) * qJD(3);
t79 = t94 * t145;
t80 = t91 * qJD(5);
t8 = -t16 + t79 + t80;
t6 = t41 * pkin(9) + t8;
t2 = -t111 * qJD(6) + t92 * t5 - t89 * t6;
t177 = 0.2e1 * qJD(2);
t176 = 0.2e1 * qJD(5);
t174 = pkin(8) - pkin(9);
t167 = t92 * t90;
t166 = t93 * t94;
t164 = t94 * t96;
t162 = t85 + t87;
t86 = t91 ^ 2;
t88 = t94 ^ 2;
t161 = t86 - t88;
t160 = t86 + t88;
t157 = qJD(3) * t28;
t156 = qJD(3) * t93;
t151 = qJD(6) * t92;
t150 = qJD(6) * t94;
t146 = qJ(2) * qJD(3);
t144 = -0.2e1 * pkin(3) * qJD(4);
t143 = pkin(4) * t82;
t142 = pkin(8) * t155;
t141 = pkin(8) * t83;
t67 = t174 * t93;
t136 = t92 * t150;
t132 = t90 * t83;
t131 = t92 * t82;
t129 = t93 * t148;
t128 = t93 * t82;
t127 = t91 * t82;
t125 = t94 * t162;
t124 = -t93 * t53 + t69;
t122 = t161 * qJD(3);
t71 = 0.2e1 * t127;
t120 = t90 * t129;
t119 = t174 * t155;
t118 = qJD(4) * t67;
t26 = -t91 * pkin(4) + t124;
t110 = t25 * t93 + t26 * t90;
t109 = t25 * t90 - t26 * t93;
t66 = t174 * t90;
t108 = t89 * t66 + t92 * t67;
t107 = t89 * t93 - t167;
t106 = t92 * qJ(5) - t175 * t89;
t1 = -t20 * t151 + t21 * t152 - t89 * t5 - t92 * t6;
t32 = t51 * t94;
t102 = qJD(3) * t51;
t18 = -t104 * t148 + t178 * t94;
t99 = -t18 + (t54 * t94 - t173) * qJD(4);
t15 = t121 - t143;
t97 = -t109 * qJD(4) + t15 * t90 + t8 * t93;
t47 = pkin(3) - t180;
t40 = t90 * t82 + t91 * t83;
t39 = t160 * t83;
t38 = -t129 - t140;
t37 = t91 * t155 - t128;
t36 = t160 * t155;
t31 = t89 * t166 - t94 * t167;
t30 = t89 * qJD(5) + t106 * qJD(6);
t29 = qJ(5) * t152 - t92 * qJD(5) + t151 * t175;
t27 = t103 * qJD(4) + t149;
t24 = t182 * t94;
t23 = t90 * t151 - t181;
t14 = t108 * qJD(6) - t92 * t118 - t89 * t119;
t13 = -t89 * t118 + t92 * t119 - t66 * t151 + t67 * t152;
t12 = t94 * t102 + (qJD(4) - qJD(6)) * t91 * t107;
t11 = t91 * t102 - t90 * t136 + t181 * t94;
t10 = qJD(4) * t32 + t89 * t129 - t92 * t133 - t93 * t136 - t150 * t169;
t9 = t89 * t128 - t90 * t131 + t183;
t7 = (t180 * qJD(4) + t147) * t94 - t182 * t148;
t3 = [0, 0, 0, 0, t177, qJ(2) * t177, -0.2e1 * t127, 0.2e1 * t122, 0, 0, 0, 0.2e1 * qJD(2) * t91 + 0.2e1 * t94 * t146, 0.2e1 * qJD(2) * t94 - 0.2e1 * t91 * t146, -0.2e1 * t87 * t127 - 0.2e1 * t88 * t132, 0.4e1 * t94 * t120 + 0.2e1 * t88 * t123, -0.2e1 * t91 * t140 - 0.2e1 * t161 * t156, 0.2e1 * t90 * t122 - 0.2e1 * t91 * t138, t71, -0.2e1 * t88 * t137 - 0.2e1 * t121 * t91 + 0.2e1 * (-t124 + 0.2e1 * t69) * t82, 0.2e1 * t88 * t139 + 0.2e1 * t16 * t91 + 0.2e1 * (-t49 + t70) * t82, 0.2e1 * (-t90 * t157 - t15) * t91 + 0.2e1 * (-qJD(3) * t26 + t18 * t90 + t28 * t83) * t94, 0.2e1 * t109 * t148 + 0.2e1 * (-t110 * qJD(4) + t15 * t93 - t8 * t90) * t94, 0.2e1 * (t28 * t156 + t8) * t91 + 0.2e1 * (qJD(3) * t25 + t28 * t155 - t18 * t93) * t94, 0.2e1 * t26 * t15 + 0.2e1 * t28 * t18 + 0.2e1 * t25 * t8, -0.2e1 * t32 * t11, 0.2e1 * t32 * t10 + 0.2e1 * t11 * t31, 0.2e1 * t11 * t91 - 0.2e1 * t32 * t82, -0.2e1 * t10 * t91 + 0.2e1 * t31 * t82, t71, -0.2e1 * t2 * t91 - 0.2e1 * (t92 * t20 - t89 * t21) * t82 + 0.2e1 * t7 * t31 - 0.2e1 * t24 * t10, -0.2e1 * t1 * t91 - 0.2e1 * t24 * t11 + 0.2e1 * t111 * t82 + 0.2e1 * t7 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t36, -t39, 0, -t36 (t110 * qJD(3) - t18) * t94 + (t97 + t157) * t91, 0, 0, 0, 0, 0, -t94 * t10 + (t9 + (t107 * t94 - t31) * qJD(3)) * t91, -t94 * t11 + t12 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t162) * t71, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t148, -t82, 0, -t96 * t148, -t126, -t94 * t123 - t120, -0.4e1 * t94 * t132 + t163 * t148, t40, -t37, 0 (-t117 * t93 - t90 * t164) * qJD(4) + (t116 * t90 - t70) * qJD(3) (t117 * t90 - t93 * t164) * qJD(4) + (-pkin(8) * t166 + (pkin(3) * t93 + t90 * t96) * t91) * qJD(3), -t179 * t90 + t99 * t93, t97, t179 * t93 + t99 * t90, pkin(8) * t97 + t18 * t54 + t28 * t33, t107 * t11 - t32 * t22, -t10 * t107 + t11 * t51 + t22 * t31 - t32 * t23, t107 * t82 + t183, t23 * t91 + t51 * t82, 0, t14 * t91 - (t92 * t66 - t89 * t67) * t82 + t27 * t31 - t47 * t10 + t7 * t51 + t24 * t23, -t107 * t7 + t108 * t82 - t47 * t11 - t13 * t91 - t24 * t22 + t27 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, -t82, 0, 0, 0, 0, 0, t38, -t41, t38, qJD(3) * t125, t41, -t165 + (pkin(8) * t125 + t170) * qJD(3), 0, 0, 0, 0, 0, -t51 * t148 + t94 * t23, t107 * t148 - t94 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t132, -0.2e1 * t123, 0, 0, 0, t90 * t144, t93 * t144, 0.2e1 * t54 * t155 - 0.2e1 * t33 * t93, 0, -0.2e1 * t33 * t90 - 0.2e1 * t54 * t83, 0.2e1 * t54 * t33, 0.2e1 * t107 * t22, 0.2e1 * t107 * t23 + 0.2e1 * t22 * t51, 0, 0, 0, 0.2e1 * t47 * t23 + 0.2e1 * t27 * t51, -0.2e1 * t107 * t27 - 0.2e1 * t47 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t41, t82, -t121, t16, -t121 + 0.2e1 * t143 (pkin(4) * t148 - qJ(5) * t154) * t93 + (t91 * t145 + (pkin(4) * qJD(4) - qJD(5)) * t94) * t90, -t16 + 0.2e1 * t79 + 0.2e1 * t80, -t15 * pkin(4) + t8 * qJ(5) + t25 * qJD(5), 0, 0, t11, -t10, t82, t30 * t91 - (-t89 * qJ(5) - t175 * t92) * t82 - t2, t106 * t82 - t29 * t91 - t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t37, -t40, 0, -t37, -t112 * t82 - t178 * t91, 0, 0, 0, 0, 0, t9, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t155, 0, -t141, t142, -t141, -t178, -t142, -t178 * pkin(8), 0, 0, t22, t23, 0, t14, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, qJ(5) * t176, 0, 0, 0, 0, 0, 0.2e1 * t30, -0.2e1 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t38, 0, t15, 0, 0, 0, 0, 0, t91 * t152 - t131, t151 * t91 + t82 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, t141, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, -t82, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, 0, -t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, -t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;

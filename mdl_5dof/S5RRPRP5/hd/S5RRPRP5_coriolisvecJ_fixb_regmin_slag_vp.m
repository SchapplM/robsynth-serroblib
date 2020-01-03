% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:56
% EndTime: 2019-12-31 19:55:00
% DurationCPUTime: 1.13s
% Computational Cost: add. (2457->203), mult. (6464->270), div. (0->0), fcn. (4676->6), ass. (0->133)
t122 = sin(qJ(4));
t120 = sin(pkin(8));
t121 = cos(pkin(8));
t123 = sin(qJ(2));
t124 = cos(qJ(2));
t101 = t120 * t124 + t121 * t123;
t136 = t101 * qJD(2);
t176 = cos(qJ(4));
t129 = t176 * t136;
t143 = t120 * t123 - t121 * t124;
t93 = t143 * qJD(1);
t94 = t101 * qJD(1);
t138 = t122 * t93 - t176 * t94;
t154 = qJD(1) * qJD(2);
t149 = t124 * t154;
t150 = t123 * t154;
t84 = -t120 * t150 + t121 * t149;
t128 = -qJD(1) * t129 + t138 * qJD(4) - t122 * t84;
t117 = qJD(2) + qJD(4);
t164 = t138 * t117;
t187 = t128 - t164;
t131 = qJD(1) * t136;
t151 = qJD(4) * t176;
t155 = qJD(4) * t122;
t137 = -t122 * t131 - t93 * t151 - t94 * t155 + t176 * t84;
t49 = -t122 * t94 - t176 * t93;
t163 = t49 * t117;
t7 = t137 - t163;
t171 = t49 ^ 2;
t185 = t138 ^ 2;
t186 = -t171 + t185;
t152 = -t124 * pkin(2) - pkin(1);
t144 = t152 * qJD(1);
t108 = qJD(3) + t144;
t58 = t93 * pkin(3) + t108;
t12 = -pkin(4) * t49 + qJ(5) * t138 + t58;
t184 = t12 * t49;
t183 = t58 * t49;
t174 = t12 * t138;
t169 = t138 * t49;
t182 = t138 * t58;
t21 = -pkin(4) * t138 - t49 * qJ(5);
t181 = -0.2e1 * t154;
t113 = t121 * pkin(2) + pkin(3);
t175 = pkin(2) * t120;
t153 = t122 * t175;
t178 = t93 * pkin(7);
t168 = -qJ(3) - pkin(6);
t109 = t168 * t123;
t105 = qJD(1) * t109;
t110 = t168 * t124;
t106 = qJD(1) * t110;
t162 = t121 * t106;
t56 = -t120 * t105 + t162;
t40 = t56 + t178;
t177 = t94 * pkin(7);
t96 = t120 * t106;
t57 = t121 * t105 + t96;
t41 = t57 - t177;
t180 = qJD(4) * t153 - t113 * t151 + t122 * t40 + t176 * t41;
t60 = t121 * t109 + t120 * t110;
t44 = -t101 * pkin(7) + t60;
t61 = t120 * t109 - t121 * t110;
t45 = -t143 * pkin(7) + t61;
t139 = -t122 * t45 + t176 * t44;
t135 = t143 * qJD(2);
t147 = qJD(2) * t168;
t90 = t124 * qJD(3) + t123 * t147;
t91 = -t123 * qJD(3) + t124 * t147;
t42 = -t120 * t90 + t121 * t91;
t32 = pkin(7) * t135 + t42;
t43 = t120 * t91 + t121 * t90;
t33 = -pkin(7) * t136 + t43;
t4 = t139 * qJD(4) + t122 * t32 + t176 * t33;
t173 = t4 * t117;
t18 = t122 * t44 + t176 * t45;
t5 = t18 * qJD(4) + t122 * t33 - t176 * t32;
t170 = t5 * t117;
t134 = t122 * t113 + t176 * t175;
t167 = -t134 * qJD(4) + t122 * t41 - t176 * t40;
t166 = -qJD(5) + t180;
t71 = t90 * qJD(1);
t72 = t91 * qJD(1);
t39 = t120 * t72 + t121 * t71;
t165 = qJD(2) * pkin(2);
t100 = t105 + t165;
t53 = t120 * t100 - t162;
t126 = qJD(1) ^ 2;
t161 = t124 * t126;
t125 = qJD(2) ^ 2;
t160 = t125 * t123;
t159 = t125 * t124;
t52 = t121 * t100 + t96;
t36 = qJD(2) * pkin(3) - t177 + t52;
t37 = t53 - t178;
t10 = -t122 * t37 + t176 * t36;
t158 = qJD(5) - t10;
t157 = t123 ^ 2 - t124 ^ 2;
t156 = qJD(1) * t123;
t115 = t123 * t165;
t65 = pkin(2) * t156 + t94 * pkin(3);
t38 = -t120 * t71 + t121 * t72;
t28 = -t84 * pkin(7) + t38;
t29 = -pkin(7) * t131 + t39;
t146 = -t122 * t28 - t36 * t151 + t37 * t155 - t176 * t29;
t2 = t122 * t29 + t37 * t151 + t36 * t155 - t176 * t28;
t145 = pkin(1) * t181;
t116 = t117 * qJD(5);
t1 = t116 - t146;
t141 = t10 * t117 + t146;
t11 = t122 * t36 + t176 * t37;
t140 = t11 * t117 - t2;
t133 = t167 * t117 - t2;
t132 = t176 * t143;
t130 = t137 + t163;
t66 = pkin(3) * t136 + t115;
t74 = t143 * pkin(3) + t152;
t112 = pkin(2) * t150;
t59 = pkin(3) * t131 + t112;
t55 = t176 * t101 - t122 * t143;
t127 = -t128 - t164;
t3 = -pkin(4) * t128 - qJ(5) * t137 + qJD(5) * t138 + t59;
t88 = -t176 * t113 - pkin(4) + t153;
t87 = qJ(5) + t134;
t54 = t122 * t101 + t132;
t23 = t55 * qJD(4) - t122 * t135 + t129;
t22 = t101 * t155 + t117 * t132 + t122 * t136;
t16 = t54 * pkin(4) - t55 * qJ(5) + t74;
t15 = t21 + t65;
t9 = t117 * qJ(5) + t11;
t8 = -t117 * pkin(4) + t158;
t6 = t23 * pkin(4) + t22 * qJ(5) - t55 * qJD(5) + t66;
t13 = [0, 0, 0, 0.2e1 * t123 * t149, t157 * t181, t159, -t160, 0, -pkin(6) * t159 + t123 * t145, pkin(6) * t160 + t124 * t145, -t43 * t93 - t39 * t143 - t42 * t94 - t60 * t84 - t38 * t101 + (-t53 * t101 + t52 * t143 - t61 * t94) * qJD(2), t38 * t60 + t39 * t61 + t52 * t42 + t53 * t43 + (t108 + t144) * t115, t137 * t55 + t138 * t22, t128 * t55 - t137 * t54 + t138 * t23 - t22 * t49, -t22 * t117, -t23 * t117, 0, -t128 * t74 + t58 * t23 - t49 * t66 + t59 * t54 - t170, t137 * t74 - t138 * t66 - t58 * t22 + t59 * t55 - t173, t12 * t23 - t128 * t16 + t3 * t54 - t49 * t6 - t170, -t1 * t54 + t128 * t18 - t137 * t139 - t138 * t5 + t2 * t55 - t8 * t22 - t9 * t23 + t4 * t49, t12 * t22 - t137 * t16 + t138 * t6 - t3 * t55 + t173, t1 * t18 + t12 * t6 - t139 * t2 + t3 * t16 + t9 * t4 + t8 * t5; 0, 0, 0, -t123 * t161, t157 * t126, 0, 0, 0, t126 * pkin(1) * t123, pkin(1) * t161, (t53 + t56) * t94 - (-t57 + t52) * t93 + (-t120 * t131 - t121 * t84) * pkin(2), -t52 * t56 - t53 * t57 + (-t108 * t156 + t39 * t120 + t121 * t38) * pkin(2), t169, t186, t7, t187, 0, t49 * t65 + t133 + t182, t180 * t117 + t138 * t65 + t146 - t183, t15 * t49 + t133 + t174, t137 * t88 + t87 * t128 + (-t166 - t8) * t49 + (t167 - t9) * t138, -t166 * t117 - t138 * t15 + t1 + t184, t1 * t87 - t12 * t15 - t166 * t9 - t167 * t8 + t2 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93 ^ 2 - t94 ^ 2, t52 * t94 + t53 * t93 + t112, 0, 0, 0, 0, 0, t127, t130, t127, -t171 - t185, -t130, t138 * t8 - t9 * t49 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t186, t7, t187, 0, t140 + t182, t141 - t183, t21 * t49 + t140 + t174, -pkin(4) * t137 + t128 * qJ(5) - (-t11 + t9) * t138 - (t8 - t158) * t49, -t138 * t21 + 0.2e1 * t116 - t141 + t184, -t2 * pkin(4) + t1 * qJ(5) - t8 * t11 - t12 * t21 + t158 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t7, -t117 ^ 2 - t185, -t9 * t117 - t174 + t2;];
tauc_reg = t13;

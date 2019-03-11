% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:37
% EndTime: 2019-03-09 01:42:41
% DurationCPUTime: 1.33s
% Computational Cost: add. (1679->220), mult. (4094->284), div. (0->0), fcn. (3000->8), ass. (0->125)
t101 = sin(qJ(4));
t164 = cos(qJ(4));
t96 = sin(pkin(10));
t98 = cos(pkin(10));
t79 = t101 * t98 + t164 * t96;
t175 = t79 * qJD(1);
t176 = qJD(6) + t175;
t100 = sin(qJ(6));
t179 = t100 * t176;
t102 = cos(qJ(6));
t140 = qJD(1) * t101;
t134 = t96 * t140;
t136 = t164 * t98;
t126 = qJD(1) * t136;
t86 = qJD(4) * t126;
t62 = qJD(4) * t134 - t86;
t55 = t102 * t62;
t114 = -t176 * t179 - t55;
t137 = t100 * qJD(4);
t69 = -t126 + t134;
t51 = -t102 * t69 + t137;
t180 = t176 * t51;
t178 = qJD(6) - t176;
t177 = t175 * qJD(4);
t153 = pkin(7) * qJD(1);
t90 = sin(pkin(9)) * pkin(1) + qJ(3);
t81 = t90 * qJD(1);
t93 = t98 * qJD(2);
t49 = t93 + (-t81 - t153) * t96;
t60 = t96 * qJD(2) + t98 * t81;
t50 = t98 * t153 + t60;
t156 = -t101 * t50 + t164 * t49;
t173 = qJD(5) - t156;
t166 = t69 * pkin(5);
t22 = t101 * t49 + t164 * t50;
t19 = -qJD(4) * qJ(5) - t22;
t10 = -t19 - t166;
t168 = pkin(4) + pkin(8);
t172 = t168 * t62 + (t10 - t22 + t166) * t176;
t138 = qJD(6) * t102;
t78 = t101 * t96 - t136;
t135 = t78 * t138;
t158 = t78 * t62;
t74 = t79 * qJD(4);
t171 = -t100 * (-t176 * t74 + t158) + t176 * t135;
t132 = qJD(4) * t164;
t133 = qJD(3) * t164;
t165 = pkin(7) + t90;
t75 = t165 * t96;
t76 = t165 * t98;
t23 = (qJD(3) * t96 + qJD(4) * t76) * t101 + t75 * t132 - t98 * t133;
t170 = t69 ^ 2;
t169 = t175 ^ 2;
t63 = qJD(1) * t74;
t167 = t63 * pkin(4);
t80 = -cos(pkin(9)) * pkin(1) - t98 * pkin(3) - pkin(2);
t107 = -t79 * qJ(5) + t80;
t25 = t168 * t78 + t107;
t163 = t25 * t62;
t68 = t80 * qJD(1) + qJD(3);
t105 = -qJ(5) * t175 + t68;
t29 = t69 * pkin(4) + t105;
t162 = t29 * t175;
t53 = t102 * qJD(4) + t100 * t69;
t161 = t53 * t69;
t160 = t69 * t51;
t159 = t175 * t69;
t31 = -qJD(6) * t137 + t100 * t63 + t69 * t138;
t139 = qJD(4) * t101;
t73 = -t98 * t132 + t96 * t139;
t157 = t31 * t79 - t53 * t73;
t154 = t96 ^ 2 + t98 ^ 2;
t151 = t102 * t176;
t150 = t31 * t102;
t149 = t62 * qJ(5);
t148 = t69 * qJ(5);
t147 = qJD(6) * t78;
t146 = t23 * qJD(4);
t113 = -t101 * t75 + t164 * t76;
t24 = t79 * qJD(3) + t113 * qJD(4);
t145 = t24 * qJD(4);
t144 = t73 * qJD(4);
t65 = t74 * qJD(4);
t142 = pkin(5) * t175 + t173;
t131 = qJD(3) * t140;
t130 = qJD(1) * t154;
t125 = qJD(1) * t133;
t129 = t98 * t125 - t96 * t131 + t49 * t132 - t50 * t139;
t9 = t96 * t125 + t98 * t131 + t50 * t132 + t49 * t139;
t35 = t101 * t76 + t164 * t75;
t56 = t102 * t63;
t32 = t53 * qJD(6) - t56;
t124 = -t79 * t32 + t73 * t51;
t123 = (-t96 * t81 + t93) * t96 - t60 * t98;
t121 = t175 * t74 - t158;
t120 = -t79 * t63 + t73 * t69;
t15 = t168 * t69 + t105;
t6 = -t168 * qJD(4) + t142;
t2 = t100 * t6 + t102 * t15;
t119 = t100 * t15 - t102 * t6;
t117 = -qJD(5) * t175 + t149;
t116 = t73 * qJ(5) - t79 * qJD(5);
t112 = t22 * qJD(4) - t9;
t111 = -t147 * t179 + t74 * t151 - t78 * t55;
t7 = -qJD(4) * qJD(5) - t129;
t3 = -t63 * pkin(5) - t7;
t110 = t3 + (t168 * t176 + t148) * t176;
t27 = t79 * pkin(5) + t35;
t109 = t10 * t74 + t27 * t62 + t3 * t78;
t108 = t100 * t62 - t151 * t176;
t66 = qJD(4) * t69;
t40 = t62 * t79;
t38 = pkin(4) * t175 + t148;
t34 = t78 * pkin(4) + t107;
t30 = t74 * pkin(4) + t116;
t28 = -t78 * pkin(5) + t113;
t20 = t117 + t167;
t18 = -qJD(4) * pkin(4) + t173;
t16 = t168 * t74 + t116;
t14 = -t73 * pkin(5) + t24;
t13 = -t74 * pkin(5) - t23;
t8 = t168 * t63 + t117;
t5 = -t62 * pkin(5) + t9;
t4 = t102 * t5;
t1 = [0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t130 (t90 * t130 - t123) * qJD(3), -t175 * t73 - t40, t120 - t121, -t144, -t65, 0, t80 * t63 + t68 * t74 - t145, -t80 * t62 - t68 * t73 + t146, -t113 * t63 + t175 * t24 - t18 * t73 + t19 * t74 + t23 * t69 - t35 * t62 + t7 * t78 + t9 * t79, -t20 * t78 - t29 * t74 - t30 * t69 - t34 * t63 + t145, -t175 * t30 - t20 * t79 + t29 * t73 + t34 * t62 - t146, -t113 * t7 + t18 * t24 + t19 * t23 + t20 * t34 + t29 * t30 + t9 * t35, t53 * t135 + (t31 * t78 + t53 * t74) * t100 (-t100 * t51 + t102 * t53) * t74 + (-t100 * t32 + t150 + (-t100 * t53 - t102 * t51) * qJD(6)) * t78, t157 + t171, t111 + t124, -t176 * t73 - t40, t119 * t73 + t13 * t51 + t28 * t32 + t4 * t79 + (-t16 * t176 - t8 * t79 + t163) * t100 + (t14 * t176 - t109) * t102 + ((-t100 * t27 - t102 * t25) * t176 - t2 * t79 + t10 * t100 * t78) * qJD(6), t13 * t53 + t2 * t73 + t28 * t31 + (-(qJD(6) * t27 + t16) * t176 + t163 - (qJD(6) * t6 + t8) * t79 + t10 * t147) * t102 + (-(-qJD(6) * t25 + t14) * t176 - (-qJD(6) * t15 + t5) * t79 + t109) * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t144, t120 + t121, t65, -t144, t18 * t74 + t19 * t73 - t7 * t79 + t9 * t78, 0, 0, 0, 0, 0, t111 - t124, t157 - t171; 0, 0, 0, 0, 0, 0, -t154 * qJD(1) ^ 2, t123 * qJD(1), 0, 0, 0, 0, 0, 0.2e1 * t177, t86 + (-t69 - t134) * qJD(4), -t169 - t170, -0.2e1 * t177, t62 + t66, t167 + t149 - t19 * t69 + (-qJD(5) - t18) * t175, 0, 0, 0, 0, 0, t108 + t160, -t114 + t161; 0, 0, 0, 0, 0, 0, 0, 0, t159, t169 - t170, t86 + (t69 - t134) * qJD(4), 0, 0, -t175 * t68 + t112, qJD(4) * t156 + t68 * t69 - t129, pkin(4) * t62 - qJ(5) * t63 + (-t19 - t22) * t175 + (t18 - t173) * t69, t38 * t69 - t112 + t162, -t29 * t69 + t38 * t175 + (0.2e1 * qJD(5) - t156) * qJD(4) + t129, -t9 * pkin(4) - t7 * qJ(5) - t173 * t19 - t18 * t22 - t29 * t38, -t179 * t53 + t150 (-t176 * t53 - t32) * t102 + (-t31 + t180) * t100, t114 + t161, t108 - t160, t176 * t69, qJ(5) * t32 + t110 * t100 + t102 * t172 - t119 * t69 + t142 * t51, qJ(5) * t31 - t100 * t172 + t110 * t102 + t142 * t53 - t2 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62 + t66, -t159, -qJD(4) ^ 2 - t169, t19 * qJD(4) + t162 + t9, 0, 0, 0, 0, 0, -qJD(4) * t51 + t114, -qJD(4) * t53 + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t51, -t51 ^ 2 + t53 ^ 2, t31 + t180, -t178 * t53 + t56, -t62, -t10 * t53 - t100 * t8 - t178 * t2 + t4, t10 * t51 - t100 * t5 - t102 * t8 + t119 * t178;];
tauc_reg  = t1;

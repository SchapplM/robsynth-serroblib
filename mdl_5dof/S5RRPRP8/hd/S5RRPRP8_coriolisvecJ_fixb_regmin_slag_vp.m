% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:19
% EndTime: 2021-01-15 20:52:25
% DurationCPUTime: 1.27s
% Computational Cost: add. (1306->201), mult. (3101->263), div. (0->0), fcn. (1861->4), ass. (0->130)
t100 = qJD(2) - qJD(4);
t106 = sin(qJ(4));
t107 = sin(qJ(2));
t108 = cos(qJ(4));
t109 = cos(qJ(2));
t124 = t107 * t106 + t109 * t108;
t47 = t124 * qJD(1);
t156 = t47 * t100;
t118 = t124 * qJD(4);
t143 = qJD(1) * qJD(2);
t138 = t109 * t143;
t139 = t107 * t143;
t18 = qJD(1) * t118 - t106 * t139 - t108 * t138;
t174 = t18 + t156;
t144 = qJD(4) * t108;
t145 = qJD(4) * t106;
t146 = qJD(2) * t109;
t176 = t106 * t146 + t107 * t144 - t109 * t145;
t175 = -0.2e1 * t143;
t148 = qJD(1) * t109;
t149 = qJD(1) * t107;
t49 = -t106 * t148 + t108 * t149;
t154 = t49 * t100;
t19 = t176 * qJD(1) - t108 * t139;
t173 = t19 + t154;
t170 = t49 ^ 2;
t46 = t47 ^ 2;
t172 = -t46 + t170;
t110 = -pkin(2) - pkin(3);
t169 = pkin(6) - pkin(7);
t141 = t110 * qJD(2);
t92 = pkin(6) * t149;
t66 = pkin(7) * t149 - t92;
t37 = qJD(3) + t141 - t66;
t102 = qJD(2) * qJ(3);
t93 = pkin(6) * t148;
t68 = -pkin(7) * t148 + t93;
t51 = t102 + t68;
t136 = -t106 * t51 + t108 * t37;
t155 = t49 * qJ(5);
t7 = t136 - t155;
t6 = -t100 * pkin(4) + t7;
t168 = t6 - t7;
t167 = t49 * t47;
t127 = -t106 * t37 - t108 * t51;
t157 = t47 * qJ(5);
t8 = -t127 - t157;
t166 = t8 * t100;
t135 = -t106 * t66 + t108 * t68;
t10 = t135 - t157;
t71 = t108 * qJ(3) + t106 * t110;
t39 = -t106 * qJD(3) - t71 * qJD(4);
t165 = -t10 + t39;
t163 = t106 * t68 + t108 * t66;
t11 = t155 + t163;
t125 = -t106 * qJ(3) + t108 * t110;
t38 = t108 * qJD(3) + t125 * qJD(4);
t164 = -t11 + t38;
t96 = t107 * qJD(3);
t160 = qJ(3) * t138 + qJD(1) * t96;
t159 = qJ(3) * t146 + t96;
t158 = qJD(2) * pkin(2);
t112 = qJD(1) ^ 2;
t153 = t109 * t112;
t111 = qJD(2) ^ 2;
t152 = t111 * t107;
t151 = t111 * t109;
t103 = t107 ^ 2;
t150 = -t109 ^ 2 + t103;
t147 = qJD(2) * t107;
t73 = -t109 * pkin(2) - t107 * qJ(3) - pkin(1);
t52 = -qJD(1) * pkin(1) - pkin(2) * t148 - qJ(3) * t149;
t76 = t169 * t109;
t142 = t107 * t153;
t60 = t109 * pkin(3) - t73;
t137 = -t47 * pkin(4) - qJD(5);
t101 = qJD(2) * qJD(3);
t67 = t169 * t147;
t40 = -qJD(1) * t67 + t101;
t85 = pkin(6) * t138;
t59 = -pkin(7) * t138 + t85;
t134 = -t106 * t59 - t108 * t40 - t37 * t144 + t51 * t145;
t133 = -t106 * t40 + t108 * t59 - t51 * t144 - t37 * t145;
t132 = qJD(1) * t73 + t52;
t131 = pkin(1) * t175;
t34 = pkin(3) * t148 - t52;
t130 = qJD(3) - t158;
t129 = t100 ^ 2;
t128 = t107 * t141;
t75 = t169 * t107;
t126 = -t106 * t75 - t108 * t76;
t1 = -t19 * qJ(5) - t47 * qJD(5) - t134;
t88 = qJ(3) * t148;
t43 = t110 * t149 + t88;
t33 = pkin(2) * t139 - t160;
t45 = pkin(2) * t147 - t159;
t123 = -pkin(6) * t111 - qJD(1) * t45 - t33;
t122 = -t34 * t49 + t133;
t121 = t34 * t47 + t134;
t120 = t18 * qJ(5) + t133;
t69 = qJD(2) * t76;
t119 = t106 * t69 - t108 * t67 + t75 * t144 - t76 * t145;
t29 = t128 + t159;
t14 = -t137 + t34;
t116 = t14 * t47 - t1;
t24 = qJD(1) * t128 + t160;
t5 = t19 * pkin(4) + t24;
t114 = t126 * qJD(4) + t106 * t67 + t108 * t69;
t70 = -pkin(6) * t139 + t101;
t72 = t130 + t92;
t74 = t93 + t102;
t113 = t70 * t109 + (t109 * t72 + (-t74 + t93) * t107) * qJD(2);
t65 = -pkin(4) + t125;
t64 = pkin(2) * t149 - t88;
t63 = -t109 * t106 + t107 * t108;
t28 = t39 * t100;
t27 = t38 * t100;
t23 = pkin(4) * t124 + t60;
t22 = t124 * qJD(2) - t118;
t21 = -t108 * t147 + t176;
t20 = -t49 * pkin(4) + t43;
t17 = -t108 * t129 - t49 * t149;
t16 = -t106 * t129 - t47 * t149;
t13 = -qJ(5) * t124 - t126;
t12 = -t63 * qJ(5) - t106 * t76 + t108 * t75;
t9 = t21 * pkin(4) + t29;
t4 = -t22 * qJ(5) - t63 * qJD(5) + t114;
t3 = -t21 * qJ(5) - qJD(5) * t124 + t119;
t2 = -t49 * qJD(5) + t120;
t15 = [0, 0, 0, 0.2e1 * t107 * t138, t150 * t175, t151, -t152, 0, -pkin(6) * t151 + t107 * t131, pkin(6) * t152 + t109 * t131, t123 * t109 + t132 * t147, t113, t123 * t107 - t132 * t146, pkin(6) * t113 + t33 * t73 + t52 * t45, -t18 * t63 + t49 * t22, t124 * t18 - t63 * t19 - t49 * t21 - t22 * t47, -t22 * t100, t21 * t100, 0, -t114 * t100 + t124 * t24 + t60 * t19 + t34 * t21 + t29 * t47, t119 * t100 - t60 * t18 + t34 * t22 + t24 * t63 + t29 * t49, -t4 * t100 + t124 * t5 + t14 * t21 + t23 * t19 + t9 * t47, t3 * t100 + t14 * t22 - t23 * t18 + t9 * t49 + t5 * t63, -t1 * t124 + t12 * t18 - t13 * t19 - t2 * t63 - t8 * t21 - t6 * t22 - t3 * t47 - t4 * t49, t1 * t13 + t2 * t12 + t14 * t9 + t5 * t23 + t8 * t3 + t6 * t4; 0, 0, 0, -t142, t150 * t112, 0, 0, 0, t112 * pkin(1) * t107, pkin(1) * t153, (-t107 * t52 + t109 * t64) * qJD(1), ((t74 - t102) * t107 + (t130 - t72) * t109) * qJD(1), 0.2e1 * t101 + (t107 * t64 + t109 * t52) * qJD(1), t70 * qJ(3) + t74 * qJD(3) - t52 * t64 + (t107 * t74 + (-t72 - t158) * t109) * qJD(1) * pkin(6), -t167, -t172, t174, t173, 0, t135 * t100 - t43 * t47 - t122 - t28, -t163 * t100 - t43 * t49 - t121 + t27, t10 * t100 - t20 * t47 - t28 + (qJD(5) + t14) * t49 - t120, -t11 * t100 - t20 * t49 - t116 + t27, t65 * t18 - t71 * t19 + (-t8 - t165) * t49 + (t6 - t164) * t47, t1 * t71 - t14 * t20 + t164 * t8 + t165 * t6 + t2 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, 0, -t103 * t112 - t111, -t74 * qJD(2) + t52 * t149 + t85, 0, 0, 0, 0, 0, t16, t17, t16, t17, -t173 * t106 + t174 * t108, -t14 * t149 + (t2 - t166) * t108 + (t100 * t6 + t1) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167, t172, -t174, -t173, 0, t127 * t100 + t122, -t136 * t100 + t121, -t166 + (t137 - t14) * t49 + t120, -t170 * pkin(4) - t7 * t100 + t116, t18 * pkin(4) - t168 * t47, t168 * t8 + (-t14 * t49 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 - t154, -t18 + t156, -t46 - t170, t8 * t47 + t6 * t49 + t5;];
tauc_reg = t15;

% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:43
% EndTime: 2019-12-05 17:47:49
% DurationCPUTime: 1.32s
% Computational Cost: add. (2560->212), mult. (5501->290), div. (0->0), fcn. (3738->6), ass. (0->134)
t127 = sin(qJ(5));
t129 = cos(qJ(5));
t157 = qJD(5) * t129;
t158 = qJD(5) * t127;
t126 = cos(pkin(8));
t130 = cos(qJ(3));
t153 = qJD(1) * qJD(3);
t148 = t130 * t153;
t107 = t126 * t148;
t125 = sin(pkin(8));
t128 = sin(qJ(3));
t162 = qJD(1) * t128;
t149 = t125 * t162;
t79 = qJD(3) * t149 - t107;
t100 = t125 * t130 + t126 * t128;
t96 = t100 * qJD(3);
t80 = qJD(1) * t96;
t91 = t100 * qJD(1);
t161 = qJD(1) * t130;
t94 = t126 * t161 - t149;
t11 = -t127 * t79 + t129 * t80 + t91 * t157 + t94 * t158;
t120 = qJD(3) + qJD(5);
t38 = t127 * t94 + t129 * t91;
t170 = t38 * t120;
t192 = -t11 + t170;
t140 = -t127 * t91 + t129 * t94;
t175 = t140 ^ 2;
t176 = t38 ^ 2;
t191 = t175 - t176;
t174 = t38 * t140;
t12 = t140 * qJD(5) - t127 * t80 - t129 * t79;
t171 = t140 * t120;
t190 = -t12 + t171;
t131 = -pkin(1) - pkin(6);
t108 = t131 * qJD(1) + qJD(2);
t154 = t130 * qJD(4);
t160 = qJD(3) * t128;
t53 = -t108 * t160 + (qJ(4) * t160 - t154) * qJD(1);
t155 = t128 * qJD(4);
t159 = qJD(3) * t130;
t54 = t108 * t159 + (-qJ(4) * t159 - t155) * qJD(1);
t25 = -t125 * t54 + t126 * t53;
t14 = t80 * pkin(7) + t25;
t26 = t125 * t53 + t126 * t54;
t15 = t79 * pkin(7) + t26;
t180 = t94 * pkin(7);
t84 = -qJ(4) * t162 + t128 * t108;
t65 = t125 * t84;
t85 = -qJ(4) * t161 + t130 * t108;
t69 = qJD(3) * pkin(3) + t85;
t29 = t126 * t69 - t65;
t21 = qJD(3) * pkin(4) - t180 + t29;
t181 = t91 * pkin(7);
t172 = t126 * t84;
t30 = t125 * t69 + t172;
t22 = t30 - t181;
t1 = (qJD(5) * t21 + t15) * t129 + t127 * t14 - t22 * t158;
t104 = pkin(3) * t162 + qJD(1) * qJ(2) + qJD(4);
t52 = t91 * pkin(4) + t104;
t189 = t52 * t38 - t1;
t101 = -t125 * t128 + t126 * t130;
t138 = -t127 * t100 + t129 * t101;
t188 = t11 * t138;
t93 = t125 * t160 - t126 * t159;
t134 = t138 * qJD(5) - t127 * t96 - t129 * t93;
t187 = t134 * t120;
t6 = t127 * t21 + t129 * t22;
t2 = -t6 * qJD(5) - t127 * t15 + t129 * t14;
t186 = -t140 * t52 + t2;
t139 = t129 * t100 + t127 * t101;
t185 = t12 * t139 + t134 * t38;
t184 = t1 * t139 + t134 * t6 + t138 * t2;
t183 = t94 ^ 2;
t121 = qJD(1) * qJD(2);
t182 = 0.2e1 * t121;
t35 = -t125 * t85 - t172;
t27 = t35 + t181;
t36 = t126 * t85 - t65;
t28 = t36 - t180;
t113 = t126 * pkin(3) + pkin(4);
t177 = pkin(3) * t125;
t88 = t127 * t113 + t129 * t177;
t179 = t88 * qJD(5) - t127 * t28 + t129 * t27;
t87 = t129 * t113 - t127 * t177;
t178 = t87 * qJD(5) - t127 * t27 - t129 * t28;
t173 = t94 * t91;
t165 = qJ(4) - t131;
t81 = t165 * t160 - t154;
t106 = t165 * t130;
t82 = -qJD(3) * t106 - t155;
t32 = t125 * t81 + t126 * t82;
t105 = t165 * t128;
t51 = -t126 * t105 - t125 * t106;
t132 = qJD(3) ^ 2;
t169 = t132 * t128;
t168 = t132 * t130;
t133 = qJD(1) ^ 2;
t167 = t133 * qJ(2);
t166 = t93 * qJD(3);
t114 = t128 * pkin(3) + qJ(2);
t103 = pkin(3) * t148 + t121;
t164 = t128 ^ 2 - t130 ^ 2;
t163 = -t132 - t133;
t156 = t104 * qJD(1);
t109 = pkin(3) * t159 + qJD(2);
t152 = 0.2e1 * qJD(1);
t151 = pkin(3) * t161;
t150 = t130 * t133 * t128;
t31 = -t125 * t82 + t126 * t81;
t145 = t127 * t93 - t129 * t96;
t50 = t125 * t105 - t126 * t106;
t46 = -t79 * pkin(4) + t103;
t143 = t128 * t148;
t142 = -t100 * t79 - t91 * t93;
t141 = -t101 * t80 - t94 * t96;
t33 = -t101 * pkin(7) + t50;
t34 = -t100 * pkin(7) + t51;
t9 = -t127 * t34 + t129 * t33;
t10 = t127 * t33 + t129 * t34;
t136 = t26 * t100 + t25 * t101 - t29 * t96 - t30 * t93;
t118 = qJ(2) * t182;
t89 = t91 ^ 2;
t83 = t96 * qJD(3);
t70 = t100 * pkin(4) + t114;
t58 = t94 * pkin(4) + t151;
t55 = -t93 * pkin(4) + t109;
t24 = t93 * pkin(7) + t32;
t23 = t96 * pkin(7) + t31;
t18 = -t139 * qJD(5) + t145;
t17 = t100 * t157 + t101 * t158 - t145;
t5 = -t127 * t22 + t129 * t21;
t4 = -t10 * qJD(5) - t127 * t24 + t129 * t23;
t3 = t9 * qJD(5) + t127 * t23 + t129 * t24;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, t118, -0.2e1 * t143, 0.2e1 * t164 * t153, -t169, 0.2e1 * t143, -t168, 0, -t131 * t169 + (qJ(2) * t159 + qJD(2) * t128) * t152, -t131 * t168 + (-qJ(2) * t160 + qJD(2) * t130) * t152, 0, t118, t141, t80 * t100 + t101 * t79 + t96 * t91 + t94 * t93, -t83, t142, t166, 0, t31 * qJD(3) + t103 * t100 - t104 * t93 + t109 * t91 - t114 * t79, -t32 * qJD(3) + t103 * t101 - t104 * t96 + t109 * t94 - t114 * t80, -t31 * t94 - t32 * t91 + t50 * t80 + t51 * t79 - t136, t103 * t114 + t104 * t109 + t25 * t50 + t26 * t51 + t29 * t31 + t30 * t32, -t140 * t17 - t188, t11 * t139 - t12 * t138 - t134 * t140 + t17 * t38, -t17 * t120, t185, -t187, 0, t70 * t12 + t4 * t120 + t134 * t52 + t139 * t46 + t55 * t38, -t70 * t11 - t3 * t120 + t138 * t46 + t140 * t55 - t52 * t17, -t10 * t12 + t9 * t11 - t140 * t4 + t5 * t17 - t3 * t38 - t184, t1 * t10 + t2 * t9 + t6 * t3 + t5 * t4 + t46 * t70 + t52 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t167, 0, 0, 0, 0, 0, 0, t163 * t128, t163 * t130, 0, -t167, 0, 0, 0, 0, 0, 0, -qJD(1) * t91 - t83, -qJD(1) * t94 + t166, -t141 - t142, t136 - t156, 0, 0, 0, 0, 0, 0, -qJD(1) * t38 + t18 * t120, -qJD(1) * t140 - t187, -t140 * t18 - t185 + t188, -t52 * qJD(1) + t5 * t18 + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, -t164 * t133, 0, -t150, 0, 0, -t130 * t167, t128 * t167, 0, 0, t173, -t89 + t183, 0, -t173, -t107 + (t94 + t149) * qJD(3), 0, -t35 * qJD(3) - t104 * t94 - t91 * t151 + t25, t36 * qJD(3) + t104 * t91 - t94 * t151 - t26, (t30 + t35) * t94 + (-t29 + t36) * t91 + (t125 * t79 + t126 * t80) * pkin(3), -t29 * t35 - t30 * t36 + (t125 * t26 + t126 * t25 - t130 * t156) * pkin(3), t174, t191, t192, -t174, t190, 0, -t179 * t120 - t58 * t38 + t186, -t178 * t120 - t140 * t58 + t189, t87 * t11 - t88 * t12 + (t179 + t6) * t140 + (-t178 - t5) * t38, t1 * t88 + t178 * t6 - t179 * t5 + t2 * t87 - t52 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107 + (t94 - t149) * qJD(3), -0.2e1 * t91 * qJD(3), -t89 - t183, t29 * t94 + t30 * t91 + t103, 0, 0, 0, 0, 0, 0, t12 + t171, -t11 - t170, -t175 - t176, t140 * t5 + t38 * t6 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t191, t192, -t174, t190, 0, t6 * t120 + t186, t5 * t120 + t189, 0, 0;];
tauc_reg = t7;

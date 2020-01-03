% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:37
% EndTime: 2019-12-31 20:30:43
% DurationCPUTime: 1.60s
% Computational Cost: add. (1651->238), mult. (3876->332), div. (0->0), fcn. (2566->6), ass. (0->145)
t101 = sin(qJ(4));
t102 = sin(qJ(2));
t104 = cos(qJ(4));
t105 = cos(qJ(2));
t56 = t102 * t101 + t105 * t104;
t189 = t56 * qJD(1);
t197 = qJD(5) + t189;
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t154 = qJD(1) * t105;
t155 = qJD(1) * t102;
t51 = -t101 * t154 + t104 * t155;
t94 = qJD(2) - qJD(4);
t28 = t100 * t51 + t103 * t94;
t202 = t197 * t28;
t121 = t100 * t94 - t103 * t51;
t201 = t121 * t197;
t150 = qJD(4) * t104;
t151 = qJD(4) * t101;
t152 = qJD(2) * t105;
t200 = t101 * t152 + t102 * t150 - t105 * t151;
t199 = qJD(5) - t197;
t86 = pkin(6) * t155;
t198 = -pkin(7) * t155 + qJD(3) + t86;
t147 = qJD(1) * qJD(2);
t139 = t102 * t147;
t21 = t200 * qJD(1) - t104 * t139;
t162 = t100 * t21;
t190 = t103 * t197;
t196 = t121 * t51 + t190 * t197 + t162;
t148 = qJD(5) * t103;
t149 = qJD(5) * t100;
t116 = t56 * qJD(4);
t138 = t105 * t147;
t169 = t101 * t139 + t104 * t138;
t20 = -qJD(1) * t116 + t169;
t8 = t103 * t20 - t94 * t148 - t51 * t149;
t172 = t8 * t100;
t195 = t121 * t190 - t172;
t9 = -t121 * qJD(5) + t100 * t20;
t194 = (-t9 + t201) * t100 - (-t8 + t202) * t103;
t193 = -0.2e1 * t147;
t106 = -pkin(2) - pkin(3);
t141 = t106 * qJD(2);
t41 = t141 + t198;
t87 = pkin(6) * t154;
t63 = -pkin(7) * t154 + t87;
t96 = qJD(2) * qJ(3);
t52 = t63 + t96;
t17 = -t101 * t52 + t104 * t41;
t14 = t94 * pkin(4) - t17;
t192 = t197 * t14;
t191 = t51 * t94 + t21;
t22 = t51 * pkin(4) + pkin(8) * t189;
t84 = qJ(3) * t154;
t46 = t106 * t155 + t84;
t120 = t104 * qJ(3) + t101 * t106;
t60 = -pkin(8) + t120;
t153 = qJD(2) * t102;
t182 = pkin(6) - pkin(7);
t62 = t182 * t153;
t95 = qJD(2) * qJD(3);
t44 = -qJD(1) * t62 + t95;
t80 = pkin(6) * t138;
t54 = -pkin(7) * t138 + t80;
t7 = t101 * t44 - t104 * t54 + t52 * t150 + t41 * t151;
t188 = (qJD(5) * t60 - t22 + t46) * t197 - t7;
t187 = t7 + (pkin(8) * qJD(5) + t22) * t197;
t160 = t103 * t21;
t184 = t100 * t197 ^ 2 - t28 * t51 - t160;
t69 = t182 * t102;
t70 = t182 * t105;
t31 = t101 * t70 - t104 * t69;
t64 = qJD(2) * t70;
t10 = -t31 * qJD(4) + t101 * t64 - t104 * t62;
t53 = -qJD(1) * pkin(1) - pkin(2) * t154 - qJ(3) * t155;
t39 = pkin(3) * t154 - t53;
t12 = pkin(4) * t189 - t51 * pkin(8) + t39;
t67 = -t105 * pkin(2) - t102 * qJ(3) - pkin(1);
t55 = t105 * pkin(3) - t67;
t57 = -t105 * t101 + t102 * t104;
t16 = t56 * pkin(4) - t57 * pkin(8) + t55;
t26 = t56 * qJD(2) - t116;
t32 = t101 * t69 + t104 * t70;
t6 = t101 * t54 + t104 * t44 + t41 * t150 - t52 * t151;
t183 = -(qJD(5) * t16 + t10) * t197 - t32 * t21 - (qJD(5) * t12 + t6) * t56 + t7 * t57 + t14 * t26;
t18 = t101 * t41 + t104 * t52;
t15 = -t94 * pkin(8) + t18;
t122 = t100 * t15 - t103 * t12;
t181 = t122 * t51;
t2 = t100 * t12 + t103 * t15;
t180 = t2 * t51;
t179 = t14 * t57;
t178 = t16 * t21;
t177 = t197 * t51;
t176 = t189 * t94;
t175 = t51 * t189;
t173 = t57 * t21;
t119 = -t101 * qJ(3) + t104 * t106;
t171 = -t119 * qJD(4) + t101 * t63 - t198 * t104;
t170 = t120 * qJD(4) + t198 * t101 + t104 * t63;
t90 = t102 * qJD(3);
t167 = qJ(3) * t138 + qJD(1) * t90;
t166 = qJ(3) * t152 + t90;
t97 = t102 ^ 2;
t165 = -t105 ^ 2 + t97;
t164 = qJD(2) * pkin(2);
t108 = qJD(1) ^ 2;
t159 = t105 * t108;
t107 = qJD(2) ^ 2;
t158 = t107 * t102;
t157 = t107 * t105;
t146 = t189 ^ 2 - t51 ^ 2;
t145 = t197 * t155;
t144 = t57 * t149;
t143 = t197 * t148;
t142 = t102 * t159;
t135 = qJD(1) * t67 + t53;
t130 = pkin(1) * t193;
t129 = qJD(3) - t164;
t128 = t94 ^ 2;
t127 = t94 * t197;
t126 = t102 * t141;
t125 = t197 * t26 + t173;
t38 = pkin(2) * t139 - t167;
t47 = pkin(2) * t153 - t166;
t118 = -pkin(6) * t107 - qJD(1) * t47 - t38;
t35 = t126 + t166;
t114 = t39 * t51 + t7;
t27 = qJD(1) * t126 + t167;
t112 = -pkin(8) * t21 + t17 * t197 + t192;
t111 = -t189 * t39 + t6;
t110 = t171 * t197 - t60 * t21 - t192;
t65 = -pkin(6) * t139 + t95;
t66 = t129 + t86;
t68 = t87 + t96;
t109 = t65 * t105 + (t105 * t66 + (-t68 + t87) * t102) * qJD(2);
t59 = pkin(4) - t119;
t58 = pkin(2) * t155 - t84;
t25 = -t104 * t153 + t200;
t11 = t32 * qJD(4) - t101 * t62 - t104 * t64;
t5 = t25 * pkin(4) - t26 * pkin(8) + t35;
t4 = t21 * pkin(4) - t20 * pkin(8) + t27;
t3 = t103 * t4;
t1 = [0, 0, 0, 0.2e1 * t102 * t138, t165 * t193, t157, -t158, 0, -pkin(6) * t157 + t102 * t130, pkin(6) * t158 + t105 * t130, t118 * t105 + t135 * t153, t109, t118 * t102 - t135 * t152, t109 * pkin(6) + t38 * t67 + t53 * t47, t20 * t57 + t51 * t26, -t189 * t26 - t20 * t56 - t51 * t25 - t173, -t26 * t94, t25 * t94, 0, t11 * t94 + t189 * t35 + t55 * t21 + t39 * t25 + t27 * t56, t10 * t94 + t55 * t20 + t39 * t26 + t27 * t57 + t35 * t51, t121 * t144 + (-t121 * t26 + t57 * t8) * t103, (t100 * t121 - t103 * t28) * t26 + (-t172 - t103 * t9 + (t100 * t28 + t103 * t121) * qJD(5)) * t57, t103 * t125 - t121 * t25 - t144 * t197 + t8 * t56, -t100 * t125 - t143 * t57 - t28 * t25 - t9 * t56, t197 * t25 + t21 * t56, -t122 * t25 + t11 * t28 + t3 * t56 + t31 * t9 + (t178 + t5 * t197 + (-t15 * t56 - t197 * t32 + t179) * qJD(5)) * t103 + t183 * t100, -t11 * t121 - t2 * t25 + t31 * t8 + (-(-qJD(5) * t32 + t5) * t197 - t178 - (-qJD(5) * t15 + t4) * t56 - qJD(5) * t179) * t100 + t183 * t103; 0, 0, 0, -t142, t165 * t108, 0, 0, 0, t108 * pkin(1) * t102, pkin(1) * t159, (-t102 * t53 + t105 * t58) * qJD(1), ((t68 - t96) * t102 + (t129 - t66) * t105) * qJD(1), 0.2e1 * t95 + (t102 * t58 + t105 * t53) * qJD(1), t65 * qJ(3) + t68 * qJD(3) - t53 * t58 + (t102 * t68 + (-t66 - t164) * t105) * qJD(1) * pkin(6), -t175, t146, qJD(4) * t189 - t169 + t176, t191, 0, t170 * t94 - t189 * t46 + t114, -t171 * t94 - t46 * t51 + t111, t195, -t194, -t196, t184, t177, t110 * t100 - t188 * t103 + t170 * t28 + t59 * t9 - t181, t188 * t100 + t110 * t103 - t121 * t170 + t59 * t8 - t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, 0, -t97 * t108 - t107, -t68 * qJD(2) + t53 * t155 + t80, 0, 0, 0, 0, 0, -t101 * t128 - t155 * t189, -t104 * t128 - t51 * t155, 0, 0, 0, 0, 0, -t103 * t145 + (t100 * t127 - t9) * t104 + (-t28 * t94 - t143 - t162) * t101, t100 * t145 + (t103 * t127 - t8) * t104 + (t121 * t94 + t149 * t197 - t160) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, -t146, t20 - t176, -t191, 0, -t18 * t94 - t114, -t17 * t94 - t111, -t195, t194, t196, -t184, -t177, -pkin(4) * t9 + t112 * t100 - t187 * t103 - t18 * t28 + t181, -pkin(4) * t8 + t187 * t100 + t112 * t103 + t121 * t18 + t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121 * t28, t121 ^ 2 - t28 ^ 2, t8 + t202, -t9 - t201, t21, -t100 * t6 + t14 * t121 - t199 * t2 + t3, -t100 * t4 - t103 * t6 + t199 * t122 + t14 * t28;];
tauc_reg = t1;

% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:36:30
% EndTime: 2019-12-05 17:36:34
% DurationCPUTime: 1.52s
% Computational Cost: add. (1502->257), mult. (3023->342), div. (0->0), fcn. (1920->10), ass. (0->159)
t89 = sin(pkin(8));
t91 = cos(pkin(8));
t107 = pkin(3) * t91 + pkin(6) * t89 + pkin(2);
t92 = cos(pkin(7));
t180 = pkin(1) * t92;
t53 = -t107 - t180;
t96 = cos(qJ(4));
t166 = t91 * t96;
t90 = sin(pkin(7));
t77 = pkin(1) * t90 + qJ(3);
t57 = t77 * t166;
t94 = sin(qJ(4));
t26 = t94 * t53 + t57;
t185 = qJD(4) * t26;
t86 = qJ(1) + pkin(7);
t83 = cos(t86);
t176 = g(2) * t83;
t82 = sin(t86);
t115 = g(3) * t82 + t176;
t138 = qJD(1) * qJD(4);
t125 = t96 * t138;
t144 = qJDD(1) * t94;
t104 = t125 + t144;
t184 = t104 * t89;
t159 = pkin(1) * qJDD(1);
t56 = qJD(1) * qJD(3) + qJDD(1) * t77;
t167 = t91 * t94;
t41 = t82 * t167 + t83 * t96;
t43 = t83 * t167 - t82 * t96;
t183 = -g(2) * t41 + g(3) * t43;
t178 = g(1) * t94;
t182 = t89 * t178 + t183;
t157 = qJD(1) * t91;
t73 = -qJD(4) + t157;
t158 = qJD(1) * t89;
t128 = qJ(5) * t158;
t38 = t53 * qJD(1) + qJD(3);
t62 = t77 * qJD(1);
t47 = qJD(2) * t89 + t62 * t91;
t14 = t96 * t38 - t47 * t94;
t8 = -t96 * t128 + t14;
t5 = -pkin(4) * t73 + t8;
t181 = -t8 + t5;
t179 = pkin(4) * t94;
t141 = t91 * qJDD(1);
t72 = -qJDD(4) + t141;
t173 = t72 * pkin(4);
t97 = cos(qJ(1));
t172 = t97 * pkin(1);
t81 = t91 * qJD(2);
t46 = t62 * t89 - t81;
t171 = t46 * t89;
t170 = t72 * t91;
t84 = t89 ^ 2;
t98 = qJD(1) ^ 2;
t169 = t84 * t98;
t168 = t89 * t56;
t149 = qJD(4) * t96;
t153 = qJD(3) * t96;
t165 = t53 * t149 + t91 * t153;
t164 = t115 * t89;
t85 = t91 ^ 2;
t163 = t84 + t85;
t87 = t94 ^ 2;
t88 = t96 ^ 2;
t162 = -t87 - t88;
t161 = t87 - t88;
t160 = qJ(5) * t89;
t156 = qJD(1) * t94;
t155 = qJD(1) * t96;
t154 = qJD(3) * t94;
t15 = t38 * t94 + t47 * t96;
t152 = qJD(4) * t15;
t151 = qJD(4) * t47;
t150 = qJD(4) * t94;
t148 = qJD(5) * t89;
t28 = qJD(5) - t81 + (pkin(4) * t156 + t62) * t89;
t147 = qJD(5) + t28;
t78 = -pkin(2) - t180;
t145 = qJDD(1) * t78;
t143 = qJDD(1) * t96;
t142 = t89 * qJDD(1);
t140 = qJ(5) * qJDD(1);
t137 = qJD(1) * qJD(5);
t136 = t77 * t167;
t135 = t94 * t169;
t134 = t96 * t160;
t133 = t89 * t156;
t132 = t77 * t150;
t131 = t73 * t150;
t130 = t46 * t158;
t129 = t91 * t154;
t95 = sin(qJ(1));
t127 = -t95 * pkin(1) + t83 * qJ(3);
t80 = t91 * qJDD(2);
t36 = -t80 + t168;
t126 = -t36 * t91 - g(1);
t35 = t53 * qJDD(1) + qJDD(3);
t37 = qJDD(2) * t89 + t56 * t91;
t3 = t38 * t149 - t47 * t150 + t94 * t35 + t96 * t37;
t50 = (t77 + t179) * t89;
t124 = qJD(1) * t50 + t28;
t59 = qJDD(3) + t145;
t123 = t59 + t145;
t122 = t72 - t141;
t121 = t72 + t141;
t120 = t89 * pkin(4) * t125 + t142 * t179 + qJDD(5) - t80;
t119 = qJD(4) * t133;
t118 = t94 * t125;
t117 = -g(2) * t43 - g(3) * t41;
t42 = t82 * t166 - t83 * t94;
t44 = -t83 * t166 - t82 * t94;
t116 = -g(2) * t44 + g(3) * t42;
t114 = g(2) * t82 - g(3) * t83;
t113 = g(2) * t97 + g(3) * t95;
t9 = -t94 * t128 + t15;
t112 = t5 * t96 + t9 * t94;
t111 = t5 * t94 - t9 * t96;
t110 = t14 * t94 - t15 * t96;
t109 = t36 * t89 + t37 * t91;
t108 = t47 * t91 + t171;
t106 = (pkin(4) * t96 + pkin(3)) * t91 - t89 * (-qJ(5) - pkin(6)) + pkin(2);
t19 = t120 + t168;
t55 = (pkin(4) * t149 + qJD(3)) * t89;
t105 = qJD(1) * t55 + qJDD(1) * t50 + t19;
t31 = t96 * t35;
t103 = t31 + qJ(5) * t119 + (-qJD(4) * t38 - t37) * t94;
t102 = g(2) * t172 - g(3) * t127;
t101 = -t73 ^ 2 - t169;
t100 = g(1) * t89 * t96 - g(2) * t42 - g(3) * t44 - t3;
t4 = -t94 * t37 - t152 + t31;
t70 = t96 * t142;
t61 = t96 * t135;
t58 = t91 * t119;
t54 = t89 * t73 * t155;
t52 = t161 * t169;
t51 = t162 * t142;
t49 = t96 * t53;
t40 = (qJDD(1) * t88 - 0.2e1 * t118) * t84;
t39 = (qJDD(1) * t87 + 0.2e1 * t118) * t84;
t27 = 0.2e1 * (t161 * t138 - t94 * t143) * t84;
t25 = t49 - t136;
t24 = -t54 - t184;
t23 = t70 + (-qJD(4) - t73) * t133;
t22 = -t94 * t160 + t26;
t21 = -t129 - t185;
t20 = -t91 * t132 + t165;
t18 = -t134 + t49 + (-t77 * t94 - pkin(4)) * t91;
t17 = t101 * t96 + t94 * t72;
t16 = t101 * t94 - t96 * t72;
t13 = (t122 * t94 + (t73 - t157) * t149) * t89;
t12 = (t121 * t94 + (t73 + t157) * t149) * t89;
t11 = t58 + (-t121 * t96 + t131) * t89;
t10 = t58 + (t122 * t96 - t131) * t89;
t7 = -t129 - t96 * t148 + (-t57 + (-t53 + t160) * t94) * qJD(4);
t6 = -t94 * t148 + (-t134 - t136) * qJD(4) + t165;
t2 = (-t104 * qJ(5) - t94 * t137) * t89 + t3;
t1 = -t173 + (-t151 + (-t137 - t140) * t89) * t96 + t103;
t29 = [0, 0, 0, 0, 0, qJDD(1), t113, -g(2) * t95 + g(3) * t97, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t92 * t159 + t115, -0.2e1 * t90 * t159 - t114, 0, (t113 + (t90 ^ 2 + t92 ^ 2) * t159) * pkin(1), t84 * qJDD(1), 0.2e1 * t89 * t141, 0, t85 * qJDD(1), 0, 0, (t115 - t123) * t91, t123 * t89 - t164, t56 * t163 + t109 + t114, t59 * t78 - g(2) * (-t83 * pkin(2) - t82 * qJ(3) - t172) - g(3) * (-pkin(2) * t82 + t127) + t109 * t77 + t108 * qJD(3), t40, t27, t11, t39, t12, t170, -t21 * t73 - t25 * t72 - t4 * t91 + (t46 * t149 + t36 * t94) * t89 + (t77 * t144 + (t77 * t149 + t154) * qJD(1)) * t84 + t116, t20 * t73 + t26 * t72 + t3 * t91 + (-t46 * t150 + t36 * t96) * t89 + (t77 * t143 + (-t132 + t153) * qJD(1)) * t84 + t117, ((-t152 - qJDD(1) * t25 - t4 + (-t21 - t185) * qJD(1)) * t96 + (qJD(4) * t14 - qJDD(1) * t26 - t3 + (qJD(4) * t25 - t20) * qJD(1)) * t94) * t89 + t164, t3 * t26 + t15 * t20 + t4 * t25 + t14 * t21 + (qJD(3) * t46 + t36 * t77) * t89 + t107 * t176 + (g(2) * qJ(3) + g(3) * t107) * t82 + t102, t40, t27, t11, t39, t12, t170, -t1 * t91 - t18 * t72 - t7 * t73 + (t105 * t94 + t124 * t149) * t89 + t116, t2 * t91 + t22 * t72 + t6 * t73 + (t105 * t96 - t124 * t150) * t89 + t117, ((-qJD(4) * t9 - qJDD(1) * t18 - t1 + (-qJD(4) * t22 - t7) * qJD(1)) * t96 + (qJD(4) * t5 - qJDD(1) * t22 - t2 + (qJD(4) * t18 - t6) * qJD(1)) * t94) * t89 + t164, t2 * t22 + t9 * t6 + t1 * t18 + t5 * t7 + t19 * t50 + t28 * t55 + (g(2) * t106 - g(3) * t179) * t83 + (-g(2) * (-qJ(3) - t179) + g(3) * t106) * t82 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t89 + t126, 0, 0, 0, 0, 0, 0, t13, t10, 0, (t3 * t96 - t4 * t94 + (-t14 * t96 - t15 * t94) * qJD(4)) * t89 + t126, 0, 0, 0, 0, 0, 0, t13, t10, 0, -t19 * t91 - g(1) + (-qJD(4) * t112 - t1 * t94 + t2 * t96) * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, t142, -t163 * t98, -t108 * qJD(1) - t115 + t59, 0, 0, 0, 0, 0, 0, t16, t17, t51, t3 * t94 + t4 * t96 - t110 * qJD(4) + (t110 * t91 - t171) * qJD(1) - t115, 0, 0, 0, 0, 0, 0, t16, t17, t51, t1 * t96 + t2 * t94 - t111 * qJD(4) + (t111 * t91 - t28 * t89) * qJD(1) - t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t52, t23, -t61, t24, -t72, -t96 * t130 - t15 * t73 + t182 + t4, t94 * t130 - t14 * t73 + t100, 0, 0, t61, -t52, t23, -t61, t24, -t72, -0.2e1 * t173 - t73 * t9 + (-pkin(4) * t135 - t151 + (-t147 * qJD(1) - t140) * t89) * t96 + t103 + t182, -pkin(4) * t88 * t169 - t73 * t8 + (t94 * t140 + (qJ(5) * t149 + t147 * t94) * qJD(1)) * t89 + t100, (-pkin(4) * t143 + (pkin(4) * qJD(4) - t181) * t156) * t89, t181 * t9 + (t1 + (-t28 * t155 + t178) * t89 + t183) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54 + t184, t70 + (-qJD(4) + t73) * t133, t162 * t169, g(1) * t91 + (qJD(1) * t112 + t114 + t56) * t89 + t120;];
tau_reg = t29;

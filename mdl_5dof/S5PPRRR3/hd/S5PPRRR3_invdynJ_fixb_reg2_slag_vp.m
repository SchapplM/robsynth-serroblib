% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:17:06
% EndTime: 2019-12-05 15:17:10
% DurationCPUTime: 2.14s
% Computational Cost: add. (1901->296), mult. (4454->420), div. (0->0), fcn. (3519->12), ass. (0->158)
t107 = sin(qJ(5));
t110 = cos(qJ(5));
t111 = cos(qJ(4));
t153 = t111 * qJDD(3);
t108 = sin(qJ(4));
t154 = t108 * qJDD(3);
t68 = t107 * t111 + t108 * t110;
t99 = qJD(4) + qJD(5);
t199 = t99 * t68;
t21 = qJD(3) * t199 + t107 * t154 - t110 * t153;
t103 = sin(pkin(9));
t198 = qJD(2) * qJD(3) + qJDD(1) * t103;
t109 = sin(qJ(3));
t173 = t110 * t111;
t174 = t107 * t108;
t67 = -t173 + t174;
t53 = t67 * t109;
t162 = qJD(5) * t107;
t189 = qJD(4) * pkin(4);
t112 = cos(qJ(3));
t168 = qJD(1) * t103;
t66 = qJD(2) * t109 + t112 * t168;
t59 = qJD(3) * pkin(6) + t66;
t140 = pkin(7) * qJD(3) + t59;
t105 = cos(pkin(9));
t167 = qJD(1) * t105;
t146 = t111 * t167;
t31 = -t140 * t108 - t146;
t29 = t31 + t189;
t79 = t108 * t167;
t32 = t140 * t111 - t79;
t155 = t105 * qJDD(1);
t143 = qJD(3) * t168;
t152 = -t109 * qJDD(2) - t198 * t112;
t44 = -t109 * t143 - t152;
t38 = qJDD(3) * pkin(6) + t44;
t74 = qJD(4) * t79;
t6 = qJDD(4) * pkin(4) + t74 + (-pkin(7) * qJDD(3) - t38) * t108 + (-t140 * qJD(4) - t155) * t111;
t132 = -t108 * t155 + t111 * t38;
t7 = pkin(7) * t153 + t31 * qJD(4) + t132;
t1 = (qJD(5) * t29 + t7) * t110 + t107 * t6 - t32 * t162;
t197 = pkin(7) + pkin(6);
t179 = t103 * t112;
t60 = -t105 * t111 - t108 * t179;
t196 = g(3) * t60;
t195 = t105 ^ 2 * qJDD(1) - g(3);
t194 = g(3) * t103;
t164 = qJD(3) * t111;
t147 = t110 * t164;
t166 = qJD(3) * t108;
t148 = t107 * t166;
t62 = -t147 + t148;
t64 = t68 * qJD(3);
t193 = t64 * t62;
t71 = t197 * t108;
t72 = t197 * t111;
t47 = -t107 * t72 - t110 * t71;
t80 = t109 * t168;
t65 = t112 * qJD(2) - t80;
t150 = qJD(4) * t197;
t69 = t108 * t150;
t70 = t111 * t150;
t192 = t47 * qJD(5) - t107 * t70 - t110 * t69 + t67 * t65;
t48 = -t107 * t71 + t110 * t72;
t191 = -t48 * qJD(5) + t107 * t69 - t110 * t70 + t68 * t65;
t190 = qJD(3) * pkin(3);
t188 = t107 * t32;
t187 = t110 * t32;
t186 = t111 * t59;
t58 = -t65 - t190;
t185 = qJD(3) * t58;
t184 = qJDD(3) * pkin(3);
t104 = sin(pkin(8));
t183 = t103 * t104;
t106 = cos(pkin(8));
t182 = t103 * t106;
t181 = t103 * t108;
t180 = t103 * t111;
t178 = t104 * t109;
t177 = t104 * t112;
t176 = t106 * t109;
t175 = t106 * t112;
t115 = qJD(3) ^ 2;
t172 = t112 * t115;
t100 = t108 ^ 2;
t101 = t111 ^ 2;
t171 = t100 - t101;
t170 = t100 + t101;
t114 = qJD(4) ^ 2;
t169 = t114 + t115;
t165 = qJD(3) * t109;
t163 = qJD(3) * t112;
t161 = qJD(5) * t110;
t159 = qJD(3) * qJD(4);
t157 = qJDD(3) * t112;
t156 = qJDD(4) * t108;
t151 = t108 * t115 * t111;
t91 = pkin(4) * t111 + pkin(3);
t149 = t103 * t165;
t145 = -g(1) * t104 + g(2) * t106;
t142 = t108 * t159;
t141 = t111 * t159;
t138 = (-qJDD(2) + t143) * t112 + t198 * t109;
t137 = -qJD(5) * t147 - t107 * t153 + (-t141 - t154) * t110;
t136 = t170 * qJDD(3);
t135 = t108 * t141;
t134 = t108 * t189 - t66;
t55 = t105 * t177 - t176;
t57 = t105 * t175 + t178;
t133 = g(1) * t57 + g(2) * t55;
t130 = t99 * t174;
t9 = t107 * t29 + t187;
t61 = -t105 * t108 + t111 * t179;
t24 = -t107 * t61 + t110 * t60;
t25 = t107 * t60 + t110 * t61;
t45 = -t108 * t59 - t146;
t46 = -t79 + t186;
t129 = t108 * t45 - t111 * t46;
t128 = -t109 * t115 + t157;
t127 = qJDD(3) * t109 + t172;
t37 = t138 - t184;
t126 = t142 - t153;
t54 = -t105 * t178 - t175;
t56 = -t105 * t176 + t177;
t125 = -g(1) * t56 - g(2) * t54 + t109 * t194;
t124 = g(3) * t179 + t133;
t123 = -pkin(6) * qJDD(4) + (t58 + t65 - t190) * qJD(4);
t2 = -t9 * qJD(5) - t107 * t7 + t110 * t6;
t122 = qJD(3) * t66 + t125;
t12 = t45 * qJD(4) + t132;
t13 = -t108 * t38 + t74 + (-qJD(4) * t59 - t155) * t111;
t121 = t12 * t111 + (-t108 * t46 - t111 * t45) * qJD(4) - t13 * t108;
t119 = -pkin(6) * t114 + t122 + t184 - t37;
t49 = -t91 * qJD(3) - t65;
t102 = qJ(4) + qJ(5);
t95 = sin(t102);
t96 = cos(t102);
t118 = -g(1) * (-t95 * t182 - t57 * t96) - g(2) * (-t95 * t183 - t55 * t96) - g(3) * (t105 * t95 - t96 * t179) + t49 * t62 - t1;
t117 = -g(1) * (t96 * t182 - t57 * t95) - g(2) * (t96 * t183 - t55 * t95) - g(3) * (-t105 * t96 - t95 * t179) - t49 * t64 + t2;
t116 = t121 - t124;
t98 = qJDD(4) + qJDD(5);
t52 = t68 * t109;
t42 = t60 * qJD(4) - t111 * t149;
t41 = -t61 * qJD(4) + t108 * t149;
t39 = -qJD(4) * t173 - t111 * t161 + t130;
t23 = t126 * pkin(4) + t37;
t22 = -t62 ^ 2 + t64 ^ 2;
t20 = qJD(3) * t130 + t137;
t17 = -t68 * t163 + t53 * t99;
t16 = -t109 * t199 - t67 * t163;
t15 = t64 * t99 - t21;
t14 = -t137 + (-t148 + t62) * t99;
t11 = t110 * t31 - t188;
t10 = -t107 * t31 - t187;
t8 = t110 * t29 - t188;
t4 = -t25 * qJD(5) - t107 * t42 + t110 * t41;
t3 = t24 * qJD(5) + t107 * t41 + t110 * t42;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) * t103 ^ 2 + t195, 0, 0, 0, 0, 0, 0, -t127 * t103, -t128 * t103, 0, (t109 * t138 + t112 * t44 + (-t109 * t66 - t112 * t65) * qJD(3)) * t103 + t195, 0, 0, 0, 0, 0, 0, qJD(4) * t41 + qJDD(4) * t60 + (t126 * t109 - t111 * t172) * t103, -qJD(4) * t42 - qJDD(4) * t61 + (t108 * t127 + t109 * t141) * t103, (-t108 * t60 + t111 * t61) * qJDD(3) + (-t108 * t41 + t111 * t42 + (-t108 * t61 - t111 * t60) * qJD(4)) * qJD(3), t12 * t61 + t13 * t60 + t41 * t45 + t42 * t46 - g(3) + (t109 * t37 + t163 * t58) * t103, 0, 0, 0, 0, 0, 0, t24 * t98 + t4 * t99 + (t109 * t21 + t163 * t62) * t103, -t25 * t98 - t3 * t99 + (-t109 * t20 + t163 * t64) * t103, t20 * t24 - t21 * t25 - t3 * t62 - t4 * t64, t1 * t25 + t2 * t24 + t3 * t9 + t4 * t8 - g(3) + (t109 * t23 + t163 * t49) * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t145, 0, 0, 0, 0, 0, 0, t128, -t127, 0, t109 * t44 - t112 * t138 + (-t109 * t65 + t112 * t66) * qJD(3) + t145, 0, 0, 0, 0, 0, 0, (-0.2e1 * t142 + t153) * t112 + (-t169 * t111 - t156) * t109, (-qJDD(4) * t109 - 0.2e1 * t112 * t159) * t111 + (t109 * t169 - t157) * t108, t109 * t136 + t170 * t172, (-qJD(3) * t129 - t37) * t112 + (t121 + t185) * t109 + t145, 0, 0, 0, 0, 0, 0, -t112 * t21 + t165 * t62 + t17 * t99 - t52 * t98, t112 * t20 - t16 * t99 + t165 * t64 + t53 * t98, -t16 * t62 - t17 * t64 - t20 * t52 + t21 * t53, -t1 * t53 - t112 * t23 + t16 * t9 + t165 * t49 + t17 * t8 - t2 * t52 + t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t122 - t138, (t65 + t80) * qJD(3) + t124 + t152, 0, 0, qJDD(3) * t100 + 0.2e1 * t135, 0.2e1 * t108 * t153 - 0.2e1 * t171 * t159, t111 * t114 + t156, qJDD(3) * t101 - 0.2e1 * t135, qJDD(4) * t111 - t108 * t114, 0, t123 * t108 + t119 * t111, -t108 * t119 + t111 * t123, -qJD(3) * t170 * t65 + pkin(6) * t136 + t116, -t58 * t66 + t129 * t65 + (t125 - t37) * pkin(3) + t116 * pkin(6), -t20 * t68 - t39 * t64, -t199 * t64 + t20 * t67 - t21 * t68 + t39 * t62, -t39 * t99 + t68 * t98, t199 * t62 + t21 * t67, -t199 * t99 - t67 * t98, 0, t125 * t96 + t134 * t62 + t191 * t99 + t199 * t49 - t21 * t91 + t23 * t67 + t47 * t98, -t125 * t95 + t134 * t64 - t192 * t99 + t20 * t91 + t23 * t68 - t39 * t49 - t48 * t98, -t1 * t67 - t191 * t64 - t192 * t62 - t199 * t9 - t2 * t68 + t20 * t47 - t21 * t48 + t39 * t8 - t124, t1 * t48 + t2 * t47 - t23 * t91 - g(1) * (t197 * t57 + t56 * t91) - g(2) * (t197 * t55 + t54 * t91) + t192 * t9 + t191 * t8 + t134 * t49 - (-t109 * t91 + t112 * t197) * t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, t171 * t115, t154, t151, t153, qJDD(4), -t196 + t74 + (-t155 + (-g(1) * t106 - g(2) * t104) * t103) * t111 + (t133 - t38 - t185) * t108 + (t46 - t186) * qJD(4), -t58 * t164 - g(1) * (-t106 * t181 - t111 * t57) - g(2) * (-t104 * t181 - t111 * t55) + g(3) * t61 - t132, 0, 0, t193, t22, t14, -t193, t15, t98, -t10 * t99 + (t110 * t98 - t162 * t99 - t166 * t62) * pkin(4) + t117, t11 * t99 + (-t107 * t98 - t161 * t99 - t166 * t64) * pkin(4) + t118, (t10 + t9) * t64 + (t11 - t8) * t62 + (-t107 * t21 + t110 * t20 + (t107 * t64 - t110 * t62) * qJD(5)) * pkin(4), -t8 * t10 - t9 * t11 + (t1 * t107 + t2 * t110 - t49 * t166 - g(1) * (t106 * t180 - t108 * t57) - g(2) * (t104 * t180 - t108 * t55) - t196 + (-t107 * t8 + t110 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t22, t14, -t193, t15, t98, t9 * t99 + t117, t8 * t99 + t118, 0, 0;];
tau_reg = t5;

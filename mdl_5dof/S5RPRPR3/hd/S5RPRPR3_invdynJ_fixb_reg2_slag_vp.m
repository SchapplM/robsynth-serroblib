% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:30
% EndTime: 2020-01-03 11:36:36
% DurationCPUTime: 1.42s
% Computational Cost: add. (2398->238), mult. (4079->325), div. (0->0), fcn. (2468->14), ass. (0->153)
t110 = sin(pkin(8));
t195 = pkin(1) * t110;
t153 = qJD(3) * t195;
t112 = cos(pkin(8));
t92 = pkin(1) * t112 + pkin(2);
t199 = qJD(1) * t153 - t92 * qJDD(1);
t105 = qJ(1) + pkin(8);
t97 = qJ(3) + t105;
t90 = sin(t97);
t91 = cos(t97);
t188 = g(2) * t91 + g(3) * t90;
t114 = sin(qJ(3));
t117 = cos(qJ(3));
t175 = pkin(1) * qJDD(1);
t152 = t110 * t175;
t70 = t92 * qJD(1);
t179 = t114 * t70;
t146 = qJD(3) * t179 + t114 * t152 + t199 * t117;
t134 = qJDD(4) + t146;
t101 = qJDD(1) + qJDD(3);
t193 = t101 * pkin(3);
t24 = t134 - t193;
t127 = -t24 - t188;
t109 = sin(pkin(9));
t181 = t109 * t90;
t111 = cos(pkin(9));
t194 = pkin(4) * t111;
t198 = pkin(7) * t181 + t90 * t194;
t180 = t109 * t91;
t197 = pkin(7) * t180 + t91 * t194;
t104 = qJD(1) + qJD(3);
t102 = t109 ^ 2;
t103 = t111 ^ 2;
t162 = t102 + t103;
t196 = t104 * t162;
t154 = qJD(1) * t195;
t41 = -t114 * t154 + t117 * t70;
t137 = qJD(4) - t41;
t113 = sin(qJ(5));
t116 = cos(qJ(5));
t159 = qJD(4) * t111;
t164 = t111 * t116;
t165 = t111 * t113;
t59 = -pkin(7) * t109 - pkin(3) - t194;
t43 = -qJ(4) * t165 + t116 * t59;
t173 = qJD(5) * t43;
t42 = t117 * t154 + t179;
t192 = -t113 * t42 + t116 * t159 - t164 * t41 + t173;
t44 = qJ(4) * t164 + t113 * t59;
t172 = qJD(5) * t44;
t191 = -t113 * t159 - t116 * t42 + t165 * t41 - t172;
t190 = g(2) * t180 + g(3) * t181;
t189 = t91 * pkin(3) + t90 * qJ(4);
t187 = -g(2) * t90 + g(3) * t91;
t51 = t114 * t92 + t117 * t195;
t115 = sin(qJ(1));
t95 = sin(t105);
t186 = t115 * pkin(1) + pkin(2) * t95;
t118 = cos(qJ(1));
t96 = cos(t105);
t185 = t118 * pkin(1) + pkin(2) * t96;
t25 = t104 * t59 + t137;
t33 = qJ(4) * t104 + t42;
t29 = qJD(2) * t109 + t111 * t33;
t8 = t113 * t25 + t116 * t29;
t184 = qJD(5) * t8;
t183 = t104 * t42;
t47 = t51 * qJD(3);
t182 = t104 * t47;
t160 = qJD(3) * t117;
t147 = -t199 * t114 + t117 * t152 + t70 * t160;
t18 = qJ(4) * t101 + qJD(4) * t104 + t147;
t14 = -t111 * qJDD(2) + t109 * t18;
t178 = t14 * t109;
t28 = -t111 * qJD(2) + t109 * t33;
t177 = t28 * t109;
t166 = t111 * t101;
t69 = -qJDD(5) + t166;
t176 = t69 * t111;
t50 = -t114 * t195 + t117 * t92;
t34 = -t50 + t59;
t48 = qJ(4) + t51;
t20 = t116 * t34 - t48 * t165;
t174 = qJD(5) * t20;
t171 = t101 * t113;
t170 = t101 * t116;
t100 = t104 ^ 2;
t169 = t102 * t100;
t168 = t104 * t111;
t167 = t109 * t101;
t163 = t113 * t116;
t106 = t113 ^ 2;
t107 = t116 ^ 2;
t161 = t106 - t107;
t158 = qJD(5) * t113;
t157 = qJD(5) * t116;
t150 = t109 * t158;
t7 = -t113 * t29 + t116 * t25;
t156 = t7 * t150 - t190;
t155 = t24 * t109 + t190;
t73 = -qJD(5) + t168;
t151 = t73 * t158;
t149 = -t14 * t111 - g(1);
t148 = t90 * pkin(3) - qJ(4) * t91;
t145 = t185 + t189;
t144 = t69 - t166;
t143 = t69 + t166;
t142 = t104 * (-qJD(5) - t73);
t141 = t162 * t101;
t139 = t163 * t169;
t138 = t104 * t113 * t157;
t136 = -g(2) * t118 - g(3) * t115;
t135 = t113 * t7 - t116 * t8;
t133 = -t183 - t193;
t49 = -pkin(3) - t50;
t132 = t101 * t49 + t182;
t131 = t148 + t186;
t15 = qJDD(2) * t109 + t111 * t18;
t130 = t15 * t111 + t178 + t187;
t129 = -t146 - t188;
t128 = -t147 - t187;
t46 = -t114 * t153 + t92 * t160;
t126 = g(1) * t109 - qJD(5) * t25 - t15;
t21 = t113 * t34 + t48 * t164;
t125 = -qJD(5) * t29 - t104 * t177;
t45 = qJD(4) + t46;
t124 = (t14 * t48 + t28 * t45) * t109;
t123 = (t14 * qJ(4) + t137 * t28) * t109;
t13 = t59 * t101 + t134;
t11 = t116 * t13;
t3 = -t113 * t15 + t11 - t184;
t38 = -t113 * t91 + t164 * t90;
t40 = t113 * t90 + t164 * t91;
t122 = -g(2) * t40 - g(3) * t38 - t3 * t111 + t113 * t178 + t157 * t177;
t121 = -t73 ^ 2 - t169;
t2 = qJD(5) * t7 + t113 * t13 + t116 * t15;
t37 = -t116 * t91 - t165 * t90;
t39 = -t90 * t116 + t165 * t91;
t120 = g(2) * t39 - g(3) * t37 + t2 * t111 + t116 * t178 - t28 * t150;
t108 = qJDD(2) - g(1);
t85 = t103 * t101;
t84 = t102 * t101;
t57 = 0.2e1 * t109 * t166;
t54 = t150 * t168;
t36 = (t101 * t107 - 0.2e1 * t138) * t102;
t35 = (t101 * t106 + 0.2e1 * t138) * t102;
t32 = -pkin(3) * t104 + t137;
t27 = 0.2e1 * (qJD(5) * t104 * t161 - t101 * t163) * t102;
t17 = (t143 * t113 + (t73 + t168) * t157) * t109;
t16 = t54 + (-t116 * t143 + t151) * t109;
t5 = -qJD(5) * t21 + t116 * t47 - t165 * t45;
t4 = t113 * t47 + t164 * t45 + t174;
t1 = [0, 0, 0, 0, 0, qJDD(1), t136, g(2) * t115 - g(3) * t118, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -g(2) * t96 - g(3) * t95 + 0.2e1 * t112 * t175, g(2) * t95 - g(3) * t96 - 0.2e1 * t152, 0, (t136 + (t110 ^ 2 + t112 ^ 2) * t175) * pkin(1), 0, 0, 0, 0, 0, t101, t101 * t50 + t129 - t182, -t101 * t51 - t104 * t46 + t128, 0, -g(2) * t185 - g(3) * t186 - t146 * t50 + t147 * t51 - t41 * t47 + t42 * t46, t84, t57, 0, t85, 0, 0, (t127 - t132) * t111, t109 * t132 + t155, t141 * t48 + t45 * t196 + t130, t24 * t49 + t32 * t47 - g(2) * t145 - g(3) * t131 + (t15 * t48 + t29 * t45) * t111 + t124, t36, t27, t16, t35, t17, t176, -t20 * t69 - t5 * t73 + (t48 * t171 + (t113 * t45 + t157 * t48) * t104) * t102 + t122, t21 * t69 + t4 * t73 + (t48 * t170 + (t116 * t45 - t158 * t48) * t104) * t102 + t120, ((-t101 * t21 - t2 + (-t4 + t174) * t104) * t113 + (-t101 * t20 - t104 * t5 - t3 + (-t104 * t21 - t8) * qJD(5)) * t116) * t109 + t156, t2 * t21 + t8 * t4 + t3 * t20 + t7 * t5 - g(2) * (t145 + t197) - g(3) * (t131 + t198) + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t15 + t149, 0, 0, 0, 0, 0, 0, (t144 * t113 + (t73 - t168) * t157) * t109, t54 + (t116 * t144 - t151) * t109, 0, (-t113 * t3 + t116 * t2 + (-t113 * t8 - t116 * t7) * qJD(5)) * t109 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t129 + t183, t104 * t41 + t128, 0, 0, t84, t57, 0, t85, 0, 0, (t127 - t133) * t111, t109 * t133 + t155, qJ(4) * t141 + t137 * t196 + t130, -t24 * pkin(3) - t32 * t42 - g(2) * t189 - g(3) * t148 + (t15 * qJ(4) + t137 * t29) * t111 + t123, t36, t27, t16, t35, t17, t176, -t43 * t69 - t191 * t73 + (qJ(4) * t171 + (qJ(4) * t157 + t113 * t137) * t104) * t102 + t122, t44 * t69 + t192 * t73 + (qJ(4) * t170 + (-qJ(4) * t158 + t116 * t137) * t104) * t102 + t120, ((-t101 * t44 - t2) * t113 + (-t101 * t43 - t184 - t3) * t116 + ((-t172 - t191) * t116 + (t173 - t192) * t113) * t104) * t109 + t156, t2 * t44 + t3 * t43 - g(2) * (t189 + t197) - g(3) * (t148 + t198) + t192 * t8 + t191 * t7 + t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, t167, -t162 * t100, (-t111 * t29 - t177) * t104 - t127, 0, 0, 0, 0, 0, 0, t113 * t121 - t116 * t69, t113 * t69 + t116 * t121, (-t106 - t107) * t167, t2 * t113 + t3 * t116 - t135 * qJD(5) + (t111 * t135 - t177) * t104 + t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t161 * t169, (t113 * t142 + t170) * t109, -t139, (t116 * t142 - t171) * t109, -t69, -g(2) * t37 - g(3) * t39 + t113 * t126 + t116 * t125 - t8 * t73 + t11, g(2) * t38 - g(3) * t40 - t7 * t73 + t126 * t116 + (-t125 - t13) * t113, 0, 0;];
tau_reg = t1;

% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:53
% EndTime: 2019-12-31 18:02:56
% DurationCPUTime: 1.97s
% Computational Cost: add. (2495->325), mult. (4471->418), div. (0->0), fcn. (2662->8), ass. (0->169)
t151 = qJ(2) * qJDD(1);
t98 = -pkin(1) - pkin(2);
t60 = t98 * qJDD(1) + qJDD(2);
t91 = sin(pkin(8));
t92 = cos(pkin(8));
t175 = t92 * t151 + t91 * t60;
t150 = qJD(1) * qJD(2);
t67 = t92 * t150;
t32 = t67 + t175;
t28 = -qJDD(1) * pkin(6) + t32;
t208 = -qJD(3) * qJD(4) - t28;
t149 = qJD(1) * qJD(4);
t97 = cos(qJ(4));
t137 = t97 * t149;
t95 = sin(qJ(4));
t152 = t95 * qJDD(1);
t207 = t137 + t152;
t138 = t95 * t149;
t80 = t97 * qJDD(1);
t206 = -t138 + t80;
t195 = sin(qJ(1));
t196 = cos(qJ(1));
t47 = -t195 * t92 + t196 * t91;
t202 = g(2) * t47;
t46 = -t195 * t91 - t196 * t92;
t204 = g(1) * t46;
t127 = t202 + t204;
t200 = g(3) * t97;
t155 = qJ(2) * qJD(1);
t61 = t98 * qJD(1) + qJD(2);
t39 = t92 * t155 + t91 * t61;
t34 = -qJD(1) * pkin(6) + t39;
t191 = t34 * t97;
t26 = qJD(3) * t95 + t191;
t79 = t97 * qJDD(3);
t10 = -qJD(4) * t26 - t28 * t95 + t79;
t6 = -qJDD(4) * pkin(4) - t10;
t107 = -t127 * t95 + t200 - t6;
t63 = qJD(1) * t97 + qJD(5);
t205 = -pkin(7) * qJD(5) * t63 + t107;
t192 = t34 * t95;
t25 = qJD(3) * t97 - t192;
t147 = -t95 * qJDD(3) + t208 * t97;
t162 = qJD(4) * t95;
t9 = -t34 * t162 - t147;
t102 = -(t25 * t97 + t26 * t95) * qJD(4) - t10 * t95 + t9 * t97;
t203 = g(2) * t46;
t201 = g(3) * t95;
t199 = t6 * t95;
t22 = qJD(4) * pkin(7) + t26;
t130 = pkin(4) * t97 + pkin(7) * t95;
t38 = -t91 * t155 + t61 * t92;
t33 = qJD(1) * pkin(3) - t38;
t23 = t130 * qJD(1) + t33;
t94 = sin(qJ(5));
t96 = cos(qJ(5));
t7 = -t22 * t94 + t23 * t96;
t198 = t7 * t63;
t8 = t22 * t96 + t23 * t94;
t197 = t8 * t63;
t161 = qJD(4) * t96;
t166 = qJD(1) * t95;
t48 = t94 * t166 + t161;
t159 = qJD(5) * t48;
t17 = qJDD(4) * t94 - t207 * t96 + t159;
t194 = t17 * t94;
t156 = t94 * qJD(4);
t157 = qJD(5) * t96;
t113 = t97 * t156 + t95 * t157;
t18 = qJD(1) * t113 - qJD(5) * t156 + qJDD(4) * t96 + t94 * t152;
t193 = t18 * t96;
t190 = t47 * t94;
t189 = t47 * t96;
t49 = t96 * t166 - t156;
t188 = t48 * t49;
t187 = t48 * t63;
t186 = t48 * t94;
t185 = t48 * t96;
t184 = t49 * t63;
t183 = t49 * t94;
t182 = t49 * t96;
t181 = t94 * t97;
t180 = t95 * t17;
t179 = t95 * t18;
t178 = t96 * t97;
t141 = t91 * t162;
t41 = -t91 * t181 - t92 * t96;
t177 = qJD(5) * t41 - t96 * t141 - (t92 * t178 + t91 * t94) * qJD(1);
t42 = t91 * t178 - t92 * t94;
t176 = -qJD(5) * t42 + t94 * t141 - (-t92 * t181 + t91 * t96) * qJD(1);
t55 = t92 * qJ(2) + t91 * t98;
t174 = t196 * pkin(1) + t195 * qJ(2);
t173 = g(1) * t195 - g(2) * t196;
t88 = t95 ^ 2;
t89 = t97 ^ 2;
t172 = t88 - t89;
t171 = t88 + t89;
t100 = qJD(1) ^ 2;
t170 = t100 * t92;
t99 = qJD(4) ^ 2;
t169 = t100 + t99;
t168 = pkin(1) * qJDD(1);
t167 = qJD(1) * t33;
t165 = qJD(2) * t92;
t164 = qJD(4) * t48;
t163 = qJD(4) * t49;
t160 = qJD(4) * t97;
t158 = qJD(5) * t94;
t154 = qJDD(4) * t95;
t153 = qJDD(4) * t97;
t146 = t196 * pkin(2) + t174;
t144 = t95 * t100 * t97;
t143 = t63 * t156;
t142 = t63 * t161;
t140 = 0.2e1 * t150;
t65 = t91 * t150;
t135 = -t91 * t151 + t60 * t92;
t54 = -t91 * qJ(2) + t92 * t98;
t134 = qJDD(1) * t171;
t133 = qJDD(2) - t168;
t50 = pkin(3) - t54;
t132 = t95 * t137;
t31 = t135 - t65;
t131 = -t195 * pkin(1) + t196 * qJ(2);
t129 = -pkin(4) * t95 + pkin(7) * t97;
t128 = -g(1) * t47 + t203;
t51 = -pkin(6) + t55;
t126 = -qJD(5) * t51 * t97 + qJD(2) * t91 + t129 * qJD(4);
t125 = t7 * t96 + t8 * t94;
t124 = -t7 * t94 + t8 * t96;
t121 = t25 * t95 - t26 * t97;
t120 = t38 * t91 - t39 * t92;
t119 = -qJD(5) * t22 + t203;
t27 = qJDD(1) * pkin(3) - t31;
t118 = -t46 * pkin(3) + pkin(6) * t47 + t146;
t45 = -qJDD(5) - t206;
t117 = t63 * t157 - t45 * t94;
t116 = t63 * t158 + t45 * t96;
t115 = g(1) * t196 + g(2) * t195;
t114 = -t91 * t160 + t92 * t166;
t112 = -t127 + t167;
t111 = -t195 * pkin(2) + t131;
t21 = -qJD(4) * pkin(4) - t25;
t110 = pkin(7) * t45 + t63 * t21;
t109 = t127 * t97 + t201;
t5 = qJDD(4) * pkin(7) + t9;
t108 = -qJD(5) * t23 - t97 * t202 - t201 - t5;
t35 = t130 + t50;
t106 = qJD(5) * t35 - t51 * t162 + t97 * t165;
t105 = t47 * pkin(3) + t46 * pkin(6) + t111;
t12 = t206 * pkin(4) + t207 * pkin(7) + t27;
t1 = qJD(5) * t7 + t94 * t12 + t96 * t5;
t11 = t96 * t12;
t2 = -qJD(5) * t8 - t94 * t5 + t11;
t104 = -t125 * qJD(5) + t1 * t96 - t2 * t94;
t103 = -qJDD(4) * t51 + (-qJD(1) * t50 - t165 - t33) * qJD(4);
t101 = qJDD(1) * t50 - t51 * t99 + t128 + t27 + t65;
t58 = -t95 * t99 + t153;
t57 = -t97 * t99 - t154;
t53 = t129 * qJD(1);
t20 = -t46 * t178 + t190;
t19 = t46 * t181 + t189;
t16 = t51 * t178 + t35 * t94;
t15 = -t51 * t181 + t35 * t96;
t14 = t25 * t96 + t53 * t94;
t13 = -t25 * t94 + t53 * t96;
t4 = -t106 * t94 + t126 * t96;
t3 = t106 * t96 + t126 * t94;
t24 = [0, 0, 0, 0, 0, qJDD(1), t173, t115, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t168 + t173, 0, -t115 + t140 + 0.2e1 * t151, -t133 * pkin(1) - g(1) * t131 - g(2) * t174 + (t140 + t151) * qJ(2), 0, 0, 0, 0, 0, qJDD(1), -qJDD(1) * t54 + t128 - t135 + 0.2e1 * t65, qJDD(1) * t55 + t127 + t175 + 0.2e1 * t67, 0, -g(1) * t111 - g(2) * t146 - t120 * qJD(2) + t31 * t54 + t32 * t55, qJDD(1) * t88 + 0.2e1 * t132, -0.2e1 * t172 * t149 + 0.2e1 * t95 * t80, t57, qJDD(1) * t89 - 0.2e1 * t132, -t58, 0, t101 * t97 + t103 * t95, -t101 * t95 + t103 * t97, -t51 * t134 - t171 * t67 - t102 - t127, t27 * t50 - g(1) * t105 - g(2) * t118 + (-t121 * t92 + t33 * t91) * qJD(2) + t102 * t51, -t96 * t180 + (-t95 * t158 + t96 * t160) * t49, (-t183 - t185) * t160 + (t194 - t193 + (-t182 + t186) * qJD(5)) * t95, (t17 - t142) * t97 + (t116 + t163) * t95, t113 * t48 + t94 * t179, (t18 + t143) * t97 + (t117 - t164) * t95, -t63 * t162 - t45 * t97, -t94 * t204 - g(2) * t20 - t15 * t45 + t4 * t63 + (-g(1) * t189 + t2 + (-t21 * t94 - t48 * t51) * qJD(4)) * t97 + (-t7 * qJD(4) - t157 * t21 - t165 * t48 - t51 * t18 - t6 * t94) * t95, -t96 * t204 - g(2) * t19 + t16 * t45 - t3 * t63 + (g(1) * t190 - t1 + (-t21 * t96 - t49 * t51) * qJD(4)) * t97 + (t8 * qJD(4) + t158 * t21 - t165 * t49 + t51 * t17 - t6 * t96) * t95, -t15 * t17 + t16 * t18 + t3 * t48 + t4 * t49 + t125 * t160 + (qJD(5) * t124 + t1 * t94 + t2 * t96 + t128) * t95, t1 * t16 + t8 * t3 + t2 * t15 + t7 * t4 + t51 * t199 - g(1) * (t130 * t47 + t105) - g(2) * (-t130 * t46 + t118) + (t51 * t160 + t165 * t95) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t100, -qJ(2) * t100 + t133 - t173, 0, 0, 0, 0, 0, 0, -qJDD(1) * t92 - t100 * t91, qJDD(1) * t91 - t170, 0, qJD(1) * t120 + t31 * t92 + t32 * t91 - t173, 0, 0, 0, 0, 0, 0, (0.2e1 * t138 - t80) * t92 + (-t169 * t97 - t154) * t91, (0.2e1 * t137 + t152) * t92 + (t169 * t95 - t153) * t91, -t91 * t134 + t171 * t170, (qJD(1) * t121 - t27) * t92 + (t102 - t167) * t91 - t173, 0, 0, 0, 0, 0, 0, t114 * t48 + t176 * t63 - t179 * t91 - t41 * t45, t114 * t49 - t177 * t63 + t180 * t91 + t42 * t45, -t17 * t41 + t176 * t49 + t177 * t48 + t18 * t42, t1 * t42 - t114 * t21 + t176 * t7 + t177 * t8 + t91 * t199 + t2 * t41 - t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + g(3), 0, 0, 0, 0, 0, 0, t58, t57, 0, -qJD(4) * t121 + t10 * t97 + t9 * t95 + g(3), 0, 0, 0, 0, 0, 0, (t18 - t143) * t97 + (-t117 - t164) * t95, (-t17 - t142) * t97 + (t116 - t163) * t95, (-t183 + t185) * t160 + (t194 + t193 + (-t182 - t186) * qJD(5)) * t95, g(3) + (qJD(4) * t124 - t6) * t97 + (qJD(4) * t21 + t104) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t172 * t100, -t152, t144, -t80, qJDD(4), t200 + t79 + (t26 - t191) * qJD(4) + (t112 + t208) * t95, -t201 + (t25 + t192) * qJD(4) + t112 * t97 + t147, 0, 0, -t63 * t182 + t194, (t17 + t187) * t96 + (t18 + t184) * t94, (t63 * t178 - t49 * t95) * qJD(1) + t117, -t63 * t186 + t193, (-t63 * t181 + t48 * t95) * qJD(1) - t116, t63 * t166, pkin(4) * t18 + t110 * t94 - t13 * t63 + t7 * t166 + t205 * t96 + t26 * t48, -pkin(4) * t17 + t110 * t96 + t14 * t63 - t8 * t166 - t205 * t94 + t26 * t49, -t13 * t49 - t14 * t48 + (t1 - t198 + (-qJD(5) * t49 + t18) * pkin(7)) * t96 + (-t2 - t197 + (t17 - t159) * pkin(7)) * t94 + t109, -t7 * t13 - t8 * t14 - t21 * t26 + t107 * pkin(4) + (t104 + t109) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, -t48 ^ 2 + t49 ^ 2, t17 - t187, -t188, t18 - t184, -t45, -g(1) * t19 + t108 * t94 + t119 * t96 + t21 * t49 + t11 + t197, g(1) * t20 - t21 * t48 + t198 + (-t119 - t12) * t94 + t108 * t96, 0, 0;];
tau_reg = t24;

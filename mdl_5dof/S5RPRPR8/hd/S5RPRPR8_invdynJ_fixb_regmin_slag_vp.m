% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR8
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
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:13
% EndTime: 2019-12-31 18:22:20
% DurationCPUTime: 2.43s
% Computational Cost: add. (1823->296), mult. (3883->424), div. (0->0), fcn. (2662->14), ass. (0->161)
t128 = cos(qJ(3));
t125 = sin(qJ(3));
t117 = qJ(1) + pkin(8);
t111 = sin(t117);
t113 = cos(t117);
t162 = g(1) * t113 + g(2) * t111;
t147 = t162 * t125;
t137 = -g(3) * t128 + t147;
t121 = sin(pkin(8));
t104 = pkin(1) * t121 + pkin(6);
t100 = t104 * qJD(1);
t183 = qJD(3) * t128;
t98 = t104 * qJDD(1);
t210 = -qJD(2) * qJD(3) - t98;
t175 = -t100 * t183 + t210 * t125;
t144 = -qJDD(3) * pkin(3) + qJDD(4) - t175;
t23 = -qJDD(2) * t128 + t144;
t211 = t23 - t137;
t187 = qJD(1) * t128;
t102 = -qJD(5) + t187;
t124 = sin(qJ(5));
t127 = cos(qJ(5));
t120 = sin(pkin(9));
t122 = cos(pkin(9));
t180 = t122 * qJD(3);
t188 = qJD(1) * t125;
t81 = t120 * t188 - t180;
t185 = qJD(3) * t120;
t83 = t122 * t188 + t185;
t152 = t124 * t81 - t127 * t83;
t209 = t102 * t152;
t68 = qJD(2) * t128 - t125 * t100;
t190 = qJDD(2) - g(3);
t208 = t190 * t128;
t114 = t128 * qJDD(1);
t179 = qJD(1) * qJD(3);
t207 = t125 * t179 - t114;
t206 = -qJD(5) - t102;
t108 = t122 * qJDD(3);
t168 = t128 * t179;
t176 = t125 * qJDD(1);
t141 = t168 + t176;
t53 = t141 * t120 - t108;
t177 = qJDD(3) * t120;
t54 = t141 * t122 + t177;
t9 = -t152 * qJD(5) + t124 * t54 + t127 * t53;
t123 = cos(pkin(8));
t205 = pkin(1) * t123;
t204 = g(3) * t125;
t202 = pkin(7) + qJ(4);
t20 = qJDD(3) * qJ(4) + qJDD(2) * t125 + t128 * t98 + (qJD(4) + t68) * qJD(3);
t156 = pkin(3) * t125 - qJ(4) * t128;
t72 = t156 * qJD(3) - qJD(4) * t125;
t150 = pkin(3) * t128 + qJ(4) * t125 + pkin(2);
t78 = -t150 - t205;
t26 = t72 * qJD(1) + t78 * qJDD(1);
t7 = t120 * t26 + t122 * t20;
t87 = t120 * t127 + t122 * t124;
t143 = t87 * t128;
t181 = qJD(5) * t127;
t182 = qJD(5) * t125;
t192 = t122 * t125;
t194 = t120 * t124;
t28 = qJD(3) * t143 + t181 * t192 - t182 * t194;
t64 = t87 * t125;
t85 = qJDD(5) + t207;
t201 = t28 * t102 - t64 * t85;
t69 = t125 * qJD(2) + t128 * t100;
t57 = qJD(3) * qJ(4) + t69;
t60 = t78 * qJD(1);
t15 = t120 * t60 + t122 * t57;
t92 = t156 * qJD(1);
t32 = t120 * t92 + t122 * t68;
t86 = -t127 * t122 + t194;
t142 = t86 * t128;
t200 = qJD(1) * t142 - t86 * qJD(5);
t199 = -qJD(1) * t143 + t87 * qJD(5);
t184 = qJD(3) * t125;
t172 = t104 * t184;
t34 = t120 * t172 + t122 * t72;
t191 = t122 * t128;
t41 = t104 * t191 + t120 * t78;
t197 = t124 * t83;
t36 = t127 * t81 + t197;
t198 = t102 * t36;
t196 = t111 * t128;
t195 = t113 * t128;
t193 = t120 * t128;
t118 = t125 ^ 2;
t189 = -t128 ^ 2 + t118;
t105 = -pkin(2) - t205;
t101 = qJD(1) * t105;
t6 = -t120 * t20 + t122 * t26;
t4 = t207 * pkin(4) - pkin(7) * t54 + t6;
t5 = -pkin(7) * t53 + t7;
t174 = -t124 * t5 + t127 * t4;
t173 = t120 * t187;
t8 = -qJD(5) * t197 - t124 * t53 + t127 * t54 - t81 * t181;
t171 = -t128 * t8 - t152 * t184;
t170 = qJ(4) * t114;
t167 = t120 * t176;
t166 = t122 * t176;
t165 = pkin(4) * t120 + t104;
t14 = -t120 * t57 + t122 * t60;
t31 = -t120 * t68 + t122 * t92;
t163 = qJD(1) * t189;
t161 = g(1) * t111 - g(2) * t113;
t126 = sin(qJ(1));
t129 = cos(qJ(1));
t160 = g(1) * t126 - g(2) * t129;
t159 = -t6 * t120 + t7 * t122;
t158 = t124 * t4 + t127 * t5;
t27 = -qJD(3) * t142 - t87 * t182;
t65 = t86 * t125;
t157 = -t102 * t27 - t65 * t85;
t10 = -pkin(4) * t187 - pkin(7) * t83 + t14;
t12 = -pkin(7) * t81 + t15;
t1 = t10 * t127 - t12 * t124;
t2 = t10 * t124 + t12 * t127;
t155 = -t120 * t14 + t122 * t15;
t67 = t122 * t78;
t25 = -pkin(7) * t192 + t67 + (-t104 * t120 - pkin(4)) * t128;
t30 = -pkin(7) * t120 * t125 + t41;
t154 = -t124 * t30 + t127 * t25;
t153 = t124 * t25 + t127 * t30;
t151 = pkin(4) * t125 - pkin(7) * t191;
t149 = t160 * pkin(1);
t148 = t128 * t9 - t36 * t184;
t97 = t202 * t122;
t146 = t151 * qJD(1) + qJD(4) * t120 + qJD(5) * t97 + t31;
t96 = t202 * t120;
t145 = pkin(7) * t173 + qJD(4) * t122 - qJD(5) * t96 - t32;
t55 = -qJD(3) * pkin(3) + qJD(4) - t68;
t140 = -qJD(1) * t101 + t162;
t139 = -qJ(4) * t184 + (qJD(4) - t55) * t128;
t138 = 0.2e1 * t101 * qJD(3) - qJDD(3) * t104;
t135 = -t162 * t128 - t204;
t130 = qJD(3) ^ 2;
t133 = -0.2e1 * qJDD(1) * t105 - t104 * t130 + t161;
t131 = qJD(1) ^ 2;
t116 = pkin(9) + qJ(5);
t112 = cos(t116);
t110 = sin(t116);
t106 = -pkin(4) * t122 - pkin(3);
t94 = qJDD(3) * t128 - t125 * t130;
t93 = qJDD(3) * t125 + t128 * t130;
t71 = t165 * t125;
t63 = t165 * t183;
t61 = t120 * t72;
t52 = t110 * t111 + t112 * t195;
t51 = -t110 * t195 + t111 * t112;
t50 = t110 * t113 - t112 * t196;
t49 = t110 * t196 + t112 * t113;
t42 = pkin(4) * t173 + t69;
t40 = -t104 * t193 + t67;
t35 = -t122 * t172 + t61;
t33 = pkin(4) * t81 + t55;
t24 = t61 + (-pkin(7) * t193 - t104 * t192) * qJD(3);
t18 = t151 * qJD(3) + t34;
t11 = pkin(4) * t53 + t23;
t3 = [qJDD(1), t160, g(1) * t129 + g(2) * t126, (t121 ^ 2 + t123 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t149, qJDD(1) * t118 + 0.2e1 * t125 * t168, -0.2e1 * qJD(3) * t163 + 0.2e1 * t125 * t114, t93, t94, 0, t138 * t125 + t133 * t128, -t133 * t125 + t138 * t128, -t162 * t120 + (t104 * t53 + t23 * t120 + (qJD(1) * t40 + t14) * qJD(3)) * t125 + (-t34 * qJD(1) - t40 * qJDD(1) - t6 + t161 * t122 + (t104 * t81 + t120 * t55) * qJD(3)) * t128, -t162 * t122 + (t104 * t54 + t23 * t122 + (-qJD(1) * t41 - t15) * qJD(3)) * t125 + (t35 * qJD(1) + t41 * qJDD(1) + t7 - t161 * t120 + (t104 * t83 + t122 * t55) * qJD(3)) * t128, -t34 * t83 - t35 * t81 - t40 * t54 - t41 * t53 + (-t120 * t15 - t122 * t14) * t183 + (-t120 * t7 - t122 * t6 + t161) * t125, t14 * t34 + t15 * t35 + t6 * t40 + t7 * t41 + t149 + (-g(1) * pkin(6) - g(2) * t150) * t113 + (-g(2) * pkin(6) + g(1) * t150) * t111 + (t125 * t23 + t55 * t183) * t104, -t152 * t27 - t65 * t8, t152 * t28 - t27 * t36 - t64 * t8 + t65 * t9, t157 + t171, t148 + t201, -t102 * t184 - t128 * t85, -(-t124 * t24 + t127 * t18) * t102 + t154 * t85 - t174 * t128 + t1 * t184 + t63 * t36 + t71 * t9 + t11 * t64 + t33 * t28 - g(1) * t50 - g(2) * t52 + (t102 * t153 + t128 * t2) * qJD(5), (t124 * t18 + t127 * t24) * t102 - t153 * t85 + t158 * t128 - t2 * t184 - t63 * t152 + t71 * t8 - t11 * t65 + t33 * t27 - g(1) * t49 - g(2) * t51 + (t1 * t128 + t102 * t154) * qJD(5); 0, 0, 0, t190, 0, 0, 0, 0, 0, t94, -t93, (-t53 + t167) * t128 + (-t120 * t163 + t125 * t81) * qJD(3), (-t54 + t166) * t128 + (-t122 * t163 + t125 * t83) * qJD(3), (t120 * t54 - t122 * t53) * t125 + (t120 * t83 - t122 * t81) * t183, -t128 * t23 - g(3) + t159 * t125 + (t125 * t55 + t155 * t128) * qJD(3), 0, 0, 0, 0, 0, -t148 + t201, -t157 + t171; 0, 0, 0, 0, -t125 * t131 * t128, t189 * t131, t176, t114, qJDD(3), qJD(3) * t69 + t140 * t125 + t175 + t208, qJD(3) * t68 + (qJD(3) * t100 - t190) * t125 + (t140 + t210) * t128, t120 * t170 - pkin(3) * t53 - t69 * t81 - t211 * t122 + (t139 * t120 - t125 * t14 + t128 * t31) * qJD(1), t122 * t170 - pkin(3) * t54 - t69 * t83 + t211 * t120 + (t139 * t122 + t125 * t15 - t128 * t32) * qJD(1), t31 * t83 + t32 * t81 + (-qJ(4) * t53 - qJD(4) * t81 + t14 * t187 + t7) * t122 + (qJ(4) * t54 + qJD(4) * t83 + t15 * t187 - t6) * t120 + t135, -t14 * t31 - t15 * t32 - t55 * t69 + t155 * qJD(4) - t211 * pkin(3) + (t135 + t159) * qJ(4), -t152 * t200 + t8 * t87, t152 * t199 - t200 * t36 - t8 * t86 - t87 * t9, -t200 * t102 + t152 * t188 + t85 * t87, t199 * t102 + t36 * t188 - t85 * t86, t102 * t188, (-t124 * t97 - t127 * t96) * t85 + t106 * t9 + t11 * t86 - t1 * t188 - t42 * t36 + t199 * t33 + (t124 * t145 + t127 * t146) * t102 + t137 * t112, -(-t124 * t96 + t127 * t97) * t85 + t106 * t8 + t11 * t87 + t2 * t188 + t42 * t152 + t200 * t33 + (-t124 * t146 + t127 * t145) * t102 - t137 * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t167 - t108 + (-t83 + t185) * t187, t166 + t177 + (t81 + t180) * t187, -t81 ^ 2 - t83 ^ 2, t14 * t83 + t15 * t81 + t144 - t147 - t208, 0, 0, 0, 0, 0, t9 + t209, t8 + t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152 * t36, t152 ^ 2 - t36 ^ 2, t8 - t198, -t9 + t209, t85, -g(1) * t51 + g(2) * t49 + t110 * t204 + t152 * t33 + t206 * t2 + t174, g(1) * t52 - g(2) * t50 + t206 * t1 + t112 * t204 + t33 * t36 - t158;];
tau_reg = t3;

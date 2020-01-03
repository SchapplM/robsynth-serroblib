% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:02:08
% EndTime: 2020-01-03 12:02:11
% DurationCPUTime: 1.08s
% Computational Cost: add. (1356->197), mult. (2183->273), div. (0->0), fcn. (1448->16), ass. (0->150)
t110 = qJD(1) + qJD(2);
t117 = sin(qJ(5));
t122 = cos(qJ(4));
t118 = sin(qJ(4));
t121 = cos(qJ(5));
t175 = t121 * t118;
t61 = t117 * t122 + t175;
t46 = t61 * t110;
t114 = qJ(1) + qJ(2);
t99 = pkin(9) + t114;
t89 = sin(t99);
t90 = cos(t99);
t144 = -g(2) * t90 - g(3) * t89;
t108 = qJDD(1) + qJDD(2);
t115 = sin(pkin(9));
t116 = cos(pkin(9));
t123 = cos(qJ(2));
t185 = pkin(1) * qJD(2);
t162 = qJD(1) * t185;
t119 = sin(qJ(2));
t166 = qJDD(1) * t119;
t199 = pkin(1) * t166 + t123 * t162;
t189 = t123 * pkin(1);
t97 = qJDD(1) * t189;
t47 = t108 * pkin(2) - t119 * t162 + t97;
t27 = -t199 * t115 + t116 * t47;
t23 = -t108 * pkin(3) - t27;
t137 = t144 - t23;
t198 = g(2) * t89 - g(3) * t90;
t186 = pkin(1) * qJD(1);
t163 = t123 * t186;
t66 = t110 * pkin(2) + t163;
t165 = t119 * t186;
t86 = t116 * t165;
t40 = t115 * t66 + t86;
t155 = t40 + (pkin(7) + pkin(8)) * t110;
t25 = t122 * qJD(3) - t155 * t118;
t109 = qJD(4) + qJD(5);
t26 = t118 * qJD(3) + t155 * t122;
t179 = t116 * t119;
t96 = pkin(2) + t189;
t187 = pkin(1) * t179 + t115 * t96;
t50 = pkin(7) + t187;
t193 = -pkin(8) - t50;
t91 = t115 * pkin(2) + pkin(7);
t192 = -pkin(8) - t91;
t191 = t116 * pkin(2);
t190 = t122 * pkin(4);
t174 = t121 * t122;
t159 = t110 * t174;
t178 = t117 * t118;
t161 = t110 * t178;
t44 = -t159 + t161;
t188 = t46 * t44;
t113 = qJ(4) + qJ(5);
t101 = sin(t113);
t184 = t101 * t89;
t183 = t101 * t90;
t182 = t121 * t26;
t181 = t110 * t118;
t180 = t115 * t119;
t176 = t118 * t108;
t173 = t122 * t108;
t172 = qJDD(3) - g(1);
t111 = t118 ^ 2;
t171 = -t122 ^ 2 + t111;
t170 = qJD(5) * t117;
t168 = t118 * qJD(4);
t167 = t122 * qJD(4);
t28 = t115 * t47 + t199 * t116;
t164 = pkin(4) * t168;
t158 = -pkin(3) - t190;
t157 = t110 * t167;
t24 = t108 * pkin(7) + t28;
t156 = pkin(8) * t108 + t24;
t102 = sin(t114);
t104 = cos(t114);
t154 = g(2) * t102 - g(3) * t104;
t153 = qJD(4) * t193;
t152 = qJD(4) * t192;
t12 = (t110 * t168 - t173) * pkin(4) + t23;
t85 = t115 * t165;
t39 = t116 * t66 - t85;
t29 = t158 * t110 - t39;
t60 = -t174 + t178;
t32 = t109 * t60;
t151 = g(2) * t183 + g(3) * t184 + t12 * t61 - t29 * t32;
t150 = -pkin(1) * t180 + t116 * t96;
t35 = -t110 * pkin(3) - t39;
t149 = -t137 * t118 + t35 * t167;
t148 = qJD(1) * (-qJD(2) + t110);
t147 = qJD(2) * (-qJD(1) - t110);
t49 = -pkin(3) - t150;
t51 = t115 * t163 + t86;
t145 = -t51 + t164;
t143 = t117 * t176 - t121 * t173;
t142 = -g(2) * t104 - g(3) * t102;
t107 = qJDD(4) + qJDD(5);
t16 = t61 * t107 - t32 * t109;
t22 = qJD(4) * pkin(4) + t25;
t140 = -t117 * t22 - t182;
t37 = t193 * t118;
t105 = t122 * pkin(8);
t38 = t122 * t50 + t105;
t139 = -t117 * t38 + t121 * t37;
t138 = t117 * t37 + t121 * t38;
t136 = t142 + t97;
t53 = t116 * t163 - t85;
t59 = t122 * t91 + t105;
t135 = qJD(5) * t59 - t118 * t53 - t122 * t152;
t58 = t192 * t118;
t134 = -qJD(5) * t58 - t118 * t152 + t122 * t53;
t125 = qJD(4) ^ 2;
t52 = (t115 * t123 + t179) * t185;
t133 = t108 * t49 + t110 * t52 + t125 * t50;
t92 = -pkin(3) - t191;
t132 = t108 * t92 - t110 * t51 + t125 * t91;
t131 = -t35 * t110 + t198 - t24;
t103 = cos(t113);
t33 = t109 * t61;
t130 = t144 * t103 + t12 * t60 + t29 * t33;
t54 = (t116 * t123 - t180) * t185;
t129 = -qJDD(4) * t50 + (t110 * t49 - t54) * qJD(4);
t128 = -qJDD(4) * t91 + (t110 * t92 + t53) * qJD(4);
t13 = qJD(5) * t159 + t108 * t175 - t109 * t161 + t117 * t173 + t121 * t157;
t98 = t122 * qJDD(3);
t4 = qJDD(4) * pkin(4) - t26 * qJD(4) - t156 * t118 + t98;
t127 = t29 * t44 + t26 * t170 + g(1) * t101 + (-t26 * t109 - t4) * t117 + t198 * t103;
t5 = t25 * qJD(4) + t118 * qJDD(3) + t156 * t122;
t126 = -g(1) * t103 + g(2) * t184 - g(3) * t183 + t140 * qJD(5) - t117 * t5 + t121 * t4 - t29 * t46;
t124 = cos(qJ(1));
t120 = sin(qJ(1));
t106 = t110 ^ 2;
t71 = qJDD(4) * t122 - t125 * t118;
t70 = qJDD(4) * t118 + t125 * t122;
t69 = t158 - t191;
t48 = t111 * t108 + 0.2e1 * t118 * t157;
t43 = t49 - t190;
t41 = t52 + t164;
t34 = -0.2e1 * t171 * t110 * qJD(4) + 0.2e1 * t118 * t173;
t30 = t35 * t168;
t21 = -t118 * t54 + t122 * t153;
t20 = t118 * t153 + t122 * t54;
t17 = -t60 * t107 - t33 * t109;
t15 = -t44 ^ 2 + t46 ^ 2;
t14 = t33 * t110 + t143;
t6 = t44 * t109 + t13;
t2 = t13 * t61 - t46 * t32;
t1 = -t13 * t60 - t61 * t14 + t32 * t44 - t46 * t33;
t3 = [qJDD(1), -g(2) * t124 - g(3) * t120, g(2) * t120 - g(3) * t124, t108, (t108 * t123 + t119 * t147) * pkin(1) + t136, ((-qJDD(1) - t108) * t119 + t123 * t147) * pkin(1) + t154, t28 * t187 + t40 * t54 + t27 * t150 - t39 * t52 - g(2) * (t124 * pkin(1) + pkin(2) * t104) - g(3) * (t120 * pkin(1) + pkin(2) * t102), t48, t34, t70, t71, 0, t30 + t129 * t118 + (-t133 + t137) * t122, t133 * t118 + t129 * t122 + t149, t2, t1, t16, t17, 0, t41 * t44 + t43 * t14 + (-t138 * qJD(5) - t117 * t20 + t121 * t21) * t109 + t139 * t107 + t130, t41 * t46 + t43 * t13 - (t139 * qJD(5) + t117 * t21 + t121 * t20) * t109 - t138 * t107 + t151; 0, 0, 0, t108, t119 * pkin(1) * t148 + t136, (t123 * t148 - t166) * pkin(1) + t154, t39 * t51 - t40 * t53 + (t115 * t28 + t116 * t27 + t142) * pkin(2), t48, t34, t70, t71, 0, t30 + t128 * t118 + (-t132 + t137) * t122, t132 * t118 + t128 * t122 + t149, t2, t1, t16, t17, 0, t69 * t14 + (-t117 * t59 + t121 * t58) * t107 + t145 * t44 + (t134 * t117 - t135 * t121) * t109 + t130, t69 * t13 - (t117 * t58 + t121 * t59) * t107 + t145 * t46 + (t135 * t117 + t134 * t121) * t109 + t151; 0, 0, 0, 0, 0, 0, t172, 0, 0, 0, 0, 0, t71, -t70, 0, 0, 0, 0, 0, t17, -t16; 0, 0, 0, 0, 0, 0, 0, -t118 * t106 * t122, t171 * t106, t176, t173, qJDD(4), -g(1) * t122 + t131 * t118 + t98, -t172 * t118 + t131 * t122, t188, t15, t6, -t143, t107, -(-t117 * t25 - t182) * t109 + (t121 * t107 - t109 * t170 - t44 * t181) * pkin(4) + t126, (-qJD(5) * t22 + t25 * t109 - t5) * t121 + (-qJD(5) * t121 * t109 - t117 * t107 - t46 * t181) * pkin(4) + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, t15, t6, -t143, t107, -t140 * t109 + t126, (-t5 + (-qJD(5) + t109) * t22) * t121 + t127;];
tau_reg = t3;

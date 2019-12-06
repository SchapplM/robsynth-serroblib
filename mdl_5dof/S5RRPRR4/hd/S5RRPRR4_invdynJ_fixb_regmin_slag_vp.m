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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:32:10
% EndTime: 2019-12-05 18:32:13
% DurationCPUTime: 1.11s
% Computational Cost: add. (1356->199), mult. (2183->274), div. (0->0), fcn. (1448->16), ass. (0->149)
t112 = qJD(1) + qJD(2);
t119 = sin(qJ(5));
t124 = cos(qJ(4));
t120 = sin(qJ(4));
t123 = cos(qJ(5));
t175 = t123 * t120;
t61 = t119 * t124 + t175;
t46 = t61 * t112;
t125 = cos(qJ(2));
t185 = pkin(1) * qJD(2);
t160 = qJD(1) * t185;
t121 = sin(qJ(2));
t166 = qJDD(1) * t121;
t202 = pkin(1) * t166 + t125 * t160;
t116 = qJ(1) + qJ(2);
t101 = pkin(9) + t116;
t89 = sin(t101);
t90 = cos(t101);
t201 = -g(2) * t89 + g(3) * t90;
t200 = g(2) * t90 + g(3) * t89;
t104 = sin(t116);
t106 = cos(t116);
t199 = g(2) * t106 + g(3) * t104;
t117 = sin(pkin(9));
t186 = pkin(1) * qJD(1);
t161 = t125 * t186;
t66 = t112 * pkin(2) + t161;
t118 = cos(pkin(9));
t163 = t121 * t186;
t86 = t118 * t163;
t40 = t117 * t66 + t86;
t153 = t40 + (pkin(7) + pkin(8)) * t112;
t25 = t124 * qJD(3) - t153 * t120;
t111 = qJD(4) + qJD(5);
t26 = t120 * qJD(3) + t153 * t124;
t179 = t118 * t121;
t190 = t125 * pkin(1);
t98 = pkin(2) + t190;
t187 = pkin(1) * t179 + t117 * t98;
t50 = pkin(7) + t187;
t194 = -pkin(8) - t50;
t91 = t117 * pkin(2) + pkin(7);
t193 = -pkin(8) - t91;
t192 = t118 * pkin(2);
t191 = t124 * pkin(4);
t174 = t123 * t124;
t157 = t112 * t174;
t178 = t119 * t120;
t159 = t112 * t178;
t44 = -t157 + t159;
t189 = t46 * t44;
t167 = t124 * qJD(4);
t110 = qJDD(1) + qJDD(2);
t99 = qJDD(1) * t190;
t47 = t110 * pkin(2) - t121 * t160 + t99;
t27 = -t202 * t117 + t118 * t47;
t23 = -t110 * pkin(3) - t27;
t85 = t117 * t163;
t39 = t118 * t66 - t85;
t35 = -t112 * pkin(3) - t39;
t188 = t23 * t120 + t35 * t167;
t115 = qJ(4) + qJ(5);
t105 = cos(t115);
t184 = t105 * t89;
t183 = t105 * t90;
t182 = t123 * t26;
t181 = t112 * t120;
t180 = t117 * t121;
t176 = t120 * t110;
t173 = t124 * t110;
t172 = qJDD(3) - g(1);
t113 = t120 ^ 2;
t171 = -t124 ^ 2 + t113;
t170 = qJD(5) * t119;
t168 = t120 * qJD(4);
t165 = t200 * t124 + t35 * t168;
t28 = t117 * t47 + t202 * t118;
t164 = t99 + t199;
t162 = pkin(4) * t168;
t156 = -pkin(3) - t191;
t155 = t112 * t167;
t24 = t110 * pkin(7) + t28;
t154 = pkin(8) * t110 + t24;
t152 = -g(2) * t104 + g(3) * t106;
t151 = qJD(4) * t194;
t150 = qJD(4) * t193;
t12 = (t112 * t168 - t173) * pkin(4) + t23;
t29 = t156 * t112 - t39;
t33 = t111 * t61;
t60 = -t174 + t178;
t149 = g(2) * t183 + g(3) * t184 + t12 * t60 + t29 * t33;
t148 = -pkin(1) * t180 + t118 * t98;
t147 = qJD(1) * (-qJD(2) + t112);
t146 = qJD(2) * (-qJD(1) - t112);
t49 = -pkin(3) - t148;
t51 = t117 * t161 + t86;
t144 = -t51 + t162;
t142 = t119 * t176 - t123 * t173;
t109 = qJDD(4) + qJDD(5);
t32 = t111 * t60;
t16 = t61 * t109 - t32 * t111;
t22 = qJD(4) * pkin(4) + t25;
t140 = -t119 * t22 - t182;
t37 = t194 * t120;
t107 = t124 * pkin(8);
t38 = t124 * t50 + t107;
t139 = -t119 * t38 + t123 * t37;
t138 = t119 * t37 + t123 * t38;
t53 = t118 * t161 - t85;
t59 = t124 * t91 + t107;
t137 = qJD(5) * t59 - t120 * t53 - t124 * t150;
t58 = t193 * t120;
t136 = -qJD(5) * t58 - t120 * t150 + t124 * t53;
t127 = qJD(4) ^ 2;
t52 = (t117 * t125 + t179) * t185;
t135 = -t110 * t49 - t112 * t52 - t127 * t50;
t92 = -pkin(3) - t192;
t134 = -t110 * t92 + t112 * t51 - t127 * t91;
t133 = -t35 * t112 + t201 - t24;
t103 = sin(t115);
t132 = -t103 * t200 + t12 * t61 - t29 * t32;
t54 = (t118 * t125 - t180) * t185;
t131 = -qJDD(4) * t50 + (t112 * t49 - t54) * qJD(4);
t130 = -qJDD(4) * t91 + (t112 * t92 + t53) * qJD(4);
t13 = qJD(5) * t157 + t110 * t175 - t111 * t159 + t119 * t173 + t123 * t155;
t100 = t124 * qJDD(3);
t4 = qJDD(4) * pkin(4) - t26 * qJD(4) - t154 * t120 + t100;
t129 = -g(2) * t184 + t29 * t44 + t26 * t170 + g(3) * t183 + g(1) * t103 + (-t26 * t111 - t4) * t119;
t5 = t25 * qJD(4) + t120 * qJDD(3) + t154 * t124;
t128 = -g(1) * t105 + t140 * qJD(5) + t201 * t103 - t119 * t5 + t123 * t4 - t29 * t46;
t126 = cos(qJ(1));
t122 = sin(qJ(1));
t108 = t112 ^ 2;
t71 = qJDD(4) * t124 - t127 * t120;
t70 = qJDD(4) * t120 + t127 * t124;
t69 = t156 - t192;
t48 = t113 * t110 + 0.2e1 * t120 * t155;
t43 = t49 - t191;
t41 = t52 + t162;
t34 = -0.2e1 * t171 * t112 * qJD(4) + 0.2e1 * t120 * t173;
t21 = -t120 * t54 + t124 * t151;
t20 = t120 * t151 + t124 * t54;
t17 = -t60 * t109 - t33 * t111;
t15 = -t44 ^ 2 + t46 ^ 2;
t14 = t33 * t112 + t142;
t6 = t44 * t111 + t13;
t2 = t13 * t61 - t46 * t32;
t1 = -t13 * t60 - t61 * t14 + t32 * t44 - t46 * t33;
t3 = [qJDD(1), g(2) * t126 + g(3) * t122, -g(2) * t122 + g(3) * t126, t110, (t110 * t125 + t121 * t146) * pkin(1) + t164, ((-qJDD(1) - t110) * t121 + t125 * t146) * pkin(1) + t152, t28 * t187 + t40 * t54 + t27 * t148 - t39 * t52 - g(2) * (-t126 * pkin(1) - pkin(2) * t106) - g(3) * (-t122 * pkin(1) - pkin(2) * t104), t48, t34, t70, t71, 0, t131 * t120 + (t135 - t23) * t124 + t165, t131 * t124 + (-t135 - t200) * t120 + t188, t2, t1, t16, t17, 0, t41 * t44 + t43 * t14 + (-t138 * qJD(5) - t119 * t20 + t123 * t21) * t111 + t139 * t109 + t149, t41 * t46 + t43 * t13 - (t139 * qJD(5) + t119 * t21 + t123 * t20) * t111 - t138 * t109 + t132; 0, 0, 0, t110, t121 * pkin(1) * t147 + t164, (t125 * t147 - t166) * pkin(1) + t152, t39 * t51 - t40 * t53 + (t117 * t28 + t118 * t27 + t199) * pkin(2), t48, t34, t70, t71, 0, t130 * t120 + (t134 - t23) * t124 + t165, t130 * t124 + (-t134 - t200) * t120 + t188, t2, t1, t16, t17, 0, t69 * t14 + (-t119 * t59 + t123 * t58) * t109 + t144 * t44 + (t136 * t119 - t137 * t123) * t111 + t149, t69 * t13 - (t119 * t58 + t123 * t59) * t109 + t144 * t46 + (t137 * t119 + t136 * t123) * t111 + t132; 0, 0, 0, 0, 0, 0, t172, 0, 0, 0, 0, 0, t71, -t70, 0, 0, 0, 0, 0, t17, -t16; 0, 0, 0, 0, 0, 0, 0, -t120 * t108 * t124, t171 * t108, t176, t173, qJDD(4), -g(1) * t124 + t133 * t120 + t100, -t172 * t120 + t133 * t124, t189, t15, t6, -t142, t109, -(-t119 * t25 - t182) * t111 + (t123 * t109 - t111 * t170 - t44 * t181) * pkin(4) + t128, (-qJD(5) * t22 + t25 * t111 - t5) * t123 + (-qJD(5) * t123 * t111 - t119 * t109 - t46 * t181) * pkin(4) + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t15, t6, -t142, t109, -t140 * t111 + t128, (-t5 + (-qJD(5) + t111) * t22) * t123 + t129;];
tau_reg = t3;

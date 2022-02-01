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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:48:32
% EndTime: 2022-01-20 10:48:37
% DurationCPUTime: 1.14s
% Computational Cost: add. (1356->201), mult. (2183->278), div. (0->0), fcn. (1448->16), ass. (0->151)
t114 = qJD(1) + qJD(2);
t121 = sin(qJ(5));
t126 = cos(qJ(4));
t122 = sin(qJ(4));
t125 = cos(qJ(5));
t176 = t125 * t122;
t61 = t121 * t126 + t176;
t46 = t61 * t114;
t112 = qJDD(1) + qJDD(2);
t119 = sin(pkin(9));
t120 = cos(pkin(9));
t127 = cos(qJ(2));
t188 = pkin(1) * qJD(2);
t162 = qJD(1) * t188;
t123 = sin(qJ(2));
t167 = qJDD(1) * t123;
t203 = pkin(1) * t167 + t127 * t162;
t194 = t127 * pkin(1);
t101 = qJDD(1) * t194;
t47 = t112 * pkin(2) - t123 * t162 + t101;
t27 = -t203 * t119 + t120 * t47;
t23 = -t112 * pkin(3) - t27;
t118 = qJ(1) + qJ(2);
t103 = pkin(9) + t118;
t92 = cos(t103);
t204 = g(2) * t92 + t23;
t106 = sin(t118);
t108 = cos(t118);
t202 = g(1) * t106 - g(2) * t108;
t189 = pkin(1) * qJD(1);
t163 = t127 * t189;
t68 = t114 * pkin(2) + t163;
t165 = t123 * t189;
t88 = t120 * t165;
t40 = t119 * t68 + t88;
t154 = t40 + (pkin(7) + pkin(8)) * t114;
t25 = t126 * qJD(3) - t154 * t122;
t113 = qJD(4) + qJD(5);
t26 = t122 * qJD(3) + t154 * t126;
t91 = sin(t103);
t201 = g(1) * t91;
t100 = pkin(2) + t194;
t180 = t120 * t123;
t191 = pkin(1) * t180 + t119 * t100;
t50 = pkin(7) + t191;
t199 = -pkin(8) - t50;
t93 = t119 * pkin(2) + pkin(7);
t198 = -pkin(8) - t93;
t196 = t120 * pkin(2);
t195 = t126 * pkin(4);
t175 = t125 * t126;
t159 = t114 * t175;
t179 = t121 * t122;
t161 = t114 * t179;
t44 = -t159 + t161;
t193 = t46 * t44;
t169 = t122 * qJD(4);
t87 = t119 * t165;
t39 = t120 * t68 - t87;
t35 = -t114 * pkin(3) - t39;
t192 = t126 * t201 + t35 * t169;
t190 = g(1) * t108 + g(2) * t106;
t117 = qJ(4) + qJ(5);
t105 = sin(t117);
t187 = t105 * t91;
t186 = t105 * t92;
t107 = cos(t117);
t185 = t107 * t91;
t184 = t107 * t92;
t183 = t125 * t26;
t182 = t114 * t122;
t181 = t119 * t123;
t177 = t122 * t112;
t174 = t126 * t112;
t173 = qJDD(3) - g(3);
t115 = t122 ^ 2;
t172 = -t126 ^ 2 + t115;
t171 = qJD(5) * t121;
t168 = t126 * qJD(4);
t166 = t204 * t122 + t35 * t168;
t28 = t119 * t47 + t203 * t120;
t164 = pkin(4) * t169;
t158 = -pkin(3) - t195;
t156 = t114 * t168;
t24 = t112 * pkin(7) + t28;
t155 = pkin(8) * t112 + t24;
t153 = qJD(4) * t199;
t152 = qJD(4) * t198;
t151 = -pkin(1) * t181 + t120 * t100;
t150 = qJD(1) * (-qJD(2) + t114);
t149 = qJD(2) * (-qJD(1) - t114);
t49 = -pkin(3) - t151;
t147 = t101 + t202;
t51 = t119 * t163 + t88;
t146 = -t51 + t164;
t145 = t121 * t177 - t125 * t174;
t111 = qJDD(4) + qJDD(5);
t60 = -t175 + t179;
t32 = t113 * t60;
t16 = t61 * t111 - t32 * t113;
t22 = qJD(4) * pkin(4) + t25;
t143 = -t121 * t22 - t183;
t37 = t199 * t122;
t109 = t126 * pkin(8);
t38 = t126 * t50 + t109;
t142 = -t121 * t38 + t125 * t37;
t141 = t121 * t37 + t125 * t38;
t12 = (t114 * t169 - t174) * pkin(4) + t23;
t29 = t158 * t114 - t39;
t140 = -g(1) * t187 + g(2) * t186 + t12 * t61 - t29 * t32;
t33 = t113 * t61;
t139 = g(1) * t185 - g(2) * t184 + t12 * t60 + t29 * t33;
t53 = t120 * t163 - t87;
t59 = t126 * t93 + t109;
t138 = qJD(5) * t59 - t122 * t53 - t126 * t152;
t58 = t198 * t122;
t137 = -qJD(5) * t58 - t122 * t152 + t126 * t53;
t129 = qJD(4) ^ 2;
t52 = (t119 * t127 + t180) * t188;
t136 = t112 * t49 + t114 * t52 + t129 * t50;
t94 = -pkin(3) - t196;
t135 = t112 * t94 - t114 * t51 + t129 * t93;
t134 = g(1) * t92 + g(2) * t91 - t35 * t114 - t24;
t54 = (t120 * t127 - t181) * t188;
t133 = -qJDD(4) * t50 + (t114 * t49 - t54) * qJD(4);
t132 = -qJDD(4) * t93 + (t114 * t94 + t53) * qJD(4);
t13 = qJD(5) * t159 + t112 * t176 - t113 * t161 + t121 * t174 + t125 * t156;
t102 = t126 * qJDD(3);
t4 = qJDD(4) * pkin(4) - t26 * qJD(4) - t155 * t122 + t102;
t131 = t29 * t44 + t26 * t171 + g(2) * t185 + g(1) * t184 + g(3) * t105 + (-t26 * t113 - t4) * t121;
t5 = t25 * qJD(4) + t122 * qJDD(3) + t155 * t126;
t130 = g(1) * t186 + g(2) * t187 - g(3) * t107 + t143 * qJD(5) - t121 * t5 + t125 * t4 - t29 * t46;
t128 = cos(qJ(1));
t124 = sin(qJ(1));
t110 = t114 ^ 2;
t73 = qJDD(4) * t126 - t129 * t122;
t72 = qJDD(4) * t122 + t129 * t126;
t71 = t158 - t196;
t48 = t115 * t112 + 0.2e1 * t122 * t156;
t43 = t49 - t195;
t41 = t52 + t164;
t34 = -0.2e1 * t172 * t114 * qJD(4) + 0.2e1 * t122 * t174;
t21 = -t122 * t54 + t126 * t153;
t20 = t122 * t153 + t126 * t54;
t17 = -t60 * t111 - t33 * t113;
t15 = -t44 ^ 2 + t46 ^ 2;
t14 = t114 * t33 + t145;
t6 = t44 * t113 + t13;
t2 = t13 * t61 - t46 * t32;
t1 = -t13 * t60 - t61 * t14 + t32 * t44 - t46 * t33;
t3 = [qJDD(1), g(1) * t124 - g(2) * t128, g(1) * t128 + g(2) * t124, t112, (t112 * t127 + t123 * t149) * pkin(1) + t147, ((-qJDD(1) - t112) * t123 + t127 * t149) * pkin(1) + t190, t28 * t191 + t40 * t54 + t27 * t151 - t39 * t52 - g(1) * (-t124 * pkin(1) - pkin(2) * t106) - g(2) * (t128 * pkin(1) + pkin(2) * t108), t48, t34, t72, t73, 0, t133 * t122 + (-t136 - t204) * t126 + t192, t133 * t126 + (t136 - t201) * t122 + t166, t2, t1, t16, t17, 0, t41 * t44 + t43 * t14 + (-t141 * qJD(5) - t121 * t20 + t125 * t21) * t113 + t142 * t111 + t139, t41 * t46 + t43 * t13 - (t142 * qJD(5) + t121 * t21 + t125 * t20) * t113 - t141 * t111 + t140; 0, 0, 0, t112, t123 * pkin(1) * t150 + t147, (t127 * t150 - t167) * pkin(1) + t190, t39 * t51 - t40 * t53 + (t119 * t28 + t120 * t27 + t202) * pkin(2), t48, t34, t72, t73, 0, t132 * t122 + (-t135 - t204) * t126 + t192, t132 * t126 + (t135 - t201) * t122 + t166, t2, t1, t16, t17, 0, t71 * t14 + (-t121 * t59 + t125 * t58) * t111 + t146 * t44 + (t137 * t121 - t138 * t125) * t113 + t139, t71 * t13 - (t121 * t58 + t125 * t59) * t111 + t146 * t46 + (t138 * t121 + t137 * t125) * t113 + t140; 0, 0, 0, 0, 0, 0, t173, 0, 0, 0, 0, 0, t73, -t72, 0, 0, 0, 0, 0, t17, -t16; 0, 0, 0, 0, 0, 0, 0, -t122 * t110 * t126, t172 * t110, t177, t174, qJDD(4), -g(3) * t126 + t134 * t122 + t102, -t173 * t122 + t134 * t126, t193, t15, t6, -t145, t111, -(-t121 * t25 - t183) * t113 + (t125 * t111 - t113 * t171 - t44 * t182) * pkin(4) + t130, (-qJD(5) * t22 + t25 * t113 - t5) * t125 + (-qJD(5) * t125 * t113 - t121 * t111 - t46 * t182) * pkin(4) + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t15, t6, -t145, t111, -t143 * t113 + t130, (-t5 + (-qJD(5) + t113) * t22) * t125 + t131;];
tau_reg = t3;

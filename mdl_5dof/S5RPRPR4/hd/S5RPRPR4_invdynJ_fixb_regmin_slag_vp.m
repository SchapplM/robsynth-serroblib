% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR4
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
% Datum: 2021-01-15 11:45
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:44:14
% EndTime: 2021-01-15 11:44:21
% DurationCPUTime: 1.41s
% Computational Cost: add. (1636->229), mult. (3529->314), div. (0->0), fcn. (2526->16), ass. (0->140)
t123 = sin(qJ(5));
t126 = cos(qJ(5));
t118 = sin(pkin(9));
t120 = cos(pkin(9));
t127 = cos(qJ(3));
t167 = t120 * t127;
t157 = qJD(1) * t167;
t124 = sin(qJ(3));
t163 = qJD(1) * t124;
t74 = t118 * t163 - t157;
t83 = t118 * t127 + t120 * t124;
t77 = t83 * qJD(1);
t143 = t123 * t74 - t126 * t77;
t159 = t127 * qJDD(1);
t160 = t124 * qJDD(1);
t148 = t118 * t160 - t120 * t159;
t76 = t83 * qJD(3);
t42 = qJD(1) * t76 + t148;
t161 = qJD(1) * qJD(3);
t156 = t124 * t161;
t137 = t83 * qJDD(1) - t118 * t156;
t155 = t127 * t161;
t43 = t120 * t155 + t137;
t133 = t143 * qJD(5) - t123 * t43 - t126 * t42;
t113 = qJD(3) + qJD(5);
t169 = t143 * t113;
t187 = t133 - t169;
t186 = 0.2e1 * qJD(3);
t66 = t126 * t74;
t35 = -t123 * t77 - t66;
t168 = t35 * t113;
t162 = qJD(5) * t123;
t7 = -qJD(5) * t66 - t123 * t42 + t126 * t43 - t77 * t162;
t185 = t7 - t168;
t184 = t143 * t35;
t119 = sin(pkin(8));
t101 = t119 * pkin(1) + pkin(6);
t166 = qJ(4) + t101;
t115 = qJ(1) + pkin(8);
t106 = sin(t115);
t108 = cos(t115);
t150 = g(2) * t106 - g(3) * t108;
t183 = t143 ^ 2 - t35 ^ 2;
t114 = qJ(3) + pkin(9);
t110 = qJ(5) + t114;
t100 = cos(t110);
t179 = t74 * pkin(7);
t153 = t166 * qJD(1);
t60 = t124 * qJD(2) + t153 * t127;
t170 = t120 * t60;
t171 = qJD(3) * pkin(3);
t59 = t127 * qJD(2) - t153 * t124;
t53 = t59 + t171;
t24 = t118 * t53 + t170;
t11 = t24 - t179;
t104 = t127 * pkin(3) + pkin(2);
t121 = cos(pkin(8));
t172 = t121 * pkin(1);
t86 = -t104 - t172;
t72 = t86 * qJD(1) + qJD(4);
t44 = t74 * pkin(4) + t72;
t99 = sin(t110);
t182 = g(1) * t99 + t150 * t100 + t11 * t162 - t44 * t35;
t109 = t127 * qJDD(2);
t89 = t101 * qJDD(1);
t134 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(3) * qJD(2) + t89;
t142 = t153 * qJD(3);
t21 = qJDD(3) * pkin(3) - t134 * t124 - t127 * t142 + t109;
t25 = (qJDD(2) - t142) * t124 + t134 * t127;
t4 = -t118 * t25 + t120 * t21;
t2 = qJDD(3) * pkin(4) - t43 * pkin(7) + t4;
t5 = t118 * t21 + t120 * t25;
t3 = -t42 * pkin(7) + t5;
t181 = -g(1) * t100 - t123 * t3 + t126 * t2 + t44 * t143 + t150 * t99;
t180 = qJD(5) - t113;
t178 = t77 * pkin(7);
t177 = pkin(3) * t118;
t176 = pkin(3) * t124;
t175 = g(1) * t127;
t49 = t118 * t60;
t27 = t120 * t59 - t49;
t152 = qJD(3) * t166;
t63 = t127 * qJD(4) - t124 * t152;
t64 = -t124 * qJD(4) - t127 * t152;
t29 = t118 * t64 + t120 * t63;
t80 = t166 * t124;
t81 = t166 * t127;
t41 = -t118 * t80 + t120 * t81;
t165 = qJDD(2) - g(1);
t116 = t124 ^ 2;
t164 = -t127 ^ 2 + t116;
t103 = -pkin(2) - t172;
t92 = qJD(1) * t103;
t158 = t124 * t171;
t23 = t120 * t53 - t49;
t26 = -t118 * t59 - t170;
t28 = -t118 * t63 + t120 * t64;
t40 = -t118 * t81 - t120 * t80;
t151 = g(2) * t108 + g(3) * t106;
t125 = sin(qJ(1));
t128 = cos(qJ(1));
t149 = -g(2) * t128 - g(3) * t125;
t10 = qJD(3) * pkin(4) - t178 + t23;
t147 = -t123 * t10 - t126 * t11;
t112 = qJDD(3) + qJDD(5);
t82 = t118 * t124 - t167;
t45 = t123 * t83 + t126 * t82;
t79 = t82 * qJD(3);
t12 = -t45 * qJD(5) - t123 * t76 - t126 * t79;
t46 = -t123 * t82 + t126 * t83;
t146 = t46 * t112 + t12 * t113;
t30 = -t83 * pkin(7) + t40;
t31 = -t82 * pkin(7) + t41;
t145 = -t123 * t31 + t126 * t30;
t144 = t123 * t30 + t126 * t31;
t102 = t120 * pkin(3) + pkin(4);
t140 = t123 * t102 + t126 * t177;
t139 = t126 * t102 - t123 * t177;
t138 = -t92 * qJD(1) + t150 - t89;
t136 = -qJDD(3) * t101 + t92 * t186;
t58 = pkin(3) * t156 + t86 * qJDD(1) + qJDD(4);
t129 = qJD(3) ^ 2;
t132 = 0.2e1 * qJDD(1) * t103 + t101 * t129 + t151;
t130 = qJD(1) ^ 2;
t122 = -qJ(4) - pkin(6);
t107 = cos(t114);
t105 = sin(t114);
t88 = qJDD(3) * t127 - t129 * t124;
t87 = qJDD(3) * t124 + t129 * t127;
t62 = t76 * pkin(4) + t158;
t61 = pkin(3) * t163 + t77 * pkin(4);
t57 = t82 * pkin(4) + t86;
t22 = t42 * pkin(4) + t58;
t17 = -t76 * pkin(7) + t29;
t16 = t79 * pkin(7) + t28;
t15 = t27 - t178;
t14 = t26 + t179;
t13 = t46 * qJD(5) - t123 * t79 + t126 * t76;
t6 = -t45 * t112 - t13 * t113;
t1 = [qJDD(1), t149, g(2) * t125 - g(3) * t128, (t149 + (t119 ^ 2 + t121 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t116 * qJDD(1) + 0.2e1 * t124 * t155, 0.2e1 * t124 * t159 - 0.2e1 * t164 * t161, t87, t88, 0, t136 * t124 - t132 * t127, t132 * t124 + t136 * t127, t40 * qJDD(3) + t86 * t42 + t58 * t82 + t72 * t76 - t151 * t107 + (t74 * t176 + t28) * qJD(3), -t41 * qJDD(3) + t86 * t43 + t58 * t83 - t72 * t79 + t151 * t105 + (t77 * t176 - t29) * qJD(3), t23 * t79 - t24 * t76 - t28 * t77 - t29 * t74 - t4 * t83 - t40 * t43 - t41 * t42 - t5 * t82 - t150, t5 * t41 + t24 * t29 + t4 * t40 + t23 * t28 + t58 * t86 + t72 * t158 - g(2) * (t128 * pkin(1) + t108 * t104 - t106 * t122) - g(3) * (t125 * pkin(1) + t106 * t104 + t108 * t122), -t12 * t143 + t7 * t46, t12 * t35 + t13 * t143 + t133 * t46 - t7 * t45, t146, t6, 0, -t62 * t35 - t57 * t133 + t22 * t45 + t44 * t13 + (-qJD(5) * t144 - t123 * t17 + t126 * t16) * t113 + t145 * t112 - t151 * t100, -t62 * t143 + t57 * t7 + t22 * t46 + t44 * t12 - (qJD(5) * t145 + t123 * t16 + t126 * t17) * t113 - t144 * t112 + t151 * t99; 0, 0, 0, t165, 0, 0, 0, 0, 0, t88, -t87, -t76 * qJD(3) - t82 * qJDD(3), t79 * qJD(3) - t83 * qJDD(3), -t83 * t42 + t82 * t43 + t79 * t74 + t76 * t77, -t23 * t76 - t24 * t79 - t4 * t82 + t5 * t83 - g(1), 0, 0, 0, 0, 0, t6, -t146; 0, 0, 0, 0, -t124 * t130 * t127, t164 * t130, t160, t159, qJDD(3), t138 * t124 + t109 - t175, -t165 * t124 + t138 * t127, -g(1) * t107 - t26 * qJD(3) - t72 * t77 + t150 * t105 + (qJDD(3) * t120 - t74 * t163) * pkin(3) + t4, g(1) * t105 + t27 * qJD(3) + t72 * t74 + t150 * t107 + (-qJDD(3) * t118 - t77 * t163) * pkin(3) - t5, (t24 + t26) * t77 + (-t23 + t27) * t74 + (-t118 * t42 - t120 * t43) * pkin(3), -t23 * t26 - t24 * t27 + (-t175 + t118 * t5 + t120 * t4 + (-qJD(1) * t72 + t150) * t124) * pkin(3), t184, t183, t185, t187, t112, t139 * t112 + t61 * t35 - (-t123 * t15 + t126 * t14) * t113 + (-t113 * t140 + t147) * qJD(5) + t181, -t140 * t112 - t126 * t3 - t123 * t2 + t61 * t143 + (t123 * t14 + t126 * t15) * t113 + (-t126 * t10 - t113 * t139) * qJD(5) + t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t186 + t148, (-t74 + t157) * qJD(3) + t137, -t74 ^ 2 - t77 ^ 2, t23 * t77 + t24 * t74 + t151 + t58, 0, 0, 0, 0, 0, -t133 - t169, t7 + t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t183, t185, t187, t112, t180 * t147 + t181, (-t11 * t113 - t2) * t123 + (-t180 * t10 - t3) * t126 + t182;];
tau_reg = t1;

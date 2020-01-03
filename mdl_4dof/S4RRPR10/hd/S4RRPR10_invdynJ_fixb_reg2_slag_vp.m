% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:58
% EndTime: 2019-12-31 17:12:01
% DurationCPUTime: 1.64s
% Computational Cost: add. (1354->287), mult. (2947->368), div. (0->0), fcn. (1646->6), ass. (0->162)
t133 = qJD(1) * qJD(2);
t86 = sin(qJ(2));
t126 = t86 * t133;
t89 = cos(qJ(2));
t134 = t89 * qJDD(1);
t188 = -t126 + t134;
t125 = t89 * t133;
t135 = t86 * qJDD(1);
t102 = t125 + t135;
t69 = pkin(5) * t135;
t127 = pkin(5) * t125 + qJDD(3) + t69;
t181 = pkin(2) + pkin(6);
t14 = t102 * pkin(3) - t181 * qJDD(2) + t127;
t77 = t86 * qJ(3);
t124 = -pkin(1) - t77;
t101 = -t181 * t89 + t124;
t63 = pkin(2) * t126;
t151 = qJ(3) * t89;
t112 = pkin(6) * t86 - t151;
t139 = t86 * qJD(3);
t98 = t112 * qJD(2) - t139;
t5 = t98 * qJD(1) + t101 * qJDD(1) + t63;
t20 = t101 * qJD(1);
t140 = t86 * qJD(1);
t71 = pkin(5) * t140;
t137 = pkin(3) * t140 + qJD(3) + t71;
t23 = -t181 * qJD(2) + t137;
t85 = sin(qJ(4));
t88 = cos(qJ(4));
t7 = t88 * t20 + t85 * t23;
t2 = -qJD(4) * t7 + t88 * t14 - t85 * t5;
t64 = qJD(4) + t140;
t187 = t7 * t64 + t2;
t180 = pkin(3) + pkin(5);
t186 = t180 * t86;
t141 = t85 * qJD(2);
t149 = qJD(1) * t89;
t38 = t88 * t149 + t141;
t109 = t38 * t64;
t11 = t38 * qJD(4) - t88 * qJDD(2) + t188 * t85;
t185 = t11 - t109;
t128 = t85 * t149;
t12 = -qJD(4) * t128 + t85 * qJDD(2) + (qJD(2) * qJD(4) + t188) * t88;
t138 = t88 * qJD(2);
t40 = -t128 + t138;
t168 = t40 * t64;
t184 = -t12 + t168;
t145 = qJD(4) * t85;
t37 = qJDD(4) + t102;
t27 = t88 * t37;
t183 = -t64 * t145 + t27;
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t115 = g(1) * t90 + g(2) * t87;
t136 = qJD(2) * qJ(3);
t72 = pkin(5) * t149;
t51 = -t72 - t136;
t73 = pkin(3) * t149;
t28 = -t51 + t73;
t182 = -t181 * t37 + t28 * t64;
t179 = g(1) * t87;
t176 = g(2) * t90;
t175 = g(3) * t86;
t82 = g(3) * t89;
t6 = -t85 * t20 + t88 * t23;
t174 = t6 * t64;
t79 = t89 * pkin(2);
t172 = t89 * pkin(6);
t171 = t12 * t85;
t169 = t40 * t38;
t167 = t64 * t86;
t166 = t85 * t37;
t165 = t86 * t87;
t164 = t86 * t90;
t93 = qJD(1) ^ 2;
t163 = t86 * t93;
t162 = t87 * t85;
t161 = t87 * t88;
t160 = t88 * t11;
t159 = t88 * t89;
t158 = t89 * t90;
t157 = t90 * t85;
t156 = t90 * t88;
t155 = t79 + t77;
t154 = t90 * pkin(1) + t87 * pkin(5);
t83 = t86 ^ 2;
t84 = t89 ^ 2;
t153 = t83 - t84;
t152 = t83 + t84;
t150 = pkin(5) * qJDD(2);
t148 = qJD(2) * t38;
t147 = qJD(2) * t40;
t146 = qJD(2) * t86;
t144 = qJD(4) * t88;
t143 = qJD(4) * t89;
t142 = qJDD(2) * pkin(2);
t132 = qJDD(2) * qJ(3);
t131 = -g(1) * t164 - g(2) * t165 + t82;
t53 = t180 * t89;
t123 = pkin(2) * t158 + t90 * t77 + t154;
t122 = -t69 - t131;
t121 = -qJD(2) * pkin(2) + qJD(3);
t70 = pkin(5) * t134;
t120 = -t70 - t132;
t119 = t155 + t172;
t118 = t86 * t125;
t117 = t152 * qJDD(1) * pkin(5);
t92 = qJD(2) ^ 2;
t116 = pkin(5) * t92 + t176;
t114 = -t6 * t85 + t7 * t88;
t34 = -pkin(1) - t119;
t17 = t186 * t85 + t88 * t34;
t16 = t186 * t88 - t85 * t34;
t46 = t121 + t71;
t111 = t46 * t89 + t51 * t86;
t110 = t64 * t85;
t108 = t124 - t79;
t24 = t127 - t142;
t105 = -t64 * t144 - t166;
t104 = -0.2e1 * pkin(1) * t133 - t150;
t103 = -t89 * t136 - t139;
t100 = 0.2e1 * qJDD(1) * pkin(1) - t116;
t29 = t108 * qJD(1);
t47 = -pkin(1) - t155;
t99 = t150 + (-qJD(1) * t47 - t29) * qJD(2);
t13 = t103 * qJD(1) + t108 * qJDD(1) + t63;
t74 = pkin(2) * t146;
t26 = t103 + t74;
t97 = qJD(1) * t26 + qJDD(1) * t47 + t116 + t13;
t21 = (-qJD(3) + t71) * qJD(2) + t120;
t96 = t111 * qJD(2) - t21 * t89 + t24 * t86;
t15 = pkin(3) * t134 + (-qJD(1) * t186 + qJD(3)) * qJD(2) - t120;
t95 = qJD(4) * t181 * t64 - t115 * t89 + t15 - t175;
t80 = t90 * pkin(5);
t75 = pkin(2) * t140;
t67 = t89 * t179;
t61 = t90 * t151;
t59 = t87 * t151;
t58 = t89 * t163;
t50 = t153 * t93;
t49 = qJDD(2) * t89 - t92 * t86;
t48 = qJDD(2) * t86 + t92 * t89;
t45 = qJD(2) * t53;
t44 = t72 + t73;
t43 = qJD(2) * t186;
t41 = -qJ(3) * t149 + t75;
t36 = t84 * qJDD(1) - 0.2e1 * t118;
t35 = t83 * qJDD(1) + 0.2e1 * t118;
t33 = -t86 * t162 + t156;
t32 = t86 * t161 + t157;
t31 = t86 * t157 + t161;
t30 = t86 * t156 - t162;
t25 = t112 * qJD(1) + t75;
t22 = t29 * t140;
t19 = t74 + t98;
t18 = -0.2e1 * t153 * t133 + 0.2e1 * t86 * t134;
t10 = t88 * t25 + t85 * t44;
t9 = -t85 * t25 + t88 * t44;
t4 = -t17 * qJD(4) - t85 * t19 + t88 * t45;
t3 = t16 * qJD(4) + t88 * t19 + t85 * t45;
t1 = qJD(4) * t6 + t85 * t14 + t88 * t5;
t8 = [0, 0, 0, 0, 0, qJDD(1), -t176 + t179, t115, 0, 0, t35, t18, t48, t36, t49, 0, t100 * t89 + t104 * t86 + t67, t104 * t89 + (-t100 - t179) * t86, -t115 + 0.2e1 * t117, -g(1) * (-t87 * pkin(1) + t80) - g(2) * t154 + (t152 * pkin(5) ^ 2 + pkin(1) ^ 2) * qJDD(1), 0, -t48, -t49, t35, t18, t36, t117 + t96 - t115, t99 * t86 + t97 * t89 - t67, t99 * t89 + (-t97 + t179) * t86, pkin(5) * t96 - g(1) * t80 - g(2) * t123 - t108 * t179 + t13 * t47 + t29 * t26, t11 * t85 * t89 + (t141 * t86 - t143 * t88) * t40, (-t38 * t85 + t40 * t88) * t146 + (t160 + t171 + (t38 * t88 + t40 * t85) * qJD(4)) * t89, (t141 * t64 - t11) * t86 + (t105 + t147) * t89, t12 * t159 + (-t138 * t86 - t143 * t85) * t38, (t138 * t64 - t12) * t86 + (-t148 - t183) * t89, t64 * qJD(2) * t89 + t37 * t86, -g(1) * t33 - g(2) * t31 + t53 * t12 + t16 * t37 - t43 * t38 + t4 * t64 + (-t138 * t28 + t2) * t86 + (qJD(2) * t6 - t145 * t28 + t15 * t88) * t89, g(1) * t32 - g(2) * t30 - t53 * t11 - t17 * t37 - t3 * t64 - t43 * t40 + (t141 * t28 - t1) * t86 + (-qJD(2) * t7 - t144 * t28 - t15 * t85) * t89, t16 * t11 - t17 * t12 - t3 * t38 - t4 * t40 + t67 + t114 * t146 + (-t176 - t1 * t88 + t2 * t85 + (t6 * t88 + t7 * t85) * qJD(4)) * t89, t1 * t17 + t7 * t3 + t2 * t16 + t6 * t4 + t15 * t53 - t28 * t43 - g(1) * (t90 * pkin(3) + t80) - g(2) * (pkin(6) * t158 + t123) + (-g(1) * (t108 - t172) - g(2) * pkin(3)) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t50, t135, t58, t134, qJDD(2), pkin(1) * t163 + t122, t175 - t70 + (pkin(1) * t93 + t115) * t89, 0, 0, qJDD(2), -t135, -t134, -t58, t50, t58, (-pkin(2) * t86 + t151) * qJDD(1) + ((-t51 - t136) * t86 + (t121 - t46) * t89) * qJD(1), -t41 * t149 + qJDD(3) - t122 - 0.2e1 * t142 + t22, 0.2e1 * t132 + 0.2e1 * qJD(2) * qJD(3) + t70 + (qJD(1) * t41 - g(3)) * t86 + (qJD(1) * t29 - t115) * t89, -t21 * qJ(3) - t51 * qJD(3) - t24 * pkin(2) - t29 * t41 - g(1) * (-pkin(2) * t164 + t61) - g(2) * (-pkin(2) * t165 + t59) - g(3) * t155 - t111 * qJD(1) * pkin(5), -t110 * t40 - t160, (-t12 - t168) * t88 + (t11 + t109) * t85, (-t167 * t85 - t40 * t89) * qJD(1) + t183, t109 * t88 + t171, (-t167 * t88 + t38 * t89) * qJD(1) + t105, -t64 * t149, qJ(3) * t12 + t137 * t38 - t6 * t149 + t182 * t88 - t9 * t64 + t95 * t85, -qJ(3) * t11 + t10 * t64 + t137 * t40 + t7 * t149 - t182 * t85 + t95 * t88, t10 * t38 + t9 * t40 + (-t7 * t140 - t11 * t181 - t2 + (t181 * t38 - t7) * qJD(4)) * t88 + (t6 * t140 + t12 * t181 - t1 + (-t181 * t40 + t6) * qJD(4)) * t85 - t131, -g(1) * t61 - g(2) * t59 - g(3) * t119 + t15 * qJ(3) - t7 * t10 + t137 * t28 - t6 * t9 + (-qJD(4) * t114 - t1 * t85 + t115 * t86 - t2 * t88) * t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, qJDD(2) + t58, -t83 * t93 - t92, t51 * qJD(2) + t131 + t22 + t24, 0, 0, 0, 0, 0, 0, -t110 * t64 - t148 + t27, -t64 ^ 2 * t88 - t147 - t166, t184 * t85 + t185 * t88, -t28 * qJD(2) + t187 * t88 + (t1 - t174) * t85 + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, -t38 ^ 2 + t40 ^ 2, -t185, -t169, t184, t37, -g(1) * t30 - g(2) * t32 + g(3) * t159 - t28 * t40 + t187, g(1) * t31 - g(2) * t33 + t28 * t38 + t174 + (-qJD(4) * t23 - t5) * t88 + (qJD(4) * t20 - t14 - t82) * t85, 0, 0;];
tau_reg = t8;

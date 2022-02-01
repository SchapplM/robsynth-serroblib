% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:30:19
% EndTime: 2022-01-23 09:30:23
% DurationCPUTime: 1.31s
% Computational Cost: add. (1801->226), mult. (3782->284), div. (0->0), fcn. (2524->12), ass. (0->146)
t186 = cos(qJ(4));
t109 = sin(qJ(4));
t110 = sin(qJ(3));
t112 = cos(qJ(3));
t107 = sin(pkin(8));
t86 = t107 * pkin(1) + pkin(6);
t187 = pkin(7) + t86;
t149 = t187 * qJD(1);
t37 = t110 * qJD(2) + t149 * t112;
t31 = t109 * t37;
t175 = qJD(3) * pkin(3);
t36 = t112 * qJD(2) - t149 * t110;
t34 = t36 + t175;
t148 = t186 * t34 - t31;
t59 = t109 * t112 + t186 * t110;
t53 = t59 * qJD(1);
t169 = t53 * qJ(5);
t7 = t148 - t169;
t102 = qJD(3) + qJD(4);
t100 = qJDD(3) + qJDD(4);
t152 = t186 * qJD(4);
t185 = pkin(3) * t102;
t195 = -t109 * pkin(3) * t100 - t152 * t185;
t108 = cos(pkin(8));
t87 = -t108 * pkin(1) - pkin(2);
t99 = t112 * pkin(3);
t194 = t87 - t99;
t71 = t86 * qJDD(1);
t143 = pkin(7) * qJDD(1) + t71;
t95 = t112 * qJDD(2);
t15 = qJDD(3) * pkin(3) - t37 * qJD(3) - t143 * t110 + t95;
t16 = t36 * qJD(3) + t110 * qJDD(2) + t143 * t112;
t193 = -t109 * t16 + t186 * t15;
t167 = t109 * t110;
t133 = t102 * t167;
t155 = t186 * t112;
t140 = qJD(1) * t155;
t145 = qJDD(1) * t186;
t158 = t112 * qJDD(1);
t144 = -t102 * t140 - t109 * t158 - t110 * t145;
t19 = qJD(1) * t133 + t144;
t173 = t19 * qJ(5);
t92 = t100 * pkin(4);
t192 = t173 + t92;
t103 = qJ(1) + pkin(8);
t94 = cos(t103);
t106 = qJ(3) + qJ(4);
t97 = sin(t106);
t181 = t94 * t97;
t93 = sin(t103);
t183 = t93 * t97;
t98 = cos(t106);
t188 = g(3) * t98;
t191 = g(1) * t181 + g(2) * t183 - t188;
t190 = t53 ^ 2;
t5 = t102 * pkin(4) + t7;
t189 = t5 - t7;
t163 = qJD(1) * t110;
t154 = t109 * t163;
t51 = -t140 + t154;
t184 = t53 * t51;
t182 = t93 * t98;
t180 = t94 * t98;
t159 = t110 * qJDD(1);
t135 = t109 * t159 - t112 * t145;
t29 = t102 * t59;
t20 = t29 * qJD(1) + t135;
t28 = -qJD(3) * t155 - t112 * t152 + t133;
t179 = -t59 * t20 + t28 * t51;
t178 = t186 * t36 - t31;
t56 = t187 * t110;
t57 = t187 * t112;
t177 = -t109 * t56 + t186 * t57;
t176 = pkin(4) * t98 + t99;
t172 = t20 * qJ(5);
t171 = t51 * qJ(5);
t170 = t51 * t102;
t74 = qJD(1) * t87;
t150 = t51 * pkin(4) + qJD(5);
t55 = t194 * qJD(1);
t26 = t150 + t55;
t166 = qJD(5) + t26;
t165 = qJDD(2) - g(3);
t104 = t110 ^ 2;
t164 = -t112 ^ 2 + t104;
t162 = qJD(4) * t109;
t160 = qJD(1) * qJD(3);
t157 = t110 * t175;
t156 = pkin(3) * t163;
t33 = t186 * t37;
t153 = qJD(3) * t187;
t151 = t110 * t160;
t147 = -t109 * t36 - t33;
t146 = -t109 * t57 - t186 * t56;
t142 = -g(1) * t183 + g(2) * t181;
t141 = g(1) * t182 - g(2) * t180;
t139 = g(1) * t94 + g(2) * t93;
t138 = -g(1) * t93 + g(2) * t94;
t111 = sin(qJ(1));
t113 = cos(qJ(1));
t137 = g(1) * t111 - g(2) * t113;
t58 = -t155 + t167;
t136 = -t58 * t19 + t53 * t29;
t132 = t59 * t100 - t28 * t102;
t130 = -t109 * t34 - t33;
t48 = t110 * t153;
t49 = t112 * t153;
t129 = -t109 * t49 - t56 * t152 - t57 * t162 - t186 * t48;
t128 = -t74 * qJD(1) + t139 - t71;
t127 = 0.2e1 * t74 * qJD(3) - qJDD(3) * t86;
t39 = pkin(3) * t151 + qJDD(1) * t194;
t126 = -t102 * t154 - t144;
t114 = qJD(3) ^ 2;
t125 = -0.2e1 * qJDD(1) * t87 - t114 * t86 - t138;
t124 = t130 * qJD(4) + t193;
t123 = -t177 * qJD(4) + t109 * t48 - t186 * t49;
t122 = t109 * t15 + t34 * t152 + t186 * t16 - t37 * t162;
t6 = t20 * pkin(4) + qJDD(5) + t39;
t121 = g(1) * t180 + g(2) * t182 + g(3) * t97 - t122;
t120 = t124 + t191;
t119 = t55 * t51 + t121;
t118 = -t55 * t53 + t120;
t117 = t166 * t51 + t121 + t172;
t115 = qJD(1) ^ 2;
t101 = -qJ(5) - pkin(7) - pkin(6);
t91 = t186 * pkin(3) + pkin(4);
t70 = qJDD(3) * t112 - t114 * t110;
t69 = qJDD(3) * t110 + t114 * t112;
t60 = pkin(2) + t176;
t50 = t51 ^ 2;
t38 = t53 * pkin(4) + t156;
t35 = t58 * pkin(4) + t194;
t25 = t29 * pkin(4) + t157;
t23 = -t50 + t190;
t22 = -t58 * qJ(5) + t177;
t21 = -t59 * qJ(5) + t146;
t18 = -t58 * t100 - t29 * t102;
t11 = t126 + t170;
t10 = -t169 + t178;
t9 = t147 + t171;
t8 = -t130 - t171;
t4 = t28 * qJ(5) - t59 * qJD(5) + t123;
t3 = -t29 * qJ(5) - t58 * qJD(5) + t129;
t2 = -t51 * qJD(5) + t122 - t172;
t1 = -t53 * qJD(5) + t124 + t192;
t12 = [qJDD(1), t137, g(1) * t113 + g(2) * t111, (t137 + (t107 ^ 2 + t108 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t104 * qJDD(1) + 0.2e1 * t112 * t151, 0.2e1 * t110 * t158 - 0.2e1 * t164 * t160, t69, t70, 0, t127 * t110 + t125 * t112, -t125 * t110 + t127 * t112, -t19 * t59 - t53 * t28, -t136 + t179, t132, t18, 0, t146 * t100 + t123 * t102 + t51 * t157 + t194 * t20 + t55 * t29 + t39 * t58 + t141, -t177 * t100 - t129 * t102 + t53 * t157 - t19 * t194 - t55 * t28 + t39 * t59 + t142, t21 * t100 + t4 * t102 + t35 * t20 + t25 * t51 + t26 * t29 + t6 * t58 + t141, -t22 * t100 - t3 * t102 - t35 * t19 + t25 * t53 - t26 * t28 + t6 * t59 + t142, -t1 * t59 + t21 * t19 - t2 * t58 - t22 * t20 + t5 * t28 - t8 * t29 - t3 * t51 - t4 * t53 - t139, t2 * t22 + t8 * t3 + t1 * t21 + t5 * t4 + t6 * t35 + t26 * t25 - g(1) * (-t111 * pkin(1) - t94 * t101 - t93 * t60) - g(2) * (t113 * pkin(1) - t93 * t101 + t94 * t60); 0, 0, 0, t165, 0, 0, 0, 0, 0, t70, -t69, 0, 0, 0, 0, 0, t18, -t132, t18, -t132, t136 + t179, -t1 * t58 + t2 * t59 - t8 * t28 - t5 * t29 - g(3); 0, 0, 0, 0, -t110 * t115 * t112, t164 * t115, t159, t158, qJDD(3), -g(3) * t112 + t128 * t110 + t95, -t165 * t110 + t128 * t112, t184, t23, t11, -t135, t100, -t147 * t102 + (t186 * t100 - t102 * t162 - t51 * t163) * pkin(3) + t118, t178 * t102 - t53 * t156 + t119 + t195, t91 * t100 - t9 * t102 - t38 * t51 - t166 * t53 + (-t33 + (-t34 - t185) * t109) * qJD(4) + t191 + t192 + t193, t10 * t102 - t38 * t53 + t117 + t195, t91 * t19 + (t8 + t9) * t53 + (t10 - t5) * t51 + (-t109 * t20 + (t109 * t53 - t186 * t51) * qJD(4)) * pkin(3), t1 * t91 - t8 * t10 - t5 * t9 - t26 * t38 - g(3) * t176 - t139 * (-t110 * pkin(3) - pkin(4) * t97) + (t2 * t109 + (-t109 * t5 + t186 * t8) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, t23, t11, -t135, t100, -t130 * t102 + t118, t148 * t102 + t119, t173 + t8 * t102 + 0.2e1 * t92 + (-t150 - t26) * t53 + t120, -t190 * pkin(4) + t7 * t102 + t117, t19 * pkin(4) - t189 * t51, t189 * t8 + (t139 * t97 - t26 * t53 + t1 - t188) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t102 + t20, t126 - t170, -t50 - t190, t5 * t53 + t8 * t51 + t138 + t6;];
tau_reg = t12;

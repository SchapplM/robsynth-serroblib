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
% Datum: 2021-01-15 12:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
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
% StartTime: 2021-01-15 12:46:06
% EndTime: 2021-01-15 12:46:12
% DurationCPUTime: 1.22s
% Computational Cost: add. (1801->225), mult. (3782->285), div. (0->0), fcn. (2524->12), ass. (0->145)
t183 = cos(qJ(4));
t107 = sin(qJ(4));
t108 = sin(qJ(3));
t110 = cos(qJ(3));
t105 = sin(pkin(8));
t84 = t105 * pkin(1) + pkin(6);
t184 = pkin(7) + t84;
t144 = t184 * qJD(1);
t37 = t108 * qJD(2) + t144 * t110;
t31 = t107 * t37;
t171 = qJD(3) * pkin(3);
t36 = t110 * qJD(2) - t144 * t108;
t34 = t36 + t171;
t143 = t183 * t34 - t31;
t59 = t107 * t110 + t183 * t108;
t53 = t59 * qJD(1);
t165 = t53 * qJ(5);
t7 = t143 - t165;
t100 = qJD(3) + qJD(4);
t147 = t183 * qJD(4);
t182 = pkin(3) * t100;
t98 = qJDD(3) + qJDD(4);
t192 = -t107 * pkin(3) * t98 - t147 * t182;
t106 = cos(pkin(8));
t85 = -t106 * pkin(1) - pkin(2);
t97 = t110 * pkin(3);
t191 = t85 - t97;
t69 = t84 * qJDD(1);
t138 = pkin(7) * qJDD(1) + t69;
t93 = t110 * qJDD(2);
t15 = qJDD(3) * pkin(3) - t37 * qJD(3) - t138 * t108 + t93;
t16 = t36 * qJD(3) + t108 * qJDD(2) + t138 * t110;
t190 = -t107 * t16 + t183 * t15;
t163 = t107 * t108;
t130 = t100 * t163;
t150 = t183 * t110;
t137 = qJD(1) * t150;
t140 = qJDD(1) * t183;
t154 = t110 * qJDD(1);
t139 = -t100 * t137 - t107 * t154 - t108 * t140;
t19 = qJD(1) * t130 + t139;
t169 = t19 * qJ(5);
t90 = t98 * pkin(4);
t189 = t169 + t90;
t101 = qJ(1) + pkin(8);
t92 = cos(t101);
t104 = qJ(3) + qJ(4);
t95 = sin(t104);
t178 = t92 * t95;
t91 = sin(t101);
t180 = t91 * t95;
t96 = cos(t104);
t185 = g(1) * t96;
t188 = g(2) * t180 - g(3) * t178 - t185;
t187 = t53 ^ 2;
t5 = t100 * pkin(4) + t7;
t186 = t5 - t7;
t159 = qJD(1) * t108;
t149 = t107 * t159;
t51 = -t137 + t149;
t181 = t53 * t51;
t179 = t91 * t96;
t177 = t92 * t96;
t155 = t108 * qJDD(1);
t132 = t107 * t155 - t110 * t140;
t29 = t100 * t59;
t20 = t29 * qJD(1) + t132;
t28 = -qJD(3) * t150 - t110 * t147 + t130;
t176 = -t59 * t20 + t28 * t51;
t175 = t183 * t36 - t31;
t56 = t184 * t108;
t57 = t184 * t110;
t174 = -t107 * t56 + t183 * t57;
t173 = g(2) * t178 + g(3) * t180;
t172 = pkin(4) * t96 + t97;
t168 = t20 * qJ(5);
t167 = t51 * qJ(5);
t166 = t51 * t100;
t72 = qJD(1) * t85;
t145 = t51 * pkin(4) + qJD(5);
t55 = t191 * qJD(1);
t26 = t145 + t55;
t162 = qJD(5) + t26;
t161 = qJDD(2) - g(1);
t102 = t108 ^ 2;
t160 = -t110 ^ 2 + t102;
t158 = qJD(4) * t107;
t156 = qJD(1) * qJD(3);
t152 = t108 * t171;
t151 = pkin(3) * t159;
t33 = t183 * t37;
t148 = qJD(3) * t184;
t146 = t108 * t156;
t142 = -t107 * t36 - t33;
t141 = -t107 * t57 - t183 * t56;
t136 = g(2) * t92 + g(3) * t91;
t135 = g(2) * t91 - g(3) * t92;
t109 = sin(qJ(1));
t111 = cos(qJ(1));
t134 = -g(2) * t111 - g(3) * t109;
t58 = -t150 + t163;
t133 = -t58 * t19 + t53 * t29;
t17 = t28 * t100 - t59 * t98;
t128 = -t107 * t34 - t33;
t48 = t108 * t148;
t49 = t110 * t148;
t127 = -t107 * t49 - t56 * t147 - t57 * t158 - t183 * t48;
t126 = -t72 * qJD(1) + t135 - t69;
t125 = 0.2e1 * t72 * qJD(3) - qJDD(3) * t84;
t39 = pkin(3) * t146 + qJDD(1) * t191;
t124 = -t100 * t149 - t139;
t112 = qJD(3) ^ 2;
t123 = 0.2e1 * qJDD(1) * t85 + t112 * t84 + t136;
t122 = t128 * qJD(4) + t190;
t121 = -t174 * qJD(4) + t107 * t48 - t183 * t49;
t120 = t107 * t15 + t34 * t147 - t37 * t158 + t183 * t16;
t6 = t20 * pkin(4) + qJDD(5) + t39;
t119 = g(1) * t95 + g(2) * t179 - g(3) * t177 - t120;
t118 = t122 + t188;
t117 = t55 * t51 + t119;
t116 = -t55 * t53 + t118;
t115 = t162 * t51 + t119 + t168;
t113 = qJD(1) ^ 2;
t99 = qJ(5) + pkin(7) + pkin(6);
t89 = t183 * pkin(3) + pkin(4);
t68 = qJDD(3) * t110 - t112 * t108;
t67 = qJDD(3) * t108 + t112 * t110;
t60 = pkin(2) + t172;
t50 = t51 ^ 2;
t38 = t53 * pkin(4) + t151;
t35 = t58 * pkin(4) + t191;
t25 = t29 * pkin(4) + t152;
t23 = -t50 + t187;
t22 = -t58 * qJ(5) + t174;
t21 = -t59 * qJ(5) + t141;
t18 = -t29 * t100 - t58 * t98;
t11 = t124 + t166;
t10 = -t165 + t175;
t9 = t142 + t167;
t8 = -t128 - t167;
t4 = t28 * qJ(5) - t59 * qJD(5) + t121;
t3 = -t29 * qJ(5) - t58 * qJD(5) + t127;
t2 = -t51 * qJD(5) + t120 - t168;
t1 = -t53 * qJD(5) + t122 + t189;
t12 = [qJDD(1), t134, g(2) * t109 - g(3) * t111, (t134 + (t105 ^ 2 + t106 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t102 * qJDD(1) + 0.2e1 * t110 * t146, 0.2e1 * t108 * t154 - 0.2e1 * t160 * t156, t67, t68, 0, t125 * t108 - t123 * t110, t123 * t108 + t125 * t110, -t19 * t59 - t53 * t28, -t133 + t176, -t17, t18, 0, -g(2) * t177 - g(3) * t179 + t121 * t100 + t141 * t98 + t51 * t152 + t191 * t20 + t55 * t29 + t39 * t58, -t127 * t100 + t53 * t152 - t174 * t98 - t19 * t191 - t55 * t28 + t39 * t59 + t173, t4 * t100 - t136 * t96 + t35 * t20 + t21 * t98 + t25 * t51 + t26 * t29 + t6 * t58, -t3 * t100 - t35 * t19 - t22 * t98 + t25 * t53 - t26 * t28 + t6 * t59 + t173, -t1 * t59 + t21 * t19 - t2 * t58 - t22 * t20 + t5 * t28 - t8 * t29 - t3 * t51 - t4 * t53 - t135, t2 * t22 + t8 * t3 + t1 * t21 + t5 * t4 + t6 * t35 + t26 * t25 - g(2) * (t111 * pkin(1) + t92 * t60 + t91 * t99) - g(3) * (t109 * pkin(1) + t91 * t60 - t92 * t99); 0, 0, 0, t161, 0, 0, 0, 0, 0, t68, -t67, 0, 0, 0, 0, 0, t18, t17, t18, t17, t133 + t176, -t1 * t58 + t2 * t59 - t8 * t28 - t5 * t29 - g(1); 0, 0, 0, 0, -t108 * t113 * t110, t160 * t113, t155, t154, qJDD(3), -g(1) * t110 + t126 * t108 + t93, -t161 * t108 + t126 * t110, t181, t23, t11, -t132, t98, -t142 * t100 + (-t100 * t158 - t51 * t159 + t183 * t98) * pkin(3) + t116, t175 * t100 - t53 * t151 + t117 + t192, -t9 * t100 - t38 * t51 + t89 * t98 - t162 * t53 + (-t33 + (-t34 - t182) * t107) * qJD(4) + t188 + t189 + t190, t10 * t100 - t38 * t53 + t115 + t192, t89 * t19 + (t8 + t9) * t53 + (t10 - t5) * t51 + (-t107 * t20 + (t107 * t53 - t183 * t51) * qJD(4)) * pkin(3), t1 * t89 - t8 * t10 - t5 * t9 - t26 * t38 - g(1) * t172 - t135 * (-t108 * pkin(3) - pkin(4) * t95) + (t2 * t107 + (-t107 * t5 + t183 * t8) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, t23, t11, -t132, t98, -t128 * t100 + t116, t143 * t100 + t117, t169 + t8 * t100 + 0.2e1 * t90 + (-t145 - t26) * t53 + t118, -t187 * pkin(4) + t7 * t100 + t115, t19 * pkin(4) - t186 * t51, t186 * t8 + (t135 * t95 - t26 * t53 + t1 - t185) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t100 + t20, t124 - t166, -t50 - t187, t5 * t53 + t8 * t51 + t136 + t6;];
tau_reg = t12;

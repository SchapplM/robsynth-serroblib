% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPR7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:47
% EndTime: 2019-12-31 17:06:50
% DurationCPUTime: 1.46s
% Computational Cost: add. (2350->302), mult. (5619->412), div. (0->0), fcn. (3891->10), ass. (0->147)
t101 = cos(qJ(1));
t98 = sin(qJ(1));
t177 = g(1) * t98 - g(2) * t101;
t91 = qJ(2) + pkin(7);
t85 = sin(t91);
t179 = t177 * t85;
t97 = sin(qJ(2));
t142 = qJD(1) * t97;
t100 = cos(qJ(2));
t144 = cos(pkin(7));
t124 = t144 * t100;
t78 = qJD(1) * t124;
t94 = sin(pkin(7));
t59 = t142 * t94 - t78;
t52 = qJD(4) + t59;
t96 = sin(qJ(4));
t125 = t52 * t96;
t128 = t144 * t97;
t69 = t100 * t94 + t128;
t62 = t69 * qJD(1);
t99 = cos(qJ(4));
t48 = qJD(2) * t96 + t62 * t99;
t178 = t48 * t125;
t119 = g(1) * t101 + g(2) * t98;
t86 = cos(t91);
t108 = -g(3) * t86 + t119 * t85;
t151 = qJ(3) + pkin(5);
t127 = qJD(2) * t151;
t109 = -t97 * qJD(3) - t100 * t127;
t131 = t151 * t97;
t33 = qJDD(2) * pkin(2) + qJD(1) * t109 - qJDD(1) * t131;
t57 = t100 * qJD(3) - t127 * t97;
t74 = t151 * t100;
t40 = qJD(1) * t57 + qJDD(1) * t74;
t13 = t144 * t33 - t94 * t40;
t11 = -qJDD(2) * pkin(3) - t13;
t81 = pkin(2) * t94 + pkin(6);
t176 = -qJD(4) * t52 * t81 + t108 - t11;
t170 = g(3) * t85;
t175 = -t119 * t86 - t170;
t174 = t62 ^ 2;
t173 = pkin(2) * t97;
t166 = pkin(2) * t100;
t84 = pkin(1) + t166;
t73 = -qJD(1) * t84 + qJD(3);
t20 = t59 * pkin(3) - t62 * pkin(6) + t73;
t72 = qJD(1) * t74;
t65 = t144 * t72;
t147 = qJD(2) * pkin(2);
t71 = qJD(1) * t131;
t67 = -t71 + t147;
t37 = t67 * t94 + t65;
t29 = qJD(2) * pkin(6) + t37;
t6 = t20 * t99 - t29 * t96;
t168 = t6 * t52;
t7 = t20 * t96 + t29 * t99;
t167 = t7 * t52;
t164 = g(3) * t100;
t139 = t99 * qJD(2);
t141 = qJD(4) * t96;
t137 = qJD(1) * qJD(2);
t130 = t97 * t137;
t105 = qJDD(1) * t69 - t130 * t94;
t39 = qJD(2) * t78 + t105;
t16 = -qJD(4) * t139 - qJDD(2) * t96 + t141 * t62 - t39 * t99;
t163 = t16 * t96;
t129 = -qJDD(2) * t99 + t96 * t39;
t17 = qJD(4) * t48 + t129;
t162 = t17 * t99;
t46 = t62 * t96 - t139;
t161 = t46 * t59;
t160 = t46 * t62;
t159 = t48 * t46;
t158 = t48 * t62;
t157 = t62 * t59;
t156 = t69 * t99;
t155 = t94 * t72;
t15 = t96 * t17;
t138 = t97 * qJDD(1);
t118 = -qJDD(1) * t124 + t138 * t94;
t61 = t69 * qJD(2);
t38 = qJD(1) * t61 + t118;
t34 = qJDD(4) + t38;
t154 = t96 * t34;
t153 = t96 * t98;
t152 = t98 * t99;
t27 = t99 * t34;
t140 = qJD(4) * t99;
t150 = -t140 * t46 - t15;
t14 = t144 * t40 + t33 * t94;
t92 = t97 ^ 2;
t93 = t100 ^ 2;
t149 = t92 - t93;
t148 = t92 + t93;
t146 = t101 * t96;
t145 = t101 * t99;
t143 = pkin(5) * qJDD(1);
t136 = t100 * qJDD(1);
t135 = t97 * t147;
t103 = qJD(1) ^ 2;
t133 = t97 * t103 * t100;
t126 = t52 * t99;
t123 = t100 * t130;
t122 = pkin(3) * t86 + pkin(6) * t85;
t121 = -t6 * t99 - t7 * t96;
t111 = -t94 * t97 + t124;
t35 = -pkin(3) * t111 - pkin(6) * t69 - t84;
t45 = -t131 * t94 + t144 * t74;
t18 = t35 * t99 - t45 * t96;
t19 = t35 * t96 + t45 * t99;
t117 = -t125 * t59 - t141 * t52 + t27;
t12 = qJDD(2) * pkin(6) + t14;
t116 = -qJD(4) * t20 - t12 + t170;
t64 = t111 * qJD(2);
t115 = t140 * t69 + t64 * t96;
t114 = -t141 * t69 + t64 * t99;
t36 = t144 * t67 - t155;
t113 = -0.2e1 * pkin(1) * t137 - pkin(5) * qJDD(2);
t28 = -qJD(2) * pkin(3) - t36;
t112 = t28 * t52 - t34 * t81;
t51 = pkin(2) * t130 - qJDD(1) * t84 + qJDD(3);
t102 = qJD(2) ^ 2;
t107 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t102 + t177;
t106 = pkin(1) * t103 + t119 - t143;
t82 = -pkin(2) * t144 - pkin(3);
t76 = t101 * t84;
t58 = t59 ^ 2;
t56 = t145 * t86 + t153;
t55 = -t146 * t86 + t152;
t54 = -t152 * t86 + t146;
t53 = t153 * t86 + t145;
t44 = t128 * t151 + t74 * t94;
t42 = -t144 * t71 - t155;
t41 = -t71 * t94 + t65;
t24 = pkin(3) * t61 - pkin(6) * t64 + t135;
t23 = pkin(2) * t142 + pkin(3) * t62 + pkin(6) * t59;
t22 = t109 * t94 + t144 * t57;
t21 = -t109 * t144 + t57 * t94;
t10 = t23 * t96 + t42 * t99;
t9 = t23 * t99 - t42 * t96;
t8 = t38 * pkin(3) - t39 * pkin(6) + t51;
t5 = t99 * t8;
t4 = -qJD(4) * t19 - t96 * t22 + t99 * t24;
t3 = qJD(4) * t18 + t99 * t22 + t96 * t24;
t2 = -qJD(4) * t7 - t96 * t12 + t5;
t1 = qJD(4) * t6 + t99 * t12 + t96 * t8;
t25 = [0, 0, 0, 0, 0, qJDD(1), t177, t119, 0, 0, qJDD(1) * t92 + 0.2e1 * t123, 0.2e1 * t136 * t97 - 0.2e1 * t137 * t149, qJDD(2) * t97 + t100 * t102, qJDD(1) * t93 - 0.2e1 * t123, qJDD(2) * t100 - t102 * t97, 0, t100 * t107 + t113 * t97, t100 * t113 - t107 * t97, 0.2e1 * t143 * t148 - t119, -g(1) * (-pkin(1) * t98 + pkin(5) * t101) - g(2) * (pkin(1) * t101 + pkin(5) * t98) + (pkin(5) ^ 2 * t148 + pkin(1) ^ 2) * qJDD(1), t39 * t69 + t62 * t64, t111 * t39 - t38 * t69 - t59 * t64 - t61 * t62, qJD(2) * t64 + qJDD(2) * t69, -t111 * t38 + t59 * t61, -qJD(2) * t61 + qJDD(2) * t111, 0, -t44 * qJDD(2) - t84 * t38 - t51 * t111 + t73 * t61 + t177 * t86 + (t173 * t59 - t21) * qJD(2), -t45 * qJDD(2) - t84 * t39 + t51 * t69 + t73 * t64 - t179 + (t173 * t62 - t22) * qJD(2), t111 * t14 - t13 * t69 + t21 * t62 - t22 * t59 - t36 * t64 - t37 * t61 - t38 * t45 + t39 * t44 - t119, t14 * t45 + t37 * t22 - t13 * t44 - t36 * t21 - t51 * t84 + t73 * t135 - g(1) * (t101 * t151 - t84 * t98) - g(2) * (t151 * t98 + t76), t114 * t48 - t156 * t16, (-t46 * t99 - t48 * t96) * t64 + (t163 - t162 + (t46 * t96 - t48 * t99) * qJD(4)) * t69, t111 * t16 + t114 * t52 + t27 * t69 + t48 * t61, t115 * t46 + t15 * t69, t111 * t17 - t115 * t52 - t154 * t69 - t46 * t61, -t111 * t34 + t52 * t61, t11 * t96 * t69 - g(1) * t54 - g(2) * t56 - t111 * t2 + t115 * t28 + t44 * t17 + t18 * t34 + t21 * t46 + t4 * t52 + t6 * t61, -g(1) * t53 - g(2) * t55 + t1 * t111 + t11 * t156 + t114 * t28 - t44 * t16 - t19 * t34 + t21 * t48 - t3 * t52 - t7 * t61, t18 * t16 - t19 * t17 - t3 * t46 - t4 * t48 + t179 + t121 * t64 + (-t1 * t96 - t2 * t99 + (t6 * t96 - t7 * t99) * qJD(4)) * t69, -g(2) * t76 + t1 * t19 + t11 * t44 + t2 * t18 + t28 * t21 + t7 * t3 + t6 * t4 + (-g(1) * t151 - g(2) * t122) * t101 + (-g(1) * (-t122 - t84) - g(2) * t151) * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, t149 * t103, t138, t133, t136, qJDD(2), t106 * t97 - t164, g(3) * t97 + t100 * t106, 0, 0, t157, -t58 + t174, (t78 + t59) * qJD(2) + t105, -t157, -t118, qJDD(2), t41 * qJD(2) - t73 * t62 + (qJDD(2) * t144 - t142 * t59) * pkin(2) + t108 + t13, t42 * qJD(2) + t73 * t59 + (-qJDD(2) * t94 - t142 * t62) * pkin(2) - t14 - t175, (t37 - t41) * t62 + (-t36 + t42) * t59 + (-t144 * t39 - t38 * t94) * pkin(2), t36 * t41 - t37 * t42 + (t144 * t13 - t164 + t14 * t94 + (-qJD(1) * t73 + t119) * t97) * pkin(2), t126 * t48 - t163, (-t16 - t161) * t99 - t178 + t150, t126 * t52 + t154 - t158, t125 * t46 - t162, t117 + t160, -t52 * t62, t112 * t96 + t82 * t17 + t176 * t99 - t41 * t46 - t9 * t52 - t6 * t62, t10 * t52 + t112 * t99 - t82 * t16 - t176 * t96 - t41 * t48 + t7 * t62, t10 * t46 + t9 * t48 + (-t17 * t81 - t59 * t6 + t1 + (t48 * t81 - t6) * qJD(4)) * t99 + (-t16 * t81 - t59 * t7 - t2 + (t46 * t81 - t7) * qJD(4)) * t96 + t175, t11 * t82 - t7 * t10 - t6 * t9 - t28 * t41 - g(3) * (t122 + t166) + (qJD(4) * t121 + t1 * t99 - t2 * t96) * t81 + t119 * (pkin(3) * t85 - pkin(6) * t86 + t173); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t62 * qJD(2) + t118, (t78 - t59) * qJD(2) + t105, -t58 - t174, t36 * t62 + t37 * t59 - t177 + t51, 0, 0, 0, 0, 0, 0, t117 - t160, -t52 ^ 2 * t99 - t154 - t158, (t16 - t161) * t99 + t178 + t150, -t28 * t62 + (t2 + t167) * t99 + (t1 - t168) * t96 - t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, -t46 ^ 2 + t48 ^ 2, t46 * t52 - t16, -t159, -t129 + (-qJD(4) + t52) * t48, t34, -g(1) * t55 + g(2) * t53 + t116 * t96 - t140 * t29 - t28 * t48 + t167 + t5, g(1) * t56 - g(2) * t54 + t28 * t46 + t168 + (qJD(4) * t29 - t8) * t96 + t116 * t99, 0, 0;];
tau_reg = t25;

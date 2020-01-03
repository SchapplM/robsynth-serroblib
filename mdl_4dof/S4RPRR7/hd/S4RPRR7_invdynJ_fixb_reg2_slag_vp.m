% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRR7
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:07
% EndTime: 2019-12-31 16:54:10
% DurationCPUTime: 1.32s
% Computational Cost: add. (2152->256), mult. (5257->342), div. (0->0), fcn. (3824->10), ass. (0->139)
t136 = qJDD(1) * pkin(1);
t94 = sin(qJ(1));
t96 = cos(qJ(1));
t168 = g(1) * t94 - g(2) * t96;
t111 = t168 - qJDD(2) + t136;
t116 = g(1) * t96 + g(2) * t94;
t88 = pkin(7) + qJ(3);
t80 = sin(t88);
t81 = cos(t88);
t101 = -g(3) * t81 + t116 * t80;
t159 = cos(qJ(3));
t132 = qJD(1) * qJD(2);
t143 = pkin(5) + qJ(2);
t166 = t143 * qJDD(1) + t132;
t89 = sin(pkin(7));
t45 = t166 * t89;
t90 = cos(pkin(7));
t46 = t166 * t90;
t93 = sin(qJ(3));
t123 = t159 * t45 + t93 * t46;
t67 = t143 * t89;
t63 = qJD(1) * t67;
t68 = t143 * t90;
t64 = qJD(1) * t68;
t35 = t159 * t64 - t93 * t63;
t12 = -t35 * qJD(3) - t123;
t9 = -qJDD(3) * pkin(3) - t12;
t100 = -t9 + t101;
t148 = t93 * t89;
t128 = qJD(1) * t148;
t127 = t159 * t90;
t74 = qJD(1) * t127;
t54 = -t74 + t128;
t48 = qJD(4) + t54;
t172 = -pkin(6) * qJD(4) * t48 + t100;
t62 = t159 * t89 + t93 * t90;
t56 = t62 * qJD(1);
t77 = t90 * pkin(2) + pkin(1);
t66 = -t77 * qJD(1) + qJD(2);
t20 = t54 * pkin(3) - t56 * pkin(6) + t66;
t28 = qJD(3) * pkin(6) + t35;
t92 = sin(qJ(4));
t95 = cos(qJ(4));
t6 = t92 * t20 + t95 * t28;
t160 = t6 * t48;
t121 = qJDD(1) * t159;
t133 = t90 * qJDD(1);
t130 = qJD(3) * t74 + t89 * t121 + t93 * t133;
t32 = qJD(3) * t128 - t130;
t134 = t89 * qJDD(1);
t113 = -t90 * t121 + t93 * t134;
t59 = t62 * qJD(3);
t33 = qJD(1) * t59 + t113;
t65 = -t77 * qJDD(1) + qJDD(2);
t10 = t33 * pkin(3) + t32 * pkin(6) + t65;
t7 = t95 * t10;
t122 = qJD(3) * t159;
t131 = t63 * t122 - t159 * t46 + t93 * t45;
t140 = qJD(3) * t93;
t11 = -t64 * t140 - t131;
t8 = qJDD(3) * pkin(6) + t11;
t2 = -qJD(4) * t6 - t92 * t8 + t7;
t171 = t2 + t160;
t170 = t168 * t80;
t120 = t48 * t92;
t41 = t92 * qJD(3) + t95 * t56;
t169 = t41 * t120;
t139 = qJD(4) * t41;
t19 = -t95 * qJDD(3) - t92 * t32 + t139;
t167 = qJ(2) * qJDD(1);
t163 = g(3) * t80;
t102 = -t116 * t81 - t163;
t165 = t56 ^ 2;
t5 = t95 * t20 - t92 * t28;
t161 = t5 * t48;
t135 = t95 * qJD(3);
t138 = qJD(4) * t92;
t18 = -qJD(4) * t135 - t92 * qJDD(3) + t56 * t138 + t95 * t32;
t158 = t18 * t92;
t157 = t19 * t95;
t39 = t92 * t56 - t135;
t156 = t39 * t54;
t155 = t41 * t39;
t154 = t41 * t56;
t153 = t56 * t39;
t152 = t56 * t54;
t151 = t62 * t95;
t15 = t92 * t19;
t26 = qJDD(4) + t33;
t150 = t92 * t26;
t149 = t93 * t64;
t147 = t94 * t92;
t146 = t94 * t95;
t23 = t95 * t26;
t145 = t96 * t92;
t144 = t96 * t95;
t137 = qJD(4) * t95;
t142 = -t39 * t137 - t15;
t86 = t89 ^ 2;
t87 = t90 ^ 2;
t141 = t86 + t87;
t125 = t141 * qJD(1) ^ 2;
t119 = t48 * t95;
t118 = 0.2e1 * t141;
t117 = t81 * pkin(3) + t80 * pkin(6);
t114 = t5 * t95 + t6 * t92;
t1 = qJD(4) * t5 + t92 * t10 + t95 * t8;
t112 = t1 - t161;
t105 = t127 - t148;
t31 = -pkin(3) * t105 - t62 * pkin(6) - t77;
t38 = t159 * t68 - t93 * t67;
t16 = t95 * t31 - t92 * t38;
t17 = t92 * t31 + t95 * t38;
t110 = -t54 * t120 - t48 * t138 + t23;
t109 = -qJD(4) * t20 + t163 - t8;
t58 = -t90 * t122 + t89 * t140;
t108 = t62 * t137 - t58 * t92;
t107 = -t62 * t138 - t58 * t95;
t34 = -t159 * t63 - t149;
t106 = -t159 * t67 - t93 * t68;
t27 = -qJD(3) * pkin(3) - t34;
t104 = -pkin(6) * t26 + t48 * t27;
t103 = t111 + t136;
t99 = t118 * t132 - t116;
t69 = t96 * t77;
t53 = t54 ^ 2;
t52 = t81 * t144 + t147;
t51 = -t81 * t145 + t146;
t50 = -t81 * t146 + t145;
t49 = t81 * t147 + t144;
t30 = t59 * pkin(3) + t58 * pkin(6);
t29 = t56 * pkin(3) + t54 * pkin(6);
t22 = t62 * qJD(2) + t38 * qJD(3);
t21 = t105 * qJD(2) + t106 * qJD(3);
t14 = t92 * t29 + t95 * t34;
t13 = t95 * t29 - t92 * t34;
t4 = -t17 * qJD(4) - t92 * t21 + t95 * t30;
t3 = t16 * qJD(4) + t95 * t21 + t92 * t30;
t24 = [0, 0, 0, 0, 0, qJDD(1), t168, t116, 0, 0, t86 * qJDD(1), 0.2e1 * t89 * t133, 0, t87 * qJDD(1), 0, 0, t103 * t90, -t103 * t89, t118 * t167 + t99, t111 * pkin(1) + (t141 * t167 + t99) * qJ(2), -t32 * t62 - t56 * t58, -t105 * t32 - t62 * t33 + t58 * t54 - t56 * t59, -t58 * qJD(3) + t62 * qJDD(3), -t105 * t33 + t54 * t59, -t59 * qJD(3) + qJDD(3) * t105, 0, -t22 * qJD(3) + qJDD(3) * t106 - t105 * t65 + t168 * t81 - t77 * t33 + t66 * t59, -t21 * qJD(3) - t38 * qJDD(3) + t77 * t32 - t66 * t58 + t65 * t62 - t170, t105 * t11 + t106 * t32 - t12 * t62 - t21 * t54 + t22 * t56 - t38 * t33 + t34 * t58 - t35 * t59 - t116, t11 * t38 + t35 * t21 + t12 * t106 - t34 * t22 - t65 * t77 - g(1) * (t143 * t96 - t94 * t77) - g(2) * (t143 * t94 + t69), t107 * t41 - t18 * t151, (t39 * t95 + t41 * t92) * t58 + (t158 - t157 + (t39 * t92 - t41 * t95) * qJD(4)) * t62, t105 * t18 + t107 * t48 + t62 * t23 + t41 * t59, t108 * t39 + t62 * t15, t105 * t19 - t108 * t48 - t62 * t150 - t39 * t59, -t105 * t26 + t48 * t59, t9 * t92 * t62 - g(1) * t50 - g(2) * t52 - t105 * t2 - t106 * t19 + t108 * t27 + t16 * t26 + t22 * t39 + t4 * t48 + t5 * t59, -g(1) * t49 - g(2) * t51 + t1 * t105 + t106 * t18 + t107 * t27 + t9 * t151 - t17 * t26 + t22 * t41 - t3 * t48 - t6 * t59, t16 * t18 - t17 * t19 - t3 * t39 - t4 * t41 + t170 + t114 * t58 + (-t1 * t92 - t2 * t95 + (t5 * t92 - t6 * t95) * qJD(4)) * t62, -g(2) * t69 + t1 * t17 + t2 * t16 + t27 * t22 + t6 * t3 - t9 * t106 + t5 * t4 + (-g(1) * t143 - g(2) * t117) * t96 + (-g(1) * (-t117 - t77) - g(2) * t143) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, t134, -t125, -qJ(2) * t125 - t111, 0, 0, 0, 0, 0, 0, 0.2e1 * t56 * qJD(3) + t113, (-t54 - t128) * qJD(3) + t130, -t53 - t165, t34 * t56 + t35 * t54 - t168 + t65, 0, 0, 0, 0, 0, 0, t110 - t153, -t48 ^ 2 * t95 - t150 - t154, (t18 - t156) * t95 + t169 + t142, t112 * t92 + t171 * t95 - t27 * t56 - t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, -t53 + t165, (t54 - t128) * qJD(3) + t130, -t152, -t113, qJDD(3), -t66 * t56 + t101 - t123, t66 * t54 + (t34 + t149) * qJD(3) + t131 - t102, 0, 0, t41 * t119 - t158, (-t18 - t156) * t95 - t169 + t142, t48 * t119 + t150 - t154, t39 * t120 - t157, t110 + t153, -t48 * t56, -pkin(3) * t19 + t104 * t92 - t13 * t48 + t172 * t95 - t35 * t39 - t5 * t56, pkin(3) * t18 + t104 * t95 + t14 * t48 - t172 * t92 - t35 * t41 + t6 * t56, t13 * t41 + t14 * t39 + ((-t19 + t139) * pkin(6) + t112) * t95 + ((qJD(4) * t39 - t18) * pkin(6) - t171) * t92 + t102, -t5 * t13 - t6 * t14 - t27 * t35 + t100 * pkin(3) + (-t114 * qJD(4) + t1 * t95 - t2 * t92 + t102) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t39 ^ 2 + t41 ^ 2, t39 * t48 - t18, -t155, t41 * t48 - t19, t26, -g(1) * t51 + g(2) * t49 + t109 * t92 - t28 * t137 - t27 * t41 + t160 + t7, g(1) * t52 - g(2) * t50 + t27 * t39 + t161 + (qJD(4) * t28 - t10) * t92 + t109 * t95, 0, 0;];
tau_reg = t24;

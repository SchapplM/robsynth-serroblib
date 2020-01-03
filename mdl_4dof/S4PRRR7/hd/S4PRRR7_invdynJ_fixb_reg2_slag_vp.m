% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRR7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:51
% EndTime: 2019-12-31 16:36:54
% DurationCPUTime: 1.75s
% Computational Cost: add. (1445->278), mult. (3652->422), div. (0->0), fcn. (2814->10), ass. (0->143)
t81 = sin(qJ(3));
t148 = qJD(1) * t81;
t78 = sin(pkin(4));
t149 = qJD(1) * t78;
t82 = sin(qJ(2));
t53 = qJD(2) * pkin(6) + t82 * t149;
t79 = cos(pkin(4));
t84 = cos(qJ(3));
t33 = t79 * t148 + t84 * t53;
t28 = qJD(3) * pkin(7) + t33;
t85 = cos(qJ(2));
t125 = t85 * t149;
t111 = t84 * pkin(3) + t81 * pkin(7);
t55 = -pkin(2) - t111;
t34 = t55 * qJD(2) - t125;
t80 = sin(qJ(4));
t83 = cos(qJ(4));
t104 = t80 * t28 - t83 * t34;
t131 = qJD(1) * qJD(2);
t122 = t82 * t131;
t161 = t78 * t85;
t107 = -qJDD(1) * t161 + t78 * t122;
t110 = pkin(3) * t81 - pkin(7) * t84;
t52 = t110 * qJD(3);
t14 = qJD(2) * t52 + t55 * qJDD(2) + t107;
t160 = t79 * t84;
t124 = qJD(1) * t160;
t133 = qJDD(1) * t79;
t121 = t85 * t131;
t138 = qJDD(2) * pkin(6);
t36 = t138 + (qJDD(1) * t82 + t121) * t78;
t128 = -qJD(3) * t124 - t81 * t133 - t84 * t36;
t145 = qJD(3) * t81;
t7 = -t53 * t145 - t128;
t5 = qJDD(3) * pkin(7) + t7;
t1 = -t104 * qJD(4) + t80 * t14 + t83 * t5;
t135 = t84 * qJD(2);
t67 = -qJD(4) + t135;
t181 = -t104 * t67 + t1;
t142 = qJD(4) * t81;
t180 = qJD(2) * t142 - qJDD(3);
t132 = t81 * qJDD(2);
t18 = ((qJD(4) + t135) * qJD(3) + t132) * t80 + t180 * t83;
t10 = t83 * t28 + t80 * t34;
t2 = -qJD(4) * t10 + t83 * t14 - t80 * t5;
t179 = t10 * t67 - t2;
t150 = cos(pkin(8));
t114 = t150 * t85;
t77 = sin(pkin(8));
t165 = t77 * t82;
t41 = -t79 * t114 + t165;
t115 = t150 * t82;
t164 = t77 * t85;
t43 = t79 * t164 + t115;
t109 = g(1) * t43 + g(2) * t41;
t139 = qJDD(2) * pkin(2);
t35 = t107 - t139;
t86 = qJD(3) ^ 2;
t178 = -pkin(6) * t86 + (-g(3) * t85 + t122) * t78 + t109 + t139 - t35;
t63 = t84 * t133;
t8 = -t33 * qJD(3) - t81 * t36 + t63;
t130 = qJD(2) * qJD(3);
t120 = t84 * t130;
t136 = t83 * qJD(3);
t17 = -qJD(4) * t136 + (-t120 - t132) * t83 + t180 * t80;
t172 = t17 * t80;
t171 = t18 * t83;
t147 = qJD(2) * t81;
t48 = t80 * t147 - t136;
t170 = t48 * t67;
t169 = t48 * t80;
t137 = t80 * qJD(3);
t50 = t83 * t147 + t137;
t168 = t50 * t48;
t167 = t50 * t67;
t166 = t50 * t83;
t163 = t78 * t82;
t162 = t78 * t84;
t159 = t80 * t84;
t158 = t81 * t53;
t157 = t83 * t84;
t156 = t84 * t85;
t140 = qJD(4) * t84;
t141 = qJD(4) * t83;
t155 = t55 * t141 + t80 * t52 + (-t81 * t136 - t80 * t140) * pkin(6) - (t83 * t156 + t80 * t82) * t149;
t143 = qJD(4) * t80;
t154 = -t55 * t143 + t83 * t52 + (t81 * t137 - t83 * t140) * pkin(6) - (-t80 * t156 + t82 * t83) * t149;
t75 = t81 ^ 2;
t76 = t84 ^ 2;
t153 = t75 - t76;
t152 = t75 + t76;
t151 = qJD(2) * pkin(2);
t146 = qJD(2) * t82;
t144 = qJD(3) * t84;
t134 = qJDD(1) - g(3);
t72 = t84 * qJDD(2);
t87 = qJD(2) ^ 2;
t129 = t81 * t87 * t84;
t127 = t78 * t146;
t126 = qJD(2) * t161;
t123 = g(3) * (pkin(2) * t161 + pkin(6) * t163);
t118 = t81 * t130;
t116 = t78 * t150;
t112 = t84 * t118;
t42 = t79 * t115 + t164;
t44 = -t79 * t165 + t114;
t108 = g(1) * t44 + g(2) * t42;
t106 = -t10 * t80 + t104 * t83;
t103 = qJDD(2) * t85 - t82 * t87;
t46 = t82 * t162 + t79 * t81;
t25 = -t83 * t161 - t46 * t80;
t101 = t80 * t161 - t46 * t83;
t45 = t81 * t163 - t160;
t32 = t124 - t158;
t47 = qJDD(4) - t72 + t118;
t100 = -t67 * t141 + t80 * t47;
t99 = t67 * t143 + t83 * t47;
t97 = -g(1) * (t77 * t162 - t44 * t81) - g(2) * (-t84 * t116 - t42 * t81) + g(3) * t45;
t20 = -t81 * t116 + t42 * t84;
t22 = t77 * t78 * t81 + t44 * t84;
t96 = -g(1) * t22 - g(2) * t20 - g(3) * t46;
t6 = -qJDD(3) * pkin(3) - t8;
t95 = -t6 + t97;
t94 = g(3) * t161 - t109;
t93 = -g(3) * t163 - t108;
t27 = -qJD(3) * pkin(3) - t32;
t92 = -pkin(7) * t47 - t67 * t27;
t90 = pkin(7) * qJD(4) * t67 + t95;
t54 = -t125 - t151;
t89 = -pkin(6) * qJDD(3) + (t125 + t54 - t151) * qJD(3);
t88 = t7 * t84 - t8 * t81 + (-t32 * t84 - t33 * t81) * qJD(3) - t108;
t51 = t110 * qJD(2);
t40 = t43 * pkin(2);
t39 = t41 * pkin(2);
t38 = pkin(6) * t157 + t80 * t55;
t37 = -pkin(6) * t159 + t83 * t55;
t24 = t46 * qJD(3) + t81 * t126;
t23 = -t45 * qJD(3) + t84 * t126;
t16 = t83 * t32 + t80 * t51;
t15 = -t80 * t32 + t83 * t51;
t4 = t25 * qJD(4) + t80 * t127 + t23 * t83;
t3 = t101 * qJD(4) + t83 * t127 - t23 * t80;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t134, 0, 0, 0, 0, 0, 0, t103 * t78, (-qJDD(2) * t82 - t85 * t87) * t78, 0, -g(3) + (t79 ^ 2 + (t82 ^ 2 + t85 ^ 2) * t78 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t24 * qJD(3) - t45 * qJDD(3) + (t103 * t84 - t85 * t118) * t78, -t23 * qJD(3) - t46 * qJDD(3) + (-t103 * t81 - t85 * t120) * t78, (t45 * t81 + t46 * t84) * qJDD(2) + (t23 * t84 + t24 * t81 + (t45 * t84 - t46 * t81) * qJD(3)) * qJD(2), t33 * t23 - t32 * t24 - t8 * t45 + t7 * t46 - g(3) + (t146 * t54 - t35 * t85) * t78, 0, 0, 0, 0, 0, 0, t45 * t18 + t24 * t48 + t25 * t47 - t3 * t67, t101 * t47 - t45 * t17 + t24 * t50 + t4 * t67, t101 * t18 + t25 * t17 - t3 * t50 - t4 * t48, -t1 * t101 + t10 * t4 - t104 * t3 + t2 * t25 + t27 * t24 + t6 * t45 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t134 * t161 + t109, -t134 * t163 + t108, 0, 0, t75 * qJDD(2) + 0.2e1 * t112, -0.2e1 * t153 * t130 + 0.2e1 * t81 * t72, qJDD(3) * t81 + t86 * t84, t76 * qJDD(2) - 0.2e1 * t112, qJDD(3) * t84 - t86 * t81, 0, t178 * t84 + t89 * t81, -t178 * t81 + t89 * t84, t152 * t138 + (-g(3) * t82 - t121 * t152) * t78 + t88, -t35 * pkin(2) + g(1) * t40 + g(2) * t39 - t123 + (-t54 * t82 + (t32 * t81 - t33 * t84) * t85) * t149 + t88 * pkin(6), -t17 * t83 * t81 + (t136 * t84 - t142 * t80) * t50, (-t48 * t83 - t50 * t80) * t144 + (t172 - t171 + (-t166 + t169) * qJD(4)) * t81, (-t136 * t67 + t17) * t84 + (qJD(3) * t50 + t99) * t81, t18 * t80 * t81 + (t137 * t84 + t141 * t81) * t48, (t137 * t67 + t18) * t84 + (-qJD(3) * t48 - t100) * t81, -t145 * t67 - t47 * t84, t37 * t47 - t154 * t67 + t93 * t80 + (-t2 + (pkin(6) * t48 + t27 * t80) * qJD(3) - t94 * t83) * t84 + (pkin(6) * t18 - qJD(3) * t104 - t125 * t48 + t141 * t27 + t6 * t80) * t81, -t38 * t47 + t155 * t67 + t93 * t83 + (t1 + (pkin(6) * t50 + t27 * t83) * qJD(3) + t94 * t80) * t84 + (-pkin(6) * t17 - t10 * qJD(3) - t125 * t50 - t143 * t27 + t6 * t83) * t81, t37 * t17 - t38 * t18 - t154 * t50 - t155 * t48 + t106 * t144 + (-t1 * t80 - t2 * t83 + (-t10 * t83 - t104 * t80) * qJD(4) - t94) * t81, t1 * t38 + t2 * t37 - g(1) * (-t111 * t43 - t40) - g(2) * (-t111 * t41 - t39) - t123 - t154 * t104 + (-g(3) * t111 - t148 * t27) * t161 + t155 * t10 + (t144 * t27 + t6 * t81 - t108) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, t153 * t87, t132, t129, t72, qJDD(3), t63 + (-qJD(2) * t54 - t36) * t81 + t97, -t54 * t135 + (t32 + t158) * qJD(3) - t96 + t128, 0, 0, -t166 * t67 - t172, (-t17 + t170) * t83 + (-t18 + t167) * t80, (t157 * t67 - t50 * t81) * qJD(2) + t100, -t169 * t67 - t171, (-t159 * t67 + t48 * t81) * qJD(2) + t99, t67 * t147, -pkin(3) * t18 + t104 * t147 + t15 * t67 - t33 * t48 + t80 * t92 + t83 * t90, pkin(3) * t17 + t10 * t147 - t16 * t67 - t33 * t50 - t80 * t90 + t83 * t92, t15 * t50 + t16 * t48 + ((qJD(4) * t50 - t18) * pkin(7) + t181) * t83 + ((qJD(4) * t48 - t17) * pkin(7) + t179) * t80 + t96, -t10 * t16 + t104 * t15 - t27 * t33 + t95 * pkin(3) + (qJD(4) * t106 + t1 * t83 - t2 * t80 + t96) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, -t48 ^ 2 + t50 ^ 2, -t17 - t170, -t168, -t167 - t18, t47, -t27 * t50 - g(1) * (-t22 * t80 + t43 * t83) - g(2) * (-t20 * t80 + t41 * t83) - g(3) * t25 - t179, t27 * t48 - g(1) * (-t22 * t83 - t43 * t80) - g(2) * (-t20 * t83 - t41 * t80) - g(3) * t101 - t181, 0, 0;];
tau_reg = t9;

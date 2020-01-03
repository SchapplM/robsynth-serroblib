% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRR4
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:38
% EndTime: 2019-12-31 16:50:40
% DurationCPUTime: 1.14s
% Computational Cost: add. (1288->224), mult. (2762->315), div. (0->0), fcn. (1677->10), ass. (0->128)
t67 = sin(pkin(7));
t51 = t67 * pkin(1) + pkin(5);
t41 = t51 * qJDD(1);
t154 = -qJD(2) * qJD(3) - t41;
t70 = sin(qJ(3));
t125 = qJD(4) * t70;
t153 = qJD(1) * t125 - qJDD(3);
t120 = t70 * qJDD(1);
t73 = cos(qJ(3));
t121 = t73 * qJD(1);
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t16 = ((qJD(4) + t121) * qJD(3) + t120) * t69 + t153 * t72;
t47 = -qJD(4) + t121;
t147 = g(3) * t73;
t43 = t51 * qJD(1);
t134 = t73 * t43;
t29 = t70 * qJD(2) + t134;
t60 = t73 * qJDD(2);
t11 = -t29 * qJD(3) - t70 * t41 + t60;
t9 = -qJDD(3) * pkin(3) - t11;
t64 = qJ(1) + pkin(7);
t55 = sin(t64);
t56 = cos(t64);
t95 = g(1) * t56 + g(2) * t55;
t81 = t95 * t70 - t147 - t9;
t152 = pkin(6) * qJD(4) * t47 + t81;
t132 = pkin(1) * qJDD(1);
t119 = qJD(1) * qJD(3);
t104 = t70 * t119;
t61 = t73 * qJDD(1);
t32 = qJDD(4) - t61 + t104;
t136 = t72 * t32;
t122 = t72 * qJD(3);
t109 = t73 * t122;
t86 = -t69 * t125 + t109;
t151 = t70 * t136 - t47 * t86;
t149 = g(1) * t55;
t148 = g(3) * t70;
t21 = qJD(3) * pkin(6) + t29;
t68 = cos(pkin(7));
t144 = t68 * pkin(1);
t97 = t73 * pkin(3) + t70 * pkin(6);
t90 = -pkin(2) - t97;
t31 = t90 - t144;
t22 = t31 * qJD(1);
t5 = -t69 * t21 + t72 * t22;
t146 = t5 * t47;
t6 = t72 * t21 + t69 * t22;
t145 = t6 * t47;
t137 = t70 * t72;
t131 = qJD(1) * t70;
t34 = t69 * t131 - t122;
t143 = -t34 * t109 - t16 * t137;
t142 = t34 * t47;
t123 = t69 * qJD(3);
t36 = t72 * t131 + t123;
t141 = t36 * t34;
t140 = t36 * t47;
t139 = t69 * t73;
t138 = t70 * t43;
t135 = t72 * t73;
t65 = t70 ^ 2;
t66 = t73 ^ 2;
t133 = t65 - t66;
t52 = -pkin(2) - t144;
t44 = qJD(1) * t52;
t130 = qJD(3) * t34;
t129 = qJD(3) * t70;
t128 = qJD(3) * t73;
t127 = qJD(4) * t34;
t126 = qJD(4) * t69;
t124 = qJD(4) * t72;
t42 = qJDD(1) * t52;
t76 = qJD(1) ^ 2;
t117 = t70 * t76 * t73;
t115 = -t70 * qJDD(2) + t154 * t73;
t74 = cos(qJ(1));
t114 = t74 * pkin(1) + t56 * pkin(2) + t55 * pkin(5);
t112 = t36 * t128;
t111 = t47 * t123;
t110 = t51 * t129;
t108 = t70 * t124;
t71 = sin(qJ(1));
t106 = -t71 * pkin(1) + t56 * pkin(5);
t15 = -qJD(1) * t109 - qJD(4) * t122 - t72 * t120 + t153 * t69;
t103 = t36 * t129 + t15 * t73;
t101 = -t15 + t127;
t100 = t36 * t108;
t98 = t73 * t104;
t96 = pkin(3) * t70 - pkin(6) * t73;
t94 = -g(2) * t56 + t149;
t93 = g(1) * t71 - g(2) * t74;
t92 = -t5 * t72 - t6 * t69;
t91 = t5 * t69 - t6 * t72;
t28 = t73 * qJD(2) - t138;
t19 = t51 * t135 + t69 * t31;
t18 = -t51 * t139 + t72 * t31;
t10 = -t43 * t129 - t115;
t8 = qJDD(3) * pkin(6) + t10;
t88 = -qJD(4) * t22 + t148 - t8;
t87 = t47 * t124 - t69 * t32;
t38 = t96 * qJD(3);
t85 = -qJD(1) * t44 + t95;
t20 = -qJD(3) * pkin(3) - t28;
t84 = -pkin(6) * t32 - t47 * t20;
t83 = 0.2e1 * t44 * qJD(3) - qJDD(3) * t51;
t82 = -t95 * t73 - t148;
t75 = qJD(3) ^ 2;
t80 = -t51 * t75 - 0.2e1 * t42 + t94;
t17 = qJD(1) * t38 + t31 * qJDD(1);
t1 = qJD(4) * t5 + t69 * t17 + t72 * t8;
t12 = t72 * t17;
t2 = -qJD(4) * t6 - t69 * t8 + t12;
t79 = t92 * qJD(4) + t1 * t72 - t2 * t69;
t78 = t10 * t73 - t11 * t70 + (-t28 * t73 - t29 * t70) * qJD(3);
t40 = qJDD(3) * t73 - t75 * t70;
t39 = qJDD(3) * t70 + t75 * t73;
t37 = t96 * qJD(1);
t27 = t56 * t135 + t55 * t69;
t26 = -t56 * t139 + t55 * t72;
t25 = -t55 * t135 + t56 * t69;
t24 = t55 * t139 + t56 * t72;
t14 = t72 * t28 + t69 * t37;
t13 = -t69 * t28 + t72 * t37;
t4 = -t19 * qJD(4) + t69 * t110 + t72 * t38;
t3 = t18 * qJD(4) - t72 * t110 + t69 * t38;
t7 = [0, 0, 0, 0, 0, qJDD(1), t93, g(1) * t74 + g(2) * t71, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t68 * t132 + t94, -0.2e1 * t67 * t132 + t95, 0, (t93 + (t67 ^ 2 + t68 ^ 2) * t132) * pkin(1), t65 * qJDD(1) + 0.2e1 * t98, -0.2e1 * t133 * t119 + 0.2e1 * t70 * t61, t39, t66 * qJDD(1) - 0.2e1 * t98, t40, 0, t83 * t70 + t80 * t73, -t80 * t70 + t83 * t73, (t65 + t66) * t41 + t78 - t95, t42 * t52 - g(1) * (-t55 * pkin(2) + t106) - g(2) * t114 + t78 * t51, -t15 * t137 + t86 * t36, -t100 + (-t112 + (t15 + t127) * t70) * t69 + t143, t103 + t151, t16 * t69 * t70 + (t123 * t73 + t108) * t34, (t16 + t111) * t73 + (t87 - t130) * t70, -t129 * t47 - t32 * t73, -g(1) * t25 - g(2) * t27 + t18 * t32 - t4 * t47 + (-t2 + (t20 * t69 + t34 * t51) * qJD(3)) * t73 + (qJD(3) * t5 + t124 * t20 + t16 * t51 + t9 * t69) * t70, -g(1) * t24 - g(2) * t26 - t19 * t32 + t3 * t47 + (t1 + (t20 * t72 + t36 * t51) * qJD(3)) * t73 + (-qJD(3) * t6 - t126 * t20 - t15 * t51 + t9 * t72) * t70, t18 * t15 - t19 * t16 - t3 * t34 - t4 * t36 + t92 * t128 + (qJD(4) * t91 - t1 * t69 - t2 * t72 + t94) * t70, t1 * t19 + t6 * t3 + t2 * t18 + t5 * t4 - g(1) * t106 - g(2) * (t56 * t97 + t114) - t90 * t149 + (t128 * t20 + t9 * t70) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t40, -t39, 0, t10 * t70 + t11 * t73 - g(3) + (-t28 * t70 + t29 * t73) * qJD(3), 0, 0, 0, 0, 0, 0, (-t16 + t111) * t73 + (t87 + t130) * t70, t103 - t151, t100 + (t101 * t70 + t112) * t69 + t143, -g(3) + (-qJD(3) * t91 - t9) * t73 + (qJD(3) * t20 + t79) * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, t133 * t76, t120, t117, t61, qJDD(3), -t147 + t60 + (t29 - t134) * qJD(3) + (t85 + t154) * t70, t148 + (t28 + t138) * qJD(3) + t85 * t73 + t115, 0, 0, -t72 * t140 - t15 * t69, (-t15 + t142) * t72 + (-t16 + t140) * t69, (t47 * t135 - t36 * t70) * qJD(1) - t87, -t69 * t142 - t16 * t72, t47 * t126 + t136 + (-t47 * t139 + t34 * t70) * qJD(1), t47 * t131, -pkin(3) * t16 + t13 * t47 - t5 * t131 + t152 * t72 - t29 * t34 + t84 * t69, pkin(3) * t15 + t6 * t131 - t14 * t47 - t152 * t69 - t29 * t36 + t84 * t72, t13 * t36 + t14 * t34 + (t1 + t146 + (qJD(4) * t36 - t16) * pkin(6)) * t72 + (pkin(6) * t101 + t145 - t2) * t69 + t82, -t5 * t13 - t6 * t14 - t20 * t29 + t81 * pkin(3) + (t79 + t82) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, -t34 ^ 2 + t36 ^ 2, -t15 - t142, -t141, -t140 - t16, t32, -g(1) * t26 + g(2) * t24 - t124 * t21 - t20 * t36 + t69 * t88 + t12 - t145, g(1) * t27 - g(2) * t25 + t20 * t34 - t146 + (qJD(4) * t21 - t17) * t69 + t88 * t72, 0, 0;];
tau_reg = t7;

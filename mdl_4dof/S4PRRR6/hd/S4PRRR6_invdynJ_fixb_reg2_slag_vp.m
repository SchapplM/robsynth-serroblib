% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRR6
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:09
% EndTime: 2019-12-31 16:35:11
% DurationCPUTime: 1.13s
% Computational Cost: add. (1062->199), mult. (2480->285), div. (0->0), fcn. (1718->10), ass. (0->121)
t76 = cos(qJ(3));
t120 = t76 * qJDD(2);
t73 = sin(qJ(3));
t121 = t73 * qJDD(2);
t72 = sin(qJ(4));
t75 = cos(qJ(4));
t43 = t72 * t76 + t75 * t73;
t66 = qJD(3) + qJD(4);
t156 = t66 * t43;
t11 = qJD(2) * t156 - t75 * t120 + t72 * t121;
t74 = sin(qJ(2));
t149 = g(3) * t74;
t70 = sin(pkin(7));
t71 = cos(pkin(7));
t102 = g(1) * t71 + g(2) * t70;
t77 = cos(qJ(2));
t96 = t102 * t77;
t87 = t96 + t149;
t148 = g(3) * t77;
t97 = t102 * t74;
t86 = t97 - t148;
t117 = qJD(2) * qJD(3);
t109 = t76 * t117;
t155 = -t109 - t121;
t67 = t73 ^ 2;
t68 = t76 ^ 2;
t136 = t67 + t68;
t111 = t136 * t77;
t125 = t74 * qJD(1);
t52 = qJD(2) * pkin(5) + t125;
t124 = t77 * qJD(1);
t134 = qJD(2) * pkin(2);
t53 = -t124 - t134;
t154 = t52 * t111 + t53 * t74;
t138 = t75 * t76;
t140 = t72 * t73;
t42 = -t138 + t140;
t35 = t42 * t74;
t118 = qJD(1) * qJD(2);
t127 = qJDD(2) * pkin(2);
t39 = -t77 * qJDD(1) + t74 * t118 - t127;
t79 = qJD(3) ^ 2;
t153 = -pkin(5) * t79 + (t102 + t118) * t74 + t127 - t148 - t39;
t129 = qJD(4) * t72;
t130 = qJD(3) * t76;
t40 = qJDD(2) * pkin(5) + t74 * qJDD(1) + t77 * t118;
t13 = qJDD(3) * pkin(3) + t155 * pkin(6) - t52 * t130 - t73 * t40;
t131 = qJD(3) * t73;
t110 = t73 * t117;
t91 = t110 - t120;
t16 = -t91 * pkin(6) - t52 * t131 + t76 * t40;
t107 = pkin(6) * qJD(2) + t52;
t32 = t107 * t73;
t25 = qJD(3) * pkin(3) - t32;
t33 = t107 * t76;
t1 = (qJD(4) * t25 + t16) * t75 - t33 * t129 + t72 * t13;
t152 = pkin(6) + pkin(5);
t46 = t152 * t73;
t47 = t152 * t76;
t22 = -t75 * t46 - t72 * t47;
t113 = qJD(3) * t152;
t44 = t73 * t113;
t45 = t76 * t113;
t92 = t42 * t77;
t147 = qJD(1) * t92 + t22 * qJD(4) - t75 * t44 - t72 * t45;
t23 = -t72 * t46 + t75 * t47;
t93 = t43 * t77;
t146 = qJD(1) * t93 - t23 * qJD(4) + t72 * t44 - t75 * t45;
t114 = qJD(2) * t138;
t133 = qJD(2) * t73;
t115 = t72 * t133;
t36 = -t114 + t115;
t38 = t43 * qJD(2);
t145 = t38 * t36;
t143 = t70 * t77;
t142 = t71 * t77;
t141 = t72 * t33;
t139 = t75 * t33;
t137 = t67 - t68;
t80 = qJD(2) ^ 2;
t135 = t79 + t80;
t132 = qJD(2) * t74;
t128 = qJD(4) * t75;
t62 = t76 * pkin(3) + pkin(2);
t41 = -t62 * qJD(2) - t124;
t126 = t41 * qJD(2);
t123 = qJDD(1) - g(3);
t122 = qJDD(3) * t73;
t119 = t77 * qJDD(2);
t116 = t73 * t80 * t76;
t112 = t136 * t40;
t106 = -qJD(4) * t114 - t72 * t120 + t155 * t75;
t104 = t136 * qJDD(2);
t103 = t73 * t109;
t101 = g(1) * t70 - g(2) * t71;
t99 = t66 * t140;
t15 = t72 * t25 + t139;
t95 = t101 * t76;
t94 = pkin(3) * t131 - t125;
t2 = -t15 * qJD(4) + t75 * t13 - t72 * t16;
t85 = -pkin(5) * qJDD(3) + (t124 + t53 - t134) * qJD(3);
t84 = -t53 * qJD(2) - t40 + t87;
t69 = qJ(3) + qJ(4);
t63 = sin(t69);
t64 = cos(t69);
t82 = -g(1) * (-t64 * t142 - t70 * t63) - g(2) * (-t64 * t143 + t71 * t63) + t41 * t36 + t64 * t149 - t1;
t81 = -g(1) * (-t63 * t142 + t70 * t64) - g(2) * (-t63 * t143 - t71 * t64) - t41 * t38 + t2 + t63 * t149;
t65 = qJDD(3) + qJDD(4);
t34 = t43 * t74;
t21 = t91 * pkin(3) + t39;
t19 = -t76 * t128 - t75 * t130 + t99;
t18 = -t75 * t32 - t141;
t17 = t72 * t32 - t139;
t14 = t75 * t25 - t141;
t12 = -t36 ^ 2 + t38 ^ 2;
t10 = qJD(2) * t99 + t106;
t6 = -qJD(2) * t93 + t66 * t35;
t5 = -qJD(2) * t92 - t156 * t74;
t4 = t38 * t66 - t11;
t3 = -t106 + (-t115 + t36) * t66;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t123, 0, 0, 0, 0, 0, 0, -t80 * t74 + t119, -qJDD(2) * t74 - t80 * t77, 0, -g(3) + (t74 ^ 2 + t77 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, (-0.2e1 * t110 + t120) * t77 + (-t135 * t76 - t122) * t74, (-qJDD(3) * t74 - 0.2e1 * t77 * t117) * t76 + (t135 * t74 - t119) * t73, t74 * t104 + t80 * t111, t154 * qJD(2) + t74 * t112 - t39 * t77 - g(3), 0, 0, 0, 0, 0, 0, -t77 * t11 + t36 * t132 - t34 * t65 + t6 * t66, t77 * t10 + t38 * t132 + t35 * t65 - t5 * t66, -t34 * t10 + t35 * t11 - t5 * t36 - t6 * t38, -t1 * t35 + t74 * t126 + t14 * t6 + t15 * t5 - t2 * t34 - t21 * t77 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t123 * t77 + t97, -t123 * t74 + t96, 0, 0, t67 * qJDD(2) + 0.2e1 * t103, -0.2e1 * t137 * t117 + 0.2e1 * t73 * t120, t79 * t76 + t122, t68 * qJDD(2) - 0.2e1 * t103, qJDD(3) * t76 - t79 * t73, 0, t153 * t76 + t85 * t73, -t153 * t73 + t85 * t76, -t149 + t112 + pkin(5) * t104 + (-t136 * t118 - t102) * t77, (-t39 + t86) * pkin(2) + (t112 - t87) * pkin(5) - t154 * qJD(1), -t10 * t43 - t38 * t19, t10 * t42 - t43 * t11 - t156 * t38 + t19 * t36, -t19 * t66 + t43 * t65, t11 * t42 + t156 * t36, -t156 * t66 - t42 * t65, 0, -t62 * t11 + t146 * t66 + t156 * t41 + t21 * t42 + t22 * t65 + t94 * t36 + t86 * t64, t62 * t10 - t147 * t66 - t41 * t19 + t21 * t43 - t23 * t65 + t94 * t38 - t86 * t63, -t1 * t42 + t22 * t10 - t23 * t11 + t14 * t19 - t146 * t38 - t147 * t36 - t15 * t156 - t2 * t43 - t87, t1 * t23 + t2 * t22 - t21 * t62 - g(3) * (t152 * t74 + t77 * t62) + t94 * t41 + t147 * t15 + t146 * t14 + t102 * (-t152 * t77 + t62 * t74); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, t137 * t80, t121, t116, t120, qJDD(3), t84 * t73 - t95, t101 * t73 + t84 * t76, 0, 0, t145, t12, t3, -t145, t4, t65, -t17 * t66 + (-t66 * t129 - t36 * t133 + t65 * t75) * pkin(3) + t81, t18 * t66 + (-t66 * t128 - t38 * t133 - t65 * t72) * pkin(3) + t82, (t15 + t17) * t38 + (-t14 + t18) * t36 + (t10 * t75 - t11 * t72 + (-t36 * t75 + t38 * t72) * qJD(4)) * pkin(3), -t14 * t17 - t15 * t18 + (t1 * t72 + t2 * t75 - t95 + (t87 - t126) * t73 + (-t14 * t72 + t15 * t75) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, t12, t3, -t145, t4, t65, t15 * t66 + t81, t14 * t66 + t82, 0, 0;];
tau_reg = t7;

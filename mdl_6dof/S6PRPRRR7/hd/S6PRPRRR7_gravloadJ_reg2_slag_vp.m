% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_gravloadJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:55:40
% EndTime: 2019-05-05 01:55:43
% DurationCPUTime: 1.07s
% Computational Cost: add. (1775->183), mult. (5156->304), div. (0->0), fcn. (6752->18), ass. (0->124)
t142 = cos(pkin(8));
t72 = sin(pkin(8));
t139 = sin(pkin(6));
t140 = cos(pkin(14));
t118 = t140 * t139;
t73 = sin(pkin(7));
t113 = t73 * t118;
t137 = sin(pkin(14));
t138 = sin(pkin(13));
t143 = cos(pkin(7));
t144 = cos(pkin(6));
t122 = t144 * t138;
t141 = cos(pkin(13));
t151 = sin(qJ(2));
t153 = cos(qJ(2));
t66 = -t151 * t122 + t141 * t153;
t95 = t153 * t122 + t141 * t151;
t93 = t95 * t140;
t81 = -t138 * t113 + t66 * t137 + t143 * t93;
t120 = t143 * t139;
t88 = -t138 * t120 - t95 * t73;
t159 = t81 * t142 + t88 * t72;
t123 = t144 * t141;
t65 = t151 * t123 + t138 * t153;
t96 = -t153 * t123 + t138 * t151;
t91 = t96 * t140;
t82 = t141 * t113 + t65 * t137 + t143 * t91;
t87 = t141 * t120 - t96 * t73;
t158 = t82 * t142 + t87 * t72;
t103 = t143 * t118;
t117 = t139 * t137;
t133 = t73 * t144;
t86 = -t153 * t103 + t151 * t117 - t140 * t133;
t125 = t153 * t139;
t97 = t73 * t125 - t144 * t143;
t157 = t86 * t142 + t97 * t72;
t132 = t73 * t142;
t121 = t143 * t140;
t90 = t96 * t137;
t48 = -t65 * t121 + t90;
t35 = t65 * t132 - t48 * t72;
t92 = t95 * t137;
t50 = -t66 * t121 + t92;
t36 = t66 * t132 - t50 * t72;
t156 = -g(1) * t36 - g(2) * t35;
t152 = cos(qJ(4));
t61 = -t151 * t103 - t153 * t117;
t150 = t61 * t72;
t149 = t72 * t73;
t74 = sin(qJ(6));
t78 = cos(qJ(5));
t148 = t74 * t78;
t77 = cos(qJ(6));
t147 = t77 * t78;
t124 = t139 * t151;
t114 = t73 * t124;
t146 = pkin(2) * t125 + qJ(3) * t114;
t145 = qJ(3) * t73;
t112 = t73 * t117;
t42 = -t141 * t112 + t65 * t140 - t143 * t90;
t76 = sin(qJ(4));
t14 = t158 * t152 + t42 * t76;
t15 = t42 * t152 - t158 * t76;
t136 = -t14 * pkin(4) + t15 * pkin(11);
t43 = t138 * t112 + t66 * t140 - t143 * t92;
t16 = t159 * t152 + t43 * t76;
t17 = t43 * t152 - t159 * t76;
t135 = -t16 * pkin(4) + t17 * pkin(11);
t102 = t143 * t117;
t57 = t153 * t102 + t151 * t118 + t137 * t133;
t27 = t157 * t152 + t57 * t76;
t28 = t57 * t152 - t157 * t76;
t134 = -t27 * pkin(4) + t28 * pkin(11);
t101 = t142 * t114;
t62 = -t151 * t102 + t153 * t118;
t131 = t62 * pkin(3) + pkin(10) * t101 + t146;
t130 = t152 * t149;
t129 = -t96 * pkin(2) + t65 * t145;
t128 = -t95 * pkin(2) + t66 * t145;
t75 = sin(qJ(5));
t127 = -pkin(5) * t78 - pkin(12) * t75;
t126 = t142 * t152;
t119 = t143 * t137;
t49 = -t65 * t119 - t91;
t116 = t49 * pkin(3) + t129;
t51 = -t66 * t119 - t93;
t115 = t51 * pkin(3) + t128;
t41 = -t97 * t142 + t86 * t72;
t10 = -t28 * t75 + t41 * t78;
t29 = -t87 * t142 + t82 * t72;
t2 = -t15 * t75 + t29 * t78;
t30 = -t88 * t142 + t81 * t72;
t4 = -t17 * t75 + t30 * t78;
t111 = g(1) * t4 + g(2) * t2 + g(3) * t10;
t11 = t28 * t78 + t41 * t75;
t3 = t15 * t78 + t29 * t75;
t5 = t17 * t78 + t30 * t75;
t110 = g(1) * t5 + g(2) * t3 + g(3) * t11;
t108 = t72 * t114;
t38 = t62 * t152 + (t142 * t61 + t108) * t76;
t53 = t101 - t150;
t24 = t38 * t75 - t53 * t78;
t21 = t49 * t152 + (t142 * t48 + t65 * t149) * t76;
t6 = t21 * t75 - t35 * t78;
t23 = t51 * t152 + (t142 * t50 + t66 * t149) * t76;
t8 = t23 * t75 - t36 * t78;
t109 = g(1) * t8 + g(2) * t6 + g(3) * t24;
t107 = g(1) * t16 + g(2) * t14 + g(3) * t27;
t106 = g(1) * t17 + g(2) * t15 + g(3) * t28;
t20 = -t48 * t126 - t65 * t130 + t49 * t76;
t22 = -t50 * t126 - t66 * t130 + t51 * t76;
t37 = -t152 * t108 - t61 * t126 + t62 * t76;
t105 = g(1) * t22 + g(2) * t20 + g(3) * t37;
t104 = t38 * pkin(4) + pkin(11) * t37 + t131;
t100 = t21 * pkin(4) + t20 * pkin(11) + t116;
t99 = t23 * pkin(4) + t22 * pkin(11) + t115;
t98 = -g(1) * t66 - g(2) * t65 - g(3) * t124;
t89 = (g(3) * t150 + t156) * pkin(10);
t31 = g(1) * t88 + g(2) * t87 + g(3) * t97;
t25 = t38 * t78 + t53 * t75;
t9 = t23 * t78 + t36 * t75;
t7 = t21 * t78 + t35 * t75;
t1 = t107 * t75;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t95 + g(2) * t96 - g(3) * t125, -t98, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t51 - g(2) * t49 - g(3) * t62, -g(1) * t50 - g(2) * t48 - g(3) * t61, t98 * t73, -g(1) * t128 - g(2) * t129 - g(3) * t146, 0, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t21 - g(3) * t38, t105, -g(3) * t53 + t156, -g(1) * t115 - g(2) * t116 - g(3) * t131 + t89, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7 - g(3) * t25, t109, -t105, -g(1) * t99 - g(2) * t100 - g(3) * t104 + t89, 0, 0, 0, 0, 0, 0, -g(1) * (t22 * t74 + t77 * t9) - g(2) * (t20 * t74 + t7 * t77) - g(3) * (t25 * t77 + t37 * t74) -g(1) * (t22 * t77 - t74 * t9) - g(2) * (t20 * t77 - t7 * t74) - g(3) * (-t25 * t74 + t37 * t77) -t109, -g(1) * (t9 * pkin(5) + t8 * pkin(12) + t99) - g(2) * (t7 * pkin(5) + t6 * pkin(12) + t100) - g(3) * (pkin(5) * t25 + pkin(12) * t24 + t104) + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t106, 0, 0, 0, 0, 0, 0, 0, 0, t107 * t78, -t1, -t106, -g(1) * t135 - g(2) * t136 - g(3) * t134, 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t147 + t17 * t74) - g(2) * (-t14 * t147 + t15 * t74) - g(3) * (-t27 * t147 + t28 * t74) -g(1) * (t16 * t148 + t17 * t77) - g(2) * (t14 * t148 + t15 * t77) - g(3) * (t27 * t148 + t28 * t77) t1, -g(1) * (t127 * t16 + t135) - g(2) * (t127 * t14 + t136) - g(3) * (t127 * t27 + t134); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, t110, 0, 0, 0, 0, 0, 0, 0, 0, -t111 * t77, t111 * t74, -t110, -g(1) * (pkin(5) * t4 + pkin(12) * t5) - g(2) * (pkin(5) * t2 + pkin(12) * t3) - g(3) * (pkin(5) * t10 + pkin(12) * t11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t16 * t77 - t5 * t74) - g(2) * (t14 * t77 - t3 * t74) - g(3) * (-t11 * t74 + t27 * t77) -g(1) * (-t16 * t74 - t5 * t77) - g(2) * (-t14 * t74 - t3 * t77) - g(3) * (-t11 * t77 - t27 * t74) 0, 0;];
taug_reg  = t12;

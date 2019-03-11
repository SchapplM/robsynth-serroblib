% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t116 = sin(pkin(7));
t69 = cos(qJ(1));
t118 = cos(pkin(12));
t120 = cos(pkin(6));
t101 = t120 * t118;
t115 = sin(pkin(12));
t128 = sin(qJ(1));
t88 = -t69 * t101 + t128 * t115;
t117 = sin(pkin(6));
t119 = cos(pkin(7));
t99 = t117 * t119;
t135 = t88 * t116 - t69 * t99;
t129 = cos(qJ(3));
t98 = t117 * t116;
t143 = t88 * t119 + t69 * t98;
t100 = t120 * t115;
t52 = t69 * t100 + t128 * t118;
t66 = sin(qJ(3));
t34 = -t52 * t129 + t143 * t66;
t65 = sin(qJ(4));
t68 = cos(qJ(4));
t16 = t135 * t68 + t34 * t65;
t31 = t143 * t129 + t52 * t66;
t64 = sin(qJ(6));
t67 = cos(qJ(6));
t147 = t16 * t64 - t31 * t67;
t146 = t16 * t67 + t31 * t64;
t17 = -t135 * t65 + t34 * t68;
t82 = t128 * t101 + t69 * t115;
t142 = t82 * t116 + t128 * t99;
t121 = qJ(5) * t65;
t127 = t31 * t68;
t139 = -pkin(4) * t127 - t31 * t121;
t136 = t82 * t119 - t128 * t98;
t53 = -t128 * t100 + t69 * t118;
t35 = t129 * t136 + t53 * t66;
t126 = t35 * t68;
t138 = -pkin(4) * t126 - t35 * t121;
t134 = t116 * t120 + t118 * t99;
t97 = t117 * t115;
t43 = -t129 * t134 + t66 * t97;
t125 = t43 * t68;
t137 = -pkin(4) * t125 - t43 * t121;
t132 = pkin(5) + pkin(10);
t131 = t31 * pkin(10);
t130 = t35 * pkin(10);
t124 = t64 * t65;
t123 = t65 * t67;
t102 = t128 * t117;
t122 = t69 * pkin(1) + qJ(2) * t102;
t25 = t31 * pkin(3);
t114 = -pkin(10) * t34 - t25;
t27 = t35 * pkin(3);
t36 = t53 * t129 - t136 * t66;
t113 = t36 * pkin(10) - t27;
t42 = t43 * pkin(3);
t44 = t129 * t97 + t134 * t66;
t112 = t44 * pkin(10) - t42;
t111 = t69 * t117;
t110 = pkin(4) * t16 - qJ(5) * t17;
t18 = -t142 * t68 + t36 * t65;
t19 = t142 * t65 + t36 * t68;
t109 = -t18 * pkin(4) + t19 * qJ(5);
t81 = -t118 * t98 + t120 * t119;
t29 = t44 * t65 - t81 * t68;
t30 = t44 * t68 + t81 * t65;
t108 = -t29 * pkin(4) + t30 * qJ(5);
t107 = -t128 * pkin(1) + qJ(2) * t111;
t106 = g(1) * t16 + g(2) * t18;
t105 = g(1) * t17 + g(2) * t19;
t104 = -g(1) * t31 + g(2) * t35;
t2 = g(1) * t18 - g(2) * t16 + g(3) * t29;
t93 = g(1) * t19 - g(2) * t17 + g(3) * t30;
t92 = g(1) * t35 + g(2) * t31 + g(3) * t43;
t91 = g(1) * t36 - g(2) * t34 + g(3) * t44;
t76 = -t52 * pkin(2) - t135 * pkin(9) + t107;
t75 = t34 * pkin(3) + t76;
t74 = t53 * pkin(2) + t142 * pkin(9) + t122;
t73 = t36 * pkin(3) + t74;
t72 = t17 * pkin(4) + t16 * qJ(5) + t75;
t70 = t19 * pkin(4) + t18 * qJ(5) + t73;
t49 = -g(1) * t102 + g(2) * t111 - g(3) * t120;
t7 = t18 * t64 + t35 * t67;
t6 = t18 * t67 - t35 * t64;
t5 = t92 * t68;
t4 = t92 * t65;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t128 - g(2) * t69, g(1) * t69 + g(2) * t128, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t52 - g(2) * t53, -g(1) * t88 + g(2) * t82, -g(1) * t111 - g(2) * t102, -g(1) * t107 - g(2) * t122, 0, 0, 0, 0, 0, 0, -g(1) * t34 - g(2) * t36, t104, g(1) * t135 - g(2) * t142, -g(1) * t76 - g(2) * t74, 0, 0, 0, 0, 0, 0, -t105, t106, -t104, -g(1) * (t75 - t131) - g(2) * (t73 + t130) 0, 0, 0, 0, 0, 0, -t104, t105, -t106, -g(1) * (t72 - t131) - g(2) * (t70 + t130) 0, 0, 0, 0, 0, 0, -g(1) * t147 - g(2) * t7, -g(1) * t146 - g(2) * t6, -t105, -g(1) * (t17 * pkin(11) - t132 * t31 + t72) - g(2) * (t19 * pkin(11) + t132 * t35 + t70); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t91, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, -t91, -g(1) * t113 - g(2) * t114 - g(3) * t112, 0, 0, 0, 0, 0, 0, -t91, -t5, t4, -g(1) * (t113 + t138) - g(2) * (t114 + t139) - g(3) * (t112 + t137) 0, 0, 0, 0, 0, 0, -g(1) * (-t35 * t124 + t36 * t67) - g(2) * (-t31 * t124 - t34 * t67) - g(3) * (-t43 * t124 + t44 * t67) -g(1) * (-t35 * t123 - t36 * t64) - g(2) * (-t31 * t123 + t34 * t64) - g(3) * (-t43 * t123 - t44 * t64) t5, -g(1) * (-pkin(11) * t126 + t132 * t36 + t138 - t27) - g(2) * (-pkin(11) * t127 - t132 * t34 + t139 - t25) - g(3) * (-pkin(11) * t125 + t132 * t44 + t137 - t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t93, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t93, -g(1) * t109 - g(2) * t110 - g(3) * t108, 0, 0, 0, 0, 0, 0, -t93 * t64, -t93 * t67, t2, -g(1) * (-t18 * pkin(11) + t109) - g(2) * (pkin(11) * t16 + t110) - g(3) * (-t29 * pkin(11) + t108); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 + g(2) * t146 - g(3) * (t29 * t67 - t43 * t64) g(1) * t7 - g(2) * t147 - g(3) * (-t29 * t64 - t43 * t67) 0, 0;];
taug_reg  = t1;

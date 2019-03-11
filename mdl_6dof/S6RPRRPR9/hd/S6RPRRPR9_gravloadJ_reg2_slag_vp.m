% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t133 = cos(qJ(3));
t116 = cos(pkin(7));
t74 = cos(qJ(1));
t115 = cos(pkin(12));
t117 = cos(pkin(6));
t101 = t117 * t115;
t112 = sin(pkin(12));
t132 = sin(qJ(1));
t87 = -t74 * t101 + t132 * t112;
t113 = sin(pkin(7));
t114 = sin(pkin(6));
t98 = t114 * t113;
t144 = t87 * t116 + t74 * t98;
t100 = t117 * t112;
t52 = t74 * t100 + t132 * t115;
t71 = sin(qJ(3));
t28 = -t52 * t133 + t144 * t71;
t99 = t114 * t116;
t41 = -t87 * t113 + t74 * t99;
t67 = qJ(4) + pkin(13);
t64 = sin(t67);
t65 = cos(t67);
t10 = t28 * t65 + t41 * t64;
t25 = t144 * t133 + t52 * t71;
t69 = sin(qJ(6));
t72 = cos(qJ(6));
t148 = t10 * t69 + t25 * t72;
t147 = t10 * t72 - t25 * t69;
t70 = sin(qJ(4));
t128 = t41 * t70;
t73 = cos(qJ(4));
t143 = t28 * t73 + t128;
t137 = t28 * t70 - t41 * t73;
t142 = t28 * t64 - t41 * t65;
t81 = t132 * t101 + t74 * t112;
t43 = t81 * t113 + t132 * t99;
t135 = t113 * t117 + t115 * t99;
t97 = t114 * t112;
t40 = t133 * t97 + t135 * t71;
t51 = -t115 * t98 + t117 * t116;
t138 = -t40 * t70 + t51 * t73;
t136 = t81 * t116 - t132 * t98;
t53 = -t132 * t100 + t74 * t115;
t30 = t53 * t133 - t136 * t71;
t13 = -t30 * t70 + t43 * t73;
t126 = t43 * t70;
t123 = t65 * t69;
t122 = t65 * t72;
t63 = t73 * pkin(4) + pkin(3);
t68 = -qJ(5) - pkin(10);
t121 = -t25 * t63 + t28 * t68;
t29 = t136 * t133 + t53 * t71;
t120 = -t29 * t63 - t30 * t68;
t39 = -t135 * t133 + t71 * t97;
t119 = -t39 * t63 - t40 * t68;
t102 = t132 * t114;
t118 = t74 * pkin(1) + qJ(2) * t102;
t111 = t74 * t114;
t110 = t137 * pkin(4);
t109 = t13 * pkin(4);
t108 = t138 * pkin(4);
t11 = t30 * t64 - t43 * t65;
t107 = g(1) * t142 + g(2) * t11;
t106 = -t132 * pkin(1) + qJ(2) * t111;
t105 = -pkin(5) * t65 - pkin(11) * t64;
t104 = -g(1) * t25 + g(2) * t29;
t19 = -t40 * t64 + t51 * t65;
t93 = g(1) * t11 - g(2) * t142 - g(3) * t19;
t12 = t30 * t65 + t43 * t64;
t20 = t40 * t65 + t51 * t64;
t92 = g(1) * t12 - g(2) * t10 + g(3) * t20;
t91 = g(1) * t29 + g(2) * t25 + g(3) * t39;
t90 = g(1) * t30 - g(2) * t28 + g(3) * t40;
t78 = -t52 * pkin(2) + t41 * pkin(9) + t106;
t77 = t53 * pkin(2) + t43 * pkin(9) + t118;
t76 = pkin(4) * t128 + t25 * t68 + t28 * t63 + t78;
t75 = pkin(4) * t126 - t29 * t68 + t30 * t63 + t77;
t48 = -g(1) * t102 + g(2) * t111 - g(3) * t117;
t14 = t30 * t73 + t126;
t3 = t12 * t72 + t29 * t69;
t2 = -t12 * t69 + t29 * t72;
t1 = t91 * t64;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t132 - g(2) * t74, g(1) * t74 + g(2) * t132, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t52 - g(2) * t53, -g(1) * t87 + g(2) * t81, -g(1) * t111 - g(2) * t102, -g(1) * t106 - g(2) * t118, 0, 0, 0, 0, 0, 0, -g(1) * t28 - g(2) * t30, t104, -g(1) * t41 - g(2) * t43, -g(1) * t78 - g(2) * t77, 0, 0, 0, 0, 0, 0, -g(1) * t143 - g(2) * t14, g(1) * t137 - g(2) * t13, -t104, -g(1) * (t28 * pkin(3) - pkin(10) * t25 + t78) - g(2) * (t30 * pkin(3) + t29 * pkin(10) + t77) 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, t107, -t104, -g(1) * t76 - g(2) * t75, 0, 0, 0, 0, 0, 0, -g(1) * t147 - g(2) * t3, g(1) * t148 - g(2) * t2, -t107, -g(1) * (t10 * pkin(5) + pkin(11) * t142 + t76) - g(2) * (t12 * pkin(5) + t11 * pkin(11) + t75); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t90, 0, 0, 0, 0, 0, 0, 0, 0, t91 * t73, -t91 * t70, -t90, -g(1) * (-t29 * pkin(3) + t30 * pkin(10)) - g(2) * (-t25 * pkin(3) - pkin(10) * t28) - g(3) * (-t39 * pkin(3) + t40 * pkin(10)) 0, 0, 0, 0, 0, 0, t91 * t65, -t1, -t90, -g(1) * t120 - g(2) * t121 - g(3) * t119, 0, 0, 0, 0, 0, 0, -g(1) * (-t29 * t122 + t30 * t69) - g(2) * (-t25 * t122 - t28 * t69) - g(3) * (-t39 * t122 + t40 * t69) -g(1) * (t29 * t123 + t30 * t72) - g(2) * (t25 * t123 - t28 * t72) - g(3) * (t39 * t123 + t40 * t72) t1, -g(1) * (t105 * t29 + t120) - g(2) * (t105 * t25 + t121) - g(3) * (t105 * t39 + t119); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t137 - g(3) * t138, g(1) * t14 - g(2) * t143 - g(3) * (-t40 * t73 - t51 * t70) 0, 0, 0, 0, 0, 0, 0, 0, t93, t92, 0, -g(1) * t109 - g(2) * t110 - g(3) * t108, 0, 0, 0, 0, 0, 0, t93 * t72, -t93 * t69, -t92, -g(1) * (-t11 * pkin(5) + t12 * pkin(11) + t109) - g(2) * (pkin(5) * t142 - pkin(11) * t10 + t110) - g(3) * (t19 * pkin(5) + t20 * pkin(11) + t108); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t148 - g(3) * (-t20 * t69 + t39 * t72) g(1) * t3 - g(2) * t147 - g(3) * (-t20 * t72 - t39 * t69) 0, 0;];
taug_reg  = t4;

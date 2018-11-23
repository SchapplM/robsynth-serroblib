% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x38]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_gravloadJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From joint_gravload_fixb_regressor_minpar_matlab.m
t101 = sin(qJ(6));
t107 = cos(qJ(6));
t103 = sin(qJ(4));
t138 = pkin(8) + qJ(4);
t125 = sin(t138) / 0.2e1;
t139 = pkin(8) - qJ(4);
t131 = sin(t139);
t117 = t125 + t131 / 0.2e1;
t112 = cos(qJ(1));
t146 = cos(pkin(7));
t98 = sin(pkin(6));
t137 = t98 * t146;
t105 = sin(qJ(2));
t106 = sin(qJ(1));
t142 = pkin(6) + qJ(2);
t130 = cos(t142) / 0.2e1;
t143 = pkin(6) - qJ(2);
t136 = cos(t143);
t118 = t136 / 0.2e1 + t130;
t70 = t106 * t105 - t112 * t118;
t97 = sin(pkin(7));
t120 = -t112 * t137 + t70 * t97;
t128 = cos(t138) / 0.2e1;
t134 = cos(t139);
t124 = t134 / 0.2e1 + t128;
t110 = cos(qJ(3));
t147 = t112 * t98;
t111 = cos(qJ(2));
t127 = sin(t142) / 0.2e1;
t133 = sin(t143);
t81 = t127 - t133 / 0.2e1;
t71 = t106 * t111 + t112 * t81;
t140 = pkin(7) + qJ(3);
t126 = sin(t140) / 0.2e1;
t141 = pkin(7) - qJ(3);
t132 = sin(t141);
t79 = t126 - t132 / 0.2e1;
t129 = cos(t140) / 0.2e1;
t135 = cos(t141);
t83 = t129 - t135 / 0.2e1;
t48 = t71 * t110 + t83 * t147 - t70 * t79;
t104 = sin(qJ(3));
t78 = t126 + t132 / 0.2e1;
t84 = t135 / 0.2e1 + t129;
t50 = t71 * t104 + t78 * t147 + t70 * t84;
t14 = t48 * t103 - t120 * t117 + t50 * t124;
t102 = sin(qJ(5));
t108 = cos(qJ(5));
t109 = cos(qJ(4));
t77 = t125 - t131 / 0.2e1;
t82 = t128 - t134 / 0.2e1;
t114 = t48 * t109 - t120 * t82 - t50 * t77;
t96 = sin(pkin(8));
t99 = cos(pkin(8));
t39 = t120 * t99 + t50 * t96;
t4 = t39 * t102 + t108 * t114;
t157 = t4 * t101 - t14 * t107;
t156 = t14 * t101 + t4 * t107;
t153 = t102 * t114 - t39 * t108;
t152 = t82 * t97;
t151 = t97 * t99;
t150 = t102 * t96;
t149 = t106 * t98;
t148 = t108 * t96;
t145 = t101 * t108;
t144 = t107 * t108;
t75 = t106 * t81 - t112 * t111;
t100 = cos(pkin(6));
t80 = t127 + t133 / 0.2e1;
t122 = -t100 * t146 + t80 * t97;
t85 = t130 - t136 / 0.2e1;
t61 = t100 * t78 + t85 * t104 + t80 * t84;
t63 = t100 * t83 + t85 * t110 - t80 * t79;
t115 = -t109 * t63 + t122 * t82 + t61 * t77;
t46 = -t122 * t99 - t61 * t96;
t73 = -t112 * t105 - t106 * t118;
t119 = -t106 * t137 + t73 * t97;
t52 = t104 * t75 + t78 * t149 + t73 * t84;
t54 = t110 * t75 + t83 * t149 - t73 * t79;
t113 = -t109 * t54 + t119 * t82 + t52 * t77;
t41 = -t119 * t99 - t52 * t96;
t6 = -t102 * t113 + t41 * t108;
t123 = g(1) * t6 - g(2) * t153 + g(3) * (-t102 * t115 + t46 * t108);
t19 = -t103 * t54 + t119 * t117 - t52 * t124;
t32 = -t103 * t63 + t122 * t117 - t61 * t124;
t121 = g(1) * t19 + g(2) * t14 + g(3) * t32;
t116 = t97 * t117;
t67 = t80 * t110 + t85 * t79;
t66 = -t80 * t104 + t85 * t84;
t60 = t73 * t110 + t75 * t79;
t59 = -t73 * t104 + t75 * t84;
t58 = -t70 * t110 - t71 * t79;
t57 = t70 * t104 - t71 * t84;
t56 = -t85 * t151 - t66 * t96;
t43 = -t75 * t151 - t59 * t96;
t42 = t151 * t71 - t57 * t96;
t38 = t67 * t109 + t85 * t152 + t66 * t77;
t37 = t67 * t103 + t85 * t116 - t66 * t124;
t36 = t61 * t109 + t63 * t77;
t35 = t61 * t103 - t63 * t124;
t31 = t52 * t109 + t54 * t77;
t30 = t52 * t103 - t54 * t124;
t29 = -t109 * t50 - t48 * t77;
t28 = -t103 * t50 + t124 * t48;
t27 = t60 * t109 + t75 * t152 + t59 * t77;
t26 = t60 * t103 + t75 * t116 - t59 * t124;
t25 = t58 * t109 - t152 * t71 + t57 * t77;
t24 = t58 * t103 - t116 * t71 - t57 * t124;
t23 = t36 * t108 - t63 * t150;
t22 = t56 * t102 + t38 * t108;
t13 = t46 * t102 + t108 * t115;
t11 = t31 * t108 - t54 * t150;
t10 = t29 * t108 + t150 * t48;
t9 = t43 * t102 + t27 * t108;
t8 = t42 * t102 + t25 * t108;
t7 = t41 * t102 + t108 * t113;
t2 = t19 * t101 + t7 * t107;
t1 = -t7 * t101 + t19 * t107;
t3 = [0, g(1) * t106 - g(2) * t112, g(1) * t112 + g(2) * t106, 0, 0, 0, 0, 0, g(1) * t71 + g(2) * t75, -g(1) * t70 - g(2) * t73, 0, 0, 0, 0, 0, g(1) * t48 + g(2) * t54, -g(1) * t50 - g(2) * t52, 0, 0, 0, 0, 0, g(1) * t114 - g(2) * t113, -g(1) * t14 + g(2) * t19, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t7, -g(1) * t153 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t156 - g(2) * t2, -g(1) * t157 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t73 + g(2) * t70 - g(3) * t80, -g(1) * t75 + g(2) * t71 - g(3) * t85, 0, 0, 0, 0, 0, -g(1) * t60 - g(2) * t58 - g(3) * t67, -g(1) * t59 - g(2) * t57 - g(3) * t66, 0, 0, 0, 0, 0, -g(1) * t27 - g(2) * t25 - g(3) * t38, g(1) * t26 + g(2) * t24 + g(3) * t37, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t8 - g(3) * t22, -g(1) * (-t27 * t102 + t43 * t108) - g(2) * (-t25 * t102 + t42 * t108) - g(3) * (-t38 * t102 + t56 * t108) 0, 0, 0, 0, 0, -g(1) * (t26 * t101 + t9 * t107) - g(2) * (t24 * t101 + t8 * t107) - g(3) * (t37 * t101 + t22 * t107) -g(1) * (-t9 * t101 + t26 * t107) - g(2) * (-t8 * t101 + t24 * t107) - g(3) * (-t22 * t101 + t37 * t107); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t52 + g(2) * t50 - g(3) * t61, -g(1) * t54 + g(2) * t48 - g(3) * t63, 0, 0, 0, 0, 0, -g(1) * t31 - g(2) * t29 - g(3) * t36, g(1) * t30 + g(2) * t28 + g(3) * t35, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t10 - g(3) * t23, -g(1) * (-t31 * t102 - t54 * t148) - g(2) * (-t29 * t102 + t148 * t48) - g(3) * (-t36 * t102 - t63 * t148) 0, 0, 0, 0, 0, -g(1) * (t30 * t101 + t11 * t107) - g(2) * (t10 * t107 + t28 * t101) - g(3) * (t35 * t101 + t23 * t107) -g(1) * (-t11 * t101 + t30 * t107) - g(2) * (-t10 * t101 + t28 * t107) - g(3) * (-t23 * t101 + t35 * t107); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, g(1) * t113 + g(2) * t114 + g(3) * t115, 0, 0, 0, 0, 0, t121 * t108, -t121 * t102, 0, 0, 0, 0, 0, -g(1) * (t101 * t113 - t19 * t144) - g(2) * (t101 * t114 - t14 * t144) - g(3) * (t101 * t115 - t32 * t144) -g(1) * (t107 * t113 + t19 * t145) - g(2) * (t107 * t114 + t14 * t145) - g(3) * (t107 * t115 + t32 * t145); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, g(1) * t7 + g(2) * t4 + g(3) * t13, 0, 0, 0, 0, 0, -t123 * t107, t123 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t157 - g(3) * (-t13 * t101 + t32 * t107) g(1) * t2 + g(2) * t156 - g(3) * (-t32 * t101 - t13 * t107);];
taug_reg  = t3;

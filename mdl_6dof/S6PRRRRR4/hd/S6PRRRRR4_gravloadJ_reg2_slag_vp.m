% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t123 = cos(pkin(7));
t79 = sin(qJ(3));
t109 = t79 * t123;
t124 = cos(pkin(6));
t75 = sin(pkin(7));
t112 = t75 * t124;
t141 = sin(qJ(2));
t142 = cos(qJ(3));
t143 = cos(qJ(2));
t76 = sin(pkin(6));
t45 = t79 * t112 + (t143 * t109 + t141 * t142) * t76;
t116 = t76 * t143;
t59 = -t75 * t116 + t124 * t123;
t78 = sin(qJ(4));
t81 = cos(qJ(4));
t152 = -t45 * t78 + t59 * t81;
t121 = sin(pkin(13));
t110 = t76 * t121;
t122 = cos(pkin(13));
t94 = t124 * t121;
t86 = t122 * t141 + t143 * t94;
t148 = t75 * t110 - t86 * t123;
t61 = t122 * t143 - t141 * t94;
t29 = t61 * t142 + t148 * t79;
t47 = t123 * t110 + t86 * t75;
t151 = -t29 * t78 + t47 * t81;
t111 = t76 * t122;
t95 = t124 * t122;
t85 = t121 * t141 - t143 * t95;
t149 = t75 * t111 + t85 * t123;
t60 = t121 * t143 + t141 * t95;
t27 = t60 * t142 - t149 * t79;
t46 = -t123 * t111 + t85 * t75;
t150 = -t27 * t78 + t46 * t81;
t147 = -g(1) * t61 - g(2) * t60;
t146 = pkin(9) * t75;
t74 = qJ(4) + qJ(5);
t72 = sin(t74);
t134 = t72 * t75;
t73 = cos(t74);
t133 = t73 * t75;
t77 = sin(qJ(6));
t132 = t73 * t77;
t80 = cos(qJ(6));
t131 = t73 * t80;
t130 = t75 * t78;
t129 = t75 * t81;
t26 = t149 * t142 + t60 * t79;
t71 = t81 * pkin(4) + pkin(3);
t82 = -pkin(11) - pkin(10);
t128 = -t26 * t71 - t27 * t82;
t28 = -t148 * t142 + t61 * t79;
t127 = -t28 * t71 - t29 * t82;
t115 = t76 * t141;
t98 = t123 * t142;
t44 = -t142 * t112 + t79 * t115 - t98 * t116;
t126 = -t44 * t71 - t45 * t82;
t108 = t75 * t115;
t125 = pkin(2) * t116 + pkin(9) * t108;
t34 = t60 * t98 - t85 * t79;
t35 = -t60 * t109 - t85 * t142;
t57 = t85 * pkin(2);
t120 = -t34 * t82 + t35 * t71 - t57;
t36 = t61 * t98 - t86 * t79;
t37 = -t61 * t109 - t86 * t142;
t58 = t86 * pkin(2);
t119 = -t36 * t82 + t37 * t71 - t58;
t11 = -t27 * t72 + t46 * t73;
t12 = t27 * t73 + t46 * t72;
t117 = t11 * pkin(5) + t12 * pkin(12);
t13 = -t29 * t72 + t47 * t73;
t14 = t29 * t73 + t47 * t72;
t114 = t13 * pkin(5) + t14 * pkin(12);
t22 = -t45 * t72 + t59 * t73;
t23 = t45 * t73 + t59 * t72;
t113 = t22 * pkin(5) + t23 * pkin(12);
t107 = t150 * pkin(4);
t106 = t151 * pkin(4);
t105 = t152 * pkin(4);
t104 = t60 * t146 - t57;
t103 = t61 * t146 - t58;
t97 = t123 * t141;
t54 = (t142 * t97 + t143 * t79) * t76;
t55 = (t143 * t142 - t79 * t97) * t76;
t96 = t78 * t108;
t100 = pkin(4) * t96 - t54 * t82 + t55 * t71 + t125;
t99 = -pkin(5) * t73 - pkin(12) * t72;
t93 = g(1) * t13 + g(2) * t11 + g(3) * t22;
t5 = g(1) * t14 + g(2) * t12 + g(3) * t23;
t15 = -t60 * t133 + t35 * t72;
t17 = -t61 * t133 + t37 * t72;
t38 = -t73 * t108 + t55 * t72;
t92 = g(1) * t17 + g(2) * t15 + g(3) * t38;
t91 = g(1) * t28 + g(2) * t26 + g(3) * t44;
t90 = g(1) * t29 + g(2) * t27 + g(3) * t45;
t89 = g(1) * t36 + g(2) * t34 + g(3) * t54;
t88 = -g(3) * t115 + t147;
t87 = t147 * (pkin(4) * t78 + pkin(9)) * t75;
t39 = t72 * t108 + t55 * t73;
t18 = t61 * t134 + t37 * t73;
t16 = t60 * t134 + t35 * t73;
t6 = t91 * t72;
t2 = t93 * t80;
t1 = t93 * t77;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t86 + g(2) * t85 - g(3) * t116, -t88, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t37 - g(2) * t35 - g(3) * t55, t89, t88 * t75, -g(1) * t103 - g(2) * t104 - g(3) * t125, 0, 0, 0, 0, 0, 0, -g(1) * (t61 * t130 + t37 * t81) - g(2) * (t60 * t130 + t35 * t81) - g(3) * (t55 * t81 + t96) -g(1) * (t61 * t129 - t37 * t78) - g(2) * (t60 * t129 - t35 * t78) - g(3) * (t81 * t108 - t55 * t78) -t89, -g(1) * (t37 * pkin(3) + t36 * pkin(10) + t103) - g(2) * (t35 * pkin(3) + t34 * pkin(10) + t104) - g(3) * (t55 * pkin(3) + t54 * pkin(10) + t125) 0, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t16 - g(3) * t39, t92, -t89, -g(1) * t119 - g(2) * t120 - g(3) * t100 + t87, 0, 0, 0, 0, 0, 0, -g(1) * (t18 * t80 + t36 * t77) - g(2) * (t16 * t80 + t34 * t77) - g(3) * (t39 * t80 + t54 * t77) -g(1) * (-t18 * t77 + t36 * t80) - g(2) * (-t16 * t77 + t34 * t80) - g(3) * (-t39 * t77 + t54 * t80) -t92, -g(1) * (t18 * pkin(5) + t17 * pkin(12) + t119) - g(2) * (t16 * pkin(5) + t15 * pkin(12) + t120) - g(3) * (t39 * pkin(5) + t38 * pkin(12) + t100) + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t90, 0, 0, 0, 0, 0, 0, 0, 0, t91 * t81, -t91 * t78, -t90, -g(1) * (-t28 * pkin(3) + t29 * pkin(10)) - g(2) * (-t26 * pkin(3) + t27 * pkin(10)) - g(3) * (-t44 * pkin(3) + t45 * pkin(10)) 0, 0, 0, 0, 0, 0, t91 * t73, -t6, -t90, -g(1) * t127 - g(2) * t128 - g(3) * t126, 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t131 + t29 * t77) - g(2) * (-t26 * t131 + t27 * t77) - g(3) * (-t44 * t131 + t45 * t77) -g(1) * (t28 * t132 + t29 * t80) - g(2) * (t26 * t132 + t27 * t80) - g(3) * (t44 * t132 + t45 * t80) t6, -g(1) * (t99 * t28 + t127) - g(2) * (t99 * t26 + t128) - g(3) * (t99 * t44 + t126); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t151 - g(2) * t150 - g(3) * t152, -g(1) * (-t29 * t81 - t47 * t78) - g(2) * (-t27 * t81 - t46 * t78) - g(3) * (-t45 * t81 - t59 * t78) 0, 0, 0, 0, 0, 0, 0, 0, -t93, t5, 0, -g(1) * t106 - g(2) * t107 - g(3) * t105, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (t106 + t114) - g(2) * (t107 + t117) - g(3) * (t105 + t113); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * t114 - g(2) * t117 - g(3) * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t77 + t28 * t80) - g(2) * (-t12 * t77 + t26 * t80) - g(3) * (-t23 * t77 + t44 * t80) -g(1) * (-t14 * t80 - t28 * t77) - g(2) * (-t12 * t80 - t26 * t77) - g(3) * (-t23 * t80 - t44 * t77) 0, 0;];
taug_reg  = t3;

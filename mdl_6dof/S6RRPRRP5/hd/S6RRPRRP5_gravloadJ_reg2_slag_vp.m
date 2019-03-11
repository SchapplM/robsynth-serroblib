% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t58 = sin(pkin(11));
t64 = sin(qJ(2));
t68 = cos(qJ(2));
t99 = cos(pkin(11));
t45 = -t68 * t58 - t64 * t99;
t59 = sin(pkin(6));
t69 = cos(qJ(1));
t112 = t59 * t69;
t60 = cos(pkin(6));
t101 = t45 * t60;
t65 = sin(qJ(1));
t81 = -t64 * t58 + t68 * t99;
t24 = t69 * t101 - t65 * t81;
t63 = sin(qJ(4));
t67 = cos(qJ(4));
t15 = -t63 * t112 - t24 * t67;
t74 = t81 * t60;
t25 = t65 * t45 + t69 * t74;
t62 = sin(qJ(5));
t66 = cos(qJ(5));
t121 = t15 * t62 + t25 * t66;
t120 = t15 * t66 - t25 * t62;
t113 = t59 * t68;
t103 = t69 * t64;
t107 = t65 * t68;
t41 = -t60 * t107 - t103;
t119 = -g(1) * t41 - g(3) * t113;
t114 = t59 * t65;
t29 = t65 * t101 + t69 * t81;
t19 = t63 * t114 + t29 * t67;
t28 = t69 * t45 - t65 * t74;
t12 = -t19 * t62 - t28 * t66;
t37 = t45 * t59;
t31 = -t37 * t67 + t60 * t63;
t36 = t81 * t59;
t1 = g(2) * t121 - g(3) * (-t31 * t62 - t36 * t66) - g(1) * t12;
t111 = t62 * t67;
t108 = t65 * t64;
t106 = t66 * t67;
t102 = t69 * t68;
t100 = pkin(2) * t113 + t36 * pkin(3);
t97 = t60 * t102;
t96 = pkin(5) * t62 + pkin(9);
t38 = t60 * t64 * pkin(2) + (-pkin(8) - qJ(3)) * t59;
t57 = t68 * pkin(2) + pkin(1);
t95 = -t65 * t38 + t69 * t57;
t92 = -t37 * pkin(9) + t100;
t91 = t29 * pkin(3) + t95;
t90 = pkin(4) * t67 + pkin(10) * t63;
t14 = t67 * t112 - t24 * t63;
t18 = -t67 * t114 + t29 * t63;
t89 = -g(1) * t14 + g(2) * t18;
t88 = g(1) * t25 - g(2) * t28;
t87 = g(1) * t69 + g(2) * t65;
t86 = g(1) * t65 - g(2) * t69;
t85 = -t69 * t38 - t65 * t57;
t56 = t66 * pkin(5) + pkin(4);
t61 = -qJ(6) - pkin(10);
t84 = t56 * t67 - t61 * t63;
t46 = pkin(2) * t97;
t83 = -pkin(2) * t108 + t25 * pkin(3) + t46;
t82 = pkin(3) * t24 + t85;
t80 = -t28 * pkin(9) + t91;
t30 = -t37 * t63 - t60 * t67;
t79 = g(1) * t18 + g(2) * t14 + g(3) * t30;
t78 = g(1) * t19 + g(2) * t15 + g(3) * t31;
t77 = -g(1) * t29 + g(2) * t24 + g(3) * t37;
t76 = g(1) * t28 + g(2) * t25 + g(3) * t36;
t75 = -t24 * pkin(9) + t83;
t73 = t25 * pkin(9) + t82;
t72 = t41 * pkin(2) + t28 * pkin(3);
t71 = pkin(9) * t29 + t72;
t43 = t87 * t59;
t42 = -t60 * t108 + t102;
t40 = -t60 * t103 - t107;
t39 = -t97 + t108;
t35 = -g(3) * t60 - t86 * t59;
t13 = t19 * t66 - t28 * t62;
t10 = t76 * t63;
t8 = t79 * t66;
t7 = t79 * t62;
t6 = g(1) * t120 - g(2) * t13;
t5 = -g(1) * t121 - g(2) * t12;
t4 = -g(1) * (t28 * t106 + t29 * t62) - g(2) * (t25 * t106 - t24 * t62) - g(3) * (t36 * t106 - t37 * t62);
t3 = -g(1) * (-t28 * t111 + t29 * t66) - g(2) * (-t25 * t111 - t24 * t66) - g(3) * (-t36 * t111 - t37 * t66);
t2 = g(1) * t13 + g(2) * t120 - g(3) * (-t31 * t66 + t36 * t62);
t9 = [0, 0, 0, 0, 0, 0, t86, t87, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t40 - g(2) * t42, -g(1) * t39 - g(2) * t41, -t43, -g(1) * (-t65 * pkin(1) + pkin(8) * t112) - g(2) * (t69 * pkin(1) + pkin(8) * t114) 0, 0, 0, 0, 0, 0, -g(1) * t24 - g(2) * t29, t88, -t43, -g(1) * t85 - g(2) * t95, 0, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t19, t89, -t88, -g(1) * t73 - g(2) * t80, 0, 0, 0, 0, 0, 0, t6, t5, -t89, -g(1) * (-pkin(4) * t15 - pkin(10) * t14 + t73) - g(2) * (t19 * pkin(4) + t18 * pkin(10) + t80) 0, 0, 0, 0, 0, 0, t6, t5, -t89, -g(1) * (t14 * t61 - t15 * t56 + t96 * t25 + t82) - g(2) * (-t18 * t61 + t19 * t56 - t96 * t28 + t91); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t39 + t119, g(3) * t59 * t64 + g(1) * t42 - g(2) * t40, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t77, 0, -g(2) * t46 + (g(2) * t108 + t119) * pkin(2), 0, 0, 0, 0, 0, 0, -t76 * t67, t10, t77, -g(1) * t71 - g(2) * t75 - g(3) * t92, 0, 0, 0, 0, 0, 0, t4, t3, -t10, -g(1) * (t90 * t28 + t71) - g(2) * (t90 * t25 + t75) - g(3) * (t90 * t36 + t92) 0, 0, 0, 0, 0, 0, t4, t3, -t10, -g(1) * (t84 * t28 + t29 * t96 + t72) - g(2) * (-t96 * t24 + t84 * t25 + t83) - g(3) * (t84 * t36 - t96 * t37 + t100); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t78, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t78, -g(1) * (-t18 * pkin(4) + t19 * pkin(10)) - g(2) * (-t14 * pkin(4) + t15 * pkin(10)) - g(3) * (-t30 * pkin(4) + t31 * pkin(10)) 0, 0, 0, 0, 0, 0, t8, -t7, -t78, -g(1) * (-t18 * t56 - t19 * t61) - g(2) * (-t14 * t56 - t15 * t61) - g(3) * (-t30 * t56 - t31 * t61); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79;];
taug_reg  = t9;

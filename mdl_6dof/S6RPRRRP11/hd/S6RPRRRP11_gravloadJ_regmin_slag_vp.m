% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t43 = cos(qJ(5));
t79 = cos(pkin(12));
t81 = cos(pkin(6));
t69 = t81 * t79;
t76 = sin(pkin(12));
t87 = sin(qJ(1));
t89 = cos(qJ(1));
t55 = -t89 * t69 + t87 * t76;
t77 = sin(pkin(7));
t78 = sin(pkin(6));
t65 = t78 * t77;
t80 = cos(pkin(7));
t103 = t55 * t80 + t89 * t65;
t67 = t81 * t76;
t28 = t89 * t67 + t87 * t79;
t42 = sin(qJ(3));
t88 = cos(qJ(3));
t15 = t103 * t42 - t28 * t88;
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t66 = t80 * t78;
t96 = t55 * t77 - t89 * t66;
t7 = t15 * t44 - t41 * t96;
t12 = t103 * t88 + t28 * t42;
t40 = sin(qJ(5));
t86 = t12 * t40;
t106 = t7 * t43 - t86;
t105 = t12 * t43 + t7 * t40;
t6 = t15 * t41 + t44 * t96;
t50 = t87 * t69 + t89 * t76;
t99 = t50 * t80 - t87 * t65;
t64 = t78 * t76;
t97 = t79 * t66 + t81 * t77;
t22 = t97 * t42 + t88 * t64;
t49 = -t79 * t65 + t81 * t80;
t11 = t22 * t44 + t41 * t49;
t29 = -t87 * t67 + t89 * t79;
t16 = t29 * t42 + t99 * t88;
t17 = t29 * t88 - t99 * t42;
t45 = -t50 * t77 - t87 * t66;
t9 = t17 * t44 - t41 * t45;
t2 = t16 * t43 - t9 * t40;
t21 = t42 * t64 - t97 * t88;
t98 = -g(1) * t2 - g(2) * t105 - g(3) * (-t11 * t40 + t21 * t43);
t60 = g(1) * t16 + g(2) * t12 + g(3) * t21;
t95 = g(1) * t17 - g(2) * t15 + g(3) * t22;
t85 = t16 * t40;
t84 = t40 * t44;
t83 = t43 * t44;
t70 = t87 * t78;
t82 = t89 * pkin(1) + qJ(2) * t70;
t8 = t17 * t41 + t44 * t45;
t74 = g(1) * t6 + g(2) * t8;
t71 = t89 * t78;
t73 = -t87 * pkin(1) + qJ(2) * t71;
t10 = t22 * t41 - t44 * t49;
t62 = g(1) * t8 - g(2) * t6 + g(3) * t10;
t61 = g(1) * t9 - g(2) * t7 + g(3) * t11;
t39 = -qJ(6) - pkin(11);
t37 = t43 * pkin(5) + pkin(4);
t27 = -g(1) * t70 + g(2) * t71 - g(3) * t81;
t3 = t9 * t43 + t85;
t1 = t60 * t41;
t4 = [0, g(1) * t87 - g(2) * t89, g(1) * t89 + g(2) * t87, g(1) * t28 - g(2) * t29, -g(1) * t55 + g(2) * t50, -g(1) * t71 - g(2) * t70, -g(1) * t73 - g(2) * t82, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t12 + g(2) * t16, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, t74, 0, 0, 0, 0, 0, -g(1) * t106 - g(2) * t3, g(1) * t105 - g(2) * t2, -t74, -g(1) * (-t28 * pkin(2) + t15 * pkin(3) - pkin(5) * t86 - pkin(10) * t12 + t7 * t37 - t6 * t39 + t73) - g(2) * (t29 * pkin(2) + t17 * pkin(3) + pkin(5) * t85 + t16 * pkin(10) + t9 * t37 - t8 * t39 + t82) + (g(1) * t96 + g(2) * t45) * pkin(9); 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t95, 0, 0, 0, 0, 0, t60 * t44, -t1, 0, 0, 0, 0, 0, -g(1) * (-t16 * t83 + t17 * t40) - g(2) * (-t12 * t83 - t15 * t40) - g(3) * (-t21 * t83 + t22 * t40) -g(1) * (t16 * t84 + t17 * t43) - g(2) * (t12 * t84 - t15 * t43) - g(3) * (t21 * t84 + t22 * t43) t1, -t95 * (t40 * pkin(5) + pkin(10)) + t60 * (t37 * t44 - t39 * t41 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t61, 0, 0, 0, 0, 0, t62 * t43, -t62 * t40, -t61, -g(1) * (-t8 * t37 - t9 * t39) - g(2) * (t37 * t6 + t39 * t7) - g(3) * (-t10 * t37 - t11 * t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, g(1) * t3 - g(2) * t106 - g(3) * (-t11 * t43 - t21 * t40) 0, t98 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62;];
taug_reg  = t4;

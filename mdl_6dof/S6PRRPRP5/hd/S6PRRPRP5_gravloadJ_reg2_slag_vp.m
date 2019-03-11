% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t64 = sin(qJ(2));
t67 = cos(qJ(2));
t92 = cos(pkin(10));
t93 = cos(pkin(6));
t78 = t93 * t92;
t91 = sin(pkin(10));
t44 = t64 * t78 + t91 * t67;
t63 = sin(qJ(3));
t66 = cos(qJ(3));
t61 = sin(pkin(6));
t86 = t61 * t92;
t23 = t44 * t66 - t63 * t86;
t77 = t93 * t91;
t46 = -t64 * t77 + t92 * t67;
t85 = t61 * t91;
t25 = t46 * t66 + t63 * t85;
t100 = t61 * t64;
t48 = t66 * t100 + t93 * t63;
t72 = g(1) * t25 + g(2) * t23 + g(3) * t48;
t22 = t44 * t63 + t66 * t86;
t105 = t22 * pkin(9);
t24 = t46 * t63 - t66 * t85;
t104 = t24 * pkin(9);
t47 = t63 * t100 - t93 * t66;
t103 = t47 * pkin(9);
t43 = t91 * t64 - t67 * t78;
t102 = t43 * t66;
t45 = t92 * t64 + t67 * t77;
t101 = t45 * t66;
t99 = t61 * t67;
t62 = sin(qJ(5));
t98 = t62 * t63;
t97 = t62 * t67;
t65 = cos(qJ(5));
t96 = t63 * t65;
t95 = pkin(2) * t99 + pkin(8) * t100;
t94 = qJ(4) * t63;
t90 = t66 * t99;
t89 = t65 * t99;
t88 = -t43 * pkin(2) + t44 * pkin(8);
t87 = -t45 * pkin(2) + t46 * pkin(8);
t20 = t22 * pkin(3);
t84 = t23 * qJ(4) - t20;
t21 = t24 * pkin(3);
t83 = t25 * qJ(4) - t21;
t42 = t47 * pkin(3);
t82 = t48 * qJ(4) - t42;
t81 = pkin(3) * t90 + t94 * t99 + t95;
t80 = -pkin(3) * t102 - t43 * t94 + t88;
t79 = -pkin(3) * t101 - t45 * t94 + t87;
t76 = pkin(4) * t100 + pkin(9) * t90 + t81;
t26 = t47 * t65 + t61 * t97;
t7 = -t22 * t65 + t43 * t62;
t9 = -t24 * t65 + t45 * t62;
t1 = g(1) * t9 + g(2) * t7 - g(3) * t26;
t10 = t24 * t62 + t45 * t65;
t27 = -t47 * t62 + t89;
t8 = t22 * t62 + t43 * t65;
t74 = g(1) * t10 + g(2) * t8 - g(3) * t27;
t14 = t43 * t96 + t44 * t62;
t16 = t45 * t96 + t46 * t62;
t30 = t62 * t100 - t63 * t89;
t73 = g(1) * t16 + g(2) * t14 + g(3) * t30;
t6 = g(1) * t24 + g(2) * t22 + g(3) * t47;
t71 = t44 * pkin(4) - pkin(9) * t102 + t80;
t70 = t46 * pkin(4) - pkin(9) * t101 + t79;
t69 = -g(1) * t45 - g(2) * t43 + g(3) * t99;
t68 = g(1) * t46 + g(2) * t44 + g(3) * t100;
t31 = (t63 * t97 + t64 * t65) * t61;
t17 = -t45 * t98 + t46 * t65;
t15 = -t43 * t98 + t44 * t65;
t12 = t69 * t66;
t11 = t69 * t63;
t4 = t72 * t65;
t3 = t72 * t62;
t2 = -g(1) * t17 - g(2) * t15 - g(3) * t31;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t68, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, -t68, -g(1) * t87 - g(2) * t88 - g(3) * t95, 0, 0, 0, 0, 0, 0, -t68, t12, -t11, -g(1) * t79 - g(2) * t80 - g(3) * t81, 0, 0, 0, 0, 0, 0, t2, t73, -t12, -g(1) * t70 - g(2) * t71 - g(3) * t76, 0, 0, 0, 0, 0, 0, t2, -t12, -t73, -g(1) * (t17 * pkin(5) + t16 * qJ(6) + t70) - g(2) * (t15 * pkin(5) + t14 * qJ(6) + t71) - g(3) * (t31 * pkin(5) + t30 * qJ(6) + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t72, -g(1) * t83 - g(2) * t84 - g(3) * t82, 0, 0, 0, 0, 0, 0, -t3, -t4, t6, -g(1) * (t83 - t104) - g(2) * (t84 - t105) - g(3) * (t82 - t103) 0, 0, 0, 0, 0, 0, -t3, t6, t4, -g(1) * (-t21 - t104) - g(2) * (-t20 - t105) - g(3) * (-t42 - t103) - t72 * (pkin(5) * t62 - qJ(6) * t65 + qJ(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t74, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t74, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - g(3) * (t26 * pkin(5) - t27 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;

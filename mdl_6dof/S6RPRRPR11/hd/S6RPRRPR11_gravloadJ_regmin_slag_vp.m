% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR11
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
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t76 = sin(pkin(12));
t81 = cos(pkin(6));
t68 = t81 * t76;
t79 = cos(pkin(12));
t87 = sin(qJ(1));
t89 = cos(qJ(1));
t28 = t89 * t68 + t87 * t79;
t44 = sin(qJ(3));
t88 = cos(qJ(3));
t70 = t81 * t79;
t55 = -t89 * t70 + t87 * t76;
t77 = sin(pkin(7));
t78 = sin(pkin(6));
t66 = t78 * t77;
t80 = cos(pkin(7));
t99 = t55 * t80 + t89 * t66;
t12 = t28 * t44 + t99 * t88;
t40 = pkin(13) + qJ(6);
t37 = sin(t40);
t38 = cos(t40);
t15 = -t28 * t88 + t99 * t44;
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t67 = t80 * t78;
t93 = t55 * t77 - t89 * t67;
t7 = t15 * t45 - t43 * t93;
t103 = t12 * t38 + t7 * t37;
t102 = -t12 * t37 + t7 * t38;
t6 = t15 * t43 + t45 * t93;
t51 = t87 * t70 + t89 * t76;
t95 = t51 * t80 - t87 * t66;
t94 = t79 * t67 + t81 * t77;
t29 = -t87 * t68 + t89 * t79;
t16 = t29 * t44 + t95 * t88;
t65 = t78 * t76;
t21 = t44 * t65 - t94 * t88;
t61 = g(1) * t16 + g(2) * t12 + g(3) * t21;
t86 = t37 * t45;
t85 = t38 * t45;
t41 = sin(pkin(13));
t84 = t41 * t45;
t42 = cos(pkin(13));
t83 = t42 * t45;
t71 = t87 * t78;
t82 = t89 * pkin(1) + qJ(2) * t71;
t17 = t29 * t88 - t95 * t44;
t46 = -t51 * t77 - t87 * t67;
t8 = t17 * t43 + t45 * t46;
t75 = g(1) * t6 + g(2) * t8;
t72 = t89 * t78;
t74 = -t87 * pkin(1) + qJ(2) * t72;
t22 = t94 * t44 + t88 * t65;
t50 = -t79 * t66 + t81 * t80;
t10 = t22 * t43 - t50 * t45;
t63 = g(1) * t8 - g(2) * t6 + g(3) * t10;
t11 = t22 * t45 + t50 * t43;
t9 = t17 * t45 - t43 * t46;
t62 = g(1) * t9 - g(2) * t7 + g(3) * t11;
t60 = g(1) * t17 - g(2) * t15 + g(3) * t22;
t27 = -g(1) * t71 + g(2) * t72 - g(3) * t81;
t3 = t61 * t43;
t2 = t16 * t37 + t9 * t38;
t1 = t16 * t38 - t9 * t37;
t4 = [0, g(1) * t87 - g(2) * t89, g(1) * t89 + g(2) * t87, g(1) * t28 - g(2) * t29, -g(1) * t55 + g(2) * t51, -g(1) * t72 - g(2) * t71, -g(1) * t74 - g(2) * t82, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t12 + g(2) * t16, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, t75, -g(1) * (-t12 * t41 + t7 * t42) - g(2) * (t16 * t41 + t9 * t42) -g(1) * (-t12 * t42 - t7 * t41) - g(2) * (t16 * t42 - t9 * t41) -t75, -g(1) * (-t28 * pkin(2) + t15 * pkin(3) + t7 * pkin(4) - pkin(10) * t12 + t6 * qJ(5) + t74) - g(2) * (t29 * pkin(2) + t17 * pkin(3) + t9 * pkin(4) + t16 * pkin(10) + t8 * qJ(5) + t82) + (g(1) * t93 + g(2) * t46) * pkin(9), 0, 0, 0, 0, 0, -g(1) * t102 - g(2) * t2, g(1) * t103 - g(2) * t1; 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t60, 0, 0, 0, 0, 0, t61 * t45, -t3, -g(1) * (-t16 * t83 + t17 * t41) - g(2) * (-t12 * t83 - t15 * t41) - g(3) * (-t21 * t83 + t22 * t41) -g(1) * (t16 * t84 + t17 * t42) - g(2) * (t12 * t84 - t15 * t42) - g(3) * (t21 * t84 + t22 * t42) t3, -t60 * pkin(10) + t61 * (pkin(4) * t45 + qJ(5) * t43 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (-t16 * t85 + t17 * t37) - g(2) * (-t12 * t85 - t15 * t37) - g(3) * (-t21 * t85 + t22 * t37) -g(1) * (t16 * t86 + t17 * t38) - g(2) * (t12 * t86 - t15 * t38) - g(3) * (t21 * t86 + t22 * t38); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t62, t63 * t42, -t63 * t41, -t62, -g(1) * (-t8 * pkin(4) + t9 * qJ(5)) - g(2) * (pkin(4) * t6 - qJ(5) * t7) - g(3) * (-t10 * pkin(4) + t11 * qJ(5)) 0, 0, 0, 0, 0, t63 * t38, -t63 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t103 - g(3) * (-t11 * t37 + t21 * t38) g(1) * t2 - g(2) * t102 - g(3) * (-t11 * t38 - t21 * t37);];
taug_reg  = t4;

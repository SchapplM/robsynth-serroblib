% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t88 = -g(1) * t42 + g(2) * t45;
t41 = sin(qJ(3));
t43 = cos(qJ(4));
t60 = t45 * t43;
t40 = sin(qJ(4));
t70 = t42 * t40;
t14 = -t41 * t70 + t60;
t61 = t45 * t40;
t69 = t42 * t43;
t16 = t41 * t61 + t69;
t44 = cos(qJ(3));
t78 = g(3) * t44;
t87 = -g(1) * t14 - g(2) * t16 + t40 * t78;
t39 = qJ(4) + qJ(5);
t29 = sin(t39);
t63 = t45 * t29;
t30 = cos(t39);
t71 = t42 * t30;
t11 = t41 * t63 + t71;
t62 = t45 * t30;
t72 = t42 * t29;
t9 = -t41 * t72 + t62;
t3 = -g(1) * t9 - g(2) * t11 + t29 * t78;
t48 = -g(3) * t41 - t44 * t88;
t85 = -pkin(1) - pkin(7);
t46 = -pkin(9) - pkin(8);
t77 = t40 * pkin(4);
t76 = t41 * t42;
t75 = t41 * t45;
t33 = qJ(6) + t39;
t26 = sin(t33);
t74 = t42 * t26;
t27 = cos(t33);
t73 = t42 * t27;
t68 = t44 * t45;
t67 = t44 * t46;
t20 = pkin(5) * t29 + t77;
t66 = t45 * t20;
t65 = t45 * t26;
t64 = t45 * t27;
t34 = t43 * pkin(4);
t21 = pkin(5) * t30 + t34;
t59 = t45 * pkin(1) + t42 * qJ(2);
t56 = t45 * pkin(7) + t59;
t55 = g(2) * t56;
t53 = t41 * pkin(3) - t44 * pkin(8);
t23 = g(1) * t45 + g(2) * t42;
t19 = pkin(3) + t21;
t38 = -pkin(10) + t46;
t51 = t41 * t19 + t44 * t38;
t28 = t34 + pkin(3);
t49 = t41 * t28 + t67;
t32 = t45 * qJ(2);
t18 = t23 * t44;
t17 = t41 * t60 - t70;
t15 = t41 * t69 + t61;
t13 = g(1) * t76 - g(2) * t75 + t78;
t12 = t41 * t62 - t72;
t10 = t41 * t71 + t63;
t8 = t41 * t64 - t74;
t7 = t41 * t65 + t73;
t6 = t41 * t73 + t65;
t5 = -t41 * t74 + t64;
t4 = g(1) * t10 - g(2) * t12 + t30 * t78;
t2 = g(1) * t6 - g(2) * t8 + t27 * t78;
t1 = -g(1) * t5 - g(2) * t7 + t26 * t78;
t22 = [0, 0, 0, 0, 0, 0, -t88, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, -t23, -g(1) * (-t42 * pkin(1) + t32) - g(2) * t59, 0, 0, 0, 0, 0, 0, -t23 * t41, -t18, -t88, -g(1) * (t85 * t42 + t32) - t55, 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t15, g(1) * t16 - g(2) * t14, t18, -g(1) * (pkin(3) * t75 - pkin(8) * t68 + t32) - t55 + (-g(1) * t85 - g(2) * t53) * t42, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t18, -g(1) * (t28 * t75 + t45 * t67 + t32) - g(2) * (pkin(4) * t61 + t56) + (-g(1) * (-t77 + t85) - g(2) * t49) * t42, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5, t18, -g(1) * (t19 * t75 + t38 * t68 + t32) - g(2) * (t56 + t66) + (-g(1) * (-t20 + t85) - g(2) * t51) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t13, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t43, t48 * t40, -t13, g(3) * t53 + t88 * (pkin(3) * t44 + pkin(8) * t41) 0, 0, 0, 0, 0, 0, -t48 * t30, t48 * t29, -t13, g(3) * t49 + t88 * (t28 * t44 - t41 * t46) 0, 0, 0, 0, 0, 0, -t48 * t27, t48 * t26, -t13, g(3) * t51 + t88 * (t19 * t44 - t38 * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, g(1) * t15 - g(2) * t17 + t43 * t78, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t87 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t20 * t76 + t45 * t21) - g(2) * (t42 * t21 + t41 * t66) + t20 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t22;

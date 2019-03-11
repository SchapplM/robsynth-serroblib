% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t44 = qJ(2) + qJ(3);
t40 = qJ(4) + t44;
t32 = sin(t40);
t33 = cos(t40);
t48 = cos(qJ(5));
t34 = t48 * pkin(5) + pkin(4);
t51 = -pkin(11) - pkin(10);
t93 = -t32 * t51 + t33 * t34;
t92 = t33 * pkin(4) + t32 * pkin(10);
t47 = sin(qJ(1));
t50 = cos(qJ(1));
t28 = g(1) * t50 + g(2) * t47;
t72 = t50 * t48;
t45 = sin(qJ(5));
t77 = t47 * t45;
t16 = t33 * t77 + t72;
t73 = t50 * t45;
t76 = t47 * t48;
t18 = -t33 * t73 + t76;
t82 = g(3) * t32;
t91 = -g(1) * t18 + g(2) * t16 + t45 * t82;
t7 = -g(3) * t33 + t28 * t32;
t52 = -pkin(8) - pkin(7);
t37 = sin(t44);
t90 = pkin(3) * t37;
t89 = pkin(4) * t32;
t88 = pkin(10) * t33;
t39 = cos(t44);
t31 = pkin(3) * t39;
t49 = cos(qJ(2));
t41 = t49 * pkin(2);
t70 = t31 + t41;
t24 = pkin(1) + t70;
t20 = t50 * t24;
t84 = g(2) * t20;
t43 = qJ(5) + qJ(6);
t36 = sin(t43);
t79 = t47 * t36;
t38 = cos(t43);
t78 = t47 * t38;
t75 = t50 * t36;
t74 = t50 * t38;
t68 = t31 + t92;
t42 = -pkin(9) + t52;
t67 = pkin(5) * t45 - t42;
t26 = t47 * t88;
t65 = -t47 * t89 + t26;
t27 = t50 * t88;
t64 = -t50 * t89 + t27;
t63 = t31 + t93;
t62 = -t89 - t90;
t60 = g(1) * t47 - g(2) * t50;
t58 = t32 * t34 + t33 * t51;
t57 = t58 * t47;
t56 = t58 * t50;
t9 = -g(3) * t39 + t28 * t37;
t46 = sin(qJ(2));
t53 = -g(3) * t49 + t28 * t46;
t35 = t41 + pkin(1);
t25 = -t46 * pkin(2) - t90;
t22 = t50 * t25;
t21 = t47 * t25;
t19 = t33 * t72 + t77;
t17 = -t33 * t76 + t73;
t15 = t60 * t32;
t14 = t33 * t74 + t79;
t13 = -t33 * t75 + t78;
t12 = -t33 * t78 + t75;
t11 = t33 * t79 + t74;
t10 = g(3) * t37 + t28 * t39;
t8 = t28 * t33 + t82;
t6 = t7 * t48;
t5 = t7 * t45;
t4 = t7 * t38;
t3 = t7 * t36;
t2 = g(1) * t14 - g(2) * t12 + t38 * t82;
t1 = -g(1) * t13 + g(2) * t11 + t36 * t82;
t23 = [0, 0, 0, 0, 0, 0, t60, t28, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t49, -t60 * t46, -t28, -g(1) * (-t47 * pkin(1) + t50 * pkin(7)) - g(2) * (t50 * pkin(1) + t47 * pkin(7)) 0, 0, 0, 0, 0, 0, t60 * t39, -t60 * t37, -t28, -g(1) * (-t47 * t35 - t50 * t52) - g(2) * (t50 * t35 - t47 * t52) 0, 0, 0, 0, 0, 0, t60 * t33, -t15, -t28, -g(1) * (-t47 * t24 - t50 * t42) - g(2) * (-t47 * t42 + t20) 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t19, -g(1) * t16 - g(2) * t18, t15, -t84 + (g(1) * t42 - g(2) * t92) * t50 + (-g(1) * (-t24 - t92) + g(2) * t42) * t47, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t15, -t84 + (-g(1) * t67 - g(2) * t93) * t50 + (-g(1) * (-t24 - t93) - g(2) * t67) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, g(3) * t46 + t28 * t49, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t53 * pkin(2), 0, 0, 0, 0, 0, 0, t7, t8, 0, -g(3) * t70 - t25 * t28, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t22 + t64) - g(2) * (t21 + t65) - g(3) * (t41 + t68) 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t22 - t56) - g(2) * (t21 - t57) - g(3) * (t41 + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t9 * pkin(3), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t50 * t62 + t27) - g(2) * (t47 * t62 + t26) - g(3) * t68, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t63 + t28 * (t58 + t90); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t64 - g(2) * t65 - g(3) * t92, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, g(1) * t56 + g(2) * t57 - g(3) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, g(1) * t19 - g(2) * t17 + t48 * t82, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t91 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t23;

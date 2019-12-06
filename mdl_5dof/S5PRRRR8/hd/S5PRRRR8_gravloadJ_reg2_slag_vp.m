% Calculate inertial parameters regressor of gravitation load for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t62 = cos(pkin(5));
t36 = sin(pkin(5));
t39 = sin(qJ(2));
t70 = t36 * t39;
t77 = -t38 * t70 + t62 * t41;
t42 = cos(qJ(2));
t35 = sin(pkin(10));
t55 = t35 * t62;
t61 = cos(pkin(10));
t24 = -t39 * t55 + t61 * t42;
t69 = t36 * t41;
t76 = -t24 * t38 + t35 * t69;
t75 = g(3) * t36;
t34 = qJ(3) + qJ(4);
t33 = cos(t34);
t37 = sin(qJ(5));
t73 = t33 * t37;
t40 = cos(qJ(5));
t72 = t33 * t40;
t71 = t35 * t36;
t68 = t36 * t42;
t67 = t37 * t42;
t43 = -pkin(8) - pkin(7);
t66 = t39 * t43;
t65 = t40 * t42;
t49 = t62 * t61;
t21 = t35 * t39 - t42 * t49;
t22 = t35 * t42 + t39 * t49;
t31 = t41 * pkin(3) + pkin(2);
t64 = -t21 * t31 - t22 * t43;
t23 = t61 * t39 + t42 * t55;
t63 = -t23 * t31 - t24 * t43;
t32 = sin(t34);
t54 = t36 * t61;
t10 = -t22 * t32 - t33 * t54;
t11 = t22 * t33 - t32 * t54;
t58 = t10 * pkin(4) + t11 * pkin(9);
t12 = -t24 * t32 + t33 * t71;
t13 = t24 * t33 + t32 * t71;
t57 = t12 * pkin(4) + t13 * pkin(9);
t17 = -t32 * t70 + t62 * t33;
t18 = t62 * t32 + t33 * t70;
t56 = t17 * pkin(4) + t18 * pkin(9);
t52 = t76 * pkin(3);
t51 = pkin(4) * t33 + pkin(9) * t32;
t50 = t77 * pkin(3);
t48 = g(1) * t12 + g(2) * t10 + g(3) * t17;
t5 = g(1) * t13 + g(2) * t11 + g(3) * t18;
t47 = -t22 * t38 - t41 * t54;
t46 = -g(1) * t23 - g(2) * t21 + g(3) * t68;
t45 = g(1) * t24 + g(2) * t22 + g(3) * t70;
t44 = t47 * pkin(3);
t25 = t31 * t68;
t6 = t46 * t32;
t2 = t48 * t40;
t1 = t48 * t37;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t45, 0, 0, 0, 0, 0, 0, 0, 0, -t46 * t41, t46 * t38, -t45, -g(1) * (-t23 * pkin(2) + t24 * pkin(7)) - g(2) * (-t21 * pkin(2) + t22 * pkin(7)) - (pkin(2) * t42 + pkin(7) * t39) * t75, 0, 0, 0, 0, 0, 0, -t46 * t33, t6, -t45, -g(1) * t63 - g(2) * t64 - g(3) * (-t36 * t66 + t25), 0, 0, 0, 0, 0, 0, -g(1) * (-t23 * t72 + t24 * t37) - g(2) * (-t21 * t72 + t22 * t37) - (t33 * t65 + t37 * t39) * t75, -g(1) * (t23 * t73 + t24 * t40) - g(2) * (t21 * t73 + t22 * t40) - (-t33 * t67 + t39 * t40) * t75, -t6, -g(1) * (-t51 * t23 + t63) - g(2) * (-t51 * t21 + t64) - g(3) * t25 - (t51 * t42 - t66) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t76 - g(2) * t47 - g(3) * t77, -g(1) * (-t24 * t41 - t38 * t71) - g(2) * (-t22 * t41 + t38 * t54) - g(3) * (-t62 * t38 - t39 * t69), 0, 0, 0, 0, 0, 0, 0, 0, -t48, t5, 0, -g(1) * t52 - g(2) * t44 - g(3) * t50, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (t52 + t57) - g(2) * (t44 + t58) - g(3) * (t50 + t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * t57 - g(2) * t58 - g(3) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t13 * t37 + t23 * t40) - g(2) * (-t11 * t37 + t21 * t40) - g(3) * (-t18 * t37 - t36 * t65), -g(1) * (-t13 * t40 - t23 * t37) - g(2) * (-t11 * t40 - t21 * t37) - g(3) * (-t18 * t40 + t36 * t67), 0, 0;];
taug_reg = t3;

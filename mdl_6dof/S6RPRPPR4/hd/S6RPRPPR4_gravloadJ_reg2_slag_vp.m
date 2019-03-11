% Calculate inertial parameters regressor of gravitation load for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t35 = pkin(9) + qJ(3);
t33 = cos(t35);
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t23 = g(1) * t44 + g(2) * t42;
t32 = sin(t35);
t77 = t23 * t32;
t7 = -g(3) * t33 + t77;
t26 = t32 * qJ(4);
t62 = t33 * pkin(3) + t26;
t40 = -pkin(7) - qJ(2);
t75 = g(2) * t40;
t73 = g(3) * t32;
t71 = t32 * t42;
t70 = t32 * t44;
t38 = cos(pkin(10));
t69 = t33 * t38;
t68 = t33 * t44;
t36 = sin(pkin(10));
t67 = t42 * t36;
t66 = t42 * t38;
t65 = t44 * t36;
t64 = t44 * t38;
t63 = t44 * t40;
t61 = qJ(4) * t33;
t60 = qJ(5) * t36;
t39 = cos(pkin(9));
t31 = pkin(2) * t39 + pkin(1);
t21 = t44 * t31;
t59 = pkin(3) * t68 + t26 * t44 + t21;
t58 = -pkin(3) - t60;
t57 = pkin(4) * t69 + t33 * t60 + t62;
t12 = t33 * t67 + t64;
t14 = t33 * t65 - t66;
t56 = g(1) * t12 - g(2) * t14;
t22 = g(1) * t42 - g(2) * t44;
t13 = t33 * t66 - t65;
t41 = sin(qJ(6));
t43 = cos(qJ(6));
t55 = t12 * t43 - t13 * t41;
t54 = t12 * t41 + t13 * t43;
t53 = t36 * t43 - t38 * t41;
t52 = t36 * t41 + t38 * t43;
t50 = -t31 - t62;
t49 = -pkin(4) * t13 - t12 * qJ(5) - t63;
t48 = g(3) * t53;
t15 = t33 * t64 + t67;
t47 = pkin(4) * t15 + qJ(5) * t14 + t59;
t45 = (-g(1) * t50 + t75) * t42;
t20 = t44 * t61;
t17 = t42 * t61;
t11 = g(1) * t71 - g(2) * t70;
t8 = t23 * t33 + t73;
t6 = t7 * t38;
t5 = t7 * t36;
t4 = g(1) * t13 - g(2) * t15;
t3 = t14 * t41 + t15 * t43;
t2 = t14 * t43 - t15 * t41;
t1 = -g(1) * t14 - g(2) * t12 - t36 * t73;
t9 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t39, -t22 * sin(pkin(9)) -t23, -g(1) * (-pkin(1) * t42 + qJ(2) * t44) - g(2) * (pkin(1) * t44 + qJ(2) * t42) 0, 0, 0, 0, 0, 0, t22 * t33, -t11, -t23, -g(1) * (-t42 * t31 - t63) - g(2) * (-t40 * t42 + t21) 0, 0, 0, 0, 0, 0, t4, -t56, t11, g(1) * t63 - g(2) * t59 + t45, 0, 0, 0, 0, 0, 0, t4, t11, t56, -g(1) * t49 - g(2) * t47 + t45, 0, 0, 0, 0, 0, 0, g(1) * t54 - g(2) * t3, g(1) * t55 - g(2) * t2, -t11, -g(1) * (-t13 * pkin(5) + t49) - g(2) * (pkin(5) * t15 - pkin(8) * t70 + t47) + (-g(1) * (pkin(8) * t32 + t50) + t75) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-pkin(3) * t70 + t20) - g(2) * (-pkin(3) * t71 + t17) - g(3) * t62, 0, 0, 0, 0, 0, 0, t6, -t8, t5, -g(1) * t20 - g(2) * t17 - g(3) * t57 + (pkin(4) * t38 - t58) * t77, 0, 0, 0, 0, 0, 0, t7 * t52, -t33 * t48 + t53 * t77, t8, -g(1) * (-pkin(8) * t68 + t20) - g(2) * (-t42 * t33 * pkin(8) + t17) - g(3) * (pkin(5) * t69 + t57) + (g(3) * pkin(8) + t23 * (-(-pkin(4) - pkin(5)) * t38 - t58)) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t55 - t32 * t48, g(1) * t3 + g(2) * t54 + t52 * t73, 0, 0;];
taug_reg  = t9;

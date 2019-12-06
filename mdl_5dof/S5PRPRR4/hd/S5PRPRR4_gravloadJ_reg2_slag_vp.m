% Calculate inertial parameters regressor of gravitation load for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t30 = sin(pkin(10));
t37 = sin(qJ(2));
t40 = cos(qJ(2));
t55 = cos(pkin(10));
t49 = -t37 * t30 + t40 * t55;
t31 = sin(pkin(9));
t64 = t31 * t37;
t32 = sin(pkin(5));
t36 = sin(qJ(4));
t63 = t32 * t36;
t39 = cos(qJ(4));
t62 = t32 * t39;
t61 = t32 * t40;
t34 = cos(pkin(5));
t60 = t34 * t37;
t59 = t34 * t40;
t35 = sin(qJ(5));
t58 = t35 * t39;
t38 = cos(qJ(5));
t56 = t38 * t39;
t33 = cos(pkin(9));
t54 = t33 * t59;
t20 = t49 * t32;
t24 = -t40 * t30 - t37 * t55;
t21 = t24 * t32;
t52 = pkin(2) * t61 + t20 * pkin(3) - t21 * pkin(7);
t51 = pkin(4) * t39 + pkin(8) * t36;
t22 = t24 * t34;
t10 = -t33 * t22 + t31 * t49;
t11 = -t31 * t22 - t33 * t49;
t50 = -t31 * t59 - t33 * t37;
t14 = t21 * t36 + t34 * t39;
t2 = -t10 * t36 - t33 * t62;
t4 = t11 * t36 + t31 * t62;
t48 = g(1) * t4 + g(2) * t2 + g(3) * t14;
t15 = -t21 * t39 + t34 * t36;
t3 = t10 * t39 - t33 * t63;
t5 = -t11 * t39 + t31 * t63;
t47 = g(1) * t5 + g(2) * t3 + g(3) * t15;
t46 = g(1) * t11 - g(2) * t10 + g(3) * t21;
t43 = t49 * t34;
t12 = t33 * t24 - t31 * t43;
t9 = t31 * t24 + t33 * t43;
t45 = g(1) * t12 + g(2) * t9 + g(3) * t20;
t25 = pkin(2) * t54;
t44 = -pkin(2) * t64 + t9 * pkin(3) + pkin(7) * t10 + t25;
t42 = t50 * pkin(2) + t12 * pkin(3) - t11 * pkin(7);
t41 = -g(1) * t50 - g(3) * t61;
t19 = -g(3) * t34 + (-g(1) * t31 + g(2) * t33) * t32;
t1 = t45 * t36;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * (t54 - t64) + t41, -g(1) * (t31 * t60 - t33 * t40) - g(2) * (-t31 * t40 - t33 * t60) + g(3) * t32 * t37, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t46, 0, -g(2) * t25 + (g(2) * t64 + t41) * pkin(2), 0, 0, 0, 0, 0, 0, -t45 * t39, t1, t46, -g(1) * t42 - g(2) * t44 - g(3) * t52, 0, 0, 0, 0, 0, 0, -g(1) * (-t11 * t35 + t12 * t56) - g(2) * (t10 * t35 + t9 * t56) - g(3) * (t20 * t56 - t21 * t35), -g(1) * (-t11 * t38 - t12 * t58) - g(2) * (t10 * t38 - t9 * t58) - g(3) * (-t20 * t58 - t21 * t38), -t1, -g(1) * (t51 * t12 + t42) - g(2) * (t51 * t9 + t44) - g(3) * (t51 * t20 + t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t38, t48 * t35, -t47, -g(1) * (t4 * pkin(4) + t5 * pkin(8)) - g(2) * (t2 * pkin(4) + t3 * pkin(8)) - g(3) * (t14 * pkin(4) + t15 * pkin(8)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t38 - t5 * t35) - g(2) * (-t3 * t35 - t9 * t38) - g(3) * (-t15 * t35 - t20 * t38), -g(1) * (t12 * t35 - t5 * t38) - g(2) * (-t3 * t38 + t9 * t35) - g(3) * (-t15 * t38 + t20 * t35), 0, 0;];
taug_reg = t6;

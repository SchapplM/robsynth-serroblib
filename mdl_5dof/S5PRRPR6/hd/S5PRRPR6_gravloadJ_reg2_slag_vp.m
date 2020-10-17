% Calculate inertial parameters regressor of gravitation load for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:47
% EndTime: 2019-12-05 16:32:49
% DurationCPUTime: 0.38s
% Computational Cost: add. (287->98), mult. (658->152), div. (0->0), fcn. (797->12), ass. (0->53)
t28 = sin(pkin(5));
t61 = g(3) * t28;
t26 = pkin(10) + qJ(5);
t24 = sin(t26);
t33 = cos(qJ(3));
t60 = t24 * t33;
t25 = cos(t26);
t59 = t25 * t33;
t27 = sin(pkin(10));
t32 = sin(qJ(2));
t58 = t27 * t32;
t57 = t27 * t33;
t56 = t28 * t32;
t34 = cos(qJ(2));
t55 = t28 * t34;
t29 = cos(pkin(10));
t54 = t29 * t33;
t53 = t33 * t34;
t52 = pkin(2) * t55 + pkin(7) * t56;
t51 = cos(pkin(5));
t50 = cos(pkin(9));
t49 = sin(pkin(9));
t48 = pkin(4) * t27 + pkin(7);
t40 = t51 * t50;
t11 = t32 * t40 + t49 * t34;
t10 = t49 * t32 - t34 * t40;
t8 = t10 * pkin(2);
t47 = t11 * pkin(7) - t8;
t39 = t51 * t49;
t13 = -t32 * t39 + t50 * t34;
t12 = t50 * t32 + t34 * t39;
t9 = t12 * pkin(2);
t46 = t13 * pkin(7) - t9;
t45 = g(3) * t52;
t44 = t28 * t50;
t43 = t28 * t49;
t31 = sin(qJ(3));
t42 = pkin(3) * t33 + qJ(4) * t31;
t23 = t29 * pkin(4) + pkin(3);
t30 = -pkin(8) - qJ(4);
t41 = t23 * t33 - t30 * t31;
t14 = t31 * t56 - t51 * t33;
t4 = t11 * t31 + t33 * t44;
t6 = t13 * t31 - t33 * t43;
t38 = g(1) * t6 + g(2) * t4 + g(3) * t14;
t15 = t51 * t31 + t33 * t56;
t5 = t11 * t33 - t31 * t44;
t7 = t13 * t33 + t31 * t43;
t37 = g(1) * t7 + g(2) * t5 + g(3) * t15;
t36 = -g(1) * t12 - g(2) * t10 + g(3) * t55;
t35 = g(1) * t13 + g(2) * t11 + g(3) * t56;
t3 = t36 * t31;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t35, 0, 0, 0, 0, 0, 0, 0, 0, -t36 * t33, t3, -t35, -g(1) * t46 - g(2) * t47 - t45, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t54 + t13 * t27) - g(2) * (-t10 * t54 + t11 * t27) - (t29 * t53 + t58) * t61, -g(1) * (t12 * t57 + t13 * t29) - g(2) * (t10 * t57 + t11 * t29) - (-t27 * t53 + t29 * t32) * t61, -t3, -g(1) * (-t42 * t12 + t46) - g(2) * (-t42 * t10 + t47) - g(3) * (t42 * t55 + t52), 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t59 + t13 * t24) - g(2) * (-t10 * t59 + t11 * t24) - (t24 * t32 + t25 * t53) * t61, -g(1) * (t12 * t60 + t13 * t25) - g(2) * (t10 * t60 + t11 * t25) - (-t24 * t53 + t25 * t32) * t61, -t3, -g(1) * (-t41 * t12 + t48 * t13 - t9) - g(2) * (-t41 * t10 + t48 * t11 - t8) - t45 - (pkin(4) * t58 + t41 * t34) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t37, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t29, -t38 * t27, -t37, -g(1) * (-t6 * pkin(3) + t7 * qJ(4)) - g(2) * (-t4 * pkin(3) + t5 * qJ(4)) - g(3) * (-t14 * pkin(3) + t15 * qJ(4)), 0, 0, 0, 0, 0, 0, t38 * t25, -t38 * t24, -t37, -g(1) * (-t6 * t23 - t7 * t30) - g(2) * (-t4 * t23 - t5 * t30) - g(3) * (-t14 * t23 - t15 * t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t25 - t7 * t24) - g(2) * (t10 * t25 - t5 * t24) - g(3) * (-t15 * t24 - t25 * t55), -g(1) * (-t12 * t24 - t7 * t25) - g(2) * (-t10 * t24 - t5 * t25) - g(3) * (-t15 * t25 + t24 * t55), 0, 0;];
taug_reg = t1;

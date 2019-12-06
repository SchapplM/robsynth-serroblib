% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t34 = cos(pkin(8));
t38 = cos(qJ(3));
t39 = cos(qJ(1));
t44 = t39 * t38;
t36 = sin(qJ(3));
t37 = sin(qJ(1));
t51 = t37 * t36;
t12 = t34 * t51 + t44;
t45 = t39 * t36;
t50 = t37 * t38;
t14 = t34 * t45 - t50;
t33 = sin(pkin(8));
t62 = g(1) * t33;
t63 = -g(2) * t12 + g(3) * t14 + t36 * t62;
t60 = g(2) * t39;
t29 = t39 * qJ(2);
t58 = g(3) * t29;
t57 = t36 * pkin(3);
t32 = qJ(3) + pkin(9);
t26 = sin(t32);
t18 = pkin(4) * t26 + t57;
t56 = t18 * t34;
t28 = qJ(5) + t32;
t23 = sin(t28);
t55 = t37 * t23;
t24 = cos(t28);
t54 = t37 * t24;
t53 = t37 * t26;
t27 = cos(t32);
t52 = t37 * t27;
t49 = t39 * t23;
t48 = t39 * t24;
t47 = t39 * t26;
t46 = t39 * t27;
t35 = -qJ(4) - pkin(6);
t30 = t38 * pkin(3);
t19 = pkin(4) * t27 + t30;
t21 = g(3) * t37 + t60;
t20 = g(2) * t37 - g(3) * t39;
t42 = pkin(2) * t34 + pkin(6) * t33 + pkin(1);
t41 = (pkin(2) + t19) * t34 - (-pkin(7) + t35) * t33 + pkin(1);
t40 = (t30 + pkin(2)) * t34 - t33 * t35 + pkin(1);
t16 = t21 * t33;
t15 = -t34 * t44 - t51;
t13 = t34 * t50 - t45;
t11 = g(1) * t34 + t20 * t33;
t10 = -t34 * t46 - t53;
t9 = t34 * t47 - t52;
t8 = t34 * t52 - t47;
t7 = t34 * t53 + t46;
t6 = -t34 * t48 - t55;
t5 = t34 * t49 - t54;
t4 = t34 * t54 - t49;
t3 = t34 * t55 + t48;
t2 = -g(2) * t4 - g(3) * t6 + t24 * t62;
t1 = -g(2) * t3 + g(3) * t5 + t23 * t62;
t17 = [0, 0, 0, 0, 0, 0, t21, -t20, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t34, -t16, t20, -g(2) * (-t39 * pkin(1) - t37 * qJ(2)) - g(3) * (-t37 * pkin(1) + t29), 0, 0, 0, 0, 0, 0, -g(2) * t15 + g(3) * t13, -g(2) * t14 - g(3) * t12, t16, -t58 + t42 * t60 + (g(2) * qJ(2) + g(3) * t42) * t37, 0, 0, 0, 0, 0, 0, -g(2) * t10 + g(3) * t8, -g(2) * t9 - g(3) * t7, t16, -t58 + (g(2) * t40 - g(3) * t57) * t39 + (-g(2) * (-qJ(2) - t57) + g(3) * t40) * t37, 0, 0, 0, 0, 0, 0, -g(2) * t6 + g(3) * t4, -g(2) * t5 - g(3) * t3, t16, -t58 + (g(2) * t41 - g(3) * t18) * t39 + (-g(2) * (-qJ(2) - t18) + g(3) * t41) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -g(2) * t13 - g(3) * t15 + t38 * t62, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t7 + g(3) * t9 + t26 * t62, -g(2) * t8 - g(3) * t10 + t27 * t62, 0, t63 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, t18 * t62 - g(2) * (t39 * t19 + t37 * t56) - g(3) * (t37 * t19 - t39 * t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t17;

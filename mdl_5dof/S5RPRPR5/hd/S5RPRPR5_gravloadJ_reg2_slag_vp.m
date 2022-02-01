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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% OptimizationMode: 2
% StartTime: 2022-01-23 09:26:00
% EndTime: 2022-01-23 09:26:01
% DurationCPUTime: 0.29s
% Computational Cost: add. (204->71), mult. (263->105), div. (0->0), fcn. (277->10), ass. (0->59)
t38 = qJ(3) + pkin(9);
t29 = cos(t38);
t43 = cos(qJ(3));
t34 = t43 * pkin(3);
t20 = pkin(4) * t29 + t34;
t39 = sin(pkin(8));
t40 = cos(pkin(8));
t47 = qJ(4) + pkin(6);
t69 = -(-pkin(7) - t47) * t39 + (pkin(2) + t20) * t40;
t44 = cos(qJ(1));
t48 = t44 * t43;
t41 = sin(qJ(3));
t42 = sin(qJ(1));
t56 = t42 * t41;
t13 = t40 * t56 + t48;
t49 = t44 * t41;
t55 = t42 * t43;
t15 = -t40 * t49 + t55;
t65 = g(3) * t39;
t68 = -g(1) * t15 + g(2) * t13 + t41 * t65;
t64 = t41 * pkin(3);
t63 = pkin(2) * t40 + pkin(1);
t30 = qJ(5) + t38;
t25 = sin(t30);
t60 = t42 * t25;
t26 = cos(t30);
t59 = t42 * t26;
t28 = sin(t38);
t58 = t42 * t28;
t57 = t42 * t29;
t19 = pkin(4) * t28 + t64;
t54 = t44 * t19;
t53 = t44 * t25;
t52 = t44 * t26;
t51 = t44 * t28;
t50 = t44 * t29;
t31 = t42 * qJ(2);
t46 = t44 * pkin(1) + t31;
t23 = g(1) * t44 + g(2) * t42;
t22 = g(1) * t42 - g(2) * t44;
t32 = t44 * qJ(2);
t27 = qJ(2) + t64;
t21 = t39 * pkin(6) + t63;
t17 = t22 * t39;
t16 = t40 * t48 + t56;
t14 = -t40 * t55 + t49;
t12 = g(3) * t40 - t23 * t39;
t11 = t40 * t34 + t47 * t39 + t63;
t10 = t40 * t50 + t58;
t9 = -t40 * t51 + t57;
t8 = -t40 * t57 + t51;
t7 = t40 * t58 + t50;
t6 = t40 * t52 + t60;
t5 = -t40 * t53 + t59;
t4 = -t40 * t59 + t53;
t3 = t40 * t60 + t52;
t2 = g(1) * t6 - g(2) * t4 + t26 * t65;
t1 = -g(1) * t5 + g(2) * t3 + t25 * t65;
t18 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t40, -t17, -t23, -g(1) * (-t42 * pkin(1) + t32) - g(2) * t46, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, t17, -g(1) * (-t21 * t42 + t32) - g(2) * (t21 * t44 + t31), 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t17, -g(1) * (-t11 * t42 + t27 * t44) - g(2) * (t11 * t44 + t27 * t42), 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t17, -g(1) * (t32 + t54) - g(2) * (t69 * t44 + t46) + (-g(1) * (-pkin(1) - t69) - g(2) * t19) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, g(1) * t16 - g(2) * t14 + t43 * t65, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t28 * t65, g(1) * t10 - g(2) * t8 + t29 * t65, 0, t68 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t42 * t20 - t40 * t54) - g(2) * (-t42 * t40 * t19 - t44 * t20) + t19 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t18;

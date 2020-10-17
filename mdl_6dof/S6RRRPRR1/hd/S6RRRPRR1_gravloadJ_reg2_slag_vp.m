% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:43:57
% EndTime: 2019-05-07 09:43:58
% DurationCPUTime: 0.36s
% Computational Cost: add. (542->89), mult. (399->114), div. (0->0), fcn. (373->12), ass. (0->62)
t41 = qJ(2) + qJ(3);
t34 = pkin(11) + t41;
t32 = qJ(5) + t34;
t27 = sin(t32);
t28 = cos(t32);
t68 = t28 * pkin(5) + t27 * pkin(10);
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t23 = g(1) * t47 + g(2) * t44;
t3 = -g(3) * t28 + t23 * t27;
t48 = -pkin(8) - pkin(7);
t35 = sin(t41);
t67 = pkin(3) * t35;
t66 = pkin(5) * t27;
t65 = pkin(10) * t28;
t64 = g(3) * t27;
t43 = sin(qJ(2));
t62 = t43 * pkin(2);
t42 = sin(qJ(6));
t61 = t44 * t42;
t45 = cos(qJ(6));
t60 = t44 * t45;
t59 = t47 * t42;
t58 = t47 * t45;
t30 = cos(t34);
t26 = pkin(4) * t30;
t36 = cos(t41);
t31 = pkin(3) * t36;
t56 = t26 + t31;
t46 = cos(qJ(2));
t38 = t46 * pkin(2);
t55 = t31 + t38;
t40 = -qJ(4) + t48;
t54 = t26 + t55;
t29 = sin(t34);
t17 = -pkin(4) * t29 - t67;
t16 = t17 - t62;
t53 = t16 - t66;
t52 = t17 - t66;
t51 = t56 + t68;
t22 = g(1) * t44 - g(2) * t47;
t7 = -g(3) * t36 + t23 * t35;
t49 = -g(3) * t46 + t23 * t43;
t37 = -pkin(9) + t40;
t33 = t38 + pkin(1);
t20 = t47 * t65;
t19 = t44 * t65;
t18 = pkin(1) + t55;
t15 = pkin(1) + t54;
t14 = t28 * t58 + t61;
t13 = -t28 * t59 + t60;
t12 = -t28 * t60 + t59;
t11 = t28 * t61 + t58;
t10 = t47 * t15;
t9 = t22 * t27;
t8 = g(3) * t35 + t23 * t36;
t6 = g(3) * t29 + t23 * t30;
t5 = -g(3) * t30 + t23 * t29;
t4 = t23 * t28 + t64;
t2 = t3 * t45;
t1 = t3 * t42;
t21 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t46, -t22 * t43, -t23, -g(1) * (-t44 * pkin(1) + t47 * pkin(7)) - g(2) * (t47 * pkin(1) + t44 * pkin(7)) 0, 0, 0, 0, 0, 0, t22 * t36, -t22 * t35, -t23, -g(1) * (-t44 * t33 - t47 * t48) - g(2) * (t47 * t33 - t44 * t48) 0, 0, 0, 0, 0, 0, t22 * t30, -t22 * t29, -t23, -g(1) * (-t44 * t18 - t47 * t40) - g(2) * (t47 * t18 - t44 * t40) 0, 0, 0, 0, 0, 0, t22 * t28, -t9, -t23, -g(1) * (-t44 * t15 - t47 * t37) - g(2) * (-t44 * t37 + t10) 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t9, -g(2) * t10 + (g(1) * t37 - g(2) * t68) * t47 + (-g(1) * (-t15 - t68) + g(2) * t37) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, g(3) * t43 + t23 * t46, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t49 * pkin(2), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(3) * t55 - t23 * (-t62 - t67) 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t54 - t23 * t16, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t53 * t47 + t20) - g(2) * (t53 * t44 + t19) - g(3) * (t38 + t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t56 - t23 * t17, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t52 * t47 + t20) - g(2) * (t52 * t44 + t19) - g(3) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t47 * t66 + t20) - g(2) * (-t44 * t66 + t19) - g(3) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 + g(2) * t11 + t42 * t64, g(1) * t14 - g(2) * t12 + t45 * t64, 0, 0;];
taug_reg  = t21;

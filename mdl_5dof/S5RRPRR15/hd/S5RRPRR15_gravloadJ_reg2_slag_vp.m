% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR15_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t15 = g(1) * t35 + g(2) * t32;
t31 = sin(qJ(2));
t66 = t15 * t31;
t23 = t31 * qJ(3);
t34 = cos(qJ(2));
t46 = t34 * pkin(2) + t23;
t30 = sin(qJ(4));
t33 = cos(qJ(4));
t50 = t32 * t33;
t11 = t30 * t35 + t31 * t50;
t55 = g(3) * t34;
t49 = t33 * t35;
t9 = -t30 * t32 + t31 * t49;
t65 = -g(1) * t9 - g(2) * t11 + t33 * t55;
t8 = g(3) * t31 + t15 * t34;
t62 = pkin(4) * t30;
t61 = pkin(7) * t34;
t60 = g(1) * t32;
t53 = t30 * t34;
t52 = t31 * t32;
t51 = t31 * t35;
t48 = t34 * t35;
t36 = -pkin(8) - pkin(7);
t47 = t34 * t36;
t45 = t35 * pkin(1) + t32 * pkin(6);
t44 = qJ(3) * t34;
t43 = pkin(4) * t53;
t41 = t30 * t51;
t40 = pkin(2) * t48 + t35 * t23 + t45;
t39 = -g(2) * t35 + t60;
t38 = -pkin(1) - t46;
t29 = qJ(4) + qJ(5);
t26 = t35 * pkin(6);
t22 = cos(t29);
t21 = sin(t29);
t20 = pkin(4) * t33 + pkin(3);
t18 = t35 * t44;
t16 = t32 * t44;
t14 = t39 * t34;
t13 = t39 * t31;
t12 = -t30 * t52 + t49;
t10 = t41 + t50;
t7 = -t55 + t66;
t6 = -t21 * t52 + t22 * t35;
t5 = t21 * t35 + t22 * t52;
t4 = t21 * t51 + t22 * t32;
t3 = -t21 * t32 + t22 * t51;
t2 = g(1) * t4 - g(2) * t6 - t21 * t55;
t1 = -g(1) * t3 - g(2) * t5 + t22 * t55;
t17 = [0, 0, 0, 0, 0, 0, t39, t15, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, -t15, -g(1) * (-t32 * pkin(1) + t26) - g(2) * t45, 0, 0, 0, 0, 0, 0, -t15, -t14, t13, -g(1) * t26 - g(2) * t40 - t38 * t60, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t14, -g(1) * (pkin(3) * t35 + t26) - g(2) * (pkin(7) * t48 + t40) + (-g(1) * (t38 - t61) - g(2) * pkin(3)) * t32, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t14, -g(1) * (t35 * t20 + t26) - g(2) * (pkin(4) * t41 - t35 * t47 + t40) + (-g(1) * (-t31 * t62 + t38 + t47) - g(2) * t20) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * (-pkin(2) * t51 + t18) - g(2) * (-pkin(2) * t52 + t16) - g(3) * t46, 0, 0, 0, 0, 0, 0, -t8 * t30, -t8 * t33, t7, -g(1) * t18 - g(2) * t16 - g(3) * (t46 + t61) + (pkin(2) + pkin(7)) * t66, 0, 0, 0, 0, 0, 0, -t8 * t21, -t8 * t22, t7, -g(1) * (t35 * t43 + t18) - g(2) * (t32 * t43 + t16) - g(3) * (t46 - t47) + (-g(3) * t62 + t15 * (pkin(2) - t36)) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, g(1) * t10 - g(2) * t12 - g(3) * t53, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t65 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t17;

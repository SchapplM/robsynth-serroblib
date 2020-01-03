% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t29 = qJ(2) + qJ(3);
t26 = qJ(4) + t29;
t21 = sin(t26);
t22 = cos(t26);
t54 = t22 * pkin(4) + t21 * pkin(9);
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t17 = g(1) * t35 + g(2) * t32;
t3 = -g(3) * t22 + t17 * t21;
t36 = -pkin(7) - pkin(6);
t24 = sin(t29);
t53 = pkin(3) * t24;
t52 = pkin(4) * t21;
t51 = pkin(9) * t22;
t50 = g(3) * t21;
t30 = sin(qJ(5));
t48 = t30 * t32;
t47 = t30 * t35;
t33 = cos(qJ(5));
t46 = t32 * t33;
t45 = t33 * t35;
t25 = cos(t29);
t20 = pkin(3) * t25;
t34 = cos(qJ(2));
t27 = t34 * pkin(2);
t43 = t20 + t27;
t42 = t20 + t54;
t31 = sin(qJ(2));
t14 = -pkin(2) * t31 - t53;
t41 = t14 - t52;
t40 = -t52 - t53;
t38 = g(1) * t32 - g(2) * t35;
t5 = -g(3) * t25 + t17 * t24;
t37 = -g(3) * t34 + t17 * t31;
t28 = -pkin(8) + t36;
t23 = t27 + pkin(1);
t16 = t35 * t51;
t15 = t32 * t51;
t13 = pkin(1) + t43;
t12 = t35 * t13;
t11 = t22 * t45 + t48;
t10 = -t22 * t47 + t46;
t9 = -t22 * t46 + t47;
t8 = t22 * t48 + t45;
t7 = t38 * t21;
t6 = g(3) * t24 + t17 * t25;
t4 = t17 * t22 + t50;
t2 = t3 * t33;
t1 = t3 * t30;
t18 = [0, 0, 0, 0, 0, 0, t38, t17, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t34, -t38 * t31, -t17, -g(1) * (-pkin(1) * t32 + pkin(6) * t35) - g(2) * (pkin(1) * t35 + pkin(6) * t32), 0, 0, 0, 0, 0, 0, t38 * t25, -t38 * t24, -t17, -g(1) * (-t32 * t23 - t35 * t36) - g(2) * (t35 * t23 - t32 * t36), 0, 0, 0, 0, 0, 0, t38 * t22, -t7, -t17, -g(1) * (-t13 * t32 - t28 * t35) - g(2) * (-t28 * t32 + t12), 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, t7, -g(2) * t12 + (g(1) * t28 - g(2) * t54) * t35 + (-g(1) * (-t13 - t54) + g(2) * t28) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, g(3) * t31 + t17 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t37 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t43 - t17 * t14, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t41 * t35 + t16) - g(2) * (t41 * t32 + t15) - g(3) * (t27 + t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t40 * t35 + t16) - g(2) * (t40 * t32 + t15) - g(3) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t35 * t52 + t16) - g(2) * (-t32 * t52 + t15) - g(3) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t8 + t30 * t50, g(1) * t11 - g(2) * t9 + t33 * t50, 0, 0;];
taug_reg = t18;

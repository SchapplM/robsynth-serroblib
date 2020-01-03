% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = qJ(2) + qJ(3);
t25 = pkin(9) + t31;
t21 = sin(t25);
t22 = cos(t25);
t40 = t22 * pkin(4) + t21 * pkin(8);
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t18 = g(1) * t37 + g(2) * t34;
t3 = -g(3) * t22 + t18 * t21;
t38 = -pkin(7) - pkin(6);
t26 = sin(t31);
t53 = pkin(3) * t26;
t52 = pkin(4) * t21;
t51 = pkin(8) * t22;
t50 = g(3) * t21;
t32 = sin(qJ(5));
t48 = t34 * t32;
t35 = cos(qJ(5));
t47 = t34 * t35;
t46 = t37 * t32;
t45 = t37 * t35;
t27 = cos(t31);
t23 = pkin(3) * t27;
t36 = cos(qJ(2));
t28 = t36 * pkin(2);
t44 = t23 + t28;
t43 = t23 + t40;
t33 = sin(qJ(2));
t14 = -t33 * pkin(2) - t53;
t42 = t14 - t52;
t41 = -t52 - t53;
t17 = g(1) * t34 - g(2) * t37;
t5 = -g(3) * t27 + t18 * t26;
t39 = -g(3) * t36 + t18 * t33;
t30 = -qJ(4) + t38;
t24 = t28 + pkin(1);
t16 = t37 * t51;
t15 = t34 * t51;
t13 = pkin(1) + t44;
t12 = t37 * t13;
t11 = t22 * t45 + t48;
t10 = -t22 * t46 + t47;
t9 = -t22 * t47 + t46;
t8 = t22 * t48 + t45;
t7 = t17 * t21;
t6 = g(3) * t26 + t18 * t27;
t4 = t18 * t22 + t50;
t2 = t3 * t35;
t1 = t3 * t32;
t19 = [0, 0, 0, 0, 0, 0, t17, t18, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t36, -t17 * t33, -t18, -g(1) * (-t34 * pkin(1) + t37 * pkin(6)) - g(2) * (t37 * pkin(1) + t34 * pkin(6)), 0, 0, 0, 0, 0, 0, t17 * t27, -t17 * t26, -t18, -g(1) * (-t34 * t24 - t37 * t38) - g(2) * (t37 * t24 - t34 * t38), 0, 0, 0, 0, 0, 0, t17 * t22, -t7, -t18, -g(1) * (-t34 * t13 - t37 * t30) - g(2) * (-t34 * t30 + t12), 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, t7, -g(2) * t12 + (g(1) * t30 - g(2) * t40) * t37 + (-g(1) * (-t13 - t40) + g(2) * t30) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, g(3) * t33 + t18 * t36, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t39 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t44 - t18 * t14, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t42 * t37 + t16) - g(2) * (t42 * t34 + t15) - g(3) * (t28 + t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t41 * t37 + t16) - g(2) * (t41 * t34 + t15) - g(3) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t8 + t32 * t50, g(1) * t11 - g(2) * t9 + t35 * t50, 0, 0;];
taug_reg = t19;

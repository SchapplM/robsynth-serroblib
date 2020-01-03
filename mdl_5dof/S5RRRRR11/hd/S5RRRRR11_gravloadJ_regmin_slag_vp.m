% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = sin(qJ(2));
t27 = sin(qJ(1));
t30 = cos(qJ(2));
t37 = cos(pkin(5));
t50 = cos(qJ(1));
t34 = t37 * t50;
t14 = t27 * t26 - t30 * t34;
t22 = qJ(4) + qJ(5);
t20 = sin(t22);
t21 = cos(t22);
t15 = t26 * t34 + t27 * t30;
t25 = sin(qJ(3));
t29 = cos(qJ(3));
t23 = sin(pkin(5));
t36 = t23 * t50;
t8 = t15 * t29 - t25 * t36;
t55 = -t14 * t21 + t8 * t20;
t54 = t14 * t20 + t8 * t21;
t24 = sin(qJ(4));
t28 = cos(qJ(4));
t53 = -t14 * t28 + t8 * t24;
t52 = t14 * t24 + t8 * t28;
t51 = g(3) * t23;
t45 = t20 * t29;
t44 = t21 * t29;
t43 = t23 * t26;
t42 = t23 * t29;
t41 = t23 * t30;
t40 = t24 * t29;
t39 = t28 * t29;
t38 = t29 * t30;
t35 = t27 * t37;
t17 = -t26 * t35 + t50 * t30;
t10 = -t17 * t25 + t27 * t42;
t32 = t15 * t25 + t29 * t36;
t33 = g(1) * t10 - g(2) * t32 + g(3) * (-t25 * t43 + t37 * t29);
t16 = t50 * t26 + t30 * t35;
t31 = -g(1) * t16 - g(2) * t14 + g(3) * t41;
t13 = t37 * t25 + t26 * t42;
t11 = t27 * t23 * t25 + t17 * t29;
t6 = t11 * t28 + t16 * t24;
t5 = -t11 * t24 + t16 * t28;
t4 = t11 * t21 + t16 * t20;
t3 = -t11 * t20 + t16 * t21;
t2 = g(1) * t4 + g(2) * t54 - g(3) * (-t13 * t21 + t20 * t41);
t1 = -g(1) * t3 + g(2) * t55 - g(3) * (-t13 * t20 - t21 * t41);
t7 = [0, g(1) * t27 - g(2) * t50, g(1) * t50 + g(2) * t27, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t17, -g(1) * t14 + g(2) * t16, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t11, -g(1) * t32 - g(2) * t10, 0, 0, 0, 0, 0, g(1) * t52 - g(2) * t6, -g(1) * t53 - g(2) * t5, 0, 0, 0, 0, 0, g(1) * t54 - g(2) * t4, -g(1) * t55 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, -t31, g(1) * t17 + g(2) * t15 + g(3) * t43, 0, 0, 0, 0, 0, -t31 * t29, t31 * t25, 0, 0, 0, 0, 0, -g(1) * (-t16 * t39 + t17 * t24) - g(2) * (-t14 * t39 + t15 * t24) - (t24 * t26 + t28 * t38) * t51, -g(1) * (t16 * t40 + t17 * t28) - g(2) * (t14 * t40 + t15 * t28) - (-t24 * t38 + t26 * t28) * t51, 0, 0, 0, 0, 0, -g(1) * (-t16 * t44 + t17 * t20) - g(2) * (-t14 * t44 + t15 * t20) - (t20 * t26 + t21 * t38) * t51, -g(1) * (t16 * t45 + t17 * t21) - g(2) * (t14 * t45 + t15 * t21) - (-t20 * t38 + t21 * t26) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, g(1) * t11 + g(2) * t8 + g(3) * t13, 0, 0, 0, 0, 0, -t33 * t28, t33 * t24, 0, 0, 0, 0, 0, -t33 * t21, t33 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t53 - g(3) * (-t13 * t24 - t28 * t41), g(1) * t6 + g(2) * t52 - g(3) * (-t13 * t28 + t24 * t41), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t7;

% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR7
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:51
% EndTime: 2019-12-31 18:19:51
% DurationCPUTime: 0.22s
% Computational Cost: add. (200->54), mult. (179->70), div. (0->0), fcn. (173->10), ass. (0->37)
t18 = qJ(3) + pkin(9);
t12 = sin(t18);
t14 = cos(t18);
t45 = t14 * pkin(4) + t12 * pkin(7);
t19 = qJ(1) + pkin(8);
t13 = sin(t19);
t15 = cos(t19);
t8 = g(1) * t15 + g(2) * t13;
t28 = -g(3) * t14 + t8 * t12;
t42 = g(3) * t12;
t23 = sin(qJ(1));
t38 = t23 * pkin(1);
t25 = cos(qJ(3));
t16 = t25 * pkin(3);
t11 = t16 + pkin(2);
t26 = cos(qJ(1));
t17 = t26 * pkin(1);
t37 = t15 * t11 + t17;
t21 = sin(qJ(5));
t36 = t13 * t21;
t24 = cos(qJ(5));
t35 = t13 * t24;
t34 = t15 * t21;
t33 = t15 * t24;
t7 = g(1) * t13 - g(2) * t15;
t31 = g(1) * t23 - g(2) * t26;
t20 = -qJ(4) - pkin(6);
t30 = -t15 * t20 - t38;
t22 = sin(qJ(3));
t27 = -g(3) * t25 + t8 * t22;
t6 = t14 * t33 + t36;
t5 = -t14 * t34 + t35;
t4 = -t14 * t35 + t34;
t3 = t14 * t36 + t33;
t2 = t7 * t12;
t1 = t8 * t14 + t42;
t9 = [0, 0, 0, 0, 0, 0, t31, g(1) * t26 + g(2) * t23, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t31 * pkin(1), 0, 0, 0, 0, 0, 0, t7 * t25, -t7 * t22, -t8, -g(1) * (-t13 * pkin(2) + t15 * pkin(6) - t38) - g(2) * (t15 * pkin(2) + t13 * pkin(6) + t17), 0, 0, 0, 0, 0, 0, t7 * t14, -t2, -t8, -g(1) * (-t13 * t11 + t30) - g(2) * (-t13 * t20 + t37), 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t2, -g(1) * t30 - g(2) * (t45 * t15 + t37) + (-g(1) * (-t11 - t45) + g(2) * t20) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, g(3) * t22 + t8 * t25, 0, 0, 0, 0, 0, 0, 0, 0, t28, t1, 0, t27 * pkin(3), 0, 0, 0, 0, 0, 0, t28 * t24, -t28 * t21, -t1, -g(3) * (t16 + t45) + t8 * (pkin(3) * t22 + pkin(4) * t12 - pkin(7) * t14); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t21 * t42, g(1) * t6 - g(2) * t4 + t24 * t42, 0, 0;];
taug_reg = t9;

% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:10:46
% EndTime: 2019-05-05 15:10:47
% DurationCPUTime: 0.28s
% Computational Cost: add. (369->76), mult. (251->91), div. (0->0), fcn. (240->12), ass. (0->49)
t31 = pkin(11) + qJ(4);
t27 = qJ(5) + t31;
t20 = sin(t27);
t21 = cos(t27);
t44 = t21 * pkin(5) + t20 * pkin(9);
t32 = qJ(1) + pkin(10);
t24 = sin(t32);
t26 = cos(t32);
t12 = g(1) * t26 + g(2) * t24;
t3 = -g(3) * t21 + t12 * t20;
t55 = pkin(5) * t20;
t54 = g(3) * t20;
t37 = sin(qJ(1));
t52 = t37 * pkin(1);
t34 = cos(pkin(11));
t22 = t34 * pkin(3) + pkin(2);
t51 = t20 * t26;
t50 = t21 * t26;
t36 = sin(qJ(6));
t49 = t24 * t36;
t38 = cos(qJ(6));
t48 = t24 * t38;
t47 = t26 * t36;
t46 = t26 * t38;
t35 = -pkin(7) - qJ(3);
t25 = cos(t31);
t19 = pkin(4) * t25;
t13 = t19 + t22;
t39 = cos(qJ(1));
t29 = t39 * pkin(1);
t45 = t26 * t13 + t29;
t23 = sin(t31);
t43 = -pkin(4) * t23 - t55;
t11 = g(1) * t24 - g(2) * t26;
t42 = g(1) * t37 - g(2) * t39;
t30 = -pkin(8) + t35;
t41 = -t26 * t30 - t52;
t40 = -g(3) * t25 + t12 * t23;
t15 = pkin(9) * t50;
t14 = t24 * t21 * pkin(9);
t9 = t21 * t46 + t49;
t8 = -t21 * t47 + t48;
t7 = -t21 * t48 + t47;
t6 = t21 * t49 + t46;
t5 = t11 * t20;
t4 = t12 * t21 + t54;
t2 = t3 * t38;
t1 = t3 * t36;
t10 = [0, 0, 0, 0, 0, 0, t42, g(1) * t39 + g(2) * t37, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, t42 * pkin(1), 0, 0, 0, 0, 0, 0, t11 * t34, -t11 * sin(pkin(11)) -t12, -g(1) * (-t24 * pkin(2) + t26 * qJ(3) - t52) - g(2) * (t26 * pkin(2) + t24 * qJ(3) + t29) 0, 0, 0, 0, 0, 0, t11 * t25, -t11 * t23, -t12, -g(1) * (-t24 * t22 - t26 * t35 - t52) - g(2) * (t26 * t22 - t24 * t35 + t29) 0, 0, 0, 0, 0, 0, t11 * t21, -t5, -t12, -g(1) * (-t24 * t13 + t41) - g(2) * (-t24 * t30 + t45) 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t5, -g(1) * t41 - g(2) * (pkin(5) * t50 + pkin(9) * t51 + t45) + (-g(1) * (-t13 - t44) + g(2) * t30) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, g(3) * t23 + t12 * t25, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t40 * pkin(4), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t43 * t26 + t15) - g(2) * (t43 * t24 + t14) - g(3) * (t19 + t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-pkin(5) * t51 + t15) - g(2) * (-t24 * t55 + t14) - g(3) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t6 + t36 * t54, g(1) * t9 - g(2) * t7 + t38 * t54, 0, 0;];
taug_reg  = t10;

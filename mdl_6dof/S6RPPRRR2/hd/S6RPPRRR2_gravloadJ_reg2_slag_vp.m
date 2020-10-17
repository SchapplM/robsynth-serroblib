% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRR2
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
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:21:25
% EndTime: 2019-05-05 15:21:26
% DurationCPUTime: 0.35s
% Computational Cost: add. (377->78), mult. (296->101), div. (0->0), fcn. (297->12), ass. (0->52)
t27 = qJ(1) + pkin(10);
t20 = sin(t27);
t22 = cos(t27);
t14 = g(1) * t22 + g(2) * t20;
t26 = pkin(11) + qJ(4);
t19 = sin(t26);
t21 = cos(t26);
t44 = pkin(4) * t21 + pkin(8) * t19;
t32 = sin(qJ(5));
t49 = t22 * t32;
t34 = cos(qJ(5));
t52 = t20 * t34;
t11 = -t21 * t49 + t52;
t58 = g(3) * t19;
t48 = t22 * t34;
t53 = t20 * t32;
t9 = t21 * t53 + t48;
t65 = -g(1) * t11 + g(2) * t9 + t32 * t58;
t38 = -g(3) * t21 + t14 * t19;
t33 = sin(qJ(1));
t56 = t33 * pkin(1);
t28 = qJ(5) + qJ(6);
t23 = sin(t28);
t55 = t20 * t23;
t24 = cos(t28);
t54 = t20 * t24;
t51 = t22 * t23;
t50 = t22 * t24;
t30 = cos(pkin(11));
t17 = pkin(3) * t30 + pkin(2);
t35 = cos(qJ(1));
t25 = t35 * pkin(1);
t47 = t22 * t17 + t25;
t31 = -pkin(7) - qJ(3);
t45 = pkin(5) * t32 - t31;
t13 = g(1) * t20 - g(2) * t22;
t42 = g(1) * t33 - g(2) * t35;
t41 = -t22 * t31 - t56;
t18 = pkin(5) * t34 + pkin(4);
t36 = -pkin(9) - pkin(8);
t40 = t21 * t18 - t19 * t36;
t12 = t21 * t48 + t53;
t10 = -t21 * t52 + t49;
t8 = t13 * t19;
t7 = t21 * t50 + t55;
t6 = -t21 * t51 + t54;
t5 = -t21 * t54 + t51;
t4 = t21 * t55 + t50;
t3 = t14 * t21 + t58;
t2 = g(1) * t7 - g(2) * t5 + t24 * t58;
t1 = -g(1) * t6 + g(2) * t4 + t23 * t58;
t15 = [0, 0, 0, 0, 0, 0, t42, g(1) * t35 + g(2) * t33, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, t42 * pkin(1), 0, 0, 0, 0, 0, 0, t13 * t30, -t13 * sin(pkin(11)) -t14, -g(1) * (-pkin(2) * t20 + qJ(3) * t22 - t56) - g(2) * (pkin(2) * t22 + qJ(3) * t20 + t25) 0, 0, 0, 0, 0, 0, t13 * t21, -t8, -t14, -g(1) * (-t20 * t17 + t41) - g(2) * (-t20 * t31 + t47) 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t8, -g(1) * t41 - g(2) * (t44 * t22 + t47) + (-g(1) * (-t17 - t44) + g(2) * t31) * t20, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, -g(1) * t4 - g(2) * t6, t8, g(1) * t56 - g(2) * t47 + (-g(1) * t45 - g(2) * t40) * t22 + (-g(1) * (-t17 - t40) - g(2) * t45) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t3, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t34, -t38 * t32, -t3, -g(3) * t44 + t14 * (pkin(4) * t19 - pkin(8) * t21) 0, 0, 0, 0, 0, 0, t38 * t24, -t38 * t23, -t3, -g(3) * t40 + t14 * (t18 * t19 + t21 * t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, g(1) * t12 - g(2) * t10 + t34 * t58, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t65 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t15;

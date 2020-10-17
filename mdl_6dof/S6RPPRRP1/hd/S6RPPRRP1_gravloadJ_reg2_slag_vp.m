% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:44:53
% EndTime: 2019-05-05 14:44:53
% DurationCPUTime: 0.30s
% Computational Cost: add. (328->68), mult. (281->81), div. (0->0), fcn. (278->10), ass. (0->44)
t25 = pkin(10) + qJ(4);
t20 = sin(t25);
t22 = cos(t25);
t40 = t22 * pkin(4) + t20 * pkin(8);
t26 = qJ(1) + pkin(9);
t21 = sin(t26);
t23 = cos(t26);
t15 = g(1) * t23 + g(2) * t21;
t33 = cos(qJ(5));
t44 = t23 * t33;
t31 = sin(qJ(5));
t47 = t21 * t31;
t10 = t22 * t47 + t44;
t45 = t23 * t31;
t46 = t21 * t33;
t12 = -t22 * t45 + t46;
t52 = g(3) * t20;
t1 = -g(1) * t12 + g(2) * t10 + t31 * t52;
t7 = -g(3) * t22 + t15 * t20;
t32 = sin(qJ(1));
t48 = t32 * pkin(1);
t28 = cos(pkin(10));
t18 = t28 * pkin(3) + pkin(2);
t34 = cos(qJ(1));
t24 = t34 * pkin(1);
t43 = t23 * t18 + t24;
t30 = -pkin(7) - qJ(3);
t41 = pkin(5) * t31 - t30;
t14 = g(1) * t21 - g(2) * t23;
t38 = g(1) * t32 - g(2) * t34;
t37 = -t23 * t30 - t48;
t19 = t33 * pkin(5) + pkin(4);
t29 = -qJ(6) - pkin(8);
t36 = t22 * t19 - t20 * t29;
t13 = t22 * t44 + t47;
t11 = -t22 * t46 + t45;
t9 = t14 * t20;
t8 = t15 * t22 + t52;
t6 = t7 * t33;
t5 = t7 * t31;
t4 = -g(1) * t11 - g(2) * t13;
t3 = -g(1) * t10 - g(2) * t12;
t2 = g(1) * t13 - g(2) * t11 + t33 * t52;
t16 = [0, 0, 0, 0, 0, 0, t38, g(1) * t34 + g(2) * t32, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, 0, t38 * pkin(1), 0, 0, 0, 0, 0, 0, t14 * t28, -t14 * sin(pkin(10)) -t15, -g(1) * (-t21 * pkin(2) + t23 * qJ(3) - t48) - g(2) * (t23 * pkin(2) + t21 * qJ(3) + t24) 0, 0, 0, 0, 0, 0, t14 * t22, -t9, -t15, -g(1) * (-t21 * t18 + t37) - g(2) * (-t21 * t30 + t43) 0, 0, 0, 0, 0, 0, t4, t3, t9, -g(1) * t37 - g(2) * (t40 * t23 + t43) + (-g(1) * (-t18 - t40) + g(2) * t30) * t21, 0, 0, 0, 0, 0, 0, t4, t3, t9, g(1) * t48 - g(2) * t43 + (-g(1) * t41 - g(2) * t36) * t23 + (-g(1) * (-t18 - t36) - g(2) * t41) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t40 + t15 * (pkin(4) * t20 - pkin(8) * t22) 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t36 + t15 * (t19 * t20 + t22 * t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg  = t16;

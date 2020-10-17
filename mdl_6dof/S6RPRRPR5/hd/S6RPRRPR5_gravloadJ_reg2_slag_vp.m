% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:40:07
% EndTime: 2019-05-05 22:40:08
% DurationCPUTime: 0.34s
% Computational Cost: add. (400->87), mult. (342->98), div. (0->0), fcn. (325->10), ass. (0->51)
t33 = pkin(10) + qJ(3);
t29 = qJ(4) + t33;
t24 = sin(t29);
t25 = cos(t29);
t49 = t25 * pkin(4) + t24 * qJ(5);
t38 = sin(qJ(1));
t40 = cos(qJ(1));
t18 = g(1) * t40 + g(2) * t38;
t63 = t18 * t24;
t4 = g(3) * t24 + t18 * t25;
t27 = sin(t33);
t61 = pkin(3) * t27;
t60 = pkin(4) * t24;
t56 = g(3) * t25;
t20 = t25 * pkin(9);
t36 = -pkin(7) - qJ(2);
t32 = -pkin(8) + t36;
t55 = pkin(5) - t32;
t35 = cos(pkin(10));
t26 = t35 * pkin(2) + pkin(1);
t37 = sin(qJ(6));
t54 = t38 * t37;
t39 = cos(qJ(6));
t53 = t38 * t39;
t52 = t40 * t32;
t51 = t40 * t37;
t50 = t40 * t39;
t48 = qJ(5) * t25;
t28 = cos(t33);
t23 = pkin(3) * t28;
t47 = t23 + t49;
t12 = t23 + t26;
t11 = t40 * t12;
t46 = g(2) * (t49 * t40 + t11);
t45 = -t60 - t61;
t17 = g(1) * t38 - g(2) * t40;
t44 = -t12 - t49;
t42 = -g(3) * t28 + t18 * t27;
t41 = (pkin(4) + pkin(9)) * t63;
t15 = t40 * t48;
t13 = t38 * t48;
t10 = -t24 * t54 + t50;
t9 = t24 * t53 + t51;
t8 = t24 * t51 + t53;
t7 = t24 * t50 - t54;
t6 = t17 * t25;
t5 = t17 * t24;
t3 = -t56 + t63;
t2 = t4 * t39;
t1 = t4 * t37;
t14 = [0, 0, 0, 0, 0, 0, t17, t18, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t35, -t17 * sin(pkin(10)) -t18, -g(1) * (-t38 * pkin(1) + t40 * qJ(2)) - g(2) * (t40 * pkin(1) + t38 * qJ(2)) 0, 0, 0, 0, 0, 0, t17 * t28, -t17 * t27, -t18, -g(1) * (-t38 * t26 - t40 * t36) - g(2) * (t40 * t26 - t38 * t36) 0, 0, 0, 0, 0, 0, t6, -t5, -t18, -g(1) * (-t38 * t12 - t52) - g(2) * (-t38 * t32 + t11) 0, 0, 0, 0, 0, 0, -t18, -t6, t5, g(1) * t52 - t46 + (-g(1) * t44 + g(2) * t32) * t38, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, t6, -t46 + (-g(1) * t55 - g(2) * t20) * t40 + (-g(1) * (t44 - t20) - g(2) * t55) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, g(3) * t27 + t18 * t28, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t42 * pkin(3), 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (t45 * t40 + t15) - g(2) * (t45 * t38 + t13) - g(3) * t47, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * (-t40 * t61 + t15) - g(2) * (-t38 * t61 + t13) - g(3) * (t20 + t47) + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (-t40 * t60 + t15) - g(2) * (-t38 * t60 + t13) - g(3) * t49, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * t15 - g(2) * t13 - g(3) * (t20 + t49) + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t39 * t56, g(1) * t8 - g(2) * t10 - t37 * t56, 0, 0;];
taug_reg  = t14;

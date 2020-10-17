% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:15:12
% EndTime: 2019-05-06 01:15:13
% DurationCPUTime: 0.42s
% Computational Cost: add. (443->82), mult. (430->117), div. (0->0), fcn. (434->10), ass. (0->55)
t39 = cos(qJ(4));
t31 = t39 * pkin(4);
t26 = t31 + pkin(3);
t37 = sin(qJ(3));
t40 = cos(qJ(3));
t42 = -pkin(9) - pkin(8);
t45 = t40 * t26 - t37 * t42;
t35 = qJ(4) + qJ(5);
t30 = cos(t35);
t21 = pkin(5) * t30 + t31;
t19 = pkin(3) + t21;
t33 = -qJ(6) + t42;
t47 = t40 * t19 - t37 * t33;
t34 = qJ(1) + pkin(10);
t27 = sin(t34);
t28 = cos(t34);
t18 = g(1) * t28 + g(2) * t27;
t36 = sin(qJ(4));
t61 = t36 * t40;
t13 = t27 * t61 + t28 * t39;
t15 = t27 * t39 - t28 * t61;
t68 = g(3) * t37;
t76 = -g(1) * t15 + g(2) * t13 + t36 * t68;
t29 = sin(t35);
t63 = t29 * t40;
t7 = t27 * t63 + t28 * t30;
t9 = t27 * t30 - t28 * t63;
t1 = -g(1) * t9 + g(2) * t7 + t29 * t68;
t11 = -g(3) * t40 + t18 * t37;
t72 = g(1) * t27;
t66 = t36 * pkin(4);
t20 = pkin(5) * t29 + t66;
t65 = t28 * t20;
t64 = t28 * t36;
t62 = t30 * t40;
t58 = t39 * t40;
t41 = cos(qJ(1));
t53 = t41 * pkin(1) + t28 * pkin(2) + t27 * pkin(7);
t38 = sin(qJ(1));
t52 = -t38 * pkin(1) + t28 * pkin(7);
t51 = t40 * pkin(3) + t37 * pkin(8);
t49 = -g(2) * t28 + t72;
t48 = g(1) * t38 - g(2) * t41;
t17 = t49 * t37;
t16 = t27 * t36 + t28 * t58;
t14 = -t27 * t58 + t64;
t12 = t18 * t40 + t68;
t10 = t27 * t29 + t28 * t62;
t8 = -t27 * t62 + t28 * t29;
t6 = t11 * t30;
t5 = t11 * t29;
t4 = -g(1) * t8 - g(2) * t10;
t3 = -g(1) * t7 - g(2) * t9;
t2 = g(1) * t10 - g(2) * t8 + t30 * t68;
t22 = [0, 0, 0, 0, 0, 0, t48, g(1) * t41 + g(2) * t38, 0, 0, 0, 0, 0, 0, 0, 0, t49, t18, 0, t48 * pkin(1), 0, 0, 0, 0, 0, 0, t49 * t40, -t17, -t18, -g(1) * (-t27 * pkin(2) + t52) - g(2) * t53, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, t17, -g(1) * t52 - g(2) * (t51 * t28 + t53) - (-pkin(2) - t51) * t72, 0, 0, 0, 0, 0, 0, t4, t3, t17, -g(1) * (pkin(4) * t64 + t52) - g(2) * (t45 * t28 + t53) + (-g(1) * (-pkin(2) - t45) - g(2) * t66) * t27, 0, 0, 0, 0, 0, 0, t4, t3, t17, -g(1) * (t52 + t65) - g(2) * (t47 * t28 + t53) + (-g(1) * (-pkin(2) - t47) - g(2) * t20) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t39, -t11 * t36, -t12, -g(3) * t51 + t18 * (pkin(3) * t37 - pkin(8) * t40) 0, 0, 0, 0, 0, 0, t6, -t5, -t12, -g(3) * t45 + t18 * (t26 * t37 + t40 * t42) 0, 0, 0, 0, 0, 0, t6, -t5, -t12, -g(3) * t47 + t18 * (t19 * t37 + t33 * t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, g(1) * t16 - g(2) * t14 + t39 * t68, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t76 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t27 * t21 - t40 * t65) - g(2) * (-t27 * t40 * t20 - t28 * t21) + t20 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11;];
taug_reg  = t22;

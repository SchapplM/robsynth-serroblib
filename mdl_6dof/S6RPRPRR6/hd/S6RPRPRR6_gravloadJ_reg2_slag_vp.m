% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:03:22
% EndTime: 2019-05-05 19:03:24
% DurationCPUTime: 0.42s
% Computational Cost: add. (403->89), mult. (383->121), div. (0->0), fcn. (385->12), ass. (0->57)
t34 = sin(pkin(11));
t39 = -pkin(7) - qJ(2);
t49 = t34 * pkin(4) - t39;
t40 = sin(qJ(1));
t41 = cos(qJ(1));
t18 = g(1) * t41 + g(2) * t40;
t33 = pkin(10) + qJ(3);
t27 = cos(t33);
t32 = pkin(11) + qJ(5);
t24 = sin(t32);
t55 = t41 * t24;
t26 = cos(t32);
t60 = t40 * t26;
t11 = -t27 * t55 + t60;
t25 = sin(t33);
t66 = g(3) * t25;
t54 = t41 * t26;
t61 = t40 * t24;
t9 = t27 * t61 + t54;
t72 = -g(1) * t11 + g(2) * t9 + t24 * t66;
t3 = -g(3) * t27 + t18 * t25;
t37 = cos(pkin(10));
t23 = t37 * pkin(2) + pkin(1);
t16 = t41 * t23;
t68 = g(2) * t16;
t36 = cos(pkin(11));
t22 = t36 * pkin(4) + pkin(3);
t28 = qJ(6) + t32;
t20 = sin(t28);
t63 = t40 * t20;
t21 = cos(t28);
t62 = t40 * t21;
t59 = t40 * t34;
t58 = t40 * t36;
t57 = t41 * t20;
t56 = t41 * t21;
t53 = t41 * t34;
t52 = t41 * t36;
t38 = -pkin(8) - qJ(4);
t51 = pkin(5) * t24 + t49;
t17 = g(1) * t40 - g(2) * t41;
t48 = t27 * pkin(3) + t25 * qJ(4);
t14 = pkin(5) * t26 + t22;
t31 = -pkin(9) + t38;
t46 = t27 * t14 - t25 * t31;
t44 = t27 * t22 - t25 * t38;
t13 = t17 * t25;
t12 = t27 * t54 + t61;
t10 = -t27 * t60 + t55;
t8 = t27 * t56 + t63;
t7 = -t27 * t57 + t62;
t6 = -t27 * t62 + t57;
t5 = t27 * t63 + t56;
t4 = t18 * t27 + t66;
t2 = g(1) * t8 - g(2) * t6 + t21 * t66;
t1 = -g(1) * t7 + g(2) * t5 + t20 * t66;
t15 = [0, 0, 0, 0, 0, 0, t17, t18, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t37, -t17 * sin(pkin(10)) -t18, -g(1) * (-t40 * pkin(1) + t41 * qJ(2)) - g(2) * (t41 * pkin(1) + t40 * qJ(2)) 0, 0, 0, 0, 0, 0, t17 * t27, -t13, -t18, -g(1) * (-t40 * t23 - t41 * t39) - g(2) * (-t40 * t39 + t16) 0, 0, 0, 0, 0, 0, -g(1) * (-t27 * t58 + t53) - g(2) * (t27 * t52 + t59) -g(1) * (t27 * t59 + t52) - g(2) * (-t27 * t53 + t58) t13, -t68 + (g(1) * t39 - g(2) * t48) * t41 + (-g(1) * (-t23 - t48) + g(2) * t39) * t40, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, -t68 + (-g(1) * t49 - g(2) * t44) * t41 + (-g(1) * (-t23 - t44) - g(2) * t49) * t40, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, t13, -t68 + (-g(1) * t51 - g(2) * t46) * t41 + (-g(1) * (-t23 - t46) - g(2) * t51) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t36, -t3 * t34, -t4, -g(3) * t48 + t18 * (pkin(3) * t25 - qJ(4) * t27) 0, 0, 0, 0, 0, 0, t3 * t26, -t3 * t24, -t4, -g(3) * t44 + t18 * (t22 * t25 + t27 * t38) 0, 0, 0, 0, 0, 0, t3 * t21, -t3 * t20, -t4, -g(3) * t46 + t18 * (t14 * t25 + t27 * t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, g(1) * t12 - g(2) * t10 + t26 * t66, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t72 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t15;

% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:32:01
% EndTime: 2019-05-06 19:32:02
% DurationCPUTime: 0.35s
% Computational Cost: add. (516->87), mult. (377->112), div. (0->0), fcn. (354->12), ass. (0->56)
t38 = qJ(2) + pkin(11);
t33 = qJ(4) + t38;
t29 = qJ(5) + t33;
t24 = sin(t29);
t25 = cos(t29);
t63 = t25 * pkin(5) + t24 * pkin(10);
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t20 = g(1) * t45 + g(2) * t42;
t3 = -g(3) * t25 + t20 * t24;
t27 = sin(t33);
t62 = pkin(4) * t27;
t61 = pkin(5) * t24;
t60 = pkin(10) * t25;
t59 = g(3) * t24;
t40 = sin(qJ(6));
t57 = t42 * t40;
t43 = cos(qJ(6));
t56 = t42 * t43;
t55 = t45 * t40;
t54 = t45 * t43;
t39 = -qJ(3) - pkin(7);
t32 = cos(t38);
t44 = cos(qJ(2));
t35 = t44 * pkin(2);
t52 = pkin(3) * t32 + t35;
t37 = -pkin(8) + t39;
t28 = cos(t33);
t23 = pkin(4) * t28;
t51 = t23 + t52;
t31 = sin(t38);
t41 = sin(qJ(2));
t49 = -t41 * pkin(2) - pkin(3) * t31;
t14 = t49 - t62;
t50 = t14 - t61;
t48 = -t61 - t62;
t19 = g(1) * t42 - g(2) * t45;
t5 = -g(3) * t28 + t20 * t27;
t46 = -g(3) * t44 + t20 * t41;
t34 = -pkin(9) + t37;
t30 = t35 + pkin(1);
t17 = t45 * t60;
t16 = t42 * t60;
t15 = pkin(1) + t52;
t13 = t25 * t54 + t57;
t12 = -t25 * t55 + t56;
t11 = -t25 * t56 + t55;
t10 = t25 * t57 + t54;
t9 = pkin(1) + t51;
t8 = t45 * t9;
t7 = t19 * t24;
t6 = g(3) * t27 + t20 * t28;
t4 = t20 * t25 + t59;
t2 = t3 * t43;
t1 = t3 * t40;
t18 = [0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t44, -t19 * t41, -t20, -g(1) * (-t42 * pkin(1) + t45 * pkin(7)) - g(2) * (t45 * pkin(1) + t42 * pkin(7)) 0, 0, 0, 0, 0, 0, t19 * t32, -t19 * t31, -t20, -g(1) * (-t42 * t30 - t45 * t39) - g(2) * (t45 * t30 - t42 * t39) 0, 0, 0, 0, 0, 0, t19 * t28, -t19 * t27, -t20, -g(1) * (-t42 * t15 - t45 * t37) - g(2) * (t45 * t15 - t42 * t37) 0, 0, 0, 0, 0, 0, t19 * t25, -t7, -t20, -g(1) * (-t45 * t34 - t42 * t9) - g(2) * (-t42 * t34 + t8) 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13, -g(1) * t10 - g(2) * t12, t7, -g(2) * t8 + (g(1) * t34 - g(2) * t63) * t45 + (-g(1) * (-t63 - t9) + g(2) * t34) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, g(3) * t41 + t20 * t44, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t32 + t20 * t31, g(3) * t31 + t20 * t32, 0, t46 * pkin(2), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(3) * t52 - t20 * t49, 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t51 - t20 * t14, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t50 * t45 + t17) - g(2) * (t50 * t42 + t16) - g(3) * (t51 + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(4), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t48 * t45 + t17) - g(2) * (t48 * t42 + t16) - g(3) * (t23 + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t45 * t61 + t17) - g(2) * (-t42 * t61 + t16) - g(3) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 + g(2) * t10 + t40 * t59, g(1) * t13 - g(2) * t11 + t43 * t59, 0, 0;];
taug_reg  = t18;

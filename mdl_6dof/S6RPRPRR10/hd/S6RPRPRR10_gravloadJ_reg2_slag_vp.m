% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:03:15
% EndTime: 2019-05-05 20:03:16
% DurationCPUTime: 0.43s
% Computational Cost: add. (279->85), mult. (377->122), div. (0->0), fcn. (379->10), ass. (0->60)
t39 = cos(qJ(1));
t37 = sin(qJ(1));
t70 = g(1) * t37;
t75 = g(2) * t39 - t70;
t32 = pkin(10) + qJ(5);
t22 = sin(t32);
t38 = cos(qJ(3));
t67 = g(3) * t38;
t36 = sin(qJ(3));
t23 = cos(t32);
t54 = t39 * t23;
t62 = t37 * t22;
t7 = -t36 * t62 + t54;
t55 = t39 * t22;
t61 = t37 * t23;
t9 = t36 * t55 + t61;
t74 = -g(1) * t7 - g(2) * t9 + t22 * t67;
t12 = -g(3) * t36 - t38 * t75;
t71 = -pkin(1) - pkin(7);
t33 = sin(pkin(10));
t66 = t33 * pkin(4);
t34 = cos(pkin(10));
t21 = t34 * pkin(4) + pkin(3);
t65 = t36 * t39;
t24 = qJ(6) + t32;
t19 = sin(t24);
t64 = t37 * t19;
t20 = cos(t24);
t63 = t37 * t20;
t60 = t37 * t33;
t59 = t37 * t34;
t58 = t38 * t39;
t57 = t39 * t19;
t56 = t39 * t20;
t53 = t39 * t33;
t52 = t39 * t34;
t35 = -pkin(8) - qJ(4);
t51 = t39 * pkin(1) + t37 * qJ(2);
t50 = t38 * qJ(4);
t48 = t39 * pkin(7) + t51;
t47 = g(2) * t48;
t17 = g(1) * t39 + g(2) * t37;
t45 = t36 * pkin(3) - t50;
t14 = pkin(5) * t23 + t21;
t31 = -pkin(9) + t35;
t43 = t36 * t14 + t38 * t31;
t41 = t36 * t21 + t38 * t35;
t27 = t39 * qJ(2);
t15 = pkin(5) * t22 + t66;
t13 = t17 * t38;
t11 = -g(2) * t65 + t36 * t70 + t67;
t10 = t36 * t54 - t62;
t8 = t36 * t61 + t55;
t6 = t36 * t56 - t64;
t5 = t36 * t57 + t63;
t4 = t36 * t63 + t57;
t3 = -t36 * t64 + t56;
t2 = g(1) * t4 - g(2) * t6 + t20 * t67;
t1 = -g(1) * t3 - g(2) * t5 + t19 * t67;
t16 = [0, 0, 0, 0, 0, 0, -t75, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t17, -g(1) * (-t37 * pkin(1) + t27) - g(2) * t51, 0, 0, 0, 0, 0, 0, -t17 * t36, -t13, -t75, -g(1) * (t71 * t37 + t27) - t47, 0, 0, 0, 0, 0, 0, -g(1) * (t36 * t52 - t60) - g(2) * (t36 * t59 + t53) -g(1) * (-t36 * t53 - t59) - g(2) * (-t36 * t60 + t52) t13, -g(1) * (pkin(3) * t65 - t39 * t50 + t27) - t47 + (-g(1) * t71 - g(2) * t45) * t37, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, t13, -g(1) * (t21 * t65 + t35 * t58 + t27) - g(2) * (pkin(4) * t53 + t48) + (-g(1) * (-t66 + t71) - g(2) * t41) * t37, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t13, -g(1) * (t14 * t65 + t31 * t58 + t27) - g(2) * (t39 * t15 + t48) + (-g(1) * (-t15 + t71) - g(2) * t43) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t34, t12 * t33, -t11, g(3) * t45 + t75 * (pkin(3) * t38 + qJ(4) * t36) 0, 0, 0, 0, 0, 0, -t12 * t23, t12 * t22, -t11, g(3) * t41 + t75 * (t21 * t38 - t35 * t36) 0, 0, 0, 0, 0, 0, -t12 * t20, t12 * t19, -t11, g(3) * t43 + t75 * (t14 * t38 - t31 * t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, g(1) * t8 - g(2) * t10 + t23 * t67, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t74 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t16;

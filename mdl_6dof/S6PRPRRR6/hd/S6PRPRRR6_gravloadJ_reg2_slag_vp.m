% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:34:27
% EndTime: 2019-05-05 01:34:28
% DurationCPUTime: 0.50s
% Computational Cost: add. (394->114), mult. (896->175), div. (0->0), fcn. (1088->12), ass. (0->64)
t30 = sin(pkin(11));
t34 = sin(qJ(2));
t37 = cos(qJ(2));
t56 = cos(pkin(11));
t57 = cos(pkin(6));
t45 = t57 * t56;
t15 = t30 * t37 + t34 * t45;
t52 = t30 * t57;
t17 = -t34 * t52 + t56 * t37;
t76 = -g(1) * t17 - g(2) * t15;
t31 = sin(pkin(6));
t73 = g(3) * t31;
t14 = t30 * t34 - t37 * t45;
t72 = t14 * pkin(8);
t16 = t56 * t34 + t37 * t52;
t71 = t16 * pkin(8);
t32 = sin(qJ(5));
t70 = t14 * t32;
t69 = t16 * t32;
t29 = qJ(5) + qJ(6);
t27 = sin(t29);
t33 = sin(qJ(4));
t68 = t27 * t33;
t28 = cos(t29);
t67 = t28 * t33;
t66 = t30 * t31;
t65 = t31 * t34;
t64 = t31 * t37;
t63 = t32 * t33;
t62 = t32 * t37;
t61 = t33 * t34;
t35 = cos(qJ(5));
t60 = t33 * t35;
t59 = t34 * t35;
t58 = pkin(2) * t64 + qJ(3) * t65;
t55 = pkin(8) * t64 + t58;
t12 = t14 * pkin(2);
t54 = -t12 - t72;
t13 = t16 * pkin(2);
t53 = -t13 - t71;
t51 = t31 * t56;
t50 = t15 * qJ(3) - t12;
t49 = t17 * qJ(3) - t13;
t48 = g(3) * t55;
t36 = cos(qJ(4));
t47 = pkin(4) * t33 - pkin(9) * t36;
t26 = t35 * pkin(5) + pkin(4);
t38 = -pkin(10) - pkin(9);
t46 = t26 * t33 + t36 * t38;
t18 = -t57 * t33 - t36 * t64;
t7 = t16 * t36 - t33 * t66;
t9 = t14 * t36 + t33 * t51;
t42 = g(1) * t7 + g(2) * t9 + g(3) * t18;
t10 = -t14 * t33 + t36 * t51;
t19 = -t33 * t64 + t57 * t36;
t8 = t16 * t33 + t36 * t66;
t41 = g(1) * t8 - g(2) * t10 + g(3) * t19;
t5 = -g(1) * t16 - g(2) * t14 + g(3) * t64;
t40 = g(3) * t65 - t76;
t39 = -g(1) * (t17 * t35 - t8 * t32) - g(2) * (t10 * t32 + t15 * t35) - g(3) * (-t19 * t32 + t31 * t59);
t4 = t40 * t36;
t2 = -g(1) * (-t17 * t27 - t8 * t28) - g(2) * (t10 * t28 - t15 * t27) - g(3) * (-t19 * t28 - t27 * t65);
t1 = -g(1) * (t17 * t28 - t8 * t27) - g(2) * (t10 * t27 + t15 * t28) - g(3) * (-t19 * t27 + t28 * t65);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t40, -g(1) * t49 - g(2) * t50 - g(3) * t58, 0, 0, 0, 0, 0, 0, -t40 * t33, -t4, -t5, -g(1) * (t49 - t71) - g(2) * (t50 - t72) - t48, 0, 0, 0, 0, 0, 0, -g(1) * (t17 * t60 - t69) - g(2) * (t15 * t60 - t70) - (t33 * t59 + t62) * t73, -g(1) * (-t16 * t35 - t17 * t63) - g(2) * (-t14 * t35 - t15 * t63) - (-t32 * t61 + t35 * t37) * t73, t4, -g(1) * t53 - g(2) * t54 - g(3) * (t47 * t65 + t55) + t76 * (qJ(3) + t47) 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t27 + t17 * t67) - g(2) * (-t14 * t27 + t15 * t67) - (t27 * t37 + t28 * t61) * t73, -g(1) * (-t16 * t28 - t17 * t68) - g(2) * (-t14 * t28 - t15 * t68) - (-t27 * t61 + t28 * t37) * t73, t4, -g(1) * (-pkin(5) * t69 + t53) - g(2) * (-pkin(5) * t70 + t54) - t48 - (pkin(5) * t62 + t46 * t34) * t73 + t76 * (qJ(3) + t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t41, 0, 0, 0, 0, 0, 0, 0, 0, -t42 * t35, t42 * t32, -t41, -g(1) * (t7 * pkin(4) + t8 * pkin(9)) - g(2) * (t9 * pkin(4) - t10 * pkin(9)) - g(3) * (t18 * pkin(4) + t19 * pkin(9)) 0, 0, 0, 0, 0, 0, -t42 * t28, t42 * t27, -t41, -g(1) * (t7 * t26 - t8 * t38) - g(2) * (t10 * t38 + t9 * t26) - g(3) * (t18 * t26 - t19 * t38); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -g(1) * (-t17 * t32 - t8 * t35) - g(2) * (t10 * t35 - t15 * t32) - g(3) * (-t19 * t35 - t32 * t65) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t39 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;

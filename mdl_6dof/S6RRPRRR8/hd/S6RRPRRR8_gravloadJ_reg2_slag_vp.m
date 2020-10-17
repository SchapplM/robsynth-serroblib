% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR8
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
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:40:42
% EndTime: 2019-05-06 22:40:44
% DurationCPUTime: 0.58s
% Computational Cost: add. (539->111), mult. (540->160), div. (0->0), fcn. (549->12), ass. (0->71)
t52 = sin(qJ(1));
t54 = cos(qJ(1));
t26 = g(1) * t54 + g(2) * t52;
t47 = pkin(11) + qJ(4);
t36 = sin(t47);
t37 = cos(t47);
t70 = t54 * t37;
t53 = cos(qJ(2));
t77 = t52 * t53;
t14 = t36 * t77 + t70;
t71 = t54 * t36;
t16 = t52 * t37 - t53 * t71;
t51 = sin(qJ(2));
t80 = g(3) * t51;
t88 = -g(1) * t16 + g(2) * t14 + t36 * t80;
t38 = qJ(5) + t47;
t33 = cos(t38);
t32 = sin(t38);
t73 = t54 * t32;
t11 = t52 * t33 - t53 * t73;
t72 = t54 * t33;
t9 = t32 * t77 + t72;
t3 = -g(1) * t11 + g(2) * t9 + t32 * t80;
t19 = -g(3) * t53 + t26 * t51;
t84 = g(1) * t52;
t48 = sin(pkin(11));
t40 = t48 * pkin(3);
t49 = cos(pkin(11));
t34 = t49 * pkin(3) + pkin(2);
t78 = t51 * t54;
t76 = t53 * t54;
t35 = qJ(6) + t38;
t28 = sin(t35);
t75 = t54 * t28;
t29 = cos(t35);
t74 = t54 * t29;
t69 = t54 * t48;
t68 = t54 * t49;
t50 = -pkin(8) - qJ(3);
t67 = t54 * pkin(1) + t52 * pkin(7);
t46 = -pkin(9) + t50;
t31 = pkin(4) * t37;
t24 = t31 + t34;
t30 = pkin(4) * t36;
t21 = -pkin(5) * t32 - t30;
t64 = -g(2) * t54 + t84;
t63 = t53 * pkin(2) + t51 * qJ(3);
t27 = pkin(5) * t33;
t13 = t27 + t24;
t39 = -pkin(10) + t46;
t61 = t53 * t13 - t51 * t39;
t59 = t53 * t24 - t51 * t46;
t57 = t53 * t34 - t51 * t50;
t43 = t54 * pkin(7);
t25 = t30 + t40;
t23 = t64 * t51;
t22 = t27 + t31;
t20 = t26 * t53 + t80;
t18 = -t21 + t40;
t17 = t52 * t36 + t53 * t70;
t15 = -t37 * t77 + t71;
t12 = t52 * t32 + t53 * t72;
t10 = -t33 * t77 + t73;
t8 = t52 * t28 + t53 * t74;
t7 = t52 * t29 - t53 * t75;
t6 = -t29 * t77 + t75;
t5 = t28 * t77 + t74;
t4 = g(1) * t12 - g(2) * t10 + t33 * t80;
t2 = g(1) * t8 - g(2) * t6 + t29 * t80;
t1 = -g(1) * t7 + g(2) * t5 + t28 * t80;
t41 = [0, 0, 0, 0, 0, 0, t64, t26, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t53, -t23, -t26, -g(1) * (-t52 * pkin(1) + t43) - g(2) * t67, 0, 0, 0, 0, 0, 0, -g(1) * (-t49 * t77 + t69) - g(2) * (t52 * t48 + t53 * t68) -g(1) * (t48 * t77 + t68) - g(2) * (t52 * t49 - t53 * t69) t23, -g(1) * t43 - g(2) * (t63 * t54 + t67) - (-pkin(1) - t63) * t84, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t23, -g(1) * (pkin(3) * t69 + t43) - g(2) * (t34 * t76 - t50 * t78 + t67) + (-g(1) * (-pkin(1) - t57) - g(2) * t40) * t52, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t23, -g(1) * (t54 * t25 + t43) - g(2) * (t24 * t76 - t46 * t78 + t67) + (-g(1) * (-pkin(1) - t59) - g(2) * t25) * t52, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, t23, -g(1) * (t54 * t18 + t43) - g(2) * (t13 * t76 - t39 * t78 + t67) + (-g(1) * (-pkin(1) - t61) - g(2) * t18) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t49, -t19 * t48, -t20, -g(3) * t63 + t26 * (pkin(2) * t51 - qJ(3) * t53) 0, 0, 0, 0, 0, 0, t19 * t37, -t19 * t36, -t20, -g(3) * t57 + t26 * (t34 * t51 + t50 * t53) 0, 0, 0, 0, 0, 0, t19 * t33, -t19 * t32, -t20, -g(3) * t59 + t26 * (t24 * t51 + t46 * t53) 0, 0, 0, 0, 0, 0, t19 * t29, -t19 * t28, -t20, -g(3) * t61 + t26 * (t13 * t51 + t39 * t53); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, g(1) * t17 - g(2) * t15 + t37 * t80, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t88 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t21 * t76 + t52 * t22) - g(2) * (t21 * t77 - t54 * t22) - t21 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t41;

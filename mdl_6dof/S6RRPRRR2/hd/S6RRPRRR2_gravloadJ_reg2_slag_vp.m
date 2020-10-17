% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR2
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:55:44
% EndTime: 2019-05-06 19:55:46
% DurationCPUTime: 0.46s
% Computational Cost: add. (504->96), mult. (440->129), div. (0->0), fcn. (429->12), ass. (0->68)
t43 = qJ(2) + pkin(11);
t37 = qJ(4) + t43;
t31 = sin(t37);
t32 = cos(t37);
t86 = -pkin(4) * t31 + pkin(9) * t32;
t49 = cos(qJ(5));
t33 = t49 * pkin(5) + pkin(4);
t52 = -pkin(10) - pkin(9);
t85 = -t31 * t52 + t32 * t33;
t84 = t32 * pkin(4) + t31 * pkin(9);
t51 = cos(qJ(1));
t67 = t51 * t49;
t46 = sin(qJ(5));
t48 = sin(qJ(1));
t72 = t48 * t46;
t14 = t32 * t72 + t67;
t68 = t51 * t46;
t71 = t48 * t49;
t16 = -t32 * t68 + t71;
t77 = g(3) * t31;
t83 = -g(1) * t16 + g(2) * t14 + t46 * t77;
t27 = g(1) * t51 + g(2) * t48;
t7 = -g(3) * t32 + t27 * t31;
t36 = cos(t43);
t50 = cos(qJ(2));
t40 = t50 * pkin(2);
t65 = pkin(3) * t36 + t40;
t22 = pkin(1) + t65;
t18 = t51 * t22;
t78 = g(2) * t18;
t44 = qJ(5) + qJ(6);
t38 = sin(t44);
t74 = t48 * t38;
t39 = cos(t44);
t73 = t48 * t39;
t70 = t51 * t38;
t69 = t51 * t39;
t45 = -qJ(3) - pkin(7);
t42 = -pkin(8) + t45;
t63 = pkin(5) * t46 - t42;
t61 = t86 * t48;
t60 = t86 * t51;
t26 = g(1) * t48 - g(2) * t51;
t57 = t31 * t33 + t32 * t52;
t56 = t57 * t48;
t55 = t57 * t51;
t47 = sin(qJ(2));
t53 = -g(3) * t50 + t27 * t47;
t35 = sin(t43);
t34 = t40 + pkin(1);
t23 = -t47 * pkin(2) - pkin(3) * t35;
t20 = t51 * t23;
t19 = t48 * t23;
t17 = t32 * t67 + t72;
t15 = -t32 * t71 + t68;
t13 = t26 * t31;
t12 = t32 * t69 + t74;
t11 = -t32 * t70 + t73;
t10 = -t32 * t73 + t70;
t9 = t32 * t74 + t69;
t8 = t27 * t32 + t77;
t6 = t7 * t49;
t5 = t7 * t46;
t4 = t7 * t39;
t3 = t7 * t38;
t2 = g(1) * t12 - g(2) * t10 + t39 * t77;
t1 = -g(1) * t11 + g(2) * t9 + t38 * t77;
t21 = [0, 0, 0, 0, 0, 0, t26, t27, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t50, -t26 * t47, -t27, -g(1) * (-t48 * pkin(1) + t51 * pkin(7)) - g(2) * (t51 * pkin(1) + t48 * pkin(7)) 0, 0, 0, 0, 0, 0, t26 * t36, -t26 * t35, -t27, -g(1) * (-t48 * t34 - t51 * t45) - g(2) * (t51 * t34 - t48 * t45) 0, 0, 0, 0, 0, 0, t26 * t32, -t13, -t27, -g(1) * (-t48 * t22 - t51 * t42) - g(2) * (-t48 * t42 + t18) 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t13, -t78 + (g(1) * t42 - g(2) * t84) * t51 + (-g(1) * (-t22 - t84) + g(2) * t42) * t48, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, -t78 + (-g(1) * t63 - g(2) * t85) * t51 + (-g(1) * (-t22 - t85) - g(2) * t63) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, g(3) * t47 + t27 * t50, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t36 + t27 * t35, g(3) * t35 + t27 * t36, 0, t53 * pkin(2), 0, 0, 0, 0, 0, 0, t7, t8, 0, -g(3) * t65 - t23 * t27, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t20 + t60) - g(2) * (t19 + t61) - g(3) * (t65 + t84) 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t20 - t55) - g(2) * (t19 - t56) - g(3) * (t85 + t65); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t60 - g(2) * t61 - g(3) * t84, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, g(1) * t55 + g(2) * t56 - g(3) * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, g(1) * t17 - g(2) * t15 + t49 * t77, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t83 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t21;

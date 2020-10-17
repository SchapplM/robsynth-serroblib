% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:06:38
% EndTime: 2019-05-05 22:06:39
% DurationCPUTime: 0.47s
% Computational Cost: add. (442->94), mult. (400->135), div. (0->0), fcn. (403->12), ass. (0->57)
t44 = cos(qJ(4));
t34 = t44 * pkin(4);
t28 = t34 + pkin(3);
t40 = -qJ(5) - pkin(8);
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t49 = t45 * t28 - t42 * t40;
t38 = qJ(4) + pkin(11);
t31 = cos(t38);
t21 = pkin(5) * t31 + t34;
t19 = pkin(3) + t21;
t37 = -pkin(9) + t40;
t51 = t45 * t19 - t42 * t37;
t39 = qJ(1) + pkin(10);
t30 = sin(t39);
t32 = cos(t39);
t18 = g(1) * t32 + g(2) * t30;
t41 = sin(qJ(4));
t64 = t41 * t45;
t13 = t30 * t64 + t32 * t44;
t15 = t30 * t44 - t32 * t64;
t70 = g(3) * t42;
t76 = -g(1) * t15 + g(2) * t13 + t41 * t70;
t11 = -g(3) * t45 + t18 * t42;
t74 = g(1) * t30;
t68 = t41 * pkin(4);
t67 = t30 * t45;
t66 = t32 * t41;
t65 = t32 * t45;
t61 = t44 * t45;
t46 = cos(qJ(1));
t57 = t46 * pkin(1) + t32 * pkin(2) + t30 * pkin(7);
t43 = sin(qJ(1));
t56 = -t43 * pkin(1) + t32 * pkin(7);
t55 = t45 * pkin(3) + t42 * pkin(8);
t53 = -g(2) * t32 + t74;
t52 = g(1) * t43 - g(2) * t46;
t33 = qJ(6) + t38;
t29 = sin(t38);
t27 = cos(t33);
t26 = sin(t33);
t20 = pkin(5) * t29 + t68;
t17 = t53 * t42;
t16 = t30 * t41 + t32 * t61;
t14 = -t30 * t61 + t66;
t12 = t18 * t45 + t70;
t10 = t30 * t29 + t31 * t65;
t9 = -t29 * t65 + t30 * t31;
t8 = t32 * t29 - t31 * t67;
t7 = t29 * t67 + t32 * t31;
t6 = t30 * t26 + t27 * t65;
t5 = -t26 * t65 + t30 * t27;
t4 = t32 * t26 - t27 * t67;
t3 = t26 * t67 + t32 * t27;
t2 = g(1) * t6 - g(2) * t4 + t27 * t70;
t1 = -g(1) * t5 + g(2) * t3 + t26 * t70;
t22 = [0, 0, 0, 0, 0, 0, t52, g(1) * t46 + g(2) * t43, 0, 0, 0, 0, 0, 0, 0, 0, t53, t18, 0, t52 * pkin(1), 0, 0, 0, 0, 0, 0, t53 * t45, -t17, -t18, -g(1) * (-t30 * pkin(2) + t56) - g(2) * t57, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, t17, -g(1) * t56 - g(2) * (t32 * t55 + t57) - (-pkin(2) - t55) * t74, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t17, -g(1) * (pkin(4) * t66 + t56) - g(2) * (t49 * t32 + t57) + (-g(1) * (-pkin(2) - t49) - g(2) * t68) * t30, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t17, -g(1) * (t32 * t20 + t56) - g(2) * (t51 * t32 + t57) + (-g(1) * (-pkin(2) - t51) - g(2) * t20) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t44, -t11 * t41, -t12, -g(3) * t55 + t18 * (pkin(3) * t42 - pkin(8) * t45) 0, 0, 0, 0, 0, 0, t11 * t31, -t11 * t29, -t12, -g(3) * t49 + t18 * (t28 * t42 + t40 * t45) 0, 0, 0, 0, 0, 0, t11 * t27, -t11 * t26, -t12, -g(3) * t51 + t18 * (t19 * t42 + t37 * t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, g(1) * t16 - g(2) * t14 + t44 * t70, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t29 * t70, g(1) * t10 - g(2) * t8 + t31 * t70, 0, t76 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t20 * t65 + t30 * t21) - g(2) * (-t20 * t67 - t32 * t21) + t20 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t22;

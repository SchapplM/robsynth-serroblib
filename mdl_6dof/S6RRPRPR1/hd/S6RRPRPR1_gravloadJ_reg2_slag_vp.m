% Calculate inertial parameters regressor of gravitation load for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:49:49
% EndTime: 2019-05-06 12:49:51
% DurationCPUTime: 0.44s
% Computational Cost: add. (460->94), mult. (395->122), div. (0->0), fcn. (379->12), ass. (0->61)
t39 = qJ(2) + pkin(10);
t34 = qJ(4) + t39;
t26 = sin(t34);
t27 = cos(t34);
t77 = -pkin(4) * t26 + qJ(5) * t27;
t41 = cos(pkin(11));
t28 = t41 * pkin(5) + pkin(4);
t43 = -pkin(9) - qJ(5);
t76 = -t26 * t43 + t27 * t28;
t75 = t27 * pkin(4) + t26 * qJ(5);
t45 = sin(qJ(1));
t47 = cos(qJ(1));
t21 = g(1) * t47 + g(2) * t45;
t5 = -g(3) * t27 + t21 * t26;
t33 = cos(t39);
t46 = cos(qJ(2));
t35 = t46 * pkin(2);
t60 = pkin(3) * t33 + t35;
t16 = pkin(1) + t60;
t12 = t47 * t16;
t73 = g(2) * t12;
t72 = g(3) * t26;
t38 = pkin(11) + qJ(6);
t30 = sin(t38);
t69 = t45 * t30;
t32 = cos(t38);
t68 = t45 * t32;
t40 = sin(pkin(11));
t67 = t45 * t40;
t66 = t45 * t41;
t65 = t47 * t30;
t64 = t47 * t32;
t63 = t47 * t40;
t62 = t47 * t41;
t42 = -qJ(3) - pkin(7);
t37 = -pkin(8) + t42;
t58 = pkin(5) * t40 - t37;
t56 = t77 * t45;
t55 = t77 * t47;
t20 = g(1) * t45 - g(2) * t47;
t52 = t26 * t28 + t27 * t43;
t51 = t52 * t45;
t50 = t52 * t47;
t44 = sin(qJ(2));
t48 = -g(3) * t46 + t21 * t44;
t31 = sin(t39);
t29 = t35 + pkin(1);
t17 = -t44 * pkin(2) - pkin(3) * t31;
t14 = t47 * t17;
t13 = t45 * t17;
t11 = t20 * t26;
t10 = t27 * t64 + t69;
t9 = -t27 * t65 + t68;
t8 = -t27 * t68 + t65;
t7 = t27 * t69 + t64;
t6 = t21 * t27 + t72;
t4 = t5 * t41;
t3 = t5 * t40;
t2 = t5 * t32;
t1 = t5 * t30;
t15 = [0, 0, 0, 0, 0, 0, t20, t21, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t46, -t20 * t44, -t21, -g(1) * (-t45 * pkin(1) + t47 * pkin(7)) - g(2) * (t47 * pkin(1) + t45 * pkin(7)) 0, 0, 0, 0, 0, 0, t20 * t33, -t20 * t31, -t21, -g(1) * (-t45 * t29 - t47 * t42) - g(2) * (t47 * t29 - t45 * t42) 0, 0, 0, 0, 0, 0, t20 * t27, -t11, -t21, -g(1) * (-t45 * t16 - t47 * t37) - g(2) * (-t45 * t37 + t12) 0, 0, 0, 0, 0, 0, -g(1) * (-t27 * t66 + t63) - g(2) * (t27 * t62 + t67) -g(1) * (t27 * t67 + t62) - g(2) * (-t27 * t63 + t66) t11, -t73 + (g(1) * t37 - g(2) * t75) * t47 + (-g(1) * (-t16 - t75) + g(2) * t37) * t45, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t11, -t73 + (-g(1) * t58 - g(2) * t76) * t47 + (-g(1) * (-t16 - t76) - g(2) * t58) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, g(3) * t44 + t21 * t46, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t33 + t21 * t31, g(3) * t31 + t21 * t33, 0, t48 * pkin(2), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(3) * t60 - t21 * t17, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t14 + t55) - g(2) * (t13 + t56) - g(3) * (t60 + t75) 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(1) * (t14 - t50) - g(2) * (t13 - t51) - g(3) * (t76 + t60); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t55 - g(2) * t56 - g(3) * t75, 0, 0, 0, 0, 0, 0, t2, -t1, -t6, g(1) * t50 + g(2) * t51 - g(3) * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t30 * t72, g(1) * t10 - g(2) * t8 + t32 * t72, 0, 0;];
taug_reg  = t15;

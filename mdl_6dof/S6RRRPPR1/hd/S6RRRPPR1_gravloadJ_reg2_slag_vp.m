% Calculate inertial parameters regressor of gravitation load for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:09:02
% EndTime: 2019-05-07 04:09:04
% DurationCPUTime: 0.47s
% Computational Cost: add. (486->99), mult. (421->126), div. (0->0), fcn. (402->12), ass. (0->67)
t41 = qJ(2) + qJ(3);
t34 = pkin(10) + t41;
t27 = sin(t34);
t28 = cos(t34);
t43 = cos(pkin(11));
t29 = t43 * pkin(5) + pkin(4);
t44 = -pkin(9) - qJ(5);
t54 = -t27 * t44 + t28 * t29;
t55 = t28 * pkin(4) + t27 * qJ(5);
t46 = sin(qJ(1));
t48 = cos(qJ(1));
t23 = g(1) * t48 + g(2) * t46;
t5 = -g(3) * t28 + t23 * t27;
t49 = -pkin(8) - pkin(7);
t35 = sin(t41);
t77 = pkin(3) * t35;
t76 = pkin(4) * t27;
t36 = cos(t41);
t30 = pkin(3) * t36;
t47 = cos(qJ(2));
t37 = t47 * pkin(2);
t61 = t30 + t37;
t18 = pkin(1) + t61;
t14 = t48 * t18;
t74 = g(2) * t14;
t72 = g(3) * t27;
t40 = pkin(11) + qJ(6);
t32 = sin(t40);
t69 = t46 * t32;
t33 = cos(t40);
t68 = t46 * t33;
t42 = sin(pkin(11));
t67 = t46 * t42;
t66 = t46 * t43;
t65 = t48 * t32;
t64 = t48 * t33;
t63 = t48 * t42;
t62 = t48 * t43;
t60 = qJ(5) * t28;
t59 = t30 + t55;
t39 = -qJ(4) + t49;
t58 = pkin(5) * t42 - t39;
t57 = t30 + t54;
t56 = -t76 - t77;
t22 = g(1) * t46 - g(2) * t48;
t53 = -t27 * t29 - t28 * t44;
t7 = -g(3) * t36 + t23 * t35;
t45 = sin(qJ(2));
t50 = -g(3) * t47 + t23 * t45;
t31 = t37 + pkin(1);
t21 = t48 * t60;
t20 = t46 * t60;
t19 = -t45 * pkin(2) - t77;
t16 = t48 * t19;
t15 = t46 * t19;
t13 = t22 * t27;
t12 = t28 * t64 + t69;
t11 = -t28 * t65 + t68;
t10 = -t28 * t68 + t65;
t9 = t28 * t69 + t64;
t8 = g(3) * t35 + t23 * t36;
t6 = t23 * t28 + t72;
t4 = t5 * t43;
t3 = t5 * t42;
t2 = t5 * t33;
t1 = t5 * t32;
t17 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t47, -t22 * t45, -t23, -g(1) * (-t46 * pkin(1) + t48 * pkin(7)) - g(2) * (t48 * pkin(1) + t46 * pkin(7)) 0, 0, 0, 0, 0, 0, t22 * t36, -t22 * t35, -t23, -g(1) * (-t46 * t31 - t48 * t49) - g(2) * (t48 * t31 - t46 * t49) 0, 0, 0, 0, 0, 0, t22 * t28, -t13, -t23, -g(1) * (-t46 * t18 - t48 * t39) - g(2) * (-t46 * t39 + t14) 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t66 + t63) - g(2) * (t28 * t62 + t67) -g(1) * (t28 * t67 + t62) - g(2) * (-t28 * t63 + t66) t13, -t74 + (g(1) * t39 - g(2) * t55) * t48 + (-g(1) * (-t18 - t55) + g(2) * t39) * t46, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, -t74 + (-g(1) * t58 - g(2) * t54) * t48 + (-g(1) * (-t18 - t54) - g(2) * t58) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, g(3) * t45 + t23 * t47, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t50 * pkin(2), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(3) * t61 - t23 * t19, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t48 * t76 + t16 + t21) - g(2) * (-t46 * t76 + t15 + t20) - g(3) * (t37 + t59) 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(1) * (t53 * t48 + t16) - g(2) * (t53 * t46 + t15) - g(3) * (t37 + t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(3), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t56 * t48 + t21) - g(2) * (t56 * t46 + t20) - g(3) * t59, 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(3) * t57 + t23 * (-t53 + t77); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 + t32 * t72, g(1) * t12 - g(2) * t10 + t33 * t72, 0, 0;];
taug_reg  = t17;

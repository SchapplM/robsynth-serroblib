% Calculate inertial parameters regressor of gravitation load for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t33 = cos(qJ(1));
t31 = sin(qJ(1));
t67 = g(2) * t31;
t10 = g(1) * t33 + t67;
t30 = sin(qJ(2));
t72 = t10 * t30;
t32 = cos(qJ(2));
t6 = g(3) * t30 + t10 * t32;
t71 = pkin(2) + pkin(3);
t23 = t33 * pkin(7);
t70 = g(1) * t23;
t69 = g(1) * t31;
t65 = g(3) * t32;
t22 = t32 * pkin(2);
t21 = t32 * pkin(3);
t64 = g(2) * qJ(4);
t28 = cos(pkin(9));
t16 = t28 * pkin(5) + pkin(4);
t63 = t30 * t16;
t62 = t30 * t33;
t26 = pkin(9) + qJ(6);
t17 = sin(t26);
t61 = t31 * t17;
t18 = cos(t26);
t60 = t31 * t18;
t27 = sin(pkin(9));
t59 = t31 * t27;
t58 = t31 * t28;
t57 = t31 * t32;
t29 = -pkin(8) - qJ(5);
t56 = t32 * t29;
t55 = t32 * t33;
t54 = t33 * t17;
t53 = t33 * t18;
t52 = t33 * t27;
t51 = t33 * t28;
t19 = t30 * qJ(3);
t50 = t22 + t19;
t49 = t33 * pkin(1) + t31 * pkin(7);
t48 = qJ(3) * t32;
t47 = t32 * qJ(5);
t46 = -t29 + t71;
t45 = qJ(5) + t71;
t44 = t21 + t50;
t43 = -pkin(1) - t19;
t42 = pkin(5) * t27 + qJ(4);
t41 = pkin(2) * t55 + t33 * t19 + t49;
t40 = g(1) * t46;
t39 = g(1) * t45;
t38 = pkin(3) * t55 + t41;
t9 = -g(2) * t33 + t69;
t37 = g(1) * (-t33 * qJ(4) + t23);
t36 = t43 - t22;
t35 = g(2) * t38;
t13 = t33 * t48;
t11 = t31 * t48;
t8 = t9 * t32;
t7 = t9 * t30;
t5 = -t65 + t72;
t4 = t30 * t53 - t61;
t3 = -t30 * t54 - t60;
t2 = -t30 * t60 - t54;
t1 = t30 * t61 - t53;
t12 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (-t31 * pkin(1) + t23) - g(2) * t49, 0, 0, 0, 0, 0, 0, t8, -t10, t7, -g(2) * t41 - t36 * t69 - t70, 0, 0, 0, 0, 0, 0, t7, -t8, t10, -t37 - t35 + (-g(1) * (t36 - t21) + t64) * t31, 0, 0, 0, 0, 0, 0, -g(1) * (-t30 * t58 - t52) - g(2) * (t30 * t51 - t59) -g(1) * (t30 * t59 - t51) - g(2) * (-t30 * t52 - t58) t8, -t37 - g(2) * (pkin(4) * t62 + t33 * t47 + t38) + (-g(1) * (-t30 * pkin(4) + t43) + t64 + t32 * t39) * t31, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3, t8, -t70 - t35 + (g(1) * t42 - g(2) * (-t56 + t63)) * t33 + (-g(1) * (t43 - t63) + g(2) * t42 + t32 * t40) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * (-pkin(2) * t62 + t13) - g(2) * (-t31 * t30 * pkin(2) + t11) - g(3) * t50, 0, 0, 0, 0, 0, 0, -t6, -t5, 0, -g(1) * t13 - g(2) * t11 - g(3) * t44 + t71 * t72, 0, 0, 0, 0, 0, 0, -t6 * t28, t6 * t27, t5, -g(1) * (pkin(4) * t55 + t13) - g(2) * (pkin(4) * t57 + t11) - g(3) * (t44 + t47) + (-g(3) * pkin(4) + t33 * t39 + t45 * t67) * t30, 0, 0, 0, 0, 0, 0, -t6 * t18, t6 * t17, t5, -g(1) * (t16 * t55 + t13) - g(2) * (t16 * t57 + t11) - g(3) * (t44 - t56) + (-g(3) * t16 + t33 * t40 + t46 * t67) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 - t17 * t65, g(1) * t4 - g(2) * t2 - t18 * t65, 0, 0;];
taug_reg  = t12;

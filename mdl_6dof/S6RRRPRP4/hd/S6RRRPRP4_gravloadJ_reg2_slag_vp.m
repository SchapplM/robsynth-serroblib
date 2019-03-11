% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t34 = qJ(2) + qJ(3);
t29 = sin(t34);
t24 = t29 * qJ(4);
t30 = cos(t34);
t62 = t30 * pkin(3) + t24;
t37 = sin(qJ(1));
t40 = cos(qJ(1));
t16 = g(1) * t40 + g(2) * t37;
t6 = g(3) * t29 + t16 * t30;
t79 = pkin(3) + pkin(9);
t36 = sin(qJ(2));
t78 = pkin(2) * t36;
t77 = pkin(3) * t29;
t41 = -pkin(8) - pkin(7);
t74 = g(2) * t41;
t72 = g(3) * t30;
t25 = t30 * pkin(9);
t35 = sin(qJ(5));
t71 = t30 * t35;
t70 = t30 * t40;
t69 = t37 * t35;
t38 = cos(qJ(5));
t68 = t37 * t38;
t67 = t40 * t35;
t66 = t40 * t38;
t65 = t40 * t41;
t61 = qJ(4) * t30;
t17 = t37 * t61;
t59 = pkin(5) * t71;
t64 = t37 * t59 + t17;
t19 = t40 * t61;
t63 = t40 * t59 + t19;
t60 = qJ(6) * t38;
t39 = cos(qJ(2));
t32 = t39 * pkin(2);
t28 = t32 + pkin(1);
t20 = t40 * t28;
t58 = pkin(3) * t70 + t40 * t24 + t20;
t57 = t25 + t62;
t56 = t32 + t62;
t55 = t30 * t60;
t54 = t40 * pkin(4) - t65;
t53 = t25 + t56;
t11 = t29 * t68 + t67;
t9 = -t29 * t66 + t69;
t52 = g(1) * t11 + g(2) * t9;
t51 = t37 * pkin(4) + pkin(9) * t70 + t58;
t50 = -t77 - t78;
t49 = g(1) * t37 - g(2) * t40;
t48 = -t28 - t62;
t1 = g(1) * t9 - g(2) * t11 + t38 * t72;
t10 = t29 * t67 + t68;
t12 = -t29 * t69 + t66;
t47 = -g(1) * t10 + g(2) * t12 + g(3) * t71;
t46 = t16 * t79;
t45 = -g(3) * t39 + t16 * t36;
t43 = t46 * t29;
t42 = (-g(1) * (t48 - t25) + t74) * t37;
t21 = t29 * t35 * pkin(5);
t8 = t49 * t30;
t7 = t49 * t29;
t5 = t16 * t29 - t72;
t4 = t6 * t38;
t3 = t6 * t35;
t2 = -g(1) * t12 - g(2) * t10;
t13 = [0, 0, 0, 0, 0, 0, t49, t16, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t39, -t49 * t36, -t16, -g(1) * (-t37 * pkin(1) + t40 * pkin(7)) - g(2) * (t40 * pkin(1) + t37 * pkin(7)) 0, 0, 0, 0, 0, 0, t8, -t7, -t16, -g(1) * (-t37 * t28 - t65) - g(2) * (-t37 * t41 + t20) 0, 0, 0, 0, 0, 0, -t16, -t8, t7, g(1) * t65 - g(2) * t58 + (-g(1) * t48 + t74) * t37, 0, 0, 0, 0, 0, 0, t2, t52, t8, -g(1) * t54 - g(2) * t51 + t42, 0, 0, 0, 0, 0, 0, t2, t8, -t52, -g(1) * (t12 * pkin(5) + t11 * qJ(6) + t54) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t51) + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, g(3) * t36 + t16 * t39, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t45 * pkin(2), 0, 0, 0, 0, 0, 0, 0, -t5, -t6, -g(1) * (t50 * t40 + t19) - g(2) * (t50 * t37 + t17) - g(3) * t56, 0, 0, 0, 0, 0, 0, -t3, -t4, t5, -g(1) * (-t40 * t78 + t19) - g(2) * (-t37 * t78 + t17) - g(3) * t53 + t43, 0, 0, 0, 0, 0, 0, -t3, t5, t4, -g(1) * t63 - g(2) * t64 - g(3) * (-t29 * t60 + t21 + t53) + t16 * (t79 * t29 + t55 + t78); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, -g(1) * (-t40 * t77 + t19) - g(2) * (-t37 * t77 + t17) - g(3) * t62, 0, 0, 0, 0, 0, 0, -t3, -t4, t5, -g(1) * t19 - g(2) * t17 - g(3) * t57 + t43, 0, 0, 0, 0, 0, 0, -t3, t5, t4, -g(1) * (-t40 * t55 + t63) - g(2) * (-t37 * t55 + t64) - g(3) * (t21 + t57) + (g(3) * t60 + t46) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t47, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, t47, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (t11 * pkin(5) - t12 * qJ(6)) - (-pkin(5) * t38 - qJ(6) * t35) * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;

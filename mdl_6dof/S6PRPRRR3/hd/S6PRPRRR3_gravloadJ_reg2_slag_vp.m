% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t41 = pkin(12) + qJ(4);
t36 = sin(t41);
t37 = cos(t41);
t68 = cos(pkin(6));
t44 = sin(pkin(6));
t48 = sin(qJ(2));
t74 = t44 * t48;
t81 = -t36 * t74 + t68 * t37;
t50 = cos(qJ(2));
t43 = sin(pkin(11));
t61 = t43 * t68;
t67 = cos(pkin(11));
t26 = -t48 * t61 + t50 * t67;
t75 = t43 * t44;
t80 = -t26 * t36 + t37 * t75;
t79 = g(3) * t44;
t45 = cos(pkin(12));
t35 = t45 * pkin(3) + pkin(2);
t38 = qJ(5) + t41;
t34 = cos(t38);
t47 = sin(qJ(6));
t77 = t34 * t47;
t49 = cos(qJ(6));
t76 = t34 * t49;
t73 = t44 * t50;
t72 = t47 * t50;
t71 = t49 * t50;
t46 = -pkin(8) - qJ(3);
t55 = t68 * t67;
t23 = t43 * t48 - t50 * t55;
t24 = t43 * t50 + t48 * t55;
t28 = pkin(4) * t37 + t35;
t40 = -pkin(9) + t46;
t70 = -t23 * t28 - t24 * t40;
t25 = t48 * t67 + t50 * t61;
t69 = -t25 * t28 - t26 * t40;
t33 = sin(t38);
t60 = t44 * t67;
t11 = -t24 * t33 - t34 * t60;
t12 = t24 * t34 - t33 * t60;
t64 = t11 * pkin(5) + t12 * pkin(10);
t13 = -t26 * t33 + t34 * t75;
t14 = t26 * t34 + t33 * t75;
t63 = t13 * pkin(5) + t14 * pkin(10);
t18 = -t33 * t74 + t34 * t68;
t19 = t33 * t68 + t34 * t74;
t62 = t18 * pkin(5) + t19 * pkin(10);
t58 = t80 * pkin(4);
t57 = pkin(5) * t34 + pkin(10) * t33;
t56 = t81 * pkin(4);
t54 = g(1) * t13 + g(2) * t11 + g(3) * t18;
t5 = g(1) * t14 + g(2) * t12 + g(3) * t19;
t53 = -t24 * t36 - t37 * t60;
t7 = -g(1) * t25 - g(2) * t23 + g(3) * t73;
t52 = g(1) * t26 + g(2) * t24 + g(3) * t74;
t51 = t53 * pkin(4);
t22 = t28 * t73;
t6 = t7 * t33;
t2 = t54 * t49;
t1 = t54 * t47;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t52, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t45, t7 * sin(pkin(12)) -t52, -g(1) * (-t25 * pkin(2) + t26 * qJ(3)) - g(2) * (-t23 * pkin(2) + t24 * qJ(3)) - (pkin(2) * t50 + qJ(3) * t48) * t79, 0, 0, 0, 0, 0, 0, -t7 * t37, t7 * t36, -t52, -g(1) * (-t25 * t35 - t26 * t46) - g(2) * (-t23 * t35 - t24 * t46) - (t35 * t50 - t46 * t48) * t79, 0, 0, 0, 0, 0, 0, -t7 * t34, t6, -t52, -g(1) * t69 - g(2) * t70 - g(3) * (-t40 * t74 + t22) 0, 0, 0, 0, 0, 0, -g(1) * (-t25 * t76 + t26 * t47) - g(2) * (-t23 * t76 + t24 * t47) - (t34 * t71 + t47 * t48) * t79, -g(1) * (t25 * t77 + t26 * t49) - g(2) * (t23 * t77 + t24 * t49) - (-t34 * t72 + t48 * t49) * t79, -t6, -g(1) * (-t25 * t57 + t69) - g(2) * (-t23 * t57 + t70) - g(3) * t22 - (-t40 * t48 + t50 * t57) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t80 - g(2) * t53 - g(3) * t81, -g(1) * (-t26 * t37 - t36 * t75) - g(2) * (-t24 * t37 + t36 * t60) - g(3) * (-t36 * t68 - t37 * t74) 0, 0, 0, 0, 0, 0, 0, 0, -t54, t5, 0, -g(1) * t58 - g(2) * t51 - g(3) * t56, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (t58 + t63) - g(2) * (t51 + t64) - g(3) * (t56 + t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * t63 - g(2) * t64 - g(3) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t47 + t25 * t49) - g(2) * (-t12 * t47 + t23 * t49) - g(3) * (-t19 * t47 - t44 * t71) -g(1) * (-t14 * t49 - t25 * t47) - g(2) * (-t12 * t49 - t23 * t47) - g(3) * (-t19 * t49 + t44 * t72) 0, 0;];
taug_reg  = t3;

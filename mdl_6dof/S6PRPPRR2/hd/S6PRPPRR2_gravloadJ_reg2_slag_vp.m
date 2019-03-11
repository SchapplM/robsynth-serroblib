% Calculate inertial parameters regressor of gravitation load for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPPRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t33 = sin(pkin(10));
t35 = cos(pkin(10));
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t61 = sin(pkin(11));
t62 = cos(pkin(11));
t25 = t39 * t61 - t42 * t62;
t36 = cos(pkin(6));
t45 = t25 * t36;
t46 = t39 * t62 + t42 * t61;
t15 = t33 * t45 - t35 * t46;
t67 = t36 * t42;
t51 = -t33 * t67 - t35 * t39;
t83 = t51 * pkin(2) + t15 * pkin(3);
t64 = t46 * t36;
t11 = t33 * t25 - t35 * t64;
t16 = -t35 * t25 - t33 * t64;
t34 = sin(pkin(6));
t24 = t46 * t34;
t48 = -g(1) * t16 + g(2) * t11 - g(3) * t24;
t12 = -t33 * t46 - t35 * t45;
t78 = t12 * pkin(8);
t77 = t15 * pkin(8);
t23 = t25 * t34;
t76 = t23 * pkin(8);
t74 = t33 * t39;
t38 = sin(qJ(5));
t73 = t34 * t38;
t41 = cos(qJ(5));
t72 = t34 * t41;
t71 = t34 * t42;
t68 = t36 * t39;
t37 = sin(qJ(6));
t66 = t37 * t38;
t40 = cos(qJ(6));
t65 = t38 * t40;
t63 = pkin(2) * t71 - t23 * pkin(3);
t59 = t35 * t67;
t54 = t24 * qJ(4) + t63;
t27 = pkin(2) * t59;
t53 = -pkin(2) * t74 + t12 * pkin(3) + t27;
t17 = t23 * t41 - t36 * t38;
t4 = -t15 * t41 - t33 * t73;
t6 = -t12 * t41 + t35 * t73;
t50 = g(1) * t4 + g(2) * t6 + g(3) * t17;
t18 = t23 * t38 + t36 * t41;
t5 = -t15 * t38 + t33 * t72;
t7 = t12 * t38 + t35 * t72;
t49 = g(1) * t5 - g(2) * t7 + g(3) * t18;
t3 = g(1) * t15 + g(2) * t12 - g(3) * t23;
t47 = -t11 * qJ(4) + t53;
t44 = -g(1) * t51 - g(3) * t71;
t43 = qJ(4) * t16 + t83;
t22 = -g(3) * t36 + (-g(1) * t33 + g(2) * t35) * t34;
t1 = t48 * t41;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * (t59 - t74) + t44, -g(1) * (t33 * t68 - t35 * t42) - g(2) * (-t33 * t42 - t35 * t68) + g(3) * t34 * t39, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t48, 0, -g(2) * t27 + (g(2) * t74 + t44) * pkin(2), 0, 0, 0, 0, 0, 0, 0, t3, t48, -g(1) * t43 - g(2) * t47 - g(3) * t54, 0, 0, 0, 0, 0, 0, t48 * t38, t1, -t3, -g(1) * (t43 + t77) - g(2) * (t47 + t78) - g(3) * (t54 - t76) 0, 0, 0, 0, 0, 0, -g(1) * (t15 * t37 + t16 * t65) - g(2) * (-t11 * t65 + t12 * t37) - g(3) * (-t23 * t37 + t24 * t65) -g(1) * (t15 * t40 - t16 * t66) - g(2) * (t11 * t66 + t12 * t40) - g(3) * (-t23 * t40 - t24 * t66) -t1, -g(1) * (t77 + t83) - g(2) * (t53 + t78) - g(3) * (t63 - t76) + t48 * (pkin(5) * t38 - pkin(9) * t41 + qJ(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t49, 0, 0, 0, 0, 0, 0, 0, 0, -t50 * t40, t50 * t37, -t49, -g(1) * (t4 * pkin(5) + t5 * pkin(9)) - g(2) * (t6 * pkin(5) - t7 * pkin(9)) - g(3) * (t17 * pkin(5) + t18 * pkin(9)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t16 * t40 - t5 * t37) - g(2) * (-t11 * t40 + t7 * t37) - g(3) * (-t18 * t37 + t24 * t40) -g(1) * (-t16 * t37 - t5 * t40) - g(2) * (t11 * t37 + t7 * t40) - g(3) * (-t18 * t40 - t24 * t37) 0, 0;];
taug_reg  = t2;

% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t40 = pkin(11) + qJ(3);
t34 = qJ(4) + t40;
t28 = sin(t34);
t29 = cos(t34);
t47 = cos(qJ(5));
t31 = t47 * pkin(5) + pkin(4);
t49 = -pkin(10) - pkin(9);
t82 = -t28 * t49 + t29 * t31;
t81 = t29 * pkin(4) + t28 * pkin(9);
t46 = sin(qJ(1));
t48 = cos(qJ(1));
t24 = g(1) * t48 + g(2) * t46;
t61 = t48 * t47;
t45 = sin(qJ(5));
t66 = t46 * t45;
t14 = t29 * t66 + t61;
t62 = t48 * t45;
t65 = t46 * t47;
t16 = -t29 * t62 + t65;
t71 = g(3) * t28;
t80 = -g(1) * t16 + g(2) * t14 + t45 * t71;
t7 = -g(3) * t29 + t24 * t28;
t32 = sin(t40);
t79 = pkin(3) * t32;
t78 = pkin(4) * t28;
t77 = pkin(9) * t29;
t33 = cos(t40);
t27 = pkin(3) * t33;
t43 = cos(pkin(11));
t30 = t43 * pkin(2) + pkin(1);
t19 = t27 + t30;
t18 = t48 * t19;
t73 = g(2) * t18;
t41 = qJ(5) + qJ(6);
t35 = sin(t41);
t68 = t46 * t35;
t36 = cos(t41);
t67 = t46 * t36;
t64 = t48 * t35;
t63 = t48 * t36;
t44 = -pkin(7) - qJ(2);
t39 = -pkin(8) + t44;
t58 = pkin(5) * t45 - t39;
t56 = -t78 - t79;
t23 = g(1) * t46 - g(2) * t48;
t53 = t28 * t31 + t29 * t49;
t50 = -g(3) * t33 + t24 * t32;
t22 = t48 * t77;
t21 = t46 * t77;
t17 = t29 * t61 + t66;
t15 = -t29 * t65 + t62;
t13 = t23 * t28;
t12 = t29 * t63 + t68;
t11 = -t29 * t64 + t67;
t10 = -t29 * t67 + t64;
t9 = t29 * t68 + t63;
t8 = t24 * t29 + t71;
t6 = t7 * t47;
t5 = t7 * t45;
t4 = t7 * t36;
t3 = t7 * t35;
t2 = g(1) * t12 - g(2) * t10 + t36 * t71;
t1 = -g(1) * t11 + g(2) * t9 + t35 * t71;
t20 = [0, 0, 0, 0, 0, 0, t23, t24, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t43, -t23 * sin(pkin(11)) -t24, -g(1) * (-t46 * pkin(1) + t48 * qJ(2)) - g(2) * (t48 * pkin(1) + t46 * qJ(2)) 0, 0, 0, 0, 0, 0, t23 * t33, -t23 * t32, -t24, -g(1) * (-t46 * t30 - t48 * t44) - g(2) * (t48 * t30 - t46 * t44) 0, 0, 0, 0, 0, 0, t23 * t29, -t13, -t24, -g(1) * (-t46 * t19 - t48 * t39) - g(2) * (-t46 * t39 + t18) 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t13, -t73 + (g(1) * t39 - g(2) * t81) * t48 + (-g(1) * (-t19 - t81) + g(2) * t39) * t46, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, -t73 + (-g(1) * t58 - g(2) * t82) * t48 + (-g(1) * (-t19 - t82) - g(2) * t58) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, g(3) * t32 + t24 * t33, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t50 * pkin(3), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t56 * t48 + t22) - g(2) * (t56 * t46 + t21) - g(3) * (t27 + t81) 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * (t27 + t82) + t24 * (t53 + t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-t48 * t78 + t22) - g(2) * (-t46 * t78 + t21) - g(3) * t81, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t82 + t24 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, g(1) * t17 - g(2) * t15 + t47 * t71, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t80 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t20;

% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t41 = sin(qJ(5));
t44 = cos(qJ(5));
t80 = pkin(5) * t44 + qJ(6) * t41;
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t26 = g(1) * t46 + g(2) * t43;
t40 = qJ(2) + qJ(3);
t37 = qJ(4) + t40;
t32 = sin(t37);
t52 = t26 * t32;
t33 = cos(t37);
t65 = t33 * pkin(4) + t32 * pkin(10);
t47 = -pkin(8) - pkin(7);
t35 = sin(t40);
t79 = pkin(3) * t35;
t78 = pkin(4) * t32;
t74 = g(3) * t32;
t73 = g(3) * t41;
t72 = t32 * t46;
t71 = t33 * t46;
t70 = t43 * t41;
t69 = t43 * t44;
t39 = -pkin(9) + t47;
t68 = t46 * t39;
t67 = t46 * t41;
t66 = t46 * t44;
t36 = cos(t40);
t29 = pkin(3) * t36;
t45 = cos(qJ(2));
t38 = t45 * pkin(2);
t64 = t29 + t38;
t17 = pkin(1) + t64;
t14 = t46 * t17;
t62 = pkin(4) * t71 + pkin(10) * t72 + t14;
t61 = t29 + t65;
t60 = t80 * t33 + t65;
t20 = t43 * t33 * pkin(10);
t59 = -t43 * t78 + t20;
t23 = pkin(10) * t71;
t58 = -pkin(4) * t72 + t23;
t57 = t29 + t60;
t56 = -t78 - t79;
t10 = t33 * t70 + t66;
t12 = t33 * t67 - t69;
t55 = g(1) * t10 - g(2) * t12;
t54 = g(1) * t43 - g(2) * t46;
t1 = g(1) * t12 + g(2) * t10 + t32 * t73;
t11 = t33 * t69 - t67;
t13 = t33 * t66 + t70;
t51 = g(1) * t13 + g(2) * t11 + t44 * t74;
t5 = -g(3) * t33 + t52;
t7 = -g(3) * t36 + t26 * t35;
t42 = sin(qJ(2));
t50 = -g(3) * t45 + t26 * t42;
t49 = (-g(1) * (-t17 - t65) + g(2) * t39) * t43;
t48 = (pkin(4) + t80) * t52;
t34 = t38 + pkin(1);
t18 = -t42 * pkin(2) - t79;
t16 = t46 * t18;
t15 = t43 * t18;
t9 = t54 * t32;
t8 = g(3) * t35 + t26 * t36;
t6 = t26 * t33 + t74;
t4 = t5 * t44;
t3 = -t33 * t73 + t41 * t52;
t2 = g(1) * t11 - g(2) * t13;
t19 = [0, 0, 0, 0, 0, 0, t54, t26, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t45, -t54 * t42, -t26, -g(1) * (-t43 * pkin(1) + t46 * pkin(7)) - g(2) * (t46 * pkin(1) + t43 * pkin(7)) 0, 0, 0, 0, 0, 0, t54 * t36, -t54 * t35, -t26, -g(1) * (-t43 * t34 - t46 * t47) - g(2) * (t46 * t34 - t43 * t47) 0, 0, 0, 0, 0, 0, t54 * t33, -t9, -t26, -g(1) * (-t43 * t17 - t68) - g(2) * (-t43 * t39 + t14) 0, 0, 0, 0, 0, 0, t2, -t55, t9, g(1) * t68 - g(2) * t62 + t49, 0, 0, 0, 0, 0, 0, t2, t9, t55, -g(1) * (-t11 * pkin(5) - t10 * qJ(6) - t68) - g(2) * (t13 * pkin(5) + t12 * qJ(6) + t62) + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, g(3) * t42 + t26 * t45, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t50 * pkin(2), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(3) * t64 - t26 * t18, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t16 + t58) - g(2) * (t15 + t59) - g(3) * (t38 + t61) 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (t16 + t23) - g(2) * (t15 + t20) - g(3) * (t38 + t57) + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(3), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t56 * t46 + t23) - g(2) * (t56 * t43 + t20) - g(3) * t61, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (-t46 * t79 + t23) - g(2) * (-t43 * t79 + t20) - g(3) * t57 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t58 - g(2) * t59 - g(3) * t65, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t23 - g(2) * t20 - g(3) * t60 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t51, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t51, -g(1) * (-t12 * pkin(5) + t13 * qJ(6)) - g(2) * (-t10 * pkin(5) + t11 * qJ(6)) - (-pkin(5) * t41 + qJ(6) * t44) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t19;

% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t80 = g(1) * t46;
t23 = g(2) * t43 + t80;
t39 = qJ(2) + qJ(3);
t36 = sin(t39);
t87 = t23 * t36;
t40 = sin(qJ(6));
t37 = cos(t39);
t41 = sin(qJ(5));
t73 = cos(qJ(5));
t62 = t36 * t73;
t20 = -t37 * t41 + t62;
t12 = t20 * t43;
t71 = t37 * t46;
t14 = t41 * t71 - t46 * t62;
t19 = t36 * t41 + t37 * t73;
t51 = g(1) * t14 - g(2) * t12 + g(3) * t19;
t1 = t51 * t40;
t44 = cos(qJ(6));
t2 = t51 * t44;
t31 = t36 * qJ(4);
t69 = t37 * pkin(3) + t31;
t86 = -t19 * pkin(5) + t20 * pkin(10);
t15 = t19 * t46;
t85 = -t14 * pkin(5) + t15 * pkin(10);
t13 = t19 * t43;
t84 = t12 * pkin(5) + t13 * pkin(10);
t42 = sin(qJ(2));
t82 = pkin(2) * t42;
t81 = pkin(3) * t36;
t78 = g(3) * t20;
t32 = t37 * pkin(4);
t47 = -pkin(8) - pkin(7);
t74 = pkin(9) + t47;
t70 = t46 * t47;
t68 = qJ(4) * t37;
t67 = t43 * t82;
t66 = t46 * t82;
t45 = cos(qJ(2));
t38 = t45 * pkin(2);
t35 = t38 + pkin(1);
t27 = t46 * t35;
t65 = pkin(3) * t71 + t46 * t31 + t27;
t64 = t32 + t69;
t63 = t38 + t69;
t61 = pkin(4) * t71 + t65;
t24 = t43 * t68;
t60 = t24 - t84;
t26 = t46 * t68;
t59 = t26 - t85;
t58 = -t81 - t82;
t57 = -g(1) * t12 - g(2) * t14;
t56 = g(1) * t43 - g(2) * t46;
t55 = t13 * t44 + t46 * t40;
t54 = t13 * t40 - t46 * t44;
t53 = -t35 - t69;
t52 = t64 - t86;
t5 = g(1) * t15 + g(2) * t13 + t78;
t50 = -g(3) * t45 + t23 * t42;
t49 = (pkin(3) + pkin(4)) * t87;
t48 = (-g(1) * (t53 - t32) + g(2) * t74) * t43;
t17 = t56 * t37;
t16 = t56 * t36;
t9 = g(3) * t36 + t23 * t37;
t8 = -g(3) * t37 + t87;
t7 = t15 * t44 - t43 * t40;
t6 = -t15 * t40 - t43 * t44;
t3 = [0, 0, 0, 0, 0, 0, t56, t23, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t45, -t56 * t42, -t23, -g(1) * (-t43 * pkin(1) + t46 * pkin(7)) - g(2) * (t46 * pkin(1) + t43 * pkin(7)) 0, 0, 0, 0, 0, 0, t17, -t16, -t23, -g(1) * (-t43 * t35 - t70) - g(2) * (-t43 * t47 + t27) 0, 0, 0, 0, 0, 0, t17, -t23, t16, g(1) * t70 - g(2) * t65 + (-g(1) * t53 + g(2) * t47) * t43, 0, 0, 0, 0, 0, 0, g(1) * t13 - g(2) * t15, -t57, t23, -g(2) * t61 + t74 * t80 + t48, 0, 0, 0, 0, 0, 0, g(1) * t55 - g(2) * t7, -g(1) * t54 - g(2) * t6, t57, -g(1) * (-t13 * pkin(5) - t46 * pkin(9) + t12 * pkin(10) - t70) - g(2) * (t15 * pkin(5) + t14 * pkin(10) + t61) + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, g(3) * t42 + t23 * t45, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t50 * pkin(2), 0, 0, 0, 0, 0, 0, t8, 0, -t9, -g(1) * (t58 * t46 + t26) - g(2) * (t58 * t43 + t24) - g(3) * t63, 0, 0, 0, 0, 0, 0, -t51, -t5, 0, -g(1) * (t26 - t66) - g(2) * (t24 - t67) - g(3) * (t32 + t63) + t49, 0, 0, 0, 0, 0, 0, -t2, t1, t5, -g(1) * (t59 - t66) - g(2) * (t60 - t67) - g(3) * (t38 + t52) + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, -t9, -g(1) * (-t46 * t81 + t26) - g(2) * (-t43 * t81 + t24) - g(3) * t69, 0, 0, 0, 0, 0, 0, -t51, -t5, 0, -g(1) * t26 - g(2) * t24 - g(3) * t64 + t49, 0, 0, 0, 0, 0, 0, -t2, t1, t5, -g(1) * t59 - g(2) * t60 - g(3) * t52 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t5, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t5, -g(1) * t85 - g(2) * t84 - g(3) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 + g(2) * t54 + t40 * t78, g(1) * t7 + g(2) * t55 + t44 * t78, 0, 0;];
taug_reg  = t3;

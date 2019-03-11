% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP4
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t47 = qJ(4) + qJ(5);
t42 = sin(t47);
t44 = cos(t47);
t101 = pkin(5) * t44 + qJ(6) * t42;
t52 = cos(qJ(4));
t40 = t52 * pkin(4) + pkin(3);
t48 = qJ(2) + qJ(3);
t45 = cos(t48);
t25 = t45 * t40;
t43 = sin(t48);
t55 = -pkin(10) - pkin(9);
t89 = t43 * t55;
t74 = t25 - t89;
t51 = sin(qJ(1));
t54 = cos(qJ(1));
t29 = g(1) * t54 + g(2) * t51;
t9 = -g(3) * t45 + t29 * t43;
t100 = t101 * t45;
t99 = t45 * pkin(3) + t43 * pkin(9);
t50 = sin(qJ(2));
t98 = pkin(2) * t50;
t97 = pkin(3) * t43;
t93 = g(3) * t43;
t91 = t42 * t43;
t90 = t43 * t44;
t88 = t45 * t54;
t87 = t45 * t55;
t86 = t51 * t42;
t85 = t51 * t44;
t49 = sin(qJ(4));
t84 = t51 * t49;
t83 = t51 * t52;
t82 = t54 * t42;
t81 = t54 * t44;
t80 = t54 * t49;
t79 = t54 * t52;
t56 = -pkin(8) - pkin(7);
t78 = t54 * t56;
t75 = t45 * t80;
t53 = cos(qJ(2));
t46 = t53 * pkin(2);
t41 = t46 + pkin(1);
t30 = t54 * t41;
t73 = -t51 * t56 + t30;
t13 = t45 * t86 + t81;
t14 = t45 * t85 - t82;
t72 = -t13 * pkin(5) + t14 * qJ(6);
t15 = t45 * t82 - t85;
t16 = t45 * t81 + t86;
t71 = -t15 * pkin(5) + t16 * qJ(6);
t70 = t46 + t74;
t69 = -t97 - t98;
t67 = g(1) * t13 - g(2) * t15;
t66 = g(1) * t51 - g(2) * t54;
t65 = t40 * t43 + t87;
t18 = t45 * t84 + t79;
t64 = t40 + t101;
t62 = t29 * t45;
t1 = g(1) * t15 + g(2) * t13 + g(3) * t91;
t3 = g(1) * t16 + g(2) * t14 + g(3) * t90;
t60 = pkin(4) * t84 + t40 * t88 - t54 * t89 + t73;
t10 = t62 + t93;
t59 = -g(3) * t53 + t29 * t50;
t58 = pkin(4) * t80 - t78 + (-t41 - t74) * t51;
t36 = pkin(4) * t83;
t32 = pkin(9) * t88;
t31 = t51 * t45 * pkin(9);
t23 = qJ(6) * t90;
t21 = t45 * t79 + t84;
t20 = -t75 + t83;
t19 = -t45 * t83 + t80;
t17 = t66 * t43;
t8 = t9 * t52;
t7 = t9 * t49;
t6 = t9 * t44;
t5 = t9 * t42;
t4 = g(1) * t14 - g(2) * t16;
t2 = [0, 0, 0, 0, 0, 0, t66, t29, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t53, -t66 * t50, -t29, -g(1) * (-t51 * pkin(1) + t54 * pkin(7)) - g(2) * (t54 * pkin(1) + t51 * pkin(7)) 0, 0, 0, 0, 0, 0, t66 * t45, -t17, -t29, -g(1) * (-t51 * t41 - t78) - g(2) * t73, 0, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t21, -g(1) * t18 - g(2) * t20, t17, -g(2) * t30 + (g(1) * t56 - g(2) * t99) * t54 + (-g(1) * (-t41 - t99) + g(2) * t56) * t51, 0, 0, 0, 0, 0, 0, t4, -t67, t17, -g(1) * t58 - g(2) * t60, 0, 0, 0, 0, 0, 0, t4, t17, t67, -g(1) * (-t14 * pkin(5) - t13 * qJ(6) + t58) - g(2) * (t16 * pkin(5) + t15 * qJ(6) + t60); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, g(3) * t50 + t29 * t53, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t59 * pkin(2), 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (t69 * t54 + t32) - g(2) * (t69 * t51 + t31) - g(3) * (t46 + t99) 0, 0, 0, 0, 0, 0, t6, -t5, -t10, -g(3) * t70 + t29 * (t65 + t98) 0, 0, 0, 0, 0, 0, t6, -t10, t5, -g(3) * (t70 + t100) + t29 * (t64 * t43 + t87 + t98); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (-t54 * t97 + t32) - g(2) * (-t51 * t97 + t31) - g(3) * t99, 0, 0, 0, 0, 0, 0, t6, -t5, -t10, -g(3) * t74 + t29 * t65, 0, 0, 0, 0, 0, 0, t6, -t10, t5, -g(3) * (t25 + t100) + t55 * t62 + (g(3) * t55 + t29 * t64) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t20 + g(2) * t18 + t49 * t93, g(1) * t21 - g(2) * t19 + t52 * t93, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t36 + (g(2) * t79 + t10 * t49) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (-pkin(4) * t75 + t36 + t71) - g(2) * (-t18 * pkin(4) + t72) - g(3) * (t23 + (-pkin(4) * t49 - pkin(5) * t42) * t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t71 - g(2) * t72 - g(3) * (-pkin(5) * t91 + t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;

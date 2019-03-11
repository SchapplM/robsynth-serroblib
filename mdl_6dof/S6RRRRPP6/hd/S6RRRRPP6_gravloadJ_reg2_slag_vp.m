% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t41 = qJ(3) + qJ(4);
t36 = sin(t41);
t37 = cos(t41);
t47 = cos(qJ(1));
t80 = t47 * t37;
t44 = sin(qJ(1));
t46 = cos(qJ(2));
t82 = t44 * t46;
t13 = t36 * t82 + t80;
t81 = t47 * t36;
t14 = t37 * t82 - t81;
t96 = -t14 * pkin(4) - t13 * qJ(5);
t91 = g(2) * t44;
t24 = g(1) * t47 + t91;
t43 = sin(qJ(2));
t50 = -g(3) * t46 + t24 * t43;
t42 = sin(qJ(3));
t95 = pkin(3) * t42;
t94 = pkin(4) * t37;
t93 = g(1) * t44;
t90 = g(3) * t43;
t48 = -pkin(9) - pkin(8);
t88 = pkin(5) - t48;
t87 = t36 * t43;
t86 = t37 * t43;
t85 = t43 * t48;
t84 = t44 * t42;
t45 = cos(qJ(3));
t83 = t44 * t45;
t35 = t45 * pkin(3) + pkin(2);
t27 = t46 * t35;
t79 = t47 * t42;
t78 = t47 * t45;
t77 = -pkin(4) - qJ(6);
t76 = t47 * pkin(1) + t44 * pkin(7);
t75 = qJ(5) * t36;
t73 = t13 * qJ(6);
t15 = -t44 * t37 + t46 * t81;
t72 = t15 * qJ(6);
t71 = t47 * t85;
t70 = t46 * t79;
t39 = t47 * pkin(7);
t69 = pkin(3) * t79 + t44 * t85 + t39;
t68 = t88 * t47;
t67 = -pkin(1) - t27;
t66 = t77 * t36;
t65 = -t13 * pkin(4) + t14 * qJ(5);
t16 = t44 * t36 + t46 * t80;
t64 = -t15 * pkin(4) + t16 * qJ(5);
t63 = -t35 - t75;
t62 = pkin(3) * t84 + t47 * t27 + t76;
t61 = g(3) * (t27 + (t75 + t94) * t46);
t60 = t46 * pkin(2) + t43 * pkin(8);
t5 = g(1) * t13 - g(2) * t15;
t6 = g(1) * t14 - g(2) * t16;
t58 = -g(2) * t47 + t93;
t18 = t42 * t82 + t78;
t55 = t24 * t46;
t54 = t16 * pkin(4) + t15 * qJ(5) + t62;
t2 = g(1) * t15 + g(2) * t13 + g(3) * t87;
t4 = g(1) * t16 + g(2) * t14 + g(3) * t86;
t52 = t67 * t44 + t69;
t33 = pkin(3) * t83;
t51 = -pkin(3) * t70 + t33 + t64;
t17 = t55 + t90;
t49 = -t18 * pkin(3) + t65;
t25 = qJ(5) * t86;
t22 = t58 * t43;
t21 = t46 * t78 + t84;
t20 = -t70 + t83;
t19 = -t45 * t82 + t79;
t8 = t50 * t37;
t7 = t50 * t36;
t1 = [0, 0, 0, 0, 0, 0, t58, t24, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t46, -t22, -t24, -g(1) * (-t44 * pkin(1) + t39) - g(2) * t76, 0, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t21, -g(1) * t18 - g(2) * t20, t22, -g(1) * t39 - g(2) * (t60 * t47 + t76) - (-pkin(1) - t60) * t93, 0, 0, 0, 0, 0, 0, t6, -t5, t22, -g(1) * t52 - g(2) * (t62 - t71) 0, 0, 0, 0, 0, 0, t22, -t6, t5, -g(1) * (t52 + t96) - g(2) * (t54 - t71) 0, 0, 0, 0, 0, 0, t22, t5, t6, -g(1) * (-t14 * qJ(6) + t69 + t96) - g(2) * (t16 * qJ(6) + t43 * t68 + t54) - (-t43 * pkin(5) + t67) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t17, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t45, -t50 * t42, -t17, -g(3) * t60 + t24 * (pkin(2) * t43 - pkin(8) * t46) 0, 0, 0, 0, 0, 0, t8, -t7, -t17, -g(3) * (t27 - t85) + t24 * (t35 * t43 + t46 * t48) 0, 0, 0, 0, 0, 0, -t17, -t8, t7, -t61 + t48 * t55 + (g(3) * t48 + t24 * (-t63 + t94)) * t43, 0, 0, 0, 0, 0, 0, -t17, t7, t8, -t61 + (-g(3) * qJ(6) * t37 - g(1) * t68 - t88 * t91) * t46 + (-g(3) * t88 + t24 * (-t77 * t37 - t63)) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t20 + g(2) * t18 + t42 * t90, g(1) * t21 - g(2) * t19 + t45 * t90, 0, 0, 0, 0, 0, 0, 0, 0, t2, t4, 0, -g(1) * t33 + (g(2) * t78 + t17 * t42) * pkin(3), 0, 0, 0, 0, 0, 0, 0, -t2, -t4, -g(1) * t51 - g(2) * t49 - g(3) * (t25 + (-pkin(4) * t36 - t95) * t43) 0, 0, 0, 0, 0, 0, 0, -t4, t2, -g(1) * (t51 - t72) - g(2) * (t49 - t73) - g(3) * t25 - (t66 - t95) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t4, -g(1) * t64 - g(2) * t65 - g(3) * (-pkin(4) * t87 + t25) 0, 0, 0, 0, 0, 0, 0, -t4, t2, -g(1) * (t64 - t72) - g(2) * (t65 - t73) - g(3) * (t43 * t66 + t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4;];
taug_reg  = t1;

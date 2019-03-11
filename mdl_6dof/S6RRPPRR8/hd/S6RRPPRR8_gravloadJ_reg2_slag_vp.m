% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t46 = cos(qJ(2));
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t20 = g(1) * t47 + g(2) * t44;
t43 = sin(qJ(2));
t95 = t20 * t43;
t97 = -g(3) * t46 + t95;
t40 = sin(pkin(10));
t41 = cos(pkin(10));
t42 = sin(qJ(5));
t45 = cos(qJ(5));
t58 = t40 * t45 - t41 * t42;
t54 = g(3) * t58;
t80 = t47 * t40;
t17 = -t44 * t41 + t46 * t80;
t79 = t47 * t41;
t18 = t44 * t40 + t46 * t79;
t6 = t17 * t45 - t18 * t42;
t82 = t44 * t46;
t15 = t40 * t82 + t79;
t16 = t41 * t82 - t80;
t62 = t15 * t45 - t16 * t42;
t96 = -g(1) * t6 - g(2) * t62 - t43 * t54;
t92 = g(1) * t44;
t89 = g(3) * t43;
t35 = t46 * pkin(2);
t87 = t15 * t42;
t86 = t40 * t42;
t85 = t41 * t46;
t84 = t43 * t44;
t83 = t43 * t47;
t81 = t46 * t47;
t48 = -pkin(9) - pkin(8);
t78 = t47 * t48;
t33 = t43 * qJ(3);
t77 = t35 + t33;
t76 = t47 * pkin(1) + t44 * pkin(7);
t75 = qJ(3) * t46;
t74 = qJ(4) * t40;
t73 = -pkin(1) - t35;
t72 = -pkin(2) - t74;
t71 = pkin(5) * t42 + qJ(4);
t70 = pkin(3) * t85 + t46 * t74 + t77;
t69 = pkin(2) * t81 + t47 * t33 + t76;
t36 = t47 * pkin(7);
t68 = -t16 * pkin(3) - t15 * qJ(4) + t36;
t67 = t18 * pkin(3) + t69;
t66 = g(1) * t15 - g(2) * t17;
t65 = -g(2) * t47 + t92;
t39 = qJ(5) + qJ(6);
t31 = sin(t39);
t32 = cos(t39);
t64 = t15 * t32 - t16 * t31;
t63 = t15 * t31 + t16 * t32;
t61 = t16 * t45 + t87;
t60 = t31 * t41 - t32 * t40;
t59 = t31 * t40 + t32 * t41;
t57 = t41 * t45 + t86;
t55 = g(3) * t60;
t52 = t17 * qJ(4) + t67;
t51 = (t73 - t33) * t92;
t30 = t45 * pkin(5) + pkin(4);
t25 = t47 * t75;
t22 = t44 * t75;
t19 = g(1) * t84 - g(2) * t83;
t12 = t20 * t46 + t89;
t10 = t97 * t41;
t9 = t97 * t40;
t8 = g(1) * t16 - g(2) * t18;
t7 = t17 * t42 + t18 * t45;
t5 = -g(1) * t17 - g(2) * t15 - t40 * t89;
t4 = t17 * t31 + t18 * t32;
t3 = t17 * t32 - t18 * t31;
t2 = g(1) * t4 + g(2) * t63 + t59 * t89;
t1 = -g(1) * t3 - g(2) * t64 + t43 * t55;
t11 = [0, 0, 0, 0, 0, 0, t65, t20, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t46, -t19, -t20, -g(1) * (-t44 * pkin(1) + t36) - g(2) * t76, 0, 0, 0, 0, 0, 0, t8, -t66, t19, -g(1) * t36 - g(2) * t69 - t51, 0, 0, 0, 0, 0, 0, t8, t19, t66, -g(1) * t68 - g(2) * t52 - t51, 0, 0, 0, 0, 0, 0, g(1) * t61 - g(2) * t7, g(1) * t62 - g(2) * t6, -t19, -g(1) * (-t16 * pkin(4) + t68) - g(2) * (t18 * pkin(4) - pkin(8) * t83 + t52) - ((pkin(8) - qJ(3)) * t43 + t73) * t92, 0, 0, 0, 0, 0, 0, g(1) * t63 - g(2) * t4, g(1) * t64 - g(2) * t3, -t19, -g(1) * (-pkin(5) * t87 - t16 * t30 + t68) - g(2) * (t71 * t17 + t18 * t30 + t43 * t78 + t67) - ((-qJ(3) - t48) * t43 + t73) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t12, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, -t12, -g(1) * (-pkin(2) * t83 + t25) - g(2) * (-pkin(2) * t84 + t22) - g(3) * t77, 0, 0, 0, 0, 0, 0, t10, -t12, t9, -g(1) * t25 - g(2) * t22 - g(3) * t70 + (pkin(3) * t41 - t72) * t95, 0, 0, 0, 0, 0, 0, t97 * t57, -t46 * t54 + t58 * t95, t12, -g(1) * (-pkin(8) * t81 + t25) - g(2) * (-pkin(8) * t82 + t22) - g(3) * (pkin(4) * t85 + t70) + (g(3) * pkin(8) + t20 * (-(-pkin(3) - pkin(4)) * t41 - t72)) * t43, 0, 0, 0, 0, 0, 0, t97 * t59, t46 * t55 - t60 * t95, t12, -g(1) * (t46 * t78 + t25) - g(2) * (t48 * t82 + t22) - g(3) * (t46 * pkin(5) * t86 + t30 * t85 + t70) + (-g(3) * t48 + t20 * (pkin(2) - (-pkin(3) - t30) * t41 + t71 * t40)) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, g(1) * t7 + g(2) * t61 + t57 * t89, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t96 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t11;

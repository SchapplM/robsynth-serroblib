% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t47 = sin(qJ(1));
t50 = cos(qJ(1));
t24 = g(1) * t50 + g(2) * t47;
t46 = sin(qJ(2));
t98 = t24 * t46;
t35 = t46 * qJ(3);
t49 = cos(qJ(2));
t62 = t49 * pkin(2) + t35;
t48 = cos(qJ(4));
t63 = t50 * t48;
t45 = sin(qJ(4));
t74 = t47 * t45;
t15 = t46 * t63 - t74;
t64 = t50 * t45;
t73 = t47 * t48;
t17 = t46 * t73 + t64;
t87 = g(3) * t49;
t97 = -g(1) * t15 - g(2) * t17 + t48 * t87;
t44 = qJ(4) + qJ(5);
t33 = sin(t44);
t66 = t50 * t33;
t34 = cos(t44);
t75 = t47 * t34;
t11 = t46 * t75 + t66;
t65 = t50 * t34;
t76 = t47 * t33;
t9 = t46 * t65 - t76;
t3 = -g(1) * t9 - g(2) * t11 + t34 * t87;
t14 = g(3) * t46 + t24 * t49;
t51 = -pkin(9) - pkin(8);
t93 = g(1) * t47;
t86 = t45 * pkin(4);
t85 = t49 * pkin(8);
t82 = t45 * t49;
t81 = t46 * t47;
t80 = t46 * t50;
t22 = pkin(5) * t33 + t86;
t79 = t47 * t22;
t36 = qJ(6) + t44;
t30 = sin(t36);
t78 = t47 * t30;
t31 = cos(t36);
t77 = t47 * t31;
t43 = -pkin(10) + t51;
t72 = t49 * t43;
t71 = t49 * t50;
t70 = t49 * t51;
t69 = t50 * t22;
t68 = t50 * t30;
t67 = t50 * t31;
t38 = t48 * pkin(4);
t23 = pkin(5) * t34 + t38;
t61 = t50 * pkin(1) + t47 * pkin(7);
t60 = qJ(3) * t49;
t59 = pkin(4) * t82;
t56 = t46 * t64;
t55 = pkin(2) * t71 + t50 * t35 + t61;
t54 = -g(2) * t50 + t93;
t53 = -pkin(1) - t62;
t40 = t50 * pkin(7);
t32 = t38 + pkin(3);
t27 = t50 * t60;
t25 = t47 * t60;
t21 = pkin(3) + t23;
t20 = t54 * t49;
t19 = t54 * t46;
t18 = -t46 * t74 + t63;
t16 = t56 + t73;
t13 = -t87 + t98;
t12 = -t46 * t76 + t65;
t10 = t46 * t66 + t75;
t8 = -t46 * t78 + t67;
t7 = t46 * t77 + t68;
t6 = t46 * t68 + t77;
t5 = t46 * t67 - t78;
t4 = g(1) * t10 - g(2) * t12 - t33 * t87;
t2 = g(1) * t6 - g(2) * t8 - t30 * t87;
t1 = -g(1) * t5 - g(2) * t7 + t31 * t87;
t26 = [0, 0, 0, 0, 0, 0, t54, t24, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, -t24, -g(1) * (-t47 * pkin(1) + t40) - g(2) * t61, 0, 0, 0, 0, 0, 0, -t24, -t20, t19, -g(1) * t40 - g(2) * t55 - t53 * t93, 0, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t16, g(1) * t17 - g(2) * t15, t20, -g(1) * (t50 * pkin(3) + t40) - g(2) * (pkin(8) * t71 + t55) + (-g(1) * (t53 - t85) - g(2) * pkin(3)) * t47, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t20, -g(1) * (t50 * t32 + t40) - g(2) * (pkin(4) * t56 - t50 * t70 + t55) + (-g(1) * (-t46 * t86 + t53 + t70) - g(2) * t32) * t47, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5, t20, -g(1) * (t50 * t21 + t40) - g(2) * (-t43 * t71 + t46 * t69 + t55) + (-g(1) * (-t46 * t22 + t53 + t72) - g(2) * t21) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, -g(1) * (-pkin(2) * t80 + t27) - g(2) * (-pkin(2) * t81 + t25) - g(3) * t62, 0, 0, 0, 0, 0, 0, -t14 * t45, -t14 * t48, t13, -g(1) * t27 - g(2) * t25 - g(3) * (t62 + t85) + (pkin(2) + pkin(8)) * t98, 0, 0, 0, 0, 0, 0, -t14 * t33, -t14 * t34, t13, -g(1) * (t50 * t59 + t27) - g(2) * (t47 * t59 + t25) - g(3) * (t62 - t70) + (-g(3) * t86 + t24 * (pkin(2) - t51)) * t46, 0, 0, 0, 0, 0, 0, -t14 * t30, -t14 * t31, t13, -g(1) * (t49 * t69 + t27) - g(2) * (t49 * t79 + t25) - g(3) * (t62 - t72) + (-g(3) * t22 + t24 * (pkin(2) - t43)) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, g(1) * t16 - g(2) * t18 - g(3) * t82, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t97 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t23 * t80 - t79) - g(2) * (t23 * t81 + t69) + t23 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t26;

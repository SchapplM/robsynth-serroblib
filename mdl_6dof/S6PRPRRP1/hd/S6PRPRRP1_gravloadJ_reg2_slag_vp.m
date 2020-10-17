% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:32:30
% EndTime: 2019-05-04 23:32:31
% DurationCPUTime: 0.53s
% Computational Cost: add. (567->107), mult. (1467->165), div. (0->0), fcn. (1871->12), ass. (0->66)
t41 = sin(pkin(11));
t48 = sin(qJ(2));
t51 = cos(qJ(2));
t76 = cos(pkin(11));
t30 = -t41 * t51 - t48 * t76;
t42 = sin(pkin(10));
t43 = cos(pkin(10));
t61 = -t41 * t48 + t51 * t76;
t44 = cos(pkin(6));
t78 = t30 * t44;
t20 = t42 * t78 + t43 * t61;
t15 = -t42 * t61 + t43 * t78;
t85 = t42 * t48;
t83 = t44 * t48;
t82 = t44 * t51;
t46 = sin(qJ(5));
t50 = cos(qJ(4));
t81 = t46 * t50;
t49 = cos(qJ(5));
t80 = t49 * t50;
t75 = sin(pkin(6));
t63 = t76 * t75;
t70 = t48 * t75;
t27 = t41 * t70 - t51 * t63;
t68 = t51 * t75;
t77 = pkin(2) * t68 - pkin(3) * t27;
t74 = t43 * t82;
t73 = pkin(5) * t46 + pkin(8);
t47 = sin(qJ(4));
t72 = t47 * t75;
t69 = t50 * t75;
t28 = t41 * t68 + t48 * t63;
t67 = t28 * pkin(8) + t77;
t66 = pkin(4) * t50 + pkin(9) * t47;
t40 = pkin(5) * t49 + pkin(4);
t45 = -qJ(6) - pkin(9);
t65 = t40 * t50 - t45 * t47;
t55 = t61 * t44;
t16 = t42 * t30 + t43 * t55;
t31 = pkin(2) * t74;
t64 = -pkin(2) * t85 + pkin(3) * t16 + t31;
t62 = -t42 * t82 - t43 * t48;
t11 = t20 * t47 - t42 * t69;
t21 = t28 * t47 - t44 * t50;
t9 = -t15 * t47 + t43 * t69;
t60 = g(1) * t11 + g(2) * t9 + g(3) * t21;
t10 = -t15 * t50 - t43 * t72;
t12 = t20 * t50 + t42 * t72;
t22 = t28 * t50 + t44 * t47;
t59 = g(1) * t12 + g(2) * t10 + g(3) * t22;
t58 = -g(1) * t20 + g(2) * t15 - g(3) * t28;
t19 = t30 * t43 - t42 * t55;
t57 = g(1) * t19 + g(2) * t16 - g(3) * t27;
t56 = -t15 * pkin(8) + t64;
t54 = pkin(2) * t62 + pkin(3) * t19;
t53 = pkin(8) * t20 + t54;
t52 = -g(1) * t62 - g(3) * t68;
t1 = -g(1) * (-t12 * t46 - t19 * t49) - g(2) * (-t10 * t46 - t16 * t49) - g(3) * (-t22 * t46 + t27 * t49);
t26 = -g(3) * t44 + (-g(1) * t42 + g(2) * t43) * t75;
t8 = t57 * t47;
t6 = t60 * t49;
t5 = t60 * t46;
t4 = -g(1) * (t19 * t80 + t20 * t46) - g(2) * (-t15 * t46 + t16 * t80) - g(3) * (-t27 * t80 + t28 * t46);
t3 = -g(1) * (-t19 * t81 + t20 * t49) - g(2) * (-t15 * t49 - t16 * t81) - g(3) * (t27 * t81 + t28 * t49);
t2 = -g(1) * (-t12 * t49 + t19 * t46) - g(2) * (-t10 * t49 + t16 * t46) - g(3) * (-t22 * t49 - t27 * t46);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * (t74 - t85) + t52, -g(1) * (t42 * t83 - t43 * t51) - g(2) * (-t42 * t51 - t43 * t83) + g(3) * t70, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t58, 0, -g(2) * t31 + (g(2) * t85 + t52) * pkin(2), 0, 0, 0, 0, 0, 0, -t57 * t50, t8, t58, -g(1) * t53 - g(2) * t56 - g(3) * t67, 0, 0, 0, 0, 0, 0, t4, t3, -t8, -g(1) * (t19 * t66 + t53) - g(2) * (t16 * t66 + t56) - g(3) * (-t27 * t66 + t67) 0, 0, 0, 0, 0, 0, t4, t3, -t8, -g(1) * (t19 * t65 + t20 * t73 + t54) - g(2) * (-t15 * t73 + t16 * t65 + t64) - g(3) * (-t27 * t65 + t28 * t73 + t77); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t59, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t59, -g(1) * (-pkin(4) * t11 + pkin(9) * t12) - g(2) * (-pkin(4) * t9 + pkin(9) * t10) - g(3) * (-pkin(4) * t21 + pkin(9) * t22) 0, 0, 0, 0, 0, 0, t6, -t5, -t59, -g(1) * (-t11 * t40 - t12 * t45) - g(2) * (-t10 * t45 - t40 * t9) - g(3) * (-t21 * t40 - t22 * t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60;];
taug_reg  = t7;

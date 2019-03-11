% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t41 = sin(qJ(1));
t43 = cos(qJ(1));
t17 = g(1) * t43 + g(2) * t41;
t40 = sin(qJ(2));
t83 = t17 * t40;
t28 = t40 * qJ(3);
t42 = cos(qJ(2));
t54 = t42 * pkin(2) + t28;
t36 = pkin(10) + qJ(5);
t26 = cos(t36);
t58 = t43 * t26;
t25 = sin(t36);
t68 = t41 * t25;
t7 = t40 * t58 - t68;
t75 = g(3) * t42;
t59 = t43 * t25;
t67 = t41 * t26;
t9 = t40 * t67 + t59;
t82 = -g(1) * t7 - g(2) * t9 + t26 * t75;
t12 = g(3) * t40 + t17 * t42;
t79 = g(1) * t41;
t37 = sin(pkin(10));
t74 = t37 * pkin(4);
t38 = cos(pkin(10));
t24 = t38 * pkin(4) + pkin(3);
t71 = t40 * t43;
t27 = qJ(6) + t36;
t22 = sin(t27);
t70 = t41 * t22;
t23 = cos(t27);
t69 = t41 * t23;
t66 = t41 * t37;
t65 = t41 * t38;
t39 = -pkin(8) - qJ(4);
t35 = -pkin(9) + t39;
t64 = t42 * t35;
t63 = t42 * t39;
t62 = t42 * t43;
t61 = t43 * t22;
t60 = t43 * t23;
t57 = t43 * t37;
t56 = t43 * t38;
t53 = t43 * pkin(1) + t41 * pkin(7);
t52 = qJ(3) * t42;
t51 = t42 * qJ(4);
t50 = t42 * t74;
t48 = t40 * t57;
t47 = pkin(2) * t62 + t43 * t28 + t53;
t46 = -g(2) * t43 + t79;
t45 = -pkin(1) - t54;
t32 = t43 * pkin(7);
t20 = t43 * t52;
t18 = t41 * t52;
t16 = pkin(5) * t25 + t74;
t15 = pkin(5) * t26 + t24;
t14 = t46 * t42;
t13 = t46 * t40;
t11 = -t75 + t83;
t10 = -t40 * t68 + t58;
t8 = t40 * t59 + t67;
t6 = -t40 * t70 + t60;
t5 = t40 * t69 + t61;
t4 = t40 * t61 + t69;
t3 = t40 * t60 - t70;
t2 = g(1) * t4 - g(2) * t6 - t22 * t75;
t1 = -g(1) * t3 - g(2) * t5 + t23 * t75;
t19 = [0, 0, 0, 0, 0, 0, t46, t17, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, -t17, -g(1) * (-t41 * pkin(1) + t32) - g(2) * t53, 0, 0, 0, 0, 0, 0, -t17, -t14, t13, -g(1) * t32 - g(2) * t47 - t45 * t79, 0, 0, 0, 0, 0, 0, -g(1) * (-t40 * t66 + t56) - g(2) * (t48 + t65) -g(1) * (-t40 * t65 - t57) - g(2) * (t40 * t56 - t66) t14, -g(1) * (t43 * pkin(3) + t32) - g(2) * (t43 * t51 + t47) + (-g(1) * (t45 - t51) - g(2) * pkin(3)) * t41, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, t14, -g(1) * (t43 * t24 + t32) - g(2) * (pkin(4) * t48 - t39 * t62 + t47) + (-g(1) * (-t40 * t74 + t45 + t63) - g(2) * t24) * t41, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t14, -g(1) * (t43 * t15 + t32) - g(2) * (t16 * t71 - t35 * t62 + t47) + (-g(1) * (-t40 * t16 + t45 + t64) - g(2) * t15) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -g(1) * (-pkin(2) * t71 + t20) - g(2) * (-t41 * t40 * pkin(2) + t18) - g(3) * t54, 0, 0, 0, 0, 0, 0, -t12 * t37, -t12 * t38, t11, -g(1) * t20 - g(2) * t18 - g(3) * (t51 + t54) + (pkin(2) + qJ(4)) * t83, 0, 0, 0, 0, 0, 0, -t12 * t25, -t12 * t26, t11, -g(1) * (t43 * t50 + t20) - g(2) * (t41 * t50 + t18) - g(3) * (t54 - t63) + (-g(3) * t74 + t17 * (pkin(2) - t39)) * t40, 0, 0, 0, 0, 0, 0, -t12 * t22, -t12 * t23, t11, -g(1) * (t16 * t62 + t20) - g(2) * (t41 * t42 * t16 + t18) - g(3) * (t54 - t64) + (-g(3) * t16 + t17 * (pkin(2) - t35)) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, g(1) * t8 - g(2) * t10 - t25 * t75, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t82 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t19;

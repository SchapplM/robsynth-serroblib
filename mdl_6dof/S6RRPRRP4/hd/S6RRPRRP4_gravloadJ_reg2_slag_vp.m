% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t45 = cos(qJ(4));
t31 = t45 * pkin(4) + pkin(3);
t39 = qJ(2) + pkin(10);
t33 = sin(t39);
t34 = cos(t39);
t48 = -pkin(9) - pkin(8);
t82 = t34 * t31 - t33 * t48;
t44 = sin(qJ(1));
t47 = cos(qJ(1));
t24 = g(1) * t47 + g(2) * t44;
t52 = -g(3) * t34 + t24 * t33;
t43 = sin(qJ(2));
t81 = pkin(2) * t43;
t78 = g(3) * t33;
t40 = qJ(4) + qJ(5);
t35 = sin(t40);
t76 = t33 * t35;
t36 = cos(t40);
t75 = t33 * t36;
t73 = t44 * t35;
t72 = t44 * t36;
t42 = sin(qJ(4));
t71 = t44 * t42;
t70 = t44 * t45;
t69 = t47 * t35;
t68 = t47 * t36;
t41 = -qJ(3) - pkin(7);
t67 = t47 * t41;
t66 = t47 * t42;
t65 = t47 * t45;
t64 = t34 * t66;
t10 = t34 * t73 + t68;
t11 = t34 * t72 - t69;
t63 = -t10 * pkin(5) + t11 * qJ(6);
t12 = t34 * t69 - t72;
t13 = t34 * t68 + t73;
t62 = -t12 * pkin(5) + t13 * qJ(6);
t46 = cos(qJ(2));
t37 = t46 * pkin(2);
t32 = t37 + pkin(1);
t25 = t47 * t32;
t61 = -t44 * t41 + t25;
t60 = t37 + t82;
t59 = t34 * pkin(3) + t33 * pkin(8);
t58 = g(1) * t10 - g(2) * t12;
t23 = g(1) * t44 - g(2) * t47;
t57 = -t34 * t48 - t81;
t56 = pkin(5) * t36 + qJ(6) * t35;
t15 = t34 * t71 + t65;
t1 = g(1) * t12 + g(2) * t10 + g(3) * t76;
t3 = g(1) * t13 + g(2) * t11 + g(3) * t75;
t53 = pkin(4) * t71 + t82 * t47 + t61;
t7 = t24 * t34 + t78;
t51 = -g(3) * t46 + t24 * t43;
t50 = pkin(4) * t66 - t67 + (-t32 - t82) * t44;
t29 = pkin(4) * t70;
t20 = qJ(6) * t75;
t18 = t34 * t65 + t71;
t17 = -t64 + t70;
t16 = -t34 * t70 + t66;
t14 = t23 * t33;
t6 = t52 * t36;
t5 = t52 * t35;
t4 = g(1) * t11 - g(2) * t13;
t2 = [0, 0, 0, 0, 0, 0, t23, t24, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t46, -t23 * t43, -t24, -g(1) * (-t44 * pkin(1) + t47 * pkin(7)) - g(2) * (t47 * pkin(1) + t44 * pkin(7)) 0, 0, 0, 0, 0, 0, t23 * t34, -t14, -t24, -g(1) * (-t44 * t32 - t67) - g(2) * t61, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, t14, -g(2) * t25 + (g(1) * t41 - g(2) * t59) * t47 + (-g(1) * (-t32 - t59) + g(2) * t41) * t44, 0, 0, 0, 0, 0, 0, t4, -t58, t14, -g(1) * t50 - g(2) * t53, 0, 0, 0, 0, 0, 0, t4, t14, t58, -g(1) * (-t11 * pkin(5) - t10 * qJ(6) + t50) - g(2) * (t13 * pkin(5) + t12 * qJ(6) + t53); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, g(3) * t43 + t24 * t46, 0, 0, 0, 0, 0, 0, 0, 0, t52, t7, 0, t51 * pkin(2), 0, 0, 0, 0, 0, 0, t52 * t45, -t52 * t42, -t7, -g(3) * (t37 + t59) + t24 * (pkin(3) * t33 - pkin(8) * t34 + t81) 0, 0, 0, 0, 0, 0, t6, -t5, -t7, -g(3) * t60 + t24 * (t31 * t33 - t57) 0, 0, 0, 0, 0, 0, t6, -t7, t5, -g(3) * (t56 * t34 + t60) + t24 * (-(-t31 - t56) * t33 - t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t17 + g(2) * t15 + t42 * t78, g(1) * t18 - g(2) * t16 + t45 * t78, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t29 + (g(2) * t65 + t7 * t42) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (-pkin(4) * t64 + t29 + t62) - g(2) * (-t15 * pkin(4) + t63) - g(3) * (t20 + (-pkin(4) * t42 - pkin(5) * t35) * t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t62 - g(2) * t63 - g(3) * (-pkin(5) * t76 + t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;

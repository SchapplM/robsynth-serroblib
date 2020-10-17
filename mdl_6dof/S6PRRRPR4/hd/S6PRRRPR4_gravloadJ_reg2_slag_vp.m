% Calculate inertial parameters regressor of gravitation load for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:46:25
% EndTime: 2019-05-05 07:46:27
% DurationCPUTime: 0.72s
% Computational Cost: add. (598->157), mult. (1165->242), div. (0->0), fcn. (1417->14), ass. (0->69)
t38 = sin(pkin(6));
t78 = g(3) * t38;
t40 = sin(qJ(4));
t77 = t40 * pkin(4);
t37 = qJ(4) + pkin(12);
t32 = sin(t37);
t19 = pkin(5) * t32 + t77;
t76 = pkin(8) + t19;
t34 = qJ(6) + t37;
t29 = sin(t34);
t44 = cos(qJ(3));
t75 = t29 * t44;
t30 = cos(t34);
t74 = t30 * t44;
t73 = t32 * t44;
t33 = cos(t37);
t72 = t33 * t44;
t42 = sin(qJ(2));
t71 = t38 * t42;
t45 = cos(qJ(2));
t70 = t38 * t45;
t69 = t40 * t42;
t68 = t40 * t44;
t43 = cos(qJ(4));
t67 = t43 * t44;
t66 = t44 * t45;
t39 = -qJ(5) - pkin(9);
t65 = pkin(2) * t70 + pkin(8) * t71;
t35 = t43 * pkin(4);
t20 = pkin(5) * t33 + t35;
t64 = cos(pkin(6));
t63 = cos(pkin(11));
t62 = sin(pkin(11));
t61 = pkin(8) + t77;
t60 = g(3) * t65;
t52 = t64 * t63;
t12 = t62 * t42 - t45 * t52;
t10 = t12 * pkin(2);
t13 = t42 * t52 + t62 * t45;
t59 = t13 * pkin(8) - t10;
t51 = t64 * t62;
t14 = t63 * t42 + t45 * t51;
t11 = t14 * pkin(2);
t15 = -t42 * t51 + t63 * t45;
t58 = t15 * pkin(8) - t11;
t57 = t38 * t63;
t56 = t38 * t62;
t41 = sin(qJ(3));
t55 = pkin(3) * t44 + pkin(9) * t41;
t18 = pkin(3) + t20;
t36 = -pkin(10) + t39;
t54 = t18 * t44 - t36 * t41;
t31 = t35 + pkin(3);
t53 = t31 * t44 - t39 * t41;
t16 = t41 * t71 - t64 * t44;
t6 = t13 * t41 + t44 * t57;
t8 = t15 * t41 - t44 * t56;
t50 = g(1) * t8 + g(2) * t6 + g(3) * t16;
t17 = t64 * t41 + t44 * t71;
t7 = t13 * t44 - t41 * t57;
t9 = t15 * t44 + t41 * t56;
t49 = g(1) * t9 + g(2) * t7 + g(3) * t17;
t48 = -g(1) * t14 - g(2) * t12 + g(3) * t70;
t47 = g(1) * t15 + g(2) * t13 + g(3) * t71;
t46 = -g(1) * (t14 * t43 - t9 * t40) - g(2) * (t12 * t43 - t7 * t40) - g(3) * (-t17 * t40 - t43 * t70);
t5 = t48 * t41;
t2 = -g(1) * (-t14 * t29 - t9 * t30) - g(2) * (-t12 * t29 - t7 * t30) - g(3) * (-t17 * t30 + t29 * t70);
t1 = -g(1) * (t14 * t30 - t9 * t29) - g(2) * (t12 * t30 - t7 * t29) - g(3) * (-t17 * t29 - t30 * t70);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t47, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * t44, t5, -t47, -g(1) * t58 - g(2) * t59 - t60, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t67 + t15 * t40) - g(2) * (-t12 * t67 + t13 * t40) - (t43 * t66 + t69) * t78, -g(1) * (t14 * t68 + t15 * t43) - g(2) * (t12 * t68 + t13 * t43) - (-t40 * t66 + t42 * t43) * t78, -t5, -g(1) * (-t55 * t14 + t58) - g(2) * (-t55 * t12 + t59) - g(3) * (t55 * t70 + t65) 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t72 + t15 * t32) - g(2) * (-t12 * t72 + t13 * t32) - (t32 * t42 + t33 * t66) * t78, -g(1) * (t14 * t73 + t15 * t33) - g(2) * (t12 * t73 + t13 * t33) - (-t32 * t66 + t33 * t42) * t78, -t5, -g(1) * (-t53 * t14 + t61 * t15 - t11) - g(2) * (-t53 * t12 + t61 * t13 - t10) - t60 - (pkin(4) * t69 + t53 * t45) * t78, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t74 + t15 * t29) - g(2) * (-t12 * t74 + t13 * t29) - (t29 * t42 + t30 * t66) * t78, -g(1) * (t14 * t75 + t15 * t30) - g(2) * (t12 * t75 + t13 * t30) - (-t29 * t66 + t30 * t42) * t78, -t5, -g(1) * (-t54 * t14 + t76 * t15 - t11) - g(2) * (-t54 * t12 + t76 * t13 - t10) - t60 - (t19 * t42 + t54 * t45) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t49, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t43, -t50 * t40, -t49, -g(1) * (-t8 * pkin(3) + t9 * pkin(9)) - g(2) * (-t6 * pkin(3) + t7 * pkin(9)) - g(3) * (-t16 * pkin(3) + t17 * pkin(9)) 0, 0, 0, 0, 0, 0, t50 * t33, -t50 * t32, -t49, -g(1) * (-t8 * t31 - t9 * t39) - g(2) * (-t6 * t31 - t7 * t39) - g(3) * (-t16 * t31 - t17 * t39) 0, 0, 0, 0, 0, 0, t50 * t30, -t50 * t29, -t49, -g(1) * (-t8 * t18 - t9 * t36) - g(2) * (-t6 * t18 - t7 * t36) - g(3) * (-t16 * t18 - t17 * t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -g(1) * (-t14 * t40 - t9 * t43) - g(2) * (-t12 * t40 - t7 * t43) - g(3) * (-t17 * t43 + t40 * t70) 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t33 - t9 * t32) - g(2) * (t12 * t33 - t7 * t32) - g(3) * (-t17 * t32 - t33 * t70) -g(1) * (-t14 * t32 - t9 * t33) - g(2) * (-t12 * t32 - t7 * t33) - g(3) * (-t17 * t33 + t32 * t70) 0, t46 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t14 * t20 - t9 * t19) - g(2) * (t12 * t20 - t7 * t19) - g(3) * (-t17 * t19 - t20 * t70); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;

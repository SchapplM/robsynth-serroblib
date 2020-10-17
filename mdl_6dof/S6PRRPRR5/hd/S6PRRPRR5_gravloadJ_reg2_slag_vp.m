% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:26:27
% EndTime: 2019-05-05 05:26:29
% DurationCPUTime: 0.64s
% Computational Cost: add. (563->142), mult. (1078->215), div. (0->0), fcn. (1311->14), ass. (0->67)
t37 = sin(pkin(6));
t76 = g(3) * t37;
t36 = sin(pkin(12));
t75 = t36 * pkin(4);
t35 = pkin(12) + qJ(5);
t30 = sin(t35);
t19 = pkin(5) * t30 + t75;
t74 = pkin(8) + t19;
t38 = cos(pkin(12));
t29 = t38 * pkin(4) + pkin(3);
t32 = qJ(6) + t35;
t27 = sin(t32);
t42 = cos(qJ(3));
t73 = t27 * t42;
t28 = cos(t32);
t72 = t28 * t42;
t71 = t30 * t42;
t31 = cos(t35);
t70 = t31 * t42;
t41 = sin(qJ(2));
t69 = t36 * t41;
t68 = t36 * t42;
t67 = t37 * t41;
t43 = cos(qJ(2));
t66 = t37 * t43;
t65 = t38 * t42;
t64 = t42 * t43;
t39 = -pkin(9) - qJ(4);
t63 = pkin(2) * t66 + pkin(8) * t67;
t62 = cos(pkin(6));
t61 = cos(pkin(11));
t60 = sin(pkin(11));
t59 = pkin(8) + t75;
t58 = g(3) * t63;
t50 = t62 * t61;
t12 = t60 * t41 - t43 * t50;
t10 = t12 * pkin(2);
t13 = t41 * t50 + t60 * t43;
t57 = t13 * pkin(8) - t10;
t49 = t62 * t60;
t14 = t61 * t41 + t43 * t49;
t11 = t14 * pkin(2);
t15 = -t41 * t49 + t61 * t43;
t56 = t15 * pkin(8) - t11;
t55 = t37 * t61;
t54 = t37 * t60;
t40 = sin(qJ(3));
t53 = pkin(3) * t42 + qJ(4) * t40;
t18 = pkin(5) * t31 + t29;
t34 = -pkin(10) + t39;
t52 = t18 * t42 - t34 * t40;
t51 = t29 * t42 - t39 * t40;
t16 = t40 * t67 - t62 * t42;
t6 = t13 * t40 + t42 * t55;
t8 = t15 * t40 - t42 * t54;
t48 = g(1) * t8 + g(2) * t6 + g(3) * t16;
t17 = t62 * t40 + t42 * t67;
t7 = t13 * t42 - t40 * t55;
t9 = t15 * t42 + t40 * t54;
t47 = g(1) * t9 + g(2) * t7 + g(3) * t17;
t46 = -g(1) * t14 - g(2) * t12 + g(3) * t66;
t45 = g(1) * t15 + g(2) * t13 + g(3) * t67;
t44 = -g(1) * (t14 * t31 - t9 * t30) - g(2) * (t12 * t31 - t7 * t30) - g(3) * (-t17 * t30 - t31 * t66);
t5 = t46 * t40;
t2 = -g(1) * (-t14 * t27 - t9 * t28) - g(2) * (-t12 * t27 - t7 * t28) - g(3) * (-t17 * t28 + t27 * t66);
t1 = -g(1) * (t14 * t28 - t9 * t27) - g(2) * (t12 * t28 - t7 * t27) - g(3) * (-t17 * t27 - t28 * t66);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t45, 0, 0, 0, 0, 0, 0, 0, 0, -t46 * t42, t5, -t45, -g(1) * t56 - g(2) * t57 - t58, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t65 + t15 * t36) - g(2) * (-t12 * t65 + t13 * t36) - (t38 * t64 + t69) * t76, -g(1) * (t14 * t68 + t15 * t38) - g(2) * (t12 * t68 + t13 * t38) - (-t36 * t64 + t38 * t41) * t76, -t5, -g(1) * (-t53 * t14 + t56) - g(2) * (-t53 * t12 + t57) - g(3) * (t53 * t66 + t63) 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t70 + t15 * t30) - g(2) * (-t12 * t70 + t13 * t30) - (t30 * t41 + t31 * t64) * t76, -g(1) * (t14 * t71 + t15 * t31) - g(2) * (t12 * t71 + t13 * t31) - (-t30 * t64 + t31 * t41) * t76, -t5, -g(1) * (-t51 * t14 + t59 * t15 - t11) - g(2) * (-t51 * t12 + t59 * t13 - t10) - t58 - (pkin(4) * t69 + t51 * t43) * t76, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t72 + t15 * t27) - g(2) * (-t12 * t72 + t13 * t27) - (t27 * t41 + t28 * t64) * t76, -g(1) * (t14 * t73 + t15 * t28) - g(2) * (t12 * t73 + t13 * t28) - (-t27 * t64 + t28 * t41) * t76, -t5, -g(1) * (-t52 * t14 + t74 * t15 - t11) - g(2) * (-t52 * t12 + t74 * t13 - t10) - t58 - (t19 * t41 + t52 * t43) * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t47, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t38, -t48 * t36, -t47, -g(1) * (-t8 * pkin(3) + t9 * qJ(4)) - g(2) * (-t6 * pkin(3) + t7 * qJ(4)) - g(3) * (-t16 * pkin(3) + t17 * qJ(4)) 0, 0, 0, 0, 0, 0, t48 * t31, -t48 * t30, -t47, -g(1) * (-t8 * t29 - t9 * t39) - g(2) * (-t6 * t29 - t7 * t39) - g(3) * (-t16 * t29 - t17 * t39) 0, 0, 0, 0, 0, 0, t48 * t28, -t48 * t27, -t47, -g(1) * (-t8 * t18 - t9 * t34) - g(2) * (-t6 * t18 - t7 * t34) - g(3) * (-t16 * t18 - t17 * t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -g(1) * (-t14 * t30 - t9 * t31) - g(2) * (-t12 * t30 - t7 * t31) - g(3) * (-t17 * t31 + t30 * t66) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t44 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;

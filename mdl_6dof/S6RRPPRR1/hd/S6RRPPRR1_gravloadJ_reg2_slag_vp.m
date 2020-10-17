% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRR1
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:33:33
% EndTime: 2019-05-06 09:33:34
% DurationCPUTime: 0.38s
% Computational Cost: add. (424->97), mult. (546->121), div. (0->0), fcn. (591->10), ass. (0->65)
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t71 = g(1) * t43;
t21 = g(2) * t40 + t71;
t35 = qJ(2) + pkin(10);
t31 = sin(t35);
t77 = t21 * t31;
t37 = sin(qJ(6));
t32 = cos(t35);
t38 = sin(qJ(5));
t67 = cos(qJ(5));
t13 = t31 * t38 + t32 * t67;
t60 = t31 * t67;
t14 = -t32 * t38 + t60;
t7 = t14 * t40;
t65 = t32 * t43;
t9 = t38 * t65 - t43 * t60;
t47 = g(1) * t9 - g(2) * t7 + g(3) * t13;
t76 = t47 * t37;
t41 = cos(qJ(6));
t75 = t47 * t41;
t26 = t31 * qJ(4);
t74 = -t32 * pkin(3) - t26;
t39 = sin(qJ(2));
t72 = pkin(2) * t39;
t69 = g(3) * t14;
t27 = t32 * pkin(4);
t36 = -qJ(3) - pkin(7);
t68 = pkin(8) + t36;
t64 = t43 * t36;
t63 = qJ(4) * t32;
t42 = cos(qJ(2));
t33 = t42 * pkin(2);
t30 = t33 + pkin(1);
t24 = t43 * t30;
t62 = pkin(3) * t65 + t43 * t26 + t24;
t61 = t33 - t74;
t59 = pkin(4) * t65 + t62;
t58 = t27 + t61;
t17 = t40 * t63;
t57 = -t40 * t72 + t17;
t19 = t43 * t63;
t56 = -t43 * t72 + t19;
t8 = t13 * t40;
t55 = t7 * pkin(5) + t8 * pkin(9);
t54 = -g(1) * t7 - g(2) * t9;
t10 = t13 * t43;
t53 = -t9 * pkin(5) + t10 * pkin(9);
t52 = -pkin(3) * t31 - t72;
t51 = -t13 * pkin(5) + t14 * pkin(9);
t20 = g(1) * t40 - g(2) * t43;
t50 = t8 * t37 - t43 * t41;
t49 = t43 * t37 + t8 * t41;
t48 = -t30 + t74;
t2 = g(1) * t10 + g(2) * t8 + t69;
t46 = -g(3) * t42 + t21 * t39;
t45 = (pkin(3) + pkin(4)) * t77;
t44 = (-g(1) * (t48 - t27) + g(2) * t68) * t40;
t12 = t20 * t32;
t11 = t20 * t31;
t6 = g(3) * t31 + t21 * t32;
t5 = -g(3) * t32 + t77;
t4 = t10 * t41 - t40 * t37;
t3 = -t10 * t37 - t40 * t41;
t1 = [0, 0, 0, 0, 0, 0, t20, t21, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t42, -t20 * t39, -t21, -g(1) * (-t40 * pkin(1) + t43 * pkin(7)) - g(2) * (t43 * pkin(1) + t40 * pkin(7)) 0, 0, 0, 0, 0, 0, t12, -t11, -t21, -g(1) * (-t40 * t30 - t64) - g(2) * (-t40 * t36 + t24) 0, 0, 0, 0, 0, 0, t12, -t21, t11, g(1) * t64 - g(2) * t62 + (-g(1) * t48 + g(2) * t36) * t40, 0, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t10, -t54, t21, -g(2) * t59 + t68 * t71 + t44, 0, 0, 0, 0, 0, 0, g(1) * t49 - g(2) * t4, -g(1) * t50 - g(2) * t3, t54, -g(1) * (-t8 * pkin(5) - t43 * pkin(8) + t7 * pkin(9) - t64) - g(2) * (t10 * pkin(5) + t9 * pkin(9) + t59) + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, g(3) * t39 + t21 * t42, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t46 * pkin(2), 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * (t52 * t43 + t19) - g(2) * (t52 * t40 + t17) - g(3) * t61, 0, 0, 0, 0, 0, 0, -t47, -t2, 0, -g(1) * t56 - g(2) * t57 - g(3) * t58 + t45, 0, 0, 0, 0, 0, 0, -t75, t76, t2, -g(1) * (-t53 + t56) - g(2) * (-t55 + t57) - g(3) * (-t51 + t58) + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t2, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t76, -t2, -g(1) * t53 - g(2) * t55 - g(3) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t50 + t37 * t69, g(1) * t4 + g(2) * t49 + t41 * t69, 0, 0;];
taug_reg  = t1;

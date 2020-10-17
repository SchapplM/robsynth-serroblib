% Calculate inertial parameters regressor of gravitation load for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:12:53
% EndTime: 2019-05-06 08:12:54
% DurationCPUTime: 0.49s
% Computational Cost: add. (366->97), mult. (511->131), div. (0->0), fcn. (543->10), ass. (0->68)
t36 = qJ(2) + pkin(9);
t33 = cos(t36);
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t22 = g(1) * t45 + g(2) * t42;
t32 = sin(t36);
t83 = t22 * t32;
t7 = -g(3) * t33 + t83;
t26 = t32 * qJ(4);
t84 = -t33 * pkin(3) - t26;
t41 = sin(qJ(2));
t82 = pkin(2) * t41;
t81 = g(1) * t42;
t39 = -qJ(3) - pkin(7);
t79 = g(2) * t39;
t77 = g(3) * t32;
t75 = t32 * t45;
t38 = cos(pkin(10));
t74 = t33 * t38;
t73 = t33 * t45;
t37 = sin(pkin(10));
t72 = t42 * t37;
t71 = t42 * t38;
t70 = t45 * t37;
t69 = t45 * t38;
t68 = t45 * t39;
t67 = qJ(4) * t33;
t66 = qJ(5) * t37;
t44 = cos(qJ(2));
t34 = t44 * pkin(2);
t31 = t34 + pkin(1);
t25 = t45 * t31;
t65 = pkin(3) * t73 + t45 * t26 + t25;
t64 = t34 - t84;
t63 = -pkin(3) - t66;
t17 = t42 * t67;
t62 = -t42 * t82 + t17;
t20 = t45 * t67;
t61 = -t45 * t82 + t20;
t60 = pkin(4) * t74 + t33 * t66 + t64;
t59 = -pkin(3) * t32 - t82;
t12 = t33 * t72 + t69;
t14 = t33 * t70 - t71;
t58 = g(1) * t12 - g(2) * t14;
t21 = -g(2) * t45 + t81;
t13 = t33 * t71 - t70;
t40 = sin(qJ(6));
t43 = cos(qJ(6));
t57 = t12 * t43 - t13 * t40;
t56 = t12 * t40 + t13 * t43;
t55 = t37 * t43 - t38 * t40;
t54 = t37 * t40 + t38 * t43;
t52 = -t31 + t84;
t51 = -t13 * pkin(4) - t12 * qJ(5) - t68;
t50 = g(3) * t55;
t15 = t33 * t69 + t72;
t49 = t15 * pkin(4) + t14 * qJ(5) + t65;
t47 = -g(3) * t44 + t22 * t41;
t46 = (-g(1) * t52 + t79) * t42;
t11 = -g(2) * t75 + t32 * t81;
t8 = t22 * t33 + t77;
t6 = t7 * t38;
t5 = t7 * t37;
t4 = g(1) * t13 - g(2) * t15;
t3 = t14 * t40 + t15 * t43;
t2 = t14 * t43 - t15 * t40;
t1 = -g(1) * t14 - g(2) * t12 - t37 * t77;
t9 = [0, 0, 0, 0, 0, 0, t21, t22, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t44, -t21 * t41, -t22, -g(1) * (-t42 * pkin(1) + t45 * pkin(7)) - g(2) * (t45 * pkin(1) + t42 * pkin(7)) 0, 0, 0, 0, 0, 0, t21 * t33, -t11, -t22, -g(1) * (-t42 * t31 - t68) - g(2) * (-t42 * t39 + t25) 0, 0, 0, 0, 0, 0, t4, -t58, t11, g(1) * t68 - g(2) * t65 + t46, 0, 0, 0, 0, 0, 0, t4, t11, t58, -g(1) * t51 - g(2) * t49 + t46, 0, 0, 0, 0, 0, 0, g(1) * t56 - g(2) * t3, g(1) * t57 - g(2) * t2, -t11, -g(1) * (-t13 * pkin(5) + t51) - g(2) * (t15 * pkin(5) - pkin(8) * t75 + t49) + (-g(1) * (t32 * pkin(8) + t52) + t79) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, g(3) * t41 + t22 * t44, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t47 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t45 * t59 + t20) - g(2) * (t42 * t59 + t17) - g(3) * t64, 0, 0, 0, 0, 0, 0, t6, -t8, t5, -g(1) * t61 - g(2) * t62 - g(3) * t60 + (pkin(4) * t38 - t63) * t83, 0, 0, 0, 0, 0, 0, t7 * t54, -t33 * t50 + t55 * t83, t8, -g(1) * (-pkin(8) * t73 + t61) - g(2) * (-t42 * t33 * pkin(8) + t62) - g(3) * (pkin(5) * t74 + t60) + (g(3) * pkin(8) + t22 * (-(-pkin(4) - pkin(5)) * t38 - t63)) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t57 - t32 * t50, g(1) * t3 + g(2) * t56 + t54 * t77, 0, 0;];
taug_reg  = t9;

% Calculate inertial parameters regressor of gravitation load for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:56:21
% EndTime: 2019-05-04 20:56:23
% DurationCPUTime: 0.53s
% Computational Cost: add. (867->119), mult. (2303->195), div. (0->0), fcn. (2987->16), ass. (0->72)
t69 = sin(pkin(13));
t70 = sin(pkin(12));
t54 = t70 * t69;
t73 = cos(pkin(13));
t74 = cos(pkin(12));
t61 = t74 * t73;
t76 = cos(pkin(6));
t45 = -t76 * t61 + t54;
t71 = sin(pkin(7));
t72 = sin(pkin(6));
t58 = t72 * t71;
t75 = cos(pkin(7));
t84 = t45 * t75 + t74 * t58;
t56 = t70 * t73;
t59 = t74 * t69;
t46 = t76 * t56 + t59;
t55 = t70 * t72;
t83 = t46 * t75 - t71 * t55;
t82 = t73 * t75 * t72 + t76 * t71;
t81 = cos(qJ(3));
t35 = qJ(5) + qJ(6);
t33 = sin(t35);
t40 = cos(qJ(4));
t80 = t33 * t40;
t34 = cos(t35);
t79 = t34 * t40;
t36 = sin(qJ(5));
t78 = t36 * t40;
t39 = cos(qJ(5));
t77 = t39 * t40;
t68 = pkin(5) * t36 + pkin(9);
t26 = t76 * t59 + t56;
t38 = sin(qJ(3));
t12 = t26 * t81 - t84 * t38;
t11 = t26 * t38 + t84 * t81;
t9 = t11 * pkin(3);
t67 = t12 * pkin(9) - t9;
t27 = -t76 * t54 + t61;
t13 = t27 * t38 + t83 * t81;
t10 = t13 * pkin(3);
t14 = t27 * t81 - t83 * t38;
t66 = t14 * pkin(9) - t10;
t57 = t72 * t69;
t18 = t38 * t57 - t82 * t81;
t17 = t18 * pkin(3);
t19 = t82 * t38 + t81 * t57;
t65 = t19 * pkin(9) - t17;
t37 = sin(qJ(4));
t64 = -pkin(4) * t40 - pkin(10) * t37;
t32 = t39 * pkin(5) + pkin(4);
t41 = -pkin(11) - pkin(10);
t63 = -t32 * t40 + t37 * t41;
t60 = t74 * t72;
t25 = -t73 * t58 + t76 * t75;
t15 = -t19 * t37 + t25 * t40;
t20 = t45 * t71 - t75 * t60;
t5 = -t12 * t37 + t20 * t40;
t21 = t46 * t71 + t75 * t55;
t7 = -t14 * t37 + t21 * t40;
t53 = g(1) * t7 + g(2) * t5 + g(3) * t15;
t16 = t19 * t40 + t25 * t37;
t6 = t12 * t40 + t20 * t37;
t8 = t14 * t40 + t21 * t37;
t52 = g(1) * t8 + g(2) * t6 + g(3) * t16;
t51 = g(1) * t13 + g(2) * t11 + g(3) * t18;
t50 = g(1) * t14 + g(2) * t12 + g(3) * t19;
t42 = -g(1) * (t13 * t39 - t8 * t36) - g(2) * (t11 * t39 - t6 * t36) - g(3) * (-t16 * t36 + t18 * t39);
t24 = -g(1) * t55 + g(2) * t60 - g(3) * t76;
t4 = t51 * t37;
t2 = -g(1) * (-t13 * t33 - t8 * t34) - g(2) * (-t11 * t33 - t6 * t34) - g(3) * (-t16 * t34 - t18 * t33);
t1 = -g(1) * (t13 * t34 - t8 * t33) - g(2) * (t11 * t34 - t6 * t33) - g(3) * (-t16 * t33 + t18 * t34);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t50, 0, 0, 0, 0, 0, 0, 0, 0, t51 * t40, -t4, -t50, -g(1) * t66 - g(2) * t67 - g(3) * t65, 0, 0, 0, 0, 0, 0, -g(1) * (-t13 * t77 + t14 * t36) - g(2) * (-t11 * t77 + t12 * t36) - g(3) * (-t18 * t77 + t19 * t36) -g(1) * (t13 * t78 + t14 * t39) - g(2) * (t11 * t78 + t12 * t39) - g(3) * (t18 * t78 + t19 * t39) t4, -g(1) * (t64 * t13 + t66) - g(2) * (t64 * t11 + t67) - g(3) * (t64 * t18 + t65) 0, 0, 0, 0, 0, 0, -g(1) * (-t13 * t79 + t14 * t33) - g(2) * (-t11 * t79 + t12 * t33) - g(3) * (-t18 * t79 + t19 * t33) -g(1) * (t13 * t80 + t14 * t34) - g(2) * (t11 * t80 + t12 * t34) - g(3) * (t18 * t80 + t19 * t34) t4, -g(1) * (t63 * t13 + t68 * t14 - t10) - g(2) * (t63 * t11 + t68 * t12 - t9) - g(3) * (t63 * t18 + t68 * t19 - t17); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t52, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * t39, t53 * t36, -t52, -g(1) * (t7 * pkin(4) + t8 * pkin(10)) - g(2) * (t5 * pkin(4) + t6 * pkin(10)) - g(3) * (t15 * pkin(4) + t16 * pkin(10)) 0, 0, 0, 0, 0, 0, -t53 * t34, t53 * t33, -t52, -g(1) * (t7 * t32 - t8 * t41) - g(2) * (t5 * t32 - t6 * t41) - g(3) * (t15 * t32 - t16 * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -g(1) * (-t13 * t36 - t8 * t39) - g(2) * (-t11 * t36 - t6 * t39) - g(3) * (-t16 * t39 - t18 * t36) 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t42 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;

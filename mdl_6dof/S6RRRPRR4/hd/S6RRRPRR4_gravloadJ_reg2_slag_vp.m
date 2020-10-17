% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:27:30
% EndTime: 2019-05-07 10:27:32
% DurationCPUTime: 0.52s
% Computational Cost: add. (540->106), mult. (527->141), div. (0->0), fcn. (517->12), ass. (0->74)
t46 = cos(pkin(11));
t31 = t46 * pkin(4) + pkin(3);
t43 = pkin(11) + qJ(5);
t36 = cos(t43);
t21 = pkin(5) * t36 + t31;
t44 = qJ(2) + qJ(3);
t38 = sin(t44);
t39 = cos(t44);
t47 = -pkin(9) - qJ(4);
t42 = -pkin(10) + t47;
t98 = t39 * t21 - t38 * t42;
t97 = t39 * t31 - t38 * t47;
t96 = t39 * pkin(3) + t38 * qJ(4);
t45 = sin(pkin(11));
t52 = -pkin(8) - pkin(7);
t66 = t45 * pkin(4) - t52;
t49 = sin(qJ(1));
t51 = cos(qJ(1));
t24 = g(1) * t51 + g(2) * t49;
t73 = t51 * t36;
t35 = sin(t43);
t80 = t49 * t35;
t15 = t39 * t80 + t73;
t74 = t51 * t35;
t79 = t49 * t36;
t17 = -t39 * t74 + t79;
t87 = g(3) * t38;
t95 = -g(1) * t17 + g(2) * t15 + t35 * t87;
t9 = -g(3) * t39 + t24 * t38;
t48 = sin(qJ(2));
t94 = pkin(2) * t48;
t93 = pkin(3) * t38;
t50 = cos(qJ(2));
t41 = t50 * pkin(2);
t34 = t41 + pkin(1);
t27 = t51 * t34;
t89 = g(2) * t27;
t37 = qJ(6) + t43;
t29 = sin(t37);
t82 = t49 * t29;
t30 = cos(t37);
t81 = t49 * t30;
t78 = t49 * t45;
t77 = t49 * t46;
t76 = t51 * t29;
t75 = t51 * t30;
t72 = t51 * t45;
t71 = t51 * t46;
t70 = pkin(5) * t35 + t66;
t68 = qJ(4) * t39;
t63 = -t93 - t94;
t62 = g(1) * t49 - g(2) * t51;
t59 = t21 * t38 + t39 * t42;
t57 = t31 * t38 + t39 * t47;
t53 = -g(3) * t50 + t24 * t48;
t26 = t51 * t68;
t25 = t49 * t68;
t19 = t62 * t38;
t18 = t39 * t73 + t80;
t16 = -t39 * t79 + t74;
t14 = t39 * t75 + t82;
t13 = -t39 * t76 + t81;
t12 = -t39 * t81 + t76;
t11 = t39 * t82 + t75;
t10 = t24 * t39 + t87;
t8 = t9 * t46;
t7 = t9 * t45;
t6 = t9 * t36;
t5 = t9 * t35;
t4 = t9 * t30;
t3 = t9 * t29;
t2 = g(1) * t14 - g(2) * t12 + t30 * t87;
t1 = -g(1) * t13 + g(2) * t11 + t29 * t87;
t20 = [0, 0, 0, 0, 0, 0, t62, t24, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t50, -t62 * t48, -t24, -g(1) * (-t49 * pkin(1) + t51 * pkin(7)) - g(2) * (t51 * pkin(1) + t49 * pkin(7)) 0, 0, 0, 0, 0, 0, t62 * t39, -t19, -t24, -g(1) * (-t49 * t34 - t51 * t52) - g(2) * (-t49 * t52 + t27) 0, 0, 0, 0, 0, 0, -g(1) * (-t39 * t77 + t72) - g(2) * (t39 * t71 + t78) -g(1) * (t39 * t78 + t71) - g(2) * (-t39 * t72 + t77) t19, -t89 + (g(1) * t52 - g(2) * t96) * t51 + (-g(1) * (-t34 - t96) + g(2) * t52) * t49, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, t19, -t89 + (-g(1) * t66 - g(2) * t97) * t51 + (-g(1) * (-t34 - t97) - g(2) * t66) * t49, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t19, -t89 + (-g(1) * t70 - g(2) * t98) * t51 + (-g(1) * (-t34 - t98) - g(2) * t70) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, g(3) * t48 + t24 * t50, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t53 * pkin(2), 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (t63 * t51 + t26) - g(2) * (t63 * t49 + t25) - g(3) * (t41 + t96) 0, 0, 0, 0, 0, 0, t6, -t5, -t10, -g(3) * (t41 + t97) + t24 * (t57 + t94) 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(3) * (t41 + t98) + t24 * (t59 + t94); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (-t51 * t93 + t26) - g(2) * (-t49 * t93 + t25) - g(3) * t96, 0, 0, 0, 0, 0, 0, t6, -t5, -t10, -g(3) * t97 + t24 * t57, 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(3) * t98 + t24 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, g(1) * t18 - g(2) * t16 + t36 * t87, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t95 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t20;

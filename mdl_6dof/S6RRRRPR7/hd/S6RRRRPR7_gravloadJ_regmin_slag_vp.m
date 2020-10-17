% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 21:20:21
% EndTime: 2019-05-07 21:20:23
% DurationCPUTime: 0.48s
% Computational Cost: add. (414->96), mult. (649->171), div. (0->0), fcn. (788->14), ass. (0->60)
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t45 = cos(qJ(2));
t46 = cos(qJ(1));
t59 = cos(pkin(6));
t55 = t46 * t59;
t18 = t41 * t42 - t45 * t55;
t39 = sin(qJ(6));
t43 = cos(qJ(6));
t19 = t41 * t55 + t42 * t45;
t37 = qJ(3) + qJ(4);
t32 = pkin(12) + t37;
t29 = sin(t32);
t30 = cos(t32);
t38 = sin(pkin(6));
t62 = t38 * t46;
t9 = -t19 * t30 + t29 * t62;
t77 = t18 * t43 + t39 * t9;
t76 = -t18 * t39 + t43 * t9;
t75 = g(1) * t46 + g(2) * t42;
t56 = t42 * t59;
t21 = -t41 * t56 + t45 * t46;
t33 = sin(t37);
t34 = cos(t37);
t64 = t38 * t42;
t12 = -t21 * t33 + t34 * t64;
t53 = t19 * t33 + t34 * t62;
t65 = t38 * t41;
t3 = -g(3) * (-t33 * t65 + t34 * t59) + g(2) * t53 - g(1) * t12;
t71 = g(3) * t38;
t67 = t30 * t39;
t66 = t30 * t43;
t44 = cos(qJ(3));
t63 = t38 * t44;
t61 = t39 * t45;
t60 = t43 * t45;
t24 = pkin(3) * t44 + pkin(4) * t34;
t58 = t19 * t34 - t33 * t62;
t40 = sin(qJ(3));
t57 = t19 * t44 - t40 * t62;
t20 = t41 * t46 + t45 * t56;
t54 = g(1) * t18 - g(2) * t20;
t52 = t19 * t40 + t44 * t62;
t51 = g(1) * (-t21 * t29 + t30 * t64) + g(2) * (-t19 * t29 - t30 * t62) + g(3) * (-t29 * t65 + t30 * t59);
t49 = -g(1) * t20 - g(2) * t18 + t45 * t71;
t48 = g(1) * t21 + g(2) * t19 + g(3) * t65;
t36 = -qJ(5) - pkin(10) - pkin(9);
t23 = pkin(3) * t40 + pkin(4) * t33;
t22 = pkin(2) + t24;
t17 = t29 * t59 + t30 * t65;
t15 = t21 * t44 + t40 * t64;
t14 = -t21 * t40 + t42 * t63;
t13 = t21 * t34 + t33 * t64;
t11 = t21 * t30 + t29 * t64;
t6 = t11 * t43 + t20 * t39;
t5 = -t11 * t39 + t20 * t43;
t4 = g(1) * t13 + g(2) * t58 - g(3) * (-t33 * t59 - t34 * t65);
t2 = t51 * t43;
t1 = t51 * t39;
t7 = [0, g(1) * t42 - g(2) * t46, t75, 0, 0, 0, 0, 0, g(1) * t19 - g(2) * t21, -t54, 0, 0, 0, 0, 0, g(1) * t57 - g(2) * t15, -g(1) * t52 - g(2) * t14, 0, 0, 0, 0, 0, g(1) * t58 - g(2) * t13, -g(1) * t53 - g(2) * t12, t54, -g(1) * (-t42 * pkin(1) + t18 * t36 - t19 * t22) - g(2) * (t46 * pkin(1) - t20 * t36 + t21 * t22) - t75 * t38 * (pkin(8) + t23) 0, 0, 0, 0, 0, -g(1) * t76 - g(2) * t6, g(1) * t77 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, -t49, t48, 0, 0, 0, 0, 0, -t49 * t44, t49 * t40, 0, 0, 0, 0, 0, -t49 * t34, t49 * t33, -t48, -g(1) * (-t20 * t22 - t21 * t36) - g(2) * (-t18 * t22 - t19 * t36) - (t22 * t45 - t36 * t41) * t71, 0, 0, 0, 0, 0, -g(1) * (-t20 * t66 + t21 * t39) - g(2) * (-t18 * t66 + t19 * t39) - (t30 * t60 + t39 * t41) * t71, -g(1) * (t20 * t67 + t21 * t43) - g(2) * (t18 * t67 + t19 * t43) - (-t30 * t61 + t41 * t43) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t14 + g(2) * t52 - g(3) * (-t40 * t65 + t44 * t59) g(1) * t15 + g(2) * t57 - g(3) * (-t40 * t59 - t41 * t63) 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (-t21 * t23 + t24 * t64) - g(2) * (-t19 * t23 - t24 * t62) - g(3) * (-t23 * t65 + t24 * t59) 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t3 * pkin(4), 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t77 - g(3) * (-t17 * t39 - t38 * t60) g(1) * t6 - g(2) * t76 - g(3) * (-t17 * t43 + t38 * t61);];
taug_reg  = t7;

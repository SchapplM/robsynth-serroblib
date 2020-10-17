% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR9
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
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:22:28
% EndTime: 2019-05-07 22:22:31
% DurationCPUTime: 0.58s
% Computational Cost: add. (601->114), mult. (978->193), div. (0->0), fcn. (1212->14), ass. (0->65)
t45 = sin(qJ(2));
t46 = sin(qJ(1));
t48 = cos(qJ(2));
t62 = cos(pkin(6));
t74 = cos(qJ(1));
t55 = t62 * t74;
t25 = t45 * t55 + t46 * t48;
t40 = qJ(3) + qJ(4);
t37 = sin(t40);
t38 = cos(t40);
t42 = sin(pkin(6));
t61 = t42 * t74;
t14 = t25 * t38 - t37 * t61;
t24 = t46 * t45 - t48 * t55;
t39 = pkin(12) + qJ(6);
t35 = sin(t39);
t36 = cos(t39);
t83 = t14 * t35 - t24 * t36;
t82 = t14 * t36 + t24 * t35;
t58 = t46 * t62;
t27 = -t45 * t58 + t74 * t48;
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t64 = t42 * t47;
t19 = -t27 * t44 + t46 * t64;
t53 = t25 * t44 + t47 * t61;
t66 = t42 * t45;
t81 = g(2) * t53 - g(3) * (-t44 * t66 + t62 * t47) - g(1) * t19;
t26 = t74 * t45 + t48 * t58;
t80 = -g(1) * t26 - g(2) * t24;
t76 = g(2) * t46;
t75 = g(3) * t42;
t71 = t35 * t38;
t70 = t36 * t38;
t41 = sin(pkin(12));
t69 = t38 * t41;
t43 = cos(pkin(12));
t68 = t38 * t43;
t67 = t38 * t48;
t65 = t42 * t46;
t63 = t42 * t48;
t60 = t44 * t74;
t59 = t25 * t47 - t42 * t60;
t13 = t25 * t37 + t38 * t61;
t17 = t27 * t37 - t38 * t65;
t57 = -g(1) * t13 + g(2) * t17;
t56 = g(1) * t27 + g(2) * t25;
t22 = t37 * t66 - t62 * t38;
t5 = g(1) * t17 + g(2) * t13 + g(3) * t22;
t18 = t27 * t38 + t37 * t65;
t23 = t62 * t37 + t38 * t66;
t7 = g(1) * t18 + g(2) * t14 + g(3) * t23;
t52 = g(3) * t63 + t80;
t51 = -g(1) * (-t17 * pkin(4) + t18 * qJ(5)) - g(2) * (-t13 * pkin(4) + t14 * qJ(5)) - g(3) * (-t22 * pkin(4) + t23 * qJ(5));
t49 = -pkin(10) - pkin(9);
t34 = t47 * pkin(3) + pkin(2);
t20 = t27 * t47 + t44 * t65;
t10 = t52 * t37;
t9 = t18 * t36 + t26 * t35;
t8 = -t18 * t35 + t26 * t36;
t4 = t5 * t43;
t3 = t5 * t41;
t2 = t5 * t36;
t1 = t5 * t35;
t6 = [0, g(1) * t46 - g(2) * t74, g(1) * t74 + t76, 0, 0, 0, 0, 0, g(1) * t25 - g(2) * t27, -g(1) * t24 + g(2) * t26, 0, 0, 0, 0, 0, g(1) * t59 - g(2) * t20, -g(1) * t53 - g(2) * t19, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t18, t57, -g(1) * (-t14 * t43 - t24 * t41) - g(2) * (t18 * t43 + t26 * t41) -g(1) * (t14 * t41 - t24 * t43) - g(2) * (-t18 * t41 + t26 * t43) -t57, -g(1) * (-t46 * pkin(1) - pkin(4) * t14 - qJ(5) * t13 + t24 * t49 - t25 * t34) - g(2) * (t74 * pkin(1) + t18 * pkin(4) + t17 * qJ(5) - t26 * t49 + t27 * t34) + (-g(1) * (pkin(3) * t60 + t74 * pkin(8)) - (pkin(3) * t44 + pkin(8)) * t76) * t42, 0, 0, 0, 0, 0, g(1) * t82 - g(2) * t9, -g(1) * t83 - g(2) * t8; 0, 0, 0, 0, 0, 0, 0, 0, -t52, g(3) * t66 + t56, 0, 0, 0, 0, 0, -t52 * t47, t52 * t44, 0, 0, 0, 0, 0, -t52 * t38, t10, -g(1) * (-t26 * t68 + t27 * t41) - g(2) * (-t24 * t68 + t25 * t41) - (t41 * t45 + t43 * t67) * t75, -g(1) * (t26 * t69 + t27 * t43) - g(2) * (t24 * t69 + t25 * t43) - (-t41 * t67 + t43 * t45) * t75, -t10 (t45 * t75 + t56) * t49 + (-t48 * t75 - t80) * (pkin(4) * t38 + qJ(5) * t37 + t34) 0, 0, 0, 0, 0, -g(1) * (-t26 * t70 + t27 * t35) - g(2) * (-t24 * t70 + t25 * t35) - (t35 * t45 + t36 * t67) * t75, -g(1) * (t26 * t71 + t27 * t36) - g(2) * (t24 * t71 + t25 * t36) - (-t35 * t67 + t36 * t45) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, g(1) * t20 + g(2) * t59 - g(3) * (-t62 * t44 - t45 * t64) 0, 0, 0, 0, 0, t5, t7, t4, -t3, -t7, t81 * pkin(3) + t51, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t7, t4, -t3, -t7, t51, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t83 - g(3) * (-t23 * t35 - t36 * t63) g(1) * t9 + g(2) * t82 - g(3) * (-t23 * t36 + t35 * t63);];
taug_reg  = t6;

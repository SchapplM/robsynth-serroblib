% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR15_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:25
% EndTime: 2024-09-27 22:26:25
% DurationCPUTime: 0.36s
% Computational Cost: add. (554->82), mult. (494->142), div. (0->0), fcn. (582->22), ass. (0->83)
t53 = sin(qJ(1));
t85 = -t53 / 0.2e1;
t48 = sin(pkin(5));
t84 = g(3) * t48;
t47 = sin(pkin(6));
t83 = t47 * t48;
t50 = cos(pkin(5));
t82 = t47 * t50;
t49 = cos(pkin(6));
t51 = sin(qJ(5));
t81 = t49 * t51;
t54 = cos(qJ(5));
t80 = t49 * t54;
t79 = t53 * t51;
t52 = sin(qJ(2));
t78 = t53 * t52;
t55 = cos(qJ(2));
t77 = t53 * t55;
t76 = t54 * t53;
t56 = cos(qJ(1));
t75 = t56 * t51;
t74 = t56 * t52;
t73 = t56 * t54;
t72 = t56 * t55;
t46 = qJ(2) + qJ(3);
t71 = t53 * t83;
t70 = t56 * t83;
t69 = t50 * t75;
t68 = t50 * t79;
t67 = t50 * t76;
t66 = t50 * t73;
t41 = pkin(5) + t46;
t37 = qJ(4) + t41;
t65 = sin(t37) / 0.2e1;
t42 = pkin(5) - t46;
t64 = -qJ(4) + t42;
t59 = sin(t64);
t21 = t65 - t59 / 0.2e1;
t45 = qJ(4) + t46;
t40 = cos(t45);
t63 = t56 * t21 + t53 * t40;
t62 = t53 * t21 - t56 * t40;
t33 = sin(t41);
t34 = sin(t42);
t27 = t33 - t34;
t44 = cos(t46);
t61 = t53 * t44 + t56 * t27 / 0.2e1;
t60 = t27 * t85 + t56 * t44;
t13 = t49 * t68 - t73;
t19 = t49 * t75 + t67;
t39 = sin(t45);
t58 = -t13 * t40 - t19 * t39 + t51 * t71;
t16 = t49 * t66 - t79;
t18 = t49 * t76 + t69;
t57 = -t16 * t40 + t18 * t39 + t54 * t70;
t43 = sin(t46);
t36 = cos(t42);
t35 = cos(t41);
t32 = cos(t64);
t31 = cos(t37) / 0.2e1;
t28 = t36 + t35;
t26 = -t50 * t78 + t72;
t25 = -t50 * t77 - t74;
t24 = -t50 * t74 - t77;
t23 = -t50 * t72 + t78;
t22 = t32 / 0.2e1 + t31;
t20 = -t49 * t73 + t68;
t17 = t49 * t79 - t66;
t15 = t49 * t69 + t76;
t14 = -t49 * t67 - t75;
t12 = t28 * t85 - t56 * t43;
t11 = t53 * t43 - t56 * t28 / 0.2e1;
t10 = -t53 * t22 - t56 * t39;
t9 = -t56 * t22 + t53 * t39;
t8 = -t15 * t40 + t17 * t39 + t51 * t70;
t7 = t14 * t40 + t20 * t39 + t54 * t71;
t6 = g(1) * t60 + g(2) * t61 - g(3) * (t35 / 0.2e1 - t36 / 0.2e1);
t5 = -g(1) * t12 + g(2) * t11 - g(3) * (t33 / 0.2e1 + t34 / 0.2e1);
t4 = -g(1) * t62 + g(2) * t63 - g(3) * (t31 - t32 / 0.2e1);
t3 = -g(1) * t10 + g(2) * t9 - g(3) * (t65 + t59 / 0.2e1);
t2 = -g(1) * (-t14 * t39 + t20 * t40) - g(2) * (-t16 * t39 - t18 * t40) - (-t39 * t80 - t40 * t51) * t84;
t1 = -g(1) * (t13 * t39 - t19 * t40) - g(2) * (-t15 * t39 - t17 * t40) - (-t39 * t81 + t40 * t54) * t84;
t29 = [0, g(1) * t53 - g(2) * t56, g(1) * t56 + g(2) * t53, 0, 0, 0, 0, 0, -g(1) * t24 - g(2) * t26, -g(1) * t23 - g(2) * t25, 0, 0, 0, 0, 0, g(1) * t61 - g(2) * t60, -g(1) * t11 - g(2) * t12, 0, 0, 0, 0, 0, g(1) * t63 + g(2) * t62, -g(1) * t9 - g(2) * t10, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t58, -g(1) * t57 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t25 + g(2) * t23 - t55 * t84, g(1) * t26 - g(2) * t24 + t52 * t84, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t57 - g(3) * (t54 * t82 + (-t39 * t51 + t40 * t80) * t48), g(1) * t58 - g(2) * t8 - g(3) * (-t51 * t82 + (-t39 * t54 - t40 * t81) * t48);];
taug_reg = t29;

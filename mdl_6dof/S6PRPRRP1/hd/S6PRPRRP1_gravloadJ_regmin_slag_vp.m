% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [6x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:24:18
% EndTime: 2021-01-16 01:24:20
% DurationCPUTime: 0.38s
% Computational Cost: add. (539->91), mult. (642->149), div. (0->0), fcn. (761->18), ass. (0->76)
t44 = qJ(2) + pkin(11);
t43 = cos(t44);
t49 = cos(pkin(10));
t33 = t49 * t43;
t42 = sin(t44);
t46 = sin(pkin(10));
t50 = cos(pkin(6));
t73 = t46 * t50;
t17 = t42 * t73 - t33;
t53 = sin(qJ(4));
t47 = sin(pkin(6));
t56 = cos(qJ(4));
t71 = t47 * t56;
t11 = t17 * t53 + t46 * t71;
t72 = t47 * t53;
t22 = t42 * t72 - t50 * t56;
t32 = t46 * t43;
t69 = t49 * t50;
t20 = t42 * t69 + t32;
t9 = t20 * t53 + t49 * t71;
t60 = g(1) * t11 - g(2) * t9 - g(3) * t22;
t81 = t46 / 0.2e1;
t80 = -t49 / 0.2e1;
t78 = g(3) * t47;
t40 = pkin(6) + t44;
t34 = sin(t40);
t41 = pkin(6) - t44;
t35 = sin(t41);
t26 = -t34 + t35;
t14 = t26 * t80 + t32;
t52 = sin(qJ(5));
t77 = t14 * t52;
t16 = t26 * t81 + t33;
t76 = t16 * t52;
t36 = cos(t40);
t37 = cos(t41);
t25 = t37 / 0.2e1 - t36 / 0.2e1;
t75 = t25 * t52;
t74 = t46 * t42;
t70 = t49 * t42;
t54 = sin(qJ(2));
t68 = t50 * t54;
t57 = cos(qJ(2));
t67 = t50 * t57;
t66 = t52 * t56;
t55 = cos(qJ(5));
t65 = t55 * t56;
t64 = t43 * t71;
t8 = t17 * t56 - t46 * t72;
t39 = t55 * pkin(5) + pkin(4);
t51 = -qJ(6) - pkin(9);
t62 = t56 * t39 - t53 * t51;
t10 = t20 * t56 - t49 * t72;
t23 = t42 * t71 + t50 * t53;
t61 = g(1) * t8 - g(2) * t10 - g(3) * t23;
t18 = t43 * t73 + t70;
t19 = t43 * t69 - t74;
t59 = g(1) * t18 - g(2) * t19 - t43 * t78;
t58 = -g(1) * (-t46 * t67 - t49 * t54) - g(2) * (-t46 * t54 + t49 * t67) - t57 * t78;
t27 = t36 + t37;
t13 = t27 * t80 + t74;
t15 = t27 * t81 + t70;
t24 = -t35 / 0.2e1 - t34 / 0.2e1;
t1 = -g(1) * (t15 * t55 + t52 * t8) - g(2) * (-t10 * t52 + t13 * t55) - g(3) * (-t23 * t52 + t24 * t55);
t48 = cos(pkin(11));
t45 = sin(pkin(11));
t29 = -t45 * pkin(3) + t48 * pkin(8);
t28 = t48 * pkin(3) + t45 * pkin(8) + pkin(2);
t21 = -g(3) * t50 + (-g(1) * t46 + g(2) * t49) * t47;
t7 = t59 * t53;
t6 = t60 * t55;
t5 = t60 * t52;
t4 = -g(1) * (-t18 * t65 + t76) - g(2) * (t19 * t65 + t77) - g(3) * (t55 * t64 + t75);
t3 = -g(1) * (t16 * t55 + t18 * t66) - g(2) * (t14 * t55 - t19 * t66) - g(3) * (t25 * t55 - t52 * t64);
t2 = -g(1) * (-t15 * t52 + t55 * t8) - g(2) * (-t10 * t55 - t13 * t52) - g(3) * (-t23 * t55 - t24 * t52);
t12 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t58, -g(1) * (t46 * t68 - t49 * t57) - g(2) * (-t46 * t57 - t49 * t68) + t54 * t78, t58 * pkin(2), 0, 0, 0, 0, 0, t59 * t56, -t7, 0, 0, 0, 0, 0, t4, t3, t4, t3, t7, -g(1) * (pkin(5) * t76 - (t49 * t28 + t29 * t73) * t54 + (-t28 * t73 + t49 * t29) * t57 - t62 * t18) - g(2) * (pkin(5) * t77 - (t46 * t28 - t29 * t69) * t54 + (t28 * t69 + t46 * t29) * t57 + t62 * t19) - g(3) * pkin(5) * t75 - (t28 * t57 + t29 * t54 + t62 * t43) * t78; 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t61, 0, 0, 0, 0, 0, -t6, t5, -t6, t5, t61, -g(1) * (t11 * t39 + t8 * t51) - g(2) * (-t10 * t51 - t9 * t39) - g(3) * (-t22 * t39 - t23 * t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60;];
taug_reg = t12;

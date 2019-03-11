% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t40 = sin(pkin(7));
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t43 = cos(pkin(12));
t68 = cos(pkin(6));
t62 = t43 * t68;
t66 = sin(pkin(12));
t51 = t66 * t46 - t48 * t62;
t67 = cos(pkin(7));
t41 = sin(pkin(6));
t71 = t41 * t43;
t81 = t40 * t71 + t51 * t67;
t57 = t68 * t66;
t52 = t43 * t46 + t48 * t57;
t63 = t41 * t66;
t80 = -t40 * t63 + t52 * t67;
t79 = pkin(9) * t40;
t78 = cos(qJ(3));
t38 = pkin(13) + qJ(5);
t36 = sin(t38);
t77 = t36 * t40;
t37 = cos(t38);
t76 = t37 * t40;
t44 = sin(qJ(6));
t75 = t37 * t44;
t47 = cos(qJ(6));
t74 = t37 * t47;
t39 = sin(pkin(13));
t73 = t39 * t40;
t42 = cos(pkin(13));
t72 = t40 * t42;
t70 = t41 * t46;
t69 = t41 * t48;
t64 = t40 * t70;
t45 = sin(qJ(3));
t61 = t45 * t67;
t60 = t68 * t40;
t58 = t67 * t78;
t29 = t46 * t62 + t66 * t48;
t10 = t29 * t78 - t81 * t45;
t30 = t43 * t48 - t46 * t57;
t12 = t30 * t78 - t80 * t45;
t19 = t45 * t60 + (t78 * t46 + t48 * t61) * t41;
t20 = t51 * t40 - t67 * t71;
t21 = t52 * t40 + t67 * t63;
t28 = -t40 * t69 + t68 * t67;
t56 = g(1) * (-t12 * t36 + t21 * t37) + g(2) * (-t10 * t36 + t20 * t37) + g(3) * (-t19 * t36 + t28 * t37);
t11 = t30 * t45 + t80 * t78;
t18 = t45 * t70 - t58 * t69 - t78 * t60;
t9 = t29 * t45 + t81 * t78;
t55 = g(1) * t11 + g(2) * t9 + g(3) * t18;
t54 = g(1) * t12 + g(2) * t10 + g(3) * t19;
t13 = t29 * t58 - t51 * t45;
t15 = t30 * t58 - t52 * t45;
t26 = (t45 * t48 + t46 * t58) * t41;
t53 = g(1) * t15 + g(2) * t13 + g(3) * t26;
t27 = (-t46 * t61 + t78 * t48) * t41;
t17 = t27 * t37 + t36 * t64;
t16 = -t30 * t61 - t52 * t78;
t14 = -t29 * t61 - t51 * t78;
t8 = t19 * t37 + t28 * t36;
t6 = t16 * t37 + t30 * t77;
t5 = t14 * t37 + t29 * t77;
t4 = t12 * t37 + t21 * t36;
t2 = t10 * t37 + t20 * t36;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t52 + g(2) * t51 - g(3) * t69, g(1) * t30 + g(2) * t29 + g(3) * t70, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14 - g(3) * t27, t53, -g(1) * (t16 * t42 + t30 * t73) - g(2) * (t14 * t42 + t29 * t73) - g(3) * (t27 * t42 + t39 * t64) -g(1) * (-t16 * t39 + t30 * t72) - g(2) * (-t14 * t39 + t29 * t72) - g(3) * (-t27 * t39 + t42 * t64) -t53, -g(1) * (-t52 * pkin(2) + t16 * pkin(3) + t15 * qJ(4) + t30 * t79) - g(2) * (-t51 * pkin(2) + t14 * pkin(3) + t13 * qJ(4) + t29 * t79) - g(3) * (t27 * pkin(3) + t26 * qJ(4) + (pkin(2) * t48 + t46 * t79) * t41) 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t5 - g(3) * t17, -g(1) * (-t16 * t36 + t30 * t76) - g(2) * (-t14 * t36 + t29 * t76) - g(3) * (-t27 * t36 + t37 * t64) 0, 0, 0, 0, 0, -g(1) * (t15 * t44 + t6 * t47) - g(2) * (t13 * t44 + t5 * t47) - g(3) * (t17 * t47 + t26 * t44) -g(1) * (t15 * t47 - t6 * t44) - g(2) * (t13 * t47 - t5 * t44) - g(3) * (-t17 * t44 + t26 * t47); 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t54, t55 * t42, -t55 * t39, -t54, -g(1) * (-t11 * pkin(3) + t12 * qJ(4)) - g(2) * (-t9 * pkin(3) + t10 * qJ(4)) - g(3) * (-t18 * pkin(3) + t19 * qJ(4)) 0, 0, 0, 0, 0, t55 * t37, -t55 * t36, 0, 0, 0, 0, 0, -g(1) * (-t11 * t74 + t12 * t44) - g(2) * (t10 * t44 - t9 * t74) - g(3) * (-t18 * t74 + t19 * t44) -g(1) * (t11 * t75 + t12 * t47) - g(2) * (t10 * t47 + t9 * t75) - g(3) * (t18 * t75 + t19 * t47); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, g(1) * t4 + g(2) * t2 + g(3) * t8, 0, 0, 0, 0, 0, -t56 * t47, t56 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t47 - t4 * t44) - g(2) * (-t2 * t44 + t9 * t47) - g(3) * (t18 * t47 - t8 * t44) -g(1) * (-t11 * t44 - t4 * t47) - g(2) * (-t2 * t47 - t9 * t44) - g(3) * (-t18 * t44 - t8 * t47);];
taug_reg  = t1;

% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRR3
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
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t44 = sin(pkin(12));
t47 = cos(pkin(12));
t53 = sin(qJ(2));
t49 = cos(pkin(6));
t57 = cos(qJ(2));
t73 = t49 * t57;
t34 = -t44 * t53 + t47 * t73;
t36 = -t44 * t73 - t47 * t53;
t45 = sin(pkin(7));
t48 = cos(pkin(7));
t46 = sin(pkin(6));
t77 = t46 * t57;
t80 = t46 * t47;
t83 = t44 * t46;
t85 = g(1) * (t36 * t48 + t45 * t83) - g(2) * (-t34 * t48 + t45 * t80) + g(3) * (t45 * t49 + t48 * t77);
t84 = g(3) * t46;
t51 = sin(qJ(5));
t82 = t45 * t51;
t55 = cos(qJ(5));
t81 = t45 * t55;
t79 = t46 * t48;
t78 = t46 * t53;
t52 = sin(qJ(3));
t76 = t48 * t52;
t56 = cos(qJ(3));
t75 = t48 * t56;
t74 = t49 * t53;
t50 = sin(qJ(6));
t72 = t50 * t55;
t54 = cos(qJ(6));
t71 = t54 * t55;
t70 = cos(pkin(13));
t69 = t45 * t78;
t43 = sin(pkin(13));
t39 = -t56 * t43 - t52 * t70;
t65 = -t52 * t43 + t56 * t70;
t25 = -t34 * t45 - t47 * t79;
t26 = -t36 * t45 + t44 * t79;
t33 = -t45 * t77 + t49 * t48;
t29 = t39 * t45;
t31 = t39 * t48;
t59 = -t49 * t29 + (-t31 * t57 + t53 * t65) * t46;
t37 = -t44 * t74 + t47 * t57;
t60 = -t29 * t83 - t36 * t31 + t37 * t65;
t35 = t44 * t57 + t47 * t74;
t61 = t29 * t80 - t34 * t31 + t35 * t65;
t64 = g(1) * (t26 * t55 - t51 * t60) + g(2) * (t25 * t55 - t51 * t61) + g(3) * (t33 * t55 - t51 * t59);
t28 = t65 * t45;
t30 = t65 * t48;
t11 = t28 * t83 + t36 * t30 + t37 * t39;
t16 = t49 * t28 + (t30 * t57 + t39 * t53) * t46;
t8 = -t28 * t80 + t34 * t30 + t35 * t39;
t63 = g(1) * t11 + g(2) * t8 + g(3) * t16;
t62 = g(1) * t37 + g(2) * t35 + g(3) * t78;
t58 = t62 * t52 - t85 * t56;
t42 = t56 * pkin(3) + pkin(2);
t32 = pkin(3) * t76 + (-pkin(9) - qJ(4)) * t45;
t24 = (t31 * t53 + t57 * t65) * t46;
t23 = (t30 * t53 - t39 * t57) * t46;
t22 = t24 * t55 + t51 * t69;
t21 = t37 * t31 + t36 * t65;
t20 = t37 * t30 - t36 * t39;
t19 = t35 * t31 + t34 * t65;
t18 = t35 * t30 - t34 * t39;
t14 = t33 * t51 + t55 * t59;
t6 = t21 * t55 + t37 * t82;
t5 = t19 * t55 + t35 * t82;
t4 = t26 * t51 + t55 * t60;
t2 = t25 * t51 + t55 * t61;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(1) * t36 - g(2) * t34 - g(3) * t77, t62, 0, 0, 0, 0, 0, -g(1) * (t36 * t56 - t37 * t76) - g(2) * (t34 * t56 - t35 * t76) - (-t53 * t76 + t56 * t57) * t84, -g(1) * (-t36 * t52 - t37 * t75) - g(2) * (-t34 * t52 - t35 * t75) - (-t52 * t57 - t53 * t75) * t84, -t62 * t45, -g(1) * (-t37 * t32 + t36 * t42) - g(2) * (-t35 * t32 + t34 * t42) - (-t32 * t53 + t42 * t57) * t84, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t5 - g(3) * t22, -g(1) * (-t21 * t51 + t37 * t81) - g(2) * (-t19 * t51 + t35 * t81) - g(3) * (-t24 * t51 + t55 * t69) 0, 0, 0, 0, 0, -g(1) * (t20 * t50 + t6 * t54) - g(2) * (t18 * t50 + t5 * t54) - g(3) * (t22 * t54 + t23 * t50) -g(1) * (t20 * t54 - t6 * t50) - g(2) * (t18 * t54 - t5 * t50) - g(3) * (-t22 * t50 + t23 * t54); 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t85 * t52 + t62 * t56, 0, t58 * pkin(3), 0, 0, 0, 0, 0, -t63 * t55, t63 * t51, 0, 0, 0, 0, 0, -g(1) * (t11 * t71 + t50 * t60) - g(2) * (t50 * t61 + t71 * t8) - g(3) * (t16 * t71 + t50 * t59) -g(1) * (-t11 * t72 + t54 * t60) - g(2) * (t54 * t61 - t8 * t72) - g(3) * (-t16 * t72 + t54 * t59); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t26 - g(2) * t25 - g(3) * t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, g(1) * t4 + g(2) * t2 + g(3) * t14, 0, 0, 0, 0, 0, -t64 * t54, t64 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t11 * t54 - t4 * t50) - g(2) * (-t2 * t50 - t8 * t54) - g(3) * (-t14 * t50 - t16 * t54) -g(1) * (t11 * t50 - t4 * t54) - g(2) * (-t2 * t54 + t8 * t50) - g(3) * (-t14 * t54 + t16 * t50);];
taug_reg  = t1;

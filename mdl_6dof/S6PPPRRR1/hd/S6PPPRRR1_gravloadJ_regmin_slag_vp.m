% Calculate minimal parameter regressor of gravitation load for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% taug_reg [6x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPPRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_gravloadJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = sin(pkin(14));
t64 = sin(pkin(12));
t68 = cos(pkin(13));
t56 = t64 * t68;
t63 = sin(pkin(13));
t69 = cos(pkin(12));
t57 = t69 * t63;
t72 = cos(pkin(6));
t47 = t72 * t57 + t56;
t67 = cos(pkin(14));
t55 = t64 * t63;
t58 = t69 * t68;
t46 = t72 * t58 - t55;
t27 = sin(pkin(6));
t66 = sin(pkin(7));
t62 = t27 * t66;
t71 = cos(pkin(7));
t77 = t46 * t71 - t69 * t62;
t35 = t47 * t26 - t77 * t67;
t40 = t69 * t27 * t71 + t46 * t66;
t65 = sin(pkin(8));
t70 = cos(pkin(8));
t80 = t35 * t70 + t40 * t65;
t49 = -t72 * t55 + t58;
t48 = -t72 * t56 - t57;
t61 = t27 * t64;
t76 = t48 * t71 + t66 * t61;
t36 = t49 * t26 - t76 * t67;
t41 = t48 * t66 - t71 * t61;
t79 = t36 * t70 + t41 * t65;
t59 = t71 * t68;
t60 = t72 * t66;
t42 = t67 * t60 + (-t63 * t26 + t67 * t59) * t27;
t50 = t68 * t62 - t72 * t71;
t78 = t42 * t70 - t50 * t65;
t75 = cos(qJ(4));
t28 = sin(qJ(6));
t32 = cos(qJ(5));
t74 = t28 * t32;
t31 = cos(qJ(6));
t73 = t31 * t32;
t18 = t76 * t26 + t49 * t67;
t30 = sin(qJ(4));
t10 = t18 * t75 - t79 * t30;
t23 = t27 * t63 * t67 + (t27 * t59 + t60) * t26;
t12 = t23 * t75 + t78 * t30;
t13 = t35 * t65 - t40 * t70;
t14 = t36 * t65 - t41 * t70;
t19 = -t42 * t65 - t50 * t70;
t29 = sin(qJ(5));
t17 = t77 * t26 + t47 * t67;
t8 = t17 * t75 - t80 * t30;
t52 = g(1) * (-t10 * t29 + t14 * t32) + g(2) * (t13 * t32 - t8 * t29) + g(3) * (-t12 * t29 + t19 * t32);
t11 = t23 * t30 - t78 * t75;
t7 = t17 * t30 + t80 * t75;
t9 = t18 * t30 + t79 * t75;
t51 = g(1) * t9 + g(2) * t7 + g(3) * t11;
t25 = -g(3) * t72 + (-t64 * g(1) + t69 * g(2)) * t27;
t6 = t12 * t32 + t19 * t29;
t4 = t10 * t32 + t14 * t29;
t2 = t13 * t29 + t8 * t32;
t1 = [-g(3), -g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, t25, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t41 + g(2) * t40 + g(3) * t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t51, g(1) * t10 + g(2) * t8 + g(3) * t12, 0, 0, 0, 0, 0, t51 * t32, -t51 * t29, 0, 0, 0, 0, 0, -g(1) * (t10 * t28 - t9 * t73) - g(2) * (t8 * t28 - t7 * t73) - g(3) * (-t11 * t73 + t12 * t28) -g(1) * (t10 * t31 + t9 * t74) - g(2) * (t8 * t31 + t7 * t74) - g(3) * (t11 * t74 + t12 * t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, g(1) * t4 + g(2) * t2 + g(3) * t6, 0, 0, 0, 0, 0, -t52 * t31, t52 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t28 + t9 * t31) - g(2) * (-t2 * t28 + t7 * t31) - g(3) * (t11 * t31 - t6 * t28) -g(1) * (-t9 * t28 - t4 * t31) - g(2) * (-t2 * t31 - t7 * t28) - g(3) * (-t11 * t28 - t6 * t31);];
taug_reg  = t1;

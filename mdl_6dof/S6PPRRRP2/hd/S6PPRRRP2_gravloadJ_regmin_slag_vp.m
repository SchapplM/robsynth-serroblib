% Calculate minimal parameter regressor of gravitation load for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [6x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t76 = sin(pkin(12));
t77 = sin(pkin(11));
t67 = t77 * t76;
t80 = cos(pkin(12));
t81 = cos(pkin(11));
t75 = t81 * t80;
t83 = cos(pkin(6));
t55 = -t83 * t75 + t67;
t78 = sin(pkin(7));
t79 = sin(pkin(6));
t72 = t79 * t78;
t82 = cos(pkin(7));
t95 = t55 * t82 + t81 * t72;
t69 = t77 * t80;
t73 = t81 * t76;
t56 = t83 * t69 + t73;
t68 = t77 * t79;
t94 = t56 * t82 - t78 * t68;
t93 = t80 * t82 * t79 + t78 * t83;
t42 = t83 * t73 + t69;
t50 = sin(qJ(3));
t86 = cos(qJ(3));
t26 = t42 * t50 + t95 * t86;
t43 = -t83 * t67 + t75;
t28 = t43 * t50 + t94 * t86;
t71 = t79 * t76;
t34 = t50 * t71 - t93 * t86;
t60 = g(1) * t28 + g(2) * t26 + g(3) * t34;
t27 = t42 * t86 - t95 * t50;
t29 = t43 * t86 - t94 * t50;
t35 = t93 * t50 + t86 * t71;
t74 = t81 * t79;
t36 = t55 * t78 - t82 * t74;
t37 = t56 * t78 + t82 * t68;
t41 = -t80 * t72 + t83 * t82;
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t62 = g(3) * (-t35 * t49 + t41 * t52) + g(2) * (-t27 * t49 + t36 * t52) + g(1) * (-t29 * t49 + t37 * t52);
t48 = sin(qJ(5));
t85 = t48 * t52;
t51 = cos(qJ(5));
t84 = t51 * t52;
t31 = t35 * t52 + t41 * t49;
t18 = t31 * t48 - t34 * t51;
t15 = t27 * t52 + t36 * t49;
t5 = t15 * t48 - t26 * t51;
t17 = t29 * t52 + t37 * t49;
t7 = t17 * t48 - t28 * t51;
t1 = g(1) * t7 + g(2) * t5 + g(3) * t18;
t19 = t31 * t51 + t34 * t48;
t6 = t15 * t51 + t26 * t48;
t8 = t17 * t51 + t28 * t48;
t64 = g(1) * t8 + g(2) * t6 + g(3) * t19;
t10 = -t26 * t85 - t27 * t51;
t12 = -t28 * t85 - t29 * t51;
t20 = -t34 * t85 - t35 * t51;
t63 = g(1) * t12 + g(2) * t10 + g(3) * t20;
t61 = g(1) * t17 + g(2) * t15 + g(3) * t31;
t40 = -g(1) * t68 + g(2) * t74 - g(3) * t83;
t21 = -t34 * t84 + t35 * t48;
t13 = -t28 * t84 + t29 * t48;
t11 = -t26 * t84 + t27 * t48;
t9 = t60 * t49;
t4 = t62 * t51;
t3 = t62 * t48;
t2 = -g(1) * t13 - g(2) * t11 - g(3) * t21;
t14 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40; 0, 0, 0, t60, g(1) * t29 + g(2) * t27 + g(3) * t35, 0, 0, 0, 0, 0, t60 * t52, -t9, 0, 0, 0, 0, 0, t2, t63, t2, t9, -t63, -g(1) * (t13 * pkin(5) + t29 * pkin(9) + t12 * qJ(6)) - g(2) * (t11 * pkin(5) + t27 * pkin(9) + t10 * qJ(6)) - g(3) * (t21 * pkin(5) + t35 * pkin(9) + t20 * qJ(6)) + t60 * (pkin(4) * t52 + pkin(10) * t49 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t61, 0, 0, 0, 0, 0, -t4, t3, -t4, -t61, -t3, -t61 * pkin(10) - t62 * (pkin(5) * t51 + qJ(6) * t48 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t64, t1, 0, -t64, -g(1) * (-t7 * pkin(5) + t8 * qJ(6)) - g(2) * (-t5 * pkin(5) + t6 * qJ(6)) - g(3) * (-t18 * pkin(5) + t19 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t14;

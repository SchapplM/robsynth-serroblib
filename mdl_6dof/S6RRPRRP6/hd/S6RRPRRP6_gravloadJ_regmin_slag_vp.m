% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t63 = sin(pkin(6));
t72 = cos(qJ(1));
t101 = t63 * t72;
t64 = cos(pkin(6));
t62 = sin(pkin(11));
t67 = sin(qJ(2));
t71 = cos(qJ(2));
t93 = cos(pkin(11));
t79 = -t71 * t62 - t67 * t93;
t46 = t79 * t64;
t68 = sin(qJ(1));
t78 = -t67 * t62 + t71 * t93;
t34 = -t72 * t46 + t68 * t78;
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t20 = -t66 * t101 + t34 * t70;
t45 = t78 * t64;
t33 = t72 * t45 + t68 * t79;
t65 = sin(qJ(5));
t69 = cos(qJ(5));
t7 = t20 * t65 + t33 * t69;
t8 = t20 * t69 - t33 * t65;
t95 = t72 * t67;
t97 = t68 * t71;
t51 = -t64 * t97 - t95;
t102 = t63 * t71;
t114 = -g(1) * t51 - g(3) * t102;
t36 = -t68 * t45 + t72 * t79;
t43 = t78 * t63;
t113 = -g(1) * t36 - g(2) * t33 - g(3) * t43;
t103 = t63 * t68;
t35 = -t68 * t46 - t72 * t78;
t23 = -t70 * t103 - t35 * t66;
t44 = t79 * t63;
t89 = -t70 * t101 - t34 * t66;
t75 = -g(3) * (t44 * t66 + t64 * t70) - g(2) * t89 + g(1) * t23;
t100 = t65 * t70;
t98 = t68 * t67;
t96 = t69 * t70;
t94 = t72 * t71;
t90 = t64 * t94;
t47 = t64 * t67 * pkin(2) + (-pkin(8) - qJ(3)) * t63;
t61 = t71 * pkin(2) + pkin(1);
t88 = -t68 * t47 + t72 * t61;
t24 = t66 * t103 - t35 * t70;
t11 = t24 * t65 + t36 * t69;
t86 = -g(1) * t7 + g(2) * t11;
t85 = g(1) * t89 + g(2) * t23;
t84 = g(1) * t72 + g(2) * t68;
t83 = g(1) * t68 - g(2) * t72;
t82 = -t72 * t47 - t68 * t61;
t39 = -t44 * t70 + t64 * t66;
t17 = t39 * t65 + t43 * t69;
t1 = g(1) * t11 + g(2) * t7 + g(3) * t17;
t12 = t24 * t69 - t36 * t65;
t18 = t39 * t69 - t43 * t65;
t77 = g(1) * t12 + g(2) * t8 + g(3) * t18;
t13 = t33 * t100 - t34 * t69;
t15 = t36 * t100 + t35 * t69;
t25 = t43 * t100 + t44 * t69;
t76 = g(1) * t15 + g(2) * t13 + g(3) * t25;
t74 = g(1) * t24 + g(2) * t20 + g(3) * t39;
t54 = pkin(2) * t90;
t52 = -t64 * t98 + t94;
t50 = -t64 * t95 - t97;
t49 = -t90 + t98;
t42 = -g(3) * t64 - t83 * t63;
t26 = t43 * t96 - t44 * t65;
t16 = -t35 * t65 + t36 * t96;
t14 = t33 * t96 + t34 * t65;
t6 = t113 * t66;
t5 = t75 * t69;
t4 = t75 * t65;
t3 = g(1) * t8 - g(2) * t12;
t2 = -g(1) * t16 - g(2) * t14 - g(3) * t26;
t9 = [0, t83, t84, 0, 0, 0, 0, 0, -g(1) * t50 - g(2) * t52, -g(1) * t49 - g(2) * t51, -t84 * t63, -g(1) * t82 - g(2) * t88, 0, 0, 0, 0, 0, g(1) * t20 - g(2) * t24, t85, 0, 0, 0, 0, 0, t3, t86, t3, -t85, -t86, -g(1) * (-t34 * pkin(3) - pkin(4) * t20 - pkin(5) * t8 + t33 * pkin(9) + pkin(10) * t89 - qJ(6) * t7 + t82) - g(2) * (-pkin(3) * t35 + t24 * pkin(4) + t12 * pkin(5) - t36 * pkin(9) + t23 * pkin(10) + t11 * qJ(6) + t88); 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t49 + t114, g(3) * t63 * t67 + g(1) * t52 - g(2) * t50, 0, -g(2) * t54 + (g(2) * t98 + t114) * pkin(2), 0, 0, 0, 0, 0, t113 * t70, -t6, 0, 0, 0, 0, 0, t2, t76, t2, t6, -t76, -g(1) * (t51 * pkin(2) + t16 * pkin(5) - t35 * pkin(9) + t15 * qJ(6)) - g(2) * (-pkin(2) * t98 + t14 * pkin(5) + pkin(9) * t34 + t13 * qJ(6) + t54) - g(3) * (pkin(2) * t102 + t26 * pkin(5) - t44 * pkin(9) + t25 * qJ(6)) + t113 * (pkin(4) * t70 + pkin(10) * t66 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, 0, 0, 0, 0, 0, t5, -t4, t5, -t74, t4, -pkin(10) * t74 + t75 * (pkin(5) * t69 + qJ(6) * t65 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t77, t1, 0, -t77, -g(1) * (-t11 * pkin(5) + t12 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - g(3) * (-t17 * pkin(5) + t18 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t9;

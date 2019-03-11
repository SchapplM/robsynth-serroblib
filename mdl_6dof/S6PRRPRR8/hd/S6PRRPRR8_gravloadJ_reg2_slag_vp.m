% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t65 = sin(qJ(2));
t68 = cos(qJ(2));
t100 = cos(pkin(6));
t61 = cos(pkin(12));
t92 = t61 * t100;
t98 = sin(pkin(12));
t47 = t65 * t92 + t98 * t68;
t82 = t100 * t98;
t49 = t61 * t68 - t65 * t82;
t120 = -g(1) * t49 - g(2) * t47;
t60 = sin(pkin(6));
t106 = t60 * t61;
t109 = cos(qJ(3));
t46 = -t98 * t65 + t68 * t92;
t59 = sin(pkin(7));
t64 = sin(qJ(3));
t99 = cos(pkin(7));
t91 = t64 * t99;
t16 = -t59 * t64 * t106 + t47 * t109 + t46 * t91;
t48 = -t61 * t65 - t68 * t82;
t93 = t60 * t98;
t86 = t59 * t93;
t18 = t49 * t109 + (t99 * t48 + t86) * t64;
t90 = t100 * t59;
t32 = t64 * t90 + (t109 * t65 + t68 * t91) * t60;
t74 = g(1) * t18 + g(2) * t16 + g(3) * t32;
t118 = pkin(9) * t59;
t83 = t99 * t109;
t94 = t60 * t109;
t15 = t61 * t59 * t94 - t46 * t83 + t47 * t64;
t112 = t15 * pkin(10);
t17 = -t109 * t86 - t48 * t83 + t49 * t64;
t111 = t17 * pkin(10);
t104 = t60 * t68;
t105 = t60 * t65;
t31 = -t83 * t104 + t64 * t105 - t109 * t90;
t110 = t31 * pkin(10);
t63 = sin(qJ(5));
t108 = t59 * t63;
t67 = cos(qJ(5));
t107 = t59 * t67;
t62 = sin(qJ(6));
t103 = t62 * t63;
t66 = cos(qJ(6));
t102 = t63 * t66;
t95 = t59 * t105;
t101 = pkin(2) * t104 + pkin(9) * t95;
t97 = t47 * t118;
t96 = t49 * t118;
t13 = t15 * pkin(3);
t89 = t16 * qJ(4) - t13;
t14 = t17 * pkin(3);
t88 = t18 * qJ(4) - t14;
t30 = t31 * pkin(3);
t87 = t32 * qJ(4) - t30;
t24 = t46 * t64 + t47 * t83;
t25 = t46 * t109 - t47 * t91;
t43 = t46 * pkin(2);
t85 = t25 * pkin(3) + t24 * qJ(4) + t43;
t26 = t48 * t64 + t49 * t83;
t27 = t48 * t109 - t49 * t91;
t44 = t48 * pkin(2);
t84 = t27 * pkin(3) + t26 * qJ(4) + t44;
t41 = (t64 * t68 + t65 * t83) * t60;
t42 = -t91 * t105 + t68 * t94;
t81 = t42 * pkin(3) + t41 * qJ(4) + t101;
t45 = t100 * t99 - t59 * t104;
t19 = t31 * t67 - t45 * t63;
t33 = -t99 * t106 - t46 * t59;
t5 = t15 * t67 - t33 * t63;
t34 = -t48 * t59 + t99 * t93;
t7 = t17 * t67 - t34 * t63;
t79 = g(1) * t7 + g(2) * t5 + g(3) * t19;
t20 = t31 * t63 + t45 * t67;
t6 = t15 * t63 + t33 * t67;
t8 = t17 * t63 + t34 * t67;
t78 = g(1) * t8 + g(2) * t6 + g(3) * t20;
t11 = -t49 * t108 + t26 * t67;
t28 = -t41 * t67 + t63 * t95;
t9 = -t47 * t108 + t24 * t67;
t77 = g(1) * t11 + g(2) * t9 - g(3) * t28;
t76 = t25 * pkin(10) + t85;
t75 = t27 * pkin(10) + t84;
t3 = g(1) * t17 + g(2) * t15 + g(3) * t31;
t73 = g(1) * t26 + g(2) * t24 + g(3) * t41;
t72 = g(1) * t27 + g(2) * t25 + g(3) * t42;
t71 = g(3) * t105 - t120;
t70 = pkin(4) * t95 + t42 * pkin(10) + t81;
t69 = t120 * (pkin(4) + pkin(9)) * t59;
t29 = t41 * t63 + t67 * t95;
t21 = t71 * t59;
t12 = t49 * t107 + t26 * t63;
t10 = t47 * t107 + t24 * t63;
t1 = t74 * t67;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t48 - g(2) * t46 - g(3) * t104, t71, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t73, -t21, -g(1) * (t44 + t96) - g(2) * (t43 + t97) - g(3) * t101, 0, 0, 0, 0, 0, 0, -t21, t72, -t73, -g(1) * (t84 + t96) - g(2) * (t85 + t97) - g(3) * t81, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10 - g(3) * t29, -t77, -t72, -g(1) * t75 - g(2) * t76 - g(3) * t70 + t69, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t66 + t27 * t62) - g(2) * (t10 * t66 + t25 * t62) - g(3) * (t29 * t66 + t42 * t62) -g(1) * (-t12 * t62 + t27 * t66) - g(2) * (-t10 * t62 + t25 * t66) - g(3) * (-t29 * t62 + t42 * t66) t77, -g(1) * (t12 * pkin(5) - t11 * pkin(11) + t75) - g(2) * (t10 * pkin(5) - t9 * pkin(11) + t76) - g(3) * (t29 * pkin(5) + t28 * pkin(11) + t70) + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t74, -g(1) * t88 - g(2) * t89 - g(3) * t87, 0, 0, 0, 0, 0, 0, -t74 * t63, -t1, t3, -g(1) * (t88 - t111) - g(2) * (t89 - t112) - g(3) * (t87 - t110) 0, 0, 0, 0, 0, 0, -g(1) * (t18 * t102 - t17 * t62) - g(2) * (t16 * t102 - t15 * t62) - g(3) * (t32 * t102 - t31 * t62) -g(1) * (-t18 * t103 - t17 * t66) - g(2) * (-t16 * t103 - t15 * t66) - g(3) * (-t32 * t103 - t31 * t66) t1, -g(1) * (-t14 - t111) - g(2) * (-t13 - t112) - g(3) * (-t30 - t110) - t74 * (pkin(5) * t63 - pkin(11) * t67 + qJ(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t78, 0, 0, 0, 0, 0, 0, 0, 0, -t79 * t66, t79 * t62, -t78, -g(1) * (t7 * pkin(5) + t8 * pkin(11)) - g(2) * (t5 * pkin(5) + t6 * pkin(11)) - g(3) * (t19 * pkin(5) + t20 * pkin(11)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t18 * t66 - t8 * t62) - g(2) * (t16 * t66 - t6 * t62) - g(3) * (-t20 * t62 + t32 * t66) -g(1) * (-t18 * t62 - t8 * t66) - g(2) * (-t16 * t62 - t6 * t66) - g(3) * (-t20 * t66 - t32 * t62) 0, 0;];
taug_reg  = t2;

% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t105 = sin(pkin(12));
t74 = sin(pkin(6));
t102 = t74 * t105;
t106 = cos(pkin(7));
t73 = sin(pkin(7));
t75 = cos(pkin(12));
t79 = sin(qJ(2));
t82 = cos(qJ(2));
t107 = cos(pkin(6));
t96 = t107 * t105;
t86 = t75 * t79 + t82 * t96;
t124 = -t102 * t73 + t106 * t86;
t112 = t74 * t75;
t101 = t75 * t107;
t85 = -t101 * t82 + t105 * t79;
t123 = t106 * t85 + t112 * t73;
t115 = cos(qJ(3));
t65 = t101 * t79 + t105 * t82;
t78 = sin(qJ(3));
t37 = t115 * t123 + t65 * t78;
t66 = t75 * t82 - t79 * t96;
t39 = t115 * t124 + t66 * t78;
t110 = t74 * t82;
t111 = t74 * t79;
t97 = t106 * t115;
t99 = t107 * t73;
t53 = -t110 * t97 + t111 * t78 - t115 * t99;
t87 = g(1) * t39 + g(2) * t37 + g(3) * t53;
t38 = t65 * t115 - t123 * t78;
t40 = t115 * t66 - t124 * t78;
t100 = t78 * t106;
t54 = t78 * t99 + (t100 * t82 + t115 * t79) * t74;
t55 = -t106 * t112 + t73 * t85;
t56 = t102 * t106 + t73 * t86;
t64 = t106 * t107 - t110 * t73;
t77 = sin(qJ(4));
t81 = cos(qJ(4));
t90 = g(3) * (-t54 * t77 + t64 * t81) + g(2) * (-t38 * t77 + t55 * t81) + g(1) * (-t40 * t77 + t56 * t81);
t122 = pkin(9) * t73;
t114 = t73 * t77;
t113 = t73 * t81;
t76 = sin(qJ(5));
t109 = t76 * t81;
t80 = cos(qJ(5));
t108 = t80 * t81;
t103 = t73 * t111;
t42 = t54 * t81 + t64 * t77;
t19 = t42 * t76 - t53 * t80;
t22 = t38 * t81 + t55 * t77;
t6 = t22 * t76 - t37 * t80;
t24 = t40 * t81 + t56 * t77;
t8 = t24 * t76 - t39 * t80;
t1 = g(1) * t8 + g(2) * t6 + g(3) * t19;
t20 = t42 * t80 + t53 * t76;
t7 = t22 * t80 + t37 * t76;
t9 = t24 * t80 + t39 * t76;
t93 = g(1) * t9 + g(2) * t7 + g(3) * t20;
t11 = -t109 * t37 - t38 * t80;
t13 = -t109 * t39 - t40 * t80;
t25 = -t109 * t53 - t54 * t80;
t92 = g(1) * t13 + g(2) * t11 + g(3) * t25;
t46 = -t100 * t65 - t115 * t85;
t28 = t114 * t65 + t46 * t81;
t45 = t65 * t97 - t78 * t85;
t15 = t28 * t76 - t45 * t80;
t48 = -t100 * t66 - t115 * t86;
t30 = t114 * t66 + t48 * t81;
t47 = t66 * t97 - t78 * t86;
t17 = t30 * t76 - t47 * t80;
t63 = (-t100 * t79 + t115 * t82) * t74;
t50 = t103 * t77 + t63 * t81;
t62 = (t78 * t82 + t79 * t97) * t74;
t31 = t50 * t76 - t62 * t80;
t91 = g(1) * t17 + g(2) * t15 + g(3) * t31;
t89 = g(1) * t24 + g(2) * t22 + g(3) * t42;
t27 = -t113 * t65 + t46 * t77;
t29 = -t113 * t66 + t48 * t77;
t49 = -t103 * t81 + t63 * t77;
t88 = g(1) * t29 + g(2) * t27 + g(3) * t49;
t32 = t50 * t80 + t62 * t76;
t26 = -t108 * t53 + t54 * t76;
t18 = t30 * t80 + t47 * t76;
t16 = t28 * t80 + t45 * t76;
t14 = -t108 * t39 + t40 * t76;
t12 = -t108 * t37 + t38 * t76;
t10 = t87 * t77;
t5 = t90 * t80;
t4 = t90 * t76;
t3 = -g(1) * t18 - g(2) * t16 - g(3) * t32;
t2 = -g(1) * t14 - g(2) * t12 - g(3) * t26;
t21 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, g(1) * t86 + g(2) * t85 - g(3) * t110, g(1) * t66 + g(2) * t65 + g(3) * t111, 0, 0, 0, 0, 0, -g(1) * t48 - g(2) * t46 - g(3) * t63, g(1) * t47 + g(2) * t45 + g(3) * t62, 0, 0, 0, 0, 0, -g(1) * t30 - g(2) * t28 - g(3) * t50, t88, 0, 0, 0, 0, 0, t3, t91, t3, -t88, -t91, -g(1) * (-pkin(2) * t86 + t48 * pkin(3) + t30 * pkin(4) + t18 * pkin(5) + t47 * pkin(10) + t29 * pkin(11) + t17 * qJ(6) + t122 * t66) - g(2) * (-pkin(2) * t85 + t46 * pkin(3) + t28 * pkin(4) + t16 * pkin(5) + t45 * pkin(10) + t27 * pkin(11) + t15 * qJ(6) + t122 * t65) - g(3) * (t63 * pkin(3) + t50 * pkin(4) + t32 * pkin(5) + t62 * pkin(10) + t49 * pkin(11) + t31 * qJ(6) + (pkin(2) * t82 + t122 * t79) * t74); 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, g(1) * t40 + g(2) * t38 + g(3) * t54, 0, 0, 0, 0, 0, t87 * t81, -t10, 0, 0, 0, 0, 0, t2, t92, t2, t10, -t92, -g(1) * (t14 * pkin(5) + t40 * pkin(10) + t13 * qJ(6)) - g(2) * (t12 * pkin(5) + t38 * pkin(10) + t11 * qJ(6)) - g(3) * (t26 * pkin(5) + t54 * pkin(10) + t25 * qJ(6)) + t87 * (pkin(4) * t81 + pkin(11) * t77 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, t89, 0, 0, 0, 0, 0, -t5, t4, -t5, -t89, -t4, -t89 * pkin(11) - t90 * (pkin(5) * t80 + qJ(6) * t76 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t93, t1, 0, -t93, -g(1) * (-pkin(5) * t8 + qJ(6) * t9) - g(2) * (-pkin(5) * t6 + qJ(6) * t7) - g(3) * (-pkin(5) * t19 + qJ(6) * t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t21;

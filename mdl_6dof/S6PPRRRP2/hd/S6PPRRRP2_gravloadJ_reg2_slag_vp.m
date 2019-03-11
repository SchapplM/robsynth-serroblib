% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t62 = sin(qJ(4));
t65 = cos(qJ(4));
t114 = -pkin(4) * t65 - pkin(10) * t62;
t103 = cos(pkin(11));
t104 = cos(pkin(7));
t105 = cos(pkin(6));
t98 = sin(pkin(12));
t99 = sin(pkin(11));
t79 = t99 * t98;
t102 = cos(pkin(12));
t86 = t103 * t102;
t68 = -t105 * t86 + t79;
t100 = sin(pkin(7));
t101 = sin(pkin(6));
t83 = t101 * t100;
t113 = t103 * t83 + t68 * t104;
t81 = t99 * t102;
t84 = t103 * t98;
t69 = t105 * t81 + t84;
t80 = t99 * t101;
t112 = -t100 * t80 + t69 * t104;
t111 = t102 * t104 * t101 + t105 * t100;
t108 = cos(qJ(3));
t61 = sin(qJ(5));
t107 = t61 * t65;
t64 = cos(qJ(5));
t106 = t64 * t65;
t55 = t105 * t84 + t81;
t63 = sin(qJ(3));
t36 = t113 * t108 + t55 * t63;
t37 = t55 * t108 - t113 * t63;
t97 = -t36 * pkin(3) + t37 * pkin(9);
t56 = -t105 * t79 + t86;
t38 = t112 * t108 + t56 * t63;
t39 = t56 * t108 - t112 * t63;
t96 = -t38 * pkin(3) + t39 * pkin(9);
t82 = t101 * t98;
t47 = -t111 * t108 + t63 * t82;
t48 = t108 * t82 + t111 * t63;
t95 = -t47 * pkin(3) + t48 * pkin(9);
t85 = t103 * t101;
t49 = t68 * t100 - t104 * t85;
t17 = -t37 * t62 + t49 * t65;
t18 = t37 * t65 + t49 * t62;
t94 = t17 * pkin(4) + t18 * pkin(10);
t50 = t69 * t100 + t104 * t80;
t19 = -t39 * t62 + t50 * t65;
t20 = t39 * t65 + t50 * t62;
t93 = t19 * pkin(4) + t20 * pkin(10);
t54 = -t102 * t83 + t105 * t104;
t40 = -t48 * t62 + t54 * t65;
t41 = t48 * t65 + t54 * t62;
t92 = t40 * pkin(4) + t41 * pkin(10);
t91 = pkin(5) * t64 + qJ(6) * t61;
t90 = t114 * t36 + t97;
t89 = t114 * t38 + t96;
t88 = t114 * t47 + t95;
t21 = t41 * t61 - t47 * t64;
t6 = t18 * t61 - t36 * t64;
t8 = t20 * t61 - t38 * t64;
t1 = g(1) * t8 + g(2) * t6 + g(3) * t21;
t22 = t41 * t64 + t47 * t61;
t7 = t18 * t64 + t36 * t61;
t9 = t20 * t64 + t38 * t61;
t78 = g(1) * t9 + g(2) * t7 + g(3) * t22;
t11 = -t36 * t107 - t37 * t64;
t13 = -t38 * t107 - t39 * t64;
t23 = -t47 * t107 - t48 * t64;
t77 = g(1) * t13 + g(2) * t11 + g(3) * t23;
t76 = g(1) * t19 + g(2) * t17 + g(3) * t40;
t75 = g(1) * t20 + g(2) * t18 + g(3) * t41;
t74 = g(1) * t38 + g(2) * t36 + g(3) * t47;
t73 = g(1) * t39 + g(2) * t37 + g(3) * t48;
t53 = -g(1) * t80 + g(2) * t85 - g(3) * t105;
t24 = -t47 * t106 + t48 * t61;
t14 = -t38 * t106 + t39 * t61;
t12 = -t36 * t106 + t37 * t61;
t10 = t74 * t62;
t4 = t76 * t64;
t3 = t76 * t61;
t2 = -g(1) * t14 - g(2) * t12 - g(3) * t24;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t73, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t65, -t10, -t73, -g(1) * t96 - g(2) * t97 - g(3) * t95, 0, 0, 0, 0, 0, 0, t2, t77, t10, -g(1) * t89 - g(2) * t90 - g(3) * t88, 0, 0, 0, 0, 0, 0, t2, t10, -t77, -g(1) * (t14 * pkin(5) + t13 * qJ(6) + t89) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t90) - g(3) * (t24 * pkin(5) + t23 * qJ(6) + t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t75, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t75, -g(1) * t93 - g(2) * t94 - g(3) * t92, 0, 0, 0, 0, 0, 0, -t4, -t75, -t3, -g(1) * (t91 * t19 + t93) - g(2) * (t91 * t17 + t94) - g(3) * (t91 * t40 + t92); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t78, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t78, -g(1) * (-t8 * pkin(5) + t9 * qJ(6)) - g(2) * (-t6 * pkin(5) + t7 * qJ(6)) - g(3) * (-t21 * pkin(5) + t22 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;

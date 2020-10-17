% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:10:49
% EndTime: 2019-05-05 05:10:51
% DurationCPUTime: 0.73s
% Computational Cost: add. (566->129), mult. (1491->182), div. (0->0), fcn. (1873->12), ass. (0->80)
t63 = sin(qJ(6));
t61 = sin(pkin(6));
t66 = sin(qJ(2));
t104 = t61 * t66;
t65 = sin(qJ(3));
t69 = cos(qJ(3));
t97 = cos(pkin(6));
t49 = t65 * t104 - t97 * t69;
t103 = t61 * t69;
t50 = t66 * t103 + t97 * t65;
t64 = sin(qJ(5));
t68 = cos(qJ(5));
t21 = t49 * t68 - t50 * t64;
t70 = cos(qJ(2));
t62 = cos(pkin(11));
t90 = t62 * t97;
t96 = sin(pkin(11));
t46 = t66 * t90 + t96 * t70;
t27 = t62 * t103 + t46 * t65;
t28 = -t62 * t61 * t65 + t46 * t69;
t5 = t27 * t68 - t28 * t64;
t81 = t97 * t96;
t48 = t62 * t70 - t66 * t81;
t91 = t61 * t96;
t29 = t48 * t65 - t69 * t91;
t30 = t48 * t69 + t65 * t91;
t9 = t29 * t68 - t30 * t64;
t78 = g(1) * t9 + g(2) * t5 + g(3) * t21;
t115 = t78 * t63;
t67 = cos(qJ(6));
t114 = t78 * t67;
t6 = t27 * t64 + t28 * t68;
t113 = t5 * pkin(5) + t6 * pkin(10);
t10 = t29 * t64 + t30 * t68;
t112 = t9 * pkin(5) + t10 * pkin(10);
t22 = t49 * t64 + t50 * t68;
t111 = t21 * pkin(5) + t22 * pkin(10);
t77 = g(1) * t10 + g(2) * t6 + g(3) * t22;
t45 = -t96 * t66 + t70 * t90;
t106 = t45 * t69;
t98 = qJ(4) * t65;
t110 = pkin(3) * t106 + t45 * t98;
t47 = -t62 * t66 - t70 * t81;
t105 = t47 * t69;
t109 = pkin(3) * t105 + t47 * t98;
t108 = t64 * t69 - t65 * t68;
t107 = pkin(8) - pkin(9);
t102 = t61 * t70;
t99 = pkin(2) * t102 + pkin(8) * t104;
t95 = t65 * t102;
t94 = t69 * t102;
t41 = t45 * pkin(2);
t93 = t46 * pkin(8) + t41;
t42 = t47 * pkin(2);
t92 = t48 * pkin(8) + t42;
t89 = -t27 * pkin(3) + t28 * qJ(4);
t88 = -t29 * pkin(3) + t30 * qJ(4);
t87 = -t49 * pkin(3) + t50 * qJ(4);
t86 = pkin(3) * t94 + qJ(4) * t95 + t99;
t85 = -t27 * pkin(4) + t89;
t84 = -t29 * pkin(4) + t88;
t83 = -t49 * pkin(4) + t87;
t82 = t64 * t65 + t68 * t69;
t13 = t108 * t45;
t15 = t108 * t47;
t31 = t64 * t94 - t68 * t95;
t76 = g(1) * t15 + g(2) * t13 + g(3) * t31;
t2 = g(1) * t29 + g(2) * t27 + g(3) * t49;
t75 = g(1) * t30 + g(2) * t28 + g(3) * t50;
t74 = pkin(4) * t106 + t107 * t46 + t110 + t41;
t73 = pkin(4) * t105 + t107 * t48 + t109 + t42;
t72 = g(1) * t47 + g(2) * t45 + g(3) * t102;
t18 = g(1) * t48 + g(2) * t46 + g(3) * t104;
t71 = pkin(4) * t94 - pkin(9) * t104 + t86;
t32 = t82 * t102;
t16 = t82 * t47;
t14 = t82 * t45;
t12 = t72 * t69;
t11 = t72 * t65;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t18, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, -t18, -g(1) * t92 - g(2) * t93 - g(3) * t99, 0, 0, 0, 0, 0, 0, -t12, -t18, -t11, -g(1) * (t92 + t109) - g(2) * (t93 + t110) - g(3) * t86, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14 - g(3) * t32, t76, t18, -g(1) * t73 - g(2) * t74 - g(3) * t71, 0, 0, 0, 0, 0, 0, -g(1) * (t16 * t67 - t48 * t63) - g(2) * (t14 * t67 - t46 * t63) - g(3) * (-t63 * t104 + t32 * t67) -g(1) * (-t16 * t63 - t48 * t67) - g(2) * (-t14 * t63 - t46 * t67) - g(3) * (-t67 * t104 - t32 * t63) -t76, -g(1) * (t16 * pkin(5) + t15 * pkin(10) + t73) - g(2) * (t14 * pkin(5) + t13 * pkin(10) + t74) - g(3) * (t32 * pkin(5) + t31 * pkin(10) + t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t75, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t75, -g(1) * t88 - g(2) * t89 - g(3) * t87, 0, 0, 0, 0, 0, 0, t78, -t77, 0, -g(1) * t84 - g(2) * t85 - g(3) * t83, 0, 0, 0, 0, 0, 0, t114, -t115, t77, -g(1) * (-t112 + t84) - g(2) * (-t113 + t85) - g(3) * (-t111 + t83); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t77, 0, 0, 0, 0, 0, 0, 0, 0, -t114, t115, -t77, -g(1) * t112 - g(2) * t113 - g(3) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t10 * t63 + t47 * t67) - g(2) * (t45 * t67 - t6 * t63) - g(3) * (t67 * t102 - t22 * t63) -g(1) * (-t10 * t67 - t47 * t63) - g(2) * (-t45 * t63 - t6 * t67) - g(3) * (-t63 * t102 - t22 * t67) 0, 0;];
taug_reg  = t1;

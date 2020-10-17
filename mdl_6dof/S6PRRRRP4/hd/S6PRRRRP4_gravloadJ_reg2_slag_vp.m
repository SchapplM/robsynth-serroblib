% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:55:42
% EndTime: 2019-05-05 09:55:45
% DurationCPUTime: 0.70s
% Computational Cost: add. (689->145), mult. (1428->212), div. (0->0), fcn. (1761->12), ass. (0->89)
t68 = sin(qJ(3));
t73 = -pkin(10) - pkin(9);
t103 = t68 * t73;
t70 = cos(qJ(4));
t61 = t70 * pkin(4) + pkin(3);
t71 = cos(qJ(3));
t119 = -t61 * t71 + t103;
t66 = sin(pkin(6));
t118 = g(3) * t66;
t65 = sin(pkin(11));
t69 = sin(qJ(2));
t72 = cos(qJ(2));
t95 = cos(pkin(11));
t96 = cos(pkin(6));
t83 = t96 * t95;
t47 = t65 * t72 + t69 * t83;
t87 = t66 * t95;
t27 = t47 * t71 - t68 * t87;
t67 = sin(qJ(4));
t117 = t27 * t67;
t88 = t65 * t96;
t49 = -t69 * t88 + t95 * t72;
t29 = t65 * t66 * t68 + t49 * t71;
t116 = t29 * t67;
t46 = t65 * t69 - t72 * t83;
t115 = t46 * t70;
t114 = t47 * t67;
t48 = t95 * t69 + t72 * t88;
t113 = t48 * t70;
t112 = t49 * t67;
t64 = qJ(4) + qJ(5);
t62 = sin(t64);
t110 = t62 * t71;
t63 = cos(t64);
t109 = t63 * t71;
t108 = t66 * t69;
t107 = t66 * t71;
t106 = t66 * t72;
t105 = t67 * t69;
t104 = t67 * t71;
t102 = t70 * t71;
t101 = t71 * t72;
t26 = -t47 * t68 - t71 * t87;
t100 = t26 * t61 - t27 * t73;
t28 = t65 * t107 - t49 * t68;
t99 = t28 * t61 - t29 * t73;
t50 = -t68 * t108 + t96 * t71;
t51 = t69 * t107 + t96 * t68;
t98 = t50 * t61 - t51 * t73;
t97 = pkin(2) * t106 + pkin(8) * t108;
t94 = t62 * t106;
t93 = t70 * t106;
t92 = -t46 * pkin(2) + t47 * pkin(8);
t91 = -t48 * pkin(2) + t49 * pkin(8);
t10 = t27 * t62 - t46 * t63;
t11 = t27 * t63 + t46 * t62;
t90 = -t10 * pkin(5) + t11 * qJ(6);
t12 = t29 * t62 - t48 * t63;
t13 = t29 * t63 + t48 * t62;
t89 = -t12 * pkin(5) + t13 * qJ(6);
t22 = t63 * t106 + t51 * t62;
t23 = t51 * t63 - t94;
t86 = -t22 * pkin(5) + t23 * qJ(6);
t85 = pkin(3) * t71 + pkin(9) * t68;
t84 = pkin(5) * t63 + qJ(6) * t62;
t82 = -t51 * t67 - t93;
t1 = g(1) * t12 + g(2) * t10 + g(3) * t22;
t3 = g(1) * t13 + g(2) * t11 + g(3) * t23;
t15 = -t46 * t110 - t47 * t63;
t17 = -t48 * t110 - t49 * t63;
t32 = -t63 * t108 + t71 * t94;
t81 = g(1) * t17 + g(2) * t15 + g(3) * t32;
t80 = g(1) * t28 + g(2) * t26 + g(3) * t50;
t79 = g(1) * t29 + g(2) * t27 + g(3) * t51;
t78 = -g(1) * t48 - g(2) * t46 + g(3) * t106;
t77 = g(1) * t49 + g(2) * t47 + g(3) * t108;
t76 = -t103 * t106 + t97 + (pkin(4) * t105 + t101 * t61) * t66;
t75 = pkin(4) * t114 + t119 * t46 + t92;
t74 = pkin(4) * t112 + t119 * t48 + t91;
t41 = pkin(4) * t113;
t39 = pkin(4) * t115;
t33 = (t63 * t101 + t62 * t69) * t66;
t18 = -t48 * t109 + t49 * t62;
t16 = -t46 * t109 + t47 * t62;
t14 = t78 * t68;
t6 = t80 * t63;
t5 = t80 * t62;
t4 = -g(1) * t18 - g(2) * t16 - g(3) * t33;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t77, 0, 0, 0, 0, 0, 0, 0, 0, -t78 * t71, t14, -t77, -g(1) * t91 - g(2) * t92 - g(3) * t97, 0, 0, 0, 0, 0, 0, -g(1) * (-t48 * t102 + t112) - g(2) * (-t46 * t102 + t114) - (t70 * t101 + t105) * t118, -g(1) * (t48 * t104 + t49 * t70) - g(2) * (t46 * t104 + t47 * t70) - (-t67 * t101 + t69 * t70) * t118, -t14, -g(1) * (-t85 * t48 + t91) - g(2) * (-t85 * t46 + t92) - g(3) * (t85 * t106 + t97) 0, 0, 0, 0, 0, 0, t4, t81, -t14, -g(1) * t74 - g(2) * t75 - g(3) * t76, 0, 0, 0, 0, 0, 0, t4, -t14, -t81, -g(1) * (t18 * pkin(5) + t17 * qJ(6) + t74) - g(2) * (t16 * pkin(5) + t15 * qJ(6) + t75) - g(3) * (t33 * pkin(5) + t32 * qJ(6) + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t79, 0, 0, 0, 0, 0, 0, 0, 0, -t80 * t70, t80 * t67, -t79, -g(1) * (t28 * pkin(3) + t29 * pkin(9)) - g(2) * (t26 * pkin(3) + t27 * pkin(9)) - g(3) * (t50 * pkin(3) + t51 * pkin(9)) 0, 0, 0, 0, 0, 0, -t6, t5, -t79, -g(1) * t99 - g(2) * t100 - g(3) * t98, 0, 0, 0, 0, 0, 0, -t6, -t79, -t5, -g(1) * (t84 * t28 + t99) - g(2) * (t84 * t26 + t100) - g(3) * (t84 * t50 + t98); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t113 - t116) - g(2) * (t115 - t117) - g(3) * t82, -g(1) * (-t29 * t70 - t48 * t67) - g(2) * (-t27 * t70 - t46 * t67) - g(3) * (t67 * t106 - t51 * t70) 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t41 - g(2) * t39 + (g(3) * t93 + t79 * t67) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (-pkin(4) * t116 + t41 + t89) - g(2) * (-pkin(4) * t117 + t39 + t90) - g(3) * (t82 * pkin(4) + t86); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t89 - g(2) * t90 - g(3) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;

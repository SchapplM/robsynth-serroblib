% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:22:08
% EndTime: 2019-05-07 08:22:11
% DurationCPUTime: 0.69s
% Computational Cost: add. (595->124), mult. (1113->196), div. (0->0), fcn. (1375->12), ass. (0->80)
t67 = sin(pkin(6));
t76 = cos(qJ(1));
t103 = t67 * t76;
t71 = sin(qJ(2));
t72 = sin(qJ(1));
t75 = cos(qJ(2));
t98 = cos(pkin(6));
t90 = t76 * t98;
t45 = t71 * t90 + t72 * t75;
t66 = qJ(3) + pkin(11);
t63 = sin(t66);
t64 = cos(t66);
t20 = -t63 * t103 + t45 * t64;
t44 = t72 * t71 - t75 * t90;
t69 = sin(qJ(5));
t73 = cos(qJ(5));
t6 = t20 * t69 - t44 * t73;
t7 = t20 * t73 + t44 * t69;
t70 = sin(qJ(3));
t74 = cos(qJ(3));
t93 = t74 * t103;
t82 = t45 * t70 + t93;
t107 = t67 * t71;
t106 = t67 * t72;
t91 = t72 * t98;
t47 = -t71 * t91 + t76 * t75;
t23 = -t64 * t106 + t47 * t63;
t83 = -t64 * t103 - t45 * t63;
t78 = -g(3) * (-t63 * t107 + t98 * t64) - g(2) * t83 + g(1) * t23;
t110 = t47 * t70;
t109 = t64 * t69;
t108 = t64 * t73;
t105 = t67 * t74;
t104 = t67 * t75;
t68 = -qJ(4) - pkin(9);
t102 = t68 * t71;
t101 = t73 * t75;
t62 = t74 * pkin(3) + pkin(2);
t100 = -t44 * t62 - t45 * t68;
t46 = t76 * t71 + t75 * t91;
t99 = -t46 * t62 - t47 * t68;
t97 = t70 * t107;
t96 = t69 * t104;
t95 = t70 * t106;
t94 = t72 * t105;
t57 = t70 * t103;
t92 = t45 * t74 - t57;
t89 = t98 * t74;
t24 = t63 * t106 + t47 * t64;
t10 = t24 * t69 - t46 * t73;
t88 = -g(1) * t6 + g(2) * t10;
t87 = t76 * pkin(1) + pkin(3) * t95 + pkin(8) * t106 - t46 * t68 + t47 * t62;
t86 = pkin(4) * t64 + pkin(10) * t63;
t85 = g(1) * t44 - g(2) * t46;
t81 = -t72 * pkin(1) + pkin(3) * t57 + pkin(8) * t103 + t44 * t68 - t45 * t62;
t34 = t64 * t107 + t98 * t63;
t17 = t67 * t101 + t34 * t69;
t1 = g(1) * t10 + g(2) * t6 + g(3) * t17;
t11 = t24 * t73 + t46 * t69;
t18 = t34 * t73 - t96;
t80 = g(1) * t11 + g(2) * t7 + g(3) * t18;
t12 = -t44 * t109 - t45 * t73;
t14 = -t46 * t109 - t47 * t73;
t27 = -t73 * t107 + t64 * t96;
t79 = g(1) * t14 + g(2) * t12 + g(3) * t27;
t16 = -g(1) * t46 - g(2) * t44 + g(3) * t104;
t77 = g(1) * t47 + g(2) * t45 + g(3) * t107;
t61 = pkin(3) * t89;
t52 = pkin(3) * t94;
t48 = t62 * t104;
t28 = (t64 * t101 + t69 * t71) * t67;
t26 = t47 * t74 + t95;
t25 = t94 - t110;
t15 = -t46 * t108 + t47 * t69;
t13 = -t44 * t108 + t45 * t69;
t5 = t78 * t73;
t4 = t78 * t69;
t3 = -g(1) * t15 - g(2) * t13 - g(3) * t28;
t2 = g(1) * t7 - g(2) * t11;
t8 = [0, g(1) * t72 - g(2) * t76, g(1) * t76 + g(2) * t72, 0, 0, 0, 0, 0, g(1) * t45 - g(2) * t47, -t85, 0, 0, 0, 0, 0, g(1) * t92 - g(2) * t26, -g(1) * t82 - g(2) * t25, t85, -g(1) * t81 - g(2) * t87, 0, 0, 0, 0, 0, t2, t88, t2, -g(1) * t83 - g(2) * t23, -t88, -g(1) * (-pkin(4) * t20 - pkin(5) * t7 + pkin(10) * t83 - qJ(6) * t6 + t81) - g(2) * (t24 * pkin(4) + t11 * pkin(5) + t23 * pkin(10) + t10 * qJ(6) + t87); 0, 0, 0, 0, 0, 0, 0, 0, -t16, t77, 0, 0, 0, 0, 0, -t16 * t74, t16 * t70, -t77, -g(1) * t99 - g(2) * t100 - g(3) * (-t67 * t102 + t48) 0, 0, 0, 0, 0, t3, t79, t3, -t16 * t63, -t79, -g(1) * (t15 * pkin(5) + t14 * qJ(6) - t86 * t46 + t99) - g(2) * (t13 * pkin(5) + t12 * qJ(6) - t86 * t44 + t100) + (-t28 * pkin(5) - t27 * qJ(6) - t48 - (t86 * t75 - t102) * t67) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t25 + g(2) * t82 - g(3) * (t89 - t97) g(1) * t26 + g(2) * t92 - g(3) * (-t71 * t105 - t98 * t70) 0, -g(1) * t52 - g(3) * t61 + (g(2) * t93 + t77 * t70) * pkin(3), 0, 0, 0, 0, 0, t5, -t4, t5, -g(1) * t24 - g(2) * t20 - g(3) * t34, t4, -g(1) * (-pkin(3) * t110 + t24 * pkin(10) + t52) - g(2) * (-t82 * pkin(3) + t20 * pkin(10)) - g(3) * (-pkin(3) * t97 + t34 * pkin(10) + t61) + t78 * (pkin(5) * t73 + qJ(6) * t69 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t80, t1, 0, -t80, -g(1) * (-t10 * pkin(5) + t11 * qJ(6)) - g(2) * (-t6 * pkin(5) + t7 * qJ(6)) - g(3) * (-t17 * pkin(5) + t18 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t8;

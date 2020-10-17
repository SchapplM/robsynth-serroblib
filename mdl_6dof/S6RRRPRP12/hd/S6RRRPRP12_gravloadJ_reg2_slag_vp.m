% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:33:46
% EndTime: 2019-05-07 09:33:50
% DurationCPUTime: 0.80s
% Computational Cost: add. (639->146), mult. (1651->200), div. (0->0), fcn. (2032->10), ass. (0->92)
t135 = cos(qJ(1));
t79 = sin(pkin(6));
t117 = t79 * t135;
t120 = cos(pkin(6));
t105 = t120 * t135;
t82 = sin(qJ(2));
t83 = sin(qJ(1));
t86 = cos(qJ(2));
t60 = t82 * t105 + t83 * t86;
t81 = sin(qJ(3));
t85 = cos(qJ(3));
t32 = t85 * t117 + t60 * t81;
t59 = -t86 * t105 + t83 * t82;
t80 = sin(qJ(5));
t84 = cos(qJ(5));
t10 = t32 * t80 + t59 * t84;
t11 = t32 * t84 - t59 * t80;
t33 = -t81 * t117 + t60 * t85;
t129 = t79 * t83;
t113 = t83 * t120;
t62 = -t82 * t113 + t135 * t86;
t37 = t81 * t129 + t62 * t85;
t128 = t79 * t85;
t58 = t120 * t81 + t82 * t128;
t95 = g(1) * t37 + g(2) * t33 + g(3) * t58;
t138 = t32 * pkin(10);
t36 = -t83 * t128 + t62 * t81;
t137 = t36 * pkin(10);
t130 = t79 * t82;
t57 = -t120 * t85 + t81 * t130;
t136 = t57 * pkin(10);
t132 = t59 * t85;
t61 = t86 * t113 + t135 * t82;
t131 = t61 * t85;
t127 = t79 * t86;
t126 = t80 * t81;
t125 = t80 * t86;
t124 = t81 * t84;
t123 = pkin(2) * t127 + pkin(9) * t130;
t122 = t135 * pkin(1) + pkin(8) * t129;
t121 = qJ(4) * t81;
t119 = t85 * t127;
t118 = t84 * t127;
t116 = -t83 * pkin(1) + pkin(8) * t117;
t115 = -t59 * pkin(2) + t60 * pkin(9);
t114 = -t61 * pkin(2) + t62 * pkin(9);
t26 = t32 * pkin(3);
t112 = t33 * qJ(4) - t26;
t28 = t36 * pkin(3);
t111 = t37 * qJ(4) - t28;
t48 = t57 * pkin(3);
t110 = t58 * qJ(4) - t48;
t109 = pkin(3) * t119 + t121 * t127 + t123;
t13 = -t36 * t84 + t61 * t80;
t108 = g(1) * t11 + g(2) * t13;
t107 = -g(1) * t32 + g(2) * t36;
t106 = -g(1) * t33 + g(2) * t37;
t22 = g(1) * t59 - g(2) * t61;
t104 = -pkin(3) * t132 - t59 * t121 + t115;
t103 = -pkin(3) * t131 - t61 * t121 + t114;
t102 = t62 * pkin(2) + t61 * pkin(9) + t122;
t101 = pkin(4) * t130 + pkin(10) * t119 + t109;
t100 = g(1) * t135 + g(2) * t83;
t98 = -t60 * pkin(2) - t59 * pkin(9) + t116;
t30 = t79 * t125 + t57 * t84;
t1 = g(1) * t13 - g(2) * t11 - g(3) * t30;
t14 = t36 * t80 + t61 * t84;
t31 = -t57 * t80 + t118;
t97 = g(1) * t14 + g(2) * t10 - g(3) * t31;
t18 = t59 * t124 + t60 * t80;
t20 = t61 * t124 + t62 * t80;
t40 = -t81 * t118 + t80 * t130;
t96 = g(1) * t20 + g(2) * t18 + g(3) * t40;
t7 = g(1) * t36 + g(2) * t32 + g(3) * t57;
t94 = t60 * pkin(4) - pkin(10) * t132 + t104;
t93 = t62 * pkin(4) - pkin(10) * t131 + t103;
t92 = -g(1) * t61 - g(2) * t59 + g(3) * t127;
t91 = g(1) * t62 + g(2) * t60 + g(3) * t130;
t90 = t37 * pkin(3) + t36 * qJ(4) + t102;
t89 = -pkin(3) * t33 - qJ(4) * t32 + t98;
t88 = t61 * pkin(4) + t37 * pkin(10) + t90;
t87 = -t59 * pkin(4) - pkin(10) * t33 + t89;
t41 = (t81 * t125 + t82 * t84) * t79;
t21 = -t61 * t126 + t62 * t84;
t19 = -t59 * t126 + t60 * t84;
t16 = t92 * t85;
t15 = t92 * t81;
t5 = t95 * t84;
t4 = t95 * t80;
t3 = g(1) * t10 - g(2) * t14;
t2 = -g(1) * t21 - g(2) * t19 - g(3) * t41;
t6 = [0, 0, 0, 0, 0, 0, g(1) * t83 - g(2) * t135, t100, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t60 - g(2) * t62, -t22, -t100 * t79, -g(1) * t116 - g(2) * t122, 0, 0, 0, 0, 0, 0, -t106, t107, t22, -g(1) * t98 - g(2) * t102, 0, 0, 0, 0, 0, 0, t22, t106, -t107, -g(1) * t89 - g(2) * t90, 0, 0, 0, 0, 0, 0, t3, t108, -t106, -g(1) * t87 - g(2) * t88, 0, 0, 0, 0, 0, 0, t3, -t106, -t108, -g(1) * (-pkin(5) * t10 + t11 * qJ(6) + t87) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t88); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t91, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, -t91, -g(1) * t114 - g(2) * t115 - g(3) * t123, 0, 0, 0, 0, 0, 0, -t91, t16, -t15, -g(1) * t103 - g(2) * t104 - g(3) * t109, 0, 0, 0, 0, 0, 0, t2, t96, -t16, -g(1) * t93 - g(2) * t94 - g(3) * t101, 0, 0, 0, 0, 0, 0, t2, -t16, -t96, -g(1) * (t21 * pkin(5) + t20 * qJ(6) + t93) - g(2) * (t19 * pkin(5) + t18 * qJ(6) + t94) - g(3) * (t41 * pkin(5) + t40 * qJ(6) + t101); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t95, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t95, -g(1) * t111 - g(2) * t112 - g(3) * t110, 0, 0, 0, 0, 0, 0, -t4, -t5, t7, -g(1) * (t111 - t137) - g(2) * (t112 - t138) - g(3) * (t110 - t136) 0, 0, 0, 0, 0, 0, -t4, t7, t5, -g(1) * (-t28 - t137) - g(2) * (-t26 - t138) - g(3) * (-t48 - t136) - t95 * (pkin(5) * t80 - qJ(6) * t84 + qJ(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t97, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t97, -g(1) * (-t13 * pkin(5) + t14 * qJ(6)) - g(2) * (pkin(5) * t11 + t10 * qJ(6)) - g(3) * (t30 * pkin(5) - t31 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t6;

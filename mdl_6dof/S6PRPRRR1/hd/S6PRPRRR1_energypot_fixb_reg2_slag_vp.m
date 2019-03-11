% Calculate inertial parameters regressor of potential energy for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:04
% EndTime: 2019-03-08 20:25:05
% DurationCPUTime: 0.27s
% Computational Cost: add. (359->99), mult. (712->156), div. (0->0), fcn. (880->14), ass. (0->57)
t100 = sin(pkin(12));
t103 = cos(pkin(12));
t108 = sin(qJ(2));
t111 = cos(qJ(2));
t84 = -t108 * t100 + t111 * t103;
t101 = sin(pkin(11));
t104 = cos(pkin(11));
t102 = sin(pkin(6));
t105 = cos(pkin(6));
t127 = t105 * t108;
t82 = pkin(2) * t127 + (-pkin(7) - qJ(3)) * t102;
t94 = t111 * pkin(2) + pkin(1);
t135 = t101 * t94 + t104 * t82;
t134 = t105 * pkin(7) + qJ(1);
t133 = t101 * t102;
t132 = t102 * t104;
t107 = sin(qJ(4));
t131 = t102 * t107;
t130 = t102 * t108;
t110 = cos(qJ(4));
t129 = t102 * t110;
t128 = t105 * t107;
t126 = t105 * t111;
t123 = t101 * t131;
t122 = t104 * t131;
t121 = -t101 * t82 + t104 * t94;
t120 = pkin(2) * t130 + t105 * qJ(3) + t134;
t119 = g(1) * t101 - g(2) * t104;
t118 = t111 * t100 + t108 * t103;
t112 = -pkin(9) - pkin(8);
t114 = t84 * t105;
t70 = -t101 * t114 - t104 * t118;
t81 = t118 * t105;
t71 = -t101 * t81 + t104 * t84;
t93 = t110 * pkin(4) + pkin(3);
t117 = pkin(4) * t123 + t70 * t112 + t71 * t93 + t121;
t79 = t84 * t102;
t80 = t118 * t102;
t116 = pkin(4) * t128 + t79 * t112 + t80 * t93 + t120;
t69 = t101 * t84 + t104 * t81;
t99 = qJ(4) + qJ(5);
t96 = sin(t99);
t97 = cos(t99);
t60 = t97 * t132 + t69 * t96;
t62 = -t97 * t133 + t71 * t96;
t72 = -t105 * t97 + t80 * t96;
t115 = g(1) * t62 + g(2) * t60 + g(3) * t72;
t68 = -t101 * t118 + t104 * t114;
t59 = g(1) * t70 + g(2) * t68 + g(3) * t79;
t113 = -pkin(4) * t122 + t68 * t112 + t69 * t93 + t135;
t109 = cos(qJ(6));
t106 = sin(qJ(6));
t78 = -g(3) * t105 - t119 * t102;
t73 = t105 * t96 + t80 * t97;
t63 = t96 * t133 + t71 * t97;
t61 = -t96 * t132 + t69 * t97;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t101, t119, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * (-t101 * t127 + t104 * t111) - g(2) * (t101 * t111 + t104 * t127) - g(3) * t130, -g(1) * (-t101 * t126 - t104 * t108) - g(2) * (-t101 * t108 + t104 * t126) - g(3) * t102 * t111, t78, -g(1) * (t104 * pkin(1) + pkin(7) * t133) - g(2) * (t101 * pkin(1) - pkin(7) * t132) - g(3) * t134, 0, 0, 0, 0, 0, 0, -g(1) * t71 - g(2) * t69 - g(3) * t80, -t59, t78, -g(1) * t121 - g(2) * t135 - g(3) * t120, 0, 0, 0, 0, 0, 0, -g(1) * (t71 * t110 + t123) - g(2) * (t69 * t110 - t122) - g(3) * (t80 * t110 + t128) -g(1) * (t101 * t129 - t71 * t107) - g(2) * (-t104 * t129 - t69 * t107) - g(3) * (t105 * t110 - t80 * t107) t59, -g(1) * (t71 * pkin(3) - t70 * pkin(8) + t121) - g(2) * (t69 * pkin(3) - t68 * pkin(8) + t135) - g(3) * (t80 * pkin(3) - t79 * pkin(8) + t120) 0, 0, 0, 0, 0, 0, -g(1) * t63 - g(2) * t61 - g(3) * t73, t115, t59, -g(1) * t117 - g(2) * t113 - g(3) * t116, 0, 0, 0, 0, 0, 0, -g(1) * (-t70 * t106 + t63 * t109) - g(2) * (-t68 * t106 + t61 * t109) - g(3) * (-t79 * t106 + t73 * t109) -g(1) * (-t63 * t106 - t70 * t109) - g(2) * (-t61 * t106 - t68 * t109) - g(3) * (-t73 * t106 - t79 * t109) -t115, -g(1) * (t63 * pkin(5) + t62 * pkin(10) + t117) - g(2) * (t61 * pkin(5) + t60 * pkin(10) + t113) - g(3) * (t73 * pkin(5) + t72 * pkin(10) + t116);];
U_reg  = t1;

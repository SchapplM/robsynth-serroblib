% Calculate inertial parameters regressor of potential energy for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:45:21
% EndTime: 2019-03-09 00:45:21
% DurationCPUTime: 0.28s
% Computational Cost: add. (331->106), mult. (533->155), div. (0->0), fcn. (633->14), ass. (0->54)
t104 = cos(pkin(12));
t102 = sin(pkin(12));
t103 = sin(pkin(6));
t130 = t102 * t103;
t134 = t104 * pkin(1) + pkin(7) * t130;
t106 = sin(qJ(5));
t108 = sin(qJ(2));
t105 = cos(pkin(6));
t111 = cos(qJ(2));
t122 = t105 * t111;
t79 = t102 * t108 - t104 * t122;
t133 = t79 * t106;
t81 = t102 * t122 + t104 * t108;
t132 = t81 * t106;
t131 = t105 * pkin(7) + qJ(1);
t129 = t103 * t104;
t107 = sin(qJ(3));
t128 = t103 * t107;
t127 = t103 * t108;
t110 = cos(qJ(3));
t126 = t103 * t110;
t125 = t103 * t111;
t124 = t105 * t107;
t123 = t105 * t108;
t121 = t102 * t128;
t120 = t106 * t125;
t97 = t102 * pkin(1);
t119 = -pkin(7) * t129 + t97;
t113 = -pkin(9) - pkin(8);
t82 = -t102 * t123 + t104 * t111;
t92 = t110 * pkin(3) + pkin(2);
t118 = pkin(3) * t121 - t81 * t113 + t82 * t92 + t134;
t117 = pkin(3) * t124 + t113 * t125 + t92 * t127 + t131;
t116 = g(1) * t102 - g(2) * t104;
t80 = t102 * t111 + t104 * t123;
t101 = qJ(3) + qJ(4);
t94 = sin(t101);
t96 = cos(t101);
t69 = t96 * t129 + t80 * t94;
t71 = -t96 * t130 + t82 * t94;
t75 = -t105 * t96 + t94 * t127;
t115 = g(1) * t71 + g(2) * t69 + g(3) * t75;
t68 = -g(1) * t81 - g(2) * t79 + g(3) * t125;
t114 = t80 * t92 - t79 * t113 + t97 + (-pkin(3) * t107 - pkin(7)) * t129;
t112 = -pkin(11) - pkin(10);
t109 = cos(qJ(5));
t100 = qJ(5) + qJ(6);
t95 = cos(t100);
t93 = sin(t100);
t91 = t109 * pkin(5) + pkin(4);
t76 = t105 * t94 + t96 * t127;
t72 = t94 * t130 + t82 * t96;
t70 = -t94 * t129 + t80 * t96;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t102, t116, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t80 - g(3) * t127, -t68, -g(3) * t105 - t116 * t103, -g(1) * t134 - g(2) * t119 - g(3) * t131, 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t110 + t121) - g(2) * (-t104 * t128 + t80 * t110) - g(3) * (t108 * t126 + t124) -g(1) * (t102 * t126 - t82 * t107) - g(2) * (-t104 * t126 - t80 * t107) - g(3) * (t105 * t110 - t107 * t127) t68, -g(1) * (t82 * pkin(2) + t81 * pkin(8) + t134) - g(2) * (t80 * pkin(2) + t79 * pkin(8) + t119) - g(3) * ((pkin(2) * t108 - pkin(8) * t111) * t103 + t131) 0, 0, 0, 0, 0, 0, -g(1) * t72 - g(2) * t70 - g(3) * t76, t115, t68, -g(1) * t118 - g(2) * t114 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t109 + t132) - g(2) * (t70 * t109 + t133) - g(3) * (t76 * t109 - t120) -g(1) * (-t72 * t106 + t81 * t109) - g(2) * (-t70 * t106 + t79 * t109) - g(3) * (-t76 * t106 - t109 * t125) -t115, -g(1) * (t72 * pkin(4) + t71 * pkin(10) + t118) - g(2) * (t70 * pkin(4) + t69 * pkin(10) + t114) - g(3) * (t76 * pkin(4) + t75 * pkin(10) + t117) 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t95 + t81 * t93) - g(2) * (t70 * t95 + t79 * t93) - g(3) * (-t93 * t125 + t76 * t95) -g(1) * (-t72 * t93 + t81 * t95) - g(2) * (-t70 * t93 + t79 * t95) - g(3) * (-t95 * t125 - t76 * t93) -t115, -g(1) * (pkin(5) * t132 - t71 * t112 + t72 * t91 + t118) - g(2) * (pkin(5) * t133 - t69 * t112 + t70 * t91 + t114) - g(3) * (-pkin(5) * t120 - t75 * t112 + t76 * t91 + t117);];
U_reg  = t1;

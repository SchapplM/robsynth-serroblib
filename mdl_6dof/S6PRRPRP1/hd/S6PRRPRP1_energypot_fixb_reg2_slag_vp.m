% Calculate inertial parameters regressor of potential energy for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:26:33
% EndTime: 2019-03-08 21:26:34
% DurationCPUTime: 0.25s
% Computational Cost: add. (319->95), mult. (533->137), div. (0->0), fcn. (633->12), ass. (0->53)
t101 = cos(pkin(10));
t100 = sin(pkin(6));
t99 = sin(pkin(10));
t130 = t100 * t99;
t131 = t101 * pkin(1) + pkin(7) * t130;
t105 = sin(qJ(5));
t107 = sin(qJ(2));
t102 = cos(pkin(6));
t110 = cos(qJ(2));
t119 = t102 * t110;
t79 = -t101 * t119 + t99 * t107;
t129 = t79 * t105;
t81 = t101 * t107 + t99 * t119;
t128 = t81 * t105;
t127 = t102 * pkin(7) + qJ(1);
t126 = t100 * t101;
t106 = sin(qJ(3));
t125 = t100 * t106;
t124 = t100 * t107;
t109 = cos(qJ(3));
t123 = t100 * t109;
t122 = t100 * t110;
t121 = t102 * t106;
t120 = t102 * t107;
t118 = t99 * t125;
t117 = t105 * t122;
t95 = t99 * pkin(1);
t116 = -pkin(7) * t126 + t95;
t104 = -qJ(4) - pkin(8);
t82 = t101 * t110 - t99 * t120;
t92 = t109 * pkin(3) + pkin(2);
t115 = pkin(3) * t118 - t81 * t104 + t82 * t92 + t131;
t114 = pkin(3) * t121 + t104 * t122 + t92 * t124 + t127;
t113 = g(1) * t99 - g(2) * t101;
t80 = t101 * t120 + t99 * t110;
t98 = qJ(3) + pkin(11);
t93 = sin(t98);
t94 = cos(t98);
t69 = t94 * t126 + t80 * t93;
t71 = -t94 * t130 + t82 * t93;
t75 = -t102 * t94 + t93 * t124;
t112 = g(1) * t71 + g(2) * t69 + g(3) * t75;
t68 = -g(1) * t81 - g(2) * t79 + g(3) * t122;
t111 = t80 * t92 - t79 * t104 + t95 + (-pkin(3) * t106 - pkin(7)) * t126;
t108 = cos(qJ(5));
t103 = -qJ(6) - pkin(9);
t91 = t108 * pkin(5) + pkin(4);
t76 = t102 * t93 + t94 * t124;
t72 = t93 * t130 + t82 * t94;
t70 = -t93 * t126 + t80 * t94;
t66 = -g(1) * (t72 * t108 + t128) - g(2) * (t70 * t108 + t129) - g(3) * (t76 * t108 - t117);
t65 = -g(1) * (-t72 * t105 + t81 * t108) - g(2) * (-t70 * t105 + t79 * t108) - g(3) * (-t76 * t105 - t108 * t122);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t101 - g(2) * t99, t113, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t80 - g(3) * t124, -t68, -g(3) * t102 - t113 * t100, -g(1) * t131 - g(2) * t116 - g(3) * t127, 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t109 + t118) - g(2) * (-t101 * t125 + t80 * t109) - g(3) * (t107 * t123 + t121) -g(1) * (-t82 * t106 + t99 * t123) - g(2) * (-t101 * t123 - t80 * t106) - g(3) * (t102 * t109 - t106 * t124) t68, -g(1) * (t82 * pkin(2) + t81 * pkin(8) + t131) - g(2) * (t80 * pkin(2) + t79 * pkin(8) + t116) - g(3) * ((pkin(2) * t107 - pkin(8) * t110) * t100 + t127) 0, 0, 0, 0, 0, 0, -g(1) * t72 - g(2) * t70 - g(3) * t76, t112, t68, -g(1) * t115 - g(2) * t111 - g(3) * t114, 0, 0, 0, 0, 0, 0, t66, t65, -t112, -g(1) * (t72 * pkin(4) + t71 * pkin(9) + t115) - g(2) * (t70 * pkin(4) + t69 * pkin(9) + t111) - g(3) * (t76 * pkin(4) + t75 * pkin(9) + t114) 0, 0, 0, 0, 0, 0, t66, t65, -t112, -g(1) * (pkin(5) * t128 - t71 * t103 + t72 * t91 + t115) - g(2) * (pkin(5) * t129 - t69 * t103 + t70 * t91 + t111) - g(3) * (-pkin(5) * t117 - t75 * t103 + t76 * t91 + t114);];
U_reg  = t1;

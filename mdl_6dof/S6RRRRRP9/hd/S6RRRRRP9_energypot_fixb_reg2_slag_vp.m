% Calculate inertial parameters regressor of potential energy for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:11:38
% EndTime: 2019-03-10 02:11:38
% DurationCPUTime: 0.22s
% Computational Cost: add. (299->97), mult. (592->136), div. (0->0), fcn. (717->12), ass. (0->51)
t108 = sin(qJ(4));
t127 = t108 * pkin(4) + pkin(9);
t115 = -pkin(11) - pkin(10);
t106 = qJ(4) + qJ(5);
t99 = sin(t106);
t139 = pkin(5) * t99 + t127;
t138 = cos(qJ(3));
t134 = cos(pkin(6));
t136 = t134 * pkin(8) + pkin(7);
t112 = cos(qJ(4));
t98 = t112 * pkin(4) + pkin(3);
t114 = cos(qJ(1));
t107 = sin(pkin(6));
t111 = sin(qJ(1));
t132 = t107 * t111;
t135 = t114 * pkin(1) + pkin(8) * t132;
t110 = sin(qJ(2));
t133 = t107 * t110;
t113 = cos(qJ(2));
t131 = t107 * t113;
t130 = t107 * t114;
t129 = pkin(2) * t133 + t136;
t125 = t111 * t134;
t90 = -t110 * t125 + t114 * t113;
t128 = t90 * pkin(2) + t135;
t126 = t107 * t138;
t124 = t114 * t134;
t123 = t111 * pkin(1) - pkin(8) * t130;
t122 = g(1) * t111 - g(2) * t114;
t89 = t114 * t110 + t113 * t125;
t121 = t89 * pkin(9) + t128;
t88 = t110 * t124 + t111 * t113;
t120 = t88 * pkin(2) + t123;
t119 = -pkin(9) * t131 + t129;
t109 = sin(qJ(3));
t79 = t88 * t109 + t114 * t126;
t81 = t90 * t109 - t111 * t126;
t85 = t109 * t133 - t134 * t138;
t118 = g(1) * t81 + g(2) * t79 + g(3) * t85;
t87 = t111 * t110 - t113 * t124;
t117 = t87 * pkin(9) + t120;
t116 = -g(1) * t89 - g(2) * t87 + g(3) * t131;
t105 = -qJ(6) + t115;
t100 = cos(t106);
t91 = pkin(5) * t100 + t98;
t86 = t134 * t109 + t110 * t126;
t82 = t109 * t132 + t90 * t138;
t80 = -t109 * t130 + t88 * t138;
t77 = -g(1) * (t82 * t100 + t89 * t99) - g(2) * (t80 * t100 + t87 * t99) - g(3) * (t86 * t100 - t99 * t131);
t76 = -g(1) * (t89 * t100 - t82 * t99) - g(2) * (t87 * t100 - t80 * t99) - g(3) * (-t100 * t131 - t86 * t99);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t114 - g(2) * t111, t122, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t90 - g(2) * t88 - g(3) * t133, -t116, -g(3) * t134 - t122 * t107, -g(1) * t135 - g(2) * t123 - g(3) * t136, 0, 0, 0, 0, 0, 0, -g(1) * t82 - g(2) * t80 - g(3) * t86, t118, t116, -g(1) * t121 - g(2) * t117 - g(3) * t119, 0, 0, 0, 0, 0, 0, -g(1) * (t89 * t108 + t82 * t112) - g(2) * (t87 * t108 + t80 * t112) - g(3) * (-t108 * t131 + t86 * t112) -g(1) * (-t82 * t108 + t89 * t112) - g(2) * (-t80 * t108 + t87 * t112) - g(3) * (-t86 * t108 - t112 * t131) -t118, -g(1) * (t82 * pkin(3) + t81 * pkin(10) + t121) - g(2) * (t80 * pkin(3) + t79 * pkin(10) + t117) - g(3) * (t86 * pkin(3) + t85 * pkin(10) + t119) 0, 0, 0, 0, 0, 0, t77, t76, -t118, -g(1) * (-t81 * t115 + t127 * t89 + t82 * t98 + t128) - g(2) * (-t79 * t115 + t127 * t87 + t80 * t98 + t120) - g(3) * (-t85 * t115 - t127 * t131 + t86 * t98 + t129) 0, 0, 0, 0, 0, 0, t77, t76, -t118, -g(1) * (-t81 * t105 + t139 * t89 + t82 * t91 + t128) - g(2) * (-t79 * t105 + t139 * t87 + t80 * t91 + t120) - g(3) * (-t85 * t105 - t139 * t131 + t86 * t91 + t129);];
U_reg  = t1;

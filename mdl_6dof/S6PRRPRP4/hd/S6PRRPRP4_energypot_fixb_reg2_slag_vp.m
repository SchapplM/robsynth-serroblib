% Calculate inertial parameters regressor of potential energy for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:43:57
% EndTime: 2019-03-08 21:43:57
% DurationCPUTime: 0.22s
% Computational Cost: add. (265->85), mult. (600->113), div. (0->0), fcn. (727->10), ass. (0->53)
t130 = pkin(4) + pkin(8);
t93 = sin(pkin(6));
t129 = pkin(7) * t93;
t100 = cos(qJ(2));
t119 = cos(pkin(6));
t112 = t100 * t119;
t92 = sin(pkin(10));
t94 = cos(pkin(10));
t98 = sin(qJ(2));
t76 = -t94 * t112 + t92 * t98;
t128 = t76 * pkin(8);
t78 = t92 * t112 + t94 * t98;
t127 = t78 * pkin(8);
t99 = cos(qJ(5));
t126 = t99 * pkin(5) + t130;
t125 = cos(qJ(3));
t97 = sin(qJ(3));
t124 = t93 * t97;
t123 = t93 * t98;
t122 = t94 * pkin(1) + t92 * t129;
t121 = t100 * t93;
t120 = t119 * pkin(7) + qJ(1);
t118 = pkin(8) * t121;
t113 = t98 * t119;
t79 = t94 * t100 - t92 * t113;
t117 = t79 * pkin(2) + t122;
t116 = pkin(2) * t123 + t120;
t115 = t93 * t125;
t96 = sin(qJ(5));
t114 = pkin(5) * t96 + qJ(4);
t72 = t92 * t124 + t79 * t125;
t111 = t72 * pkin(3) + t117;
t81 = t98 * t115 + t119 * t97;
t110 = t81 * pkin(3) + t116;
t109 = t92 * pkin(1) - t94 * t129;
t108 = g(1) * t92 - g(2) * t94;
t77 = t92 * t100 + t94 * t113;
t107 = t77 * pkin(2) + t109;
t70 = -t94 * t124 + t77 * t125;
t106 = t70 * pkin(3) + t107;
t71 = -t92 * t115 + t79 * t97;
t105 = t71 * qJ(4) + t111;
t80 = -t119 * t125 + t97 * t123;
t104 = t80 * qJ(4) + t110;
t69 = t94 * t115 + t77 * t97;
t103 = g(1) * t71 + g(2) * t69 + g(3) * t80;
t102 = g(1) * t72 + g(2) * t70 + g(3) * t81;
t66 = -g(1) * t78 - g(2) * t76 + g(3) * t121;
t101 = t69 * qJ(4) + t106;
t95 = -qJ(6) - pkin(9);
t64 = -g(1) * (t71 * t96 + t78 * t99) - g(2) * (t69 * t96 + t76 * t99) - g(3) * (-t99 * t121 + t80 * t96);
t63 = -g(1) * (t71 * t99 - t78 * t96) - g(2) * (t69 * t99 - t76 * t96) - g(3) * (t96 * t121 + t80 * t99);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t94 - g(2) * t92, t108, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t79 - g(2) * t77 - g(3) * t123, -t66, -g(3) * t119 - t108 * t93, -g(1) * t122 - g(2) * t109 - g(3) * t120, 0, 0, 0, 0, 0, 0, -t102, t103, t66, -g(1) * (t117 + t127) - g(2) * (t107 + t128) - g(3) * (t116 - t118) 0, 0, 0, 0, 0, 0, t66, t102, -t103, -g(1) * (t105 + t127) - g(2) * (t101 + t128) - g(3) * (t104 - t118) 0, 0, 0, 0, 0, 0, t64, t63, -t102, -g(1) * (t72 * pkin(9) + t130 * t78 + t105) - g(2) * (t70 * pkin(9) + t130 * t76 + t101) - g(3) * (t81 * pkin(9) - t130 * t121 + t104) 0, 0, 0, 0, 0, 0, t64, t63, -t102, -g(1) * (t114 * t71 + t126 * t78 - t72 * t95 + t111) - g(2) * (t114 * t69 + t126 * t76 - t70 * t95 + t106) - g(3) * (t114 * t80 - t126 * t121 - t81 * t95 + t110);];
U_reg  = t1;

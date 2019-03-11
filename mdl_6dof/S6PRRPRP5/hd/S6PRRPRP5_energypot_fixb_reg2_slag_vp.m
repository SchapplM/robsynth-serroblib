% Calculate inertial parameters regressor of potential energy for
% S6PRRPRP5
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:15
% EndTime: 2019-03-08 21:49:15
% DurationCPUTime: 0.18s
% Computational Cost: add. (277->78), mult. (638->108), div. (0->0), fcn. (780->10), ass. (0->53)
t132 = cos(qJ(3));
t100 = sin(pkin(6));
t131 = pkin(7) * t100;
t101 = cos(pkin(10));
t99 = sin(pkin(10));
t130 = t101 * pkin(1) + t99 * t131;
t128 = cos(pkin(6));
t129 = t128 * pkin(7) + qJ(1);
t103 = sin(qJ(3));
t127 = t100 * t103;
t104 = sin(qJ(2));
t126 = t100 * t104;
t106 = cos(qJ(2));
t125 = t100 * t106;
t124 = pkin(8) * t125;
t123 = pkin(2) * t126 + t129;
t122 = t100 * t132;
t121 = t104 * t128;
t120 = t106 * t128;
t119 = t99 * pkin(1) - t101 * t131;
t118 = g(1) * t99 - g(2) * t101;
t86 = t101 * t104 + t99 * t120;
t87 = t101 * t106 - t99 * t121;
t117 = t87 * pkin(2) + t86 * pkin(8) + t130;
t88 = t103 * t126 - t128 * t132;
t89 = t128 * t103 + t104 * t122;
t116 = t89 * pkin(3) + t88 * qJ(4) + t123;
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t85 = t101 * t121 + t99 * t106;
t73 = t101 * t122 + t85 * t103;
t84 = -t101 * t120 + t99 * t104;
t64 = t84 * t102 - t73 * t105;
t75 = t87 * t103 - t99 * t122;
t66 = t86 * t102 - t75 * t105;
t77 = t102 * t125 + t88 * t105;
t115 = g(1) * t66 + g(2) * t64 - g(3) * t77;
t114 = g(1) * t75 + g(2) * t73 + g(3) * t88;
t74 = -t101 * t127 + t85 * t132;
t76 = t99 * t127 + t87 * t132;
t113 = g(1) * t76 + g(2) * t74 + g(3) * t89;
t112 = t85 * pkin(2) + t84 * pkin(8) + t119;
t68 = -g(1) * t86 - g(2) * t84 + g(3) * t125;
t111 = t76 * pkin(3) + t75 * qJ(4) + t117;
t110 = t74 * pkin(3) + t73 * qJ(4) + t112;
t109 = t86 * pkin(4) + t76 * pkin(9) + t111;
t108 = t89 * pkin(9) + (-pkin(4) - pkin(8)) * t125 + t116;
t107 = t84 * pkin(4) + t74 * pkin(9) + t110;
t78 = t88 * t102 - t105 * t125;
t67 = t75 * t102 + t86 * t105;
t65 = t73 * t102 + t84 * t105;
t62 = -g(1) * t67 - g(2) * t65 - g(3) * t78;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t101 - g(2) * t99, t118, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t87 - g(2) * t85 - g(3) * t126, -t68, -g(3) * t128 - t118 * t100, -g(1) * t130 - g(2) * t119 - g(3) * t129, 0, 0, 0, 0, 0, 0, -t113, t114, t68, -g(1) * t117 - g(2) * t112 - g(3) * (t123 - t124) 0, 0, 0, 0, 0, 0, t68, t113, -t114, -g(1) * t111 - g(2) * t110 - g(3) * (t116 - t124) 0, 0, 0, 0, 0, 0, t62, t115, -t113, -g(1) * t109 - g(2) * t107 - g(3) * t108, 0, 0, 0, 0, 0, 0, t62, -t113, -t115, -g(1) * (t67 * pkin(5) + t66 * qJ(6) + t109) - g(2) * (t65 * pkin(5) + t64 * qJ(6) + t107) - g(3) * (t78 * pkin(5) - t77 * qJ(6) + t108);];
U_reg  = t1;

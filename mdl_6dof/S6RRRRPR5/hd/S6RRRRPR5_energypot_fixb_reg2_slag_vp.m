% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:16:16
% EndTime: 2019-03-09 22:16:16
% DurationCPUTime: 0.16s
% Computational Cost: add. (215->70), mult. (262->92), div. (0->0), fcn. (278->10), ass. (0->44)
t105 = cos(qJ(4));
t99 = qJ(2) + qJ(3);
t95 = sin(t99);
t123 = t105 * t95;
t101 = sin(qJ(4));
t125 = t101 * t95;
t131 = pkin(4) * t123 + qJ(5) * t125;
t130 = g(3) * pkin(6);
t96 = cos(t99);
t129 = pkin(3) * t96;
t128 = g(3) * t95;
t102 = sin(qJ(2));
t127 = t102 * pkin(2) + pkin(6);
t103 = sin(qJ(1));
t107 = cos(qJ(1));
t108 = -pkin(8) - pkin(7);
t106 = cos(qJ(2));
t93 = t106 * pkin(2) + pkin(1);
t126 = t103 * t93 + t107 * t108;
t124 = t103 * t95;
t122 = t107 * t95;
t121 = t103 * t101;
t120 = t103 * t105;
t119 = t107 * t101;
t118 = t107 * t105;
t117 = t95 * pkin(3) + t127;
t116 = pkin(9) * t124 + t103 * t129 + t126;
t115 = -t103 * t108 + t107 * t93;
t114 = g(1) * t107 + g(2) * t103;
t113 = -t96 * pkin(9) + t117;
t112 = pkin(9) * t122 + t107 * t129 + t115;
t77 = t96 * t121 + t118;
t78 = t96 * t120 - t119;
t111 = t78 * pkin(4) + t77 * qJ(5) + t116;
t79 = t96 * t119 - t120;
t110 = g(1) * t79 + g(2) * t77 + g(3) * t125;
t80 = t96 * t118 + t121;
t109 = t80 * pkin(4) + t79 * qJ(5) + t112;
t104 = cos(qJ(6));
t100 = sin(qJ(6));
t81 = g(1) * t103 - g(2) * t107;
t74 = -g(3) * t96 + t114 * t95;
t73 = -g(1) * t80 - g(2) * t78 - g(3) * t123;
t1 = [0, 0, 0, 0, 0, 0, -t114, t81, -g(3), -t130, 0, 0, 0, 0, 0, 0, -g(3) * t102 - t114 * t106, -g(3) * t106 + t114 * t102, -t81, -g(1) * (t107 * pkin(1) + t103 * pkin(7)) - g(2) * (t103 * pkin(1) - t107 * pkin(7)) - t130, 0, 0, 0, 0, 0, 0, -t114 * t96 - t128, t74, -t81, -g(1) * t115 - g(2) * t126 - g(3) * t127, 0, 0, 0, 0, 0, 0, t73, t110, -t74, -g(1) * t112 - g(2) * t116 - g(3) * t113, 0, 0, 0, 0, 0, 0, t73, -t74, -t110, -g(1) * t109 - g(2) * t111 - g(3) * (t113 + t131) 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t100 + t80 * t104) - g(2) * (t77 * t100 + t78 * t104) - (t100 * t101 + t104 * t105) * t128, -g(1) * (-t80 * t100 + t79 * t104) - g(2) * (-t78 * t100 + t77 * t104) - (-t100 * t105 + t101 * t104) * t128, t74, -g(1) * (t80 * pkin(5) - pkin(10) * t122 + t109) - g(2) * (t78 * pkin(5) - pkin(10) * t124 + t111) - g(3) * (pkin(5) * t123 + (-pkin(9) + pkin(10)) * t96 + t117 + t131);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:20:57
% EndTime: 2019-03-09 09:20:58
% DurationCPUTime: 0.18s
% Computational Cost: add. (208->80), mult. (457->106), div. (0->0), fcn. (530->10), ass. (0->46)
t100 = cos(pkin(6));
t131 = t100 * pkin(8) + pkin(7);
t108 = cos(qJ(1));
t104 = sin(qJ(1));
t99 = sin(pkin(6));
t128 = t104 * t99;
t130 = t108 * pkin(1) + pkin(8) * t128;
t103 = sin(qJ(2));
t129 = t103 * t99;
t107 = cos(qJ(2));
t127 = t107 * t99;
t126 = t108 * t99;
t125 = t104 * t103;
t124 = t104 * t107;
t123 = t108 * t103;
t122 = t108 * t107;
t121 = pkin(2) * t129 + t131;
t120 = qJ(3) * t127;
t119 = t104 * pkin(1) - pkin(8) * t126;
t85 = t100 * t124 + t123;
t86 = -t100 * t125 + t122;
t118 = t86 * pkin(2) + t85 * qJ(3) + t130;
t117 = pkin(3) * t129 - t100 * qJ(4) + t121;
t102 = sin(qJ(5));
t106 = cos(qJ(5));
t83 = -t100 * t122 + t125;
t70 = t83 * t102 - t106 * t126;
t72 = t85 * t102 + t106 * t128;
t81 = -t100 * t106 + t102 * t127;
t116 = g(1) * t72 + g(2) * t70 - g(3) * t81;
t84 = t100 * t123 + t124;
t115 = t84 * pkin(2) + t83 * qJ(3) + t119;
t114 = g(1) * t86 + g(2) * t84 + g(3) * t129;
t68 = -g(1) * t85 - g(2) * t83 + g(3) * t127;
t113 = t84 * pkin(3) + qJ(4) * t126 + t115;
t112 = t86 * pkin(3) - qJ(4) * t128 + t118;
t111 = pkin(9) * t129 + (-pkin(4) - qJ(3)) * t127 + t117;
t110 = t83 * pkin(4) + t84 * pkin(9) + t113;
t109 = t85 * pkin(4) + t86 * pkin(9) + t112;
t105 = cos(qJ(6));
t101 = sin(qJ(6));
t82 = -t100 * t102 - t106 * t127;
t74 = g(1) * t128 - g(2) * t126 + g(3) * t100;
t73 = -t102 * t128 + t85 * t106;
t71 = t102 * t126 + t83 * t106;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t108 - g(2) * t104, g(1) * t104 - g(2) * t108, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -t114, -t68, -t74, -g(1) * t130 - g(2) * t119 - g(3) * t131, 0, 0, 0, 0, 0, 0, -t114, -t74, t68, -g(1) * t118 - g(2) * t115 - g(3) * (-t120 + t121) 0, 0, 0, 0, 0, 0, t68, t114, t74, -g(1) * t112 - g(2) * t113 - g(3) * (t117 - t120) 0, 0, 0, 0, 0, 0, -g(1) * t73 - g(2) * t71 - g(3) * t82, t116, -t114, -g(1) * t109 - g(2) * t110 - g(3) * t111, 0, 0, 0, 0, 0, 0, -g(1) * (t86 * t101 + t73 * t105) - g(2) * (t84 * t101 + t71 * t105) - g(3) * (t101 * t129 + t82 * t105) -g(1) * (-t73 * t101 + t86 * t105) - g(2) * (-t71 * t101 + t84 * t105) - g(3) * (-t82 * t101 + t105 * t129) -t116, -g(1) * (t73 * pkin(5) + t72 * pkin(10) + t109) - g(2) * (t71 * pkin(5) + t70 * pkin(10) + t110) - g(3) * (t82 * pkin(5) - t81 * pkin(10) + t111);];
U_reg  = t1;

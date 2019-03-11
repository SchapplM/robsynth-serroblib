% Calculate inertial parameters regressor of potential energy for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:27:14
% EndTime: 2019-03-09 17:27:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (179->66), mult. (380->90), div. (0->0), fcn. (430->8), ass. (0->45)
t132 = g(3) * pkin(6);
t106 = sin(qJ(1));
t110 = cos(qJ(1));
t131 = t110 * pkin(1) + t106 * pkin(7);
t104 = sin(qJ(3));
t105 = sin(qJ(2));
t130 = t104 * t105;
t129 = t105 * t106;
t108 = cos(qJ(3));
t128 = t105 * t108;
t127 = t105 * t110;
t109 = cos(qJ(2));
t126 = t106 * t109;
t125 = t110 * t104;
t124 = t110 * t108;
t123 = t106 * pkin(1) - t110 * pkin(7);
t122 = t110 * t109 * pkin(2) + pkin(8) * t127 + t131;
t121 = t105 * pkin(2) - t109 * pkin(8) + pkin(6);
t120 = g(1) * t110 + g(2) * t106;
t119 = pkin(2) * t126 + pkin(8) * t129 + t123;
t118 = pkin(3) * t128 + qJ(4) * t130 + t121;
t103 = sin(qJ(5));
t107 = cos(qJ(5));
t83 = t104 * t126 + t124;
t84 = t108 * t126 - t125;
t70 = t84 * t103 - t83 * t107;
t85 = -t106 * t108 + t109 * t125;
t86 = t106 * t104 + t109 * t124;
t72 = t86 * t103 - t85 * t107;
t77 = t103 * t128 - t107 * t130;
t117 = g(1) * t72 + g(2) * t70 + g(3) * t77;
t116 = t86 * pkin(3) + t85 * qJ(4) + t122;
t115 = pkin(4) * t128 + t109 * pkin(9) + t118;
t114 = g(1) * t85 + g(2) * t83 + g(3) * t130;
t113 = t84 * pkin(3) + t83 * qJ(4) + t119;
t112 = t86 * pkin(4) - pkin(9) * t127 + t116;
t111 = t84 * pkin(4) - pkin(9) * t129 + t113;
t87 = g(1) * t106 - g(2) * t110;
t78 = (t103 * t104 + t107 * t108) * t105;
t74 = -g(3) * t109 + t120 * t105;
t73 = t85 * t103 + t86 * t107;
t71 = t83 * t103 + t84 * t107;
t69 = -g(1) * t86 - g(2) * t84 - g(3) * t128;
t68 = -g(1) * t73 - g(2) * t71 - g(3) * t78;
t1 = [0, 0, 0, 0, 0, 0, -t120, t87, -g(3), -t132, 0, 0, 0, 0, 0, 0, -g(3) * t105 - t120 * t109, t74, -t87, -g(1) * t131 - g(2) * t123 - t132, 0, 0, 0, 0, 0, 0, t69, t114, -t74, -g(1) * t122 - g(2) * t119 - g(3) * t121, 0, 0, 0, 0, 0, 0, t69, -t74, -t114, -g(1) * t116 - g(2) * t113 - g(3) * t118, 0, 0, 0, 0, 0, 0, t68, t117, t74, -g(1) * t112 - g(2) * t111 - g(3) * t115, 0, 0, 0, 0, 0, 0, t68, t74, -t117, -g(1) * (t73 * pkin(5) + t72 * qJ(6) + t112) - g(2) * (t71 * pkin(5) + t70 * qJ(6) + t111) - g(3) * (t78 * pkin(5) + t77 * qJ(6) + t115);];
U_reg  = t1;

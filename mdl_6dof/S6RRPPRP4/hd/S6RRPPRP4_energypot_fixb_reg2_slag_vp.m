% Calculate inertial parameters regressor of potential energy for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:40
% EndTime: 2019-03-09 08:39:40
% DurationCPUTime: 0.15s
% Computational Cost: add. (179->66), mult. (380->90), div. (0->0), fcn. (430->8), ass. (0->45)
t129 = g(3) * pkin(6);
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t128 = t107 * pkin(1) + t104 * pkin(7);
t100 = sin(pkin(9));
t103 = sin(qJ(2));
t127 = t100 * t103;
t101 = cos(pkin(9));
t126 = t101 * t103;
t125 = t103 * t104;
t124 = t103 * t107;
t106 = cos(qJ(2));
t123 = t104 * t106;
t122 = t107 * t100;
t121 = t107 * t101;
t120 = t104 * pkin(1) - t107 * pkin(7);
t119 = t107 * t106 * pkin(2) + qJ(3) * t124 + t128;
t118 = t103 * pkin(2) - t106 * qJ(3) + pkin(6);
t117 = g(1) * t107 + g(2) * t104;
t116 = pkin(2) * t123 + qJ(3) * t125 + t120;
t115 = pkin(3) * t126 + qJ(4) * t127 + t118;
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t80 = t100 * t123 + t121;
t81 = t101 * t123 - t122;
t67 = t81 * t102 - t80 * t105;
t82 = -t104 * t101 + t106 * t122;
t83 = t104 * t100 + t106 * t121;
t69 = t83 * t102 - t82 * t105;
t74 = t102 * t126 - t105 * t127;
t114 = g(1) * t69 + g(2) * t67 + g(3) * t74;
t113 = t83 * pkin(3) + t82 * qJ(4) + t119;
t112 = g(1) * t82 + g(2) * t80 + g(3) * t127;
t111 = pkin(4) * t126 + t106 * pkin(8) + t115;
t110 = t81 * pkin(3) + t80 * qJ(4) + t116;
t109 = t83 * pkin(4) - pkin(8) * t124 + t113;
t108 = t81 * pkin(4) - pkin(8) * t125 + t110;
t84 = g(1) * t104 - g(2) * t107;
t75 = (t100 * t102 + t101 * t105) * t103;
t73 = -g(3) * t106 + t117 * t103;
t70 = t82 * t102 + t83 * t105;
t68 = t80 * t102 + t81 * t105;
t66 = -g(1) * t83 - g(2) * t81 - g(3) * t126;
t65 = -g(1) * t70 - g(2) * t68 - g(3) * t75;
t1 = [0, 0, 0, 0, 0, 0, -t117, t84, -g(3), -t129, 0, 0, 0, 0, 0, 0, -g(3) * t103 - t117 * t106, t73, -t84, -g(1) * t128 - g(2) * t120 - t129, 0, 0, 0, 0, 0, 0, t66, t112, -t73, -g(1) * t119 - g(2) * t116 - g(3) * t118, 0, 0, 0, 0, 0, 0, t66, -t73, -t112, -g(1) * t113 - g(2) * t110 - g(3) * t115, 0, 0, 0, 0, 0, 0, t65, t114, t73, -g(1) * t109 - g(2) * t108 - g(3) * t111, 0, 0, 0, 0, 0, 0, t65, t73, -t114, -g(1) * (t70 * pkin(5) + t69 * qJ(6) + t109) - g(2) * (t68 * pkin(5) + t67 * qJ(6) + t108) - g(3) * (t75 * pkin(5) + t74 * qJ(6) + t111);];
U_reg  = t1;

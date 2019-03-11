% Calculate inertial parameters regressor of potential energy for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:22
% EndTime: 2019-03-09 02:48:22
% DurationCPUTime: 0.16s
% Computational Cost: add. (215->70), mult. (262->92), div. (0->0), fcn. (278->10), ass. (0->44)
t94 = pkin(9) + qJ(3);
t90 = sin(t94);
t97 = cos(pkin(10));
t120 = t90 * t97;
t95 = sin(pkin(10));
t121 = t90 * t95;
t126 = pkin(4) * t120 + qJ(5) * t121;
t125 = g(3) * pkin(6);
t91 = cos(t94);
t124 = pkin(3) * t91;
t123 = g(3) * t90;
t96 = sin(pkin(9));
t122 = t96 * pkin(2) + pkin(6);
t101 = sin(qJ(1));
t103 = cos(qJ(1));
t98 = cos(pkin(9));
t88 = t98 * pkin(2) + pkin(1);
t99 = -pkin(7) - qJ(2);
t119 = t101 * t88 + t103 * t99;
t118 = t101 * t90;
t117 = t101 * t95;
t116 = t101 * t97;
t115 = t103 * t90;
t114 = t103 * t95;
t113 = t103 * t97;
t112 = t90 * pkin(3) + t122;
t111 = -t101 * t99 + t103 * t88;
t110 = qJ(4) * t118 + t101 * t124 + t119;
t109 = g(1) * t103 + g(2) * t101;
t108 = -t91 * qJ(4) + t112;
t107 = qJ(4) * t115 + t103 * t124 + t111;
t72 = t91 * t117 + t113;
t73 = t91 * t116 - t114;
t106 = t73 * pkin(4) + t72 * qJ(5) + t110;
t74 = t91 * t114 - t116;
t105 = g(1) * t74 + g(2) * t72 + g(3) * t121;
t75 = t91 * t113 + t117;
t104 = t75 * pkin(4) + t74 * qJ(5) + t107;
t102 = cos(qJ(6));
t100 = sin(qJ(6));
t82 = g(1) * t101 - g(2) * t103;
t69 = -g(3) * t91 + t109 * t90;
t68 = -g(1) * t75 - g(2) * t73 - g(3) * t120;
t1 = [0, 0, 0, 0, 0, 0, -t109, t82, -g(3), -t125, 0, 0, 0, 0, 0, 0, -g(3) * t96 - t109 * t98, -g(3) * t98 + t109 * t96, -t82, -g(1) * (t103 * pkin(1) + t101 * qJ(2)) - g(2) * (t101 * pkin(1) - t103 * qJ(2)) - t125, 0, 0, 0, 0, 0, 0, -t109 * t91 - t123, t69, -t82, -g(1) * t111 - g(2) * t119 - g(3) * t122, 0, 0, 0, 0, 0, 0, t68, t105, -t69, -g(1) * t107 - g(2) * t110 - g(3) * t108, 0, 0, 0, 0, 0, 0, t68, -t69, -t105, -g(1) * t104 - g(2) * t106 - g(3) * (t108 + t126) 0, 0, 0, 0, 0, 0, -g(1) * (t74 * t100 + t75 * t102) - g(2) * (t72 * t100 + t73 * t102) - (t100 * t95 + t102 * t97) * t123, -g(1) * (-t75 * t100 + t74 * t102) - g(2) * (-t73 * t100 + t72 * t102) - (-t100 * t97 + t102 * t95) * t123, t69, -g(1) * (t75 * pkin(5) - pkin(8) * t115 + t104) - g(2) * (t73 * pkin(5) - pkin(8) * t118 + t106) - g(3) * (pkin(5) * t120 + (pkin(8) - qJ(4)) * t91 + t112 + t126);];
U_reg  = t1;

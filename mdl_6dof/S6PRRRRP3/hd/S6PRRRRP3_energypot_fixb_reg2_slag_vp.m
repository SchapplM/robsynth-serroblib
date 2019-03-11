% Calculate inertial parameters regressor of potential energy for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:11:06
% EndTime: 2019-03-09 00:11:06
% DurationCPUTime: 0.25s
% Computational Cost: add. (299->97), mult. (592->136), div. (0->0), fcn. (717->12), ass. (0->51)
t100 = sin(qJ(4));
t116 = pkin(4) * t100 + pkin(8);
t105 = -pkin(10) - pkin(9);
t98 = sin(pkin(6));
t129 = pkin(7) * t98;
t96 = qJ(4) + qJ(5);
t89 = sin(t96);
t128 = pkin(5) * t89 + t116;
t103 = cos(qJ(4));
t88 = t103 * pkin(4) + pkin(3);
t127 = cos(qJ(3));
t97 = sin(pkin(11));
t99 = cos(pkin(11));
t125 = t99 * pkin(1) + t97 * t129;
t101 = sin(qJ(3));
t124 = t101 * t98;
t102 = sin(qJ(2));
t123 = t102 * t98;
t104 = cos(qJ(2));
t122 = t104 * t98;
t120 = cos(pkin(6));
t121 = t120 * pkin(7) + qJ(1);
t115 = t102 * t120;
t78 = t99 * t104 - t97 * t115;
t119 = t78 * pkin(2) + t125;
t118 = pkin(2) * t123 + t121;
t117 = t98 * t127;
t114 = t104 * t120;
t113 = t97 * pkin(1) - t99 * t129;
t112 = g(1) * t97 - g(2) * t99;
t76 = t97 * t104 + t99 * t115;
t111 = t76 * pkin(2) + t113;
t77 = t99 * t102 + t97 * t114;
t110 = t77 * pkin(8) + t119;
t109 = -pkin(8) * t122 + t118;
t69 = t76 * t101 + t99 * t117;
t71 = t101 * t78 - t97 * t117;
t79 = t101 * t123 - t120 * t127;
t108 = g(1) * t71 + g(2) * t69 + g(3) * t79;
t75 = t102 * t97 - t99 * t114;
t107 = t75 * pkin(8) + t111;
t106 = -g(1) * t77 - g(2) * t75 + g(3) * t122;
t95 = -qJ(6) + t105;
t90 = cos(t96);
t81 = pkin(5) * t90 + t88;
t80 = t120 * t101 + t102 * t117;
t72 = t97 * t124 + t78 * t127;
t70 = -t99 * t124 + t76 * t127;
t67 = -g(1) * (t72 * t90 + t77 * t89) - g(2) * (t70 * t90 + t75 * t89) - g(3) * (-t89 * t122 + t80 * t90);
t66 = -g(1) * (-t72 * t89 + t77 * t90) - g(2) * (-t70 * t89 + t75 * t90) - g(3) * (-t90 * t122 - t80 * t89);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t99 - g(2) * t97, t112, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t76 - g(3) * t123, -t106, -g(3) * t120 - t112 * t98, -g(1) * t125 - g(2) * t113 - g(3) * t121, 0, 0, 0, 0, 0, 0, -g(1) * t72 - g(2) * t70 - g(3) * t80, t108, t106, -g(1) * t110 - g(2) * t107 - g(3) * t109, 0, 0, 0, 0, 0, 0, -g(1) * (t100 * t77 + t103 * t72) - g(2) * (t100 * t75 + t103 * t70) - g(3) * (-t100 * t122 + t103 * t80) -g(1) * (-t100 * t72 + t103 * t77) - g(2) * (-t100 * t70 + t103 * t75) - g(3) * (-t100 * t80 - t103 * t122) -t108, -g(1) * (pkin(3) * t72 + pkin(9) * t71 + t110) - g(2) * (pkin(3) * t70 + pkin(9) * t69 + t107) - g(3) * (pkin(3) * t80 + pkin(9) * t79 + t109) 0, 0, 0, 0, 0, 0, t67, t66, -t108, -g(1) * (-t105 * t71 + t116 * t77 + t72 * t88 + t119) - g(2) * (-t105 * t69 + t116 * t75 + t70 * t88 + t111) - g(3) * (-t105 * t79 - t116 * t122 + t80 * t88 + t118) 0, 0, 0, 0, 0, 0, t67, t66, -t108, -g(1) * (t128 * t77 - t71 * t95 + t72 * t81 + t119) - g(2) * (t128 * t75 - t69 * t95 + t70 * t81 + t111) - g(3) * (-t128 * t122 - t79 * t95 + t80 * t81 + t118);];
U_reg  = t1;

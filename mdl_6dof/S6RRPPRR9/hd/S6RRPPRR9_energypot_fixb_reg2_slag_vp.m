% Calculate inertial parameters regressor of potential energy for
% S6RRPPRR9
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:32:06
% EndTime: 2019-03-09 09:32:06
% DurationCPUTime: 0.22s
% Computational Cost: add. (204->81), mult. (447->109), div. (0->0), fcn. (516->10), ass. (0->49)
t132 = -pkin(3) - pkin(8);
t123 = cos(pkin(6));
t131 = t123 * pkin(8) + pkin(7);
t106 = cos(qJ(1));
t102 = sin(qJ(1));
t98 = sin(pkin(6));
t128 = t102 * t98;
t130 = t106 * pkin(1) + pkin(8) * t128;
t101 = sin(qJ(2));
t129 = t101 * t98;
t104 = cos(qJ(5));
t127 = t104 * t98;
t105 = cos(qJ(2));
t126 = t105 * t98;
t125 = t106 * t98;
t120 = t102 * t123;
t82 = t106 * t101 + t105 * t120;
t124 = t82 * qJ(3);
t122 = pkin(8) * t125;
t83 = -t101 * t120 + t106 * t105;
t121 = t83 * pkin(2) + t130;
t119 = t106 * t123;
t80 = t102 * t101 - t105 * t119;
t81 = t101 * t119 + t102 * t105;
t96 = t102 * pkin(1);
t118 = t81 * pkin(2) + t80 * qJ(3) + t96;
t117 = g(1) * t102 - g(2) * t106;
t116 = g(2) * (-pkin(4) + t132) * t125;
t115 = pkin(2) * t129 - qJ(3) * t126 + t131;
t114 = pkin(3) * t128 + t83 * qJ(4) + t121;
t100 = sin(qJ(5));
t70 = t100 * t128 - t83 * t104;
t72 = t100 * t125 + t81 * t104;
t78 = t123 * t100 - t101 * t127;
t113 = g(1) * t70 - g(2) * t72 + g(3) * t78;
t112 = t81 * qJ(4) + t118;
t111 = g(1) * t83 + g(2) * t81 + g(3) * t129;
t67 = -g(1) * t82 - g(2) * t80 + g(3) * t126;
t110 = t123 * pkin(3) + qJ(4) * t129 + t115;
t109 = -t80 * pkin(9) + t112;
t108 = t123 * pkin(4) + pkin(9) * t126 + t110;
t107 = pkin(4) * t128 + (-pkin(9) + qJ(3)) * t82 + t114;
t103 = cos(qJ(6));
t99 = sin(qJ(6));
t79 = t100 * t129 + t123 * t104;
t74 = -g(3) * t123 - t117 * t98;
t73 = t81 * t100 - t104 * t125;
t71 = t83 * t100 + t102 * t127;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t106 - g(2) * t102, t117, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -t111, -t67, t74, -g(1) * t130 - g(2) * (t96 - t122) - g(3) * t131, 0, 0, 0, 0, 0, 0, t74, t111, t67, -g(1) * (t121 + t124) - g(2) * (t118 - t122) - g(3) * t115, 0, 0, 0, 0, 0, 0, t74, t67, -t111, -g(1) * (t114 + t124) - g(2) * (t132 * t125 + t112) - g(3) * t110, 0, 0, 0, 0, 0, 0, -g(1) * t71 - g(2) * t73 - g(3) * t79, t113, -t67, -g(1) * t107 - g(2) * t109 - g(3) * t108 - t116, 0, 0, 0, 0, 0, 0, -g(1) * (t71 * t103 - t82 * t99) - g(2) * (t73 * t103 - t80 * t99) - g(3) * (t79 * t103 + t99 * t126) -g(1) * (-t82 * t103 - t71 * t99) - g(2) * (-t80 * t103 - t73 * t99) - g(3) * (t103 * t126 - t79 * t99) -t113, -g(1) * (t71 * pkin(5) + t70 * pkin(10) + t107) - g(2) * (t73 * pkin(5) - t72 * pkin(10) + t109) - g(3) * (t79 * pkin(5) + t78 * pkin(10) + t108) - t116;];
U_reg  = t1;

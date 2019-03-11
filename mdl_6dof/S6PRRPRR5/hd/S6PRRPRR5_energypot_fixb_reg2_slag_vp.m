% Calculate inertial parameters regressor of potential energy for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:20:46
% EndTime: 2019-03-08 22:20:46
% DurationCPUTime: 0.28s
% Computational Cost: add. (311->108), mult. (592->154), div. (0->0), fcn. (717->14), ass. (0->52)
t101 = sin(pkin(12));
t121 = t101 * pkin(4) + pkin(8);
t100 = pkin(12) + qJ(5);
t92 = sin(t100);
t133 = pkin(5) * t92 + t121;
t104 = cos(pkin(12));
t91 = t104 * pkin(4) + pkin(3);
t132 = cos(qJ(3));
t103 = sin(pkin(6));
t131 = pkin(7) * t103;
t106 = -pkin(9) - qJ(4);
t102 = sin(pkin(11));
t105 = cos(pkin(11));
t129 = t105 * pkin(1) + t102 * t131;
t127 = cos(pkin(6));
t128 = t127 * pkin(7) + qJ(1);
t107 = sin(qJ(3));
t126 = t103 * t107;
t108 = sin(qJ(2));
t125 = t103 * t108;
t109 = cos(qJ(2));
t124 = t103 * t109;
t119 = t108 * t127;
t79 = -t102 * t119 + t105 * t109;
t123 = t79 * pkin(2) + t129;
t122 = pkin(2) * t125 + t128;
t120 = t103 * t132;
t118 = t109 * t127;
t117 = t102 * pkin(1) - t105 * t131;
t116 = g(1) * t102 - g(2) * t105;
t78 = t102 * t118 + t105 * t108;
t115 = t78 * pkin(8) + t123;
t77 = t102 * t109 + t105 * t119;
t114 = t77 * pkin(2) + t117;
t113 = -pkin(8) * t124 + t122;
t70 = t105 * t120 + t77 * t107;
t72 = -t102 * t120 + t79 * t107;
t80 = t107 * t125 - t127 * t132;
t112 = g(1) * t72 + g(2) * t70 + g(3) * t80;
t76 = t102 * t108 - t105 * t118;
t111 = t76 * pkin(8) + t114;
t110 = -g(1) * t78 - g(2) * t76 + g(3) * t124;
t99 = -pkin(10) + t106;
t94 = qJ(6) + t100;
t93 = cos(t100);
t90 = cos(t94);
t89 = sin(t94);
t82 = pkin(5) * t93 + t91;
t81 = t127 * t107 + t108 * t120;
t73 = t102 * t126 + t79 * t132;
t71 = -t105 * t126 + t77 * t132;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t105 - g(2) * t102, t116, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t79 - g(2) * t77 - g(3) * t125, -t110, -g(3) * t127 - t116 * t103, -g(1) * t129 - g(2) * t117 - g(3) * t128, 0, 0, 0, 0, 0, 0, -g(1) * t73 - g(2) * t71 - g(3) * t81, t112, t110, -g(1) * t115 - g(2) * t111 - g(3) * t113, 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t101 + t73 * t104) - g(2) * (t76 * t101 + t71 * t104) - g(3) * (-t101 * t124 + t81 * t104) -g(1) * (-t73 * t101 + t78 * t104) - g(2) * (-t71 * t101 + t76 * t104) - g(3) * (-t81 * t101 - t104 * t124) -t112, -g(1) * (t73 * pkin(3) + t72 * qJ(4) + t115) - g(2) * (t71 * pkin(3) + t70 * qJ(4) + t111) - g(3) * (t81 * pkin(3) + t80 * qJ(4) + t113) 0, 0, 0, 0, 0, 0, -g(1) * (t73 * t93 + t78 * t92) - g(2) * (t71 * t93 + t76 * t92) - g(3) * (-t92 * t124 + t81 * t93) -g(1) * (-t73 * t92 + t78 * t93) - g(2) * (-t71 * t92 + t76 * t93) - g(3) * (-t93 * t124 - t81 * t92) -t112, -g(1) * (-t72 * t106 + t121 * t78 + t73 * t91 + t123) - g(2) * (-t70 * t106 + t121 * t76 + t71 * t91 + t114) - g(3) * (-t80 * t106 - t121 * t124 + t81 * t91 + t122) 0, 0, 0, 0, 0, 0, -g(1) * (t73 * t90 + t78 * t89) - g(2) * (t71 * t90 + t76 * t89) - g(3) * (-t89 * t124 + t81 * t90) -g(1) * (-t73 * t89 + t78 * t90) - g(2) * (-t71 * t89 + t76 * t90) - g(3) * (-t90 * t124 - t81 * t89) -t112, -g(1) * (t133 * t78 - t72 * t99 + t73 * t82 + t123) - g(2) * (t133 * t76 - t70 * t99 + t71 * t82 + t114) - g(3) * (-t133 * t124 - t80 * t99 + t81 * t82 + t122);];
U_reg  = t1;

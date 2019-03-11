% Calculate inertial parameters regressor of potential energy for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:51:31
% EndTime: 2019-03-09 00:51:31
% DurationCPUTime: 0.28s
% Computational Cost: add. (311->108), mult. (592->154), div. (0->0), fcn. (717->14), ass. (0->52)
t105 = sin(qJ(4));
t122 = t105 * pkin(4) + pkin(8);
t110 = -pkin(10) - pkin(9);
t101 = qJ(4) + qJ(5);
t93 = sin(t101);
t134 = pkin(5) * t93 + t122;
t108 = cos(qJ(4));
t92 = t108 * pkin(4) + pkin(3);
t133 = cos(qJ(3));
t103 = sin(pkin(6));
t132 = pkin(7) * t103;
t102 = sin(pkin(12));
t104 = cos(pkin(12));
t130 = t104 * pkin(1) + t102 * t132;
t128 = cos(pkin(6));
t129 = t128 * pkin(7) + qJ(1);
t106 = sin(qJ(3));
t127 = t103 * t106;
t107 = sin(qJ(2));
t126 = t103 * t107;
t109 = cos(qJ(2));
t125 = t103 * t109;
t120 = t107 * t128;
t80 = -t102 * t120 + t104 * t109;
t124 = t80 * pkin(2) + t130;
t123 = pkin(2) * t126 + t129;
t121 = t103 * t133;
t119 = t109 * t128;
t118 = t102 * pkin(1) - t104 * t132;
t117 = g(1) * t102 - g(2) * t104;
t79 = t102 * t119 + t104 * t107;
t116 = t79 * pkin(8) + t124;
t78 = t102 * t109 + t104 * t120;
t115 = t78 * pkin(2) + t118;
t114 = -pkin(8) * t125 + t123;
t71 = t104 * t121 + t78 * t106;
t73 = -t102 * t121 + t80 * t106;
t81 = t106 * t126 - t128 * t133;
t113 = g(1) * t73 + g(2) * t71 + g(3) * t81;
t77 = t102 * t107 - t104 * t119;
t112 = t77 * pkin(8) + t115;
t111 = -g(1) * t79 - g(2) * t77 + g(3) * t125;
t100 = -pkin(11) + t110;
t98 = qJ(6) + t101;
t94 = cos(t101);
t91 = cos(t98);
t90 = sin(t98);
t83 = pkin(5) * t94 + t92;
t82 = t128 * t106 + t107 * t121;
t74 = t102 * t127 + t80 * t133;
t72 = -t104 * t127 + t78 * t133;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t102, t117, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t80 - g(2) * t78 - g(3) * t126, -t111, -g(3) * t128 - t117 * t103, -g(1) * t130 - g(2) * t118 - g(3) * t129, 0, 0, 0, 0, 0, 0, -g(1) * t74 - g(2) * t72 - g(3) * t82, t113, t111, -g(1) * t116 - g(2) * t112 - g(3) * t114, 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t105 + t74 * t108) - g(2) * (t77 * t105 + t72 * t108) - g(3) * (-t105 * t125 + t82 * t108) -g(1) * (-t74 * t105 + t79 * t108) - g(2) * (-t72 * t105 + t77 * t108) - g(3) * (-t82 * t105 - t108 * t125) -t113, -g(1) * (t74 * pkin(3) + t73 * pkin(9) + t116) - g(2) * (t72 * pkin(3) + t71 * pkin(9) + t112) - g(3) * (t82 * pkin(3) + t81 * pkin(9) + t114) 0, 0, 0, 0, 0, 0, -g(1) * (t74 * t94 + t79 * t93) - g(2) * (t72 * t94 + t77 * t93) - g(3) * (-t93 * t125 + t82 * t94) -g(1) * (-t74 * t93 + t79 * t94) - g(2) * (-t72 * t93 + t77 * t94) - g(3) * (-t94 * t125 - t82 * t93) -t113, -g(1) * (-t73 * t110 + t122 * t79 + t74 * t92 + t124) - g(2) * (-t71 * t110 + t122 * t77 + t72 * t92 + t115) - g(3) * (-t81 * t110 - t122 * t125 + t82 * t92 + t123) 0, 0, 0, 0, 0, 0, -g(1) * (t74 * t91 + t79 * t90) - g(2) * (t72 * t91 + t77 * t90) - g(3) * (-t90 * t125 + t82 * t91) -g(1) * (-t74 * t90 + t79 * t91) - g(2) * (-t72 * t90 + t77 * t91) - g(3) * (-t91 * t125 - t82 * t90) -t113, -g(1) * (-t73 * t100 + t134 * t79 + t74 * t83 + t124) - g(2) * (-t71 * t100 + t134 * t77 + t72 * t83 + t115) - g(3) * (-t81 * t100 - t134 * t125 + t82 * t83 + t123);];
U_reg  = t1;

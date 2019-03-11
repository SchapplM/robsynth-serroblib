% Calculate inertial parameters regressor of potential energy for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:38:55
% EndTime: 2019-03-08 20:38:55
% DurationCPUTime: 0.28s
% Computational Cost: add. (331->106), mult. (533->153), div. (0->0), fcn. (633->14), ass. (0->52)
t104 = cos(pkin(11));
t101 = sin(pkin(11));
t102 = sin(pkin(6));
t126 = t101 * t102;
t130 = t104 * pkin(1) + pkin(7) * t126;
t107 = sin(qJ(5));
t108 = sin(qJ(2));
t105 = cos(pkin(6));
t110 = cos(qJ(2));
t120 = t105 * t110;
t77 = t101 * t108 - t104 * t120;
t129 = t77 * t107;
t79 = t101 * t120 + t104 * t108;
t128 = t79 * t107;
t127 = t105 * pkin(7) + qJ(1);
t125 = t102 * t104;
t124 = t102 * t108;
t123 = t102 * t110;
t100 = sin(pkin(12));
t122 = t105 * t100;
t121 = t105 * t108;
t119 = t100 * t126;
t118 = t107 * t123;
t95 = t101 * pkin(1);
t117 = -pkin(7) * t125 + t95;
t106 = -pkin(8) - qJ(3);
t80 = -t101 * t121 + t104 * t110;
t103 = cos(pkin(12));
t89 = t103 * pkin(3) + pkin(2);
t116 = pkin(3) * t119 - t79 * t106 + t80 * t89 + t130;
t115 = pkin(3) * t122 + t106 * t123 + t89 * t124 + t127;
t114 = g(1) * t101 - g(2) * t104;
t78 = t101 * t110 + t104 * t121;
t98 = pkin(12) + qJ(4);
t91 = sin(t98);
t92 = cos(t98);
t67 = t92 * t125 + t78 * t91;
t69 = -t92 * t126 + t80 * t91;
t73 = -t105 * t92 + t91 * t124;
t113 = g(1) * t69 + g(2) * t67 + g(3) * t73;
t66 = -g(1) * t79 - g(2) * t77 + g(3) * t123;
t112 = t78 * t89 - t77 * t106 + t95 + (-pkin(3) * t100 - pkin(7)) * t125;
t111 = -pkin(10) - pkin(9);
t109 = cos(qJ(5));
t99 = qJ(5) + qJ(6);
t94 = cos(t99);
t93 = sin(t99);
t90 = t109 * pkin(5) + pkin(4);
t74 = t105 * t91 + t92 * t124;
t70 = t91 * t126 + t80 * t92;
t68 = -t91 * t125 + t78 * t92;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t104 - g(2) * t101, t114, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t80 - g(2) * t78 - g(3) * t124, -t66, -g(3) * t105 - t114 * t102, -g(1) * t130 - g(2) * t117 - g(3) * t127, 0, 0, 0, 0, 0, 0, -g(1) * (t80 * t103 + t119) - g(2) * (-t100 * t125 + t78 * t103) - g(3) * (t103 * t124 + t122) -g(1) * (-t80 * t100 + t103 * t126) - g(2) * (-t78 * t100 - t103 * t125) - g(3) * (-t100 * t124 + t105 * t103) t66, -g(1) * (t80 * pkin(2) + t79 * qJ(3) + t130) - g(2) * (t78 * pkin(2) + t77 * qJ(3) + t117) - g(3) * ((pkin(2) * t108 - qJ(3) * t110) * t102 + t127) 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68 - g(3) * t74, t113, t66, -g(1) * t116 - g(2) * t112 - g(3) * t115, 0, 0, 0, 0, 0, 0, -g(1) * (t70 * t109 + t128) - g(2) * (t68 * t109 + t129) - g(3) * (t74 * t109 - t118) -g(1) * (-t70 * t107 + t79 * t109) - g(2) * (-t68 * t107 + t77 * t109) - g(3) * (-t74 * t107 - t109 * t123) -t113, -g(1) * (t70 * pkin(4) + t69 * pkin(9) + t116) - g(2) * (t68 * pkin(4) + t67 * pkin(9) + t112) - g(3) * (t74 * pkin(4) + t73 * pkin(9) + t115) 0, 0, 0, 0, 0, 0, -g(1) * (t70 * t94 + t79 * t93) - g(2) * (t68 * t94 + t77 * t93) - g(3) * (-t93 * t123 + t74 * t94) -g(1) * (-t70 * t93 + t79 * t94) - g(2) * (-t68 * t93 + t77 * t94) - g(3) * (-t94 * t123 - t74 * t93) -t113, -g(1) * (pkin(5) * t128 - t69 * t111 + t70 * t90 + t116) - g(2) * (pkin(5) * t129 - t67 * t111 + t68 * t90 + t112) - g(3) * (-pkin(5) * t118 - t73 * t111 + t74 * t90 + t115);];
U_reg  = t1;

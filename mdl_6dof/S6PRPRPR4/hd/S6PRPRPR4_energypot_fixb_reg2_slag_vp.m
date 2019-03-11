% Calculate inertial parameters regressor of potential energy for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_energypot_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:16
% EndTime: 2019-03-08 19:41:16
% DurationCPUTime: 0.28s
% Computational Cost: add. (331->106), mult. (533->153), div. (0->0), fcn. (633->14), ass. (0->52)
t106 = cos(pkin(10));
t102 = sin(pkin(10));
t103 = sin(pkin(6));
t126 = t102 * t103;
t130 = t106 * pkin(1) + pkin(7) * t126;
t100 = sin(pkin(12));
t110 = sin(qJ(2));
t107 = cos(pkin(6));
t111 = cos(qJ(2));
t120 = t107 * t111;
t77 = t102 * t110 - t106 * t120;
t129 = t77 * t100;
t79 = t102 * t120 + t106 * t110;
t128 = t79 * t100;
t127 = t107 * pkin(7) + qJ(1);
t125 = t103 * t106;
t124 = t103 * t110;
t123 = t103 * t111;
t101 = sin(pkin(11));
t122 = t107 * t101;
t121 = t107 * t110;
t119 = t101 * t126;
t118 = t100 * t123;
t95 = t102 * pkin(1);
t117 = -pkin(7) * t125 + t95;
t109 = -pkin(8) - qJ(3);
t80 = -t102 * t121 + t106 * t111;
t105 = cos(pkin(11));
t90 = t105 * pkin(3) + pkin(2);
t116 = pkin(3) * t119 - t79 * t109 + t80 * t90 + t130;
t115 = pkin(3) * t122 + t109 * t123 + t90 * t124 + t127;
t114 = g(1) * t102 - g(2) * t106;
t78 = t102 * t111 + t106 * t121;
t99 = pkin(11) + qJ(4);
t92 = sin(t99);
t94 = cos(t99);
t67 = t94 * t125 + t78 * t92;
t69 = -t94 * t126 + t80 * t92;
t73 = -t107 * t94 + t92 * t124;
t113 = g(1) * t69 + g(2) * t67 + g(3) * t73;
t66 = -g(1) * t79 - g(2) * t77 + g(3) * t123;
t112 = t78 * t90 - t77 * t109 + t95 + (-pkin(3) * t101 - pkin(7)) * t125;
t108 = -pkin(9) - qJ(5);
t104 = cos(pkin(12));
t98 = pkin(12) + qJ(6);
t93 = cos(t98);
t91 = sin(t98);
t89 = t104 * pkin(5) + pkin(4);
t74 = t107 * t92 + t94 * t124;
t70 = t92 * t126 + t80 * t94;
t68 = -t92 * t125 + t78 * t94;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t106 - g(2) * t102, t114, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t80 - g(2) * t78 - g(3) * t124, -t66, -g(3) * t107 - t114 * t103, -g(1) * t130 - g(2) * t117 - g(3) * t127, 0, 0, 0, 0, 0, 0, -g(1) * (t80 * t105 + t119) - g(2) * (-t101 * t125 + t78 * t105) - g(3) * (t105 * t124 + t122) -g(1) * (-t80 * t101 + t105 * t126) - g(2) * (-t78 * t101 - t105 * t125) - g(3) * (-t101 * t124 + t107 * t105) t66, -g(1) * (t80 * pkin(2) + t79 * qJ(3) + t130) - g(2) * (t78 * pkin(2) + t77 * qJ(3) + t117) - g(3) * ((pkin(2) * t110 - qJ(3) * t111) * t103 + t127) 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68 - g(3) * t74, t113, t66, -g(1) * t116 - g(2) * t112 - g(3) * t115, 0, 0, 0, 0, 0, 0, -g(1) * (t70 * t104 + t128) - g(2) * (t68 * t104 + t129) - g(3) * (t74 * t104 - t118) -g(1) * (-t70 * t100 + t79 * t104) - g(2) * (-t68 * t100 + t77 * t104) - g(3) * (-t74 * t100 - t104 * t123) -t113, -g(1) * (t70 * pkin(4) + t69 * qJ(5) + t116) - g(2) * (t68 * pkin(4) + t67 * qJ(5) + t112) - g(3) * (t74 * pkin(4) + t73 * qJ(5) + t115) 0, 0, 0, 0, 0, 0, -g(1) * (t70 * t93 + t79 * t91) - g(2) * (t68 * t93 + t77 * t91) - g(3) * (-t91 * t123 + t74 * t93) -g(1) * (-t70 * t91 + t79 * t93) - g(2) * (-t68 * t91 + t77 * t93) - g(3) * (-t93 * t123 - t74 * t91) -t113, -g(1) * (pkin(5) * t128 - t69 * t108 + t70 * t89 + t116) - g(2) * (pkin(5) * t129 - t67 * t108 + t68 * t89 + t112) - g(3) * (-pkin(5) * t118 - t73 * t108 + t74 * t89 + t115);];
U_reg  = t1;

% Calculate inertial parameters regressor of potential energy for
% S6PRRRRP1
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
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:16
% EndTime: 2019-03-08 23:59:17
% DurationCPUTime: 0.24s
% Computational Cost: add. (319->95), mult. (533->137), div. (0->0), fcn. (633->12), ass. (0->53)
t98 = sin(pkin(11));
t99 = sin(pkin(6));
t130 = t98 * t99;
t100 = cos(pkin(11));
t129 = pkin(1) * t100 + pkin(7) * t130;
t128 = t100 * t99;
t103 = sin(qJ(5));
t105 = sin(qJ(2));
t101 = cos(pkin(6));
t108 = cos(qJ(2));
t118 = t101 * t108;
t78 = -t100 * t118 + t105 * t98;
t127 = t103 * t78;
t80 = t100 * t105 + t118 * t98;
t126 = t103 * t80;
t104 = sin(qJ(3));
t125 = t104 * t99;
t124 = t105 * t99;
t107 = cos(qJ(3));
t123 = t107 * t99;
t122 = t108 * t99;
t121 = pkin(7) * t101 + qJ(1);
t120 = t101 * t104;
t119 = t101 * t105;
t117 = t98 * t125;
t116 = t103 * t122;
t94 = t98 * pkin(1);
t115 = -pkin(7) * t128 + t94;
t109 = -pkin(9) - pkin(8);
t81 = t100 * t108 - t119 * t98;
t91 = pkin(3) * t107 + pkin(2);
t114 = pkin(3) * t117 - t109 * t80 + t81 * t91 + t129;
t113 = pkin(3) * t120 + t109 * t122 + t124 * t91 + t121;
t112 = g(1) * t98 - g(2) * t100;
t79 = t100 * t119 + t108 * t98;
t97 = qJ(3) + qJ(4);
t92 = sin(t97);
t93 = cos(t97);
t68 = t128 * t93 + t79 * t92;
t70 = -t130 * t93 + t81 * t92;
t74 = -t101 * t93 + t124 * t92;
t111 = g(1) * t70 + g(2) * t68 + g(3) * t74;
t67 = -g(1) * t80 - g(2) * t78 + g(3) * t122;
t110 = t79 * t91 - t78 * t109 + t94 + (-pkin(3) * t104 - pkin(7)) * t128;
t106 = cos(qJ(5));
t102 = -qJ(6) - pkin(10);
t90 = pkin(5) * t106 + pkin(4);
t75 = t101 * t92 + t124 * t93;
t71 = t130 * t92 + t81 * t93;
t69 = -t128 * t92 + t79 * t93;
t65 = -g(1) * (t106 * t71 + t126) - g(2) * (t106 * t69 + t127) - g(3) * (t106 * t75 - t116);
t64 = -g(1) * (-t103 * t71 + t106 * t80) - g(2) * (-t103 * t69 + t106 * t78) - g(3) * (-t103 * t75 - t106 * t122);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t100 - g(2) * t98, t112, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t81 - g(2) * t79 - g(3) * t124, -t67, -g(3) * t101 - t112 * t99, -g(1) * t129 - g(2) * t115 - g(3) * t121, 0, 0, 0, 0, 0, 0, -g(1) * (t107 * t81 + t117) - g(2) * (-t100 * t125 + t107 * t79) - g(3) * (t105 * t123 + t120) -g(1) * (-t104 * t81 + t123 * t98) - g(2) * (-t100 * t123 - t104 * t79) - g(3) * (t101 * t107 - t104 * t124) t67, -g(1) * (pkin(2) * t81 + pkin(8) * t80 + t129) - g(2) * (pkin(2) * t79 + pkin(8) * t78 + t115) - g(3) * ((pkin(2) * t105 - pkin(8) * t108) * t99 + t121) 0, 0, 0, 0, 0, 0, -g(1) * t71 - g(2) * t69 - g(3) * t75, t111, t67, -g(1) * t114 - g(2) * t110 - g(3) * t113, 0, 0, 0, 0, 0, 0, t65, t64, -t111, -g(1) * (pkin(4) * t71 + pkin(10) * t70 + t114) - g(2) * (pkin(4) * t69 + pkin(10) * t68 + t110) - g(3) * (pkin(4) * t75 + pkin(10) * t74 + t113) 0, 0, 0, 0, 0, 0, t65, t64, -t111, -g(1) * (pkin(5) * t126 - t102 * t70 + t71 * t90 + t114) - g(2) * (pkin(5) * t127 - t102 * t68 + t69 * t90 + t110) - g(3) * (-pkin(5) * t116 - t102 * t74 + t75 * t90 + t113);];
U_reg  = t1;

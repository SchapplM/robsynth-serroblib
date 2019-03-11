% Calculate inertial parameters regressor of potential energy for
% S6PRRRRP4
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:17:46
% EndTime: 2019-03-09 00:17:46
% DurationCPUTime: 0.23s
% Computational Cost: add. (322->92), mult. (648->131), div. (0->0), fcn. (793->12), ass. (0->56)
t108 = sin(pkin(6));
t140 = pkin(7) * t108;
t111 = sin(qJ(4));
t107 = sin(pkin(11));
t109 = cos(pkin(11));
t113 = sin(qJ(2));
t110 = cos(pkin(6));
t116 = cos(qJ(2));
t131 = t110 * t116;
t89 = t107 * t113 - t109 * t131;
t139 = t111 * t89;
t91 = t107 * t131 + t109 * t113;
t138 = t111 * t91;
t137 = t109 * pkin(1) + t107 * t140;
t112 = sin(qJ(3));
t136 = t108 * t112;
t135 = t108 * t113;
t115 = cos(qJ(3));
t134 = t108 * t115;
t133 = t108 * t116;
t132 = t110 * t113;
t130 = t110 * pkin(7) + qJ(1);
t129 = pkin(2) * t135 + t130;
t128 = t107 * pkin(1) - t109 * t140;
t127 = g(1) * t107 - g(2) * t109;
t92 = -t107 * t132 + t109 * t116;
t126 = t92 * pkin(2) + t91 * pkin(8) + t137;
t125 = -pkin(8) * t133 + t129;
t106 = qJ(4) + qJ(5);
t101 = sin(t106);
t102 = cos(t106);
t90 = t107 * t116 + t109 * t132;
t78 = -t109 * t136 + t115 * t90;
t67 = t101 * t78 - t89 * t102;
t80 = t107 * t136 + t115 * t92;
t69 = t101 * t80 - t91 * t102;
t94 = t110 * t112 + t113 * t134;
t73 = t101 * t94 + t102 * t133;
t124 = g(1) * t69 + g(2) * t67 + g(3) * t73;
t77 = t109 * t134 + t112 * t90;
t79 = -t107 * t134 + t112 * t92;
t93 = -t110 * t115 + t112 * t135;
t123 = g(1) * t79 + g(2) * t77 + g(3) * t93;
t122 = t90 * pkin(2) + t89 * pkin(8) + t128;
t114 = cos(qJ(4));
t100 = pkin(4) * t114 + pkin(3);
t117 = -pkin(10) - pkin(9);
t121 = pkin(4) * t138 + t80 * t100 - t79 * t117 + t126;
t120 = -g(1) * t91 - g(2) * t89 + g(3) * t133;
t119 = pkin(4) * t139 + t78 * t100 - t77 * t117 + t122;
t118 = t94 * t100 - t93 * t117 + (-pkin(4) * t111 - pkin(8)) * t133 + t129;
t74 = -t101 * t133 + t102 * t94;
t70 = t101 * t91 + t102 * t80;
t68 = t101 * t89 + t102 * t78;
t65 = -g(1) * t70 - g(2) * t68 - g(3) * t74;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t109 - g(2) * t107, t127, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t92 - g(2) * t90 - g(3) * t135, -t120, -g(3) * t110 - t127 * t108, -g(1) * t137 - g(2) * t128 - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * t80 - g(2) * t78 - g(3) * t94, t123, t120, -g(1) * t126 - g(2) * t122 - g(3) * t125, 0, 0, 0, 0, 0, 0, -g(1) * (t114 * t80 + t138) - g(2) * (t114 * t78 + t139) - g(3) * (-t111 * t133 + t114 * t94) -g(1) * (-t111 * t80 + t114 * t91) - g(2) * (-t111 * t78 + t114 * t89) - g(3) * (-t111 * t94 - t114 * t133) -t123, -g(1) * (t80 * pkin(3) + t79 * pkin(9) + t126) - g(2) * (t78 * pkin(3) + t77 * pkin(9) + t122) - g(3) * (pkin(3) * t94 + pkin(9) * t93 + t125) 0, 0, 0, 0, 0, 0, t65, t124, -t123, -g(1) * t121 - g(2) * t119 - g(3) * t118, 0, 0, 0, 0, 0, 0, t65, -t123, -t124, -g(1) * (t70 * pkin(5) + t69 * qJ(6) + t121) - g(2) * (t68 * pkin(5) + t67 * qJ(6) + t119) - g(3) * (t74 * pkin(5) + t73 * qJ(6) + t118);];
U_reg  = t1;

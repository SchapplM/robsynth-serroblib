% Calculate inertial parameters regressor of potential energy for
% S6PRRRRP2
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:04:18
% EndTime: 2019-03-09 00:04:19
% DurationCPUTime: 0.23s
% Computational Cost: add. (346->91), mult. (579->133), div. (0->0), fcn. (697->12), ass. (0->57)
t108 = cos(pkin(11));
t106 = sin(pkin(11));
t107 = sin(pkin(6));
t137 = t106 * t107;
t138 = pkin(1) * t108 + pkin(7) * t137;
t136 = t107 * t108;
t111 = sin(qJ(3));
t135 = t107 * t111;
t112 = sin(qJ(2));
t134 = t107 * t112;
t114 = cos(qJ(3));
t133 = t107 * t114;
t115 = cos(qJ(2));
t132 = t107 * t115;
t109 = cos(pkin(6));
t131 = t109 * t111;
t130 = t109 * t112;
t129 = t109 * t115;
t128 = pkin(7) * t109 + qJ(1);
t127 = t106 * t135;
t102 = t106 * pkin(1);
t126 = -pkin(7) * t136 + t102;
t116 = -pkin(9) - pkin(8);
t89 = t106 * t129 + t108 * t112;
t90 = -t106 * t130 + t108 * t115;
t99 = pkin(3) * t114 + pkin(2);
t125 = pkin(3) * t127 - t116 * t89 + t90 * t99 + t138;
t124 = pkin(3) * t131 + t116 * t132 + t134 * t99 + t128;
t123 = g(1) * t106 - g(2) * t108;
t110 = sin(qJ(5));
t113 = cos(qJ(5));
t105 = qJ(3) + qJ(4);
t100 = sin(t105);
t101 = cos(t105);
t88 = t106 * t115 + t108 * t130;
t73 = -t100 * t136 + t101 * t88;
t87 = t106 * t112 - t108 * t129;
t65 = t110 * t73 - t113 * t87;
t75 = t100 * t137 + t101 * t90;
t67 = t110 * t75 - t113 * t89;
t82 = t100 * t109 + t101 * t134;
t76 = t110 * t82 + t113 * t132;
t122 = g(1) * t67 + g(2) * t65 + g(3) * t76;
t72 = t100 * t88 + t101 * t136;
t74 = t100 * t90 - t101 * t137;
t81 = t100 * t134 - t101 * t109;
t121 = g(1) * t74 + g(2) * t72 + g(3) * t81;
t120 = pkin(4) * t75 + pkin(10) * t74 + t125;
t69 = -g(1) * t89 - g(2) * t87 + g(3) * t132;
t119 = pkin(4) * t82 + pkin(10) * t81 + t124;
t118 = t102 + t88 * t99 - t87 * t116 + (-pkin(3) * t111 - pkin(7)) * t136;
t117 = pkin(4) * t73 + t72 * pkin(10) + t118;
t77 = -t110 * t132 + t113 * t82;
t68 = t110 * t89 + t113 * t75;
t66 = t110 * t87 + t113 * t73;
t63 = -g(1) * t68 - g(2) * t66 - g(3) * t77;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t108 - g(2) * t106, t123, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t90 - g(2) * t88 - g(3) * t134, -t69, -g(3) * t109 - t107 * t123, -g(1) * t138 - g(2) * t126 - g(3) * t128, 0, 0, 0, 0, 0, 0, -g(1) * (t114 * t90 + t127) - g(2) * (-t108 * t135 + t114 * t88) - g(3) * (t112 * t133 + t131) -g(1) * (t106 * t133 - t111 * t90) - g(2) * (-t108 * t133 - t111 * t88) - g(3) * (t109 * t114 - t111 * t134) t69, -g(1) * (pkin(2) * t90 + pkin(8) * t89 + t138) - g(2) * (pkin(2) * t88 + pkin(8) * t87 + t126) - g(3) * ((pkin(2) * t112 - pkin(8) * t115) * t107 + t128) 0, 0, 0, 0, 0, 0, -g(1) * t75 - g(2) * t73 - g(3) * t82, t121, t69, -g(1) * t125 - g(2) * t118 - g(3) * t124, 0, 0, 0, 0, 0, 0, t63, t122, -t121, -g(1) * t120 - g(2) * t117 - g(3) * t119, 0, 0, 0, 0, 0, 0, t63, -t121, -t122, -g(1) * (pkin(5) * t68 + qJ(6) * t67 + t120) - g(2) * (t66 * pkin(5) + t65 * qJ(6) + t117) - g(3) * (pkin(5) * t77 + qJ(6) * t76 + t119);];
U_reg  = t1;

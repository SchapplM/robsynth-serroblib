% Calculate inertial parameters regressor of potential energy for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:46:52
% EndTime: 2019-03-08 22:46:53
% DurationCPUTime: 0.23s
% Computational Cost: add. (322->92), mult. (648->131), div. (0->0), fcn. (793->12), ass. (0->56)
t106 = sin(pkin(6));
t138 = pkin(7) * t106;
t110 = sin(qJ(4));
t105 = sin(pkin(10));
t107 = cos(pkin(10));
t112 = sin(qJ(2));
t108 = cos(pkin(6));
t115 = cos(qJ(2));
t129 = t108 * t115;
t87 = t105 * t112 - t107 * t129;
t137 = t110 * t87;
t89 = t105 * t129 + t107 * t112;
t136 = t110 * t89;
t135 = t107 * pkin(1) + t105 * t138;
t111 = sin(qJ(3));
t134 = t106 * t111;
t133 = t106 * t112;
t114 = cos(qJ(3));
t132 = t106 * t114;
t131 = t106 * t115;
t130 = t108 * t112;
t128 = t108 * pkin(7) + qJ(1);
t127 = pkin(2) * t133 + t128;
t126 = t105 * pkin(1) - t107 * t138;
t125 = g(1) * t105 - g(2) * t107;
t90 = -t105 * t130 + t107 * t115;
t124 = t90 * pkin(2) + t89 * pkin(8) + t135;
t123 = -pkin(8) * t131 + t127;
t104 = qJ(4) + pkin(11);
t100 = cos(t104);
t88 = t105 * t115 + t107 * t130;
t76 = -t107 * t134 + t114 * t88;
t99 = sin(t104);
t65 = -t87 * t100 + t76 * t99;
t78 = t105 * t134 + t114 * t90;
t67 = -t89 * t100 + t78 * t99;
t92 = t108 * t111 + t112 * t132;
t71 = t100 * t131 + t92 * t99;
t122 = g(1) * t67 + g(2) * t65 + g(3) * t71;
t75 = t107 * t132 + t111 * t88;
t77 = -t105 * t132 + t111 * t90;
t91 = -t108 * t114 + t111 * t133;
t121 = g(1) * t77 + g(2) * t75 + g(3) * t91;
t120 = t88 * pkin(2) + t87 * pkin(8) + t126;
t109 = -qJ(5) - pkin(9);
t113 = cos(qJ(4));
t98 = pkin(4) * t113 + pkin(3);
t119 = pkin(4) * t136 - t77 * t109 + t78 * t98 + t124;
t118 = -g(1) * t89 - g(2) * t87 + g(3) * t131;
t117 = pkin(4) * t137 - t75 * t109 + t76 * t98 + t120;
t116 = t92 * t98 - t91 * t109 + (-pkin(4) * t110 - pkin(8)) * t131 + t127;
t72 = t92 * t100 - t99 * t131;
t68 = t100 * t78 + t89 * t99;
t66 = t100 * t76 + t87 * t99;
t63 = -g(1) * t68 - g(2) * t66 - g(3) * t72;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t107 - g(2) * t105, t125, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t90 - g(2) * t88 - g(3) * t133, -t118, -g(3) * t108 - t125 * t106, -g(1) * t135 - g(2) * t126 - g(3) * t128, 0, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t76 - g(3) * t92, t121, t118, -g(1) * t124 - g(2) * t120 - g(3) * t123, 0, 0, 0, 0, 0, 0, -g(1) * (t113 * t78 + t136) - g(2) * (t113 * t76 + t137) - g(3) * (-t110 * t131 + t92 * t113) -g(1) * (-t110 * t78 + t113 * t89) - g(2) * (-t110 * t76 + t113 * t87) - g(3) * (-t92 * t110 - t113 * t131) -t121, -g(1) * (t78 * pkin(3) + t77 * pkin(9) + t124) - g(2) * (t76 * pkin(3) + t75 * pkin(9) + t120) - g(3) * (t92 * pkin(3) + t91 * pkin(9) + t123) 0, 0, 0, 0, 0, 0, t63, t122, -t121, -g(1) * t119 - g(2) * t117 - g(3) * t116, 0, 0, 0, 0, 0, 0, t63, -t121, -t122, -g(1) * (t68 * pkin(5) + t67 * qJ(6) + t119) - g(2) * (t66 * pkin(5) + t65 * qJ(6) + t117) - g(3) * (t72 * pkin(5) + t71 * qJ(6) + t116);];
U_reg  = t1;

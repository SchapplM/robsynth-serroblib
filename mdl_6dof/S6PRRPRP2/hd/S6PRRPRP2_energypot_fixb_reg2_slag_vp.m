% Calculate inertial parameters regressor of potential energy for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:32:14
% EndTime: 2019-03-08 21:32:14
% DurationCPUTime: 0.24s
% Computational Cost: add. (346->91), mult. (579->133), div. (0->0), fcn. (697->12), ass. (0->57)
t109 = cos(pkin(10));
t107 = sin(pkin(10));
t108 = sin(pkin(6));
t138 = t107 * t108;
t139 = t109 * pkin(1) + pkin(7) * t138;
t137 = t108 * t109;
t113 = sin(qJ(3));
t136 = t108 * t113;
t114 = sin(qJ(2));
t135 = t108 * t114;
t116 = cos(qJ(3));
t134 = t108 * t116;
t117 = cos(qJ(2));
t133 = t108 * t117;
t110 = cos(pkin(6));
t132 = t110 * t113;
t131 = t110 * t114;
t130 = t110 * t117;
t129 = t110 * pkin(7) + qJ(1);
t128 = t107 * t136;
t103 = t107 * pkin(1);
t127 = -pkin(7) * t137 + t103;
t100 = t116 * pkin(3) + pkin(2);
t111 = -qJ(4) - pkin(8);
t90 = t107 * t130 + t109 * t114;
t91 = -t107 * t131 + t109 * t117;
t126 = pkin(3) * t128 + t91 * t100 - t90 * t111 + t139;
t125 = pkin(3) * t132 + t100 * t135 + t111 * t133 + t129;
t124 = g(1) * t107 - g(2) * t109;
t112 = sin(qJ(5));
t115 = cos(qJ(5));
t106 = qJ(3) + pkin(11);
t101 = sin(t106);
t102 = cos(t106);
t89 = t107 * t117 + t109 * t131;
t74 = -t101 * t137 + t89 * t102;
t88 = t107 * t114 - t109 * t130;
t66 = t74 * t112 - t88 * t115;
t76 = t101 * t138 + t91 * t102;
t68 = t76 * t112 - t90 * t115;
t83 = t110 * t101 + t102 * t135;
t77 = t83 * t112 + t115 * t133;
t123 = g(1) * t68 + g(2) * t66 + g(3) * t77;
t73 = t89 * t101 + t102 * t137;
t75 = t91 * t101 - t102 * t138;
t82 = t101 * t135 - t110 * t102;
t122 = g(1) * t75 + g(2) * t73 + g(3) * t82;
t121 = t76 * pkin(4) + t75 * pkin(9) + t126;
t70 = -g(1) * t90 - g(2) * t88 + g(3) * t133;
t120 = t83 * pkin(4) + t82 * pkin(9) + t125;
t119 = t103 + t89 * t100 - t88 * t111 + (-pkin(3) * t113 - pkin(7)) * t137;
t118 = t74 * pkin(4) + t73 * pkin(9) + t119;
t78 = -t112 * t133 + t83 * t115;
t69 = t90 * t112 + t76 * t115;
t67 = t88 * t112 + t74 * t115;
t64 = -g(1) * t69 - g(2) * t67 - g(3) * t78;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t109 - g(2) * t107, t124, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t91 - g(2) * t89 - g(3) * t135, -t70, -g(3) * t110 - t124 * t108, -g(1) * t139 - g(2) * t127 - g(3) * t129, 0, 0, 0, 0, 0, 0, -g(1) * (t91 * t116 + t128) - g(2) * (-t109 * t136 + t89 * t116) - g(3) * (t114 * t134 + t132) -g(1) * (t107 * t134 - t91 * t113) - g(2) * (-t109 * t134 - t89 * t113) - g(3) * (t110 * t116 - t113 * t135) t70, -g(1) * (t91 * pkin(2) + t90 * pkin(8) + t139) - g(2) * (t89 * pkin(2) + t88 * pkin(8) + t127) - g(3) * ((pkin(2) * t114 - pkin(8) * t117) * t108 + t129) 0, 0, 0, 0, 0, 0, -g(1) * t76 - g(2) * t74 - g(3) * t83, t122, t70, -g(1) * t126 - g(2) * t119 - g(3) * t125, 0, 0, 0, 0, 0, 0, t64, t123, -t122, -g(1) * t121 - g(2) * t118 - g(3) * t120, 0, 0, 0, 0, 0, 0, t64, -t122, -t123, -g(1) * (t69 * pkin(5) + t68 * qJ(6) + t121) - g(2) * (t67 * pkin(5) + t66 * qJ(6) + t118) - g(3) * (t78 * pkin(5) + t77 * qJ(6) + t120);];
U_reg  = t1;

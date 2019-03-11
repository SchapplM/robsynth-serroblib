% Calculate inertial parameters regressor of potential energy for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:36:45
% EndTime: 2019-03-08 23:36:45
% DurationCPUTime: 0.23s
% Computational Cost: add. (337->88), mult. (795->126), div. (0->0), fcn. (997->12), ass. (0->57)
t142 = pkin(9) - pkin(10);
t107 = cos(pkin(11));
t110 = sin(qJ(3));
t106 = sin(pkin(6));
t138 = cos(qJ(3));
t130 = t106 * t138;
t105 = sin(pkin(11));
t114 = cos(qJ(2));
t111 = sin(qJ(2));
t135 = cos(pkin(6));
t129 = t111 * t135;
t92 = t105 * t114 + t107 * t129;
t80 = t107 * t130 + t110 * t92;
t141 = t80 * pkin(9);
t94 = -t105 * t129 + t107 * t114;
t82 = -t105 * t130 + t110 * t94;
t140 = t82 * pkin(9);
t133 = t106 * t111;
t95 = t110 * t133 - t135 * t138;
t139 = t95 * pkin(9);
t137 = pkin(7) * t106;
t136 = pkin(1) * t107 + t105 * t137;
t134 = t106 * t110;
t132 = t106 * t114;
t131 = pkin(7) * t135 + qJ(1);
t128 = t114 * t135;
t127 = pkin(1) * t105 - t107 * t137;
t126 = g(1) * t105 - g(2) * t107;
t93 = t105 * t128 + t107 * t111;
t125 = pkin(2) * t94 + t93 * pkin(8) + t136;
t83 = t105 * t134 + t138 * t94;
t124 = pkin(3) * t83 + t125;
t123 = pkin(2) * t133 - pkin(8) * t132 + t131;
t109 = sin(qJ(4));
t113 = cos(qJ(4));
t81 = -t107 * t134 + t138 * t92;
t91 = t105 * t111 - t107 * t128;
t73 = t109 * t81 - t113 * t91;
t75 = t109 * t83 - t113 * t93;
t96 = t110 * t135 + t111 * t130;
t84 = t109 * t96 + t113 * t132;
t122 = g(1) * t75 + g(2) * t73 + g(3) * t84;
t70 = g(1) * t82 + g(2) * t80 + g(3) * t95;
t121 = pkin(3) * t96 + t123;
t120 = pkin(2) * t92 + t91 * pkin(8) + t127;
t119 = -g(1) * t93 - g(2) * t91 + g(3) * t132;
t118 = pkin(3) * t81 + t120;
t76 = t109 * t93 + t113 * t83;
t117 = pkin(4) * t76 + t75 * qJ(5) + t124;
t85 = -t109 * t132 + t113 * t96;
t116 = pkin(4) * t85 + t84 * qJ(5) + t121;
t74 = t109 * t91 + t113 * t81;
t115 = pkin(4) * t74 + t73 * qJ(5) + t118;
t112 = cos(qJ(6));
t108 = sin(qJ(6));
t68 = -g(1) * t76 - g(2) * t74 - g(3) * t85;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t107 - g(2) * t105, t126, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t94 - g(2) * t92 - g(3) * t133, -t119, -g(3) * t135 - t106 * t126, -g(1) * t136 - g(2) * t127 - g(3) * t131, 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t81 - g(3) * t96, t70, t119, -g(1) * t125 - g(2) * t120 - g(3) * t123, 0, 0, 0, 0, 0, 0, t68, t122, -t70, -g(1) * (t124 + t140) - g(2) * (t118 + t141) - g(3) * (t121 + t139) 0, 0, 0, 0, 0, 0, t68, -t70, -t122, -g(1) * (t117 + t140) - g(2) * (t115 + t141) - g(3) * (t116 + t139) 0, 0, 0, 0, 0, 0, -g(1) * (t108 * t75 + t112 * t76) - g(2) * (t108 * t73 + t112 * t74) - g(3) * (t108 * t84 + t112 * t85) -g(1) * (-t108 * t76 + t112 * t75) - g(2) * (-t108 * t74 + t112 * t73) - g(3) * (-t108 * t85 + t112 * t84) t70, -g(1) * (t76 * pkin(5) + t142 * t82 + t117) - g(2) * (t74 * pkin(5) + t142 * t80 + t115) - g(3) * (t85 * pkin(5) + t142 * t95 + t116);];
U_reg  = t1;

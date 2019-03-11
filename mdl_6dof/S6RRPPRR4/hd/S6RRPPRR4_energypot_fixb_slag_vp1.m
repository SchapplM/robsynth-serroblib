% Calculate potential energy for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:05
% EndTime: 2019-03-09 09:01:05
% DurationCPUTime: 0.51s
% Computational Cost: add. (336->120), mult. (724->162), div. (0->0), fcn. (881->12), ass. (0->52)
t138 = rSges(6,3) + pkin(9);
t137 = pkin(10) + rSges(7,3);
t107 = cos(pkin(6));
t136 = t107 * pkin(8) + pkin(7);
t111 = sin(qJ(1));
t115 = cos(qJ(1));
t104 = sin(pkin(11));
t106 = cos(pkin(11));
t110 = sin(qJ(2));
t114 = cos(qJ(2));
t121 = t104 * t114 + t106 * t110;
t120 = t121 * t107;
t129 = t106 * t114;
t94 = -t104 * t110 + t129;
t82 = -t111 * t120 + t115 * t94;
t100 = pkin(2) * t114 + pkin(1);
t97 = t115 * t100;
t135 = t82 * pkin(3) + t97;
t105 = sin(pkin(6));
t92 = pkin(2) * t107 * t110 + (-pkin(8) - qJ(3)) * t105;
t134 = t111 * t100 + t115 * t92;
t133 = rSges(5,3) + qJ(4);
t132 = t105 * t110;
t131 = t105 * t111;
t130 = t105 * t115;
t128 = t110 * t115;
t127 = t111 * t110;
t126 = t111 * t114;
t125 = t114 * t115;
t80 = t111 * t94 + t115 * t120;
t124 = t80 * pkin(3) + t134;
t123 = pkin(2) * t132 + t107 * qJ(3) + t136;
t91 = t121 * t105;
t122 = t91 * pkin(3) + t123;
t119 = t94 * t107;
t81 = -t111 * t119 - t115 * t121;
t118 = pkin(4) * t131 - t81 * qJ(4) - t111 * t92 + t135;
t90 = t104 * t132 - t105 * t129;
t117 = t107 * pkin(4) + t90 * qJ(4) + t122;
t79 = -t111 * t121 + t115 * t119;
t116 = -pkin(4) * t130 - t79 * qJ(4) + t124;
t113 = cos(qJ(5));
t112 = cos(qJ(6));
t109 = sin(qJ(5));
t108 = sin(qJ(6));
t84 = t107 * t113 + t109 * t90;
t83 = t107 * t109 - t90 * t113;
t75 = -t79 * t109 - t113 * t130;
t74 = t109 * t130 - t79 * t113;
t73 = -t109 * t81 + t113 * t131;
t72 = t109 * t131 + t81 * t113;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t115 - t111 * rSges(2,2)) + g(2) * (t111 * rSges(2,1) + rSges(2,2) * t115) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t115 * pkin(1) + (-t107 * t127 + t125) * rSges(3,1) + (-t107 * t126 - t128) * rSges(3,2)) + g(2) * (t111 * pkin(1) + (t107 * t128 + t126) * rSges(3,1) + (t107 * t125 - t127) * rSges(3,2)) + g(3) * (rSges(3,3) * t107 + t136) + (g(3) * (rSges(3,1) * t110 + rSges(3,2) * t114) + (g(1) * t111 - g(2) * t115) * (rSges(3,3) + pkin(8))) * t105) - m(4) * (g(1) * (rSges(4,1) * t82 + rSges(4,2) * t81 + t97 + (rSges(4,3) * t105 - t92) * t111) + g(2) * (t80 * rSges(4,1) + t79 * rSges(4,2) - rSges(4,3) * t130 + t134) + g(3) * (rSges(4,1) * t91 - rSges(4,2) * t90 + rSges(4,3) * t107 + t123)) - m(5) * (g(1) * (-rSges(5,2) * t82 - t133 * t81 + (rSges(5,1) * t105 - t92) * t111 + t135) + g(2) * (-rSges(5,1) * t130 - t80 * rSges(5,2) - t133 * t79 + t124) + g(3) * (rSges(5,1) * t107 - rSges(5,2) * t91 + t133 * t90 + t122)) - m(6) * (g(1) * (rSges(6,1) * t73 - rSges(6,2) * t72 + t138 * t82 + t118) + g(2) * (t75 * rSges(6,1) + t74 * rSges(6,2) + t138 * t80 + t116) + g(3) * (rSges(6,1) * t84 - rSges(6,2) * t83 + t138 * t91 + t117)) - m(7) * (g(1) * (t73 * pkin(5) + t82 * pkin(9) + (t108 * t82 + t112 * t73) * rSges(7,1) + (-t108 * t73 + t112 * t82) * rSges(7,2) + t137 * t72 + t118) + g(2) * (t75 * pkin(5) + t80 * pkin(9) + (t108 * t80 + t112 * t75) * rSges(7,1) + (-t108 * t75 + t112 * t80) * rSges(7,2) - t137 * t74 + t116) + g(3) * (t84 * pkin(5) + t91 * pkin(9) + (t108 * t91 + t112 * t84) * rSges(7,1) + (-t108 * t84 + t112 * t91) * rSges(7,2) + t137 * t83 + t117));
U  = t1;

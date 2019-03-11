% Calculate potential energy for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:01
% EndTime: 2019-03-09 10:38:02
% DurationCPUTime: 0.49s
% Computational Cost: add. (372->109), mult. (818->142), div. (0->0), fcn. (1011->12), ass. (0->54)
t109 = sin(pkin(11));
t114 = sin(qJ(2));
t118 = cos(qJ(2));
t138 = cos(pkin(11));
t98 = -t114 * t109 + t118 * t138;
t146 = rSges(6,1) + pkin(9);
t144 = rSges(5,3) + pkin(9);
t143 = pkin(10) + rSges(7,3);
t142 = pkin(2) * t114;
t111 = cos(pkin(6));
t141 = t111 * pkin(8) + pkin(7);
t140 = rSges(6,3) + qJ(5);
t106 = pkin(2) * t118 + pkin(1);
t115 = sin(qJ(1));
t119 = cos(qJ(1));
t110 = sin(pkin(6));
t96 = t111 * t142 + (-pkin(8) - qJ(3)) * t110;
t139 = t115 * t106 + t119 * t96;
t137 = t110 * t115;
t136 = t110 * t119;
t134 = t114 * t119;
t133 = t115 * t114;
t132 = t115 * t118;
t131 = t118 * t119;
t97 = -t118 * t109 - t114 * t138;
t95 = t97 * t111;
t84 = t115 * t98 - t119 * t95;
t130 = t84 * pkin(3) + t139;
t113 = sin(qJ(4));
t117 = cos(qJ(4));
t78 = -t113 * t136 + t117 * t84;
t128 = t78 * pkin(4) + t130;
t127 = t111 * qJ(3) + t110 * t142 + t141;
t101 = t119 * t106;
t86 = t115 * t95 + t119 * t98;
t126 = t86 * pkin(3) - t115 * t96 + t101;
t94 = t97 * t110;
t125 = -t94 * pkin(3) + t127;
t80 = t113 * t137 + t117 * t86;
t124 = t80 * pkin(4) + t126;
t89 = t111 * t113 - t117 * t94;
t123 = t89 * pkin(4) + t125;
t112 = sin(qJ(6));
t116 = cos(qJ(6));
t122 = t112 * rSges(7,1) + t116 * rSges(7,2) + qJ(5);
t121 = t116 * rSges(7,1) - t112 * rSges(7,2) + pkin(5) + pkin(9);
t120 = t111 * t98;
t93 = t98 * t110;
t88 = -t111 * t117 - t113 * t94;
t85 = -t115 * t120 + t119 * t97;
t83 = t115 * t97 + t119 * t120;
t79 = t113 * t86 - t117 * t137;
t77 = t84 * t113 + t117 * t136;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t119 - t115 * rSges(2,2)) + g(2) * (t115 * rSges(2,1) + rSges(2,2) * t119) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t119 * pkin(1) + (-t111 * t133 + t131) * rSges(3,1) + (-t111 * t132 - t134) * rSges(3,2)) + g(2) * (t115 * pkin(1) + (t111 * t134 + t132) * rSges(3,1) + (t111 * t131 - t133) * rSges(3,2)) + g(3) * (rSges(3,3) * t111 + t141) + (g(3) * (rSges(3,1) * t114 + rSges(3,2) * t118) + (g(1) * t115 - g(2) * t119) * (rSges(3,3) + pkin(8))) * t110) - m(4) * (g(1) * (rSges(4,1) * t86 + rSges(4,2) * t85 + t101 + (rSges(4,3) * t110 - t96) * t115) + g(2) * (t84 * rSges(4,1) + t83 * rSges(4,2) - rSges(4,3) * t136 + t139) + g(3) * (-rSges(4,1) * t94 + rSges(4,2) * t93 + rSges(4,3) * t111 + t127)) - m(5) * (g(1) * (rSges(5,1) * t80 - rSges(5,2) * t79 - t144 * t85 + t126) + g(2) * (rSges(5,1) * t78 - rSges(5,2) * t77 - t144 * t83 + t130) + g(3) * (rSges(5,1) * t89 - rSges(5,2) * t88 - t144 * t93 + t125)) - m(6) * (g(1) * (-rSges(6,2) * t80 + t140 * t79 - t146 * t85 + t124) + g(2) * (-rSges(6,2) * t78 + t140 * t77 - t146 * t83 + t128) + g(3) * (-rSges(6,2) * t89 + t140 * t88 - t146 * t93 + t123)) - m(7) * (g(1) * (-t121 * t85 + t122 * t79 + t143 * t80 + t124) + g(2) * (-t121 * t83 + t122 * t77 + t143 * t78 + t128) + g(3) * (-t121 * t93 + t122 * t88 + t143 * t89 + t123));
U  = t1;

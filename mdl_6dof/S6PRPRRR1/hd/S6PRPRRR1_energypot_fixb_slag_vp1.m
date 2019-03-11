% Calculate potential energy for
% S6PRPRRR1
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:14
% EndTime: 2019-03-08 20:22:15
% DurationCPUTime: 0.55s
% Computational Cost: add. (383->127), mult. (724->177), div. (0->0), fcn. (880->14), ass. (0->57)
t149 = pkin(8) + rSges(5,3);
t148 = pkin(10) + rSges(7,3);
t126 = cos(qJ(2));
t109 = pkin(2) * t126 + pkin(1);
t116 = sin(pkin(11));
t119 = cos(pkin(11));
t117 = sin(pkin(6));
t120 = cos(pkin(6));
t123 = sin(qJ(2));
t139 = t120 * t123;
t97 = pkin(2) * t139 + (-pkin(7) - qJ(3)) * t117;
t147 = t116 * t109 + t119 * t97;
t146 = t116 * t117;
t145 = t117 * t119;
t122 = sin(qJ(4));
t144 = t117 * t122;
t143 = t117 * t123;
t125 = cos(qJ(4));
t142 = t117 * t125;
t118 = cos(pkin(12));
t141 = t118 * t126;
t140 = t120 * t122;
t138 = t120 * t126;
t137 = t120 * pkin(7) + qJ(1);
t136 = t116 * t144;
t135 = t119 * t144;
t103 = t119 * t109;
t134 = -t116 * t97 + t103;
t133 = pkin(2) * t143 + t120 * qJ(3) + t137;
t115 = sin(pkin(12));
t132 = t115 * t126 + t118 * t123;
t99 = -t115 * t123 + t141;
t108 = pkin(4) * t125 + pkin(3);
t127 = -pkin(9) - pkin(8);
t130 = t99 * t120;
t86 = -t116 * t130 - t119 * t132;
t96 = t132 * t120;
t87 = -t116 * t96 + t119 * t99;
t131 = pkin(4) * t136 + t87 * t108 + t86 * t127 + t134;
t94 = t115 * t143 - t117 * t141;
t95 = t132 * t117;
t129 = pkin(4) * t140 + t95 * t108 - t94 * t127 + t133;
t84 = -t116 * t132 + t119 * t130;
t85 = t116 * t99 + t119 * t96;
t128 = -pkin(4) * t135 + t85 * t108 + t84 * t127 + t147;
t124 = cos(qJ(6));
t121 = sin(qJ(6));
t114 = qJ(4) + qJ(5);
t112 = cos(t114);
t111 = sin(t114);
t89 = t111 * t120 + t112 * t95;
t88 = t111 * t95 - t120 * t112;
t79 = t111 * t146 + t112 * t87;
t78 = t111 * t87 - t112 * t146;
t77 = -t111 * t145 + t112 * t85;
t76 = t111 * t85 + t112 * t145;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t119 - rSges(2,2) * t116) + g(2) * (rSges(2,1) * t116 + rSges(2,2) * t119) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t119 * pkin(1) + (-t116 * t139 + t119 * t126) * rSges(3,1) + (-t116 * t138 - t119 * t123) * rSges(3,2)) + g(2) * (t116 * pkin(1) + (t116 * t126 + t119 * t139) * rSges(3,1) + (-t116 * t123 + t119 * t138) * rSges(3,2)) + g(3) * (rSges(3,3) * t120 + t137) + (g(3) * (rSges(3,1) * t123 + rSges(3,2) * t126) + (g(1) * t116 - g(2) * t119) * (rSges(3,3) + pkin(7))) * t117) - m(4) * (g(1) * (rSges(4,1) * t87 + rSges(4,2) * t86 + t103 + (rSges(4,3) * t117 - t97) * t116) + g(2) * (rSges(4,1) * t85 + rSges(4,2) * t84 - rSges(4,3) * t145 + t147) + g(3) * (rSges(4,1) * t95 - rSges(4,2) * t94 + rSges(4,3) * t120 + t133)) - m(5) * (g(1) * (t87 * pkin(3) + (t125 * t87 + t136) * rSges(5,1) + (t116 * t142 - t122 * t87) * rSges(5,2) - t149 * t86 + t134) + g(2) * (t85 * pkin(3) + (t125 * t85 - t135) * rSges(5,1) + (-t119 * t142 - t122 * t85) * rSges(5,2) - t149 * t84 + t147) + g(3) * (t95 * pkin(3) + (t125 * t95 + t140) * rSges(5,1) + (t120 * t125 - t122 * t95) * rSges(5,2) + t149 * t94 + t133)) - m(6) * (g(1) * (rSges(6,1) * t79 - rSges(6,2) * t78 - t86 * rSges(6,3) + t131) + g(2) * (rSges(6,1) * t77 - rSges(6,2) * t76 - rSges(6,3) * t84 + t128) + g(3) * (rSges(6,1) * t89 - rSges(6,2) * t88 + rSges(6,3) * t94 + t129)) - m(7) * (g(1) * (t79 * pkin(5) + (-t121 * t86 + t124 * t79) * rSges(7,1) + (-t121 * t79 - t124 * t86) * rSges(7,2) + t148 * t78 + t131) + g(2) * (t77 * pkin(5) + (-t121 * t84 + t124 * t77) * rSges(7,1) + (-t121 * t77 - t124 * t84) * rSges(7,2) + t148 * t76 + t128) + g(3) * (t89 * pkin(5) + (t121 * t94 + t124 * t89) * rSges(7,1) + (-t121 * t89 + t124 * t94) * rSges(7,2) + t148 * t88 + t129));
U  = t1;

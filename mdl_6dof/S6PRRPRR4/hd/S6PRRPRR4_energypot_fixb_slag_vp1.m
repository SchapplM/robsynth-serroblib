% Calculate potential energy for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:10
% EndTime: 2019-03-08 22:10:10
% DurationCPUTime: 0.41s
% Computational Cost: add. (334->116), mult. (737->148), div. (0->0), fcn. (900->12), ass. (0->56)
t150 = rSges(5,2) + pkin(8);
t149 = rSges(4,3) + pkin(8);
t148 = rSges(6,3) - pkin(8);
t147 = pkin(10) + rSges(7,3);
t146 = cos(qJ(3));
t114 = sin(pkin(6));
t145 = pkin(7) * t114;
t144 = rSges(5,3) + qJ(4);
t143 = cos(pkin(6));
t118 = sin(qJ(3));
t142 = t114 * t118;
t119 = sin(qJ(2));
t141 = t114 * t119;
t122 = cos(qJ(2));
t140 = t114 * t122;
t139 = pkin(7) * t143 + qJ(1);
t113 = sin(pkin(11));
t115 = cos(pkin(11));
t138 = pkin(1) * t115 + t113 * t145;
t137 = -pkin(9) - t148;
t133 = t119 * t143;
t101 = -t113 * t133 + t115 * t122;
t136 = pkin(2) * t101 + t138;
t135 = pkin(2) * t141 + t139;
t134 = t114 * t146;
t132 = t122 * t143;
t92 = t101 * t146 + t113 * t142;
t131 = pkin(3) * t92 + t136;
t103 = t118 * t143 + t119 * t134;
t130 = pkin(3) * t103 + t135;
t110 = t113 * pkin(1);
t99 = t113 * t122 + t115 * t133;
t129 = pkin(2) * t99 - t115 * t145 + t110;
t116 = sin(qJ(6));
t120 = cos(qJ(6));
t128 = rSges(7,1) * t120 - rSges(7,2) * t116 + pkin(5);
t90 = -t115 * t142 + t146 * t99;
t127 = pkin(3) * t90 + t129;
t126 = -rSges(7,1) * t116 - rSges(7,2) * t120 + pkin(8) - pkin(9);
t91 = t101 * t118 - t113 * t134;
t125 = pkin(4) * t92 + t91 * qJ(4) + t131;
t102 = t118 * t141 - t143 * t146;
t124 = pkin(4) * t103 + pkin(9) * t140 + t102 * qJ(4) + t130;
t89 = t115 * t134 + t118 * t99;
t123 = pkin(4) * t90 + t89 * qJ(4) + t127;
t121 = cos(qJ(5));
t117 = sin(qJ(5));
t100 = t113 * t132 + t115 * t119;
t98 = t113 * t119 - t115 * t132;
t82 = t102 * t117 + t103 * t121;
t81 = -t102 * t121 + t103 * t117;
t80 = t117 * t91 + t121 * t92;
t79 = t117 * t92 - t121 * t91;
t78 = t117 * t89 + t121 * t90;
t77 = t117 * t90 - t121 * t89;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t115 - rSges(2,2) * t113) + g(2) * (rSges(2,1) * t113 + rSges(2,2) * t115) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t101 - rSges(3,2) * t100 + t138) + g(2) * (t99 * rSges(3,1) - t98 * rSges(3,2) + t110) + g(3) * (t143 * rSges(3,3) + t139) + (g(1) * rSges(3,3) * t113 + g(3) * (rSges(3,1) * t119 + rSges(3,2) * t122) + g(2) * (-rSges(3,3) - pkin(7)) * t115) * t114) - m(4) * (g(1) * (t92 * rSges(4,1) - t91 * rSges(4,2) + t100 * t149 + t136) + g(2) * (t90 * rSges(4,1) - t89 * rSges(4,2) + t149 * t98 + t129) + g(3) * (t103 * rSges(4,1) - t102 * rSges(4,2) - t140 * t149 + t135)) - m(5) * (g(1) * (t92 * rSges(5,1) + t100 * t150 + t144 * t91 + t131) + g(2) * (t90 * rSges(5,1) + t144 * t89 + t150 * t98 + t127) + g(3) * (t103 * rSges(5,1) + t102 * t144 - t140 * t150 + t130)) - m(6) * (g(3) * (t82 * rSges(6,1) - t81 * rSges(6,2) + t140 * t148 + t124) + (t78 * rSges(6,1) - t77 * rSges(6,2) + t137 * t98 + t123) * g(2) + (t80 * rSges(6,1) - t79 * rSges(6,2) + t100 * t137 + t125) * g(1)) - m(7) * (g(1) * (t100 * t126 + t128 * t80 + t147 * t79 + t125) + g(2) * (t126 * t98 + t128 * t78 + t147 * t77 + t123) + g(3) * (t82 * pkin(5) - pkin(8) * t140 + (t116 * t140 + t120 * t82) * rSges(7,1) + (-t116 * t82 + t120 * t140) * rSges(7,2) + t147 * t81 + t124));
U  = t1;

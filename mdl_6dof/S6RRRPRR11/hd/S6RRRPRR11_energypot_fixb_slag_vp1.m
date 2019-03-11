% Calculate potential energy for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR11_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:24:43
% EndTime: 2019-03-09 19:24:44
% DurationCPUTime: 0.41s
% Computational Cost: add. (334->116), mult. (737->148), div. (0->0), fcn. (900->12), ass. (0->56)
t150 = rSges(5,2) + pkin(9);
t149 = rSges(4,3) + pkin(9);
t148 = rSges(6,3) - pkin(9);
t147 = pkin(11) + rSges(7,3);
t146 = cos(qJ(3));
t143 = cos(pkin(6));
t145 = t143 * pkin(8) + pkin(7);
t144 = rSges(5,3) + qJ(4);
t113 = sin(pkin(6));
t117 = sin(qJ(2));
t142 = t113 * t117;
t118 = sin(qJ(1));
t141 = t113 * t118;
t121 = cos(qJ(2));
t140 = t113 * t121;
t122 = cos(qJ(1));
t139 = t113 * t122;
t138 = t122 * pkin(1) + pkin(8) * t141;
t137 = -pkin(10) - t148;
t136 = pkin(2) * t142 + t145;
t133 = t118 * t143;
t103 = -t117 * t133 + t122 * t121;
t135 = t103 * pkin(2) + t138;
t134 = t113 * t146;
t132 = t122 * t143;
t116 = sin(qJ(3));
t99 = t143 * t116 + t117 * t134;
t131 = t99 * pkin(3) + t136;
t92 = t103 * t146 + t116 * t141;
t130 = t92 * pkin(3) + t135;
t101 = t117 * t132 + t118 * t121;
t111 = t118 * pkin(1);
t129 = t101 * pkin(2) - pkin(8) * t139 + t111;
t114 = sin(qJ(6));
t119 = cos(qJ(6));
t128 = rSges(7,1) * t119 - rSges(7,2) * t114 + pkin(5);
t90 = t101 * t146 - t116 * t139;
t127 = t90 * pkin(3) + t129;
t126 = -rSges(7,1) * t114 - rSges(7,2) * t119 + pkin(9) - pkin(10);
t91 = t103 * t116 - t118 * t134;
t125 = t92 * pkin(4) + t91 * qJ(4) + t130;
t98 = t116 * t142 - t143 * t146;
t124 = t99 * pkin(4) + pkin(10) * t140 + t98 * qJ(4) + t131;
t89 = t101 * t116 + t122 * t134;
t123 = t90 * pkin(4) + t89 * qJ(4) + t127;
t120 = cos(qJ(5));
t115 = sin(qJ(5));
t102 = t122 * t117 + t121 * t133;
t100 = t118 * t117 - t121 * t132;
t82 = t115 * t98 + t120 * t99;
t81 = t115 * t99 - t98 * t120;
t80 = t115 * t91 + t120 * t92;
t79 = t115 * t92 - t91 * t120;
t78 = t115 * t89 + t120 * t90;
t77 = t115 * t90 - t89 * t120;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t122 - t118 * rSges(2,2)) + g(2) * (t118 * rSges(2,1) + rSges(2,2) * t122) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t103 - rSges(3,2) * t102 + t138) + g(2) * (t101 * rSges(3,1) - t100 * rSges(3,2) + t111) + g(3) * (t143 * rSges(3,3) + t145) + (g(1) * rSges(3,3) * t118 + g(3) * (rSges(3,1) * t117 + rSges(3,2) * t121) + g(2) * (-rSges(3,3) - pkin(8)) * t122) * t113) - m(4) * (g(1) * (t92 * rSges(4,1) - t91 * rSges(4,2) + t149 * t102 + t135) + g(2) * (t90 * rSges(4,1) - t89 * rSges(4,2) + t149 * t100 + t129) + g(3) * (t99 * rSges(4,1) - t98 * rSges(4,2) - t149 * t140 + t136)) - m(5) * (g(1) * (t92 * rSges(5,1) + t150 * t102 + t144 * t91 + t130) + g(2) * (t90 * rSges(5,1) + t150 * t100 + t144 * t89 + t127) + g(3) * (t99 * rSges(5,1) - t150 * t140 + t144 * t98 + t131)) - m(6) * (g(3) * (t82 * rSges(6,1) - t81 * rSges(6,2) + t148 * t140 + t124) + (t78 * rSges(6,1) - t77 * rSges(6,2) + t137 * t100 + t123) * g(2) + (t80 * rSges(6,1) - t79 * rSges(6,2) + t137 * t102 + t125) * g(1)) - m(7) * (g(1) * (t126 * t102 + t128 * t80 + t147 * t79 + t125) + g(2) * (t126 * t100 + t128 * t78 + t147 * t77 + t123) + g(3) * (t82 * pkin(5) - pkin(9) * t140 + (t114 * t140 + t119 * t82) * rSges(7,1) + (-t114 * t82 + t119 * t140) * rSges(7,2) + t147 * t81 + t124));
U  = t1;

% Calculate potential energy for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:07:39
% EndTime: 2019-03-09 12:07:40
% DurationCPUTime: 0.46s
% Computational Cost: add. (418->115), mult. (938->155), div. (0->0), fcn. (1177->12), ass. (0->60)
t161 = rSges(7,1) + pkin(5);
t160 = rSges(7,2) + pkin(10);
t158 = rSges(5,3) + pkin(9);
t157 = rSges(6,3) + pkin(10);
t128 = cos(pkin(6));
t156 = t128 * pkin(8) + pkin(7);
t155 = rSges(7,3) + qJ(6);
t126 = sin(pkin(6));
t131 = sin(qJ(2));
t154 = t126 * t131;
t132 = sin(qJ(1));
t153 = t126 * t132;
t136 = cos(qJ(1));
t152 = t126 * t136;
t127 = cos(pkin(11));
t135 = cos(qJ(2));
t151 = t127 * t135;
t150 = t131 * t136;
t149 = t132 * t131;
t148 = t132 * t135;
t147 = t135 * t136;
t113 = pkin(2) * t128 * t131 + (-pkin(8) - qJ(3)) * t126;
t122 = pkin(2) * t135 + pkin(1);
t146 = t136 * t113 + t132 * t122;
t125 = sin(pkin(11));
t141 = t125 * t135 + t127 * t131;
t112 = t141 * t128;
t115 = -t125 * t131 + t151;
t100 = t112 * t136 + t132 * t115;
t145 = t100 * pkin(3) + t146;
t144 = pkin(2) * t154 + t128 * qJ(3) + t156;
t102 = -t132 * t112 + t115 * t136;
t118 = t136 * t122;
t143 = t102 * pkin(3) - t113 * t132 + t118;
t111 = t141 * t126;
t142 = t111 * pkin(3) + t144;
t130 = sin(qJ(4));
t134 = cos(qJ(4));
t92 = t100 * t134 - t130 * t152;
t139 = t115 * t128;
t99 = -t132 * t141 + t136 * t139;
t140 = t92 * pkin(4) - pkin(9) * t99 + t145;
t101 = -t132 * t139 - t136 * t141;
t94 = t102 * t134 + t130 * t153;
t138 = t94 * pkin(4) - pkin(9) * t101 + t143;
t105 = t111 * t134 + t128 * t130;
t110 = t125 * t154 - t126 * t151;
t137 = t105 * pkin(4) + pkin(9) * t110 + t142;
t133 = cos(qJ(5));
t129 = sin(qJ(5));
t104 = t111 * t130 - t128 * t134;
t93 = t102 * t130 - t134 * t153;
t91 = t100 * t130 + t134 * t152;
t88 = t105 * t133 + t110 * t129;
t87 = t105 * t129 - t110 * t133;
t86 = -t101 * t129 + t133 * t94;
t85 = t101 * t133 + t129 * t94;
t84 = -t129 * t99 + t133 * t92;
t83 = t129 * t92 + t99 * t133;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t136 - t132 * rSges(2,2)) + g(2) * (t132 * rSges(2,1) + rSges(2,2) * t136) + g(3) * (pkin(7) + rSges(2,3))) - m(3) * (g(1) * (t136 * pkin(1) + (-t128 * t149 + t147) * rSges(3,1) + (-t128 * t148 - t150) * rSges(3,2)) + g(2) * (t132 * pkin(1) + (t128 * t150 + t148) * rSges(3,1) + (t128 * t147 - t149) * rSges(3,2)) + g(3) * (rSges(3,3) * t128 + t156) + (g(3) * (rSges(3,1) * t131 + rSges(3,2) * t135) + (g(1) * t132 - g(2) * t136) * (rSges(3,3) + pkin(8))) * t126) - m(4) * (g(1) * (rSges(4,1) * t102 + rSges(4,2) * t101 + t118 + (rSges(4,3) * t126 - t113) * t132) + g(2) * (t100 * rSges(4,1) + t99 * rSges(4,2) - rSges(4,3) * t152 + t146) + g(3) * (rSges(4,1) * t111 - rSges(4,2) * t110 + rSges(4,3) * t128 + t144)) - m(5) * (g(1) * (rSges(5,1) * t94 - rSges(5,2) * t93 - t158 * t101 + t143) + g(2) * (rSges(5,1) * t92 - rSges(5,2) * t91 - t158 * t99 + t145) + g(3) * (rSges(5,1) * t105 - rSges(5,2) * t104 + t158 * t110 + t142)) - m(6) * (g(1) * (rSges(6,1) * t86 - rSges(6,2) * t85 + t157 * t93 + t138) + g(2) * (rSges(6,1) * t84 - rSges(6,2) * t83 + t157 * t91 + t140) + g(3) * (rSges(6,1) * t88 - rSges(6,2) * t87 + t157 * t104 + t137)) - m(7) * (g(1) * (t155 * t85 + t160 * t93 + t161 * t86 + t138) + g(2) * (t155 * t83 + t160 * t91 + t161 * t84 + t140) + g(3) * (t160 * t104 + t155 * t87 + t161 * t88 + t137));
U  = t1;

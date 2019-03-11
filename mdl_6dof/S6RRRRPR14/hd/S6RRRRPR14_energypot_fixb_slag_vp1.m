% Calculate potential energy for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR14_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_energypot_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR14_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:08:35
% EndTime: 2019-03-10 00:08:35
% DurationCPUTime: 0.54s
% Computational Cost: add. (563->128), mult. (1388->175), div. (0->0), fcn. (1753->16), ass. (0->62)
t134 = cos(pkin(6));
t139 = sin(qJ(1));
t140 = cos(qJ(2));
t160 = t139 * t140;
t138 = sin(qJ(2));
t141 = cos(qJ(1));
t162 = t138 * t141;
t113 = -t134 * t160 - t162;
t130 = sin(pkin(7));
t133 = cos(pkin(7));
t131 = sin(pkin(6));
t166 = t131 * t139;
t150 = -t113 * t130 + t133 * t166;
t165 = t131 * t140;
t149 = -t130 * t165 + t133 * t134;
t174 = rSges(5,3) + pkin(11);
t173 = cos(qJ(3));
t172 = cos(qJ(4));
t171 = t134 * pkin(9) + pkin(8);
t170 = qJ(5) + rSges(6,3);
t169 = pkin(12) + qJ(5) + rSges(7,3);
t167 = t131 * t138;
t164 = t131 * t141;
t161 = t139 * t138;
t159 = t140 * t141;
t158 = t141 * pkin(1) + pkin(9) * t166;
t155 = t130 * t173;
t154 = t133 * t173;
t153 = t131 * t155;
t128 = pkin(13) + qJ(6);
t123 = sin(t128);
t124 = cos(t128);
t132 = cos(pkin(13));
t152 = rSges(7,1) * t124 - rSges(7,2) * t123 + pkin(5) * t132 + pkin(4);
t111 = t134 * t159 - t161;
t151 = -t111 * t130 - t133 * t164;
t114 = -t134 * t161 + t159;
t148 = t114 * pkin(2) + t150 * pkin(10) + t158;
t129 = sin(pkin(13));
t147 = rSges(7,1) * t123 + rSges(7,2) * t124 + pkin(5) * t129 + pkin(11);
t146 = pkin(2) * t167 + t149 * pkin(10) + t171;
t137 = sin(qJ(3));
t100 = t114 * t173 + (t113 * t133 + t130 * t166) * t137;
t145 = t100 * pkin(3) + t148;
t105 = t134 * t130 * t137 + (t133 * t137 * t140 + t173 * t138) * t131;
t144 = t105 * pkin(3) + t146;
t112 = t134 * t162 + t160;
t126 = t139 * pkin(1);
t143 = t112 * pkin(2) - pkin(9) * t164 + t151 * pkin(10) + t126;
t98 = t112 * t173 + (t111 * t133 - t130 * t164) * t137;
t142 = t98 * pkin(3) + t143;
t136 = sin(qJ(4));
t104 = -t134 * t155 + t137 * t167 - t154 * t165;
t99 = -t113 * t154 + t114 * t137 - t139 * t153;
t97 = -t111 * t154 + t112 * t137 + t141 * t153;
t96 = t105 * t172 + t149 * t136;
t95 = t105 * t136 - t149 * t172;
t92 = t100 * t172 + t150 * t136;
t91 = t100 * t136 - t150 * t172;
t90 = t151 * t136 + t98 * t172;
t89 = t136 * t98 - t151 * t172;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t141 - t139 * rSges(2,2)) + g(2) * (t139 * rSges(2,1) + rSges(2,2) * t141) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t114 + rSges(3,2) * t113 + t158) + g(2) * (t112 * rSges(3,1) + t111 * rSges(3,2) + t126) + g(3) * (rSges(3,3) * t134 + t171) + (g(1) * rSges(3,3) * t139 + g(3) * (rSges(3,1) * t138 + rSges(3,2) * t140) + g(2) * (-rSges(3,3) - pkin(9)) * t141) * t131) - m(4) * (g(1) * (t100 * rSges(4,1) - t99 * rSges(4,2) + t150 * rSges(4,3) + t148) + g(2) * (t98 * rSges(4,1) - t97 * rSges(4,2) + t151 * rSges(4,3) + t143) + g(3) * (t105 * rSges(4,1) - t104 * rSges(4,2) + t149 * rSges(4,3) + t146)) - m(5) * (g(1) * (rSges(5,1) * t92 - rSges(5,2) * t91 + t174 * t99 + t145) + g(2) * (t90 * rSges(5,1) - t89 * rSges(5,2) + t174 * t97 + t142) + g(3) * (rSges(5,1) * t96 - rSges(5,2) * t95 + t174 * t104 + t144)) - m(6) * (g(1) * (t92 * pkin(4) + t99 * pkin(11) + (t129 * t99 + t132 * t92) * rSges(6,1) + (-t129 * t92 + t132 * t99) * rSges(6,2) + t170 * t91 + t145) + g(2) * (t90 * pkin(4) + t97 * pkin(11) + (t129 * t97 + t132 * t90) * rSges(6,1) + (-t129 * t90 + t132 * t97) * rSges(6,2) + t170 * t89 + t142) + g(3) * (t96 * pkin(4) + t104 * pkin(11) + (t104 * t129 + t132 * t96) * rSges(6,1) + (t104 * t132 - t129 * t96) * rSges(6,2) + t170 * t95 + t144)) - m(7) * (g(1) * (t147 * t99 + t152 * t92 + t169 * t91 + t145) + g(2) * (t147 * t97 + t152 * t90 + t169 * t89 + t142) + g(3) * (t147 * t104 + t152 * t96 + t169 * t95 + t144));
U  = t1;

% Calculate potential energy for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR9_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_energypot_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR9_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:37
% EndTime: 2019-03-10 05:16:37
% DurationCPUTime: 0.54s
% Computational Cost: add. (563->128), mult. (1388->175), div. (0->0), fcn. (1753->16), ass. (0->62)
t132 = cos(pkin(6));
t137 = sin(qJ(1));
t139 = cos(qJ(2));
t160 = t137 * t139;
t136 = sin(qJ(2));
t140 = cos(qJ(1));
t161 = t136 * t140;
t113 = -t132 * t160 - t161;
t129 = sin(pkin(7));
t131 = cos(pkin(7));
t130 = sin(pkin(6));
t166 = t130 * t137;
t150 = -t113 * t129 + t131 * t166;
t165 = t130 * t139;
t149 = -t129 * t165 + t131 * t132;
t174 = rSges(5,3) + pkin(11);
t173 = pkin(12) + rSges(6,3);
t172 = cos(qJ(3));
t171 = cos(qJ(4));
t170 = t132 * pkin(9) + pkin(8);
t169 = pkin(13) + pkin(12) + rSges(7,3);
t167 = t130 * t136;
t164 = t130 * t140;
t162 = t136 * t137;
t159 = t139 * t140;
t158 = t140 * pkin(1) + pkin(9) * t166;
t155 = t129 * t172;
t154 = t131 * t172;
t153 = t130 * t155;
t128 = qJ(5) + qJ(6);
t123 = sin(t128);
t124 = cos(t128);
t138 = cos(qJ(5));
t152 = rSges(7,1) * t124 - rSges(7,2) * t123 + pkin(5) * t138 + pkin(4);
t111 = t132 * t159 - t162;
t151 = -t111 * t129 - t131 * t164;
t114 = -t132 * t162 + t159;
t148 = t114 * pkin(2) + pkin(10) * t150 + t158;
t133 = sin(qJ(5));
t147 = rSges(7,1) * t123 + rSges(7,2) * t124 + pkin(5) * t133 + pkin(11);
t146 = pkin(2) * t167 + pkin(10) * t149 + t170;
t135 = sin(qJ(3));
t100 = t114 * t172 + (t113 * t131 + t129 * t166) * t135;
t145 = t100 * pkin(3) + t148;
t105 = t132 * t129 * t135 + (t131 * t135 * t139 + t172 * t136) * t130;
t144 = t105 * pkin(3) + t146;
t112 = t132 * t161 + t160;
t126 = t137 * pkin(1);
t143 = t112 * pkin(2) - pkin(9) * t164 + t151 * pkin(10) + t126;
t98 = t112 * t172 + (t111 * t131 - t129 * t164) * t135;
t142 = t98 * pkin(3) + t143;
t134 = sin(qJ(4));
t104 = -t132 * t155 + t135 * t167 - t154 * t165;
t99 = -t113 * t154 + t114 * t135 - t137 * t153;
t97 = -t111 * t154 + t112 * t135 + t140 * t153;
t96 = t105 * t171 + t149 * t134;
t95 = t105 * t134 - t149 * t171;
t92 = t100 * t171 + t150 * t134;
t91 = t100 * t134 - t150 * t171;
t90 = t151 * t134 + t98 * t171;
t89 = t134 * t98 - t151 * t171;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t140 - rSges(2,2) * t137) + g(2) * (rSges(2,1) * t137 + rSges(2,2) * t140) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t114 + rSges(3,2) * t113 + t158) + g(2) * (t112 * rSges(3,1) + t111 * rSges(3,2) + t126) + g(3) * (t132 * rSges(3,3) + t170) + (g(1) * rSges(3,3) * t137 + g(3) * (rSges(3,1) * t136 + rSges(3,2) * t139) + g(2) * (-rSges(3,3) - pkin(9)) * t140) * t130) - m(4) * (g(1) * (t100 * rSges(4,1) - t99 * rSges(4,2) + t150 * rSges(4,3) + t148) + g(2) * (t98 * rSges(4,1) - t97 * rSges(4,2) + t151 * rSges(4,3) + t143) + g(3) * (t105 * rSges(4,1) - t104 * rSges(4,2) + t149 * rSges(4,3) + t146)) - m(5) * (g(1) * (t92 * rSges(5,1) - t91 * rSges(5,2) + t174 * t99 + t145) + g(2) * (rSges(5,1) * t90 - rSges(5,2) * t89 + t174 * t97 + t142) + g(3) * (t96 * rSges(5,1) - t95 * rSges(5,2) + t174 * t104 + t144)) - m(6) * (g(1) * (t92 * pkin(4) + t99 * pkin(11) + (t133 * t99 + t138 * t92) * rSges(6,1) + (-t133 * t92 + t138 * t99) * rSges(6,2) + t173 * t91 + t145) + g(2) * (t90 * pkin(4) + t97 * pkin(11) + (t133 * t97 + t138 * t90) * rSges(6,1) + (-t133 * t90 + t138 * t97) * rSges(6,2) + t173 * t89 + t142) + g(3) * (t96 * pkin(4) + t104 * pkin(11) + (t104 * t133 + t138 * t96) * rSges(6,1) + (t104 * t138 - t133 * t96) * rSges(6,2) + t173 * t95 + t144)) - m(7) * (g(1) * (t147 * t99 + t152 * t92 + t169 * t91 + t145) + g(2) * (t147 * t97 + t152 * t90 + t169 * t89 + t142) + g(3) * (t147 * t104 + t152 * t96 + t169 * t95 + t144));
U  = t1;

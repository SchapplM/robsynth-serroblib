% Calculate potential energy for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_energypot_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:10
% EndTime: 2019-03-08 19:03:10
% DurationCPUTime: 0.54s
% Computational Cost: add. (563->128), mult. (1388->178), div. (0->0), fcn. (1753->16), ass. (0->61)
t130 = sin(pkin(7));
t134 = cos(pkin(7));
t135 = cos(pkin(6));
t131 = sin(pkin(6));
t132 = cos(pkin(13));
t163 = t131 * t132;
t148 = -t130 * t163 + t134 * t135;
t128 = sin(pkin(13));
t133 = cos(pkin(12));
t129 = sin(pkin(12));
t164 = t129 * t135;
t112 = -t128 * t133 - t132 * t164;
t161 = t131 * t134;
t149 = -t112 * t130 + t129 * t161;
t172 = rSges(5,3) + pkin(9);
t171 = pkin(10) + rSges(6,3);
t170 = cos(qJ(3));
t169 = cos(qJ(4));
t168 = pkin(11) + pkin(10) + rSges(7,3);
t166 = t128 * t131;
t165 = t129 * t131;
t162 = t131 * t133;
t160 = t133 * t135;
t158 = t135 * qJ(2) + qJ(1);
t157 = t133 * pkin(1) + qJ(2) * t165;
t154 = t130 * t170;
t153 = t134 * t170;
t152 = t131 * t154;
t127 = qJ(5) + qJ(6);
t123 = sin(t127);
t124 = cos(t127);
t139 = cos(qJ(5));
t151 = rSges(7,1) * t124 - rSges(7,2) * t123 + pkin(5) * t139 + pkin(4);
t110 = -t128 * t129 + t132 * t160;
t150 = -t110 * t130 - t133 * t161;
t113 = -t128 * t164 + t132 * t133;
t147 = t113 * pkin(2) + t149 * pkin(8) + t157;
t136 = sin(qJ(5));
t146 = t123 * rSges(7,1) + t124 * rSges(7,2) + t136 * pkin(5) + pkin(9);
t138 = sin(qJ(3));
t97 = t113 * t170 + (t112 * t134 + t130 * t165) * t138;
t145 = t97 * pkin(3) + t147;
t144 = pkin(2) * t166 + t148 * pkin(8) + t158;
t104 = t135 * t130 * t138 + (t132 * t134 * t138 + t170 * t128) * t131;
t143 = t104 * pkin(3) + t144;
t111 = t128 * t160 + t129 * t132;
t125 = t129 * pkin(1);
t142 = t111 * pkin(2) + t150 * pkin(8) - qJ(2) * t162 + t125;
t95 = t111 * t170 + (t110 * t134 - t130 * t162) * t138;
t141 = t95 * pkin(3) + t142;
t137 = sin(qJ(4));
t103 = -t135 * t154 + t138 * t166 - t153 * t163;
t99 = t104 * t169 + t148 * t137;
t98 = t104 * t137 - t148 * t169;
t96 = -t112 * t153 + t113 * t138 - t129 * t152;
t94 = -t110 * t153 + t111 * t138 + t133 * t152;
t91 = t149 * t137 + t97 * t169;
t90 = t137 * t97 - t149 * t169;
t89 = t150 * t137 + t95 * t169;
t88 = t137 * t95 - t150 * t169;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t133 - rSges(2,2) * t129) + g(2) * (rSges(2,1) * t129 + rSges(2,2) * t133) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t113 + rSges(3,2) * t112 + t157) + g(2) * (rSges(3,1) * t111 + rSges(3,2) * t110 + t125) + g(3) * (rSges(3,3) * t135 + t158) + (g(1) * rSges(3,3) * t129 + g(3) * (rSges(3,1) * t128 + rSges(3,2) * t132) + g(2) * (-rSges(3,3) - qJ(2)) * t133) * t131) - m(4) * (g(1) * (t97 * rSges(4,1) - t96 * rSges(4,2) + t149 * rSges(4,3) + t147) + g(2) * (t95 * rSges(4,1) - t94 * rSges(4,2) + t150 * rSges(4,3) + t142) + g(3) * (t104 * rSges(4,1) - t103 * rSges(4,2) + t148 * rSges(4,3) + t144)) - m(5) * (g(1) * (rSges(5,1) * t91 - rSges(5,2) * t90 + t172 * t96 + t145) + g(2) * (rSges(5,1) * t89 - rSges(5,2) * t88 + t172 * t94 + t141) + g(3) * (rSges(5,1) * t99 - rSges(5,2) * t98 + t172 * t103 + t143)) - m(6) * (g(1) * (t91 * pkin(4) + t96 * pkin(9) + (t136 * t96 + t139 * t91) * rSges(6,1) + (-t136 * t91 + t139 * t96) * rSges(6,2) + t171 * t90 + t145) + g(2) * (t89 * pkin(4) + t94 * pkin(9) + (t136 * t94 + t139 * t89) * rSges(6,1) + (-t136 * t89 + t139 * t94) * rSges(6,2) + t171 * t88 + t141) + g(3) * (t99 * pkin(4) + t103 * pkin(9) + (t103 * t136 + t139 * t99) * rSges(6,1) + (t103 * t139 - t136 * t99) * rSges(6,2) + t171 * t98 + t143)) - m(7) * (g(1) * (t146 * t96 + t151 * t91 + t168 * t90 + t145) + g(2) * (t146 * t94 + t151 * t89 + t168 * t88 + t141) + g(3) * (t146 * t103 + t151 * t99 + t168 * t98 + t143));
U  = t1;

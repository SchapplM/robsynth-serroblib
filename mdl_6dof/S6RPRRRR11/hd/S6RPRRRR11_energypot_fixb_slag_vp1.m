% Calculate potential energy for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_energypot_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR11_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:10
% EndTime: 2019-03-09 07:38:11
% DurationCPUTime: 0.53s
% Computational Cost: add. (563->128), mult. (1388->177), div. (0->0), fcn. (1753->16), ass. (0->60)
t128 = sin(pkin(13));
t131 = cos(pkin(13));
t139 = cos(qJ(1));
t133 = cos(pkin(6));
t137 = sin(qJ(1));
t159 = t133 * t137;
t112 = -t128 * t139 - t131 * t159;
t129 = sin(pkin(7));
t132 = cos(pkin(7));
t130 = sin(pkin(6));
t162 = t130 * t137;
t149 = -t112 * t129 + t132 * t162;
t163 = t130 * t131;
t148 = -t129 * t163 + t132 * t133;
t171 = rSges(5,3) + pkin(10);
t170 = pkin(11) + rSges(6,3);
t169 = cos(qJ(3));
t168 = cos(qJ(4));
t167 = t133 * qJ(2) + pkin(8);
t166 = pkin(12) + pkin(11) + rSges(7,3);
t164 = t128 * t130;
t161 = t130 * t139;
t158 = t133 * t139;
t157 = t139 * pkin(1) + qJ(2) * t162;
t154 = t129 * t169;
t153 = t132 * t169;
t152 = t130 * t154;
t127 = qJ(5) + qJ(6);
t123 = sin(t127);
t124 = cos(t127);
t138 = cos(qJ(5));
t151 = t124 * rSges(7,1) - t123 * rSges(7,2) + pkin(5) * t138 + pkin(4);
t110 = -t128 * t137 + t131 * t158;
t150 = -t110 * t129 - t132 * t161;
t113 = -t128 * t159 + t131 * t139;
t147 = t113 * pkin(2) + t149 * pkin(9) + t157;
t134 = sin(qJ(5));
t146 = t123 * rSges(7,1) + t124 * rSges(7,2) + t134 * pkin(5) + pkin(10);
t145 = pkin(2) * t164 + t148 * pkin(9) + t167;
t136 = sin(qJ(3));
t99 = t113 * t169 + (t112 * t132 + t129 * t162) * t136;
t144 = t99 * pkin(3) + t147;
t104 = t133 * t129 * t136 + (t131 * t132 * t136 + t169 * t128) * t130;
t143 = t104 * pkin(3) + t145;
t111 = t128 * t158 + t131 * t137;
t125 = t137 * pkin(1);
t142 = t111 * pkin(2) + t150 * pkin(9) - qJ(2) * t161 + t125;
t97 = t111 * t169 + (t110 * t132 - t129 * t161) * t136;
t141 = t97 * pkin(3) + t142;
t135 = sin(qJ(4));
t103 = -t133 * t154 + t136 * t164 - t153 * t163;
t98 = -t112 * t153 + t113 * t136 - t137 * t152;
t96 = -t110 * t153 + t111 * t136 + t139 * t152;
t95 = t104 * t168 + t148 * t135;
t94 = t104 * t135 - t148 * t168;
t91 = t149 * t135 + t99 * t168;
t90 = t135 * t99 - t149 * t168;
t89 = t150 * t135 + t97 * t168;
t88 = t135 * t97 - t150 * t168;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t139 - rSges(2,2) * t137) + g(2) * (rSges(2,1) * t137 + rSges(2,2) * t139) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t113 + rSges(3,2) * t112 + t157) + g(2) * (rSges(3,1) * t111 + rSges(3,2) * t110 + t125) + g(3) * (rSges(3,3) * t133 + t167) + (g(1) * rSges(3,3) * t137 + g(3) * (rSges(3,1) * t128 + rSges(3,2) * t131) + g(2) * (-rSges(3,3) - qJ(2)) * t139) * t130) - m(4) * (g(1) * (t99 * rSges(4,1) - t98 * rSges(4,2) + t149 * rSges(4,3) + t147) + g(2) * (t97 * rSges(4,1) - t96 * rSges(4,2) + t150 * rSges(4,3) + t142) + g(3) * (t104 * rSges(4,1) - t103 * rSges(4,2) + t148 * rSges(4,3) + t145)) - m(5) * (g(1) * (rSges(5,1) * t91 - rSges(5,2) * t90 + t171 * t98 + t144) + g(2) * (rSges(5,1) * t89 - rSges(5,2) * t88 + t171 * t96 + t141) + g(3) * (rSges(5,1) * t95 - rSges(5,2) * t94 + t171 * t103 + t143)) - m(6) * (g(1) * (t91 * pkin(4) + t98 * pkin(10) + (t134 * t98 + t138 * t91) * rSges(6,1) + (-t134 * t91 + t138 * t98) * rSges(6,2) + t170 * t90 + t144) + g(2) * (t89 * pkin(4) + t96 * pkin(10) + (t134 * t96 + t138 * t89) * rSges(6,1) + (-t134 * t89 + t138 * t96) * rSges(6,2) + t170 * t88 + t141) + g(3) * (t95 * pkin(4) + t103 * pkin(10) + (t103 * t134 + t138 * t95) * rSges(6,1) + (t103 * t138 - t134 * t95) * rSges(6,2) + t170 * t94 + t143)) - m(7) * (g(1) * (t146 * t98 + t151 * t91 + t166 * t90 + t144) + g(2) * (t146 * t96 + t151 * t89 + t166 * t88 + t141) + g(3) * (t146 * t103 + t151 * t95 + t166 * t94 + t143));
U  = t1;

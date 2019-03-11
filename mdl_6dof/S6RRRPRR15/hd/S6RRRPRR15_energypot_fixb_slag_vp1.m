% Calculate potential energy for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR15_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR15_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:26
% EndTime: 2019-03-09 20:24:26
% DurationCPUTime: 0.51s
% Computational Cost: add. (477->127), mult. (1189->174), div. (0->0), fcn. (1483->14), ass. (0->62)
t137 = cos(pkin(6));
t142 = sin(qJ(1));
t145 = cos(qJ(2));
t163 = t142 * t145;
t141 = sin(qJ(2));
t146 = cos(qJ(1));
t165 = t141 * t146;
t121 = -t137 * t163 - t165;
t134 = sin(pkin(7));
t136 = cos(pkin(7));
t135 = sin(pkin(6));
t170 = t135 * t142;
t112 = -t121 * t134 + t136 * t170;
t169 = t135 * t145;
t118 = -t134 * t169 + t136 * t137;
t178 = rSges(6,3) + pkin(11);
t177 = pkin(12) + rSges(7,3);
t176 = cos(qJ(3));
t175 = t137 * pkin(9) + pkin(8);
t174 = rSges(5,3) + qJ(4);
t140 = sin(qJ(3));
t172 = t134 * t140;
t171 = t135 * t141;
t168 = t135 * t146;
t166 = t136 * t140;
t164 = t142 * t141;
t162 = t145 * t146;
t161 = t146 * pkin(1) + pkin(9) * t170;
t158 = t134 * t176;
t157 = t136 * t176;
t156 = t135 * t158;
t119 = t137 * t162 - t164;
t111 = -t119 * t134 - t136 * t168;
t122 = -t137 * t164 + t162;
t155 = t122 * pkin(2) + t112 * pkin(10) + t161;
t154 = pkin(2) * t171 + t118 * pkin(10) + t175;
t104 = t122 * t176 + (t121 * t136 + t134 * t170) * t140;
t153 = t104 * pkin(3) + t155;
t108 = t137 * t172 + (t141 * t176 + t145 * t166) * t135;
t152 = t108 * pkin(3) + t154;
t103 = -t121 * t157 + t122 * t140 - t142 * t156;
t151 = t112 * pkin(4) + t103 * qJ(4) + t153;
t107 = -t137 * t158 + t140 * t171 - t157 * t169;
t150 = t118 * pkin(4) + t107 * qJ(4) + t152;
t120 = t137 * t165 + t163;
t132 = t142 * pkin(1);
t149 = t120 * pkin(2) - pkin(9) * t168 + pkin(10) * t111 + t132;
t102 = t119 * t166 + t120 * t176 - t168 * t172;
t148 = t102 * pkin(3) + t149;
t101 = -t119 * t157 + t120 * t140 + t146 * t156;
t147 = t111 * pkin(4) + t101 * qJ(4) + t148;
t144 = cos(qJ(5));
t143 = cos(qJ(6));
t139 = sin(qJ(5));
t138 = sin(qJ(6));
t100 = t107 * t139 + t118 * t144;
t99 = -t107 * t144 + t118 * t139;
t94 = t103 * t139 + t112 * t144;
t93 = -t103 * t144 + t112 * t139;
t92 = t101 * t139 + t111 * t144;
t91 = -t101 * t144 + t111 * t139;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t146 - t142 * rSges(2,2)) + g(2) * (t142 * rSges(2,1) + rSges(2,2) * t146) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t122 + rSges(3,2) * t121 + t161) + g(2) * (t120 * rSges(3,1) + t119 * rSges(3,2) + t132) + g(3) * (rSges(3,3) * t137 + t175) + (g(1) * rSges(3,3) * t142 + g(3) * (rSges(3,1) * t141 + rSges(3,2) * t145) + g(2) * (-rSges(3,3) - pkin(9)) * t146) * t135) - m(4) * (g(1) * (rSges(4,1) * t104 - rSges(4,2) * t103 + rSges(4,3) * t112 + t155) + g(2) * (t102 * rSges(4,1) - t101 * rSges(4,2) + t111 * rSges(4,3) + t149) + g(3) * (rSges(4,1) * t108 - rSges(4,2) * t107 + rSges(4,3) * t118 + t154)) - m(5) * (g(1) * (rSges(5,1) * t112 - rSges(5,2) * t104 + t103 * t174 + t153) + g(2) * (t111 * rSges(5,1) - t102 * rSges(5,2) + t101 * t174 + t148) + g(3) * (rSges(5,1) * t118 - rSges(5,2) * t108 + t107 * t174 + t152)) - m(6) * (g(1) * (rSges(6,1) * t94 - rSges(6,2) * t93 + t104 * t178 + t151) + g(2) * (t92 * rSges(6,1) - t91 * rSges(6,2) + t102 * t178 + t147) + g(3) * (rSges(6,1) * t100 - rSges(6,2) * t99 + t108 * t178 + t150)) - m(7) * (g(1) * (t94 * pkin(5) + t104 * pkin(11) + (t104 * t138 + t143 * t94) * rSges(7,1) + (t104 * t143 - t138 * t94) * rSges(7,2) + t177 * t93 + t151) + g(2) * (t92 * pkin(5) + t102 * pkin(11) + (t102 * t138 + t143 * t92) * rSges(7,1) + (t102 * t143 - t138 * t92) * rSges(7,2) + t177 * t91 + t147) + g(3) * (t100 * pkin(5) + t108 * pkin(11) + (t100 * t143 + t108 * t138) * rSges(7,1) + (-t100 * t138 + t108 * t143) * rSges(7,2) + t177 * t99 + t150));
U  = t1;

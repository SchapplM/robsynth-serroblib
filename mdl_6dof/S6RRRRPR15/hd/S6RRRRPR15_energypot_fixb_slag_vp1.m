% Calculate potential energy for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR15_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR15_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:05
% EndTime: 2019-03-10 00:34:06
% DurationCPUTime: 0.47s
% Computational Cost: add. (532->118), mult. (1342->158), div. (0->0), fcn. (1691->14), ass. (0->63)
t133 = cos(pkin(6));
t138 = sin(qJ(1));
t140 = cos(qJ(2));
t163 = t138 * t140;
t137 = sin(qJ(2));
t141 = cos(qJ(1));
t165 = t137 * t141;
t118 = -t133 * t163 - t165;
t130 = sin(pkin(7));
t132 = cos(pkin(7));
t131 = sin(pkin(6));
t169 = t131 * t138;
t153 = -t118 * t130 + t132 * t169;
t168 = t131 * t140;
t152 = -t130 * t168 + t132 * t133;
t178 = rSges(6,1) + pkin(11);
t177 = rSges(5,3) + pkin(11);
t176 = pkin(12) + rSges(7,3);
t175 = cos(qJ(3));
t174 = cos(qJ(4));
t173 = t133 * pkin(9) + pkin(8);
t172 = rSges(6,3) + qJ(5);
t170 = t131 * t137;
t167 = t131 * t141;
t164 = t138 * t137;
t162 = t140 * t141;
t161 = t141 * pkin(1) + pkin(9) * t169;
t158 = t130 * t175;
t157 = t132 * t175;
t156 = t131 * t158;
t134 = sin(qJ(6));
t139 = cos(qJ(6));
t155 = t134 * rSges(7,1) + t139 * rSges(7,2) + qJ(5);
t116 = t133 * t162 - t164;
t154 = -t116 * t130 - t132 * t167;
t151 = t139 * rSges(7,1) - t134 * rSges(7,2) + pkin(5) + pkin(11);
t119 = -t133 * t164 + t162;
t150 = t119 * pkin(2) + t153 * pkin(10) + t161;
t149 = pkin(2) * t170 + t152 * pkin(10) + t173;
t136 = sin(qJ(3));
t105 = t119 * t175 + (t118 * t132 + t130 * t169) * t136;
t148 = t105 * pkin(3) + t150;
t110 = t133 * t130 * t136 + (t132 * t136 * t140 + t175 * t137) * t131;
t147 = t110 * pkin(3) + t149;
t135 = sin(qJ(4));
t96 = t105 * t174 + t153 * t135;
t146 = t96 * pkin(4) + t148;
t101 = t110 * t174 + t152 * t135;
t145 = t101 * pkin(4) + t147;
t117 = t133 * t165 + t163;
t128 = t138 * pkin(1);
t144 = t117 * pkin(2) - pkin(9) * t167 + t154 * pkin(10) + t128;
t103 = t117 * t175 + (t116 * t132 - t130 * t167) * t136;
t143 = t103 * pkin(3) + t144;
t94 = t103 * t174 + t154 * t135;
t142 = t94 * pkin(4) + t143;
t109 = -t133 * t158 + t136 * t170 - t157 * t168;
t104 = -t118 * t157 + t119 * t136 - t138 * t156;
t102 = -t116 * t157 + t117 * t136 + t141 * t156;
t100 = t110 * t135 - t152 * t174;
t95 = t105 * t135 - t153 * t174;
t93 = t103 * t135 - t154 * t174;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t141 - t138 * rSges(2,2)) + g(2) * (t138 * rSges(2,1) + rSges(2,2) * t141) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t119 + rSges(3,2) * t118 + t161) + g(2) * (t117 * rSges(3,1) + t116 * rSges(3,2) + t128) + g(3) * (rSges(3,3) * t133 + t173) + (g(1) * rSges(3,3) * t138 + g(3) * (rSges(3,1) * t137 + rSges(3,2) * t140) + g(2) * (-rSges(3,3) - pkin(9)) * t141) * t131) - m(4) * (g(1) * (t105 * rSges(4,1) - t104 * rSges(4,2) + t153 * rSges(4,3) + t150) + g(2) * (t103 * rSges(4,1) - t102 * rSges(4,2) + t154 * rSges(4,3) + t144) + g(3) * (t110 * rSges(4,1) - t109 * rSges(4,2) + t152 * rSges(4,3) + t149)) - m(5) * (g(1) * (rSges(5,1) * t96 - rSges(5,2) * t95 + t177 * t104 + t148) + g(2) * (t94 * rSges(5,1) - t93 * rSges(5,2) + t177 * t102 + t143) + g(3) * (rSges(5,1) * t101 - rSges(5,2) * t100 + t177 * t109 + t147)) - m(6) * (g(1) * (-rSges(6,2) * t96 + t178 * t104 + t172 * t95 + t146) + g(2) * (-t94 * rSges(6,2) + t178 * t102 + t172 * t93 + t142) + g(3) * (-rSges(6,2) * t101 + t172 * t100 + t178 * t109 + t145)) - m(7) * (g(1) * (t151 * t104 + t155 * t95 + t176 * t96 + t146) + g(2) * (t151 * t102 + t155 * t93 + t176 * t94 + t142) + g(3) * (t155 * t100 + t176 * t101 + t151 * t109 + t145));
U  = t1;

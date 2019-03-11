% Calculate potential energy for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:05
% EndTime: 2019-03-08 23:48:06
% DurationCPUTime: 0.47s
% Computational Cost: add. (532->118), mult. (1342->161), div. (0->0), fcn. (1691->14), ass. (0->62)
t130 = sin(pkin(7));
t133 = cos(pkin(7));
t134 = cos(pkin(6));
t131 = sin(pkin(6));
t140 = cos(qJ(2));
t165 = t131 * t140;
t151 = -t130 * t165 + t134 * t133;
t129 = sin(pkin(12));
t132 = cos(pkin(12));
t138 = sin(qJ(2));
t162 = t134 * t140;
t117 = -t129 * t162 - t132 * t138;
t167 = t131 * t133;
t152 = -t117 * t130 + t129 * t167;
t176 = rSges(6,1) + pkin(10);
t175 = rSges(5,3) + pkin(10);
t174 = pkin(11) + rSges(7,3);
t173 = cos(qJ(3));
t172 = cos(qJ(4));
t171 = rSges(6,3) + qJ(5);
t169 = t129 * t131;
t168 = t131 * t132;
t166 = t131 * t138;
t163 = t134 * t138;
t161 = t134 * pkin(8) + qJ(1);
t160 = t132 * pkin(1) + pkin(8) * t169;
t157 = t130 * t173;
t156 = t133 * t173;
t155 = t131 * t157;
t135 = sin(qJ(6));
t139 = cos(qJ(6));
t154 = t135 * rSges(7,1) + t139 * rSges(7,2) + qJ(5);
t115 = -t129 * t138 + t132 * t162;
t153 = -t115 * t130 - t132 * t167;
t150 = t139 * rSges(7,1) - t135 * rSges(7,2) + pkin(5) + pkin(10);
t118 = -t129 * t163 + t132 * t140;
t149 = t118 * pkin(2) + t152 * pkin(9) + t160;
t137 = sin(qJ(3));
t102 = t118 * t173 + (t117 * t133 + t130 * t169) * t137;
t148 = t102 * pkin(3) + t149;
t147 = pkin(2) * t166 + t151 * pkin(9) + t161;
t136 = sin(qJ(4));
t95 = t102 * t172 + t152 * t136;
t146 = t95 * pkin(4) + t148;
t109 = t134 * t130 * t137 + (t133 * t137 * t140 + t173 * t138) * t131;
t145 = t109 * pkin(3) + t147;
t104 = t109 * t172 + t151 * t136;
t144 = t104 * pkin(4) + t145;
t116 = t129 * t140 + t132 * t163;
t126 = t129 * pkin(1);
t143 = t116 * pkin(2) - pkin(8) * t168 + t153 * pkin(9) + t126;
t100 = t116 * t173 + (t115 * t133 - t130 * t168) * t137;
t142 = t100 * pkin(3) + t143;
t93 = t100 * t172 + t153 * t136;
t141 = t93 * pkin(4) + t142;
t108 = -t134 * t157 + t137 * t166 - t156 * t165;
t103 = t109 * t136 - t151 * t172;
t101 = -t117 * t156 + t118 * t137 - t129 * t155;
t99 = -t115 * t156 + t116 * t137 + t132 * t155;
t94 = t102 * t136 - t152 * t172;
t92 = t100 * t136 - t153 * t172;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t132 - rSges(2,2) * t129) + g(2) * (rSges(2,1) * t129 + rSges(2,2) * t132) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t118 + rSges(3,2) * t117 + t160) + g(2) * (rSges(3,1) * t116 + rSges(3,2) * t115 + t126) + g(3) * (t134 * rSges(3,3) + t161) + (g(1) * rSges(3,3) * t129 + g(3) * (rSges(3,1) * t138 + rSges(3,2) * t140) + g(2) * (-rSges(3,3) - pkin(8)) * t132) * t131) - m(4) * (g(1) * (t102 * rSges(4,1) - t101 * rSges(4,2) + t152 * rSges(4,3) + t149) + g(2) * (t100 * rSges(4,1) - t99 * rSges(4,2) + t153 * rSges(4,3) + t143) + g(3) * (t109 * rSges(4,1) - t108 * rSges(4,2) + t151 * rSges(4,3) + t147)) - m(5) * (g(1) * (rSges(5,1) * t95 - rSges(5,2) * t94 + t175 * t101 + t148) + g(2) * (rSges(5,1) * t93 - rSges(5,2) * t92 + t175 * t99 + t142) + g(3) * (t104 * rSges(5,1) - t103 * rSges(5,2) + t175 * t108 + t145)) - m(6) * (g(1) * (-rSges(6,2) * t95 + t176 * t101 + t171 * t94 + t146) + g(2) * (-rSges(6,2) * t93 + t171 * t92 + t176 * t99 + t141) + g(3) * (-t104 * rSges(6,2) + t171 * t103 + t176 * t108 + t144)) - m(7) * (g(1) * (t150 * t101 + t154 * t94 + t174 * t95 + t146) + g(2) * (t150 * t99 + t154 * t92 + t174 * t93 + t141) + g(3) * (t154 * t103 + t174 * t104 + t150 * t108 + t144));
U  = t1;

% Calculate potential energy for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:48:56
% EndTime: 2019-03-08 18:48:56
% DurationCPUTime: 0.48s
% Computational Cost: add. (532->118), mult. (1342->161), div. (0->0), fcn. (1691->14), ass. (0->62)
t131 = sin(pkin(7));
t135 = cos(pkin(7));
t136 = cos(pkin(6));
t132 = sin(pkin(6));
t133 = cos(pkin(12));
t166 = t132 * t133;
t151 = -t131 * t166 + t136 * t135;
t129 = sin(pkin(12));
t134 = cos(pkin(11));
t130 = sin(pkin(11));
t167 = t130 * t136;
t117 = -t129 * t134 - t133 * t167;
t164 = t132 * t135;
t152 = -t117 * t131 + t130 * t164;
t176 = rSges(6,1) + pkin(9);
t175 = rSges(5,3) + pkin(9);
t174 = pkin(10) + rSges(7,3);
t173 = cos(qJ(3));
t172 = cos(qJ(4));
t171 = rSges(6,3) + qJ(5);
t169 = t129 * t132;
t168 = t130 * t132;
t165 = t132 * t134;
t163 = t134 * t136;
t161 = t136 * qJ(2) + qJ(1);
t160 = t134 * pkin(1) + qJ(2) * t168;
t157 = t131 * t173;
t156 = t135 * t173;
t155 = t132 * t157;
t137 = sin(qJ(6));
t140 = cos(qJ(6));
t154 = t137 * rSges(7,1) + t140 * rSges(7,2) + qJ(5);
t115 = -t129 * t130 + t133 * t163;
t153 = -t115 * t131 - t134 * t164;
t150 = t140 * rSges(7,1) - t137 * rSges(7,2) + pkin(5) + pkin(9);
t118 = -t129 * t167 + t133 * t134;
t149 = t118 * pkin(2) + t152 * pkin(8) + t160;
t139 = sin(qJ(3));
t102 = t118 * t173 + (t117 * t135 + t131 * t168) * t139;
t148 = t102 * pkin(3) + t149;
t147 = pkin(2) * t169 + t151 * pkin(8) + t161;
t138 = sin(qJ(4));
t95 = t102 * t172 + t152 * t138;
t146 = t95 * pkin(4) + t148;
t109 = t136 * t131 * t139 + (t133 * t135 * t139 + t173 * t129) * t132;
t145 = t109 * pkin(3) + t147;
t104 = t109 * t172 + t151 * t138;
t144 = t104 * pkin(4) + t145;
t116 = t129 * t163 + t130 * t133;
t127 = t130 * pkin(1);
t143 = t116 * pkin(2) + t153 * pkin(8) - qJ(2) * t165 + t127;
t100 = t116 * t173 + (t115 * t135 - t131 * t165) * t139;
t142 = t100 * pkin(3) + t143;
t93 = t100 * t172 + t153 * t138;
t141 = t93 * pkin(4) + t142;
t108 = -t136 * t157 + t139 * t169 - t156 * t166;
t103 = t109 * t138 - t151 * t172;
t101 = -t117 * t156 + t118 * t139 - t130 * t155;
t99 = -t115 * t156 + t116 * t139 + t134 * t155;
t94 = t102 * t138 - t152 * t172;
t92 = t100 * t138 - t153 * t172;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t134 - rSges(2,2) * t130) + g(2) * (rSges(2,1) * t130 + rSges(2,2) * t134) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t118 + rSges(3,2) * t117 + t160) + g(2) * (rSges(3,1) * t116 + rSges(3,2) * t115 + t127) + g(3) * (rSges(3,3) * t136 + t161) + (g(1) * rSges(3,3) * t130 + g(3) * (rSges(3,1) * t129 + rSges(3,2) * t133) + g(2) * (-rSges(3,3) - qJ(2)) * t134) * t132) - m(4) * (g(1) * (t102 * rSges(4,1) - t101 * rSges(4,2) + t152 * rSges(4,3) + t149) + g(2) * (t100 * rSges(4,1) - t99 * rSges(4,2) + t153 * rSges(4,3) + t143) + g(3) * (t109 * rSges(4,1) - t108 * rSges(4,2) + t151 * rSges(4,3) + t147)) - m(5) * (g(1) * (rSges(5,1) * t95 - rSges(5,2) * t94 + t175 * t101 + t148) + g(2) * (rSges(5,1) * t93 - rSges(5,2) * t92 + t175 * t99 + t142) + g(3) * (rSges(5,1) * t104 - rSges(5,2) * t103 + t175 * t108 + t145)) - m(6) * (g(1) * (-rSges(6,2) * t95 + t176 * t101 + t171 * t94 + t146) + g(2) * (-rSges(6,2) * t93 + t171 * t92 + t176 * t99 + t141) + g(3) * (-rSges(6,2) * t104 + t171 * t103 + t176 * t108 + t144)) - m(7) * (g(1) * (t150 * t101 + t154 * t94 + t174 * t95 + t146) + g(2) * (t150 * t99 + t154 * t92 + t174 * t93 + t141) + g(3) * (t154 * t103 + t174 * t104 + t150 * t108 + t144));
U  = t1;

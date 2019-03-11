% Calculate potential energy for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_energypot_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR12_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:47:58
% EndTime: 2019-03-09 05:47:58
% DurationCPUTime: 0.48s
% Computational Cost: add. (532->118), mult. (1342->159), div. (0->0), fcn. (1691->14), ass. (0->62)
t131 = sin(pkin(12));
t136 = cos(pkin(6));
t142 = cos(qJ(1));
t134 = cos(pkin(12));
t140 = sin(qJ(1));
t163 = t140 * t134;
t119 = -t131 * t142 - t136 * t163;
t132 = sin(pkin(7));
t135 = cos(pkin(7));
t133 = sin(pkin(6));
t168 = t133 * t140;
t154 = -t119 * t132 + t135 * t168;
t169 = t133 * t134;
t153 = -t132 * t169 + t136 * t135;
t178 = rSges(6,1) + pkin(10);
t177 = rSges(5,3) + pkin(10);
t176 = pkin(11) + rSges(7,3);
t175 = cos(qJ(3));
t174 = cos(qJ(4));
t173 = t136 * qJ(2) + pkin(8);
t172 = rSges(6,3) + qJ(5);
t170 = t131 * t133;
t167 = t133 * t142;
t165 = t136 * t142;
t164 = t140 * t131;
t162 = t142 * pkin(1) + qJ(2) * t168;
t159 = t132 * t175;
t158 = t135 * t175;
t157 = t133 * t159;
t137 = sin(qJ(6));
t141 = cos(qJ(6));
t156 = t137 * rSges(7,1) + t141 * rSges(7,2) + qJ(5);
t117 = t134 * t165 - t164;
t155 = -t117 * t132 - t135 * t167;
t152 = t141 * rSges(7,1) - t137 * rSges(7,2) + pkin(5) + pkin(10);
t120 = t134 * t142 - t136 * t164;
t151 = t120 * pkin(2) + t154 * pkin(9) + t162;
t150 = pkin(2) * t170 + t153 * pkin(9) + t173;
t139 = sin(qJ(3));
t106 = t120 * t175 + (t119 * t135 + t132 * t168) * t139;
t149 = t106 * pkin(3) + t151;
t111 = t136 * t132 * t139 + (t134 * t135 * t139 + t175 * t131) * t133;
t148 = t111 * pkin(3) + t150;
t138 = sin(qJ(4));
t97 = t106 * t174 + t154 * t138;
t147 = t97 * pkin(4) + t149;
t102 = t111 * t174 + t153 * t138;
t146 = t102 * pkin(4) + t148;
t118 = t131 * t165 + t163;
t129 = t140 * pkin(1);
t145 = t118 * pkin(2) + t155 * pkin(9) - qJ(2) * t167 + t129;
t104 = t118 * t175 + (t117 * t135 - t132 * t167) * t139;
t144 = t104 * pkin(3) + t145;
t95 = t104 * t174 + t155 * t138;
t143 = t95 * pkin(4) + t144;
t110 = -t136 * t159 + t139 * t170 - t158 * t169;
t105 = -t119 * t158 + t120 * t139 - t140 * t157;
t103 = -t117 * t158 + t118 * t139 + t142 * t157;
t101 = t111 * t138 - t153 * t174;
t96 = t106 * t138 - t154 * t174;
t94 = t104 * t138 - t155 * t174;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t142 - t140 * rSges(2,2)) + g(2) * (t140 * rSges(2,1) + rSges(2,2) * t142) + g(3) * (pkin(8) + rSges(2,3))) - m(3) * (g(1) * (t120 * rSges(3,1) + t119 * rSges(3,2) + t162) + g(2) * (t118 * rSges(3,1) + t117 * rSges(3,2) + t129) + g(3) * (rSges(3,3) * t136 + t173) + (g(1) * rSges(3,3) * t140 + g(3) * (rSges(3,1) * t131 + rSges(3,2) * t134) + g(2) * (-rSges(3,3) - qJ(2)) * t142) * t133) - m(4) * (g(1) * (t106 * rSges(4,1) - t105 * rSges(4,2) + t154 * rSges(4,3) + t151) + g(2) * (t104 * rSges(4,1) - t103 * rSges(4,2) + t155 * rSges(4,3) + t145) + g(3) * (t111 * rSges(4,1) - t110 * rSges(4,2) + t153 * rSges(4,3) + t150)) - m(5) * (g(1) * (rSges(5,1) * t97 - rSges(5,2) * t96 + t177 * t105 + t149) + g(2) * (t95 * rSges(5,1) - t94 * rSges(5,2) + t177 * t103 + t144) + g(3) * (rSges(5,1) * t102 - rSges(5,2) * t101 + t177 * t110 + t148)) - m(6) * (g(1) * (-rSges(6,2) * t97 + t178 * t105 + t172 * t96 + t147) + g(2) * (-t95 * rSges(6,2) + t178 * t103 + t172 * t94 + t143) + g(3) * (-rSges(6,2) * t102 + t172 * t101 + t178 * t110 + t146)) - m(7) * (g(1) * (t152 * t105 + t156 * t96 + t176 * t97 + t147) + g(2) * (t152 * t103 + t156 * t94 + t176 * t95 + t143) + g(3) * (t156 * t101 + t176 * t102 + t152 * t110 + t146));
U  = t1;

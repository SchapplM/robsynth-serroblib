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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:05
% EndTime: 2019-03-08 23:48:06
% DurationCPUTime: 0.74s
% Computational Cost: add. (532->108), mult. (1351->136), div. (0->0), fcn. (1691->14), ass. (0->62)
t132 = sin(pkin(7));
t135 = cos(pkin(7));
t136 = cos(pkin(6));
t133 = sin(pkin(6));
t142 = cos(qJ(2));
t166 = t133 * t142;
t152 = -t132 * t166 + t136 * t135;
t131 = sin(pkin(12));
t134 = cos(pkin(12));
t140 = sin(qJ(2));
t163 = t136 * t142;
t119 = -t131 * t163 - t134 * t140;
t168 = t133 * t135;
t153 = -t119 * t132 + t131 * t168;
t180 = -m(7) * pkin(11) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t137 = sin(qJ(6));
t141 = cos(qJ(6));
t179 = -t137 * mrSges(7,1) - t141 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t178 = -m(7) * (pkin(5) + pkin(10)) - t141 * mrSges(7,1) + t137 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t176 = cos(qJ(3));
t175 = cos(qJ(4));
t117 = -t131 * t140 + t134 * t163;
t164 = t136 * t140;
t118 = t131 * t142 + t134 * t164;
t139 = sin(qJ(3));
t158 = t132 * t176;
t156 = t133 * t158;
t157 = t135 * t176;
t101 = -t117 * t157 + t118 * t139 + t134 * t156;
t174 = pkin(10) * t101;
t120 = -t131 * t164 + t134 * t142;
t103 = -t119 * t157 + t120 * t139 - t131 * t156;
t173 = pkin(10) * t103;
t167 = t133 * t140;
t110 = -t136 * t158 + t139 * t167 - t157 * t166;
t172 = t110 * pkin(10);
t170 = t131 * t133;
t169 = t133 * t134;
t162 = t134 * pkin(1) + pkin(8) * t170;
t161 = t136 * pkin(8) + qJ(1);
t155 = t131 * pkin(1) - pkin(8) * t169;
t154 = -t117 * t132 - t134 * t168;
t151 = t120 * pkin(2) + t153 * pkin(9) + t162;
t104 = t120 * t176 + (t119 * t135 + t132 * t170) * t139;
t150 = t104 * pkin(3) + t151;
t149 = pkin(2) * t167 + t152 * pkin(9) + t161;
t111 = t136 * t132 * t139 + (t135 * t139 * t142 + t140 * t176) * t133;
t148 = t111 * pkin(3) + t149;
t138 = sin(qJ(4));
t96 = t104 * t138 - t153 * t175;
t97 = t104 * t175 + t138 * t153;
t147 = t97 * pkin(4) + qJ(5) * t96 + t150;
t105 = t111 * t138 - t152 * t175;
t106 = t111 * t175 + t138 * t152;
t146 = t106 * pkin(4) + t105 * qJ(5) + t148;
t145 = t118 * pkin(2) + pkin(9) * t154 + t155;
t102 = t118 * t176 + (t117 * t135 - t132 * t169) * t139;
t144 = t102 * pkin(3) + t145;
t94 = t102 * t138 - t154 * t175;
t95 = t102 * t175 + t138 * t154;
t143 = t95 * pkin(4) + qJ(5) * t94 + t144;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t161 - t136 * mrSges(3,3) - (t140 * mrSges(3,1) + t142 * mrSges(3,2)) * t133 - m(4) * t149 - t111 * mrSges(4,1) - t152 * mrSges(4,3) - m(5) * (t148 + t172) - m(6) * (t146 + t172) - m(7) * t146 + t179 * t105 + t178 * t110 + t180 * t106) * g(3) + (-mrSges(1,2) - t131 * mrSges(2,1) - t134 * mrSges(2,2) - m(3) * t155 - t118 * mrSges(3,1) - t117 * mrSges(3,2) + mrSges(3,3) * t169 - m(4) * t145 - t102 * mrSges(4,1) - t154 * mrSges(4,3) - m(5) * (t144 + t174) - m(6) * (t143 + t174) - m(7) * t143 + t180 * t95 + t179 * t94 + t178 * t101) * g(2) + (-mrSges(1,1) - t134 * mrSges(2,1) + t131 * mrSges(2,2) - m(3) * t162 - t120 * mrSges(3,1) - t119 * mrSges(3,2) - mrSges(3,3) * t170 - m(4) * t151 - t104 * mrSges(4,1) - t153 * mrSges(4,3) - m(5) * (t150 + t173) - m(6) * (t147 + t173) - m(7) * t147 + t180 * t97 + t179 * t96 + t178 * t103) * g(1);
U  = t1;

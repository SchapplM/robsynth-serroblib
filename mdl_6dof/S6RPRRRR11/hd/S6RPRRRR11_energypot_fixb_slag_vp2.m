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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:11
% EndTime: 2019-03-09 07:38:11
% DurationCPUTime: 0.80s
% Computational Cost: add. (563->105), mult. (1397->135), div. (0->0), fcn. (1753->16), ass. (0->53)
t132 = sin(pkin(13));
t135 = cos(pkin(13));
t143 = cos(qJ(1));
t137 = cos(pkin(6));
t141 = sin(qJ(1));
t166 = t137 * t141;
t116 = -t132 * t143 - t135 * t166;
t133 = sin(pkin(7));
t136 = cos(pkin(7));
t134 = sin(pkin(6));
t169 = t134 * t141;
t155 = -t116 * t133 + t136 * t169;
t170 = t134 * t135;
t154 = -t133 * t170 + t136 * t137;
t179 = -m(5) - m(6);
t178 = -m(6) * pkin(11) + m(7) * (-pkin(12) - pkin(11)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t131 = qJ(5) + qJ(6);
t127 = sin(t131);
t128 = cos(t131);
t138 = sin(qJ(5));
t142 = cos(qJ(5));
t177 = -m(7) * (pkin(5) * t138 + pkin(10)) - t138 * mrSges(6,1) - t127 * mrSges(7,1) - t142 * mrSges(6,2) - t128 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t176 = -m(6) * pkin(4) - m(7) * (pkin(5) * t142 + pkin(4)) - t142 * mrSges(6,1) - t128 * mrSges(7,1) + t138 * mrSges(6,2) + t127 * mrSges(7,2) - mrSges(5,1);
t175 = cos(qJ(3));
t174 = cos(qJ(4));
t173 = t137 * qJ(2) + pkin(8);
t171 = t132 * t134;
t168 = t134 * t143;
t165 = t137 * t143;
t164 = t143 * pkin(1) + qJ(2) * t169;
t160 = t133 * t175;
t159 = t136 * t175;
t158 = t134 * t160;
t157 = t141 * pkin(1) - qJ(2) * t168;
t114 = -t132 * t141 + t135 * t165;
t156 = -t114 * t133 - t136 * t168;
t117 = -t132 * t166 + t135 * t143;
t153 = t117 * pkin(2) + t155 * pkin(9) + t164;
t152 = pkin(2) * t171 + t154 * pkin(9) + t173;
t140 = sin(qJ(3));
t103 = t117 * t175 + (t116 * t136 + t133 * t169) * t140;
t151 = t103 * pkin(3) + t153;
t108 = t137 * t133 * t140 + (t135 * t136 * t140 + t132 * t175) * t134;
t150 = t108 * pkin(3) + t152;
t115 = t132 * t165 + t135 * t141;
t147 = t115 * pkin(2) + pkin(9) * t156 + t157;
t101 = t115 * t175 + (t114 * t136 - t133 * t168) * t140;
t146 = t101 * pkin(3) + t147;
t139 = sin(qJ(4));
t107 = -t137 * t160 + t140 * t171 - t159 * t170;
t102 = -t116 * t159 + t117 * t140 - t141 * t158;
t100 = -t114 * t159 + t115 * t140 + t143 * t158;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t173 - t137 * mrSges(3,3) - (t132 * mrSges(3,1) + t135 * mrSges(3,2)) * t134 - m(4) * t152 - t108 * mrSges(4,1) - t154 * mrSges(4,3) - m(7) * t150 + t179 * (t107 * pkin(10) + t150) + t176 * (t108 * t174 + t139 * t154) + t177 * t107 + t178 * (t108 * t139 - t154 * t174)) * g(3) + (-m(3) * t157 - m(4) * t147 - m(7) * t146 - t141 * mrSges(2,1) - t115 * mrSges(3,1) - t101 * mrSges(4,1) - t143 * mrSges(2,2) - t114 * mrSges(3,2) + mrSges(3,3) * t168 - t156 * mrSges(4,3) - mrSges(1,2) + t179 * (t100 * pkin(10) + t146) + t176 * (t101 * t174 + t139 * t156) + t177 * t100 + t178 * (t101 * t139 - t156 * t174)) * g(2) + (-m(3) * t164 - m(4) * t153 - m(7) * t151 - t143 * mrSges(2,1) - t117 * mrSges(3,1) - t103 * mrSges(4,1) + t141 * mrSges(2,2) - t116 * mrSges(3,2) - mrSges(3,3) * t169 - t155 * mrSges(4,3) - mrSges(1,1) + t179 * (t102 * pkin(10) + t151) + t176 * (t103 * t174 + t139 * t155) + t177 * t102 + t178 * (t103 * t139 - t155 * t174)) * g(1);
U  = t1;

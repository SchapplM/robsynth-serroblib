% Calculate potential energy for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:38:43
% EndTime: 2019-03-08 23:38:43
% DurationCPUTime: 0.81s
% Computational Cost: add. (563->105), mult. (1397->136), div. (0->0), fcn. (1753->16), ass. (0->54)
t134 = sin(pkin(7));
t138 = cos(pkin(7));
t139 = cos(pkin(6));
t135 = sin(pkin(6));
t144 = cos(qJ(2));
t169 = t135 * t144;
t154 = -t134 * t169 + t139 * t138;
t133 = sin(pkin(12));
t137 = cos(pkin(12));
t143 = sin(qJ(2));
t166 = t139 * t144;
t116 = -t133 * t166 - t137 * t143;
t171 = t135 * t138;
t155 = -t116 * t134 + t133 * t171;
t180 = -m(5) - m(6);
t179 = -m(6) * qJ(5) + m(7) * (-pkin(11) - qJ(5)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t131 = pkin(13) + qJ(6);
t126 = sin(t131);
t127 = cos(t131);
t132 = sin(pkin(13));
t136 = cos(pkin(13));
t178 = -m(7) * (pkin(5) * t132 + pkin(10)) - t132 * mrSges(6,1) - t126 * mrSges(7,1) - t136 * mrSges(6,2) - t127 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t177 = -m(6) * pkin(4) - m(7) * (pkin(5) * t136 + pkin(4)) - t136 * mrSges(6,1) - t127 * mrSges(7,1) + t132 * mrSges(6,2) + t126 * mrSges(7,2) - mrSges(5,1);
t176 = cos(qJ(3));
t175 = cos(qJ(4));
t173 = t133 * t135;
t172 = t135 * t137;
t170 = t135 * t143;
t167 = t139 * t143;
t165 = t137 * pkin(1) + pkin(8) * t173;
t164 = t139 * pkin(8) + qJ(1);
t160 = t134 * t176;
t159 = t138 * t176;
t158 = t135 * t160;
t157 = t133 * pkin(1) - pkin(8) * t172;
t114 = -t133 * t143 + t137 * t166;
t156 = -t114 * t134 - t137 * t171;
t117 = -t133 * t167 + t137 * t144;
t153 = t117 * pkin(2) + t155 * pkin(9) + t165;
t142 = sin(qJ(3));
t101 = t117 * t176 + (t116 * t138 + t134 * t173) * t142;
t152 = t101 * pkin(3) + t153;
t151 = pkin(2) * t170 + t154 * pkin(9) + t164;
t108 = t139 * t134 * t142 + (t138 * t142 * t144 + t143 * t176) * t135;
t150 = t108 * pkin(3) + t151;
t115 = t133 * t144 + t137 * t167;
t147 = t115 * pkin(2) + pkin(9) * t156 + t157;
t99 = t115 * t176 + (t114 * t138 - t134 * t172) * t142;
t146 = t99 * pkin(3) + t147;
t141 = sin(qJ(4));
t107 = -t139 * t160 + t142 * t170 - t159 * t169;
t100 = -t116 * t159 + t117 * t142 - t133 * t158;
t98 = -t114 * t159 + t115 * t142 + t137 * t158;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t164 - t139 * mrSges(3,3) - (t143 * mrSges(3,1) + t144 * mrSges(3,2)) * t135 - m(4) * t151 - t108 * mrSges(4,1) - t154 * mrSges(4,3) - m(7) * t150 + t180 * (t107 * pkin(10) + t150) + t177 * (t108 * t175 + t141 * t154) + t178 * t107 + t179 * (t108 * t141 - t154 * t175)) * g(3) + (-m(3) * t157 - m(4) * t147 - m(7) * t146 - mrSges(2,1) * t133 - t115 * mrSges(3,1) - t99 * mrSges(4,1) - mrSges(2,2) * t137 - t114 * mrSges(3,2) + mrSges(3,3) * t172 - t156 * mrSges(4,3) - mrSges(1,2) + t180 * (pkin(10) * t98 + t146) + t177 * (t141 * t156 + t175 * t99) + t178 * t98 + t179 * (t141 * t99 - t156 * t175)) * g(2) + (-m(3) * t165 - m(4) * t153 - m(7) * t152 - mrSges(2,1) * t137 - t117 * mrSges(3,1) - t101 * mrSges(4,1) + mrSges(2,2) * t133 - t116 * mrSges(3,2) - mrSges(3,3) * t173 - t155 * mrSges(4,3) - mrSges(1,1) + t180 * (pkin(10) * t100 + t152) + t177 * (t101 * t175 + t141 * t155) + t178 * t100 + t179 * (t101 * t141 - t155 * t175)) * g(1);
U  = t1;

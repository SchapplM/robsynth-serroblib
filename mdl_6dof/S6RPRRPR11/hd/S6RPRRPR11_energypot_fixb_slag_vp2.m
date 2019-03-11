% Calculate potential energy for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:39:09
% EndTime: 2019-03-09 05:39:10
% DurationCPUTime: 0.73s
% Computational Cost: add. (563->105), mult. (1397->134), div. (0->0), fcn. (1753->16), ass. (0->54)
t134 = sin(pkin(12));
t140 = cos(pkin(6));
t145 = cos(qJ(1));
t138 = cos(pkin(12));
t144 = sin(qJ(1));
t166 = t144 * t138;
t117 = -t134 * t145 - t140 * t166;
t135 = sin(pkin(7));
t139 = cos(pkin(7));
t136 = sin(pkin(6));
t171 = t136 * t144;
t156 = -t117 * t135 + t139 * t171;
t172 = t136 * t138;
t155 = -t135 * t172 + t139 * t140;
t181 = -m(5) - m(6);
t180 = -m(6) * qJ(5) + m(7) * (-pkin(11) - qJ(5)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t132 = pkin(13) + qJ(6);
t127 = sin(t132);
t128 = cos(t132);
t133 = sin(pkin(13));
t137 = cos(pkin(13));
t179 = -m(7) * (pkin(5) * t133 + pkin(10)) - t133 * mrSges(6,1) - t127 * mrSges(7,1) - t137 * mrSges(6,2) - t128 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t178 = -m(6) * pkin(4) - m(7) * (pkin(5) * t137 + pkin(4)) - t137 * mrSges(6,1) - t128 * mrSges(7,1) + t133 * mrSges(6,2) + t127 * mrSges(7,2) - mrSges(5,1);
t177 = cos(qJ(3));
t176 = cos(qJ(4));
t175 = t140 * qJ(2) + pkin(8);
t173 = t136 * t134;
t170 = t136 * t145;
t168 = t140 * t145;
t167 = t144 * t134;
t165 = t145 * pkin(1) + qJ(2) * t171;
t161 = t135 * t177;
t160 = t139 * t177;
t159 = t136 * t161;
t158 = t144 * pkin(1) - qJ(2) * t170;
t115 = t138 * t168 - t167;
t157 = -t115 * t135 - t139 * t170;
t118 = t138 * t145 - t140 * t167;
t154 = t118 * pkin(2) + t156 * pkin(9) + t165;
t153 = pkin(2) * t173 + t155 * pkin(9) + t175;
t143 = sin(qJ(3));
t104 = t118 * t177 + (t117 * t139 + t135 * t171) * t143;
t152 = t104 * pkin(3) + t154;
t109 = t140 * t135 * t143 + (t138 * t139 * t143 + t177 * t134) * t136;
t151 = t109 * pkin(3) + t153;
t116 = t134 * t168 + t166;
t148 = t116 * pkin(2) + t157 * pkin(9) + t158;
t102 = t116 * t177 + (t115 * t139 - t135 * t170) * t143;
t147 = t102 * pkin(3) + t148;
t142 = sin(qJ(4));
t108 = -t140 * t161 + t143 * t173 - t160 * t172;
t103 = -t117 * t160 + t118 * t143 - t144 * t159;
t101 = -t115 * t160 + t116 * t143 + t145 * t159;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t175 - t140 * mrSges(3,3) - (t134 * mrSges(3,1) + t138 * mrSges(3,2)) * t136 - m(4) * t153 - t109 * mrSges(4,1) - t155 * mrSges(4,3) - m(7) * t151 + t181 * (pkin(10) * t108 + t151) + t178 * (t109 * t176 + t155 * t142) + t179 * t108 + t180 * (t109 * t142 - t155 * t176)) * g(3) + (-m(3) * t158 - m(4) * t148 - m(7) * t147 - t144 * mrSges(2,1) - t116 * mrSges(3,1) - t102 * mrSges(4,1) - t145 * mrSges(2,2) - t115 * mrSges(3,2) + mrSges(3,3) * t170 - t157 * mrSges(4,3) - mrSges(1,2) + t181 * (t101 * pkin(10) + t147) + t178 * (t102 * t176 + t157 * t142) + t179 * t101 + t180 * (t102 * t142 - t157 * t176)) * g(2) + (-m(3) * t165 - m(4) * t154 - m(7) * t152 - t145 * mrSges(2,1) - t118 * mrSges(3,1) - t104 * mrSges(4,1) + t144 * mrSges(2,2) - t117 * mrSges(3,2) - mrSges(3,3) * t171 - t156 * mrSges(4,3) - mrSges(1,1) + t181 * (pkin(10) * t103 + t152) + t178 * (t104 * t176 + t156 * t142) + t179 * t103 + t180 * (t104 * t142 - t156 * t176)) * g(1);
U  = t1;

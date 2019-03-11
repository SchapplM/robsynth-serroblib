% Calculate potential energy for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRPR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_energypot_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:12
% EndTime: 2019-03-08 18:45:13
% DurationCPUTime: 0.79s
% Computational Cost: add. (563->105), mult. (1397->136), div. (0->0), fcn. (1753->16), ass. (0->54)
t136 = sin(pkin(7));
t141 = cos(pkin(7));
t142 = cos(pkin(6));
t137 = sin(pkin(6));
t139 = cos(pkin(12));
t171 = t137 * t139;
t155 = -t136 * t171 + t141 * t142;
t134 = sin(pkin(12));
t140 = cos(pkin(11));
t135 = sin(pkin(11));
t172 = t135 * t142;
t117 = -t134 * t140 - t139 * t172;
t169 = t137 * t141;
t156 = -t117 * t136 + t135 * t169;
t181 = -m(5) - m(6);
t180 = -m(6) * qJ(5) + m(7) * (-pkin(10) - qJ(5)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t132 = pkin(13) + qJ(6);
t127 = sin(t132);
t128 = cos(t132);
t133 = sin(pkin(13));
t138 = cos(pkin(13));
t179 = -m(7) * (pkin(5) * t133 + pkin(9)) - t133 * mrSges(6,1) - t127 * mrSges(7,1) - t138 * mrSges(6,2) - t128 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t178 = -m(6) * pkin(4) - m(7) * (pkin(5) * t138 + pkin(4)) - t138 * mrSges(6,1) - t128 * mrSges(7,1) + t133 * mrSges(6,2) + t127 * mrSges(7,2) - mrSges(5,1);
t177 = cos(qJ(3));
t176 = cos(qJ(4));
t174 = t134 * t137;
t173 = t135 * t137;
t170 = t137 * t140;
t168 = t140 * t142;
t166 = t142 * qJ(2) + qJ(1);
t165 = t140 * pkin(1) + qJ(2) * t173;
t161 = t136 * t177;
t160 = t141 * t177;
t159 = t137 * t161;
t158 = t135 * pkin(1) - qJ(2) * t170;
t115 = -t134 * t135 + t139 * t168;
t157 = -t115 * t136 - t140 * t169;
t118 = -t134 * t172 + t139 * t140;
t154 = t118 * pkin(2) + t156 * pkin(8) + t165;
t145 = sin(qJ(3));
t102 = t118 * t177 + (t117 * t141 + t136 * t173) * t145;
t153 = t102 * pkin(3) + t154;
t152 = pkin(2) * t174 + t155 * pkin(8) + t166;
t109 = t142 * t136 * t145 + (t139 * t141 * t145 + t177 * t134) * t137;
t151 = t109 * pkin(3) + t152;
t116 = t134 * t168 + t135 * t139;
t148 = t116 * pkin(2) + t157 * pkin(8) + t158;
t100 = t116 * t177 + (t115 * t141 - t136 * t170) * t145;
t147 = t100 * pkin(3) + t148;
t144 = sin(qJ(4));
t108 = -t142 * t161 + t145 * t174 - t160 * t171;
t101 = -t117 * t160 + t118 * t145 - t135 * t159;
t99 = -t115 * t160 + t116 * t145 + t140 * t159;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t166 - t142 * mrSges(3,3) - (t134 * mrSges(3,1) + t139 * mrSges(3,2)) * t137 - m(4) * t152 - t109 * mrSges(4,1) - t155 * mrSges(4,3) - m(7) * t151 + t181 * (pkin(9) * t108 + t151) + t178 * (t109 * t176 + t155 * t144) + t179 * t108 + t180 * (t109 * t144 - t155 * t176)) * g(3) + (-m(3) * t158 - m(4) * t148 - m(7) * t147 - t135 * mrSges(2,1) - t116 * mrSges(3,1) - t100 * mrSges(4,1) - t140 * mrSges(2,2) - t115 * mrSges(3,2) + mrSges(3,3) * t170 - t157 * mrSges(4,3) - mrSges(1,2) + t181 * (pkin(9) * t99 + t147) + t178 * (t100 * t176 + t157 * t144) + t179 * t99 + t180 * (t100 * t144 - t157 * t176)) * g(2) + (-m(3) * t165 - m(4) * t154 - m(7) * t153 - t140 * mrSges(2,1) - t118 * mrSges(3,1) - t102 * mrSges(4,1) + t135 * mrSges(2,2) - t117 * mrSges(3,2) - mrSges(3,3) * t173 - t156 * mrSges(4,3) - mrSges(1,1) + t181 * (pkin(9) * t101 + t153) + t178 * (t102 * t176 + t156 * t144) + t179 * t101 + t180 * (t102 * t144 - t156 * t176)) * g(1);
U  = t1;

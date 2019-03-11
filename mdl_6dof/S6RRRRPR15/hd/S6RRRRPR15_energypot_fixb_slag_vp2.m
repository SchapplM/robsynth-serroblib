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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR15_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:06
% EndTime: 2019-03-10 00:34:06
% DurationCPUTime: 0.75s
% Computational Cost: add. (532->108), mult. (1351->133), div. (0->0), fcn. (1691->14), ass. (0->63)
t135 = cos(pkin(6));
t140 = sin(qJ(1));
t142 = cos(qJ(2));
t164 = t140 * t142;
t139 = sin(qJ(2));
t143 = cos(qJ(1));
t166 = t139 * t143;
t120 = -t135 * t164 - t166;
t132 = sin(pkin(7));
t134 = cos(pkin(7));
t133 = sin(pkin(6));
t170 = t133 * t140;
t154 = -t120 * t132 + t134 * t170;
t169 = t133 * t142;
t153 = -t132 * t169 + t134 * t135;
t182 = -m(7) * pkin(12) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t136 = sin(qJ(6));
t141 = cos(qJ(6));
t181 = -t136 * mrSges(7,1) - t141 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t180 = -m(7) * (pkin(5) + pkin(11)) - t141 * mrSges(7,1) + t136 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t178 = cos(qJ(3));
t177 = cos(qJ(4));
t163 = t142 * t143;
t165 = t140 * t139;
t121 = -t135 * t165 + t163;
t138 = sin(qJ(3));
t159 = t132 * t178;
t157 = t133 * t159;
t158 = t134 * t178;
t106 = -t120 * t158 + t121 * t138 - t140 * t157;
t176 = pkin(11) * t106;
t171 = t133 * t139;
t111 = -t135 * t159 + t138 * t171 - t158 * t169;
t175 = pkin(11) * t111;
t118 = t135 * t163 - t165;
t119 = t135 * t166 + t164;
t104 = -t118 * t158 + t119 * t138 + t143 * t157;
t174 = t104 * pkin(11);
t173 = t135 * pkin(9) + pkin(8);
t168 = t133 * t143;
t162 = t143 * pkin(1) + pkin(9) * t170;
t156 = t140 * pkin(1) - pkin(9) * t168;
t155 = -t118 * t132 - t134 * t168;
t152 = t121 * pkin(2) + t154 * pkin(10) + t162;
t151 = pkin(2) * t171 + t153 * pkin(10) + t173;
t107 = t121 * t178 + (t120 * t134 + t132 * t170) * t138;
t150 = t107 * pkin(3) + t152;
t112 = t135 * t132 * t138 + (t134 * t138 * t142 + t139 * t178) * t133;
t149 = t112 * pkin(3) + t151;
t137 = sin(qJ(4));
t97 = t107 * t137 - t154 * t177;
t98 = t107 * t177 + t137 * t154;
t148 = t98 * pkin(4) + qJ(5) * t97 + t150;
t102 = t112 * t137 - t153 * t177;
t103 = t112 * t177 + t137 * t153;
t147 = t103 * pkin(4) + qJ(5) * t102 + t149;
t146 = t119 * pkin(2) + pkin(10) * t155 + t156;
t105 = t119 * t178 + (t118 * t134 - t132 * t168) * t138;
t145 = t105 * pkin(3) + t146;
t95 = t105 * t137 - t155 * t177;
t96 = t105 * t177 + t137 * t155;
t144 = t96 * pkin(4) + t95 * qJ(5) + t145;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t173 - t135 * mrSges(3,3) - (t139 * mrSges(3,1) + t142 * mrSges(3,2)) * t133 - m(4) * t151 - t112 * mrSges(4,1) - t153 * mrSges(4,3) - m(5) * (t149 + t175) - m(6) * (t147 + t175) - m(7) * t147 + t181 * t102 + t180 * t111 + t182 * t103) * g(3) + (-mrSges(1,2) - t140 * mrSges(2,1) - t143 * mrSges(2,2) - m(3) * t156 - t119 * mrSges(3,1) - t118 * mrSges(3,2) + mrSges(3,3) * t168 - m(4) * t146 - t105 * mrSges(4,1) - t155 * mrSges(4,3) - m(5) * (t145 + t174) - m(6) * (t144 + t174) - m(7) * t144 + t182 * t96 + t181 * t95 + t180 * t104) * g(2) + (-mrSges(1,1) - t143 * mrSges(2,1) + t140 * mrSges(2,2) - m(3) * t162 - t121 * mrSges(3,1) - t120 * mrSges(3,2) - mrSges(3,3) * t170 - m(4) * t152 - t107 * mrSges(4,1) - t154 * mrSges(4,3) - m(5) * (t150 + t176) - m(6) * (t148 + t176) - m(7) * t148 + t182 * t98 + t181 * t97 + t180 * t106) * g(1);
U  = t1;

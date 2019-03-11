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
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:48:08
% EndTime: 2019-03-09 05:48:09
% DurationCPUTime: 0.71s
% Computational Cost: add. (532->108), mult. (1351->134), div. (0->0), fcn. (1691->14), ass. (0->62)
t132 = sin(pkin(12));
t137 = cos(pkin(6));
t143 = cos(qJ(1));
t135 = cos(pkin(12));
t141 = sin(qJ(1));
t163 = t141 * t135;
t120 = -t132 * t143 - t137 * t163;
t133 = sin(pkin(7));
t136 = cos(pkin(7));
t134 = sin(pkin(6));
t168 = t134 * t141;
t154 = -t120 * t133 + t136 * t168;
t169 = t134 * t135;
t153 = -t133 * t169 + t136 * t137;
t181 = -m(7) * pkin(11) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t138 = sin(qJ(6));
t142 = cos(qJ(6));
t180 = -t138 * mrSges(7,1) - t142 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t179 = -m(7) * (pkin(5) + pkin(10)) - t142 * mrSges(7,1) + t138 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t177 = cos(qJ(3));
t176 = cos(qJ(4));
t164 = t141 * t132;
t121 = t135 * t143 - t137 * t164;
t140 = sin(qJ(3));
t159 = t133 * t177;
t157 = t134 * t159;
t158 = t136 * t177;
t106 = -t120 * t158 + t121 * t140 - t141 * t157;
t175 = pkin(10) * t106;
t170 = t132 * t134;
t111 = -t137 * t159 + t140 * t170 - t158 * t169;
t174 = pkin(10) * t111;
t165 = t137 * t143;
t118 = t135 * t165 - t164;
t119 = t132 * t165 + t163;
t104 = -t118 * t158 + t119 * t140 + t143 * t157;
t173 = t104 * pkin(10);
t172 = t137 * qJ(2) + pkin(8);
t167 = t134 * t143;
t162 = t143 * pkin(1) + qJ(2) * t168;
t156 = t141 * pkin(1) - qJ(2) * t167;
t155 = -t118 * t133 - t136 * t167;
t152 = t121 * pkin(2) + t154 * pkin(9) + t162;
t151 = pkin(2) * t170 + t153 * pkin(9) + t172;
t107 = t121 * t177 + (t120 * t136 + t133 * t168) * t140;
t150 = t107 * pkin(3) + t152;
t112 = t137 * t133 * t140 + (t135 * t136 * t140 + t132 * t177) * t134;
t149 = t112 * pkin(3) + t151;
t139 = sin(qJ(4));
t97 = t107 * t139 - t154 * t176;
t98 = t107 * t176 + t139 * t154;
t148 = t98 * pkin(4) + qJ(5) * t97 + t150;
t102 = t112 * t139 - t153 * t176;
t103 = t112 * t176 + t139 * t153;
t147 = t103 * pkin(4) + qJ(5) * t102 + t149;
t146 = t119 * pkin(2) + pkin(9) * t155 + t156;
t105 = t119 * t177 + (t118 * t136 - t133 * t167) * t140;
t145 = t105 * pkin(3) + t146;
t95 = t105 * t139 - t155 * t176;
t96 = t105 * t176 + t139 * t155;
t144 = t96 * pkin(4) + t95 * qJ(5) + t145;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t172 - t137 * mrSges(3,3) - (t132 * mrSges(3,1) + t135 * mrSges(3,2)) * t134 - m(4) * t151 - t112 * mrSges(4,1) - t153 * mrSges(4,3) - m(5) * (t149 + t174) - m(6) * (t147 + t174) - m(7) * t147 + t180 * t102 + t179 * t111 + t181 * t103) * g(3) + (-mrSges(1,2) - t141 * mrSges(2,1) - t143 * mrSges(2,2) - m(3) * t156 - t119 * mrSges(3,1) - t118 * mrSges(3,2) + mrSges(3,3) * t167 - m(4) * t146 - t105 * mrSges(4,1) - t155 * mrSges(4,3) - m(5) * (t145 + t173) - m(6) * (t144 + t173) - m(7) * t144 + t181 * t96 + t180 * t95 + t179 * t104) * g(2) + (-mrSges(1,1) - t143 * mrSges(2,1) + t141 * mrSges(2,2) - m(3) * t162 - t121 * mrSges(3,1) - t120 * mrSges(3,2) - mrSges(3,3) * t168 - m(4) * t152 - t107 * mrSges(4,1) - t154 * mrSges(4,3) - m(5) * (t150 + t175) - m(6) * (t148 + t175) - m(7) * t148 + t181 * t98 + t180 * t97 + t179 * t106) * g(1);
U  = t1;

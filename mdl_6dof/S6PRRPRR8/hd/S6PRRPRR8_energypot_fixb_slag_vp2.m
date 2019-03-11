% Calculate potential energy for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:35:59
% EndTime: 2019-03-08 22:36:00
% DurationCPUTime: 0.74s
% Computational Cost: add. (477->102), mult. (1198->132), div. (0->0), fcn. (1483->14), ass. (0->55)
t138 = sin(pkin(7));
t141 = cos(pkin(7));
t142 = cos(pkin(6));
t139 = sin(pkin(6));
t149 = cos(qJ(2));
t171 = t139 * t149;
t121 = -t138 * t171 + t142 * t141;
t137 = sin(pkin(12));
t140 = cos(pkin(12));
t146 = sin(qJ(2));
t167 = t142 * t149;
t124 = -t137 * t167 - t140 * t146;
t173 = t139 * t141;
t115 = -t124 * t138 + t137 * t173;
t184 = -m(6) - m(7);
t183 = mrSges(4,2) - mrSges(5,3);
t182 = -mrSges(4,3) - mrSges(5,1);
t181 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t143 = sin(qJ(6));
t147 = cos(qJ(6));
t180 = -t143 * mrSges(7,1) - t147 * mrSges(7,2) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t179 = -m(7) * pkin(5) - t147 * mrSges(7,1) + t143 * mrSges(7,2) - mrSges(6,1);
t178 = cos(qJ(3));
t176 = t137 * t139;
t145 = sin(qJ(3));
t175 = t138 * t145;
t174 = t139 * t140;
t172 = t139 * t146;
t170 = t141 * t145;
t168 = t142 * t146;
t166 = t140 * pkin(1) + pkin(8) * t176;
t165 = t142 * pkin(8) + qJ(1);
t162 = t138 * t178;
t161 = t141 * t178;
t160 = t139 * t162;
t159 = t137 * pkin(1) - pkin(8) * t174;
t122 = -t137 * t146 + t140 * t167;
t114 = -t122 * t138 - t140 * t173;
t125 = -t137 * t168 + t140 * t149;
t158 = t125 * pkin(2) + t115 * pkin(9) + t166;
t157 = pkin(2) * t172 + t121 * pkin(9) + t165;
t104 = -t124 * t161 + t125 * t145 - t137 * t160;
t105 = t125 * t178 + (t124 * t141 + t138 * t176) * t145;
t156 = t105 * pkin(3) + qJ(4) * t104 + t158;
t112 = -t142 * t162 + t145 * t172 - t161 * t171;
t113 = t142 * t175 + (t178 * t146 + t149 * t170) * t139;
t155 = t113 * pkin(3) + t112 * qJ(4) + t157;
t123 = t137 * t149 + t140 * t168;
t154 = t123 * pkin(2) + t114 * pkin(9) + t159;
t102 = -t122 * t161 + t123 * t145 + t140 * t160;
t103 = t122 * t170 + t123 * t178 - t174 * t175;
t151 = t103 * pkin(3) + qJ(4) * t102 + t154;
t148 = cos(qJ(5));
t144 = sin(qJ(5));
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t165 - t142 * mrSges(3,3) - (t146 * mrSges(3,1) + t149 * mrSges(3,2)) * t139 - m(4) * t157 - m(5) * t155 + t184 * (t121 * pkin(4) + t113 * pkin(10) + t155) + t182 * t121 + t183 * t112 + t181 * (-t112 * t148 + t121 * t144) + t179 * (t112 * t144 + t121 * t148) + t180 * t113) * g(3) + (-m(3) * t159 - m(4) * t154 - m(5) * t151 - t137 * mrSges(2,1) - t123 * mrSges(3,1) - t140 * mrSges(2,2) - t122 * mrSges(3,2) + mrSges(3,3) * t174 - mrSges(1,2) + t184 * (t114 * pkin(4) + pkin(10) * t103 + t151) + t181 * (-t102 * t148 + t114 * t144) + t182 * t114 + t183 * t102 + t179 * (t102 * t144 + t114 * t148) + t180 * t103) * g(2) + (-m(3) * t166 - m(4) * t158 - m(5) * t156 - t140 * mrSges(2,1) - t125 * mrSges(3,1) + t137 * mrSges(2,2) - t124 * mrSges(3,2) - mrSges(3,3) * t176 - mrSges(1,1) + t184 * (t115 * pkin(4) + pkin(10) * t105 + t156) + t181 * (-t104 * t148 + t115 * t144) + t182 * t115 + t183 * t104 + t179 * (t104 * t144 + t115 * t148) + t180 * t105) * g(1);
U  = t1;

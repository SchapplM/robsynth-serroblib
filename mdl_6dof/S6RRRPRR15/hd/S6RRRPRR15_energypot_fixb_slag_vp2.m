% Calculate potential energy for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR15_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:26
% EndTime: 2019-03-09 20:24:27
% DurationCPUTime: 0.75s
% Computational Cost: add. (477->102), mult. (1198->129), div. (0->0), fcn. (1483->14), ass. (0->56)
t141 = cos(pkin(6));
t146 = sin(qJ(1));
t149 = cos(qJ(2));
t168 = t146 * t149;
t145 = sin(qJ(2));
t150 = cos(qJ(1));
t170 = t145 * t150;
t125 = -t141 * t168 - t170;
t138 = sin(pkin(7));
t140 = cos(pkin(7));
t139 = sin(pkin(6));
t175 = t139 * t146;
t116 = -t125 * t138 + t140 * t175;
t174 = t139 * t149;
t122 = -t138 * t174 + t140 * t141;
t186 = -m(6) - m(7);
t185 = mrSges(4,2) - mrSges(5,3);
t184 = -mrSges(4,3) - mrSges(5,1);
t183 = -m(7) * pkin(12) + mrSges(6,2) - mrSges(7,3);
t142 = sin(qJ(6));
t147 = cos(qJ(6));
t182 = -t142 * mrSges(7,1) - t147 * mrSges(7,2) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t181 = -m(7) * pkin(5) - t147 * mrSges(7,1) + t142 * mrSges(7,2) - mrSges(6,1);
t180 = cos(qJ(3));
t179 = t141 * pkin(9) + pkin(8);
t144 = sin(qJ(3));
t177 = t138 * t144;
t176 = t139 * t145;
t173 = t139 * t150;
t171 = t140 * t144;
t169 = t146 * t145;
t167 = t149 * t150;
t166 = t150 * pkin(1) + pkin(9) * t175;
t163 = t138 * t180;
t162 = t140 * t180;
t161 = t139 * t163;
t160 = t146 * pkin(1) - pkin(9) * t173;
t123 = t141 * t167 - t169;
t115 = -t123 * t138 - t140 * t173;
t126 = -t141 * t169 + t167;
t159 = t126 * pkin(2) + t116 * pkin(10) + t166;
t158 = pkin(2) * t176 + t122 * pkin(10) + t179;
t107 = -t125 * t162 + t126 * t144 - t146 * t161;
t108 = t126 * t180 + (t125 * t140 + t138 * t175) * t144;
t157 = t108 * pkin(3) + qJ(4) * t107 + t159;
t111 = -t141 * t163 + t144 * t176 - t162 * t174;
t112 = t141 * t177 + (t180 * t145 + t149 * t171) * t139;
t156 = t112 * pkin(3) + qJ(4) * t111 + t158;
t124 = t141 * t170 + t168;
t155 = t124 * pkin(2) + t115 * pkin(10) + t160;
t105 = -t123 * t162 + t124 * t144 + t150 * t161;
t106 = t123 * t171 + t124 * t180 - t173 * t177;
t152 = t106 * pkin(3) + t105 * qJ(4) + t155;
t148 = cos(qJ(5));
t143 = sin(qJ(5));
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t179 - t141 * mrSges(3,3) - (t145 * mrSges(3,1) + t149 * mrSges(3,2)) * t139 - m(4) * t158 - m(5) * t156 + t186 * (t122 * pkin(4) + pkin(11) * t112 + t156) + t184 * t122 + t185 * t111 + t183 * (-t111 * t148 + t122 * t143) + t181 * (t111 * t143 + t122 * t148) + t182 * t112) * g(3) + (-m(3) * t160 - m(4) * t155 - m(5) * t152 - t146 * mrSges(2,1) - t124 * mrSges(3,1) - t150 * mrSges(2,2) - t123 * mrSges(3,2) + mrSges(3,3) * t173 - mrSges(1,2) + t186 * (t115 * pkin(4) + t106 * pkin(11) + t152) + t183 * (-t105 * t148 + t115 * t143) + t184 * t115 + t185 * t105 + t181 * (t105 * t143 + t115 * t148) + t182 * t106) * g(2) + (-m(3) * t166 - m(4) * t159 - m(5) * t157 - t150 * mrSges(2,1) - t126 * mrSges(3,1) + t146 * mrSges(2,2) - t125 * mrSges(3,2) - mrSges(3,3) * t175 - mrSges(1,1) + t186 * (t116 * pkin(4) + pkin(11) * t108 + t157) + t183 * (-t107 * t148 + t116 * t143) + t184 * t116 + t185 * t107 + t181 * (t107 * t143 + t116 * t148) + t182 * t108) * g(1);
U  = t1;

% Calculate potential energy for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR13_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:31
% EndTime: 2019-03-09 04:21:32
% DurationCPUTime: 0.75s
% Computational Cost: add. (477->102), mult. (1198->130), div. (0->0), fcn. (1483->14), ass. (0->55)
t138 = sin(pkin(12));
t143 = cos(pkin(6));
t150 = cos(qJ(1));
t141 = cos(pkin(12));
t147 = sin(qJ(1));
t167 = t147 * t141;
t125 = -t138 * t150 - t143 * t167;
t139 = sin(pkin(7));
t142 = cos(pkin(7));
t140 = sin(pkin(6));
t173 = t140 * t147;
t116 = -t125 * t139 + t142 * t173;
t174 = t140 * t141;
t122 = -t139 * t174 + t142 * t143;
t185 = -m(6) - m(7);
t184 = mrSges(4,2) - mrSges(5,3);
t183 = -mrSges(4,3) - mrSges(5,1);
t182 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t144 = sin(qJ(6));
t148 = cos(qJ(6));
t181 = -t144 * mrSges(7,1) - t148 * mrSges(7,2) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t180 = -m(7) * pkin(5) - t148 * mrSges(7,1) + t144 * mrSges(7,2) - mrSges(6,1);
t179 = cos(qJ(3));
t178 = t143 * qJ(2) + pkin(8);
t176 = t138 * t140;
t146 = sin(qJ(3));
t175 = t139 * t146;
t172 = t140 * t150;
t170 = t142 * t146;
t169 = t143 * t150;
t168 = t147 * t138;
t166 = t150 * pkin(1) + qJ(2) * t173;
t163 = t139 * t179;
t162 = t142 * t179;
t161 = t140 * t163;
t160 = t147 * pkin(1) - qJ(2) * t172;
t123 = t141 * t169 - t168;
t115 = -t123 * t139 - t142 * t172;
t126 = t141 * t150 - t143 * t168;
t159 = t126 * pkin(2) + t116 * pkin(9) + t166;
t158 = pkin(2) * t176 + t122 * pkin(9) + t178;
t107 = -t125 * t162 + t126 * t146 - t147 * t161;
t108 = t126 * t179 + (t125 * t142 + t139 * t173) * t146;
t157 = t108 * pkin(3) + qJ(4) * t107 + t159;
t111 = -t143 * t163 + t146 * t176 - t162 * t174;
t112 = t143 * t175 + (t179 * t138 + t141 * t170) * t140;
t156 = t112 * pkin(3) + qJ(4) * t111 + t158;
t124 = t138 * t169 + t167;
t155 = t124 * pkin(2) + t115 * pkin(9) + t160;
t105 = -t123 * t162 + t124 * t146 + t150 * t161;
t106 = t123 * t170 + t124 * t179 - t172 * t175;
t152 = t106 * pkin(3) + t105 * qJ(4) + t155;
t149 = cos(qJ(5));
t145 = sin(qJ(5));
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t178 - t143 * mrSges(3,3) - (t138 * mrSges(3,1) + t141 * mrSges(3,2)) * t140 - m(4) * t158 - m(5) * t156 + t185 * (t122 * pkin(4) + pkin(10) * t112 + t156) + t183 * t122 + t184 * t111 + t182 * (-t111 * t149 + t122 * t145) + t180 * (t111 * t145 + t122 * t149) + t181 * t112) * g(3) + (-m(3) * t160 - m(4) * t155 - m(5) * t152 - t147 * mrSges(2,1) - t124 * mrSges(3,1) - t150 * mrSges(2,2) - t123 * mrSges(3,2) + mrSges(3,3) * t172 - mrSges(1,2) + t185 * (t115 * pkin(4) + t106 * pkin(10) + t152) + t182 * (-t105 * t149 + t115 * t145) + t183 * t115 + t184 * t105 + t180 * (t105 * t145 + t115 * t149) + t181 * t106) * g(2) + (-m(3) * t166 - m(4) * t159 - m(5) * t157 - t150 * mrSges(2,1) - t126 * mrSges(3,1) + t147 * mrSges(2,2) - t125 * mrSges(3,2) - mrSges(3,3) * t173 - mrSges(1,1) + t185 * (t116 * pkin(4) + pkin(10) * t108 + t157) + t182 * (-t107 * t149 + t116 * t145) + t183 * t116 + t184 * t107 + t180 * (t107 * t145 + t116 * t149) + t181 * t108) * g(1);
U  = t1;

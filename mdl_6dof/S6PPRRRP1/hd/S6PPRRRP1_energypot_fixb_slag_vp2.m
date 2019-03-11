% Calculate potential energy for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:14
% EndTime: 2019-03-08 18:52:15
% DurationCPUTime: 0.77s
% Computational Cost: add. (551->116), mult. (1397->153), div. (0->0), fcn. (1753->14), ass. (0->58)
t141 = sin(pkin(7));
t145 = cos(pkin(7));
t146 = cos(pkin(6));
t142 = sin(pkin(6));
t143 = cos(pkin(12));
t177 = t142 * t143;
t161 = -t141 * t177 + t145 * t146;
t139 = sin(pkin(12));
t144 = cos(pkin(11));
t140 = sin(pkin(11));
t179 = t140 * t146;
t126 = -t139 * t144 - t143 * t179;
t175 = t142 * t145;
t162 = -t126 * t141 + t140 * t175;
t187 = -mrSges(6,1) - mrSges(7,1);
t186 = -mrSges(6,2) - mrSges(7,2);
t148 = sin(qJ(5));
t185 = -m(7) * (pkin(5) * t148 + pkin(9)) + mrSges(4,2) - mrSges(5,3);
t184 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t183 = cos(qJ(3));
t182 = cos(qJ(4));
t180 = t140 * t142;
t178 = t142 * t139;
t176 = t142 * t144;
t174 = t144 * t146;
t172 = t146 * qJ(2) + qJ(1);
t171 = t144 * pkin(1) + qJ(2) * t180;
t167 = t141 * t183;
t166 = t145 * t183;
t165 = t142 * t167;
t164 = t140 * pkin(1) - qJ(2) * t176;
t124 = -t139 * t140 + t143 * t174;
t163 = -t124 * t141 - t144 * t175;
t127 = -t139 * t179 + t143 * t144;
t160 = t127 * pkin(2) + t162 * pkin(8) + t171;
t159 = pkin(2) * t178 + t161 * pkin(8) + t172;
t150 = sin(qJ(3));
t111 = t127 * t183 + (t126 * t145 + t141 * t180) * t150;
t158 = t111 * pkin(3) + t160;
t118 = t146 * t141 * t150 + (t143 * t145 * t150 + t139 * t183) * t142;
t157 = t118 * pkin(3) + t159;
t110 = -t126 * t166 + t127 * t150 - t140 * t165;
t156 = pkin(9) * t110 + t158;
t117 = -t146 * t167 + t150 * t178 - t166 * t177;
t155 = pkin(9) * t117 + t157;
t125 = t139 * t174 + t140 * t143;
t154 = t125 * pkin(2) + pkin(8) * t163 + t164;
t109 = t125 * t183 + (t124 * t145 - t141 * t176) * t150;
t153 = t109 * pkin(3) + t154;
t108 = -t124 * t166 + t125 * t150 + t144 * t165;
t152 = pkin(9) * t108 + t153;
t151 = cos(qJ(5));
t149 = sin(qJ(4));
t135 = pkin(5) * t151 + pkin(4);
t113 = t118 * t182 + t149 * t161;
t103 = t111 * t182 + t149 * t162;
t101 = t109 * t182 + t149 * t163;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t172 - t146 * mrSges(3,3) - (t139 * mrSges(3,1) + t143 * mrSges(3,2)) * t142 - m(4) * t159 - t118 * mrSges(4,1) - t161 * mrSges(4,3) - m(5) * t155 - t113 * mrSges(5,1) - m(6) * (pkin(4) * t113 + t155) - m(7) * (t113 * t135 + t157) + t185 * t117 + t187 * (t113 * t151 + t117 * t148) + t186 * (-t113 * t148 + t117 * t151) + t184 * (t118 * t149 - t161 * t182)) * g(3) + (-mrSges(1,2) - t140 * mrSges(2,1) - t144 * mrSges(2,2) - m(3) * t164 - t125 * mrSges(3,1) - t124 * mrSges(3,2) + mrSges(3,3) * t176 - m(4) * t154 - t109 * mrSges(4,1) - t163 * mrSges(4,3) - m(5) * t152 - t101 * mrSges(5,1) - m(6) * (pkin(4) * t101 + t152) - m(7) * (t101 * t135 + t153) + t187 * (t101 * t151 + t108 * t148) + t186 * (-t101 * t148 + t108 * t151) + t185 * t108 + t184 * (t109 * t149 - t163 * t182)) * g(2) + (-mrSges(1,1) - t144 * mrSges(2,1) + t140 * mrSges(2,2) - m(3) * t171 - t127 * mrSges(3,1) - t126 * mrSges(3,2) - mrSges(3,3) * t180 - m(4) * t160 - t111 * mrSges(4,1) - t162 * mrSges(4,3) - m(5) * t156 - t103 * mrSges(5,1) - m(6) * (pkin(4) * t103 + t156) - m(7) * (t103 * t135 + t158) + t187 * (t103 * t151 + t110 * t148) + t186 * (-t103 * t148 + t110 * t151) + t185 * t110 + t184 * (t111 * t149 - t162 * t182)) * g(1);
U  = t1;

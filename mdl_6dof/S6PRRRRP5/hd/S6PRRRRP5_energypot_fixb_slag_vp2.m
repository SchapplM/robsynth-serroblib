% Calculate potential energy for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:20:05
% EndTime: 2019-03-09 00:20:06
% DurationCPUTime: 0.79s
% Computational Cost: add. (551->116), mult. (1397->153), div. (0->0), fcn. (1753->14), ass. (0->58)
t139 = sin(pkin(7));
t142 = cos(pkin(7));
t143 = cos(pkin(6));
t140 = sin(pkin(6));
t150 = cos(qJ(2));
t175 = t140 * t150;
t160 = -t139 * t175 + t143 * t142;
t138 = sin(pkin(12));
t141 = cos(pkin(12));
t148 = sin(qJ(2));
t172 = t143 * t150;
t125 = -t138 * t172 - t141 * t148;
t177 = t140 * t142;
t161 = -t125 * t139 + t138 * t177;
t186 = -mrSges(6,1) - mrSges(7,1);
t185 = -mrSges(6,2) - mrSges(7,2);
t145 = sin(qJ(5));
t184 = -m(7) * (pkin(5) * t145 + pkin(10)) + mrSges(4,2) - mrSges(5,3);
t183 = -m(6) * pkin(11) + m(7) * (-qJ(6) - pkin(11)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t182 = cos(qJ(3));
t181 = cos(qJ(4));
t179 = t138 * t140;
t178 = t140 * t141;
t176 = t140 * t148;
t173 = t143 * t148;
t171 = t141 * pkin(1) + pkin(8) * t179;
t170 = t143 * pkin(8) + qJ(1);
t166 = t139 * t182;
t165 = t142 * t182;
t164 = t140 * t166;
t163 = t138 * pkin(1) - pkin(8) * t178;
t123 = -t138 * t148 + t141 * t172;
t162 = -t123 * t139 - t141 * t177;
t126 = -t138 * t173 + t141 * t150;
t159 = t126 * pkin(2) + t161 * pkin(9) + t171;
t158 = pkin(2) * t176 + t160 * pkin(9) + t170;
t147 = sin(qJ(3));
t110 = t126 * t182 + (t125 * t142 + t139 * t179) * t147;
t157 = t110 * pkin(3) + t159;
t117 = t143 * t139 * t147 + (t142 * t147 * t150 + t182 * t148) * t140;
t156 = t117 * pkin(3) + t158;
t109 = -t125 * t165 + t126 * t147 - t138 * t164;
t155 = pkin(10) * t109 + t157;
t116 = -t143 * t166 + t147 * t176 - t165 * t175;
t154 = t116 * pkin(10) + t156;
t124 = t138 * t150 + t141 * t173;
t153 = t124 * pkin(2) + t162 * pkin(9) + t163;
t108 = t124 * t182 + (t123 * t142 - t139 * t178) * t147;
t152 = t108 * pkin(3) + t153;
t107 = -t123 * t165 + t124 * t147 + t141 * t164;
t151 = pkin(10) * t107 + t152;
t149 = cos(qJ(5));
t146 = sin(qJ(4));
t134 = pkin(5) * t149 + pkin(4);
t112 = t117 * t181 + t160 * t146;
t104 = t110 * t181 + t161 * t146;
t102 = t108 * t181 + t162 * t146;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t170 - t143 * mrSges(3,3) - (mrSges(3,1) * t148 + mrSges(3,2) * t150) * t140 - m(4) * t158 - t117 * mrSges(4,1) - t160 * mrSges(4,3) - m(5) * t154 - t112 * mrSges(5,1) - m(6) * (t112 * pkin(4) + t154) - m(7) * (t112 * t134 + t156) + t185 * (-t112 * t145 + t116 * t149) + t184 * t116 + t186 * (t112 * t149 + t116 * t145) + t183 * (t117 * t146 - t160 * t181)) * g(3) + (-mrSges(1,2) - t138 * mrSges(2,1) - t141 * mrSges(2,2) - m(3) * t163 - t124 * mrSges(3,1) - t123 * mrSges(3,2) + mrSges(3,3) * t178 - m(4) * t153 - t108 * mrSges(4,1) - t162 * mrSges(4,3) - m(5) * t151 - t102 * mrSges(5,1) - m(6) * (pkin(4) * t102 + t151) - m(7) * (t102 * t134 + t152) + t186 * (t102 * t149 + t107 * t145) + t185 * (-t102 * t145 + t107 * t149) + t184 * t107 + t183 * (t108 * t146 - t162 * t181)) * g(2) + (-mrSges(1,1) - t141 * mrSges(2,1) + t138 * mrSges(2,2) - m(3) * t171 - t126 * mrSges(3,1) - t125 * mrSges(3,2) - mrSges(3,3) * t179 - m(4) * t159 - t110 * mrSges(4,1) - t161 * mrSges(4,3) - m(5) * t155 - t104 * mrSges(5,1) - m(6) * (pkin(4) * t104 + t155) - m(7) * (t104 * t134 + t157) + t186 * (t104 * t149 + t109 * t145) + t185 * (-t104 * t145 + t109 * t149) + t184 * t109 + t183 * (t110 * t146 - t161 * t181)) * g(1);
U  = t1;

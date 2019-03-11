% Calculate potential energy for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:34:22
% EndTime: 2019-03-09 06:34:22
% DurationCPUTime: 0.79s
% Computational Cost: add. (551->116), mult. (1397->151), div. (0->0), fcn. (1753->14), ass. (0->58)
t139 = sin(pkin(12));
t144 = cos(pkin(6));
t151 = cos(qJ(1));
t142 = cos(pkin(12));
t149 = sin(qJ(1));
t172 = t149 * t142;
t126 = -t139 * t151 - t144 * t172;
t140 = sin(pkin(7));
t143 = cos(pkin(7));
t141 = sin(pkin(6));
t177 = t141 * t149;
t162 = -t126 * t140 + t143 * t177;
t178 = t141 * t142;
t161 = -t140 * t178 + t143 * t144;
t187 = -mrSges(6,1) - mrSges(7,1);
t186 = -mrSges(6,2) - mrSges(7,2);
t146 = sin(qJ(5));
t185 = -m(7) * (pkin(5) * t146 + pkin(10)) + mrSges(4,2) - mrSges(5,3);
t184 = -m(6) * pkin(11) + m(7) * (-qJ(6) - pkin(11)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t183 = cos(qJ(3));
t182 = cos(qJ(4));
t181 = t144 * qJ(2) + pkin(8);
t179 = t139 * t141;
t176 = t141 * t151;
t174 = t144 * t151;
t173 = t149 * t139;
t171 = t151 * pkin(1) + qJ(2) * t177;
t167 = t140 * t183;
t166 = t143 * t183;
t165 = t141 * t167;
t164 = t149 * pkin(1) - qJ(2) * t176;
t124 = t142 * t174 - t173;
t163 = -t124 * t140 - t143 * t176;
t127 = t142 * t151 - t144 * t173;
t160 = t127 * pkin(2) + t162 * pkin(9) + t171;
t159 = pkin(2) * t179 + t161 * pkin(9) + t181;
t148 = sin(qJ(3));
t113 = t127 * t183 + (t126 * t143 + t140 * t177) * t148;
t158 = t113 * pkin(3) + t160;
t118 = t144 * t140 * t148 + (t142 * t143 * t148 + t183 * t139) * t141;
t157 = t118 * pkin(3) + t159;
t112 = -t126 * t166 + t127 * t148 - t149 * t165;
t156 = pkin(10) * t112 + t158;
t117 = -t144 * t167 + t148 * t179 - t166 * t178;
t155 = pkin(10) * t117 + t157;
t125 = t139 * t174 + t172;
t154 = t125 * pkin(2) + t163 * pkin(9) + t164;
t111 = t125 * t183 + (t124 * t143 - t140 * t176) * t148;
t153 = t111 * pkin(3) + t154;
t110 = -t124 * t166 + t125 * t148 + t151 * t165;
t152 = t110 * pkin(10) + t153;
t150 = cos(qJ(5));
t147 = sin(qJ(4));
t135 = pkin(5) * t150 + pkin(4);
t109 = t118 * t182 + t161 * t147;
t105 = t113 * t182 + t162 * t147;
t103 = t111 * t182 + t163 * t147;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t181 - t144 * mrSges(3,3) - (t139 * mrSges(3,1) + t142 * mrSges(3,2)) * t141 - m(4) * t159 - t118 * mrSges(4,1) - t161 * mrSges(4,3) - m(5) * t155 - t109 * mrSges(5,1) - m(6) * (pkin(4) * t109 + t155) - m(7) * (t109 * t135 + t157) + t185 * t117 + t187 * (t109 * t150 + t117 * t146) + t186 * (-t109 * t146 + t117 * t150) + t184 * (t118 * t147 - t161 * t182)) * g(3) + (-mrSges(1,2) - t149 * mrSges(2,1) - t151 * mrSges(2,2) - m(3) * t164 - t125 * mrSges(3,1) - t124 * mrSges(3,2) + mrSges(3,3) * t176 - m(4) * t154 - t111 * mrSges(4,1) - t163 * mrSges(4,3) - m(5) * t152 - t103 * mrSges(5,1) - m(6) * (t103 * pkin(4) + t152) - m(7) * (t103 * t135 + t153) + t187 * (t103 * t150 + t110 * t146) + t186 * (-t103 * t146 + t110 * t150) + t185 * t110 + t184 * (t111 * t147 - t163 * t182)) * g(2) + (-mrSges(1,1) - t151 * mrSges(2,1) + t149 * mrSges(2,2) - m(3) * t171 - t127 * mrSges(3,1) - t126 * mrSges(3,2) - mrSges(3,3) * t177 - m(4) * t160 - t113 * mrSges(4,1) - t162 * mrSges(4,3) - m(5) * t156 - t105 * mrSges(5,1) - m(6) * (pkin(4) * t105 + t156) - m(7) * (t105 * t135 + t158) + t187 * (t105 * t150 + t112 * t146) + t186 * (-t105 * t146 + t112 * t150) + t185 * t112 + t184 * (t113 * t147 - t162 * t182)) * g(1);
U  = t1;

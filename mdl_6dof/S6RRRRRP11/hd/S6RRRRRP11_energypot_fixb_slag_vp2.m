% Calculate potential energy for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:32:39
% EndTime: 2019-03-10 02:32:40
% DurationCPUTime: 0.75s
% Computational Cost: add. (551->116), mult. (1397->150), div. (0->0), fcn. (1753->14), ass. (0->59)
t142 = cos(pkin(6));
t148 = sin(qJ(1));
t150 = cos(qJ(2));
t173 = t148 * t150;
t147 = sin(qJ(2));
t151 = cos(qJ(1));
t175 = t147 * t151;
t126 = -t142 * t173 - t175;
t139 = sin(pkin(7));
t141 = cos(pkin(7));
t140 = sin(pkin(6));
t179 = t140 * t148;
t162 = -t126 * t139 + t141 * t179;
t178 = t140 * t150;
t161 = -t139 * t178 + t141 * t142;
t188 = -mrSges(6,1) - mrSges(7,1);
t187 = -mrSges(6,2) - mrSges(7,2);
t144 = sin(qJ(5));
t186 = -m(7) * (pkin(5) * t144 + pkin(11)) + mrSges(4,2) - mrSges(5,3);
t185 = -m(6) * pkin(12) + m(7) * (-qJ(6) - pkin(12)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t184 = cos(qJ(3));
t183 = cos(qJ(4));
t182 = t142 * pkin(9) + pkin(8);
t180 = t140 * t147;
t177 = t140 * t151;
t174 = t148 * t147;
t172 = t150 * t151;
t171 = t151 * pkin(1) + pkin(9) * t179;
t167 = t139 * t184;
t166 = t141 * t184;
t165 = t140 * t167;
t164 = t148 * pkin(1) - pkin(9) * t177;
t124 = t142 * t172 - t174;
t163 = -t124 * t139 - t141 * t177;
t127 = -t142 * t174 + t172;
t160 = t127 * pkin(2) + t162 * pkin(10) + t171;
t159 = pkin(2) * t180 + t161 * pkin(10) + t182;
t146 = sin(qJ(3));
t113 = t127 * t184 + (t126 * t141 + t139 * t179) * t146;
t158 = t113 * pkin(3) + t160;
t118 = t142 * t139 * t146 + (t141 * t146 * t150 + t184 * t147) * t140;
t157 = t118 * pkin(3) + t159;
t112 = -t126 * t166 + t127 * t146 - t148 * t165;
t156 = pkin(11) * t112 + t158;
t117 = -t142 * t167 + t146 * t180 - t166 * t178;
t155 = pkin(11) * t117 + t157;
t125 = t142 * t175 + t173;
t154 = t125 * pkin(2) + t163 * pkin(10) + t164;
t111 = t125 * t184 + (t124 * t141 - t139 * t177) * t146;
t153 = t111 * pkin(3) + t154;
t110 = -t124 * t166 + t125 * t146 + t151 * t165;
t152 = t110 * pkin(11) + t153;
t149 = cos(qJ(5));
t145 = sin(qJ(4));
t135 = pkin(5) * t149 + pkin(4);
t109 = t118 * t183 + t161 * t145;
t105 = t113 * t183 + t162 * t145;
t103 = t111 * t183 + t163 * t145;
t1 = (-mrSges(1,3) - m(2) * pkin(8) - mrSges(2,3) - m(3) * t182 - t142 * mrSges(3,3) - (t147 * mrSges(3,1) + t150 * mrSges(3,2)) * t140 - m(4) * t159 - t118 * mrSges(4,1) - t161 * mrSges(4,3) - m(5) * t155 - t109 * mrSges(5,1) - m(6) * (pkin(4) * t109 + t155) - m(7) * (t109 * t135 + t157) + t186 * t117 + t188 * (t109 * t149 + t117 * t144) + t187 * (-t109 * t144 + t117 * t149) + t185 * (t118 * t145 - t161 * t183)) * g(3) + (-mrSges(1,2) - t148 * mrSges(2,1) - mrSges(2,2) * t151 - m(3) * t164 - t125 * mrSges(3,1) - t124 * mrSges(3,2) + mrSges(3,3) * t177 - m(4) * t154 - t111 * mrSges(4,1) - t163 * mrSges(4,3) - m(5) * t152 - t103 * mrSges(5,1) - m(6) * (t103 * pkin(4) + t152) - m(7) * (t103 * t135 + t153) + t188 * (t103 * t149 + t110 * t144) + t187 * (-t103 * t144 + t110 * t149) + t186 * t110 + t185 * (t111 * t145 - t163 * t183)) * g(2) + (-mrSges(1,1) - mrSges(2,1) * t151 + t148 * mrSges(2,2) - m(3) * t171 - t127 * mrSges(3,1) - t126 * mrSges(3,2) - mrSges(3,3) * t179 - m(4) * t160 - t113 * mrSges(4,1) - t162 * mrSges(4,3) - m(5) * t156 - t105 * mrSges(5,1) - m(6) * (pkin(4) * t105 + t156) - m(7) * (t105 * t135 + t158) + t188 * (t105 * t149 + t112 * t144) + t187 * (-t105 * t144 + t112 * t149) + t186 * t112 + t185 * (t113 * t145 - t162 * t183)) * g(1);
U  = t1;

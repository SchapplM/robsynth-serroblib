% Calculate potential energy for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR15_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:23:52
% EndTime: 2024-09-27 22:23:53
% DurationCPUTime: 0.14s
% Computational Cost: add. (316->120), mult. (367->150), div. (0->0), fcn. (365->26), ass. (0->69)
t123 = qJ(2) + qJ(3);
t116 = pkin(5) + t123;
t110 = qJ(4) + t116;
t101 = sin(t110) / 0.2e1;
t117 = pkin(5) - t123;
t111 = -qJ(4) + t117;
t102 = cos(t111) / 0.2e1;
t138 = pkin(8) + pkin(9);
t122 = pkin(10) + t138;
t126 = cos(pkin(6));
t103 = pkin(11) * t126 + t122;
t104 = sin(t111);
t105 = cos(t110);
t106 = sin(t116);
t107 = sin(t117);
t108 = cos(t116);
t109 = cos(t117);
t125 = sin(pkin(5));
t127 = cos(pkin(5));
t130 = sin(qJ(3));
t131 = sin(qJ(2));
t135 = cos(qJ(3));
t136 = cos(qJ(2));
t139 = pkin(3) * t130 * t136 + (pkin(3) * t135 + pkin(2)) * t131;
t129 = sin(qJ(4));
t134 = cos(qJ(4));
t124 = sin(pkin(6));
t162 = pkin(11) * t124;
t96 = pkin(4) * t134 + t129 * t162 + pkin(3);
t99 = -pkin(4) * t129 + t134 * t162;
t90 = t130 * t99 + t135 * t96 + pkin(2);
t91 = -t130 * t96 + t135 * t99;
t140 = t131 * t90 - t136 * t91;
t163 = pkin(2) * t131;
t171 = -(t106 - t107) * mrSges(4,1) / 0.2e1 - (t109 + t108) * mrSges(4,2) / 0.2e1 - m(4) * (-t125 * t138 + t127 * t163) - m(5) * (-t122 * t125 + t127 * t139) + m(6) * (t103 * t125 - t127 * t140) - (t101 - t104 / 0.2e1) * mrSges(5,1) - (t102 + t105 / 0.2e1) * mrSges(5,2) - mrSges(2,2);
t170 = t124 * mrSges(6,3);
t168 = m(3) * pkin(8) + t126 * mrSges(6,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t115 = t136 * pkin(2) + pkin(1);
t119 = cos(t123);
t166 = -m(3) * pkin(1) - m(4) * t115 - m(5) * (pkin(3) * t119 + t115) - m(6) * (t131 * t91 + t136 * t90 + pkin(1)) - t119 * mrSges(4,1) + sin(t123) * mrSges(4,2) - mrSges(2,1);
t120 = qJ(4) + t123;
t112 = sin(t120);
t132 = sin(qJ(1));
t161 = t112 * t132;
t137 = cos(qJ(1));
t160 = t112 * t137;
t113 = cos(t120);
t159 = t113 * t126;
t158 = t113 * t132;
t157 = t113 * t137;
t156 = t125 * t132;
t155 = t125 * t137;
t128 = sin(qJ(5));
t154 = t128 * t132;
t153 = t128 * t137;
t152 = t131 * t132;
t151 = t131 * t137;
t133 = cos(qJ(5));
t150 = t132 * t133;
t149 = t132 * t136;
t148 = t133 * t137;
t147 = t136 * t137;
t146 = t124 * t155;
t145 = t127 * t154;
t144 = t127 * t153;
t143 = t127 * t150;
t142 = t127 * t148;
t141 = t124 * t156;
t1 = (-mrSges(1,3) - mrSges(2,3) - (t109 / 0.2e1 - t108 / 0.2e1) * mrSges(4,1) - (t106 / 0.2e1 + t107 / 0.2e1) * mrSges(4,2) - (t102 - t105 / 0.2e1) * mrSges(5,1) - (t101 + t104 / 0.2e1) * mrSges(5,2) + (-t131 * mrSges(3,1) - t136 * mrSges(3,2) - m(4) * t163 - m(5) * t139 - m(6) * t140 - (t112 * t133 + t128 * t159) * mrSges(6,1) - (-t112 * t128 + t133 * t159) * mrSges(6,2) + t113 * t170) * t125 + (-m(2) - m(3) - m(4) - m(5) - m(6)) * pkin(7) + (-m(4) * t138 - m(5) * t122 - m(6) * t103 - (mrSges(6,1) * t128 + mrSges(6,2) * t133) * t124 - t168) * t127) * g(3) + (-mrSges(1,2) - (t127 * t151 + t149) * mrSges(3,1) - (t127 * t147 - t152) * mrSges(3,2) - t158 * mrSges(5,1) + t161 * mrSges(5,2) - ((t126 * t144 + t150) * t113 + (-t126 * t154 + t142) * t112 - t128 * t146) * mrSges(6,1) - ((t126 * t142 - t154) * t113 + (-t126 * t150 - t144) * t112 - t133 * t146) * mrSges(6,2) - (-t127 * t157 + t161) * t170 + t166 * t132 + t168 * t155 + t171 * t137) * g(2) + (-mrSges(1,1) - (-t127 * t152 + t147) * mrSges(3,1) - (-t127 * t149 - t151) * mrSges(3,2) - t157 * mrSges(5,1) + t160 * mrSges(5,2) - ((-t126 * t145 + t148) * t113 + (-t126 * t153 - t143) * t112 + t128 * t141) * mrSges(6,1) - ((-t126 * t143 - t153) * t113 + (-t126 * t148 + t145) * t112 + t133 * t141) * mrSges(6,2) - (t127 * t158 + t160) * t170 + t166 * t137 - t168 * t156 - t171 * t132) * g(1);
U = t1;

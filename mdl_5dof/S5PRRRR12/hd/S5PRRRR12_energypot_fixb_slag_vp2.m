% Calculate potential energy for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:09
% EndTime: 2024-09-28 18:07:10
% DurationCPUTime: 0.50s
% Computational Cost: add. (316->127), mult. (405->162), div. (0->0), fcn. (399->26), ass. (0->67)
t134 = qJ(2) + qJ(3);
t131 = qJ(4) + t134;
t122 = cos(t131);
t139 = cos(pkin(6));
t181 = t122 * t139;
t125 = pkin(5) + t134;
t119 = qJ(4) + t125;
t109 = sin(t119) / 0.2e1;
t126 = pkin(5) - t134;
t120 = -qJ(4) + t126;
t110 = cos(t120) / 0.2e1;
t112 = sin(t125) / 0.2e1;
t113 = cos(t126) / 0.2e1;
t115 = sin(t120);
t116 = cos(t119);
t117 = sin(t126);
t118 = cos(t125);
t121 = sin(t131);
t149 = pkin(7) + pkin(8);
t133 = pkin(9) + t149;
t137 = sin(pkin(5));
t140 = cos(pkin(5));
t141 = sin(qJ(5));
t145 = cos(qJ(5));
t148 = cos(qJ(2));
t143 = sin(qJ(3));
t144 = sin(qJ(2));
t147 = cos(qJ(3));
t152 = pkin(3) * t143 * t148 + (pkin(3) * t147 + pkin(2)) * t144;
t160 = t140 * t145;
t161 = t140 * t144;
t162 = t140 * t141;
t136 = sin(pkin(6));
t170 = t136 * t137;
t175 = m(3) * pkin(7) + m(6) * (pkin(10) * t139 + t133) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t180 = (t160 * t121 - t141 * t170 + t162 * t181) * mrSges(6,1) - (t162 * t121 + t145 * t170 - t160 * t181) * mrSges(6,2) + m(4) * (pkin(2) * t161 - t137 * t149) + m(5) * (-t133 * t137 + t140 * t152) + t161 * mrSges(3,1) + (t112 - t117 / 0.2e1) * mrSges(4,1) + (t109 - t115 / 0.2e1) * mrSges(5,1) + t140 * t148 * mrSges(3,2) + (t113 + t118 / 0.2e1) * mrSges(4,2) + (t110 + t116 / 0.2e1) * mrSges(5,2) - t175 * t137 + mrSges(2,2);
t124 = t148 * pkin(2) + pkin(1);
t128 = cos(t134);
t165 = t139 * t145;
t166 = t139 * t141;
t176 = -mrSges(2,1) + t144 * mrSges(3,2) - m(4) * t124 - t128 * mrSges(4,1) + sin(t134) * mrSges(4,2) - m(5) * (pkin(3) * t128 + t124) - t122 * mrSges(5,1) - (-t166 * t121 + t145 * t122) * mrSges(6,1) - (-t165 * t121 - t141 * t122) * mrSges(6,2) + (-m(6) * pkin(2) - mrSges(3,1)) * t148;
t174 = pkin(10) * t136;
t135 = sin(pkin(11));
t173 = t121 * t135;
t172 = t122 * t136;
t138 = cos(pkin(11));
t169 = t136 * t138;
t167 = t137 * t139;
t164 = t140 * t135;
t163 = t140 * t138;
t158 = t135 * t174;
t157 = pkin(10) * t169;
t101 = pkin(4) * t164 - t157;
t142 = sin(qJ(4));
t146 = cos(qJ(4));
t99 = pkin(4) * t138 + t140 * t158;
t151 = pkin(3) * t138 - t101 * t142 + t146 * t99;
t100 = -pkin(4) * t135 + t140 * t157;
t102 = pkin(4) * t163 + t158;
t150 = pkin(3) * t135 - t100 * t146 + t102 * t142;
t130 = t138 * pkin(1);
t129 = t135 * pkin(1);
t107 = -pkin(4) * t142 + t146 * t174;
t106 = pkin(4) * t146 + t142 * t174 + pkin(3);
t95 = pkin(3) * t163 + t100 * t142 + t102 * t146;
t94 = -pkin(3) * t164 - t101 * t146 - t142 * t99;
t1 = (-mrSges(1,3) - mrSges(2,3) - (t113 - t118 / 0.2e1) * mrSges(4,1) - (t112 + t117 / 0.2e1) * mrSges(4,2) - (t110 - t116 / 0.2e1) * mrSges(5,1) - (t109 + t115 / 0.2e1) * mrSges(5,2) + (-m(5) * t152 - (t121 * t145 + t122 * t166) * mrSges(6,1) - (-t121 * t141 + t122 * t165) * mrSges(6,2) + mrSges(6,3) * t172 + (-mrSges(3,2) + m(6) * (-t106 * t143 + t107 * t147)) * t148 + (-mrSges(3,1) - m(4) * pkin(2) - m(6) * (t106 * t147 + t107 * t143 + pkin(2))) * t144) * t137 + (-m(2) - m(3) - m(4) - m(5) - m(6)) * qJ(1) + (-m(4) * t149 - m(5) * t133 - t139 * mrSges(6,3) - (mrSges(6,1) * t141 + mrSges(6,2) * t145) * t136 - t175) * t140) * g(3) + (-mrSges(1,2) - m(3) * t129 + t173 * mrSges(5,2) - m(6) * ((t95 * t143 + t147 * t150) * t148 + (pkin(2) * t163 - t143 * t150 + t95 * t147) * t144 + t129) - (-t122 * t163 + t173) * t136 * mrSges(6,3) + t176 * t135 + (t167 * mrSges(6,3) - t180) * t138) * g(2) + (-mrSges(1,1) - m(3) * t130 - m(6) * ((t94 * t143 + t147 * t151) * t148 + (-pkin(2) * t164 - t143 * t151 + t94 * t147) * t144 + t130) - t121 * t169 * mrSges(6,3) + (t121 * mrSges(5,2) + t176) * t138 + (-(t140 * t172 + t167) * mrSges(6,3) + t180) * t135) * g(1);
U = t1;

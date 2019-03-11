% Calculate potential energy for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:27:40
% EndTime: 2019-03-09 13:27:40
% DurationCPUTime: 0.75s
% Computational Cost: add. (383->95), mult. (733->117), div. (0->0), fcn. (880->14), ass. (0->54)
t164 = -m(4) - m(5);
t163 = -m(6) - m(7);
t162 = -m(3) * pkin(1) - mrSges(2,1);
t161 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t132 = cos(qJ(4));
t160 = -t132 * mrSges(5,2) - mrSges(3,3) - mrSges(4,3);
t159 = m(3) * pkin(8) - t160;
t128 = sin(qJ(4));
t158 = -m(5) * pkin(3) - t132 * mrSges(5,1) + t128 * mrSges(5,2) - mrSges(4,1);
t127 = sin(qJ(6));
t131 = cos(qJ(6));
t157 = -m(7) * pkin(5) - t131 * mrSges(7,1) + t127 * mrSges(7,2) - mrSges(6,1);
t156 = m(5) * pkin(9) + t127 * mrSges(7,1) + t131 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t126 = cos(pkin(6));
t155 = t126 * pkin(8) + pkin(7);
t124 = sin(pkin(6));
t129 = sin(qJ(2));
t154 = t124 * t129;
t130 = sin(qJ(1));
t153 = t124 * t130;
t134 = cos(qJ(1));
t152 = t124 * t134;
t125 = cos(pkin(12));
t133 = cos(qJ(2));
t151 = t125 * t133;
t150 = t126 * t128;
t149 = t129 * t130;
t148 = t129 * t134;
t147 = t130 * t133;
t146 = t133 * t134;
t105 = pkin(2) * t126 * t129 + (-pkin(8) - qJ(3)) * t124;
t117 = pkin(2) * t133 + pkin(1);
t145 = t134 * t105 + t130 * t117;
t144 = t128 * t153;
t143 = t128 * t152;
t142 = -t105 * t130 + t134 * t117;
t141 = pkin(2) * t154 + t126 * qJ(3) + t155;
t123 = sin(pkin(12));
t140 = t123 * t133 + t125 * t129;
t107 = -t123 * t129 + t151;
t137 = t107 * t126;
t135 = -pkin(10) - pkin(9);
t122 = qJ(4) + qJ(5);
t120 = cos(t122);
t119 = sin(t122);
t116 = pkin(4) * t132 + pkin(3);
t104 = t140 * t126;
t103 = t140 * t124;
t102 = t123 * t154 - t124 * t151;
t95 = -t104 * t130 + t107 * t134;
t94 = -t130 * t137 - t134 * t140;
t93 = t104 * t134 + t107 * t130;
t92 = -t130 * t140 + t134 * t137;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t155 - (mrSges(3,1) * t129 + mrSges(3,2) * t133) * t124 - t150 * mrSges(5,1) + t163 * (pkin(4) * t150 - t102 * t135 + t103 * t116 + t141) + t161 * (t103 * t119 - t126 * t120) + t164 * t141 + t158 * t103 + t160 * t126 + t157 * (t103 * t120 + t119 * t126) - t156 * t102) * g(3) + (-mrSges(1,2) - t134 * mrSges(2,2) - (t126 * t148 + t147) * mrSges(3,1) - (t126 * t146 - t149) * mrSges(3,2) + t143 * mrSges(5,1) + t164 * t145 + t158 * t93 + t163 * (-pkin(4) * t143 + t93 * t116 + t92 * t135 + t145) + t161 * (t119 * t93 + t120 * t152) + t162 * t130 + t159 * t152 + t157 * (-t119 * t152 + t120 * t93) + t156 * t92) * g(2) + (-mrSges(1,1) + t130 * mrSges(2,2) - (-t126 * t149 + t146) * mrSges(3,1) - (-t126 * t147 - t148) * mrSges(3,2) - t144 * mrSges(5,1) + t164 * t142 + t158 * t95 + t163 * (pkin(4) * t144 + t95 * t116 + t94 * t135 + t142) + t161 * (t119 * t95 - t120 * t153) + t162 * t134 - t159 * t153 + t157 * (t119 * t153 + t120 * t95) + t156 * t94) * g(1);
U  = t1;

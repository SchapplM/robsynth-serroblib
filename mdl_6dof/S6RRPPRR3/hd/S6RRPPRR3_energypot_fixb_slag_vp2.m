% Calculate potential energy for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:25
% EndTime: 2019-03-09 08:54:25
% DurationCPUTime: 0.76s
% Computational Cost: add. (383->95), mult. (733->117), div. (0->0), fcn. (880->14), ass. (0->54)
t164 = -m(4) - m(5);
t163 = -m(6) - m(7);
t162 = -m(3) * pkin(1) - mrSges(2,1);
t161 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t126 = cos(pkin(12));
t160 = -t126 * mrSges(5,2) - mrSges(3,3) - mrSges(4,3);
t159 = m(3) * pkin(8) - t160;
t123 = sin(pkin(12));
t158 = -m(5) * pkin(3) - t126 * mrSges(5,1) + t123 * mrSges(5,2) - mrSges(4,1);
t130 = sin(qJ(6));
t133 = cos(qJ(6));
t157 = -m(7) * pkin(5) - t133 * mrSges(7,1) + t130 * mrSges(7,2) - mrSges(6,1);
t156 = m(5) * qJ(4) + t130 * mrSges(7,1) + t133 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t128 = cos(pkin(6));
t155 = t128 * pkin(8) + pkin(7);
t154 = t123 * t128;
t125 = sin(pkin(6));
t131 = sin(qJ(2));
t153 = t125 * t131;
t132 = sin(qJ(1));
t152 = t125 * t132;
t135 = cos(qJ(1));
t151 = t125 * t135;
t127 = cos(pkin(11));
t134 = cos(qJ(2));
t150 = t127 * t134;
t149 = t131 * t135;
t148 = t132 * t131;
t147 = t132 * t134;
t146 = t134 * t135;
t105 = pkin(2) * t128 * t131 + (-pkin(8) - qJ(3)) * t125;
t117 = pkin(2) * t134 + pkin(1);
t145 = t135 * t105 + t132 * t117;
t144 = t123 * t152;
t143 = t123 * t151;
t142 = -t105 * t132 + t135 * t117;
t141 = pkin(2) * t153 + t128 * qJ(3) + t155;
t124 = sin(pkin(11));
t140 = t124 * t134 + t127 * t131;
t107 = -t124 * t131 + t150;
t137 = t107 * t128;
t129 = -pkin(9) - qJ(4);
t122 = pkin(12) + qJ(5);
t119 = cos(t122);
t118 = sin(t122);
t116 = pkin(4) * t126 + pkin(3);
t104 = t140 * t128;
t103 = t140 * t125;
t102 = t124 * t153 - t125 * t150;
t95 = -t132 * t104 + t107 * t135;
t94 = -t132 * t137 - t135 * t140;
t93 = t104 * t135 + t132 * t107;
t92 = -t132 * t140 + t135 * t137;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t155 - (mrSges(3,1) * t131 + mrSges(3,2) * t134) * t125 - t154 * mrSges(5,1) + t163 * (pkin(4) * t154 - t102 * t129 + t103 * t116 + t141) + t161 * (t103 * t118 - t128 * t119) + t164 * t141 + t158 * t103 + t160 * t128 + t157 * (t103 * t119 + t118 * t128) - t156 * t102) * g(3) + (-mrSges(1,2) - t135 * mrSges(2,2) - (t128 * t149 + t147) * mrSges(3,1) - (t128 * t146 - t148) * mrSges(3,2) + t143 * mrSges(5,1) + t164 * t145 + t158 * t93 + t163 * (-pkin(4) * t143 + t93 * t116 + t92 * t129 + t145) + t161 * (t93 * t118 + t119 * t151) + t162 * t132 + t159 * t151 + t157 * (-t118 * t151 + t93 * t119) + t156 * t92) * g(2) + (-mrSges(1,1) + t132 * mrSges(2,2) - (-t128 * t148 + t146) * mrSges(3,1) - (-t128 * t147 - t149) * mrSges(3,2) - t144 * mrSges(5,1) + t164 * t142 + t158 * t95 + t163 * (pkin(4) * t144 + t95 * t116 + t94 * t129 + t142) + t161 * (t118 * t95 - t119 * t152) + t162 * t135 - t159 * t152 + t157 * (t118 * t152 + t119 * t95) + t156 * t94) * g(1);
U  = t1;

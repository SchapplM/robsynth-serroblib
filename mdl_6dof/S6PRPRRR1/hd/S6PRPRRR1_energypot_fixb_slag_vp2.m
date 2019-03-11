% Calculate potential energy for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:14
% EndTime: 2019-03-08 20:22:15
% DurationCPUTime: 0.79s
% Computational Cost: add. (383->92), mult. (733->111), div. (0->0), fcn. (880->14), ass. (0->55)
t136 = cos(qJ(4));
t170 = t136 * mrSges(5,2);
t137 = cos(qJ(2));
t169 = t137 * mrSges(3,2);
t168 = -m(4) - m(5);
t167 = -m(6) - m(7);
t166 = mrSges(3,3) + mrSges(4,3);
t165 = m(3) * pkin(7) + t166;
t164 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t134 = sin(qJ(2));
t163 = -m(3) * pkin(1) - t137 * mrSges(3,1) + t134 * mrSges(3,2) - mrSges(2,1);
t133 = sin(qJ(4));
t162 = -m(5) * pkin(3) - mrSges(5,1) * t136 + mrSges(5,2) * t133 - mrSges(4,1);
t132 = sin(qJ(6));
t135 = cos(qJ(6));
t161 = -m(7) * pkin(5) - t135 * mrSges(7,1) + t132 * mrSges(7,2) - mrSges(6,1);
t128 = sin(pkin(6));
t131 = cos(pkin(6));
t151 = t131 * t134;
t160 = t151 * mrSges(3,1) - t128 * t170 + t131 * t169 + mrSges(2,2);
t159 = m(5) * pkin(8) + t132 * mrSges(7,1) + t135 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t127 = sin(pkin(11));
t158 = t127 * t128;
t130 = cos(pkin(11));
t157 = t128 * t130;
t156 = t128 * t133;
t155 = t128 * t134;
t129 = cos(pkin(12));
t153 = t129 * t137;
t152 = t131 * t133;
t108 = pkin(2) * t151 + (-pkin(7) - qJ(3)) * t128;
t120 = pkin(2) * t137 + pkin(1);
t149 = t130 * t108 + t127 * t120;
t148 = t131 * pkin(7) + qJ(1);
t147 = t127 * t156;
t146 = t130 * t156;
t145 = -t108 * t127 + t130 * t120;
t144 = pkin(2) * t155 + t131 * qJ(3) + t148;
t126 = sin(pkin(12));
t143 = t126 * t137 + t129 * t134;
t110 = -t126 * t134 + t153;
t141 = t110 * t131;
t138 = -pkin(9) - pkin(8);
t125 = qJ(4) + qJ(5);
t123 = cos(t125);
t122 = sin(t125);
t119 = pkin(4) * t136 + pkin(3);
t107 = t143 * t131;
t106 = t143 * t128;
t105 = t126 * t155 - t128 * t153;
t98 = -t107 * t127 + t110 * t130;
t97 = -t127 * t141 - t130 * t143;
t96 = t107 * t130 + t110 * t127;
t95 = -t127 * t143 + t130 * t141;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t148 - (t134 * mrSges(3,1) + t169) * t128 - t152 * mrSges(5,1) + t167 * (pkin(4) * t152 - t105 * t138 + t106 * t119 + t144) + t164 * (t106 * t122 - t131 * t123) + t168 * t144 + t162 * t106 + (-t166 - t170) * t131 + t161 * (t106 * t123 + t122 * t131) - t159 * t105) * g(3) + (t146 * mrSges(5,1) - mrSges(1,2) + t168 * t149 + t162 * t96 + t167 * (-pkin(4) * t146 + t96 * t119 + t95 * t138 + t149) + t164 * (t122 * t96 + t123 * t157) - t160 * t130 + t163 * t127 + t165 * t157 + t161 * (-t122 * t157 + t123 * t96) + t159 * t95) * g(2) + (-t147 * mrSges(5,1) - mrSges(1,1) + t168 * t145 + t162 * t98 + t167 * (pkin(4) * t147 + t98 * t119 + t97 * t138 + t145) + t164 * (t122 * t98 - t123 * t158) + t160 * t127 + t163 * t130 - t165 * t158 + t161 * (t122 * t158 + t123 * t98) + t159 * t97) * g(1);
U  = t1;

% Calculate potential energy for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:22:33
% EndTime: 2019-03-08 21:22:34
% DurationCPUTime: 0.75s
% Computational Cost: add. (343->112), mult. (554->131), div. (0->0), fcn. (633->12), ass. (0->47)
t145 = -mrSges(6,1) - mrSges(7,1);
t144 = -mrSges(6,2) - mrSges(7,2);
t143 = m(4) * pkin(8) - mrSges(3,2) + mrSges(4,3);
t142 = -mrSges(5,3) - t143;
t141 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t119 = sin(qJ(3));
t140 = pkin(3) * t119;
t118 = sin(qJ(5));
t112 = sin(pkin(10));
t114 = cos(pkin(10));
t120 = sin(qJ(2));
t115 = cos(pkin(6));
t123 = cos(qJ(2));
t132 = t115 * t123;
t92 = t112 * t120 - t114 * t132;
t139 = t118 * t92;
t94 = t112 * t132 + t114 * t120;
t138 = t118 * t94;
t113 = sin(pkin(6));
t137 = t112 * t113;
t136 = t113 * t120;
t135 = t113 * t123;
t134 = t114 * t113;
t133 = t115 * t120;
t131 = t114 * pkin(1) + pkin(7) * t137;
t130 = t115 * pkin(7) + qJ(1);
t129 = t119 * t137;
t128 = t118 * t135;
t108 = t112 * pkin(1);
t127 = -pkin(7) * t134 + t108;
t122 = cos(qJ(3));
t105 = pkin(3) * t122 + pkin(2);
t117 = -qJ(4) - pkin(8);
t95 = -t112 * t133 + t114 * t123;
t126 = pkin(3) * t129 + t95 * t105 - t94 * t117 + t131;
t125 = t105 * t136 + t115 * t140 + t117 * t135 + t130;
t93 = t112 * t123 + t114 * t133;
t124 = t108 + t93 * t105 - t92 * t117 + (-pkin(7) - t140) * t134;
t121 = cos(qJ(5));
t111 = qJ(3) + pkin(11);
t107 = cos(t111);
t106 = sin(t111);
t104 = pkin(5) * t121 + pkin(4);
t89 = t106 * t115 + t107 * t136;
t83 = t106 * t137 + t107 * t95;
t81 = -t106 * t134 + t107 * t93;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(5) * t125 - t89 * mrSges(5,1) + mrSges(5,3) * t135 - m(6) * (pkin(4) * t89 + t125) - m(7) * (-pkin(5) * t128 + t89 * t104 + t125) + t145 * (t89 * t121 - t128) + t144 * (-t89 * t118 - t121 * t135) + (-m(3) - m(4)) * t130 + (-t119 * mrSges(4,1) - t122 * mrSges(4,2) - mrSges(3,3)) * t115 + (t143 * t123 + (-m(4) * pkin(2) - t122 * mrSges(4,1) + t119 * mrSges(4,2) - mrSges(3,1)) * t120) * t113 + t141 * (t106 * t136 - t115 * t107)) * g(3) + (-mrSges(1,2) - t112 * mrSges(2,1) - t114 * mrSges(2,2) - m(3) * t127 - t93 * mrSges(3,1) + mrSges(3,3) * t134 - m(4) * (pkin(2) * t93 + t127) - (-t119 * t134 + t122 * t93) * mrSges(4,1) - (-t119 * t93 - t122 * t134) * mrSges(4,2) - m(5) * t124 - t81 * mrSges(5,1) - m(6) * (pkin(4) * t81 + t124) - m(7) * (pkin(5) * t139 + t104 * t81 + t124) + t145 * (t121 * t81 + t139) + t144 * (-t118 * t81 + t121 * t92) + t142 * t92 + t141 * (t106 * t93 + t107 * t134)) * g(2) + (-mrSges(1,1) - t114 * mrSges(2,1) + t112 * mrSges(2,2) - m(3) * t131 - t95 * mrSges(3,1) - mrSges(3,3) * t137 - m(4) * (pkin(2) * t95 + t131) - (t122 * t95 + t129) * mrSges(4,1) - (-t119 * t95 + t122 * t137) * mrSges(4,2) - m(5) * t126 - t83 * mrSges(5,1) - m(6) * (pkin(4) * t83 + t126) - m(7) * (pkin(5) * t138 + t104 * t83 + t126) + t145 * (t121 * t83 + t138) + t144 * (-t118 * t83 + t121 * t94) + t142 * t94 + t141 * (t106 * t95 - t107 * t137)) * g(1);
U  = t1;

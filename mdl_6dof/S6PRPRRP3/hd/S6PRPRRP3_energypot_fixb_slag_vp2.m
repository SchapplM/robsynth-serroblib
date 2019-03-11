% Calculate potential energy for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:36
% EndTime: 2019-03-08 20:04:37
% DurationCPUTime: 0.74s
% Computational Cost: add. (343->112), mult. (554->131), div. (0->0), fcn. (633->12), ass. (0->47)
t143 = -mrSges(6,1) - mrSges(7,1);
t142 = -mrSges(6,2) - mrSges(7,2);
t141 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t140 = -mrSges(5,3) - t141;
t139 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t110 = sin(pkin(11));
t138 = pkin(3) * t110;
t118 = sin(qJ(5));
t111 = sin(pkin(10));
t114 = cos(pkin(10));
t119 = sin(qJ(2));
t115 = cos(pkin(6));
t121 = cos(qJ(2));
t130 = t115 * t121;
t90 = t111 * t119 - t114 * t130;
t137 = t118 * t90;
t92 = t111 * t130 + t114 * t119;
t136 = t118 * t92;
t112 = sin(pkin(6));
t135 = t111 * t112;
t134 = t112 * t114;
t133 = t112 * t119;
t132 = t112 * t121;
t131 = t115 * t119;
t129 = t114 * pkin(1) + pkin(7) * t135;
t128 = t115 * pkin(7) + qJ(1);
t127 = t110 * t135;
t126 = t118 * t132;
t106 = t111 * pkin(1);
t125 = -pkin(7) * t134 + t106;
t113 = cos(pkin(11));
t102 = pkin(3) * t113 + pkin(2);
t117 = -pkin(8) - qJ(3);
t93 = -t111 * t131 + t114 * t121;
t124 = pkin(3) * t127 + t93 * t102 - t92 * t117 + t129;
t123 = t102 * t133 + t115 * t138 + t117 * t132 + t128;
t91 = t111 * t121 + t114 * t131;
t122 = t106 + t91 * t102 - t90 * t117 + (-pkin(7) - t138) * t134;
t120 = cos(qJ(5));
t109 = pkin(11) + qJ(4);
t105 = cos(t109);
t104 = sin(t109);
t103 = pkin(5) * t120 + pkin(4);
t87 = t104 * t115 + t105 * t133;
t81 = t104 * t135 + t105 * t93;
t79 = -t104 * t134 + t105 * t91;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(5) * t123 - t87 * mrSges(5,1) + mrSges(5,3) * t132 - m(6) * (pkin(4) * t87 + t123) - m(7) * (-pkin(5) * t126 + t87 * t103 + t123) + t143 * (t87 * t120 - t126) + t142 * (-t87 * t118 - t120 * t132) + (-m(3) - m(4)) * t128 + (-t110 * mrSges(4,1) - t113 * mrSges(4,2) - mrSges(3,3)) * t115 + (t141 * t121 + (-m(4) * pkin(2) - mrSges(4,1) * t113 + mrSges(4,2) * t110 - mrSges(3,1)) * t119) * t112 + t139 * (t104 * t133 - t115 * t105)) * g(3) + (-mrSges(1,2) - t111 * mrSges(2,1) - t114 * mrSges(2,2) - m(3) * t125 - t91 * mrSges(3,1) + mrSges(3,3) * t134 - m(4) * (pkin(2) * t91 + t125) - (-t110 * t134 + t113 * t91) * mrSges(4,1) - (-t110 * t91 - t113 * t134) * mrSges(4,2) - m(5) * t122 - t79 * mrSges(5,1) - m(6) * (pkin(4) * t79 + t122) - m(7) * (pkin(5) * t137 + t103 * t79 + t122) + t143 * (t120 * t79 + t137) + t142 * (-t118 * t79 + t120 * t90) + t140 * t90 + t139 * (t104 * t91 + t105 * t134)) * g(2) + (-mrSges(1,1) - t114 * mrSges(2,1) + t111 * mrSges(2,2) - m(3) * t129 - t93 * mrSges(3,1) - mrSges(3,3) * t135 - m(4) * (pkin(2) * t93 + t129) - (t113 * t93 + t127) * mrSges(4,1) - (-t110 * t93 + t113 * t135) * mrSges(4,2) - m(5) * t124 - t81 * mrSges(5,1) - m(6) * (pkin(4) * t81 + t124) - m(7) * (pkin(5) * t136 + t103 * t81 + t124) + t143 * (t120 * t81 + t136) + t142 * (-t118 * t81 + t120 * t92) + t140 * t92 + t139 * (t104 * t93 - t105 * t135)) * g(1);
U  = t1;

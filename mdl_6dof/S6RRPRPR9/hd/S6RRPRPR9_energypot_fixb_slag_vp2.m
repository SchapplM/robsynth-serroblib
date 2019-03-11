% Calculate potential energy for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_energypot_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:55:19
% EndTime: 2019-03-09 10:55:20
% DurationCPUTime: 0.90s
% Computational Cost: add. (355->108), mult. (554->123), div. (0->0), fcn. (633->14), ass. (0->44)
t106 = pkin(12) + qJ(6);
t101 = cos(t106);
t108 = sin(pkin(12));
t111 = cos(pkin(12));
t115 = -pkin(9) - qJ(3);
t143 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t144 = -m(5) - m(6) - m(7);
t99 = sin(t106);
t145 = -t99 * mrSges(7,1) - t111 * mrSges(6,2) - t101 * mrSges(7,2) - t144 * t115 - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t108 - t143;
t141 = -m(6) * qJ(5) + m(7) * (-pkin(10) - qJ(5)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t97 = pkin(5) * t111 + pkin(4);
t139 = -m(6) * pkin(4) - m(7) * t97 - mrSges(6,1) * t111 - mrSges(7,1) * t101 + mrSges(6,2) * t108 + mrSges(7,2) * t99 - mrSges(5,1);
t109 = sin(pkin(11));
t138 = pkin(3) * t109;
t113 = cos(pkin(6));
t137 = pkin(8) * t113 + pkin(7);
t119 = cos(qJ(1));
t110 = sin(pkin(6));
t117 = sin(qJ(1));
t132 = t110 * t117;
t134 = pkin(1) * t119 + pkin(8) * t132;
t116 = sin(qJ(2));
t133 = t110 * t116;
t118 = cos(qJ(2));
t131 = t110 * t118;
t130 = t110 * t119;
t129 = t116 * t117;
t128 = t116 * t119;
t127 = t117 * t118;
t126 = t118 * t119;
t125 = t108 * t131;
t124 = t109 * t132;
t112 = cos(pkin(11));
t98 = pkin(3) * t112 + pkin(2);
t123 = t113 * t138 + t115 * t131 + t133 * t98 + t137;
t104 = t117 * pkin(1);
t122 = -pkin(8) * t130 + t104;
t107 = pkin(11) + qJ(4);
t102 = cos(t107);
t100 = sin(t107);
t88 = -t113 * t129 + t126;
t86 = t113 * t128 + t127;
t82 = t100 * t113 + t102 * t133;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(5) * t123 - t82 * mrSges(5,1) + mrSges(5,3) * t131 - m(6) * (pkin(4) * t82 + t123) - (t111 * t82 - t125) * mrSges(6,1) - (-t108 * t82 - t111 * t131) * mrSges(6,2) - m(7) * (-pkin(5) * t125 + t82 * t97 + t123) - (t101 * t82 - t131 * t99) * mrSges(7,1) - (-t101 * t131 - t82 * t99) * mrSges(7,2) + (-m(3) - m(4)) * t137 + (-t109 * mrSges(4,1) - t112 * mrSges(4,2) - mrSges(3,3)) * t113 + (t143 * t118 + (-m(4) * pkin(2) - mrSges(4,1) * t112 + mrSges(4,2) * t109 - mrSges(3,1)) * t116) * t110 + t141 * (t100 * t133 - t102 * t113)) * g(3) + (-mrSges(1,2) - t117 * mrSges(2,1) - t119 * mrSges(2,2) - m(3) * t122 - t86 * mrSges(3,1) + mrSges(3,3) * t130 - m(4) * (pkin(2) * t86 + t122) - (-t109 * t130 + t112 * t86) * mrSges(4,1) - (-t109 * t86 - t112 * t130) * mrSges(4,2) + t144 * (t104 + t86 * t98 + (-pkin(8) - t138) * t130) + t139 * (-t100 * t130 + t102 * t86) + t141 * (t100 * t86 + t102 * t130) + t145 * (-t113 * t126 + t129)) * g(2) + (-mrSges(1,1) - t119 * mrSges(2,1) + t117 * mrSges(2,2) - m(3) * t134 - t88 * mrSges(3,1) - mrSges(3,3) * t132 - m(4) * (pkin(2) * t88 + t134) - (t112 * t88 + t124) * mrSges(4,1) - (-t109 * t88 + t112 * t132) * mrSges(4,2) + t144 * (pkin(3) * t124 + t88 * t98 + t134) + t139 * (t100 * t132 + t102 * t88) + t141 * (t100 * t88 - t102 * t132) + t145 * (t113 * t127 + t128)) * g(1);
U  = t1;

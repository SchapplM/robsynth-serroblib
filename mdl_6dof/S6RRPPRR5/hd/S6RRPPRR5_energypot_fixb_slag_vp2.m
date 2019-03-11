% Calculate potential energy for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:07
% EndTime: 2019-03-09 09:07:08
% DurationCPUTime: 0.64s
% Computational Cost: add. (232->88), mult. (478->97), div. (0->0), fcn. (530->10), ass. (0->44)
t141 = -m(6) - m(7);
t142 = pkin(9) - qJ(3);
t101 = sin(pkin(6));
t110 = cos(qJ(1));
t124 = t110 * t101;
t102 = cos(pkin(6));
t105 = sin(qJ(2));
t123 = t110 * t105;
t106 = sin(qJ(1));
t109 = cos(qJ(2));
t125 = t106 * t109;
t87 = t102 * t123 + t125;
t140 = t87 * pkin(3) + qJ(4) * t124;
t139 = -mrSges(3,1) - mrSges(4,1) - mrSges(5,1);
t138 = -mrSges(3,2) + mrSges(5,2) + mrSges(4,3);
t137 = -mrSges(5,3) + mrSges(4,2) + mrSges(3,3);
t136 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t103 = sin(qJ(6));
t107 = cos(qJ(6));
t135 = -m(7) * pkin(5) - t107 * mrSges(7,1) + t103 * mrSges(7,2) - mrSges(6,1);
t134 = t103 * mrSges(7,1) + t107 * mrSges(7,2) - t141 * t142 + mrSges(6,3) - t138;
t133 = t102 * pkin(8) + pkin(7);
t88 = t102 * t125 + t123;
t131 = t88 * qJ(3);
t127 = t106 * t101;
t130 = t110 * pkin(1) + pkin(8) * t127;
t129 = t101 * t105;
t128 = t101 * t109;
t126 = t106 * t105;
t122 = t110 * t109;
t121 = pkin(2) * t129 + t133;
t89 = -t102 * t126 + t122;
t120 = t89 * pkin(2) + t130;
t119 = t106 * pkin(1) - pkin(8) * t124;
t117 = t87 * pkin(2) + t119;
t116 = pkin(3) * t129 - t102 * qJ(4) + t121;
t115 = t89 * pkin(3) - qJ(4) * t127 + t120;
t86 = -t102 * t122 + t126;
t114 = t86 * qJ(3) + t117;
t111 = pkin(4) * t129 + t142 * t128 + t116;
t108 = cos(qJ(5));
t104 = sin(qJ(5));
t85 = -t102 * t104 + t108 * t129;
t1 = (-mrSges(1,3) - m(2) * pkin(7) - mrSges(2,3) - m(3) * t133 - m(4) * t121 - m(5) * t116 - m(6) * t111 - t85 * mrSges(6,1) - mrSges(6,3) * t128 - m(7) * (t85 * pkin(5) + t111) - (t103 * t128 + t85 * t107) * mrSges(7,1) - (-t85 * t103 + t107 * t128) * mrSges(7,2) + t136 * (t102 * t108 + t104 * t129) - t137 * t102 + (((m(4) + m(5)) * qJ(3) + t138) * t109 + t139 * t105) * t101) * g(3) + (-mrSges(1,2) - t106 * mrSges(2,1) - t110 * mrSges(2,2) - m(3) * t119 - m(4) * t114 - m(5) * (t114 + t140) + t141 * (t87 * pkin(4) + t117 + t140) + t136 * (t87 * t104 - t108 * t124) + t139 * t87 + t135 * (t104 * t124 + t87 * t108) + t137 * t124 + t134 * t86) * g(2) + (-mrSges(1,1) - t110 * mrSges(2,1) + t106 * mrSges(2,2) - m(3) * t130 - m(4) * (t120 + t131) - m(5) * (t115 + t131) + t141 * (t89 * pkin(4) + t115) + t136 * (t89 * t104 + t108 * t127) + t139 * t89 + t135 * (-t104 * t127 + t89 * t108) - t137 * t127 + t134 * t88) * g(1);
U  = t1;

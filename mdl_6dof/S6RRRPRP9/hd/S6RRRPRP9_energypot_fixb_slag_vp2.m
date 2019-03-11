% Calculate potential energy for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP9_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:47
% EndTime: 2019-03-09 17:21:47
% DurationCPUTime: 0.51s
% Computational Cost: add. (203->78), mult. (401->87), div. (0->0), fcn. (430->8), ass. (0->35)
t130 = -m(6) - m(7);
t129 = -mrSges(4,1) - mrSges(5,1);
t128 = mrSges(2,2) - mrSges(3,3);
t127 = mrSges(4,2) - mrSges(5,3);
t103 = cos(qJ(2));
t99 = sin(qJ(2));
t126 = -t103 * mrSges(3,1) + t99 * mrSges(3,2) - mrSges(2,1);
t125 = -mrSges(4,3) - mrSges(5,2) + mrSges(6,3) + mrSges(7,2);
t124 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t123 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t98 = sin(qJ(3));
t122 = t98 * t99;
t100 = sin(qJ(1));
t104 = cos(qJ(1));
t121 = t104 * pkin(1) + t100 * pkin(7);
t120 = t100 * t99;
t102 = cos(qJ(3));
t119 = t102 * t99;
t118 = t104 * t99;
t117 = t100 * t103;
t116 = t103 * t104;
t115 = t100 * pkin(1) - pkin(7) * t104;
t114 = pkin(2) * t116 + pkin(8) * t118 + t121;
t113 = t99 * pkin(2) - pkin(8) * t103 + pkin(6);
t111 = pkin(2) * t117 + pkin(8) * t120 + t115;
t110 = pkin(3) * t119 + qJ(4) * t122 + t113;
t81 = -t100 * t102 + t98 * t116;
t82 = t100 * t98 + t102 * t116;
t109 = t82 * pkin(3) + t81 * qJ(4) + t114;
t79 = t102 * t104 + t98 * t117;
t80 = t102 * t117 - t104 * t98;
t107 = t80 * pkin(3) + t79 * qJ(4) + t111;
t101 = cos(qJ(5));
t97 = sin(qJ(5));
t1 = (-m(4) * t113 - m(5) * t110 - mrSges(1,3) - mrSges(2,3) + t130 * (pkin(4) * t119 + t103 * pkin(9) + t110) + t123 * (-t101 * t122 + t97 * t119) + (-m(2) - m(3)) * pkin(6) + (-mrSges(3,2) - t125) * t103 + (t124 * (t101 * t102 + t97 * t98) + t129 * t102 + t127 * t98 - mrSges(3,1)) * t99) * g(3) + (-m(3) * t115 - m(4) * t111 - m(5) * t107 - mrSges(1,2) + t129 * t80 + t127 * t79 + t130 * (t80 * pkin(4) - pkin(9) * t120 + t107) + t124 * (t101 * t80 + t79 * t97) + t123 * (-t79 * t101 + t80 * t97) - t128 * t104 + t126 * t100 + t125 * t120) * g(2) + (-m(3) * t121 - m(4) * t114 - m(5) * t109 - mrSges(1,1) + t129 * t82 + t127 * t81 + t130 * (t82 * pkin(4) - pkin(9) * t118 + t109) + t124 * (t101 * t82 + t81 * t97) + t123 * (-t81 * t101 + t82 * t97) + t126 * t104 + t128 * t100 + t125 * t118) * g(1);
U  = t1;

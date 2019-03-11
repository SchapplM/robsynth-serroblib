% Calculate potential energy for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:20
% EndTime: 2019-03-09 18:20:20
% DurationCPUTime: 0.47s
% Computational Cost: add. (206->77), mult. (213->71), div. (0->0), fcn. (184->10), ass. (0->35)
t113 = -mrSges(4,1) + mrSges(5,2) + m(7) * (-pkin(10) - pkin(9)) - mrSges(7,3);
t74 = qJ(5) + qJ(6);
t69 = sin(t74);
t71 = cos(t74);
t76 = sin(qJ(5));
t112 = -m(7) * pkin(5) * t76 - t69 * mrSges(7,1) - t71 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t75 = qJ(2) + qJ(3);
t72 = cos(t75);
t81 = cos(qJ(1));
t101 = t72 * t81;
t70 = sin(t75);
t95 = qJ(4) * t70;
t111 = pkin(3) * t101 + t81 * t95;
t110 = m(5) + m(6) + m(7);
t109 = -m(6) * pkin(9) - mrSges(6,3);
t77 = sin(qJ(2));
t80 = cos(qJ(2));
t106 = -m(3) * pkin(1) - t80 * mrSges(3,1) + t77 * mrSges(3,2) + t112 * t70 + t113 * t72 - mrSges(2,1);
t105 = -m(3) * pkin(7) - t71 * mrSges(7,1) + t69 * mrSges(7,2) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t103 = t77 * pkin(2) + pkin(6);
t78 = sin(qJ(1));
t102 = t72 * t78;
t100 = t76 * t78;
t99 = t76 * t81;
t79 = cos(qJ(5));
t98 = t78 * t79;
t97 = t79 * t81;
t67 = pkin(2) * t80 + pkin(1);
t83 = -pkin(8) - pkin(7);
t96 = t78 * t67 + t81 * t83;
t62 = t81 * t67;
t93 = t62 + t111;
t91 = -t78 * t83 + t62;
t66 = pkin(5) * t79 + pkin(4);
t1 = (-m(4) * t103 - t77 * mrSges(3,1) - t80 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t110 * (t70 * pkin(3) + t103) + (t76 * mrSges(6,1) + t79 * mrSges(6,2) + t110 * qJ(4) - t112) * t72 + (t109 + t113) * t70) * g(3) + (-mrSges(1,2) - m(4) * t96 - (t100 * t70 - t97) * mrSges(6,1) - (t70 * t98 + t99) * mrSges(6,2) + t109 * t102 - t110 * (pkin(3) * t102 + t78 * t95 + t96) + (m(6) * pkin(4) + m(7) * t66 - t105) * t81 + t106 * t78) * g(2) + (-mrSges(1,1) - m(4) * t91 - m(5) * (t91 + t111) - m(6) * (pkin(9) * t101 + t93) - (t70 * t99 + t98) * mrSges(6,1) - (t70 * t97 - t100) * mrSges(6,2) - mrSges(6,3) * t101 - m(7) * t93 + t106 * t81 + (-m(6) * (pkin(4) - t83) - m(7) * (t66 - t83) + t105) * t78) * g(1);
U  = t1;

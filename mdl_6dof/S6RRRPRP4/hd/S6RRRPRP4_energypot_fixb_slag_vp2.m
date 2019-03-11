% Calculate potential energy for
% S6RRRPRP4
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:29
% EndTime: 2019-03-09 16:43:29
% DurationCPUTime: 0.45s
% Computational Cost: add. (201->70), mult. (223->66), div. (0->0), fcn. (198->8), ass. (0->32)
t115 = -mrSges(4,1) + mrSges(5,2);
t114 = mrSges(4,2) - mrSges(5,3);
t113 = -m(6) - m(7);
t112 = -mrSges(6,3) - mrSges(7,2);
t81 = qJ(2) + qJ(3);
t77 = sin(t81);
t78 = cos(t81);
t83 = sin(qJ(2));
t86 = cos(qJ(2));
t111 = -m(3) * pkin(1) - t86 * mrSges(3,1) + t83 * mrSges(3,2) + t114 * t77 + t115 * t78 - mrSges(2,1);
t110 = -m(3) * pkin(7) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t109 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t108 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t107 = t83 * pkin(2) + pkin(6);
t82 = sin(qJ(5));
t84 = sin(qJ(1));
t106 = t82 * t84;
t105 = t84 * t78;
t85 = cos(qJ(5));
t104 = t84 * t85;
t87 = cos(qJ(1));
t103 = t87 * t77;
t102 = t87 * t78;
t75 = pkin(2) * t86 + pkin(1);
t88 = -pkin(8) - pkin(7);
t101 = t84 * t75 + t87 * t88;
t100 = qJ(4) * t77;
t99 = t77 * pkin(3) + t107;
t97 = t87 * t75 - t84 * t88;
t95 = pkin(3) * t105 + t84 * t100 + t101;
t92 = pkin(3) * t102 + t87 * t100 + t97;
t1 = (-m(4) * t107 - m(5) * t99 - mrSges(3,1) * t83 - mrSges(3,2) * t86 - mrSges(1,3) - mrSges(2,3) + t113 * (t77 * pkin(9) + t99) + (-m(2) - m(3)) * pkin(6) + (-t108 * t85 + t109 * t82 + (m(5) - t113) * qJ(4) - t114) * t78 + (t112 + t115) * t77) * g(3) + (-m(4) * t101 - m(5) * t95 - mrSges(1,2) + t113 * (-pkin(4) * t87 + pkin(9) * t105 + t95) - t109 * (t77 * t106 - t85 * t87) + t108 * (t77 * t104 + t82 * t87) + t112 * t105 - t110 * t87 + t111 * t84) * g(2) + (-m(4) * t97 - m(5) * t92 - mrSges(1,1) + t113 * (t84 * pkin(4) + pkin(9) * t102 + t92) - t109 * (t82 * t103 + t104) - t108 * (-t85 * t103 + t106) + t112 * t102 + t111 * t87 + t110 * t84) * g(1);
U  = t1;

% Calculate potential energy for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:10
% EndTime: 2019-03-09 08:17:11
% DurationCPUTime: 0.60s
% Computational Cost: add. (165->74), mult. (302->73), div. (0->0), fcn. (297->8), ass. (0->37)
t117 = -mrSges(3,1) + mrSges(4,2);
t116 = mrSges(3,2) - mrSges(4,3);
t115 = -m(6) - m(7);
t84 = sin(qJ(2));
t101 = qJ(3) * t84;
t85 = sin(qJ(1));
t87 = cos(qJ(2));
t104 = t85 * t87;
t114 = pkin(2) * t104 + t85 * t101;
t113 = t116 * t84 + t117 * t87 - mrSges(2,1);
t112 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t111 = m(7) * pkin(8) - mrSges(6,2) - mrSges(5,3) + mrSges(7,3);
t83 = sin(qJ(6));
t86 = cos(qJ(6));
t110 = t83 * mrSges(7,1) + t86 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t109 = m(7) * pkin(5) + t86 * mrSges(7,1) - t83 * mrSges(7,2) + mrSges(5,1) + mrSges(6,1);
t108 = t84 * pkin(2) + pkin(6);
t81 = sin(pkin(9));
t107 = t81 * t85;
t88 = cos(qJ(1));
t106 = t84 * t88;
t82 = cos(pkin(9));
t105 = t85 * t82;
t103 = t87 * t88;
t102 = t88 * pkin(1) + t85 * pkin(7);
t100 = qJ(4) * t87;
t99 = t84 * qJ(4) + t108;
t79 = t85 * pkin(1);
t98 = -pkin(7) * t88 + t79;
t95 = pkin(2) * t103 + t88 * t101 + t102;
t92 = t85 * pkin(3) + t88 * t100 + t95;
t91 = t85 * t100 + t79 + (-pkin(3) - pkin(7)) * t88 + t114;
t66 = t84 * t107 - t82 * t88;
t65 = t84 * t105 + t81 * t88;
t64 = t81 * t106 + t105;
t63 = -t82 * t106 + t107;
t1 = (-m(4) * t108 - m(5) * t99 - mrSges(1,3) - mrSges(2,3) + t115 * (t87 * t82 * qJ(5) + t99) + (-m(2) - m(3)) * pkin(6) + (t115 * (-pkin(4) * t81 - qJ(3)) - t110 * t82 + t109 * t81 + (m(4) + m(5)) * qJ(3) - t116) * t87 + (t111 + t117) * t84) * g(3) + (-mrSges(1,2) - m(3) * t98 - m(4) * (t98 + t114) - m(5) * t91 + t115 * (t66 * pkin(4) - t65 * qJ(5) + t91) - t109 * t66 + t110 * t65 - t112 * t88 + t113 * t85 + t111 * t104) * g(2) + (-m(3) * t102 - m(4) * t95 - m(5) * t92 - mrSges(1,1) + t115 * (t64 * pkin(4) + t63 * qJ(5) + t92) - t109 * t64 - t110 * t63 + t113 * t88 + t112 * t85 + t111 * t103) * g(1);
U  = t1;

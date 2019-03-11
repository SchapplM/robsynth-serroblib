% Calculate potential energy for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:25
% EndTime: 2019-03-09 10:03:25
% DurationCPUTime: 0.54s
% Computational Cost: add. (155->71), mult. (276->67), div. (0->0), fcn. (261->6), ass. (0->36)
t117 = mrSges(3,2) - mrSges(4,3);
t116 = m(7) * qJ(6) - mrSges(3,1) + mrSges(4,2);
t114 = -m(6) - m(7);
t83 = sin(qJ(1));
t85 = cos(qJ(2));
t102 = t83 * t85;
t82 = sin(qJ(2));
t98 = qJ(3) * t82;
t113 = pkin(2) * t102 + t83 * t98;
t112 = mrSges(5,1) + mrSges(6,1) + mrSges(7,1);
t111 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t110 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t109 = -mrSges(5,3) - mrSges(6,2) + mrSges(7,3);
t108 = -m(7) * pkin(5) - t112;
t107 = t116 * t85 + t117 * t82 - mrSges(2,1);
t106 = pkin(2) * t82 + pkin(6);
t81 = sin(qJ(4));
t105 = t81 * t83;
t86 = cos(qJ(1));
t104 = t81 * t86;
t84 = cos(qJ(4));
t103 = t83 * t84;
t101 = t84 * t86;
t100 = t85 * t86;
t99 = pkin(1) * t86 + pkin(7) * t83;
t96 = pkin(8) * t82 + t106;
t79 = t83 * pkin(1);
t95 = -pkin(7) * t86 + t79;
t93 = pkin(2) * t100 + t86 * t98 + t99;
t90 = pkin(3) * t83 + pkin(8) * t100 + t93;
t89 = pkin(8) * t102 + t79 + (-pkin(3) - pkin(7)) * t86 + t113;
t66 = t105 * t82 - t101;
t65 = t103 * t82 + t104;
t64 = t104 * t82 + t103;
t63 = -t101 * t82 + t105;
t1 = (-m(4) * t106 - m(5) * t96 - mrSges(1,3) - mrSges(2,3) + t114 * (qJ(5) * t84 * t85 + t96) + (-m(2) - m(3)) * pkin(6) + (t110 * t84 + (m(4) + m(5) - t114) * qJ(3) + (m(6) * pkin(4) - m(7) * (-pkin(4) - pkin(5)) + t112) * t81 - t117) * t85 + (t109 + t116) * t82) * g(3) + (-mrSges(1,2) - m(3) * t95 - m(4) * (t95 + t113) - m(5) * t89 + t114 * (pkin(4) * t66 - t65 * qJ(5) + t89) - t111 * t86 + t107 * t83 + t108 * t66 - t110 * t65 + t109 * t102) * g(2) + (-m(3) * t99 - m(4) * t93 - m(5) * t90 - mrSges(1,1) + t114 * (pkin(4) * t64 + qJ(5) * t63 + t90) + t107 * t86 + t111 * t83 + t108 * t64 + t110 * t63 + t109 * t100) * g(1);
U  = t1;

% Calculate potential energy for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:35
% EndTime: 2019-03-09 13:55:36
% DurationCPUTime: 0.55s
% Computational Cost: add. (188->74), mult. (331->67), div. (0->0), fcn. (338->10), ass. (0->35)
t125 = -mrSges(3,1) - mrSges(4,1);
t124 = mrSges(3,2) - mrSges(4,3);
t117 = -m(6) * pkin(9) + m(7) * (-pkin(10) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t90 = sin(qJ(4));
t123 = t117 * t90;
t91 = sin(qJ(2));
t94 = cos(qJ(4));
t112 = t91 * t94;
t88 = qJ(5) + qJ(6);
t80 = sin(t88);
t81 = cos(t88);
t89 = sin(qJ(5));
t93 = cos(qJ(5));
t115 = -m(6) * pkin(4) - m(7) * (pkin(5) * t93 + pkin(4)) - t93 * mrSges(6,1) - t81 * mrSges(7,1) + t89 * mrSges(6,2) + t80 * mrSges(7,2) - mrSges(5,1);
t95 = cos(qJ(2));
t69 = t90 * t91 + t94 * t95;
t122 = -t117 * t112 + t69 * t115 + t124 * t91 + t125 * t95 - mrSges(2,1);
t121 = -m(5) - m(6);
t108 = qJ(3) * t91;
t92 = sin(qJ(1));
t111 = t92 * t95;
t120 = pkin(2) * t111 + t92 * t108;
t96 = cos(qJ(1));
t119 = pkin(3) * t111 + t96 * pkin(8);
t116 = -t89 * mrSges(6,1) - t80 * mrSges(7,1) - t93 * mrSges(6,2) - t81 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3);
t113 = pkin(5) * t89;
t110 = t95 * t96;
t109 = t96 * pkin(1) + t92 * pkin(7);
t85 = t92 * pkin(1);
t107 = -t96 * pkin(7) + t85;
t106 = pkin(2) * t110 + t96 * t108 + t109;
t105 = t91 * pkin(2) - t95 * qJ(3) + pkin(6);
t104 = pkin(3) * t110 + t106;
t101 = t107 + t120;
t1 = (-m(4) * t105 - mrSges(1,3) - mrSges(2,3) - t124 * t95 + t125 * t91 + (-m(2) - m(3)) * pkin(6) + t115 * (-t90 * t95 + t112) + (-m(7) + t121) * (t91 * pkin(3) + t105) + t117 * t69) * g(3) + (-mrSges(1,2) - m(3) * t107 - m(4) * t101 - m(7) * (t85 + t119 + t120) + t121 * (t101 + t119) + t111 * t123 + (-m(7) * (-pkin(7) + t113) + t116) * t96 + t122 * t92) * g(2) + (-m(3) * t109 - m(4) * t106 - m(7) * t104 - mrSges(1,1) + t121 * (-t92 * pkin(8) + t104) + t110 * t123 + (-m(7) * (-pkin(8) - t113) - t116) * t92 + t122 * t96) * g(1);
U  = t1;

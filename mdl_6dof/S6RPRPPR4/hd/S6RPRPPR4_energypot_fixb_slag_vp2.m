% Calculate potential energy for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:41
% EndTime: 2019-03-09 02:46:41
% DurationCPUTime: 0.55s
% Computational Cost: add. (239->71), mult. (283->68), div. (0->0), fcn. (278->10), ass. (0->37)
t126 = -mrSges(4,2) + mrSges(6,2) + mrSges(5,3) - mrSges(7,3);
t86 = pkin(9) + qJ(3);
t83 = sin(t86);
t84 = cos(t86);
t125 = pkin(3) * t84 + qJ(4) * t83;
t88 = sin(pkin(9));
t90 = cos(pkin(9));
t124 = -m(3) * pkin(1) - t90 * mrSges(3,1) - t84 * mrSges(4,1) + t88 * mrSges(3,2) - mrSges(2,1) + (m(7) * pkin(8) - t126) * t83;
t123 = -m(6) - m(7);
t87 = sin(pkin(10));
t89 = cos(pkin(10));
t122 = (pkin(4) * t89 + qJ(5) * t87) * t83;
t119 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t92 = sin(qJ(6));
t94 = cos(qJ(6));
t117 = -t92 * mrSges(7,1) - t94 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t116 = -m(7) * pkin(5) - t94 * mrSges(7,1) + t92 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t114 = t88 * pkin(2) + pkin(6);
t93 = sin(qJ(1));
t111 = t93 * t87;
t110 = t93 * t89;
t95 = cos(qJ(1));
t109 = t95 * t87;
t108 = t95 * t89;
t81 = t90 * pkin(2) + pkin(1);
t91 = -pkin(7) - qJ(2);
t107 = t93 * t81 + t95 * t91;
t105 = t83 * pkin(3) + t114;
t103 = t95 * t81 - t93 * t91;
t102 = t125 * t93 + t107;
t100 = -t84 * qJ(4) + t105;
t99 = t125 * t95 + t103;
t70 = t84 * t108 + t111;
t69 = t84 * t109 - t110;
t68 = t84 * t110 - t109;
t67 = t84 * t111 + t108;
t1 = (-mrSges(1,3) - mrSges(2,3) - t88 * mrSges(3,1) - t90 * mrSges(3,2) - m(4) * t114 - m(5) * t100 - m(6) * (t100 + t122) - m(7) * (t105 + t122) + (-m(2) - m(3)) * pkin(6) + (-m(7) * (pkin(8) - qJ(4)) + t126) * t84 + (t116 * t89 + t117 * t87 - mrSges(4,1)) * t83) * g(3) + (-m(4) * t107 - m(5) * t102 - mrSges(1,2) + t123 * (t68 * pkin(4) + t67 * qJ(5) + t102) + t116 * t68 + t117 * t67 - t119 * t95 + t124 * t93) * g(2) + (-m(4) * t103 - m(5) * t99 - mrSges(1,1) + t123 * (t70 * pkin(4) + t69 * qJ(5) + t99) + t116 * t70 + t117 * t69 + t119 * t93 + t124 * t95) * g(1);
U  = t1;

% Calculate potential energy for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:35
% EndTime: 2019-03-09 11:11:36
% DurationCPUTime: 0.53s
% Computational Cost: add. (186->78), mult. (246->70), div. (0->0), fcn. (221->10), ass. (0->35)
t82 = -qJ(5) - pkin(8);
t121 = -mrSges(3,1) + mrSges(4,2) + m(7) * (-pkin(9) + t82) - mrSges(7,3) + m(6) * t82 - mrSges(6,3);
t83 = sin(qJ(4));
t109 = t83 * pkin(4);
t81 = qJ(4) + pkin(10);
t74 = qJ(6) + t81;
t69 = sin(t74);
t70 = cos(t74);
t72 = sin(t81);
t73 = cos(t81);
t120 = -m(6) * t109 - t72 * mrSges(6,1) - t73 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - m(7) * (pkin(5) * t72 + t109) - t69 * mrSges(7,1) - t70 * mrSges(7,2);
t84 = sin(qJ(2));
t100 = qJ(3) * t84;
t85 = sin(qJ(1));
t87 = cos(qJ(2));
t105 = t85 * t87;
t119 = pkin(2) * t105 + t85 * t100;
t118 = -m(5) - m(6) - m(7);
t117 = m(4) - t118;
t116 = -m(5) * pkin(8) - mrSges(5,3);
t111 = t120 * t84 + t121 * t87 - mrSges(2,1);
t110 = t73 * mrSges(6,1) + t70 * mrSges(7,1) - t72 * mrSges(6,2) - t69 * mrSges(7,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t86 = cos(qJ(4));
t71 = t86 * pkin(4) + pkin(3);
t107 = t85 * t83;
t106 = t85 * t86;
t88 = cos(qJ(1));
t104 = t88 * t83;
t103 = t88 * t86;
t102 = t88 * t87;
t101 = t88 * pkin(1) + t85 * pkin(7);
t77 = t85 * pkin(1);
t98 = -t88 * pkin(7) + t77;
t63 = pkin(5) * t73 + t71;
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t117 * (t84 * pkin(2) + pkin(6)) + (t83 * mrSges(5,1) + t86 * mrSges(5,2) + t117 * qJ(3) - t120) * t87 + (t116 + t121) * t84) * g(3) + (-mrSges(1,2) - m(3) * t98 - m(4) * (t98 + t119) - (t84 * t107 - t103) * mrSges(5,1) - (t84 * t106 + t104) * mrSges(5,2) + t116 * t105 + t118 * (t77 + t119) + (-m(5) * (-pkin(3) - pkin(7)) - m(6) * (-pkin(7) - t71) - m(7) * (-pkin(7) - t63) + t110) * t88 + t111 * t85) * g(2) + (-mrSges(1,1) - m(3) * t101 - (t84 * t104 + t106) * mrSges(5,1) - (t84 * t103 - t107) * mrSges(5,2) + t116 * t102 - t117 * (pkin(2) * t102 + t88 * t100 + t101) + t111 * t88 + (-m(5) * pkin(3) - m(6) * t71 - m(7) * t63 - t110) * t85) * g(1);
U  = t1;

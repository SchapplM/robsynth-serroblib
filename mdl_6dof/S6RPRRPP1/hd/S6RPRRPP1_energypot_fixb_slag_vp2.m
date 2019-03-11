% Calculate potential energy for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:27:49
% EndTime: 2019-03-09 04:27:49
% DurationCPUTime: 0.48s
% Computational Cost: add. (250->70), mult. (222->65), div. (0->0), fcn. (201->10), ass. (0->32)
t83 = sin(qJ(4));
t86 = cos(qJ(4));
t112 = -m(5) * pkin(3) - t86 * mrSges(5,1) + t83 * mrSges(5,2) - mrSges(4,1);
t111 = -m(5) * pkin(8) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t110 = m(4) + m(5);
t109 = -m(6) - m(7);
t108 = -t83 * mrSges(5,1) - t86 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t106 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t105 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t84 = sin(qJ(3));
t87 = cos(qJ(3));
t104 = t111 * t84 + t112 * t87 - mrSges(3,1);
t103 = pkin(4) * t83;
t85 = sin(qJ(1));
t77 = t85 * pkin(1);
t88 = cos(qJ(1));
t78 = t88 * pkin(1);
t80 = qJ(1) + pkin(9);
t74 = sin(t80);
t102 = t74 * t87;
t76 = cos(t80);
t101 = t76 * t87;
t81 = -qJ(5) - pkin(8);
t100 = t81 * t84;
t82 = qJ(2) + pkin(6);
t97 = t74 * pkin(2) + t77;
t96 = t76 * pkin(2) + t74 * pkin(7) + t78;
t79 = qJ(4) + pkin(10);
t75 = cos(t79);
t73 = sin(t79);
t72 = pkin(4) * t86 + pkin(3);
t1 = (-m(2) * pkin(6) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t109 * (t84 * t72 + t87 * t81 + t82) + (-m(3) - t110) * t82 - t111 * t87 + (t105 * t73 + t106 * t75 + t112) * t84) * g(3) + (-m(3) * t77 - t85 * mrSges(2,1) - mrSges(2,2) * t88 - mrSges(1,2) - t110 * t97 + t109 * (-t74 * t100 + t72 * t102 + (-pkin(7) - t103) * t76 + t97) + t106 * (t75 * t102 - t73 * t76) + t105 * (t73 * t102 + t75 * t76) + (t110 * pkin(7) - t108) * t76 + t104 * t74) * g(2) + (-m(3) * t78 - mrSges(2,1) * t88 + t85 * mrSges(2,2) - mrSges(1,1) - t110 * t96 + t109 * (-t76 * t100 + t72 * t101 + t74 * t103 + t96) + t106 * (t75 * t101 + t73 * t74) + t105 * (t73 * t101 - t74 * t75) + t108 * t74 + t104 * t76) * g(1);
U  = t1;

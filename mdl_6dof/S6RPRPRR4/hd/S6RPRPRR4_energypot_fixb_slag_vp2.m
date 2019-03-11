% Calculate potential energy for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR4_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:20
% EndTime: 2019-03-09 03:44:20
% DurationCPUTime: 0.47s
% Computational Cost: add. (211->70), mult. (199->59), div. (0->0), fcn. (170->10), ass. (0->31)
t112 = -mrSges(4,1) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(9) - pkin(8)) - mrSges(7,3);
t75 = qJ(5) + qJ(6);
t69 = sin(t75);
t70 = cos(t75);
t77 = sin(qJ(5));
t80 = cos(qJ(5));
t111 = -t69 * mrSges(7,1) - t80 * mrSges(6,2) - t70 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t77;
t108 = m(6) * pkin(8);
t74 = qJ(1) + pkin(10);
t67 = sin(t74);
t81 = cos(qJ(3));
t100 = t67 * t81;
t78 = sin(qJ(3));
t94 = qJ(4) * t78;
t107 = pkin(3) * t100 + t67 * t94;
t106 = m(5) + m(6) + m(7);
t103 = t80 * mrSges(6,1) + t70 * mrSges(7,1) - t77 * mrSges(6,2) - t69 * mrSges(7,2) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t102 = t111 * t78 + t112 * t81 - mrSges(3,1);
t79 = sin(qJ(1));
t72 = t79 * pkin(1);
t82 = cos(qJ(1));
t73 = t82 * pkin(1);
t68 = cos(t74);
t99 = t68 * t81;
t76 = qJ(2) + pkin(6);
t95 = t67 * pkin(2) + t72;
t92 = t68 * pkin(2) + t67 * pkin(7) + t73;
t91 = t95 + t107;
t90 = -t68 * pkin(7) + t95;
t66 = t80 * pkin(5) + pkin(4);
t1 = (-m(2) * pkin(6) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - m(4)) * t76 - t106 * (t78 * pkin(3) + t76) + (t106 * qJ(4) - t111) * t81 + (-t108 + t112) * t78) * g(3) + (-mrSges(1,2) - t79 * mrSges(2,1) - t82 * mrSges(2,2) - m(3) * t72 - m(4) * t90 - m(5) * (t90 + t107) - m(6) * (pkin(8) * t100 + t91) - m(7) * t91 + (-m(6) * (-pkin(4) - pkin(7)) - m(7) * (-pkin(7) - t66) + t103) * t68 + t102 * t67) * g(2) + (-t99 * t108 - m(3) * t73 - m(4) * t92 - t82 * mrSges(2,1) + t79 * mrSges(2,2) - mrSges(1,1) - t106 * (pkin(3) * t99 + t68 * t94 + t92) + (-m(6) * pkin(4) - m(7) * t66 - t103) * t67 + t102 * t68) * g(1);
U  = t1;

% Calculate potential energy for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPP8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:53
% EndTime: 2019-03-09 04:53:54
% DurationCPUTime: 0.46s
% Computational Cost: add. (143->70), mult. (245->64), div. (0->0), fcn. (230->6), ass. (0->34)
t111 = mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t109 = -m(6) - m(7);
t110 = -m(5) + t109;
t108 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t107 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t106 = -m(7) * pkin(5) - mrSges(7,1);
t78 = sin(qJ(3));
t81 = cos(qJ(3));
t105 = -t78 * mrSges(4,1) + t111 * t81 + mrSges(2,2) - mrSges(3,3);
t104 = m(7) * qJ(6) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t103 = pkin(2) + pkin(6);
t102 = pkin(3) * t78;
t77 = sin(qJ(4));
t82 = cos(qJ(1));
t101 = t77 * t82;
t79 = sin(qJ(1));
t100 = t79 * t77;
t80 = cos(qJ(4));
t99 = t79 * t80;
t98 = t79 * t81;
t95 = t82 * t80;
t73 = t79 * pkin(1);
t94 = t79 * pkin(7) + t73;
t93 = t82 * pkin(1) + t79 * qJ(2);
t92 = pkin(8) * t98;
t91 = t82 * t81 * pkin(8) + t94;
t90 = t82 * pkin(7) + t93;
t87 = t79 * t102 + t90;
t60 = t78 * t100 - t95;
t61 = t78 * t99 + t101;
t83 = t61 * pkin(4) + qJ(5) * t60 + t87;
t63 = t78 * t95 - t100;
t62 = t78 * t101 + t99;
t1 = (-m(4) * t103 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) + (t106 - t111) * t78 + t110 * (t81 * pkin(3) + t78 * pkin(8) + t103) + (t109 * (pkin(4) * t80 + qJ(5) * t77) - t104 * t80 + t107 * t77 - mrSges(4,1)) * t81) * g(3) + (-m(3) * t73 - m(4) * t94 - m(5) * t91 - mrSges(1,2) + t109 * (-t63 * pkin(4) - t62 * qJ(5) + t91) + t108 * t79 + t104 * t63 - t107 * t62 + (t106 * t81 + t110 * (-qJ(2) - t102) + (m(3) + m(4)) * qJ(2) - t105) * t82) * g(2) + (-mrSges(1,1) - m(3) * t93 - m(4) * t90 - m(5) * (t87 - t92) - m(6) * (t83 - t92) - m(7) * t83 - (m(7) * (-pkin(5) - pkin(8)) - mrSges(7,1)) * t98 + t108 * t82 - t104 * t61 + t107 * t60 + t105 * t79) * g(1);
U  = t1;

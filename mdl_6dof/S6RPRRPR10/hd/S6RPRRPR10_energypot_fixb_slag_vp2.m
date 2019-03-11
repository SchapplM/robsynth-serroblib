% Calculate potential energy for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:17
% EndTime: 2019-03-09 05:34:18
% DurationCPUTime: 0.54s
% Computational Cost: add. (153->73), mult. (271->68), div. (0->0), fcn. (266->8), ass. (0->36)
t113 = -mrSges(4,2) + mrSges(6,2) + mrSges(5,3);
t111 = -m(6) - m(7);
t112 = -m(5) + t111;
t110 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t109 = m(7) * pkin(9) + mrSges(7,3);
t79 = sin(qJ(3));
t83 = cos(qJ(3));
t108 = -t79 * mrSges(4,1) + t113 * t83 + mrSges(2,2) - mrSges(3,3);
t77 = sin(qJ(6));
t81 = cos(qJ(6));
t107 = t77 * mrSges(7,1) + t81 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t106 = -m(7) * pkin(5) - mrSges(7,1) * t81 + mrSges(7,2) * t77 - mrSges(5,1) - mrSges(6,1);
t105 = pkin(2) + pkin(6);
t104 = pkin(3) * t79;
t78 = sin(qJ(4));
t84 = cos(qJ(1));
t103 = t78 * t84;
t80 = sin(qJ(1));
t102 = t80 * t78;
t82 = cos(qJ(4));
t101 = t80 * t82;
t100 = t80 * t83;
t97 = t84 * t82;
t73 = t80 * pkin(1);
t96 = pkin(7) * t80 + t73;
t95 = pkin(1) * t84 + qJ(2) * t80;
t94 = pkin(8) * t100;
t93 = pkin(8) * t83 * t84 + t96;
t92 = pkin(7) * t84 + t95;
t89 = t104 * t80 + t92;
t61 = t102 * t79 - t97;
t62 = t101 * t79 + t103;
t85 = pkin(4) * t62 + qJ(5) * t61 + t89;
t64 = -t79 * t97 + t102;
t63 = t103 * t79 + t101;
t1 = (-m(4) * t105 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) + (t109 - t113) * t79 + t112 * (pkin(3) * t83 + pkin(8) * t79 + t105) + (t111 * (pkin(4) * t82 + qJ(5) * t78) + t106 * t82 - t107 * t78 - mrSges(4,1)) * t83) * g(3) + (-m(3) * t73 - m(4) * t96 - m(5) * t93 - mrSges(1,2) + t111 * (t64 * pkin(4) - t63 * qJ(5) + t93) + t106 * t64 + t107 * t63 + t110 * t80 + (t109 * t83 + t112 * (-qJ(2) - t104) + (m(3) + m(4)) * qJ(2) - t108) * t84) * g(2) + (-mrSges(1,1) - m(3) * t95 - m(4) * t92 - m(5) * (t89 - t94) - m(6) * (t85 - t94) - m(7) * t85 - (m(7) * (-pkin(8) + pkin(9)) + mrSges(7,3)) * t100 + t106 * t62 - t107 * t61 + t110 * t84 + t108 * t80) * g(1);
U  = t1;

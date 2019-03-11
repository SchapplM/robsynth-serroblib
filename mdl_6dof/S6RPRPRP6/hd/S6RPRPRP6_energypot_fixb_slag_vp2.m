% Calculate potential energy for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:18:03
% EndTime: 2019-03-09 03:18:03
% DurationCPUTime: 0.42s
% Computational Cost: add. (196->72), mult. (213->65), div. (0->0), fcn. (184->8), ass. (0->33)
t81 = sin(qJ(5));
t113 = -m(7) * pkin(5) * t81 + mrSges(4,2) - mrSges(5,3);
t112 = m(7) * (-qJ(6) - pkin(8)) - mrSges(4,1) + mrSges(5,2) - mrSges(7,3);
t111 = -m(5) - m(7);
t110 = mrSges(6,1) + mrSges(7,1);
t109 = mrSges(6,2) + mrSges(7,2);
t76 = pkin(9) + qJ(3);
t74 = cos(t76);
t84 = cos(qJ(1));
t100 = t74 * t84;
t73 = sin(t76);
t94 = qJ(4) * t73;
t108 = pkin(3) * t100 + t84 * t94;
t107 = m(6) - t111;
t106 = -m(6) * pkin(8) - mrSges(6,3);
t77 = sin(pkin(9));
t78 = cos(pkin(9));
t105 = -m(3) * pkin(1) - t78 * mrSges(3,1) + t77 * mrSges(3,2) + t112 * t74 + t113 * t73 - mrSges(2,1);
t83 = cos(qJ(5));
t104 = -m(7) * (pkin(5) * t83 + pkin(4)) - mrSges(5,1) + mrSges(2,2) - mrSges(4,3) - m(3) * qJ(2) - mrSges(3,3);
t102 = t77 * pkin(2) + pkin(6);
t82 = sin(qJ(1));
t101 = t74 * t82;
t99 = t81 * t84;
t98 = t82 * t81;
t97 = t82 * t83;
t96 = t83 * t84;
t70 = pkin(2) * t78 + pkin(1);
t80 = -pkin(7) - qJ(2);
t95 = t82 * t70 + t84 * t80;
t66 = t84 * t70;
t91 = -t82 * t80 + t66;
t1 = (-m(4) * t102 - t77 * mrSges(3,1) - t78 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t107 * (t73 * pkin(3) + t102) + (t107 * qJ(4) + t109 * t83 + t110 * t81 - t113) * t74 + (t106 + t112) * t73) * g(3) + (-m(4) * t95 - mrSges(1,2) + t106 * t101 - t107 * (pkin(3) * t101 + t82 * t94 + t95) - t110 * (t73 * t98 - t96) - t109 * (t73 * t97 + t99) + (m(6) * pkin(4) - t104) * t84 + t105 * t82) * g(2) + (-mrSges(1,1) - m(4) * t91 - m(6) * (pkin(8) * t100 + t108 + t66) - mrSges(6,3) * t100 + t111 * (t91 + t108) - t110 * (t73 * t99 + t97) - t109 * (t73 * t96 - t98) + (-m(6) * (pkin(4) - t80) + t104) * t82 + t105 * t84) * g(1);
U  = t1;

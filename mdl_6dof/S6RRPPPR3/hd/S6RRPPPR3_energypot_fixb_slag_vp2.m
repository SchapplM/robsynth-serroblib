% Calculate potential energy for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:46
% EndTime: 2019-03-09 08:13:46
% DurationCPUTime: 0.43s
% Computational Cost: add. (152->73), mult. (235->59), div. (0->0), fcn. (206->8), ass. (0->31)
t74 = pkin(9) + qJ(6);
t66 = sin(t74);
t67 = cos(t74);
t75 = sin(pkin(9));
t76 = cos(pkin(9));
t113 = t76 * mrSges(6,1) + t67 * mrSges(7,1) - t75 * mrSges(6,2) - t66 * mrSges(7,2) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t112 = -m(6) * qJ(5) - mrSges(3,1) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(8) - qJ(5)) - mrSges(7,3);
t111 = -m(6) - m(7);
t79 = sin(qJ(1));
t80 = cos(qJ(2));
t101 = t79 * t80;
t78 = sin(qJ(2));
t98 = qJ(3) * t78;
t110 = pkin(2) * t101 + t79 * t98;
t81 = cos(qJ(1));
t109 = pkin(3) * t101 + t81 * qJ(4);
t108 = m(5) - t111;
t65 = pkin(5) * t76 + pkin(4);
t105 = -mrSges(2,1) + t112 * t80 + (-m(6) * pkin(4) - m(7) * t65 - t113) * t78;
t104 = -t75 * mrSges(6,1) - t66 * mrSges(7,1) - t76 * mrSges(6,2) - t67 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3);
t103 = pkin(5) * t75;
t102 = t78 * pkin(2) + pkin(6);
t100 = t80 * t81;
t99 = t81 * pkin(1) + t79 * pkin(7);
t72 = t79 * pkin(1);
t96 = -pkin(7) * t81 + t72;
t95 = pkin(2) * t100 + t81 * t98 + t99;
t94 = -qJ(3) * t80 + t102;
t84 = t96 + t110;
t69 = t78 * pkin(3);
t1 = (-mrSges(1,3) - mrSges(2,3) - m(4) * t94 - m(5) * (t69 + t94) + t111 * (t69 + t102) + (-m(2) - m(3)) * pkin(6) + (-m(6) * (-pkin(4) - qJ(3)) - m(7) * (-qJ(3) - t65) + t113) * t80 + t112 * t78) * g(3) + (-mrSges(1,2) - m(3) * t96 - m(4) * t84 - m(5) * (t84 + t109) + t111 * (t72 + t109 + t110) + (m(6) * pkin(7) - m(7) * (-pkin(7) + t103) + t104) * t81 + t105 * t79) * g(2) + (-m(3) * t99 - m(4) * t95 - mrSges(1,1) - t108 * (pkin(3) * t100 + t95) + t105 * t81 + (m(7) * t103 + t108 * qJ(4) - t104) * t79) * g(1);
U  = t1;

% Calculate potential energy for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:31
% EndTime: 2019-03-09 05:55:31
% DurationCPUTime: 0.37s
% Computational Cost: add. (238->70), mult. (194->68), div. (0->0), fcn. (169->10), ass. (0->33)
t113 = -m(3) - m(4);
t112 = -m(6) - m(7);
t111 = -mrSges(6,3) - mrSges(7,2);
t85 = qJ(3) + qJ(4);
t79 = sin(t85);
t80 = cos(t85);
t88 = sin(qJ(3));
t91 = cos(qJ(3));
t110 = -m(4) * pkin(2) - mrSges(4,1) * t91 - mrSges(5,1) * t80 + mrSges(4,2) * t88 + mrSges(5,2) * t79 - mrSges(3,1);
t109 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t108 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t107 = m(4) * pkin(7) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t106 = pkin(4) * t80;
t89 = sin(qJ(1));
t82 = t89 * pkin(1);
t92 = cos(qJ(1));
t83 = t92 * pkin(1);
t84 = qJ(1) + pkin(10);
t77 = sin(t84);
t105 = t77 * t79;
t78 = cos(t84);
t104 = t78 * t79;
t87 = sin(qJ(5));
t103 = t80 * t87;
t90 = cos(qJ(5));
t102 = t80 * t90;
t86 = qJ(2) + pkin(6);
t76 = pkin(3) * t91 + pkin(2);
t93 = -pkin(8) - pkin(7);
t101 = t77 * t76 + t78 * t93 + t82;
t100 = t88 * pkin(3) + t86;
t99 = t78 * t76 - t77 * t93 + t83;
t1 = (-m(2) * pkin(6) - m(5) * t100 - t88 * mrSges(4,1) - t91 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t112 * (t79 * pkin(4) - pkin(9) * t80 + t100) + t113 * t86 + (-mrSges(5,2) - t111) * t80 + (t108 * t87 + t109 * t90 - mrSges(5,1)) * t79) * g(3) + (-m(5) * t101 - t89 * mrSges(2,1) - t92 * mrSges(2,2) - mrSges(1,2) + t112 * (pkin(9) * t105 + t77 * t106 + t101) + t113 * t82 + t109 * (t77 * t102 - t78 * t87) + t108 * (t77 * t103 + t78 * t90) + t111 * t105 + t107 * t78 + t110 * t77) * g(2) + (-m(5) * t99 - t92 * mrSges(2,1) + t89 * mrSges(2,2) - mrSges(1,1) + t112 * (pkin(9) * t104 + t78 * t106 + t99) + t113 * t83 + t109 * (t78 * t102 + t77 * t87) + t108 * (t78 * t103 - t77 * t90) + t111 * t104 + t110 * t78 - t107 * t77) * g(1);
U  = t1;

% Calculate potential energy for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:25:20
% EndTime: 2019-03-09 08:25:20
% DurationCPUTime: 0.48s
% Computational Cost: add. (239->73), mult. (236->70), div. (0->0), fcn. (215->10), ass. (0->36)
t84 = sin(pkin(10));
t85 = cos(pkin(10));
t119 = -m(5) * pkin(3) - t85 * mrSges(5,1) + t84 * mrSges(5,2) - mrSges(4,1);
t118 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3);
t117 = -m(4) - m(5);
t116 = -m(6) - m(7);
t115 = -mrSges(6,3) - mrSges(7,2);
t83 = qJ(2) + pkin(9);
t78 = sin(t83);
t80 = cos(t83);
t88 = sin(qJ(2));
t90 = cos(qJ(2));
t113 = -m(3) * pkin(1) - t90 * mrSges(3,1) + t88 * mrSges(3,2) + t118 * t78 + t119 * t80 - mrSges(2,1);
t112 = m(3) * pkin(7) + t84 * mrSges(5,1) + t85 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3);
t111 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t110 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t109 = pkin(4) * t84;
t108 = pkin(2) * t88 + pkin(6);
t91 = cos(qJ(1));
t107 = t80 * t91;
t82 = pkin(10) + qJ(5);
t77 = sin(t82);
t89 = sin(qJ(1));
t106 = t89 * t77;
t105 = t89 * t78;
t79 = cos(t82);
t104 = t89 * t79;
t103 = t91 * t78;
t76 = pkin(2) * t90 + pkin(1);
t86 = -qJ(3) - pkin(7);
t102 = t76 * t89 + t86 * t91;
t70 = t91 * t76;
t100 = -t89 * t86 + t70;
t87 = -pkin(8) - qJ(4);
t74 = pkin(4) * t85 + pkin(3);
t1 = (-mrSges(3,1) * t88 - mrSges(3,2) * t90 - mrSges(1,3) - mrSges(2,3) + t116 * (t74 * t78 + t80 * t87 + t108) + t117 * t108 + (-m(2) - m(3)) * pkin(6) + (-t115 - t118) * t80 + (t110 * t77 + t111 * t79 + t119) * t78) * g(3) + (-mrSges(1,2) + t116 * (t74 * t80 * t89 - t105 * t87 - t109 * t91 + t102) + t111 * (t104 * t80 - t77 * t91) + t110 * (t106 * t80 + t79 * t91) + t115 * t105 + t117 * t102 + t112 * t91 + t113 * t89) * g(2) + (-m(4) * t100 - m(5) * t70 - mrSges(1,1) + t116 * (-t103 * t87 + t107 * t74 + t109 * t89 + t100) + t111 * (t107 * t79 + t106) + t110 * (t107 * t77 - t104) + t115 * t103 + t113 * t91 + (m(5) * t86 - t112) * t89) * g(1);
U  = t1;

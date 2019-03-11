% Calculate potential energy for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:49:49
% EndTime: 2019-03-09 09:49:49
% DurationCPUTime: 0.38s
% Computational Cost: add. (223->72), mult. (257->69), div. (0->0), fcn. (242->8), ass. (0->37)
t114 = -m(6) - m(7);
t80 = qJ(2) + pkin(9);
t77 = sin(t80);
t82 = sin(qJ(4));
t85 = cos(qJ(4));
t113 = (pkin(4) * t85 + qJ(5) * t82) * t77;
t112 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t111 = mrSges(5,2) - mrSges(7,2) - mrSges(6,3);
t110 = -mrSges(5,3) - mrSges(6,2) + mrSges(7,3);
t109 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t78 = cos(t80);
t83 = sin(qJ(2));
t86 = cos(qJ(2));
t108 = -m(3) * pkin(1) - t86 * mrSges(3,1) - t78 * mrSges(4,1) + t83 * mrSges(3,2) - mrSges(2,1) + (m(7) * qJ(6) + mrSges(4,2)) * t77;
t107 = pkin(3) * t78;
t106 = t83 * pkin(2) + pkin(6);
t87 = cos(qJ(1));
t105 = t77 * t87;
t104 = t82 * t87;
t84 = sin(qJ(1));
t103 = t84 * t77;
t102 = t84 * t82;
t101 = t84 * t85;
t100 = t85 * t87;
t76 = pkin(2) * t86 + pkin(1);
t81 = -qJ(3) - pkin(7);
t99 = t84 * t76 + t87 * t81;
t97 = t77 * pkin(3) + t106;
t95 = t87 * t76 - t84 * t81;
t94 = pkin(8) * t103 + t84 * t107 + t99;
t92 = -pkin(8) * t78 + t97;
t91 = pkin(8) * t105 + t87 * t107 + t95;
t64 = t78 * t100 + t102;
t63 = t78 * t104 - t101;
t62 = t78 * t101 - t104;
t61 = t78 * t102 + t100;
t1 = (-mrSges(1,3) - mrSges(2,3) - t83 * mrSges(3,1) - t86 * mrSges(3,2) - m(4) * t106 - m(5) * t92 - m(6) * (t92 + t113) - m(7) * (t97 + t113) + (-m(2) - m(3)) * pkin(6) + (-mrSges(4,2) - m(7) * (-pkin(8) + qJ(6)) - t110) * t78 + (t109 * t85 + t111 * t82 - mrSges(4,1)) * t77) * g(3) + (-m(4) * t99 - m(5) * t94 - mrSges(1,2) + t114 * (t62 * pkin(4) + t61 * qJ(5) + t94) - t112 * t87 + t108 * t84 + t109 * t62 + t111 * t61 + t110 * t103) * g(2) + (-m(4) * t95 - m(5) * t91 - mrSges(1,1) + t114 * (t64 * pkin(4) + t63 * qJ(5) + t91) + t108 * t87 + t112 * t84 + t109 * t64 + t111 * t63 + t110 * t105) * g(1);
U  = t1;

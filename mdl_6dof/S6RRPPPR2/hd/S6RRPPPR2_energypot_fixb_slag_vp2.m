% Calculate potential energy for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:09:50
% EndTime: 2019-03-09 08:09:50
% DurationCPUTime: 0.50s
% Computational Cost: add. (206->76), mult. (213->69), div. (0->0), fcn. (184->10), ass. (0->34)
t76 = pkin(10) + qJ(6);
t71 = sin(t76);
t73 = cos(t76);
t78 = sin(pkin(10));
t117 = -m(7) * pkin(5) * t78 - t71 * mrSges(7,1) - t73 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t116 = -mrSges(4,1) + mrSges(5,2) + m(7) * (-pkin(8) - qJ(5)) - mrSges(7,3) - m(6) * qJ(5);
t114 = -m(6) - m(7);
t77 = qJ(2) + pkin(9);
t74 = cos(t77);
t85 = cos(qJ(1));
t104 = t74 * t85;
t72 = sin(t77);
t98 = qJ(4) * t72;
t113 = pkin(3) * t104 + t85 * t98;
t112 = m(5) - t114;
t109 = -m(3) * pkin(7) - t73 * mrSges(7,1) + t71 * mrSges(7,2) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t82 = sin(qJ(2));
t84 = cos(qJ(2));
t108 = -m(3) * pkin(1) - t84 * mrSges(3,1) + t82 * mrSges(3,2) + t116 * t74 + t117 * t72 - mrSges(2,1);
t106 = t82 * pkin(2) + pkin(6);
t83 = sin(qJ(1));
t105 = t74 * t83;
t103 = t78 * t85;
t79 = cos(pkin(10));
t102 = t79 * t85;
t101 = t83 * t78;
t100 = t83 * t79;
t70 = pkin(2) * t84 + pkin(1);
t80 = -qJ(3) - pkin(7);
t99 = t83 * t70 + t85 * t80;
t66 = t85 * t70;
t93 = -t83 * t80 + t66;
t68 = pkin(5) * t79 + pkin(4);
t1 = (-m(4) * t106 - t82 * mrSges(3,1) - t84 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t112 * (t72 * pkin(3) + t106) + (t78 * mrSges(6,1) + t79 * mrSges(6,2) + t112 * qJ(4) - t117) * t74 + (-mrSges(6,3) + t116) * t72) * g(3) + (-mrSges(1,2) - m(4) * t99 - (t72 * t101 - t102) * mrSges(6,1) - (t72 * t100 + t103) * mrSges(6,2) - mrSges(6,3) * t105 - t112 * (pkin(3) * t105 + t83 * t98 + t99) + (m(6) * pkin(4) + m(7) * t68 - t109) * t85 + t108 * t83) * g(2) + (-mrSges(1,1) - m(4) * t93 - m(5) * (t93 + t113) - (t72 * t103 + t100) * mrSges(6,1) - (t72 * t102 - t101) * mrSges(6,2) - mrSges(6,3) * t104 + t114 * (t66 + t113) + t108 * t85 + (-m(6) * (pkin(4) - t80) - m(7) * (t68 - t80) + t109) * t83) * g(1);
U  = t1;

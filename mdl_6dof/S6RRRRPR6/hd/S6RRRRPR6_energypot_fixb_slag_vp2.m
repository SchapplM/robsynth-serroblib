% Calculate potential energy for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:18:33
% EndTime: 2019-03-09 22:18:33
% DurationCPUTime: 0.41s
% Computational Cost: add. (245->61), mult. (249->50), div. (0->0), fcn. (228->12), ass. (0->27)
t83 = cos(qJ(3));
t68 = t83 * pkin(3) + pkin(2);
t79 = qJ(3) + qJ(4);
t71 = cos(t79);
t59 = pkin(4) * t71 + t68;
t69 = pkin(11) + t79;
t67 = qJ(6) + t69;
t61 = sin(t67);
t62 = cos(t67);
t63 = sin(t69);
t64 = cos(t69);
t70 = sin(t79);
t80 = sin(qJ(3));
t109 = -mrSges(3,1) - m(4) * pkin(2) - t83 * mrSges(4,1) + t80 * mrSges(4,2) - m(5) * t68 - t71 * mrSges(5,1) + t70 * mrSges(5,2) - m(6) * t59 - t64 * mrSges(6,1) + t63 * mrSges(6,2) - m(7) * (pkin(5) * t64 + t59) - t62 * mrSges(7,1) + t61 * mrSges(7,2);
t86 = -pkin(9) - pkin(8);
t78 = -qJ(5) + t86;
t108 = mrSges(3,2) + m(7) * (-pkin(10) + t78) - mrSges(7,3) + m(6) * t78 - mrSges(6,3) + m(5) * t86 - mrSges(5,3) - m(4) * pkin(8) - mrSges(4,3);
t103 = m(3) + m(4) + m(6) + m(7) + m(5);
t81 = sin(qJ(2));
t84 = cos(qJ(2));
t102 = t108 * t81 + t109 * t84 - mrSges(2,1);
t73 = t80 * pkin(3);
t60 = pkin(4) * t70 + t73;
t101 = m(5) * t73 + m(6) * t60 + m(7) * (pkin(5) * t63 + t60) - mrSges(2,2) + mrSges(3,3) + t61 * mrSges(7,1) + t62 * mrSges(7,2) + t63 * mrSges(6,1) + t64 * mrSges(6,2) + t70 * mrSges(5,1) + t71 * mrSges(5,2) + t80 * mrSges(4,1) + t83 * mrSges(4,2);
t85 = cos(qJ(1));
t82 = sin(qJ(1));
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) - t103) * pkin(6) - t108 * t84 + t109 * t81) * g(3) + (-mrSges(1,2) + (t103 * pkin(7) + t101) * t85 + (-t103 * pkin(1) + t102) * t82) * g(2) + (-mrSges(1,1) - t103 * (t85 * pkin(1) + t82 * pkin(7)) + t102 * t85 - t101 * t82) * g(1);
U  = t1;

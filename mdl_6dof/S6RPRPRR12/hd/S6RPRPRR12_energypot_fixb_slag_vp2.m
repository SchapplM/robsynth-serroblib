% Calculate potential energy for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRR12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:06
% EndTime: 2019-03-09 04:18:06
% DurationCPUTime: 0.43s
% Computational Cost: add. (138->68), mult. (201->54), div. (0->0), fcn. (172->8), ass. (0->24)
t89 = -mrSges(4,1) + mrSges(5,2);
t88 = mrSges(6,3) + mrSges(7,3);
t87 = m(5) + m(6) + m(7);
t59 = qJ(5) + qJ(6);
t50 = sin(t59);
t51 = cos(t59);
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t86 = -t50 * mrSges(7,1) - t63 * mrSges(6,2) - t51 * mrSges(7,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t60;
t61 = sin(qJ(3));
t64 = cos(qJ(3));
t85 = -t64 * mrSges(4,2) + t89 * t61 + mrSges(2,2) - mrSges(3,3);
t84 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,1) - m(6) * pkin(4) - t63 * mrSges(6,1) + t60 * mrSges(6,2) - m(7) * (pkin(5) * t63 + pkin(4)) - t51 * mrSges(7,1) + t50 * mrSges(7,2);
t66 = -pkin(9) - pkin(8);
t83 = -m(6) * pkin(8) + m(7) * t66 - t88;
t82 = pkin(2) + pkin(6);
t81 = pkin(3) * t61;
t62 = sin(qJ(1));
t55 = t62 * pkin(1);
t78 = t62 * pkin(7) + t55;
t65 = cos(qJ(1));
t77 = t65 * pkin(1) + t62 * qJ(2);
t75 = t65 * pkin(7) + t77;
t1 = (-m(4) * t82 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t87 * (t64 * pkin(3) + t61 * qJ(4) + t82) + (t83 + t89) * t64 + (mrSges(4,2) + t86) * t61) * g(3) + (-m(3) * t55 - m(4) * t78 - mrSges(1,2) - t87 * (t65 * t64 * qJ(4) + t78) + (m(5) * t81 + (-m(6) * (-pkin(3) - pkin(8)) - m(7) * (-pkin(3) + t66) + t88) * t61 + t86 * t64 + (m(3) + m(4) + t87) * qJ(2) - t85) * t65 + t84 * t62) * g(2) + (-m(3) * t77 - m(4) * t75 - mrSges(1,1) - t87 * (t62 * t81 + t75) + t84 * t65 + (t83 * t61 + (t87 * qJ(4) - t86) * t64 + t85) * t62) * g(1);
U  = t1;

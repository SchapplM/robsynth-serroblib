% Calculate potential energy for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:41
% EndTime: 2019-03-09 02:58:42
% DurationCPUTime: 0.47s
% Computational Cost: add. (121->66), mult. (188->51), div. (0->0), fcn. (155->6), ass. (0->24)
t87 = m(6) + m(7);
t86 = -mrSges(4,1) - mrSges(5,1);
t84 = mrSges(7,3) - mrSges(6,2);
t83 = m(5) + t87;
t57 = sin(qJ(6));
t60 = cos(qJ(6));
t82 = -m(7) * pkin(5) - t60 * mrSges(7,1) + t57 * mrSges(7,2) - mrSges(6,1) - mrSges(5,3);
t58 = sin(qJ(3));
t61 = cos(qJ(3));
t81 = -t61 * mrSges(4,2) + t58 * t86 + mrSges(2,2) - mrSges(3,3);
t80 = -m(7) * pkin(8) - t84;
t79 = t57 * mrSges(7,1) + t60 * mrSges(7,2) + qJ(5) * t87 - mrSges(2,1) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t78 = pkin(2) + pkin(6);
t77 = -pkin(3) - pkin(4);
t59 = sin(qJ(1));
t75 = t59 * t58;
t52 = t59 * pkin(1);
t74 = t59 * pkin(7) + t52;
t62 = cos(qJ(1));
t73 = t62 * pkin(1) + t59 * qJ(2);
t71 = t62 * pkin(7) + t73;
t70 = t61 * pkin(3) + t58 * qJ(4) + t78;
t69 = pkin(3) * t75 + t71;
t1 = (-m(4) * t78 - m(5) * t70 - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) - t87 * (t61 * pkin(4) + t70) + (-m(2) - m(3)) * pkin(6) + (t80 + t86) * t61 + (mrSges(4,2) + t82) * t58) * g(3) + (-m(3) * t52 - m(4) * t74 - mrSges(1,2) - t83 * (t62 * t61 * qJ(4) + t74) + (t82 * t61 + (m(5) * pkin(3) - m(6) * t77 - m(7) * (-pkin(8) + t77) + t84) * t58 + (m(3) + m(4) + t83) * qJ(2) - t81) * t62 + t79 * t59) * g(2) + (-m(3) * t73 - m(4) * t71 - m(5) * t69 - mrSges(1,1) - t87 * (pkin(4) * t75 + t69) + t79 * t62 + (t80 * t58 + (t83 * qJ(4) - t82) * t61 + t81) * t59) * g(1);
U  = t1;

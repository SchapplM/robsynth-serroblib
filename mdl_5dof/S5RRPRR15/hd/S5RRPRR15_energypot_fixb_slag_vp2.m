% Calculate potential energy for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR15_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:04
% EndTime: 2019-12-31 20:41:05
% DurationCPUTime: 0.42s
% Computational Cost: add. (116->65), mult. (179->62), div. (0->0), fcn. (156->8), ass. (0->29)
t95 = -mrSges(3,1) + mrSges(4,2) + m(6) * (-pkin(8) - pkin(7)) - mrSges(6,3);
t61 = qJ(4) + qJ(5);
t55 = sin(t61);
t56 = cos(t61);
t62 = sin(qJ(4));
t94 = -m(6) * pkin(4) * t62 - t55 * mrSges(6,1) - t56 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3);
t64 = sin(qJ(1));
t63 = sin(qJ(2));
t77 = qJ(3) * t63;
t66 = cos(qJ(2));
t82 = t64 * t66;
t93 = pkin(2) * t82 + t64 * t77;
t92 = m(4) + m(5) + m(6);
t91 = -m(5) * pkin(7) - mrSges(5,3);
t88 = t94 * t63 + t95 * t66 - mrSges(2,1);
t87 = t56 * mrSges(6,1) - t55 * mrSges(6,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t84 = t64 * t62;
t65 = cos(qJ(4));
t83 = t64 * t65;
t67 = cos(qJ(1));
t81 = t67 * t62;
t80 = t67 * t65;
t79 = t67 * t66;
t78 = t67 * pkin(1) + t64 * pkin(6);
t59 = t64 * pkin(1);
t76 = t59 + t93;
t75 = -t67 * pkin(6) + t59;
t54 = t65 * pkin(4) + pkin(3);
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(5) - t92 * (t63 * pkin(2) + pkin(5)) + (t62 * mrSges(5,1) + t65 * mrSges(5,2) + t92 * qJ(3) - t94) * t66 + (t91 + t95) * t63) * g(3) + (-mrSges(1,2) - m(3) * t75 - m(4) * (t75 + t93) - m(5) * (pkin(7) * t82 + t76) - (t63 * t84 - t80) * mrSges(5,1) - (t63 * t83 + t81) * mrSges(5,2) - mrSges(5,3) * t82 - m(6) * t76 + (-m(5) * (-pkin(3) - pkin(6)) - m(6) * (-pkin(6) - t54) + t87) * t67 + t88 * t64) * g(2) + (-mrSges(1,1) - m(3) * t78 - (t63 * t81 + t83) * mrSges(5,1) - (t63 * t80 - t84) * mrSges(5,2) + t91 * t79 - t92 * (pkin(2) * t79 + t67 * t77 + t78) + t88 * t67 + (-m(5) * pkin(3) - m(6) * t54 - t87) * t64) * g(1);
U = t1;

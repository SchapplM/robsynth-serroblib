% Calculate potential energy for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:39
% EndTime: 2019-12-05 15:59:39
% DurationCPUTime: 0.41s
% Computational Cost: add. (116->58), mult. (179->50), div. (0->0), fcn. (156->8), ass. (0->25)
t97 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3) + m(6) * (-pkin(7) - pkin(6)) - mrSges(6,3);
t62 = qJ(4) + qJ(5);
t56 = sin(t62);
t57 = cos(t62);
t65 = sin(qJ(4));
t67 = cos(qJ(4));
t96 = -t56 * mrSges(6,1) - t67 * mrSges(5,2) - t57 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t65;
t93 = m(5) * pkin(6);
t63 = sin(pkin(8));
t66 = sin(qJ(2));
t79 = qJ(3) * t66;
t68 = cos(qJ(2));
t85 = t63 * t68;
t92 = pkin(2) * t85 + t63 * t79;
t91 = m(4) + m(5) + m(6);
t88 = t67 * mrSges(5,1) + t57 * mrSges(6,1) - t65 * mrSges(5,2) - t56 * mrSges(6,2) + mrSges(4,1) - mrSges(2,2) + mrSges(3,3);
t87 = t96 * t66 + t97 * t68 - mrSges(2,1);
t64 = cos(pkin(8));
t84 = t64 * t68;
t80 = t64 * pkin(1) + t63 * pkin(5);
t59 = t63 * pkin(1);
t77 = t59 + t92;
t76 = -t64 * pkin(5) + t59;
t55 = t67 * pkin(4) + pkin(3);
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * qJ(1) - t91 * (t66 * pkin(2) + qJ(1)) + (t91 * qJ(3) - t96) * t68 + (-t93 + t97) * t66) * g(3) + (-mrSges(1,2) - m(3) * t76 - m(4) * (t76 + t92) - m(5) * (pkin(6) * t85 + t77) - m(6) * t77 + (-m(5) * (-pkin(3) - pkin(5)) - m(6) * (-pkin(5) - t55) + t88) * t64 + t87 * t63) * g(2) + (-t84 * t93 - m(3) * t80 - mrSges(1,1) - t91 * (pkin(2) * t84 + t64 * t79 + t80) + (-m(5) * pkin(3) - m(6) * t55 - t88) * t63 + t87 * t64) * g(1);
U = t1;

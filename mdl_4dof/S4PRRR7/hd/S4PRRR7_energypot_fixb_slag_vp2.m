% Calculate potential energy for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_energypot_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_energypot_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_energypot_fixb_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:00
% EndTime: 2019-12-31 16:36:00
% DurationCPUTime: 0.35s
% Computational Cost: add. (130->64), mult. (277->83), div. (0->0), fcn. (311->10), ass. (0->31)
t96 = -m(4) - m(5);
t95 = -m(5) * pkin(7) + mrSges(4,2) - mrSges(5,3);
t74 = sin(qJ(4));
t77 = cos(qJ(4));
t94 = -t74 * mrSges(5,1) - t77 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t93 = -m(5) * pkin(3) - t77 * mrSges(5,1) + t74 * mrSges(5,2) - mrSges(4,1);
t70 = sin(pkin(8));
t71 = sin(pkin(4));
t92 = t70 * t71;
t76 = sin(qJ(2));
t91 = t71 * t76;
t78 = cos(qJ(3));
t90 = t71 * t78;
t79 = cos(qJ(2));
t89 = t71 * t79;
t72 = cos(pkin(8));
t88 = t72 * t71;
t73 = cos(pkin(4));
t87 = t73 * t76;
t86 = t73 * t79;
t85 = t72 * pkin(1) + pkin(5) * t92;
t84 = t73 * pkin(5) + qJ(1);
t83 = t70 * pkin(1) - pkin(5) * t88;
t81 = pkin(2) * t91 - pkin(6) * t89 + t84;
t75 = sin(qJ(3));
t61 = t73 * t75 + t76 * t90;
t59 = -t70 * t87 + t72 * t79;
t58 = t70 * t86 + t72 * t76;
t57 = t70 * t79 + t72 * t87;
t56 = t70 * t76 - t72 * t86;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t84 - t73 * mrSges(3,3) - (t76 * mrSges(3,1) + t79 * mrSges(3,2)) * t71 - m(4) * t81 - t61 * mrSges(4,1) + mrSges(4,3) * t89 - m(5) * (t61 * pkin(3) + t81) - (t61 * t77 - t74 * t89) * mrSges(5,1) - (-t61 * t74 - t77 * t89) * mrSges(5,2) + t95 * (-t73 * t78 + t75 * t91)) * g(3) + (-m(3) * t83 - t70 * mrSges(2,1) - t57 * mrSges(3,1) - t72 * mrSges(2,2) + mrSges(3,3) * t88 - mrSges(1,2) + t96 * (t57 * pkin(2) + pkin(6) * t56 + t83) + t93 * (t57 * t78 - t75 * t88) + t94 * t56 + t95 * (t57 * t75 + t78 * t88)) * g(2) + (-m(3) * t85 - t72 * mrSges(2,1) - t59 * mrSges(3,1) + t70 * mrSges(2,2) - mrSges(3,3) * t92 - mrSges(1,1) + t96 * (t59 * pkin(2) + pkin(6) * t58 + t85) + t93 * (t59 * t78 + t75 * t92) + t94 * t58 + t95 * (t59 * t75 - t70 * t90)) * g(1);
U = t1;

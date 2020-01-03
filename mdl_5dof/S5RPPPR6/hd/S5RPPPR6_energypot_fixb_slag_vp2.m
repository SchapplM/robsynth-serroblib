% Calculate potential energy for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:33
% EndTime: 2019-12-31 17:47:33
% DurationCPUTime: 0.48s
% Computational Cost: add. (116->73), mult. (210->77), div. (0->0), fcn. (199->8), ass. (0->33)
t101 = -mrSges(3,1) + mrSges(4,2);
t100 = mrSges(3,2) - mrSges(4,3);
t99 = -m(5) - m(6);
t75 = sin(qJ(1));
t71 = sin(pkin(7));
t86 = qJ(3) * t71;
t73 = cos(pkin(7));
t93 = t73 * t75;
t98 = pkin(2) * t93 + t75 * t86;
t97 = t100 * t71 + t101 * t73 - mrSges(2,1);
t96 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t95 = m(6) * pkin(6) - mrSges(5,2) + mrSges(6,3);
t94 = t71 * pkin(2) + pkin(5);
t77 = cos(qJ(1));
t92 = t73 * t77;
t70 = sin(pkin(8));
t91 = t75 * t70;
t72 = cos(pkin(8));
t90 = t75 * t72;
t89 = t77 * t70;
t88 = t77 * t72;
t87 = t77 * pkin(1) + t75 * qJ(2);
t85 = qJ(4) * t73;
t68 = t75 * pkin(1);
t83 = -t77 * qJ(2) + t68;
t82 = pkin(2) * t92 + t77 * t86 + t87;
t79 = t75 * pkin(3) + t77 * t85 + t82;
t78 = t75 * t85 + t68 + (-pkin(3) - qJ(2)) * t77 + t98;
t76 = cos(qJ(5));
t74 = sin(qJ(5));
t56 = t71 * t91 - t88;
t54 = t71 * t89 + t90;
t1 = (-m(4) * t94 - mrSges(1,3) - mrSges(2,3) + t99 * (t71 * qJ(4) + t94) + (-m(2) - m(3)) * pkin(5) + ((m(6) * pkin(4) + t76 * mrSges(6,1) - t74 * mrSges(6,2) + mrSges(5,1)) * t70 - t95 * t72 + (m(4) - t99) * qJ(3) - t100) * t73 + (-t74 * mrSges(6,1) - t76 * mrSges(6,2) - mrSges(5,3) + t101) * t71) * g(3) + (-mrSges(1,2) - m(3) * t83 - m(4) * (t83 + t98) - m(5) * t78 - t56 * mrSges(5,1) - mrSges(5,3) * t93 - m(6) * (t56 * pkin(4) + t78) - (t56 * t76 + t74 * t93) * mrSges(6,1) - (-t56 * t74 + t76 * t93) * mrSges(6,2) + t95 * (t71 * t90 + t89) - t96 * t77 + t97 * t75) * g(2) + (-mrSges(1,1) - m(3) * t87 - m(4) * t82 - m(5) * t79 - t54 * mrSges(5,1) - mrSges(5,3) * t92 - m(6) * (t54 * pkin(4) + t79) - (t54 * t76 + t74 * t92) * mrSges(6,1) - (-t54 * t74 + t76 * t92) * mrSges(6,2) - t95 * (-t71 * t88 + t91) + t97 * t77 + t96 * t75) * g(1);
U = t1;

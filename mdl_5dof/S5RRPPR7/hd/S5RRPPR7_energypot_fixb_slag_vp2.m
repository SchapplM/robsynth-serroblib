% Calculate potential energy for
% S5RRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR7_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR7_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:34:51
% EndTime: 2019-12-31 19:34:51
% DurationCPUTime: 0.35s
% Computational Cost: add. (133->65), mult. (148->62), div. (0->0), fcn. (121->8), ass. (0->30)
t87 = -mrSges(4,1) + mrSges(5,2);
t86 = mrSges(4,2) - mrSges(5,3);
t85 = m(5) + m(6);
t65 = cos(qJ(1));
t58 = qJ(2) + pkin(8);
t55 = sin(t58);
t73 = qJ(4) * t55;
t56 = cos(t58);
t77 = t65 * t56;
t84 = pkin(3) * t77 + t65 * t73;
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t83 = -m(3) * pkin(1) - t64 * mrSges(3,1) + t61 * mrSges(3,2) + t86 * t55 + t87 * t56 - mrSges(2,1);
t82 = -m(3) * pkin(6) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t81 = t61 * pkin(2) + pkin(5);
t62 = sin(qJ(1));
t80 = t62 * t56;
t60 = sin(qJ(5));
t79 = t62 * t60;
t63 = cos(qJ(5));
t78 = t62 * t63;
t76 = t65 * t60;
t75 = t65 * t63;
t54 = t64 * pkin(2) + pkin(1);
t59 = -qJ(3) - pkin(6);
t74 = t62 * t54 + t65 * t59;
t51 = t65 * t54;
t70 = -t62 * t59 + t51;
t69 = pkin(3) * t80 + t62 * t73 + t74;
t1 = (-m(4) * t81 - t61 * mrSges(3,1) - t64 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) - t85 * (t55 * pkin(3) + t81) + (-m(2) - m(3)) * pkin(5) + (t60 * mrSges(6,1) + t63 * mrSges(6,2) + t85 * qJ(4) - t86) * t56 + (-m(6) * pkin(7) - mrSges(6,3) + t87) * t55) * g(3) + (-mrSges(1,2) - m(4) * t74 - m(5) * t69 - m(6) * (pkin(7) * t80 + t69) - (t55 * t79 - t75) * mrSges(6,1) - (t55 * t78 + t76) * mrSges(6,2) - mrSges(6,3) * t80 + (m(6) * pkin(4) - t82) * t65 + t83 * t62) * g(2) + (-mrSges(1,1) - m(4) * t70 - m(5) * (t70 + t84) - m(6) * (pkin(7) * t77 + t51 + t84) - (t55 * t76 + t78) * mrSges(6,1) - (t55 * t75 - t79) * mrSges(6,2) - mrSges(6,3) * t77 + t83 * t65 + (-m(6) * (pkin(4) - t59) + t82) * t62) * g(1);
U = t1;

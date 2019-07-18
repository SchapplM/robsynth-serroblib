% Calculate potential energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energypot_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:04
% EndTime: 2019-07-18 17:17:05
% DurationCPUTime: 0.33s
% Computational Cost: add. (126->65), mult. (149->72), div. (0->0), fcn. (128->10), ass. (0->31)
t84 = -mrSges(5,3) - mrSges(6,3);
t58 = qJ(2) + qJ(3);
t53 = sin(t58);
t55 = cos(t58);
t60 = sin(qJ(2));
t63 = cos(qJ(2));
t81 = pkin(1) * t63;
t83 = -m(4) * t81 - mrSges(3,1) * t63 - mrSges(4,1) * t55 + mrSges(3,2) * t60 + mrSges(4,2) * t53 - mrSges(2,1);
t82 = mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t80 = t60 * pkin(1) + pkin(4);
t61 = sin(qJ(1));
t79 = t53 * t61;
t64 = cos(qJ(1));
t78 = t53 * t64;
t77 = t55 * t61;
t76 = t55 * t64;
t59 = sin(qJ(4));
t75 = t59 * t64;
t57 = qJ(4) + qJ(5);
t52 = sin(t57);
t74 = t61 * t52;
t54 = cos(t57);
t73 = t61 * t54;
t72 = t61 * t59;
t62 = cos(qJ(4));
t71 = t61 * t62;
t70 = t62 * t64;
t69 = pkin(5) * t79 + t61 * t81;
t68 = pkin(5) * t78 + t64 * t81;
t51 = pkin(3) * t62 + pkin(2);
t1 = (-m(4) * t80 - t60 * mrSges(3,1) - t63 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + (-m(5) - m(6)) * (-pkin(5) * t55 + t80) + (-m(2) - m(3)) * pkin(4) + (-mrSges(4,2) - t84) * t55 + (-m(5) * pkin(2) - m(6) * t51 - mrSges(5,1) * t62 - mrSges(6,1) * t54 + mrSges(5,2) * t59 + mrSges(6,2) * t52 - mrSges(4,1)) * t53) * g(3) + (-mrSges(1,2) - m(5) * (pkin(2) * t77 + t69) + t75 * mrSges(5,1) + t70 * mrSges(5,2) - m(6) * (-pkin(3) * t75 + t51 * t77 + t69) + t84 * t79 + (-t71 * mrSges(5,1) - t73 * mrSges(6,1) + t72 * mrSges(5,2) + t74 * mrSges(6,2)) * t55 + (t52 * mrSges(6,1) + t54 * mrSges(6,2) - t82) * t64 + t83 * t61) * g(2) + (-mrSges(1,1) - m(5) * (pkin(2) * t76 + t68) - (t55 * t70 + t72) * mrSges(5,1) - (-t55 * t75 + t71) * mrSges(5,2) - m(6) * (pkin(3) * t72 + t51 * t76 + t68) - (t54 * t76 + t74) * mrSges(6,1) - (-t52 * t76 + t73) * mrSges(6,2) + t84 * t78 + t83 * t64 + t82 * t61) * g(1);
U  = t1;

% Calculate potential energy for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:22
% EndTime: 2019-03-09 13:12:22
% DurationCPUTime: 0.31s
% Computational Cost: add. (231->62), mult. (171->50), div. (0->0), fcn. (138->12), ass. (0->29)
t74 = sin(qJ(6));
t77 = cos(qJ(6));
t99 = -m(7) * pkin(5) - t77 * mrSges(7,1) + t74 * mrSges(7,2) - mrSges(6,1);
t98 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t97 = m(6) + m(7);
t78 = cos(qJ(2));
t64 = t78 * pkin(2) + pkin(1);
t72 = qJ(2) + pkin(11);
t66 = cos(t72);
t54 = pkin(3) * t66 + t64;
t67 = qJ(4) + t72;
t63 = qJ(5) + t67;
t57 = sin(t63);
t58 = cos(t63);
t61 = sin(t67);
t62 = cos(t67);
t65 = sin(t72);
t75 = sin(qJ(2));
t95 = -m(3) * pkin(1) - m(4) * t64 - m(5) * t54 - t78 * mrSges(3,1) - t66 * mrSges(4,1) - t62 * mrSges(5,1) + t75 * mrSges(3,2) + t65 * mrSges(4,2) + t61 * mrSges(5,2) + t98 * t57 + t99 * t58 - mrSges(2,1);
t73 = -qJ(3) - pkin(7);
t71 = -pkin(8) + t73;
t94 = m(3) * pkin(7) - m(4) * t73 - m(5) * t71 + t74 * mrSges(7,1) + t77 * mrSges(7,2) - mrSges(2,2) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t93 = t75 * pkin(2) + pkin(6);
t91 = pkin(3) * t65 + t93;
t79 = cos(qJ(1));
t76 = sin(qJ(1));
t68 = -pkin(9) + t71;
t53 = pkin(4) * t62 + t54;
t1 = (-m(4) * t93 - m(5) * t91 - t75 * mrSges(3,1) - t65 * mrSges(4,1) - t61 * mrSges(5,1) - t78 * mrSges(3,2) - t66 * mrSges(4,2) - t62 * mrSges(5,2) - mrSges(1,3) - mrSges(2,3) - t97 * (pkin(4) * t61 + t91) - t98 * t58 + t99 * t57 + (-m(2) - m(3)) * pkin(6)) * g(3) + (-mrSges(1,2) - t97 * (t76 * t53 + t79 * t68) + t94 * t79 + t95 * t76) * g(2) + (-mrSges(1,1) + (t97 * t68 - t94) * t76 + (-t97 * t53 + t95) * t79) * g(1);
U  = t1;

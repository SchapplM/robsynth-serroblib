% Calculate potential energy for
% S6RPRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPRP7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:21:07
% EndTime: 2019-03-09 03:21:08
% DurationCPUTime: 0.45s
% Computational Cost: add. (163->71), mult. (189->63), div. (0->0), fcn. (160->8), ass. (0->32)
t95 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t94 = m(7) * pkin(5);
t93 = m(3) + m(4);
t92 = -mrSges(6,1) - mrSges(7,1);
t91 = mrSges(6,2) + mrSges(7,2);
t90 = -m(5) - m(6) - m(7);
t62 = qJ(3) + pkin(9);
t56 = sin(t62);
t57 = cos(t62);
t66 = sin(qJ(3));
t69 = cos(qJ(3));
t89 = t66 * mrSges(4,1) + t56 * mrSges(5,1) + t69 * mrSges(4,2) + t95 * t57 - mrSges(2,2) + mrSges(3,3);
t88 = -m(4) * pkin(7) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t87 = pkin(2) + pkin(6);
t86 = pkin(3) * t66;
t65 = sin(qJ(5));
t67 = sin(qJ(1));
t83 = t67 * t65;
t68 = cos(qJ(5));
t82 = t67 * t68;
t70 = cos(qJ(1));
t81 = t70 * t65;
t80 = t70 * t68;
t79 = t70 * pkin(1) + t67 * qJ(2);
t77 = -qJ(2) - t86;
t75 = pkin(4) * t56 - pkin(8) * t57;
t55 = t68 * pkin(5) + pkin(4);
t63 = -qJ(6) - pkin(8);
t72 = t55 * t56 + t57 * t63;
t64 = -qJ(4) - pkin(7);
t59 = t67 * pkin(1);
t1 = (-m(4) * t87 - t69 * mrSges(4,1) + t66 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) + t90 * (t69 * pkin(3) + t87) + (-m(6) * pkin(4) - m(7) * t55 + t91 * t65 + t92 * t68 - mrSges(5,1)) * t57 + (-m(6) * pkin(8) + m(7) * t63 + t95) * t56) * g(3) + (-t83 * t94 - mrSges(1,2) + t90 * (-t67 * t64 + t59) - t93 * t59 + t92 * (-t56 * t80 + t83) - t91 * (t56 * t81 + t82) + t88 * t67 + (-m(5) * t77 - m(6) * (-t75 + t77) - m(7) * (-t72 + t77) + t93 * qJ(2) + t89) * t70) * g(2) + (-t81 * t94 - mrSges(1,1) - t93 * t79 + t90 * (-t70 * t64 + t67 * t86 + t79) + t92 * (t56 * t82 + t81) - t91 * (-t56 * t83 + t80) + t88 * t70 + (-m(6) * t75 - m(7) * t72 - t89) * t67) * g(1);
U  = t1;

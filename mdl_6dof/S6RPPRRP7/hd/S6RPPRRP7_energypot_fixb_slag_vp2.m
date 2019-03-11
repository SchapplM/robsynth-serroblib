% Calculate potential energy for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:41
% EndTime: 2019-03-09 02:12:41
% DurationCPUTime: 0.45s
% Computational Cost: add. (163->71), mult. (189->64), div. (0->0), fcn. (160->8), ass. (0->32)
t95 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t94 = m(7) * pkin(5);
t93 = m(3) + m(4);
t92 = -mrSges(6,1) - mrSges(7,1);
t91 = mrSges(6,2) + mrSges(7,2);
t90 = -m(5) - m(6) - m(7);
t62 = pkin(9) + qJ(4);
t56 = sin(t62);
t57 = cos(t62);
t63 = sin(pkin(9));
t64 = cos(pkin(9));
t89 = t63 * mrSges(4,1) + t56 * mrSges(5,1) + t64 * mrSges(4,2) + t95 * t57 - mrSges(2,2) + mrSges(3,3);
t88 = -m(4) * qJ(3) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t87 = pkin(2) + pkin(6);
t86 = pkin(3) * t63;
t67 = sin(qJ(5));
t70 = cos(qJ(1));
t83 = t67 * t70;
t68 = sin(qJ(1));
t82 = t68 * t67;
t69 = cos(qJ(5));
t81 = t68 * t69;
t80 = t70 * t56;
t79 = t70 * pkin(1) + t68 * qJ(2);
t77 = -qJ(2) - t86;
t75 = pkin(4) * t56 - pkin(8) * t57;
t55 = pkin(5) * t69 + pkin(4);
t65 = -qJ(6) - pkin(8);
t72 = t55 * t56 + t57 * t65;
t66 = -pkin(7) - qJ(3);
t60 = t68 * pkin(1);
t1 = (-m(4) * t87 - t64 * mrSges(4,1) + t63 * mrSges(4,2) - mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) + t90 * (t64 * pkin(3) + t87) + (-m(6) * pkin(4) - m(7) * t55 + t91 * t67 + t92 * t69 - mrSges(5,1)) * t57 + (-m(6) * pkin(8) + m(7) * t65 + t95) * t56) * g(3) + (-t82 * t94 - mrSges(1,2) + t90 * (-t68 * t66 + t60) - t93 * t60 + t92 * (-t69 * t80 + t82) - t91 * (t67 * t80 + t81) + t88 * t68 + (-m(5) * t77 - m(6) * (-t75 + t77) - m(7) * (-t72 + t77) + t93 * qJ(2) + t89) * t70) * g(2) + (-t83 * t94 - mrSges(1,1) - t93 * t79 + t90 * (-t66 * t70 + t68 * t86 + t79) + t92 * (t56 * t81 + t83) - t91 * (-t56 * t82 + t69 * t70) + t88 * t70 + (-m(6) * t75 - m(7) * t72 - t89) * t68) * g(1);
U  = t1;

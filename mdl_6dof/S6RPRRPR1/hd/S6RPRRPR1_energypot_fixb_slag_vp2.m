% Calculate potential energy for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:17
% EndTime: 2019-03-09 04:57:17
% DurationCPUTime: 0.33s
% Computational Cost: add. (223->63), mult. (157->51), div. (0->0), fcn. (124->12), ass. (0->31)
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t96 = -m(7) * pkin(5) - t75 * mrSges(7,1) + t72 * mrSges(7,2) - mrSges(6,1);
t95 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t94 = -m(3) - m(4);
t93 = m(6) + m(7);
t92 = -m(5) + t94;
t70 = qJ(3) + qJ(4);
t61 = pkin(11) + t70;
t55 = sin(t61);
t56 = cos(t61);
t76 = cos(qJ(3));
t58 = t76 * pkin(3) + pkin(2);
t62 = sin(t70);
t63 = cos(t70);
t73 = sin(qJ(3));
t90 = -m(4) * pkin(2) - m(5) * t58 - t76 * mrSges(4,1) - t63 * mrSges(5,1) + t73 * mrSges(4,2) + t62 * mrSges(5,2) + t95 * t55 + t96 * t56 - mrSges(3,1);
t78 = -pkin(8) - pkin(7);
t89 = m(4) * pkin(7) - m(5) * t78 + t72 * mrSges(7,1) + t75 * mrSges(7,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t74 = sin(qJ(1));
t65 = t74 * pkin(1);
t77 = cos(qJ(1));
t67 = t77 * pkin(1);
t71 = qJ(2) + pkin(6);
t86 = t73 * pkin(3) + t71;
t69 = qJ(1) + pkin(10);
t68 = -qJ(5) + t78;
t60 = cos(t69);
t59 = sin(t69);
t53 = pkin(4) * t63 + t58;
t1 = (-m(2) * pkin(6) - m(5) * t86 - mrSges(4,1) * t73 - t62 * mrSges(5,1) - mrSges(4,2) * t76 - t63 * mrSges(5,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t93 * (pkin(4) * t62 + t86) + t94 * t71 - t95 * t56 + t96 * t55) * g(3) + (-mrSges(2,1) * t74 - mrSges(2,2) * t77 - mrSges(1,2) - t93 * (t59 * t53 + t60 * t68 + t65) + t92 * t65 + t89 * t60 + t90 * t59) * g(2) + (-mrSges(2,1) * t77 + mrSges(2,2) * t74 - mrSges(1,1) - t93 * (t60 * t53 + t67) + t92 * t67 + t90 * t60 + (t93 * t68 - t89) * t59) * g(1);
U  = t1;

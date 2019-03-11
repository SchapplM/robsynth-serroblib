% Calculate potential energy for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:17:42
% EndTime: 2019-03-09 02:17:43
% DurationCPUTime: 0.33s
% Computational Cost: add. (223->63), mult. (157->51), div. (0->0), fcn. (124->12), ass. (0->31)
t75 = sin(qJ(6));
t77 = cos(qJ(6));
t96 = -m(7) * pkin(5) - t77 * mrSges(7,1) + t75 * mrSges(7,2) - mrSges(6,1);
t95 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t94 = -m(3) - m(4);
t93 = m(6) + m(7);
t92 = -m(5) + t94;
t69 = pkin(11) + qJ(4);
t63 = qJ(5) + t69;
t56 = sin(t63);
t57 = cos(t63);
t72 = cos(pkin(11));
t58 = t72 * pkin(3) + pkin(2);
t59 = sin(t69);
t61 = cos(t69);
t71 = sin(pkin(11));
t90 = -m(4) * pkin(2) - m(5) * t58 - t72 * mrSges(4,1) - t61 * mrSges(5,1) + t71 * mrSges(4,2) + t59 * mrSges(5,2) + t95 * t56 + t96 * t57 - mrSges(3,1);
t74 = -pkin(7) - qJ(3);
t89 = m(4) * qJ(3) - m(5) * t74 + t75 * mrSges(7,1) + t77 * mrSges(7,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t76 = sin(qJ(1));
t66 = t76 * pkin(1);
t78 = cos(qJ(1));
t67 = t78 * pkin(1);
t73 = qJ(2) + pkin(6);
t86 = t71 * pkin(3) + t73;
t70 = qJ(1) + pkin(10);
t68 = -pkin(8) + t74;
t62 = cos(t70);
t60 = sin(t70);
t53 = pkin(4) * t61 + t58;
t1 = (-m(2) * pkin(6) - m(5) * t86 - mrSges(4,1) * t71 - t59 * mrSges(5,1) - mrSges(4,2) * t72 - t61 * mrSges(5,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t93 * (pkin(4) * t59 + t86) + t94 * t73 - t95 * t57 + t96 * t56) * g(3) + (-t76 * mrSges(2,1) - mrSges(2,2) * t78 - mrSges(1,2) - t93 * (t60 * t53 + t62 * t68 + t66) + t92 * t66 + t89 * t62 + t90 * t60) * g(2) + (-mrSges(2,1) * t78 + t76 * mrSges(2,2) - mrSges(1,1) - t93 * (t62 * t53 + t67) + t92 * t67 + t90 * t62 + (t93 * t68 - t89) * t60) * g(1);
U  = t1;

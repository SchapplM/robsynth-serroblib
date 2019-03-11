% Calculate potential energy for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_energypot_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:42
% EndTime: 2019-03-09 01:50:42
% DurationCPUTime: 0.31s
% Computational Cost: add. (110->64), mult. (159->46), div. (0->0), fcn. (126->6), ass. (0->26)
t77 = m(6) + m(7);
t76 = -m(7) * pkin(8) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t51 = sin(qJ(6));
t54 = cos(qJ(6));
t75 = -t51 * mrSges(7,1) - t54 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t73 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3) + (t77 * qJ(5) - t75) * t55 + t76 * t52;
t72 = -t54 * mrSges(7,1) + t51 * mrSges(7,2) - mrSges(6,1) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3);
t71 = pkin(2) + pkin(6);
t70 = pkin(4) * t52;
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t69 = t56 * pkin(1) + t53 * qJ(2);
t68 = pkin(3) + t71;
t67 = t56 * qJ(3) + t69;
t47 = t53 * pkin(1);
t66 = -t56 * qJ(2) + t47;
t44 = t53 * qJ(3);
t64 = t44 + t66;
t60 = -t53 * pkin(7) + t67;
t49 = t56 * pkin(7);
t59 = t49 + t64;
t42 = t56 * t70;
t41 = t53 * t70;
t1 = (-m(4) * t71 - m(5) * t68 - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - t77 * (t55 * pkin(4) + t52 * qJ(5) + t68) + (-m(2) - m(3)) * pkin(6) + t76 * t55 + t75 * t52) * g(3) + (-mrSges(1,2) - m(3) * t66 - m(4) * t64 - m(5) * t59 - m(6) * (t41 + t59) - m(7) * (t41 + t44 + t47 + t49) + (-m(7) * (pkin(5) - qJ(2)) + t72) * t56 + t73 * t53) * g(2) + (-mrSges(1,1) - m(3) * t69 - m(4) * t67 - m(5) * t60 - m(6) * (t42 + t60) - m(7) * (t42 + t67) + t73 * t56 + (-m(7) * (-pkin(5) - pkin(7)) - t72) * t53) * g(1);
U  = t1;

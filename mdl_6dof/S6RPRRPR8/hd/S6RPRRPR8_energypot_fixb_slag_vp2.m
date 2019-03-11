% Calculate potential energy for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:40
% EndTime: 2019-03-09 05:22:40
% DurationCPUTime: 0.39s
% Computational Cost: add. (174->59), mult. (209->47), div. (0->0), fcn. (184->10), ass. (0->23)
t62 = qJ(4) + pkin(10);
t54 = qJ(6) + t62;
t49 = sin(t54);
t50 = cos(t54);
t67 = cos(qJ(4));
t51 = t67 * pkin(4) + pkin(3);
t52 = sin(t62);
t53 = cos(t62);
t64 = sin(qJ(4));
t89 = mrSges(4,1) + m(5) * pkin(3) + t67 * mrSges(5,1) - t64 * mrSges(5,2) + m(6) * t51 + t53 * mrSges(6,1) - t52 * mrSges(6,2) + m(7) * (pkin(5) * t53 + t51) + t50 * mrSges(7,1) - t49 * mrSges(7,2);
t63 = -qJ(5) - pkin(8);
t88 = mrSges(4,2) - m(5) * pkin(8) + m(6) * t63 + m(7) * (-pkin(9) + t63) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t65 = sin(qJ(3));
t68 = cos(qJ(3));
t87 = t89 * t65 + t88 * t68 - mrSges(2,2) + mrSges(3,3);
t84 = -m(4) - m(5) - m(6) - m(7);
t80 = t64 * pkin(4);
t83 = -m(6) * t80 - t52 * mrSges(6,1) - t53 * mrSges(6,2) - mrSges(2,1) + mrSges(3,2) - mrSges(4,3) - m(7) * (pkin(5) * t52 + t80) - t49 * mrSges(7,1) - t50 * mrSges(7,2) - t64 * mrSges(5,1) - t67 * mrSges(5,2);
t66 = sin(qJ(1));
t69 = cos(qJ(1));
t78 = t69 * pkin(1) + t66 * qJ(2);
t57 = t66 * pkin(1);
t1 = (-mrSges(3,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) + t84 * (pkin(2) + pkin(6)) - t89 * t68 + t88 * t65) * g(3) + (-m(3) * t57 - mrSges(1,2) + t84 * (t66 * pkin(7) + t57) + ((m(3) - t84) * qJ(2) + t87) * t69 + t83 * t66) * g(2) + (-m(3) * t78 - mrSges(1,1) + t84 * (t69 * pkin(7) + t78) + t83 * t69 - t87 * t66) * g(1);
U  = t1;

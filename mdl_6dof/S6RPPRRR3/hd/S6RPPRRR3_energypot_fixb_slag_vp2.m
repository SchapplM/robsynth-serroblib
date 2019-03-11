% Calculate potential energy for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:44
% EndTime: 2019-03-09 02:22:44
% DurationCPUTime: 0.38s
% Computational Cost: add. (193->59), mult. (170->46), div. (0->0), fcn. (141->10), ass. (0->23)
t58 = qJ(5) + qJ(6);
t53 = sin(t58);
t54 = cos(t58);
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t84 = mrSges(5,1) + m(6) * pkin(4) + t63 * mrSges(6,1) - t60 * mrSges(6,2) + m(7) * (t63 * pkin(5) + pkin(4)) + t54 * mrSges(7,1) - t53 * mrSges(7,2);
t83 = mrSges(5,2) - m(6) * pkin(8) + m(7) * (-pkin(9) - pkin(8)) - mrSges(6,3) - mrSges(7,3);
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t82 = t84 * t61 + t83 * t64 - mrSges(3,2) + mrSges(4,3);
t80 = -m(5) - m(6) - m(7);
t78 = -t53 * mrSges(7,1) - t63 * mrSges(6,2) - t54 * mrSges(7,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t60;
t62 = sin(qJ(1));
t55 = t62 * pkin(1);
t65 = cos(qJ(1));
t56 = t65 * pkin(1);
t59 = qJ(2) + pkin(6);
t57 = qJ(1) + pkin(10);
t51 = sin(t57);
t76 = t51 * pkin(2) + t55;
t52 = cos(t57);
t73 = t52 * pkin(2) + t51 * qJ(3) + t56;
t1 = (-m(2) * pkin(6) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - m(4)) * t59 + t80 * (pkin(3) + t59) - t84 * t64 + t83 * t61) * g(3) + (-m(3) * t55 - m(4) * t76 - t62 * mrSges(2,1) - t65 * mrSges(2,2) - mrSges(1,2) + t80 * (t51 * pkin(7) + t76) + ((m(4) - t80) * qJ(3) + t82) * t52 + t78 * t51) * g(2) + (-m(3) * t56 - m(4) * t73 - t65 * mrSges(2,1) + t62 * mrSges(2,2) - mrSges(1,1) + t80 * (t52 * pkin(7) + t73) + t78 * t52 - t82 * t51) * g(1);
U  = t1;

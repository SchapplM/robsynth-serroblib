% Calculate potential energy for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:37
% EndTime: 2019-03-09 01:31:37
% DurationCPUTime: 0.38s
% Computational Cost: add. (189->64), mult. (150->50), div. (0->0), fcn. (117->10), ass. (0->26)
t82 = m(6) + m(7);
t61 = sin(qJ(6));
t63 = cos(qJ(6));
t81 = -m(7) * pkin(5) - t63 * mrSges(7,1) + t61 * mrSges(7,2) - mrSges(6,1);
t80 = m(4) + m(5);
t78 = -m(7) * pkin(8) - mrSges(7,3);
t55 = pkin(10) + qJ(5);
t48 = sin(t55);
t50 = cos(t55);
t57 = sin(pkin(10));
t58 = cos(pkin(10));
t77 = -t57 * mrSges(5,1) - t58 * mrSges(5,2) - t50 * mrSges(6,2) + t81 * t48 + mrSges(3,2) - mrSges(4,3);
t76 = -m(5) * qJ(4) - t61 * mrSges(7,1) - t63 * mrSges(7,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) + t82 * (-pkin(7) - qJ(4));
t75 = pkin(4) * t57;
t62 = sin(qJ(1));
t53 = t62 * pkin(1);
t64 = cos(qJ(1));
t54 = t64 * pkin(1);
t59 = qJ(2) + pkin(6);
t73 = pkin(3) + t59;
t56 = qJ(1) + pkin(9);
t49 = sin(t56);
t51 = cos(t56);
t72 = t51 * pkin(2) + t49 * qJ(3) + t54;
t71 = -qJ(3) - t75;
t1 = (-m(2) * pkin(6) - m(5) * t73 - t58 * mrSges(5,1) + t57 * mrSges(5,2) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t82 * (t58 * pkin(4) + t73) + (-m(3) - m(4)) * t59 + t81 * t50 + (mrSges(6,2) + t78) * t48) * g(3) + (-m(3) * t53 - t62 * mrSges(2,1) - t64 * mrSges(2,2) - mrSges(1,2) + (-t82 - t80) * (t49 * pkin(2) + t53) + (-m(6) * t71 - m(7) * (pkin(8) * t50 + t71) - t50 * mrSges(7,3) + t80 * qJ(3) - t77) * t51 + t76 * t49) * g(2) + (-m(3) * t54 - t64 * mrSges(2,1) + t62 * mrSges(2,2) - mrSges(1,1) - t80 * t72 - t82 * (t49 * t75 + t72) + t76 * t51 + (-t78 * t50 + t77) * t49) * g(1);
U  = t1;

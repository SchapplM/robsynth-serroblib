% Calculate potential energy for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:11
% EndTime: 2019-03-09 02:28:11
% DurationCPUTime: 0.30s
% Computational Cost: add. (130->62), mult. (152->46), div. (0->0), fcn. (119->8), ass. (0->25)
t57 = sin(qJ(6));
t60 = cos(qJ(6));
t85 = -m(7) * pkin(5) - t60 * mrSges(7,1) + t57 * mrSges(7,2) - mrSges(6,1);
t84 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t83 = -m(6) - m(7);
t56 = qJ(4) + qJ(5);
t48 = sin(t56);
t49 = cos(t56);
t58 = sin(qJ(4));
t61 = cos(qJ(4));
t81 = -t58 * mrSges(5,1) - t61 * mrSges(5,2) + t85 * t48 - t84 * t49 - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t80 = -t57 * mrSges(7,1) - t60 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3) - mrSges(6,3);
t79 = pkin(2) + pkin(6);
t78 = pkin(4) * t58;
t59 = sin(qJ(1));
t50 = t59 * qJ(3);
t53 = t59 * pkin(1);
t77 = t50 + t53;
t62 = cos(qJ(1));
t76 = t62 * pkin(1) + t59 * qJ(2);
t74 = pkin(3) + t79;
t72 = t62 * qJ(3) + t76;
t70 = -t62 * qJ(2) + t53;
t63 = -pkin(8) - pkin(7);
t1 = (-m(4) * t79 - m(5) * t74 - t61 * mrSges(5,1) + t58 * mrSges(5,2) - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) + t83 * (t61 * pkin(4) + t74) + t85 * t49 + t84 * t48 + (-m(2) - m(3)) * pkin(6)) * g(3) + (-mrSges(1,2) - m(3) * t70 - m(4) * (t50 + t70) - m(5) * t77 + t83 * (t59 * t78 + t77) + (-m(5) * (pkin(7) - qJ(2)) + t83 * (-qJ(2) - t63) + t80) * t62 + t81 * t59) * g(2) + (-m(3) * t76 - mrSges(1,1) + (-m(4) - m(5)) * t72 + t83 * (t59 * t63 + t62 * t78 + t72) + t81 * t62 + (m(5) * pkin(7) - t80) * t59) * g(1);
U  = t1;

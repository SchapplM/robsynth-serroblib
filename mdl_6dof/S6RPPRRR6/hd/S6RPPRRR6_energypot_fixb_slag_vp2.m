% Calculate potential energy for
% S6RPPRRR6
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:32
% EndTime: 2019-03-09 02:30:32
% DurationCPUTime: 0.32s
% Computational Cost: add. (128->62), mult. (172->45), div. (0->0), fcn. (143->8), ass. (0->24)
t53 = qJ(5) + qJ(6);
t45 = sin(t53);
t46 = cos(t53);
t54 = sin(qJ(5));
t57 = cos(qJ(5));
t83 = -mrSges(5,1) - m(6) * pkin(4) - t57 * mrSges(6,1) + t54 * mrSges(6,2) - m(7) * (pkin(5) * t57 + pkin(4)) - t46 * mrSges(7,1) + t45 * mrSges(7,2);
t82 = mrSges(5,2) + m(7) * (-pkin(9) - pkin(8)) - mrSges(7,3) - m(6) * pkin(8) - mrSges(6,3);
t81 = -m(6) - m(7);
t80 = m(5) - t81;
t55 = sin(qJ(4));
t58 = cos(qJ(4));
t77 = t83 * t55 - t82 * t58 - mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t76 = -t54 * mrSges(6,1) - t45 * mrSges(7,1) - t57 * mrSges(6,2) - t46 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3);
t75 = pkin(2) + pkin(6);
t74 = pkin(5) * t54;
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t73 = t59 * pkin(1) + t56 * qJ(2);
t50 = t56 * pkin(1);
t69 = -t59 * qJ(2) + t50;
t47 = t56 * qJ(3);
t68 = t47 + t69;
t51 = t59 * pkin(7);
t1 = (-m(4) * t75 - mrSges(3,1) - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t80 * (pkin(3) + t75) + t83 * t58 + t82 * t55) * g(3) + (-mrSges(1,2) - m(3) * t69 - m(4) * t68 - m(5) * (t51 + t68) + t81 * (t47 + t50 + t51) + (m(6) * qJ(2) - m(7) * (-qJ(2) + t74) + t76) * t59 + t77 * t56) * g(2) + (-m(3) * t73 - mrSges(1,1) + (-m(4) - t80) * (t59 * qJ(3) + t73) + t77 * t59 + (m(7) * t74 + t80 * pkin(7) - t76) * t56) * g(1);
U  = t1;

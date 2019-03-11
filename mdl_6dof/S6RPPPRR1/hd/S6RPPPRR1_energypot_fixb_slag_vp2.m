% Calculate potential energy for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:47
% EndTime: 2019-03-09 01:29:48
% DurationCPUTime: 0.31s
% Computational Cost: add. (159->60), mult. (136->44), div. (0->0), fcn. (103->8), ass. (0->25)
t52 = sin(qJ(6));
t55 = cos(qJ(6));
t74 = -m(7) * pkin(5) - mrSges(7,1) * t55 + mrSges(7,2) * t52 - mrSges(6,1);
t73 = -m(7) * pkin(8) + mrSges(6,2) - mrSges(7,3);
t72 = m(6) + m(7);
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t70 = t74 * t53 - t73 * t56 - mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t69 = t52 * mrSges(7,1) + t55 * mrSges(7,2) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t54 = sin(qJ(1));
t48 = t54 * pkin(1);
t57 = cos(qJ(1));
t49 = t57 * pkin(1);
t51 = qJ(2) + pkin(6);
t50 = qJ(1) + pkin(9);
t46 = sin(t50);
t68 = t46 * pkin(2) + t48;
t67 = pkin(3) + t51;
t47 = cos(t50);
t66 = t47 * pkin(2) + t46 * qJ(3) + t49;
t63 = -qJ(3) * t47 + t68;
t40 = t46 * qJ(4);
t60 = t40 + t63;
t44 = t47 * pkin(7);
t1 = (-m(2) * pkin(6) - m(5) * t67 - mrSges(4,1) - mrSges(5,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) - t72 * (pkin(4) + t67) + t74 * t56 + t73 * t53 + (-m(3) - m(4)) * t51) * g(3) + (-mrSges(1,2) - t54 * mrSges(2,1) - t57 * mrSges(2,2) - m(3) * t48 - m(4) * t63 - m(5) * t60 - m(6) * (t44 + t60) - m(7) * (t40 + t44 + t68) + (m(7) * qJ(3) - t69) * t47 + t70 * t46) * g(2) + (-m(3) * t49 - m(4) * t66 - t57 * mrSges(2,1) + t54 * mrSges(2,2) - mrSges(1,1) + (-m(5) - t72) * (t47 * qJ(4) + t66) + t70 * t47 + (t72 * pkin(7) + t69) * t46) * g(1);
U  = t1;

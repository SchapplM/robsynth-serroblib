% Calculate potential energy for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:56:56
% EndTime: 2019-03-09 06:56:56
% DurationCPUTime: 0.37s
% Computational Cost: add. (234->61), mult. (181->50), div. (0->0), fcn. (152->12), ass. (0->27)
t66 = qJ(5) + qJ(6);
t58 = sin(t66);
t60 = cos(t66);
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t96 = -mrSges(5,1) - m(6) * pkin(4) - t72 * mrSges(6,1) + t69 * mrSges(6,2) - m(7) * (pkin(5) * t72 + pkin(4)) - t60 * mrSges(7,1) + t58 * mrSges(7,2);
t95 = mrSges(5,2) + m(7) * (-pkin(10) - pkin(9)) - mrSges(7,3) - m(6) * pkin(9) - mrSges(6,3);
t94 = -m(3) - m(4);
t93 = m(5) + m(6) + m(7);
t67 = qJ(3) + qJ(4);
t59 = sin(t67);
t61 = cos(t67);
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t90 = -m(4) * pkin(2) - t73 * mrSges(4,1) + t70 * mrSges(4,2) + t95 * t59 + t96 * t61 - mrSges(3,1);
t89 = m(4) * pkin(7) + t58 * mrSges(7,1) + t72 * mrSges(6,2) + t60 * mrSges(7,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t69;
t71 = sin(qJ(1));
t63 = t71 * pkin(1);
t74 = cos(qJ(1));
t64 = t74 * pkin(1);
t68 = qJ(2) + pkin(6);
t76 = -pkin(8) - pkin(7);
t65 = qJ(1) + pkin(11);
t57 = cos(t65);
t56 = sin(t65);
t55 = pkin(3) * t73 + pkin(2);
t1 = (-m(2) * pkin(6) - mrSges(4,1) * t70 - mrSges(4,2) * t73 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t94 * t68 - t93 * (t70 * pkin(3) + t68) - t95 * t61 + t96 * t59) * g(3) + (-mrSges(2,1) * t71 - mrSges(2,2) * t74 - mrSges(1,2) + t94 * t63 - t93 * (t56 * t55 + t57 * t76 + t63) + t89 * t57 + t90 * t56) * g(2) + (-mrSges(2,1) * t74 + t71 * mrSges(2,2) - mrSges(1,1) + t94 * t64 - t93 * (t57 * t55 + t64) + t90 * t57 + (t93 * t76 - t89) * t56) * g(1);
U  = t1;

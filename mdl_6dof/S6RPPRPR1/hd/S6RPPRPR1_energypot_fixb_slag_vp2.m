% Calculate potential energy for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:38:56
% EndTime: 2019-03-09 01:38:56
% DurationCPUTime: 0.36s
% Computational Cost: add. (234->61), mult. (181->50), div. (0->0), fcn. (152->12), ass. (0->27)
t65 = pkin(11) + qJ(6);
t56 = sin(t65);
t59 = cos(t65);
t68 = sin(pkin(11));
t70 = cos(pkin(11));
t96 = -mrSges(5,1) - m(6) * pkin(4) - t70 * mrSges(6,1) + t68 * mrSges(6,2) - m(7) * (t70 * pkin(5) + pkin(4)) - t59 * mrSges(7,1) + t56 * mrSges(7,2);
t95 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(8) - qJ(5)) - mrSges(7,3);
t94 = -m(3) - m(4);
t93 = m(5) + m(6) + m(7);
t66 = pkin(10) + qJ(4);
t57 = sin(t66);
t60 = cos(t66);
t69 = sin(pkin(10));
t71 = cos(pkin(10));
t90 = -m(4) * pkin(2) - t71 * mrSges(4,1) + t69 * mrSges(4,2) + t95 * t57 + t96 * t60 - mrSges(3,1);
t89 = m(4) * qJ(3) + t56 * mrSges(7,1) + t70 * mrSges(6,2) + t59 * mrSges(7,2) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + (m(7) * pkin(5) + mrSges(6,1)) * t68;
t75 = sin(qJ(1));
t63 = t75 * pkin(1);
t76 = cos(qJ(1));
t64 = t76 * pkin(1);
t72 = qJ(2) + pkin(6);
t74 = -pkin(7) - qJ(3);
t67 = qJ(1) + pkin(9);
t61 = cos(t67);
t58 = sin(t67);
t55 = t71 * pkin(3) + pkin(2);
t1 = (-m(2) * pkin(6) - t69 * mrSges(4,1) - t71 * mrSges(4,2) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t94 * t72 - t93 * (t69 * pkin(3) + t72) - t95 * t60 + t96 * t57) * g(3) + (-t75 * mrSges(2,1) - t76 * mrSges(2,2) - mrSges(1,2) + t94 * t63 - t93 * (t58 * t55 + t61 * t74 + t63) + t89 * t61 + t90 * t58) * g(2) + (-t76 * mrSges(2,1) + t75 * mrSges(2,2) - mrSges(1,1) + t94 * t64 - t93 * (t61 * t55 + t64) + t90 * t61 + (t93 * t74 - t89) * t58) * g(1);
U  = t1;

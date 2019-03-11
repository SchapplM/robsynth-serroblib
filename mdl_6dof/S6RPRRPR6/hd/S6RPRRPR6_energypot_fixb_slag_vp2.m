% Calculate potential energy for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:25
% EndTime: 2019-03-09 05:15:25
% DurationCPUTime: 0.35s
% Computational Cost: add. (236->61), mult. (221->50), div. (0->0), fcn. (196->12), ass. (0->26)
t72 = qJ(4) + pkin(11);
t67 = qJ(6) + t72;
t58 = sin(t67);
t59 = cos(t67);
t79 = cos(qJ(4));
t62 = t79 * pkin(4) + pkin(3);
t64 = sin(t72);
t66 = cos(t72);
t77 = sin(qJ(4));
t103 = -mrSges(4,1) - m(5) * pkin(3) - mrSges(5,1) * t79 + mrSges(5,2) * t77 - m(6) * t62 - mrSges(6,1) * t66 + mrSges(6,2) * t64 - m(7) * (pkin(5) * t66 + t62) - mrSges(7,1) * t59 + mrSges(7,2) * t58;
t75 = -qJ(5) - pkin(8);
t102 = mrSges(4,2) - m(5) * pkin(8) - mrSges(5,3) + m(7) * (-pkin(9) + t75) - mrSges(7,3) + m(6) * t75 - mrSges(6,3);
t101 = m(4) + m(5) + m(6) + m(7);
t71 = pkin(10) + qJ(3);
t63 = sin(t71);
t65 = cos(t71);
t73 = sin(pkin(10));
t74 = cos(pkin(10));
t97 = -m(3) * pkin(1) - t74 * mrSges(3,1) + t73 * mrSges(3,2) + t102 * t63 + t103 * t65 - mrSges(2,1);
t95 = t77 * pkin(4);
t96 = m(6) * t95 + m(7) * (pkin(5) * t64 + t95) - mrSges(2,2) + mrSges(4,3) + t58 * mrSges(7,1) + t59 * mrSges(7,2) + t64 * mrSges(6,1) + t66 * mrSges(6,2) + t77 * mrSges(5,1) + t79 * mrSges(5,2) + m(3) * qJ(2) + mrSges(3,3);
t80 = cos(qJ(1));
t78 = sin(qJ(1));
t76 = -pkin(7) - qJ(2);
t60 = t74 * pkin(2) + pkin(1);
t1 = (-t73 * mrSges(3,1) - t74 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t101 * (t73 * pkin(2) + pkin(6)) - t102 * t65 + t103 * t63) * g(3) + (-mrSges(1,2) - t101 * (t78 * t60 + t80 * t76) + t96 * t80 + t97 * t78) * g(2) + (-mrSges(1,1) + (t101 * t76 - t96) * t78 + (-t101 * t60 + t97) * t80) * g(1);
U  = t1;

% Calculate potential energy for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:58:42
% EndTime: 2019-03-09 05:58:42
% DurationCPUTime: 0.43s
% Computational Cost: add. (235->67), mult. (207->61), div. (0->0), fcn. (182->10), ass. (0->29)
t78 = cos(qJ(4));
t63 = t78 * pkin(4) + pkin(3);
t73 = qJ(4) + qJ(5);
t67 = cos(t73);
t75 = sin(qJ(4));
t102 = -m(6) * t63 - m(7) * (pkin(5) * t67 + t63) - mrSges(4,1) - m(5) * pkin(3) - t78 * mrSges(5,1) + t75 * mrSges(5,2);
t81 = -pkin(9) - pkin(8);
t101 = m(6) * t81 + m(7) * (-qJ(6) + t81) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(8) - mrSges(5,3);
t100 = -mrSges(6,1) - mrSges(7,1);
t99 = mrSges(6,2) + mrSges(7,2);
t98 = -m(4) - m(6) - m(7);
t97 = -m(5) + t98;
t76 = sin(qJ(3));
t79 = cos(qJ(3));
t95 = t101 * t76 + t102 * t79 - mrSges(3,1);
t66 = sin(t73);
t93 = t75 * pkin(4);
t94 = -m(6) * t93 - m(7) * (pkin(5) * t66 + t93) + mrSges(3,2) - mrSges(4,3) - t75 * mrSges(5,1) - t78 * mrSges(5,2);
t77 = sin(qJ(1));
t68 = t77 * pkin(1);
t80 = cos(qJ(1));
t70 = t80 * pkin(1);
t92 = t66 * t79;
t91 = t67 * t79;
t72 = qJ(1) + pkin(10);
t64 = sin(t72);
t90 = t64 * pkin(2) + t68;
t65 = cos(t72);
t1 = (-m(2) * pkin(6) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) + t97) * (qJ(2) + pkin(6)) - t101 * t79 + (t100 * t67 + t99 * t66 + t102) * t76) * g(3) + (-m(3) * t68 - m(5) * t90 - t77 * mrSges(2,1) - t80 * mrSges(2,2) - mrSges(1,2) + t100 * (t64 * t91 - t65 * t66) - t99 * (-t64 * t92 - t65 * t67) + t98 * (-t65 * pkin(7) + t90) + (m(5) * pkin(7) - t94) * t65 + t95 * t64) * g(2) + (-m(3) * t70 - t80 * mrSges(2,1) + t77 * mrSges(2,2) - mrSges(1,1) + t100 * (t64 * t66 + t65 * t91) - t99 * (t64 * t67 - t65 * t92) + t97 * (t65 * pkin(2) + t64 * pkin(7) + t70) + t94 * t64 + t95 * t65) * g(1);
U  = t1;

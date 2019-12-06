% Calculate potential energy for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRR3_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:18
% EndTime: 2019-12-05 15:16:18
% DurationCPUTime: 0.60s
% Computational Cost: add. (150->58), mult. (261->58), div. (0->0), fcn. (266->10), ass. (0->27)
t80 = sin(qJ(4));
t82 = cos(qJ(4));
t114 = -(m(6) * pkin(4) + mrSges(5,1)) * t80 - t82 * mrSges(5,2) + mrSges(3,2);
t108 = m(4) + m(5) + m(6);
t76 = sin(pkin(9));
t78 = cos(pkin(9));
t110 = -mrSges(2,1) + (-t108 * pkin(2) - mrSges(3,1)) * t78 + t114 * t76;
t109 = mrSges(2,2) - mrSges(3,3);
t75 = qJ(4) + qJ(5);
t70 = sin(t75);
t71 = cos(t75);
t104 = t70 * mrSges(6,1) + t71 * mrSges(6,2) + mrSges(4,3);
t103 = -m(5) * pkin(6) + m(6) * (-pkin(7) - pkin(6)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t102 = -m(5) * pkin(3) - m(6) * (t82 * pkin(4) + pkin(3)) - t82 * mrSges(5,1) - t71 * mrSges(6,1) + t80 * mrSges(5,2) + t70 * mrSges(6,2) - mrSges(4,1);
t77 = sin(pkin(8));
t100 = t76 * t77;
t79 = cos(pkin(8));
t99 = t76 * t79;
t81 = sin(qJ(3));
t96 = t77 * t81;
t83 = cos(qJ(3));
t95 = t77 * t83;
t94 = t79 * t81;
t93 = t79 * t83;
t92 = t79 * pkin(1) + t77 * qJ(2);
t88 = t77 * pkin(1) - t79 * qJ(2);
t1 = (-mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * qJ(1) - t108 * (t76 * pkin(2) + qJ(1)) + (pkin(5) * t108 + t104 - t114) * t78 + (t102 * t83 + t103 * t81 - mrSges(3,1)) * t76) * g(3) + (-m(3) * t88 - mrSges(1,2) - t108 * (pkin(5) * t100 + t88) - t109 * t79 - t104 * t100 + t102 * (t78 * t95 - t94) + t103 * (t78 * t96 + t93) + t110 * t77) * g(2) + (-m(3) * t92 - mrSges(1,1) - t108 * (pkin(5) * t99 + t92) - t104 * t99 + t102 * (t78 * t93 + t96) + t109 * t77 + t103 * (t78 * t94 - t95) + t110 * t79) * g(1);
U = t1;

% Calculate potential energy for
% S6RPRRRP6
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:09
% EndTime: 2019-03-09 06:14:09
% DurationCPUTime: 0.40s
% Computational Cost: add. (226->67), mult. (221->59), div. (0->0), fcn. (196->10), ass. (0->31)
t82 = cos(qJ(4));
t67 = t82 * pkin(4) + pkin(3);
t76 = qJ(4) + qJ(5);
t71 = cos(t76);
t80 = sin(qJ(4));
t109 = -m(6) * t67 - m(7) * (pkin(5) * t71 + t67) - mrSges(4,1) - m(5) * pkin(3) - t82 * mrSges(5,1) + t80 * mrSges(5,2);
t84 = -pkin(9) - pkin(8);
t108 = m(6) * t84 + m(7) * (-qJ(6) + t84) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(8) - mrSges(5,3);
t107 = -mrSges(6,1) - mrSges(7,1);
t106 = mrSges(6,2) + mrSges(7,2);
t105 = -m(4) - m(6) - m(7);
t104 = -m(5) + t105;
t75 = pkin(10) + qJ(3);
t68 = sin(t75);
t69 = cos(t75);
t77 = sin(pkin(10));
t78 = cos(pkin(10));
t102 = -m(3) * pkin(1) - t78 * mrSges(3,1) + t77 * mrSges(3,2) + t108 * t68 + t109 * t69 - mrSges(2,1);
t100 = pkin(4) * t80;
t70 = sin(t76);
t101 = m(6) * t100 + m(7) * (pkin(5) * t70 + t100) - mrSges(2,2) + mrSges(4,3) + t80 * mrSges(5,1) + t82 * mrSges(5,2) + m(3) * qJ(2) + mrSges(3,3);
t81 = sin(qJ(1));
t98 = t70 * t81;
t83 = cos(qJ(1));
t97 = t70 * t83;
t96 = t71 * t81;
t95 = t71 * t83;
t79 = -pkin(7) - qJ(2);
t65 = pkin(2) * t78 + pkin(1);
t64 = t83 * t65;
t1 = (-mrSges(3,1) * t77 - mrSges(3,2) * t78 - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) + t104 * (t77 * pkin(2) + pkin(6)) - t108 * t69 + (t106 * t70 + t107 * t71 + t109) * t68) * g(3) + (-mrSges(1,2) + t107 * (t69 * t96 - t97) - t106 * (-t69 * t98 - t95) + t104 * (t81 * t65 + t83 * t79) + t101 * t83 + t102 * t81) * g(2) + (-m(5) * t64 - mrSges(1,1) + t107 * (t69 * t95 + t98) - t106 * (-t69 * t97 + t96) + t105 * (-t81 * t79 + t64) + (m(5) * t79 - t101) * t81 + t102 * t83) * g(1);
U  = t1;

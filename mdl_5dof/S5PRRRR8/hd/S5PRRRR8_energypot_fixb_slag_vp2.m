% Calculate potential energy for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:42
% EndTime: 2019-12-05 17:14:43
% DurationCPUTime: 0.69s
% Computational Cost: add. (227->91), mult. (378->111), div. (0->0), fcn. (422->12), ass. (0->41)
t125 = -m(5) - m(6);
t124 = m(4) * pkin(7) - mrSges(3,2) + mrSges(4,3);
t123 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t101 = cos(qJ(5));
t98 = sin(qJ(5));
t122 = -m(6) * pkin(4) - mrSges(6,1) * t101 + mrSges(6,2) * t98 - mrSges(5,1);
t121 = -mrSges(6,1) * t98 - mrSges(6,2) * t101 - mrSges(5,3) - t124;
t99 = sin(qJ(3));
t120 = pkin(3) * t99;
t94 = sin(pkin(10));
t95 = sin(pkin(5));
t119 = t94 * t95;
t96 = cos(pkin(10));
t118 = t95 * t96;
t117 = t95 * t99;
t116 = t96 * pkin(1) + pkin(6) * t119;
t100 = sin(qJ(2));
t115 = t100 * t95;
t97 = cos(pkin(5));
t114 = t100 * t97;
t102 = cos(qJ(3));
t113 = t102 * t95;
t103 = cos(qJ(2));
t112 = t103 * t95;
t111 = t103 * t97;
t110 = t97 * pkin(6) + qJ(1);
t109 = t94 * t117;
t90 = t94 * pkin(1);
t108 = -pkin(6) * t118 + t90;
t104 = -pkin(8) - pkin(7);
t87 = pkin(3) * t102 + pkin(2);
t106 = t104 * t112 + t87 * t115 + t97 * t120 + t110;
t93 = qJ(3) + qJ(4);
t89 = cos(t93);
t88 = sin(t93);
t78 = t103 * t96 - t114 * t94;
t77 = t100 * t96 + t111 * t94;
t76 = t103 * t94 + t114 * t96;
t75 = t100 * t94 - t111 * t96;
t72 = t115 * t89 + t88 * t97;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(5) * t106 - t72 * mrSges(5,1) + mrSges(5,3) * t112 - m(6) * (pkin(4) * t72 + t106) - (t101 * t72 - t112 * t98) * mrSges(6,1) - (-t101 * t112 - t72 * t98) * mrSges(6,2) + (-t99 * mrSges(4,1) - t102 * mrSges(4,2) - mrSges(3,3)) * t97 + (t124 * t103 + (-m(4) * pkin(2) - mrSges(4,1) * t102 + mrSges(4,2) * t99 - mrSges(3,1)) * t100) * t95 + t123 * (t115 * t88 - t97 * t89) + (-m(3) - m(4)) * t110) * g(3) + (-mrSges(1,2) - t94 * mrSges(2,1) - t96 * mrSges(2,2) - m(3) * t108 - t76 * mrSges(3,1) + mrSges(3,3) * t118 - m(4) * (pkin(2) * t76 + t108) - (t102 * t76 - t117 * t96) * mrSges(4,1) - (-t113 * t96 - t76 * t99) * mrSges(4,2) + t125 * (t76 * t87 - t75 * t104 + t90 + (-pkin(6) - t120) * t118) + t123 * (t118 * t89 + t76 * t88) + t122 * (-t118 * t88 + t76 * t89) + t121 * t75) * g(2) + (-mrSges(1,1) - t96 * mrSges(2,1) + t94 * mrSges(2,2) - m(3) * t116 - t78 * mrSges(3,1) - mrSges(3,3) * t119 - m(4) * (pkin(2) * t78 + t116) - (t102 * t78 + t109) * mrSges(4,1) - (t113 * t94 - t78 * t99) * mrSges(4,2) + t125 * (pkin(3) * t109 - t77 * t104 + t78 * t87 + t116) + t123 * (-t119 * t89 + t78 * t88) + t122 * (t119 * t88 + t78 * t89) + t121 * t77) * g(1);
U = t1;

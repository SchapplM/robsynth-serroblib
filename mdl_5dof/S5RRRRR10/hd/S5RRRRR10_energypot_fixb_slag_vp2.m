% Calculate potential energy for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:31:54
% EndTime: 2019-12-31 22:31:55
% DurationCPUTime: 0.66s
% Computational Cost: add. (227->91), mult. (378->107), div. (0->0), fcn. (422->12), ass. (0->41)
t125 = -m(5) - m(6);
t124 = m(4) * pkin(8) - mrSges(3,2) + mrSges(4,3);
t123 = -m(6) * pkin(10) + mrSges(5,2) - mrSges(6,3);
t100 = cos(qJ(5));
t96 = sin(qJ(5));
t122 = -m(6) * pkin(4) - mrSges(6,1) * t100 + mrSges(6,2) * t96 - mrSges(5,1);
t121 = -mrSges(6,1) * t96 - mrSges(6,2) * t100 - mrSges(5,3) - t124;
t97 = sin(qJ(3));
t120 = pkin(3) * t97;
t95 = cos(pkin(5));
t119 = t95 * pkin(7) + pkin(6);
t94 = sin(pkin(5));
t98 = sin(qJ(2));
t118 = t94 * t98;
t99 = sin(qJ(1));
t117 = t94 * t99;
t116 = t98 * t99;
t103 = cos(qJ(1));
t115 = t103 * pkin(1) + pkin(7) * t117;
t102 = cos(qJ(2));
t114 = t102 * t94;
t113 = t102 * t99;
t112 = t103 * t94;
t111 = t103 * t98;
t110 = t102 * t103;
t109 = t97 * t117;
t91 = t99 * pkin(1);
t108 = -pkin(7) * t112 + t91;
t104 = -pkin(9) - pkin(8);
t101 = cos(qJ(3));
t87 = pkin(3) * t101 + pkin(2);
t107 = t104 * t114 + t87 * t118 + t95 * t120 + t119;
t93 = qJ(3) + qJ(4);
t89 = cos(t93);
t88 = sin(t93);
t78 = -t95 * t116 + t110;
t77 = t95 * t113 + t111;
t76 = t95 * t111 + t113;
t75 = -t95 * t110 + t116;
t72 = t89 * t118 + t88 * t95;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - m(5) * t107 - t72 * mrSges(5,1) + mrSges(5,3) * t114 - m(6) * (pkin(4) * t72 + t107) - (t100 * t72 - t96 * t114) * mrSges(6,1) - (-t100 * t114 - t72 * t96) * mrSges(6,2) + (-t97 * mrSges(4,1) - t101 * mrSges(4,2) - mrSges(3,3)) * t95 + ((-m(4) * pkin(2) - t101 * mrSges(4,1) + t97 * mrSges(4,2) - mrSges(3,1)) * t98 + t124 * t102) * t94 + t123 * (t88 * t118 - t95 * t89) + (-m(3) - m(4)) * t119) * g(3) + (-mrSges(1,2) - t99 * mrSges(2,1) - t103 * mrSges(2,2) - m(3) * t108 - t76 * mrSges(3,1) + mrSges(3,3) * t112 - m(4) * (pkin(2) * t76 + t108) - (t101 * t76 - t112 * t97) * mrSges(4,1) - (-t101 * t112 - t76 * t97) * mrSges(4,2) + t125 * (t76 * t87 - t75 * t104 + t91 + (-pkin(7) - t120) * t112) + t123 * (t112 * t89 + t76 * t88) + t122 * (-t112 * t88 + t76 * t89) + t121 * t75) * g(2) + (-mrSges(1,1) - t103 * mrSges(2,1) + t99 * mrSges(2,2) - m(3) * t115 - t78 * mrSges(3,1) - mrSges(3,3) * t117 - m(4) * (pkin(2) * t78 + t115) - (t101 * t78 + t109) * mrSges(4,1) - (t101 * t117 - t78 * t97) * mrSges(4,2) + t125 * (pkin(3) * t109 - t77 * t104 + t78 * t87 + t115) + t123 * (-t117 * t89 + t78 * t88) + t122 * (t117 * t88 + t78 * t89) + t121 * t77) * g(1);
U = t1;

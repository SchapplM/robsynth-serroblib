% Calculate potential energy for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR14_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR14_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:41
% EndTime: 2019-12-31 20:35:42
% DurationCPUTime: 0.66s
% Computational Cost: add. (227->91), mult. (378->107), div. (0->0), fcn. (422->12), ass. (0->41)
t125 = -m(5) - m(6);
t124 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t123 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t102 = cos(qJ(5));
t99 = sin(qJ(5));
t122 = -m(6) * pkin(4) - mrSges(6,1) * t102 + mrSges(6,2) * t99 - mrSges(5,1);
t121 = -mrSges(6,1) * t99 - mrSges(6,2) * t102 - mrSges(5,3) - t123;
t94 = sin(pkin(10));
t120 = pkin(3) * t94;
t97 = cos(pkin(5));
t119 = t97 * pkin(7) + pkin(6);
t104 = cos(qJ(1));
t101 = sin(qJ(1));
t95 = sin(pkin(5));
t116 = t101 * t95;
t118 = t104 * pkin(1) + pkin(7) * t116;
t100 = sin(qJ(2));
t117 = t100 * t95;
t103 = cos(qJ(2));
t115 = t103 * t95;
t114 = t104 * t95;
t113 = t100 * t101;
t112 = t100 * t104;
t111 = t101 * t103;
t110 = t103 * t104;
t109 = t94 * t116;
t91 = t101 * pkin(1);
t108 = -pkin(7) * t114 + t91;
t96 = cos(pkin(10));
t87 = pkin(3) * t96 + pkin(2);
t98 = -pkin(8) - qJ(3);
t107 = t98 * t115 + t87 * t117 + t97 * t120 + t119;
t93 = pkin(10) + qJ(4);
t89 = cos(t93);
t88 = sin(t93);
t78 = -t113 * t97 + t110;
t77 = t111 * t97 + t112;
t76 = t112 * t97 + t111;
t75 = -t110 * t97 + t113;
t72 = t117 * t89 + t88 * t97;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - m(5) * t107 - t72 * mrSges(5,1) + mrSges(5,3) * t115 - m(6) * (pkin(4) * t72 + t107) - (t102 * t72 - t115 * t99) * mrSges(6,1) - (-t102 * t115 - t72 * t99) * mrSges(6,2) + (-t94 * mrSges(4,1) - t96 * mrSges(4,2) - mrSges(3,3)) * t97 + (t123 * t103 + (-m(4) * pkin(2) - mrSges(4,1) * t96 + mrSges(4,2) * t94 - mrSges(3,1)) * t100) * t95 + t124 * (t117 * t88 - t97 * t89) + (-m(3) - m(4)) * t119) * g(3) + (-mrSges(1,2) - t101 * mrSges(2,1) - t104 * mrSges(2,2) - m(3) * t108 - t76 * mrSges(3,1) + mrSges(3,3) * t114 - m(4) * (pkin(2) * t76 + t108) - (-t114 * t94 + t76 * t96) * mrSges(4,1) - (-t114 * t96 - t76 * t94) * mrSges(4,2) + t125 * (t76 * t87 - t75 * t98 + t91 + (-pkin(7) - t120) * t114) + t124 * (t114 * t89 + t76 * t88) + t122 * (-t114 * t88 + t76 * t89) + t121 * t75) * g(2) + (-mrSges(1,1) - t104 * mrSges(2,1) + t101 * mrSges(2,2) - m(3) * t118 - t78 * mrSges(3,1) - mrSges(3,3) * t116 - m(4) * (pkin(2) * t78 + t118) - (t78 * t96 + t109) * mrSges(4,1) - (t116 * t96 - t78 * t94) * mrSges(4,2) + t125 * (pkin(3) * t109 - t77 * t98 + t78 * t87 + t118) + t124 * (-t116 * t89 + t78 * t88) + t122 * (t116 * t88 + t78 * t89) + t121 * t77) * g(1);
U = t1;

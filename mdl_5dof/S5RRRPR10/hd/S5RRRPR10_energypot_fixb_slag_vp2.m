% Calculate potential energy for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:32
% EndTime: 2019-12-31 21:26:33
% DurationCPUTime: 0.64s
% Computational Cost: add. (227->91), mult. (378->107), div. (0->0), fcn. (422->12), ass. (0->41)
t125 = -m(5) - m(6);
t124 = m(4) * pkin(8) - mrSges(3,2) + mrSges(4,3);
t123 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t101 = cos(qJ(5));
t97 = sin(qJ(5));
t122 = -m(6) * pkin(4) - mrSges(6,1) * t101 + mrSges(6,2) * t97 - mrSges(5,1);
t121 = -mrSges(6,1) * t97 - mrSges(6,2) * t101 - mrSges(5,3) - t124;
t98 = sin(qJ(3));
t120 = pkin(3) * t98;
t95 = cos(pkin(5));
t119 = t95 * pkin(7) + pkin(6);
t94 = sin(pkin(5));
t99 = sin(qJ(2));
t118 = t94 * t99;
t104 = cos(qJ(1));
t100 = sin(qJ(1));
t116 = t100 * t94;
t117 = t104 * pkin(1) + pkin(7) * t116;
t115 = t100 * t99;
t103 = cos(qJ(2));
t114 = t103 * t94;
t113 = t104 * t94;
t112 = t104 * t99;
t111 = t100 * t103;
t110 = t103 * t104;
t109 = t98 * t116;
t91 = t100 * pkin(1);
t108 = -pkin(7) * t113 + t91;
t102 = cos(qJ(3));
t87 = pkin(3) * t102 + pkin(2);
t96 = -qJ(4) - pkin(8);
t107 = t96 * t114 + t87 * t118 + t95 * t120 + t119;
t93 = qJ(3) + pkin(10);
t89 = cos(t93);
t88 = sin(t93);
t78 = -t115 * t95 + t110;
t77 = t111 * t95 + t112;
t76 = t112 * t95 + t111;
t75 = -t110 * t95 + t115;
t72 = t118 * t89 + t88 * t95;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - m(5) * t107 - t72 * mrSges(5,1) + mrSges(5,3) * t114 - m(6) * (pkin(4) * t72 + t107) - (t101 * t72 - t114 * t97) * mrSges(6,1) - (-t101 * t114 - t72 * t97) * mrSges(6,2) + (-t98 * mrSges(4,1) - t102 * mrSges(4,2) - mrSges(3,3)) * t95 + ((-m(4) * pkin(2) - t102 * mrSges(4,1) + t98 * mrSges(4,2) - mrSges(3,1)) * t99 + t124 * t103) * t94 + t123 * (t118 * t88 - t95 * t89) + (-m(3) - m(4)) * t119) * g(3) + (-mrSges(1,2) - t100 * mrSges(2,1) - t104 * mrSges(2,2) - m(3) * t108 - t76 * mrSges(3,1) + mrSges(3,3) * t113 - m(4) * (t76 * pkin(2) + t108) - (t76 * t102 - t113 * t98) * mrSges(4,1) - (-t102 * t113 - t76 * t98) * mrSges(4,2) + t125 * (t76 * t87 - t75 * t96 + t91 + (-pkin(7) - t120) * t113) + t123 * (t113 * t89 + t76 * t88) + t122 * (-t113 * t88 + t76 * t89) + t121 * t75) * g(2) + (-mrSges(1,1) - t104 * mrSges(2,1) + t100 * mrSges(2,2) - m(3) * t117 - t78 * mrSges(3,1) - mrSges(3,3) * t116 - m(4) * (pkin(2) * t78 + t117) - (t102 * t78 + t109) * mrSges(4,1) - (t102 * t116 - t78 * t98) * mrSges(4,2) + t125 * (pkin(3) * t109 - t77 * t96 + t78 * t87 + t117) + t123 * (-t116 * t89 + t78 * t88) + t122 * (t116 * t88 + t78 * t89) + t121 * t77) * g(1);
U = t1;

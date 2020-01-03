% Calculate potential energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:04
% EndTime: 2019-12-31 22:08:04
% DurationCPUTime: 0.49s
% Computational Cost: add. (207->85), mult. (442->103), div. (0->0), fcn. (511->10), ass. (0->42)
t124 = -mrSges(5,1) - mrSges(6,1);
t123 = -mrSges(5,2) - mrSges(6,2);
t94 = sin(qJ(4));
t107 = pkin(4) * t94 + pkin(8);
t122 = -m(6) * t107 + mrSges(3,2) - mrSges(4,3);
t121 = -m(5) * pkin(9) + m(6) * (-qJ(5) - pkin(9)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t92 = cos(pkin(5));
t120 = t92 * pkin(7) + pkin(6);
t91 = sin(pkin(5));
t96 = sin(qJ(2));
t119 = t91 * t96;
t97 = sin(qJ(1));
t118 = t91 * t97;
t99 = cos(qJ(3));
t117 = t91 * t99;
t116 = t96 * t97;
t101 = cos(qJ(1));
t115 = t101 * pkin(1) + pkin(7) * t118;
t100 = cos(qJ(2));
t114 = t100 * t91;
t113 = t101 * t91;
t112 = t101 * t96;
t111 = t97 * t100;
t110 = t100 * t101;
t109 = pkin(2) * t119 + t120;
t81 = -t92 * t116 + t110;
t108 = t81 * pkin(2) + t115;
t106 = t97 * pkin(1) - pkin(7) * t113;
t80 = t92 * t111 + t112;
t105 = pkin(8) * t80 + t108;
t79 = t92 * t112 + t111;
t104 = t79 * pkin(2) + t106;
t103 = -pkin(8) * t114 + t109;
t78 = -t92 * t110 + t116;
t102 = t78 * pkin(8) + t104;
t98 = cos(qJ(4));
t95 = sin(qJ(3));
t87 = pkin(4) * t98 + pkin(3);
t77 = t96 * t117 + t92 * t95;
t73 = t95 * t118 + t81 * t99;
t71 = -t95 * t113 + t79 * t99;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - m(3) * t120 - t92 * mrSges(3,3) - (t96 * mrSges(3,1) + t100 * mrSges(3,2)) * t91 - m(4) * t103 - t77 * mrSges(4,1) + mrSges(4,3) * t114 - m(5) * (pkin(3) * t77 + t103) - m(6) * (-t107 * t114 + t77 * t87 + t109) + t124 * (-t94 * t114 + t77 * t98) + t123 * (-t98 * t114 - t77 * t94) + t121 * (t95 * t119 - t92 * t99)) * g(3) + (-mrSges(1,2) - t97 * mrSges(2,1) - t101 * mrSges(2,2) - m(3) * t106 - t79 * mrSges(3,1) + mrSges(3,3) * t113 - m(4) * t102 - t71 * mrSges(4,1) - m(5) * (t71 * pkin(3) + t102) - m(6) * (t71 * t87 + t104) + t122 * t78 + t124 * (t71 * t98 + t78 * t94) + t123 * (-t71 * t94 + t78 * t98) + t121 * (t99 * t113 + t79 * t95)) * g(2) + (-mrSges(1,1) - t101 * mrSges(2,1) + t97 * mrSges(2,2) - m(3) * t115 - t81 * mrSges(3,1) - mrSges(3,3) * t118 - m(4) * t105 - t73 * mrSges(4,1) - m(5) * (pkin(3) * t73 + t105) - m(6) * (t73 * t87 + t108) + t122 * t80 + t124 * (t73 * t98 + t80 * t94) + t123 * (-t73 * t94 + t80 * t98) + t121 * (-t97 * t117 + t81 * t95)) * g(1);
U = t1;

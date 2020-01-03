% Calculate potential energy for
% S5RRRRP11
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
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRP11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:21
% EndTime: 2019-12-31 22:14:21
% DurationCPUTime: 0.42s
% Computational Cost: add. (222->77), mult. (488->95), div. (0->0), fcn. (575->10), ass. (0->40)
t129 = -m(5) - m(6);
t128 = mrSges(3,2) - mrSges(4,3);
t127 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t126 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t125 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t99 = cos(pkin(5));
t124 = t99 * pkin(7) + pkin(6);
t107 = cos(qJ(1));
t103 = sin(qJ(1));
t98 = sin(pkin(5));
t122 = t103 * t98;
t123 = t107 * pkin(1) + pkin(7) * t122;
t121 = t107 * t98;
t102 = sin(qJ(2));
t120 = t98 * t102;
t106 = cos(qJ(2));
t119 = t98 * t106;
t118 = t102 * t103;
t117 = t102 * t107;
t116 = t103 * t106;
t115 = t106 * t107;
t114 = t103 * pkin(1) - pkin(7) * t121;
t88 = t99 * t116 + t117;
t89 = -t99 * t118 + t115;
t113 = t89 * pkin(2) + pkin(8) * t88 + t123;
t112 = pkin(2) * t120 - pkin(8) * t119 + t124;
t86 = -t99 * t115 + t118;
t87 = t99 * t117 + t116;
t111 = t87 * pkin(2) + t86 * pkin(8) + t114;
t105 = cos(qJ(3));
t104 = cos(qJ(4));
t101 = sin(qJ(3));
t100 = sin(qJ(4));
t85 = t101 * t99 + t105 * t120;
t84 = t101 * t120 - t99 * t105;
t78 = t101 * t122 + t105 * t89;
t77 = t101 * t89 - t105 * t122;
t76 = -t101 * t121 + t87 * t105;
t75 = t87 * t101 + t105 * t121;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - m(3) * t124 - t99 * mrSges(3,3) - (t102 * mrSges(3,1) + t106 * mrSges(3,2)) * t98 - m(4) * t112 - t85 * mrSges(4,1) + mrSges(4,3) * t119 + t129 * (t85 * pkin(3) + pkin(9) * t84 + t112) + t126 * (-t100 * t119 + t104 * t85) + t125 * (t100 * t85 + t104 * t119) + t127 * t84) * g(3) + (-m(3) * t114 - m(4) * t111 - t103 * mrSges(2,1) - t87 * mrSges(3,1) - t76 * mrSges(4,1) - mrSges(2,2) * t107 + mrSges(3,3) * t121 - mrSges(1,2) + t128 * t86 + t129 * (t76 * pkin(3) + t75 * pkin(9) + t111) + t126 * (t100 * t86 + t104 * t76) + t125 * (t100 * t76 - t86 * t104) + t127 * t75) * g(2) + (-m(3) * t123 - m(4) * t113 - mrSges(2,1) * t107 - t89 * mrSges(3,1) - t78 * mrSges(4,1) + t103 * mrSges(2,2) - mrSges(3,3) * t122 - mrSges(1,1) + t128 * t88 + t129 * (t78 * pkin(3) + pkin(9) * t77 + t113) + t126 * (t100 * t88 + t104 * t78) + t125 * (t100 * t78 - t88 * t104) + t127 * t77) * g(1);
U = t1;

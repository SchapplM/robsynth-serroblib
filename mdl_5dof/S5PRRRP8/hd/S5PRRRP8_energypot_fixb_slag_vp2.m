% Calculate potential energy for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP8_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:08
% EndTime: 2019-12-05 16:58:09
% DurationCPUTime: 0.45s
% Computational Cost: add. (222->77), mult. (488->99), div. (0->0), fcn. (575->10), ass. (0->40)
t128 = -m(5) - m(6);
t127 = mrSges(3,2) - mrSges(4,3);
t126 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t125 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t124 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t97 = sin(pkin(9));
t98 = sin(pkin(5));
t123 = t97 * t98;
t99 = cos(pkin(9));
t122 = t99 * t98;
t121 = t99 * pkin(1) + pkin(6) * t123;
t102 = sin(qJ(3));
t120 = t102 * t98;
t105 = cos(qJ(3));
t119 = t105 * t98;
t103 = sin(qJ(2));
t118 = t98 * t103;
t106 = cos(qJ(2));
t117 = t98 * t106;
t100 = cos(pkin(5));
t116 = t100 * pkin(6) + qJ(1);
t115 = t100 * t103;
t114 = t100 * t106;
t113 = t97 * pkin(1) - pkin(6) * t122;
t85 = t99 * t103 + t97 * t114;
t86 = t106 * t99 - t97 * t115;
t112 = t86 * pkin(2) + pkin(7) * t85 + t121;
t111 = pkin(2) * t118 - pkin(7) * t117 + t116;
t83 = t103 * t97 - t99 * t114;
t84 = t106 * t97 + t99 * t115;
t110 = t84 * pkin(2) + pkin(7) * t83 + t113;
t104 = cos(qJ(4));
t101 = sin(qJ(4));
t88 = t100 * t102 + t105 * t118;
t87 = -t100 * t105 + t102 * t118;
t75 = t105 * t86 + t97 * t120;
t74 = t102 * t86 - t97 * t119;
t73 = t105 * t84 - t99 * t120;
t72 = t102 * t84 + t99 * t119;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t116 - t100 * mrSges(3,3) - (t103 * mrSges(3,1) + t106 * mrSges(3,2)) * t98 - m(4) * t111 - t88 * mrSges(4,1) + mrSges(4,3) * t117 + t128 * (t88 * pkin(3) + t87 * pkin(8) + t111) + t125 * (-t101 * t117 + t88 * t104) + t124 * (t88 * t101 + t104 * t117) + t126 * t87) * g(3) + (-m(3) * t113 - m(4) * t110 - mrSges(2,1) * t97 - t84 * mrSges(3,1) - t73 * mrSges(4,1) - mrSges(2,2) * t99 + mrSges(3,3) * t122 - mrSges(1,2) + t127 * t83 + t128 * (t73 * pkin(3) + pkin(8) * t72 + t110) + t125 * (t101 * t83 + t104 * t73) + t124 * (t101 * t73 - t83 * t104) + t126 * t72) * g(2) + (-m(3) * t121 - m(4) * t112 - mrSges(2,1) * t99 - t86 * mrSges(3,1) - t75 * mrSges(4,1) + mrSges(2,2) * t97 - mrSges(3,3) * t123 - mrSges(1,1) + t127 * t85 + t128 * (t75 * pkin(3) + pkin(8) * t74 + t112) + t125 * (t101 * t85 + t104 * t75) + t124 * (t101 * t75 - t85 * t104) + t126 * t74) * g(1);
U = t1;

% Calculate potential energy for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:34:54
% EndTime: 2019-12-05 16:34:54
% DurationCPUTime: 0.49s
% Computational Cost: add. (243->80), mult. (543->102), div. (0->0), fcn. (651->12), ass. (0->41)
t131 = -m(5) - m(6);
t130 = mrSges(3,2) - mrSges(4,3);
t129 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t128 = -t105 * mrSges(6,1) - t108 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t127 = -m(6) * pkin(4) - t108 * mrSges(6,1) + t105 * mrSges(6,2) - mrSges(5,1);
t103 = cos(pkin(9));
t100 = sin(pkin(9));
t101 = sin(pkin(5));
t124 = t100 * t101;
t126 = t103 * pkin(1) + pkin(6) * t124;
t104 = cos(pkin(5));
t125 = t104 * pkin(6) + qJ(1);
t107 = sin(qJ(2));
t123 = t101 * t107;
t109 = cos(qJ(3));
t122 = t101 * t109;
t110 = cos(qJ(2));
t121 = t101 * t110;
t120 = t103 * t101;
t119 = t104 * t107;
t118 = t104 * t110;
t117 = t100 * pkin(1) - pkin(6) * t120;
t87 = t100 * t118 + t103 * t107;
t88 = -t100 * t119 + t103 * t110;
t116 = t88 * pkin(2) + pkin(7) * t87 + t126;
t115 = pkin(2) * t123 - pkin(7) * t121 + t125;
t85 = t100 * t107 - t103 * t118;
t86 = t100 * t110 + t103 * t119;
t114 = t86 * pkin(2) + pkin(7) * t85 + t117;
t106 = sin(qJ(3));
t102 = cos(pkin(10));
t99 = sin(pkin(10));
t90 = t104 * t106 + t107 * t122;
t89 = -t104 * t109 + t106 * t123;
t79 = t106 * t124 + t109 * t88;
t78 = -t100 * t122 + t106 * t88;
t77 = -t106 * t120 + t109 * t86;
t76 = t106 * t86 + t109 * t120;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t125 - t104 * mrSges(3,3) - (t107 * mrSges(3,1) + t110 * mrSges(3,2)) * t101 - m(4) * t115 - t90 * mrSges(4,1) + mrSges(4,3) * t121 + t131 * (t90 * pkin(3) + t89 * qJ(4) + t115) + t127 * (t90 * t102 - t99 * t121) + t128 * t89 + t129 * (t102 * t121 + t90 * t99)) * g(3) + (-m(3) * t117 - m(4) * t114 - t100 * mrSges(2,1) - t86 * mrSges(3,1) - t77 * mrSges(4,1) - t103 * mrSges(2,2) + mrSges(3,3) * t120 - mrSges(1,2) + t130 * t85 + t131 * (t77 * pkin(3) + qJ(4) * t76 + t114) + t127 * (t102 * t77 + t85 * t99) + t128 * t76 + t129 * (-t85 * t102 + t77 * t99)) * g(2) + (-m(3) * t126 - m(4) * t116 - t103 * mrSges(2,1) - t88 * mrSges(3,1) - t79 * mrSges(4,1) + t100 * mrSges(2,2) - mrSges(3,3) * t124 - mrSges(1,1) + t130 * t87 + t131 * (t79 * pkin(3) + qJ(4) * t78 + t116) + t127 * (t102 * t79 + t87 * t99) + t128 * t78 + t129 * (-t87 * t102 + t79 * t99)) * g(1);
U = t1;

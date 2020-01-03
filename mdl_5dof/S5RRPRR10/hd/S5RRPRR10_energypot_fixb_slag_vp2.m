% Calculate potential energy for
% S5RRPRR10
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPRR10_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:42
% EndTime: 2019-12-31 20:23:42
% DurationCPUTime: 0.56s
% Computational Cost: add. (253->79), mult. (556->102), div. (0->0), fcn. (668->12), ass. (0->44)
t134 = -m(5) - m(6);
t133 = -mrSges(3,3) - mrSges(4,3);
t132 = -m(3) * pkin(1) - mrSges(2,1);
t131 = m(3) * pkin(7) - t133;
t130 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t103 = sin(qJ(5));
t107 = cos(qJ(5));
t129 = t103 * mrSges(6,1) + t107 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t128 = -m(6) * pkin(4) - t107 * mrSges(6,1) + t103 * mrSges(6,2) - mrSges(5,1);
t102 = cos(pkin(5));
t127 = t102 * pkin(7) + pkin(6);
t106 = sin(qJ(1));
t110 = cos(qJ(1));
t100 = sin(pkin(5));
t105 = sin(qJ(2));
t87 = pkin(2) * t102 * t105 + (-pkin(7) - qJ(3)) * t100;
t109 = cos(qJ(2));
t96 = pkin(2) * t109 + pkin(1);
t126 = t106 * t96 + t110 * t87;
t125 = t100 * t105;
t124 = t100 * t106;
t123 = t100 * t110;
t101 = cos(pkin(10));
t122 = t101 * t109;
t121 = t105 * t110;
t120 = t106 * t105;
t119 = t106 * t109;
t118 = t109 * t110;
t117 = pkin(2) * t125 + t102 * qJ(3) + t127;
t116 = -t106 * t87 + t110 * t96;
t99 = sin(pkin(10));
t89 = -t105 * t99 + t122;
t114 = t101 * t105 + t109 * t99;
t111 = t89 * t102;
t108 = cos(qJ(4));
t104 = sin(qJ(4));
t86 = t114 * t102;
t85 = t114 * t100;
t84 = -t100 * t122 + t99 * t125;
t78 = -t106 * t86 + t110 * t89;
t77 = -t106 * t111 - t110 * t114;
t76 = t106 * t89 + t110 * t86;
t75 = -t106 * t114 + t110 * t111;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - m(3) * t127 - (t105 * mrSges(3,1) + t109 * mrSges(3,2)) * t100 - m(4) * t117 - t85 * mrSges(4,1) + t134 * (t85 * pkin(3) + pkin(8) * t84 + t117) + t128 * (t102 * t104 + t108 * t85) - t129 * t84 + t130 * (-t102 * t108 + t104 * t85) + t133 * t102) * g(3) + (-mrSges(1,2) - mrSges(2,2) * t110 - (t102 * t121 + t119) * mrSges(3,1) - (t102 * t118 - t120) * mrSges(3,2) - m(4) * t126 - t76 * mrSges(4,1) + t134 * (t76 * pkin(3) - pkin(8) * t75 + t126) + t128 * (-t104 * t123 + t76 * t108) + t129 * t75 + t130 * (t76 * t104 + t108 * t123) + t132 * t106 + t131 * t123) * g(2) + (-mrSges(1,1) + t106 * mrSges(2,2) - (-t102 * t120 + t118) * mrSges(3,1) - (-t102 * t119 - t121) * mrSges(3,2) - m(4) * t116 - t78 * mrSges(4,1) + t134 * (t78 * pkin(3) - pkin(8) * t77 + t116) + t128 * (t104 * t124 + t108 * t78) + t129 * t77 + t130 * (t104 * t78 - t108 * t124) + t132 * t110 - t131 * t124) * g(1);
U = t1;

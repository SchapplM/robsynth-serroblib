% Calculate potential energy for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_energypot_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:42:52
% EndTime: 2019-03-08 19:42:53
% DurationCPUTime: 0.70s
% Computational Cost: add. (332->97), mult. (532->110), div. (0->0), fcn. (603->12), ass. (0->47)
t137 = -m(6) - m(7);
t136 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t135 = -m(7) * pkin(9) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t111 = sin(qJ(6));
t113 = cos(qJ(6));
t134 = -t111 * mrSges(7,1) - t113 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t133 = m(7) * pkin(5) + t113 * mrSges(7,1) - t111 * mrSges(7,2) + mrSges(6,1) + mrSges(5,3);
t132 = -t133 - t136;
t104 = sin(pkin(11));
t131 = pkin(3) * t104;
t108 = cos(pkin(10));
t105 = sin(pkin(10));
t106 = sin(pkin(6));
t129 = t105 * t106;
t130 = t108 * pkin(1) + pkin(7) * t129;
t128 = t106 * t108;
t112 = sin(qJ(2));
t127 = t106 * t112;
t114 = cos(qJ(2));
t126 = t106 * t114;
t109 = cos(pkin(6));
t125 = t109 * t112;
t124 = t109 * t114;
t123 = t109 * pkin(7) + qJ(1);
t122 = t104 * t129;
t100 = t105 * pkin(1);
t121 = -pkin(7) * t128 + t100;
t110 = -pkin(8) - qJ(3);
t86 = t105 * t124 + t108 * t112;
t87 = -t105 * t125 + t108 * t114;
t107 = cos(pkin(11));
t97 = pkin(3) * t107 + pkin(2);
t120 = pkin(3) * t122 - t86 * t110 + t87 * t97 + t130;
t119 = t109 * t131 + t110 * t126 + t97 * t127 + t123;
t84 = t105 * t112 - t108 * t124;
t85 = t105 * t114 + t108 * t125;
t116 = t100 + t85 * t97 - t84 * t110 + (-pkin(7) - t131) * t128;
t103 = pkin(11) + qJ(4);
t99 = cos(t103);
t98 = sin(t103);
t81 = t109 * t98 + t127 * t99;
t80 = -t109 * t99 + t127 * t98;
t76 = t129 * t98 + t87 * t99;
t75 = -t129 * t99 + t87 * t98;
t74 = -t128 * t98 + t85 * t99;
t73 = t128 * t99 + t85 * t98;
t1 = (-m(2) * qJ(1) - m(5) * t119 - mrSges(1,3) - mrSges(2,3) + t137 * (t81 * pkin(4) + t80 * qJ(5) + t119) + t134 * t80 + t133 * t126 + (-m(3) - m(4)) * t123 + (-t104 * mrSges(4,1) - t107 * mrSges(4,2) - mrSges(3,3)) * t109 + (t136 * t114 + (-m(4) * pkin(2) - t107 * mrSges(4,1) + t104 * mrSges(4,2) - mrSges(3,1)) * t112) * t106 + t135 * t81) * g(3) + (-mrSges(1,2) - t105 * mrSges(2,1) - t108 * mrSges(2,2) - m(3) * t121 - t85 * mrSges(3,1) + mrSges(3,3) * t128 - m(4) * (pkin(2) * t85 + t121) - (-t104 * t128 + t85 * t107) * mrSges(4,1) - (-t104 * t85 - t107 * t128) * mrSges(4,2) - m(5) * t116 + t137 * (t74 * pkin(4) + qJ(5) * t73 + t116) + t135 * t74 + t134 * t73 + t132 * t84) * g(2) + (-mrSges(1,1) - t108 * mrSges(2,1) + t105 * mrSges(2,2) - m(3) * t130 - t87 * mrSges(3,1) - mrSges(3,3) * t129 - m(4) * (pkin(2) * t87 + t130) - (t107 * t87 + t122) * mrSges(4,1) - (-t104 * t87 + t107 * t129) * mrSges(4,2) - m(5) * t120 + t137 * (t76 * pkin(4) + qJ(5) * t75 + t120) + t135 * t76 + t134 * t75 + t132 * t86) * g(1);
U  = t1;

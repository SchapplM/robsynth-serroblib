% Calculate potential energy for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% mrSges [8x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U = S7RRRRRRR1_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1),zeros(8,1),zeros(8,3)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_energypot_fixb_slag_vp2: qJ has to be [7x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_energypot_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_energypot_fixb_slag_vp2: m has to be [8x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [8,3]), ...
  'S7RRRRRRR1_energypot_fixb_slag_vp2: mrSges has to be [8x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 19:14:16
% EndTime: 2018-11-26 19:14:16
% DurationCPUTime: 0.38s
% Computational Cost: add. (294->85), mult. (698->104), div. (0->0), fcn. (857->14), ass. (0->46)
t128 = -m(4) - m(5);
t127 = -m(6) - m(7);
t110 = cos(qJ(2));
t97 = t110 * pkin(2);
t126 = t97 + pkin(1);
t125 = mrSges(2,2) - mrSges(3,3);
t124 = -mrSges(3,2) - mrSges(4,3);
t123 = mrSges(4,2) + mrSges(5,3);
t102 = sin(qJ(3));
t103 = sin(qJ(2));
t122 = t102 * t103;
t104 = sin(qJ(1));
t121 = t103 * t104;
t109 = cos(qJ(3));
t120 = t103 * t109;
t111 = cos(qJ(1));
t119 = t103 * t111;
t118 = t104 * t110;
t117 = t110 * t111;
t116 = -m(8) * pkin(3) + mrSges(5,2) - mrSges(6,3);
t115 = -m(8) * pkin(4) - mrSges(7,2) - mrSges(8,3);
t105 = cos(qJ(7));
t98 = sin(qJ(7));
t114 = -t105 * mrSges(8,1) + t98 * mrSges(8,2) - mrSges(7,1);
t113 = t98 * mrSges(8,1) + t105 * mrSges(8,2) + mrSges(6,2) - mrSges(7,3);
t112 = -mrSges(3,1) * t110 - mrSges(2,1) + ((m(8) - t128) * pkin(2) - t124) * t103;
t108 = cos(qJ(4));
t107 = cos(qJ(5));
t106 = cos(qJ(6));
t101 = sin(qJ(4));
t100 = sin(qJ(5));
t99 = sin(qJ(6));
t93 = -t104 * t102 + t109 * t117;
t92 = -t102 * t117 - t104 * t109;
t91 = t102 * t111 + t109 * t118;
t90 = -t102 * t118 + t109 * t111;
t89 = -t101 * t110 + t108 * t120;
t88 = t101 * t120 + t108 * t110;
t86 = t101 * t119 + t93 * t108;
t85 = t101 * t93 - t108 * t119;
t84 = t101 * t121 + t108 * t91;
t83 = t101 * t91 - t108 * t121;
t82 = -t100 * t122 + t107 * t89;
t78 = t100 * t92 + t107 * t86;
t76 = t100 * t90 + t107 * t84;
t1 = (-m(8) * t97 - t89 * mrSges(5,1) - t82 * mrSges(6,1) - mrSges(1,3) - mrSges(2,3) + t128 * t126 + t116 * t88 + t127 * (pkin(3) * t88 + t126) + t114 * (t106 * t82 + t88 * t99) + t113 * (t100 * t89 + t107 * t122) + t115 * (t106 * t88 - t82 * t99) + t124 * t110 + (-m(8) - m(2) - m(3)) * pkin(1) + (-mrSges(4,1) * t109 + t123 * t102 - mrSges(3,1)) * t103) * g(3) + (-t91 * mrSges(4,1) - t84 * mrSges(5,1) - t76 * mrSges(6,1) - mrSges(1,2) - t123 * t90 + t127 * (-pkin(2) * t121 + pkin(3) * t83) + t114 * (t106 * t76 + t83 * t99) + t113 * (t100 * t84 - t107 * t90) + t116 * t83 + t115 * (t106 * t83 - t76 * t99) - t125 * t111 + t112 * t104) * g(2) + (-t93 * mrSges(4,1) - t86 * mrSges(5,1) - t78 * mrSges(6,1) - mrSges(1,1) - t123 * t92 + t127 * (-pkin(2) * t119 + t85 * pkin(3)) + t114 * (t106 * t78 + t85 * t99) + t113 * (t100 * t86 - t107 * t92) + t116 * t85 + t115 * (t106 * t85 - t78 * t99) + t125 * t104 + t112 * t111) * g(1);
U  = t1;

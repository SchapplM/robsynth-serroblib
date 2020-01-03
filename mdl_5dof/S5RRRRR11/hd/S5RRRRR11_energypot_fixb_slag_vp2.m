% Calculate potential energy for
% S5RRRRR11
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR11_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:43
% EndTime: 2019-12-31 22:38:44
% DurationCPUTime: 0.53s
% Computational Cost: add. (219->85), mult. (442->100), div. (0->0), fcn. (511->12), ass. (0->40)
t117 = -m(4) - m(5);
t116 = -m(5) * pkin(9) + m(6) * (-pkin(10) - pkin(9)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t87 = sin(qJ(4));
t101 = pkin(4) * t87 + pkin(8);
t84 = qJ(4) + qJ(5);
t79 = sin(t84);
t80 = cos(t84);
t91 = cos(qJ(4));
t115 = -m(6) * t101 - t87 * mrSges(5,1) - t79 * mrSges(6,1) - t91 * mrSges(5,2) - t80 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3);
t78 = t91 * pkin(4) + pkin(3);
t114 = -m(5) * pkin(3) - m(6) * t78 - t91 * mrSges(5,1) - t80 * mrSges(6,1) + t87 * mrSges(5,2) + t79 * mrSges(6,2) - mrSges(4,1);
t86 = cos(pkin(5));
t113 = t86 * pkin(7) + pkin(6);
t85 = sin(pkin(5));
t89 = sin(qJ(2));
t112 = t85 * t89;
t93 = cos(qJ(2));
t111 = t85 * t93;
t90 = sin(qJ(1));
t110 = t90 * t85;
t109 = t90 * t89;
t108 = t90 * t93;
t94 = cos(qJ(1));
t107 = t94 * t85;
t106 = t94 * t89;
t105 = t94 * t93;
t104 = t94 * pkin(1) + pkin(7) * t110;
t103 = pkin(2) * t112 + t113;
t72 = -t86 * t109 + t105;
t102 = t72 * pkin(2) + t104;
t100 = t90 * pkin(1) - pkin(7) * t107;
t70 = t86 * t106 + t108;
t99 = t70 * pkin(2) + t100;
t97 = -pkin(8) * t111 + t103;
t92 = cos(qJ(3));
t88 = sin(qJ(3));
t71 = t86 * t108 + t106;
t69 = -t86 * t105 + t109;
t68 = t92 * t112 + t86 * t88;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - m(3) * t113 - t86 * mrSges(3,3) - (t89 * mrSges(3,1) + t93 * mrSges(3,2)) * t85 - m(4) * t97 - t68 * mrSges(4,1) + mrSges(4,3) * t111 - m(5) * (t68 * pkin(3) + t97) - (-t87 * t111 + t68 * t91) * mrSges(5,1) - (-t91 * t111 - t68 * t87) * mrSges(5,2) - m(6) * (-t101 * t111 + t68 * t78 + t103) - (-t79 * t111 + t68 * t80) * mrSges(6,1) - (-t80 * t111 - t68 * t79) * mrSges(6,2) + t116 * (t88 * t112 - t86 * t92)) * g(3) + (-m(3) * t100 - m(6) * t99 - t90 * mrSges(2,1) - t70 * mrSges(3,1) - t94 * mrSges(2,2) + mrSges(3,3) * t107 - mrSges(1,2) + t117 * (t69 * pkin(8) + t99) + t114 * (-t88 * t107 + t70 * t92) + t115 * t69 + t116 * (t92 * t107 + t70 * t88)) * g(2) + (-m(3) * t104 - m(6) * t102 - t94 * mrSges(2,1) - t72 * mrSges(3,1) + t90 * mrSges(2,2) - mrSges(3,3) * t110 - mrSges(1,1) + t117 * (t71 * pkin(8) + t102) + t114 * (t88 * t110 + t72 * t92) + t115 * t71 + t116 * (-t92 * t110 + t72 * t88)) * g(1);
U = t1;

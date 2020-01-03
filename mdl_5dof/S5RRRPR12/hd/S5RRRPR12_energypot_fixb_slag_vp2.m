% Calculate potential energy for
% S5RRRPR12
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR12_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:36:59
% EndTime: 2019-12-31 21:36:59
% DurationCPUTime: 0.53s
% Computational Cost: add. (219->85), mult. (442->100), div. (0->0), fcn. (511->12), ass. (0->40)
t117 = -m(4) - m(5);
t116 = -m(5) * qJ(4) + m(6) * (-pkin(9) - qJ(4)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t85 = sin(pkin(10));
t101 = pkin(4) * t85 + pkin(8);
t84 = pkin(10) + qJ(5);
t79 = sin(t84);
t80 = cos(t84);
t87 = cos(pkin(10));
t115 = -m(6) * t101 - t85 * mrSges(5,1) - t79 * mrSges(6,1) - t87 * mrSges(5,2) - t80 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3);
t78 = t87 * pkin(4) + pkin(3);
t114 = -m(5) * pkin(3) - m(6) * t78 - t87 * mrSges(5,1) - t80 * mrSges(6,1) + t85 * mrSges(5,2) + t79 * mrSges(6,2) - mrSges(4,1);
t88 = cos(pkin(5));
t113 = t88 * pkin(7) + pkin(6);
t86 = sin(pkin(5));
t91 = sin(qJ(2));
t112 = t86 * t91;
t94 = cos(qJ(2));
t111 = t86 * t94;
t92 = sin(qJ(1));
t110 = t92 * t86;
t109 = t92 * t91;
t108 = t92 * t94;
t95 = cos(qJ(1));
t107 = t95 * t86;
t106 = t95 * t91;
t105 = t95 * t94;
t104 = t95 * pkin(1) + pkin(7) * t110;
t103 = pkin(2) * t112 + t113;
t72 = -t88 * t109 + t105;
t102 = t72 * pkin(2) + t104;
t100 = t92 * pkin(1) - pkin(7) * t107;
t70 = t88 * t106 + t108;
t99 = t70 * pkin(2) + t100;
t97 = -pkin(8) * t111 + t103;
t93 = cos(qJ(3));
t90 = sin(qJ(3));
t71 = t88 * t108 + t106;
t69 = -t88 * t105 + t109;
t68 = t93 * t112 + t88 * t90;
t1 = (-mrSges(1,3) - m(2) * pkin(6) - mrSges(2,3) - m(3) * t113 - t88 * mrSges(3,3) - (t91 * mrSges(3,1) + t94 * mrSges(3,2)) * t86 - m(4) * t97 - t68 * mrSges(4,1) + mrSges(4,3) * t111 - m(5) * (t68 * pkin(3) + t97) - (-t85 * t111 + t68 * t87) * mrSges(5,1) - (-t87 * t111 - t68 * t85) * mrSges(5,2) - m(6) * (-t101 * t111 + t68 * t78 + t103) - (-t79 * t111 + t68 * t80) * mrSges(6,1) - (-t80 * t111 - t68 * t79) * mrSges(6,2) + t116 * (t90 * t112 - t88 * t93)) * g(3) + (-m(3) * t100 - m(6) * t99 - t92 * mrSges(2,1) - t70 * mrSges(3,1) - t95 * mrSges(2,2) + mrSges(3,3) * t107 - mrSges(1,2) + t117 * (t69 * pkin(8) + t99) + t114 * (-t90 * t107 + t70 * t93) + t115 * t69 + t116 * (t93 * t107 + t70 * t90)) * g(2) + (-m(3) * t104 - m(6) * t102 - t95 * mrSges(2,1) - t72 * mrSges(3,1) + t92 * mrSges(2,2) - mrSges(3,3) * t110 - mrSges(1,1) + t117 * (t71 * pkin(8) + t102) + t114 * (t90 * t110 + t72 * t93) + t115 * t71 + t116 * (-t93 * t110 + t72 * t90)) * g(1);
U = t1;

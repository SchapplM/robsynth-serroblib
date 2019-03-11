% Calculate potential energy for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:11:56
% EndTime: 2019-03-09 22:11:56
% DurationCPUTime: 0.47s
% Computational Cost: add. (239->75), mult. (283->73), div. (0->0), fcn. (278->10), ass. (0->40)
t118 = -m(6) - m(7);
t82 = qJ(2) + qJ(3);
t79 = sin(t82);
t84 = sin(qJ(4));
t88 = cos(qJ(4));
t117 = (pkin(4) * t88 + qJ(5) * t84) * t79;
t80 = cos(t82);
t85 = sin(qJ(2));
t89 = cos(qJ(2));
t116 = -m(3) * pkin(1) - t89 * mrSges(3,1) - t80 * mrSges(4,1) + t85 * mrSges(3,2) + t79 * mrSges(4,2) - mrSges(2,1);
t115 = mrSges(6,2) + mrSges(5,3) - mrSges(7,3);
t114 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t113 = m(7) * pkin(10) - t115;
t83 = sin(qJ(6));
t87 = cos(qJ(6));
t112 = -t83 * mrSges(7,1) - t87 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t111 = -m(7) * pkin(5) - t87 * mrSges(7,1) + t83 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t110 = pkin(3) * t80;
t109 = t85 * pkin(2) + pkin(6);
t86 = sin(qJ(1));
t108 = t86 * t79;
t107 = t86 * t84;
t106 = t86 * t88;
t90 = cos(qJ(1));
t105 = t90 * t79;
t104 = t90 * t84;
t103 = t90 * t88;
t77 = t89 * pkin(2) + pkin(1);
t91 = -pkin(8) - pkin(7);
t102 = t86 * t77 + t90 * t91;
t101 = t79 * pkin(3) + t109;
t99 = t90 * t77 - t86 * t91;
t98 = pkin(9) * t108 + t86 * t110 + t102;
t96 = -t80 * pkin(9) + t101;
t95 = pkin(9) * t105 + t90 * t110 + t99;
t66 = t80 * t103 + t107;
t65 = t80 * t104 - t106;
t64 = t80 * t106 - t104;
t63 = t80 * t107 + t103;
t1 = (-mrSges(1,3) - mrSges(2,3) - t85 * mrSges(3,1) - t89 * mrSges(3,2) - m(4) * t109 - m(5) * t96 - m(6) * (t96 + t117) - m(7) * (t101 + t117) + (-m(2) - m(3)) * pkin(6) + (-mrSges(4,2) - m(7) * (-pkin(9) + pkin(10)) + t115) * t80 + (t111 * t88 + t112 * t84 - mrSges(4,1)) * t79) * g(3) + (-m(4) * t102 - m(5) * t98 - mrSges(1,2) + t118 * (t64 * pkin(4) + t63 * qJ(5) + t98) + t111 * t64 + t112 * t63 - t114 * t90 + t116 * t86 + t113 * t108) * g(2) + (-m(4) * t99 - m(5) * t95 - mrSges(1,1) + t118 * (t66 * pkin(4) + t65 * qJ(5) + t95) + t111 * t66 + t112 * t65 + t116 * t90 + t114 * t86 + t113 * t105) * g(1);
U  = t1;

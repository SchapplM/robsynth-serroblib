% Calculate potential energy for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR5_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_energypot_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_energypot_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_energypot_fixb_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:02
% EndTime: 2019-03-09 02:50:03
% DurationCPUTime: 0.50s
% Computational Cost: add. (206->76), mult. (213->69), div. (0->0), fcn. (184->10), ass. (0->34)
t76 = pkin(10) + qJ(6);
t71 = sin(t76);
t73 = cos(t76);
t78 = sin(pkin(10));
t117 = -m(7) * pkin(5) * t78 - t71 * mrSges(7,1) - t73 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t116 = -mrSges(4,1) + mrSges(5,2) + m(7) * (-pkin(8) - qJ(5)) - mrSges(7,3) - m(6) * qJ(5);
t114 = -m(6) - m(7);
t77 = pkin(9) + qJ(3);
t74 = cos(t77);
t85 = cos(qJ(1));
t104 = t74 * t85;
t72 = sin(t77);
t98 = qJ(4) * t72;
t113 = pkin(3) * t104 + t85 * t98;
t112 = m(5) - t114;
t109 = -m(3) * qJ(2) - t73 * mrSges(7,1) + t71 * mrSges(7,2) - mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t79 = sin(pkin(9));
t81 = cos(pkin(9));
t108 = -m(3) * pkin(1) - t81 * mrSges(3,1) + t79 * mrSges(3,2) + t116 * t74 + t117 * t72 - mrSges(2,1);
t106 = pkin(2) * t79 + pkin(6);
t84 = sin(qJ(1));
t105 = t74 * t84;
t103 = t78 * t85;
t80 = cos(pkin(10));
t102 = t80 * t85;
t101 = t84 * t78;
t100 = t84 * t80;
t69 = pkin(2) * t81 + pkin(1);
t83 = -pkin(7) - qJ(2);
t99 = t69 * t84 + t83 * t85;
t64 = t85 * t69;
t93 = -t84 * t83 + t64;
t68 = pkin(5) * t80 + pkin(4);
t1 = (-m(4) * t106 - t79 * mrSges(3,1) - t81 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * pkin(6) - t112 * (pkin(3) * t72 + t106) + (t78 * mrSges(6,1) + t80 * mrSges(6,2) + qJ(4) * t112 - t117) * t74 + (-mrSges(6,3) + t116) * t72) * g(3) + (-mrSges(1,2) - m(4) * t99 - (t101 * t72 - t102) * mrSges(6,1) - (t100 * t72 + t103) * mrSges(6,2) - mrSges(6,3) * t105 - t112 * (pkin(3) * t105 + t84 * t98 + t99) + (m(6) * pkin(4) + m(7) * t68 - t109) * t85 + t108 * t84) * g(2) + (-mrSges(1,1) - m(4) * t93 - m(5) * (t93 + t113) - (t103 * t72 + t100) * mrSges(6,1) - (t102 * t72 - t101) * mrSges(6,2) - mrSges(6,3) * t104 + t114 * (t64 + t113) + t108 * t85 + (-m(6) * (pkin(4) - t83) - m(7) * (t68 - t83) + t109) * t84) * g(1);
U  = t1;

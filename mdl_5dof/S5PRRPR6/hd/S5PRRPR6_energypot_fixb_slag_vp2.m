% Calculate potential energy for
% S5PRRPR6
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:00
% EndTime: 2019-12-05 16:30:01
% DurationCPUTime: 0.60s
% Computational Cost: add. (219->85), mult. (442->103), div. (0->0), fcn. (511->12), ass. (0->39)
t115 = -m(4) - m(5);
t114 = -m(5) * qJ(4) + m(6) * (-pkin(8) - qJ(4)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t84 = sin(pkin(10));
t100 = pkin(4) * t84 + pkin(7);
t83 = pkin(10) + qJ(5);
t78 = sin(t83);
t79 = cos(t83);
t87 = cos(pkin(10));
t113 = -m(6) * t100 - t84 * mrSges(5,1) - t78 * mrSges(6,1) - t87 * mrSges(5,2) - t79 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3);
t77 = t87 * pkin(4) + pkin(3);
t112 = -m(5) * pkin(3) - m(6) * t77 - t87 * mrSges(5,1) - t79 * mrSges(6,1) + t84 * mrSges(5,2) + t78 * mrSges(6,2) - mrSges(4,1);
t85 = sin(pkin(9));
t86 = sin(pkin(5));
t111 = t85 * t86;
t92 = sin(qJ(2));
t110 = t86 * t92;
t93 = cos(qJ(3));
t109 = t86 * t93;
t94 = cos(qJ(2));
t108 = t86 * t94;
t88 = cos(pkin(9));
t107 = t88 * t86;
t89 = cos(pkin(5));
t106 = t89 * t92;
t105 = t89 * t94;
t104 = t88 * pkin(1) + pkin(6) * t111;
t103 = t89 * pkin(6) + qJ(1);
t69 = -t85 * t106 + t88 * t94;
t102 = t69 * pkin(2) + t104;
t101 = pkin(2) * t110 + t103;
t99 = t85 * pkin(1) - pkin(6) * t107;
t67 = t88 * t106 + t85 * t94;
t98 = t67 * pkin(2) + t99;
t96 = -pkin(7) * t108 + t101;
t91 = sin(qJ(3));
t71 = t92 * t109 + t89 * t91;
t68 = t85 * t105 + t88 * t92;
t66 = -t88 * t105 + t85 * t92;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t103 - t89 * mrSges(3,3) - (t92 * mrSges(3,1) + t94 * mrSges(3,2)) * t86 - m(4) * t96 - t71 * mrSges(4,1) + mrSges(4,3) * t108 - m(5) * (t71 * pkin(3) + t96) - (-t84 * t108 + t71 * t87) * mrSges(5,1) - (-t87 * t108 - t71 * t84) * mrSges(5,2) - m(6) * (-t100 * t108 + t71 * t77 + t101) - (-t78 * t108 + t71 * t79) * mrSges(6,1) - (-t79 * t108 - t71 * t78) * mrSges(6,2) + t114 * (t91 * t110 - t89 * t93)) * g(3) + (-m(3) * t99 - m(6) * t98 - t85 * mrSges(2,1) - t67 * mrSges(3,1) - t88 * mrSges(2,2) + mrSges(3,3) * t107 - mrSges(1,2) + t115 * (t66 * pkin(7) + t98) + t112 * (-t91 * t107 + t67 * t93) + t113 * t66 + t114 * (t93 * t107 + t67 * t91)) * g(2) + (-m(3) * t104 - m(6) * t102 - t88 * mrSges(2,1) - t69 * mrSges(3,1) + t85 * mrSges(2,2) - mrSges(3,3) * t111 - mrSges(1,1) + t115 * (t68 * pkin(7) + t102) + t112 * (t91 * t111 + t69 * t93) + t113 * t68 + t114 * (-t85 * t109 + t69 * t91)) * g(1);
U = t1;

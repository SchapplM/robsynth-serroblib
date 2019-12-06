% Calculate potential energy for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR9_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR9_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:18:57
% EndTime: 2019-12-05 17:18:58
% DurationCPUTime: 0.58s
% Computational Cost: add. (219->85), mult. (442->103), div. (0->0), fcn. (511->12), ass. (0->39)
t115 = -m(4) - m(5);
t114 = -m(5) * pkin(8) + m(6) * (-pkin(9) - pkin(8)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t88 = sin(qJ(4));
t100 = pkin(4) * t88 + pkin(7);
t83 = qJ(4) + qJ(5);
t78 = sin(t83);
t79 = cos(t83);
t91 = cos(qJ(4));
t113 = -m(6) * t100 - t88 * mrSges(5,1) - t78 * mrSges(6,1) - t91 * mrSges(5,2) - t79 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3);
t77 = t91 * pkin(4) + pkin(3);
t112 = -m(5) * pkin(3) - m(6) * t77 - t91 * mrSges(5,1) - t79 * mrSges(6,1) + t88 * mrSges(5,2) + t78 * mrSges(6,2) - mrSges(4,1);
t84 = sin(pkin(10));
t85 = sin(pkin(5));
t111 = t84 * t85;
t90 = sin(qJ(2));
t110 = t85 * t90;
t92 = cos(qJ(3));
t109 = t85 * t92;
t93 = cos(qJ(2));
t108 = t85 * t93;
t86 = cos(pkin(10));
t107 = t86 * t85;
t87 = cos(pkin(5));
t106 = t87 * t90;
t105 = t87 * t93;
t104 = t86 * pkin(1) + pkin(6) * t111;
t103 = t87 * pkin(6) + qJ(1);
t69 = -t84 * t106 + t86 * t93;
t102 = t69 * pkin(2) + t104;
t101 = pkin(2) * t110 + t103;
t99 = t84 * pkin(1) - pkin(6) * t107;
t67 = t86 * t106 + t84 * t93;
t98 = t67 * pkin(2) + t99;
t96 = -pkin(7) * t108 + t101;
t89 = sin(qJ(3));
t71 = t90 * t109 + t87 * t89;
t68 = t84 * t105 + t86 * t90;
t66 = -t86 * t105 + t84 * t90;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t103 - t87 * mrSges(3,3) - (t90 * mrSges(3,1) + t93 * mrSges(3,2)) * t85 - m(4) * t96 - t71 * mrSges(4,1) + mrSges(4,3) * t108 - m(5) * (t71 * pkin(3) + t96) - (-t88 * t108 + t71 * t91) * mrSges(5,1) - (-t91 * t108 - t71 * t88) * mrSges(5,2) - m(6) * (-t100 * t108 + t71 * t77 + t101) - (-t78 * t108 + t71 * t79) * mrSges(6,1) - (-t79 * t108 - t71 * t78) * mrSges(6,2) + t114 * (t89 * t110 - t87 * t92)) * g(3) + (-m(3) * t99 - m(6) * t98 - t84 * mrSges(2,1) - t67 * mrSges(3,1) - t86 * mrSges(2,2) + mrSges(3,3) * t107 - mrSges(1,2) + t115 * (t66 * pkin(7) + t98) + t112 * (-t89 * t107 + t67 * t92) + t113 * t66 + t114 * (t92 * t107 + t67 * t89)) * g(2) + (-m(3) * t104 - m(6) * t102 - t86 * mrSges(2,1) - t69 * mrSges(3,1) + t84 * mrSges(2,2) - mrSges(3,3) * t111 - mrSges(1,1) + t115 * (t68 * pkin(7) + t102) + t112 * (t89 * t111 + t69 * t92) + t113 * t68 + t114 * (-t84 * t109 + t69 * t89)) * g(1);
U = t1;

% Calculate potential energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:53
% EndTime: 2019-12-05 17:30:54
% DurationCPUTime: 0.60s
% Computational Cost: add. (158->79), mult. (320->92), div. (0->0), fcn. (348->10), ass. (0->35)
t114 = -m(5) - m(6);
t84 = sin(pkin(7));
t87 = cos(pkin(7));
t113 = -t87 * mrSges(3,1) + t84 * mrSges(3,2) - mrSges(2,1);
t112 = mrSges(2,2) - mrSges(3,3);
t111 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t88 = sin(qJ(5));
t90 = cos(qJ(5));
t110 = t88 * mrSges(6,1) + t90 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t109 = -m(6) * pkin(4) - t90 * mrSges(6,1) + t88 * mrSges(6,2) - mrSges(5,1);
t83 = sin(pkin(8));
t107 = t83 * t84;
t86 = cos(pkin(8));
t106 = t84 * t86;
t89 = sin(qJ(1));
t105 = t89 * t84;
t104 = t89 * t87;
t91 = cos(qJ(1));
t103 = t91 * t84;
t102 = t91 * t87;
t101 = t91 * pkin(1) + t89 * qJ(2);
t100 = qJ(3) * t84;
t99 = pkin(2) * t102 + t91 * t100 + t101;
t98 = t84 * pkin(2) - t87 * qJ(3) + pkin(5);
t95 = -pkin(2) * t87 - pkin(1) - t100;
t94 = pkin(3) * t106 + qJ(4) * t107 + t98;
t85 = cos(pkin(9));
t82 = sin(pkin(9));
t80 = t91 * qJ(2);
t71 = t102 * t86 + t89 * t83;
t70 = t102 * t83 - t89 * t86;
t69 = -t104 * t86 + t91 * t83;
t68 = t104 * t83 + t91 * t86;
t67 = t106 * t85 - t87 * t82;
t1 = (-m(3) * t101 - m(4) * t99 - t71 * mrSges(4,1) - mrSges(4,3) * t103 - mrSges(1,3) + t114 * (t71 * pkin(3) + t70 * qJ(4) + t99) + t113 * t91 + t112 * t89 + t109 * (t103 * t82 + t71 * t85) - t110 * t70 + t111 * (-t103 * t85 + t71 * t82)) * g(3) + (-t69 * mrSges(4,1) - mrSges(1,2) + t112 * t91 + (-m(3) - m(4)) * t80 + t114 * (t69 * pkin(3) - t68 * qJ(4) + t95 * t89 + t80) + t109 * (-t105 * t82 + t69 * t85) + t110 * t68 + t111 * (t105 * t85 + t69 * t82) + (m(3) * pkin(1) - m(4) * t95 + t84 * mrSges(4,3) - t113) * t89) * g(2) + (-mrSges(1,1) - mrSges(2,3) - m(4) * t98 - m(5) * t94 - t67 * mrSges(5,1) - mrSges(5,3) * t107 - m(6) * (t67 * pkin(4) + t94) - (t107 * t88 + t67 * t90) * mrSges(6,1) - (t107 * t90 - t67 * t88) * mrSges(6,2) + (-mrSges(3,2) + mrSges(4,3)) * t87 + (-t86 * mrSges(4,1) + t83 * mrSges(4,2) - mrSges(3,1)) * t84 + t111 * (t106 * t82 + t87 * t85) + (-m(2) - m(3)) * pkin(5)) * g(1);
U = t1;

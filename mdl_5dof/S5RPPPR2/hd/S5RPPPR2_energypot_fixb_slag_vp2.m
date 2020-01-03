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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:22:15
% EndTime: 2020-01-03 11:22:16
% DurationCPUTime: 0.67s
% Computational Cost: add. (158->78), mult. (320->91), div. (0->0), fcn. (348->10), ass. (0->36)
t116 = -m(5) - m(6);
t82 = sin(pkin(7));
t85 = cos(pkin(7));
t115 = -t85 * mrSges(3,1) + t82 * mrSges(3,2) - mrSges(2,1);
t114 = -mrSges(2,2) + mrSges(3,3);
t113 = pkin(2) * t85 + qJ(3) * t82;
t112 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t86 = sin(qJ(5));
t88 = cos(qJ(5));
t111 = t86 * mrSges(6,1) + t88 * mrSges(6,2) - mrSges(4,2) + mrSges(5,3);
t110 = -m(6) * pkin(4) - t88 * mrSges(6,1) + t86 * mrSges(6,2) - mrSges(5,1);
t81 = sin(pkin(8));
t107 = t81 * t82;
t84 = cos(pkin(8));
t106 = t82 * t84;
t87 = sin(qJ(1));
t105 = t82 * t87;
t89 = cos(qJ(1));
t104 = t82 * t89;
t103 = t85 * t89;
t102 = t87 * t81;
t101 = t87 * t84;
t99 = t87 * qJ(2);
t98 = t87 * pkin(1) - qJ(2) * t89;
t97 = t82 * pkin(2) - qJ(3) * t85 + pkin(5);
t95 = t113 * t87 + t98;
t94 = -pkin(1) - t113;
t92 = pkin(3) * t106 + qJ(4) * t107 + t97;
t83 = cos(pkin(9));
t80 = sin(pkin(9));
t71 = -t103 * t84 - t102;
t70 = t103 * t81 - t101;
t69 = t101 * t85 - t81 * t89;
t68 = t102 * t85 + t84 * t89;
t67 = t106 * t83 - t80 * t85;
t1 = (m(4) * t99 - t71 * mrSges(4,1) - mrSges(1,3) + (m(3) * qJ(2) + t114) * t87 + t116 * (t71 * pkin(3) - t70 * qJ(4) + t94 * t89 - t99) + t110 * (-t104 * t80 + t71 * t83) + t111 * t70 + t112 * (t104 * t83 + t71 * t80) + (m(3) * pkin(1) - m(4) * t94 + t82 * mrSges(4,3) - t115) * t89) * g(3) + (-m(3) * t98 - m(4) * t95 - t69 * mrSges(4,1) - mrSges(4,3) * t105 - mrSges(1,2) + t116 * (t69 * pkin(3) + t68 * qJ(4) + t95) + t114 * t89 + t115 * t87 + t110 * (t105 * t80 + t69 * t83) - t111 * t68 + t112 * (-t105 * t83 + t69 * t80)) * g(2) + (-mrSges(1,1) - mrSges(2,3) - m(4) * t97 - m(5) * t92 - t67 * mrSges(5,1) - mrSges(5,3) * t107 - m(6) * (pkin(4) * t67 + t92) - (t107 * t86 + t67 * t88) * mrSges(6,1) - (t107 * t88 - t67 * t86) * mrSges(6,2) + (-mrSges(3,2) + mrSges(4,3)) * t85 + (-t84 * mrSges(4,1) + t81 * mrSges(4,2) - mrSges(3,1)) * t82 + t112 * (t106 * t80 + t83 * t85) + (-m(2) - m(3)) * pkin(5)) * g(1);
U = t1;

% Calculate potential energy for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPPRR2_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:22
% EndTime: 2019-12-05 14:59:23
% DurationCPUTime: 0.43s
% Computational Cost: add. (158->77), mult. (320->92), div. (0->0), fcn. (348->10), ass. (0->35)
t117 = -m(5) - m(6);
t116 = mrSges(2,2) - mrSges(3,3);
t87 = sin(pkin(8));
t90 = cos(pkin(8));
t115 = -t90 * mrSges(3,1) + t87 * mrSges(3,2) - mrSges(2,1);
t114 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t92 = sin(qJ(5));
t94 = cos(qJ(5));
t113 = -t92 * mrSges(6,1) - t94 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3);
t112 = -m(6) * pkin(4) - t94 * mrSges(6,1) + t92 * mrSges(6,2) - mrSges(5,1);
t86 = sin(pkin(9));
t111 = t86 * t87;
t89 = cos(pkin(9));
t110 = t87 * t89;
t88 = sin(pkin(7));
t109 = t88 * t87;
t108 = t88 * t90;
t91 = cos(pkin(7));
t107 = t91 * t87;
t106 = t91 * t90;
t105 = t91 * pkin(1) + t88 * qJ(2);
t104 = qJ(3) * t87;
t103 = t88 * pkin(1) - t91 * qJ(2);
t102 = pkin(2) * t106 + t91 * t104 + t105;
t101 = t87 * pkin(2) - t90 * qJ(3) + qJ(1);
t99 = pkin(2) * t108 + t88 * t104 + t103;
t98 = pkin(3) * t110 + pkin(5) * t111 + t101;
t95 = cos(qJ(4));
t93 = sin(qJ(4));
t72 = t95 * t110 - t90 * t93;
t70 = t89 * t106 + t88 * t86;
t69 = t86 * t106 - t88 * t89;
t68 = t89 * t108 - t91 * t86;
t67 = t86 * t108 + t91 * t89;
t1 = (-mrSges(1,3) - mrSges(2,3) - m(4) * t101 - m(5) * t98 - t72 * mrSges(5,1) - mrSges(5,3) * t111 - m(6) * (t72 * pkin(4) + t98) - (t92 * t111 + t72 * t94) * mrSges(6,1) - (t94 * t111 - t72 * t92) * mrSges(6,2) + (-mrSges(3,2) + mrSges(4,3)) * t90 + (-t89 * mrSges(4,1) + t86 * mrSges(4,2) - mrSges(3,1)) * t87 + t114 * (t93 * t110 + t90 * t95) + (-m(2) - m(3)) * qJ(1)) * g(3) + (-m(3) * t103 - m(4) * t99 - t68 * mrSges(4,1) - mrSges(4,3) * t109 - mrSges(1,2) + t117 * (t68 * pkin(3) + t67 * pkin(5) + t99) - t116 * t91 + t115 * t88 + t112 * (t93 * t109 + t68 * t95) + t113 * t67 + t114 * (-t95 * t109 + t68 * t93)) * g(2) + (-m(3) * t105 - m(4) * t102 - t70 * mrSges(4,1) - mrSges(4,3) * t107 - mrSges(1,1) + t117 * (t70 * pkin(3) + t69 * pkin(5) + t102) + t115 * t91 + t116 * t88 + t112 * (t93 * t107 + t70 * t95) + t113 * t69 + t114 * (-t95 * t107 + t70 * t93)) * g(1);
U = t1;

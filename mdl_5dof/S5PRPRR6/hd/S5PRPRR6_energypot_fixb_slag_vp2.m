% Calculate potential energy for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRPRR6_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR6_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_energypot_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:22
% EndTime: 2019-12-05 15:56:22
% DurationCPUTime: 0.66s
% Computational Cost: add. (227->91), mult. (378->108), div. (0->0), fcn. (422->12), ass. (0->40)
t122 = -m(5) - m(6);
t121 = -m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3);
t120 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t101 = cos(qJ(5));
t99 = sin(qJ(5));
t119 = -m(6) * pkin(4) - t101 * mrSges(6,1) + t99 * mrSges(6,2) - mrSges(5,1);
t118 = -t99 * mrSges(6,1) - t101 * mrSges(6,2) - mrSges(5,3) - t120;
t92 = sin(pkin(10));
t117 = pkin(3) * t92;
t93 = sin(pkin(9));
t94 = sin(pkin(5));
t116 = t93 * t94;
t96 = cos(pkin(9));
t115 = t94 * t96;
t114 = t96 * pkin(1) + pkin(6) * t116;
t100 = sin(qJ(2));
t113 = t100 * t93;
t112 = t100 * t94;
t102 = cos(qJ(2));
t111 = t102 * t94;
t97 = cos(pkin(5));
t110 = t102 * t97;
t109 = t96 * t100;
t108 = t97 * pkin(6) + qJ(1);
t107 = t92 * t116;
t88 = t93 * pkin(1);
t106 = -pkin(6) * t115 + t88;
t95 = cos(pkin(10));
t85 = pkin(3) * t95 + pkin(2);
t98 = -pkin(7) - qJ(3);
t104 = t98 * t111 + t85 * t112 + t97 * t117 + t108;
t91 = pkin(10) + qJ(4);
t87 = cos(t91);
t86 = sin(t91);
t76 = t102 * t96 - t97 * t113;
t75 = t93 * t110 + t109;
t74 = t102 * t93 + t97 * t109;
t73 = -t96 * t110 + t113;
t70 = t87 * t112 + t86 * t97;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(5) * t104 - t70 * mrSges(5,1) + mrSges(5,3) * t111 - m(6) * (pkin(4) * t70 + t104) - (t70 * t101 - t99 * t111) * mrSges(6,1) - (-t101 * t111 - t70 * t99) * mrSges(6,2) + (-mrSges(4,1) * t92 - mrSges(4,2) * t95 - mrSges(3,3)) * t97 + (t120 * t102 + (-m(4) * pkin(2) - t95 * mrSges(4,1) + t92 * mrSges(4,2) - mrSges(3,1)) * t100) * t94 + t121 * (t86 * t112 - t97 * t87) + (-m(3) - m(4)) * t108) * g(3) + (-mrSges(1,2) - t93 * mrSges(2,1) - t96 * mrSges(2,2) - m(3) * t106 - t74 * mrSges(3,1) + mrSges(3,3) * t115 - m(4) * (pkin(2) * t74 + t106) - (-t92 * t115 + t74 * t95) * mrSges(4,1) - (-t95 * t115 - t74 * t92) * mrSges(4,2) + t122 * (t74 * t85 - t73 * t98 + t88 + (-pkin(6) - t117) * t115) + t121 * (t87 * t115 + t74 * t86) + t119 * (-t86 * t115 + t74 * t87) + t118 * t73) * g(2) + (-mrSges(1,1) - t96 * mrSges(2,1) + t93 * mrSges(2,2) - m(3) * t114 - t76 * mrSges(3,1) - mrSges(3,3) * t116 - m(4) * (pkin(2) * t76 + t114) - (t76 * t95 + t107) * mrSges(4,1) - (t95 * t116 - t76 * t92) * mrSges(4,2) + t122 * (pkin(3) * t107 - t75 * t98 + t76 * t85 + t114) + t121 * (-t87 * t116 + t76 * t86) + t119 * (t86 * t116 + t76 * t87) + t118 * t75) * g(1);
U = t1;

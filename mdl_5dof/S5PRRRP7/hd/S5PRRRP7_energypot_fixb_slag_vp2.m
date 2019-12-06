% Calculate potential energy for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRP7_energypot_fixb_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_energypot_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_energypot_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_energypot_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_energypot_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP7_energypot_fixb_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:02
% EndTime: 2019-12-05 16:54:02
% DurationCPUTime: 0.48s
% Computational Cost: add. (207->85), mult. (442->106), div. (0->0), fcn. (511->10), ass. (0->41)
t122 = -mrSges(5,1) - mrSges(6,1);
t121 = -mrSges(5,2) - mrSges(6,2);
t95 = sin(qJ(4));
t106 = pkin(4) * t95 + pkin(7);
t120 = -m(6) * t106 + mrSges(3,2) - mrSges(4,3);
t119 = -m(5) * pkin(8) + m(6) * (-qJ(5) - pkin(8)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t90 = sin(pkin(9));
t91 = sin(pkin(5));
t118 = t90 * t91;
t92 = cos(pkin(9));
t117 = t91 * t92;
t96 = sin(qJ(3));
t116 = t91 * t96;
t97 = sin(qJ(2));
t115 = t91 * t97;
t99 = cos(qJ(3));
t114 = t91 * t99;
t93 = cos(pkin(5));
t113 = t93 * t97;
t112 = t92 * pkin(1) + pkin(6) * t118;
t100 = cos(qJ(2));
t111 = t100 * t91;
t110 = t100 * t93;
t109 = t93 * pkin(6) + qJ(1);
t78 = t100 * t92 - t113 * t90;
t108 = t78 * pkin(2) + t112;
t107 = pkin(2) * t115 + t109;
t105 = t90 * pkin(1) - pkin(6) * t117;
t76 = t100 * t90 + t113 * t92;
t104 = t76 * pkin(2) + t105;
t77 = t110 * t90 + t92 * t97;
t103 = pkin(7) * t77 + t108;
t102 = -pkin(7) * t111 + t107;
t75 = -t110 * t92 + t90 * t97;
t101 = pkin(7) * t75 + t104;
t98 = cos(qJ(4));
t86 = pkin(4) * t98 + pkin(3);
t80 = t114 * t97 + t93 * t96;
t70 = t116 * t90 + t78 * t99;
t68 = -t116 * t92 + t76 * t99;
t1 = (-mrSges(1,3) - m(2) * qJ(1) - mrSges(2,3) - m(3) * t109 - t93 * mrSges(3,3) - (t97 * mrSges(3,1) + t100 * mrSges(3,2)) * t91 - m(4) * t102 - t80 * mrSges(4,1) + mrSges(4,3) * t111 - m(5) * (t80 * pkin(3) + t102) - m(6) * (-t106 * t111 + t80 * t86 + t107) + t122 * (-t111 * t95 + t80 * t98) + t121 * (-t111 * t98 - t80 * t95) + t119 * (t115 * t96 - t93 * t99)) * g(3) + (-mrSges(1,2) - t90 * mrSges(2,1) - t92 * mrSges(2,2) - m(3) * t105 - t76 * mrSges(3,1) + mrSges(3,3) * t117 - m(4) * t101 - t68 * mrSges(4,1) - m(5) * (pkin(3) * t68 + t101) - m(6) * (t68 * t86 + t104) + t120 * t75 + t122 * (t68 * t98 + t75 * t95) + t121 * (-t68 * t95 + t75 * t98) + t119 * (t114 * t92 + t76 * t96)) * g(2) + (-mrSges(1,1) - t92 * mrSges(2,1) + t90 * mrSges(2,2) - m(3) * t112 - t78 * mrSges(3,1) - mrSges(3,3) * t118 - m(4) * t103 - t70 * mrSges(4,1) - m(5) * (pkin(3) * t70 + t103) - m(6) * (t70 * t86 + t108) + t120 * t77 + t122 * (t70 * t98 + t77 * t95) + t121 * (-t70 * t95 + t77 * t98) + t119 * (-t114 * t90 + t78 * t96)) * g(1);
U = t1;

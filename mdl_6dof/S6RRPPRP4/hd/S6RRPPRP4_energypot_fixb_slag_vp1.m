% Calculate potential energy for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRP4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:36:56
% EndTime: 2019-03-09 08:36:56
% DurationCPUTime: 0.38s
% Computational Cost: add. (203->103), mult. (392->127), div. (0->0), fcn. (430->8), ass. (0->40)
t128 = rSges(7,1) + pkin(5);
t127 = rSges(7,3) + qJ(6);
t100 = sin(qJ(2));
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t126 = (g(1) * t104 + g(2) * t101) * t100;
t123 = t100 * pkin(2) + pkin(6);
t120 = t104 * pkin(1) + t101 * pkin(7);
t97 = sin(pkin(9));
t119 = t100 * t97;
t98 = cos(pkin(9));
t118 = t100 * t98;
t117 = rSges(5,3) + qJ(4);
t116 = t100 * t101;
t115 = t100 * t104;
t103 = cos(qJ(2));
t114 = t101 * t103;
t113 = t103 * t104;
t112 = pkin(3) * t118 + qJ(4) * t119 + t123;
t111 = pkin(2) * t113 + qJ(3) * t115 + t120;
t82 = t101 * t97 + t98 * t113;
t110 = t82 * pkin(3) + t111;
t94 = t101 * pkin(1);
t109 = pkin(2) * t114 - t104 * pkin(7) + qJ(3) * t116 + t94;
t108 = pkin(4) * t118 + t103 * pkin(8) + t112;
t80 = -t104 * t97 + t98 * t114;
t107 = t80 * pkin(3) + t109;
t81 = -t101 * t98 + t97 * t113;
t106 = t82 * pkin(4) + t81 * qJ(4) + t110;
t79 = t104 * t98 + t97 * t114;
t105 = t80 * pkin(4) + t79 * qJ(4) + t107;
t102 = cos(qJ(5));
t99 = sin(qJ(5));
t74 = (t102 * t98 + t97 * t99) * t100;
t73 = -t102 * t119 + t99 * t118;
t70 = t102 * t82 + t81 * t99;
t69 = -t81 * t102 + t82 * t99;
t68 = t102 * t80 + t79 * t99;
t67 = -t79 * t102 + t80 * t99;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t104 * rSges(2,1) - t101 * rSges(2,2)) + g(2) * (t101 * rSges(2,1) + t104 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t101 * rSges(3,3) + t120) + g(2) * (rSges(3,1) * t114 - rSges(3,2) * t116 + t94) + g(3) * (t100 * rSges(3,1) + t103 * rSges(3,2) + pkin(6)) + (g(1) * (rSges(3,1) * t103 - rSges(3,2) * t100) + g(2) * (-rSges(3,3) - pkin(7))) * t104) - m(4) * (g(1) * (t82 * rSges(4,1) - t81 * rSges(4,2) + rSges(4,3) * t115 + t111) + g(2) * (t80 * rSges(4,1) - t79 * rSges(4,2) + rSges(4,3) * t116 + t109) + g(3) * ((-rSges(4,3) - qJ(3)) * t103 + (rSges(4,1) * t98 - rSges(4,2) * t97) * t100 + t123)) - m(5) * (g(1) * (t82 * rSges(5,1) + rSges(5,2) * t115 + t117 * t81 + t110) + g(2) * (t80 * rSges(5,1) + rSges(5,2) * t116 + t117 * t79 + t107) + g(3) * ((-rSges(5,2) - qJ(3)) * t103 + (rSges(5,1) * t98 + rSges(5,3) * t97) * t100 + t112)) - m(6) * (g(1) * (t70 * rSges(6,1) - t69 * rSges(6,2) + t106) + g(2) * (t68 * rSges(6,1) - t67 * rSges(6,2) + t105) + g(3) * (t74 * rSges(6,1) - t73 * rSges(6,2) + (rSges(6,3) - qJ(3)) * t103 + t108) + (-rSges(6,3) - pkin(8)) * t126) - m(7) * (g(1) * (t127 * t69 + t128 * t70 + t106) + g(2) * (t127 * t67 + t128 * t68 + t105) + g(3) * (t128 * t74 + t127 * t73 + (rSges(7,2) - qJ(3)) * t103 + t108) + (-rSges(7,2) - pkin(8)) * t126);
U  = t1;

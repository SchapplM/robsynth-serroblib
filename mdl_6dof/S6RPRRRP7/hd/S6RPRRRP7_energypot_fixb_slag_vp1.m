% Calculate potential energy for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRP7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP7_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:33
% EndTime: 2019-03-09 06:18:33
% DurationCPUTime: 0.44s
% Computational Cost: add. (239->95), mult. (227->115), div. (0->0), fcn. (215->10), ass. (0->45)
t118 = rSges(7,1) + pkin(5);
t117 = rSges(5,3) + pkin(8);
t116 = rSges(7,3) + qJ(6);
t87 = sin(qJ(1));
t89 = cos(qJ(1));
t115 = g(1) * t89 + g(2) * t87;
t83 = sin(pkin(10));
t111 = t83 * pkin(2) + pkin(6);
t81 = pkin(10) + qJ(3);
t76 = sin(t81);
t110 = rSges(4,2) * t76;
t77 = cos(t81);
t109 = t87 * t77;
t82 = qJ(4) + qJ(5);
t78 = sin(t82);
t108 = t87 * t78;
t79 = cos(t82);
t107 = t87 * t79;
t86 = sin(qJ(4));
t106 = t87 * t86;
t88 = cos(qJ(4));
t105 = t87 * t88;
t104 = t89 * t77;
t103 = t89 * t78;
t102 = t89 * t79;
t101 = t89 * t86;
t100 = t89 * t88;
t84 = cos(pkin(10));
t73 = t84 * pkin(2) + pkin(1);
t85 = -pkin(7) - qJ(2);
t97 = t87 * t73 + t89 * t85;
t96 = rSges(3,3) + qJ(2);
t69 = t89 * t73;
t95 = -t87 * t85 + t69;
t75 = t88 * pkin(4) + pkin(3);
t90 = -pkin(9) - pkin(8);
t94 = t76 * t75 + t77 * t90 + t111;
t93 = pkin(4) * t106 + t75 * t104 + t95;
t92 = rSges(3,1) * t84 - rSges(3,2) * t83 + pkin(1);
t91 = -pkin(4) * t101 + t75 * t109 + t97;
t64 = t77 * t102 + t108;
t63 = t77 * t103 - t107;
t62 = t77 * t107 - t103;
t61 = t77 * t108 + t102;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t89 * rSges(2,1) - t87 * rSges(2,2)) + g(2) * (t87 * rSges(2,1) + t89 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t83 * rSges(3,1) + t84 * rSges(3,2) + pkin(6)) + (g(1) * t92 - g(2) * t96) * t89 + (g(1) * t96 + g(2) * t92) * t87) - m(4) * (g(1) * (rSges(4,1) * t104 - t89 * t110 + t69) + g(2) * (-t89 * rSges(4,3) + t97) + g(3) * (t76 * rSges(4,1) + t77 * rSges(4,2) + t111) + (g(1) * (rSges(4,3) - t85) + g(2) * (rSges(4,1) * t77 - t110)) * t87) - m(5) * (g(1) * (pkin(3) * t104 + (t77 * t100 + t106) * rSges(5,1) + (-t77 * t101 + t105) * rSges(5,2) + t95) + g(2) * (pkin(3) * t109 + (t77 * t105 - t101) * rSges(5,1) + (-t77 * t106 - t100) * rSges(5,2) + t97) + g(3) * (-t117 * t77 + t111) + (g(3) * (rSges(5,1) * t88 - rSges(5,2) * t86 + pkin(3)) + t115 * t117) * t76) - m(6) * (g(1) * (t64 * rSges(6,1) - t63 * rSges(6,2) + t93) + g(2) * (t62 * rSges(6,1) - t61 * rSges(6,2) + t91) + g(3) * (-t77 * rSges(6,3) + t94) + (g(3) * (rSges(6,1) * t79 - rSges(6,2) * t78) + t115 * (rSges(6,3) - t90)) * t76) - m(7) * (g(1) * (t116 * t63 + t118 * t64 + t93) + g(2) * (t116 * t61 + t118 * t62 + t91) + g(3) * (-t77 * rSges(7,2) + t94) + (g(3) * (t116 * t78 + t118 * t79) + t115 * (rSges(7,2) - t90)) * t76);
U  = t1;

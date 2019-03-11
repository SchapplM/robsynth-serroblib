% Calculate potential energy for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:39:02
% EndTime: 2019-03-09 16:39:02
% DurationCPUTime: 0.44s
% Computational Cost: add. (239->95), mult. (227->115), div. (0->0), fcn. (215->10), ass. (0->41)
t114 = rSges(7,1) + pkin(5);
t113 = rSges(5,3) + qJ(4);
t112 = rSges(7,3) + qJ(6);
t87 = sin(qJ(1));
t89 = cos(qJ(1));
t111 = g(1) * t89 + g(2) * t87;
t108 = rSges(3,3) + pkin(7);
t86 = sin(qJ(2));
t107 = t86 * pkin(2) + pkin(6);
t82 = qJ(2) + qJ(3);
t78 = sin(t82);
t106 = rSges(4,2) * t78;
t79 = cos(t82);
t105 = t87 * t79;
t83 = sin(pkin(10));
t104 = t87 * t83;
t84 = cos(pkin(10));
t103 = t87 * t84;
t102 = t89 * t79;
t101 = t89 * t83;
t100 = t89 * t84;
t88 = cos(qJ(2));
t74 = t88 * pkin(2) + pkin(1);
t90 = -pkin(8) - pkin(7);
t97 = t87 * t74 + t89 * t90;
t69 = t89 * t74;
t95 = -t87 * t90 + t69;
t73 = t84 * pkin(4) + pkin(3);
t85 = -pkin(9) - qJ(4);
t94 = t78 * t73 + t79 * t85 + t107;
t93 = pkin(4) * t104 + t73 * t102 + t95;
t92 = rSges(3,1) * t88 - rSges(3,2) * t86 + pkin(1);
t91 = -pkin(4) * t101 + t73 * t105 + t97;
t81 = pkin(10) + qJ(5);
t77 = cos(t81);
t76 = sin(t81);
t64 = t77 * t102 + t87 * t76;
t63 = t76 * t102 - t87 * t77;
t62 = t77 * t105 - t89 * t76;
t61 = t76 * t105 + t89 * t77;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t89 * rSges(2,1) - t87 * rSges(2,2)) + g(2) * (t87 * rSges(2,1) + t89 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t86 * rSges(3,1) + t88 * rSges(3,2) + pkin(6)) + (g(1) * t92 - g(2) * t108) * t89 + (g(1) * t108 + g(2) * t92) * t87) - m(4) * (g(1) * (rSges(4,1) * t102 - t89 * t106 + t69) + g(2) * (-t89 * rSges(4,3) + t97) + g(3) * (t78 * rSges(4,1) + t79 * rSges(4,2) + t107) + (g(1) * (rSges(4,3) - t90) + g(2) * (rSges(4,1) * t79 - t106)) * t87) - m(5) * (g(1) * (pkin(3) * t102 + (t79 * t100 + t104) * rSges(5,1) + (-t79 * t101 + t103) * rSges(5,2) + t95) + g(2) * (pkin(3) * t105 + (t79 * t103 - t101) * rSges(5,1) + (-t79 * t104 - t100) * rSges(5,2) + t97) + g(3) * (-t113 * t79 + t107) + (g(3) * (rSges(5,1) * t84 - rSges(5,2) * t83 + pkin(3)) + t111 * t113) * t78) - m(6) * (g(1) * (t64 * rSges(6,1) - t63 * rSges(6,2) + t93) + g(2) * (t62 * rSges(6,1) - t61 * rSges(6,2) + t91) + g(3) * (-t79 * rSges(6,3) + t94) + (g(3) * (rSges(6,1) * t77 - rSges(6,2) * t76) + t111 * (rSges(6,3) - t85)) * t78) - m(7) * (g(1) * (t112 * t63 + t114 * t64 + t93) + g(2) * (t112 * t61 + t114 * t62 + t91) + g(3) * (-t79 * rSges(7,2) + t94) + (g(3) * (t112 * t76 + t114 * t77) + t111 * (rSges(7,2) - t85)) * t78);
U  = t1;

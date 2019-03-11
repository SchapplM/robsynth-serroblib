% Calculate potential energy for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRPR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:08:58
% EndTime: 2019-03-09 05:08:58
% DurationCPUTime: 0.43s
% Computational Cost: add. (235->93), mult. (186->113), div. (0->0), fcn. (166->12), ass. (0->44)
t109 = rSges(7,3) + pkin(9) + qJ(5);
t108 = rSges(6,3) + qJ(5);
t82 = sin(qJ(1));
t83 = cos(qJ(1));
t107 = g(1) * t83 + g(2) * t82;
t77 = sin(pkin(10));
t104 = t77 * pkin(2) + pkin(6);
t79 = cos(pkin(10));
t65 = t79 * pkin(2) + pkin(1);
t75 = pkin(10) + qJ(3);
t70 = qJ(4) + t75;
t62 = sin(t70);
t103 = rSges(5,2) * t62;
t63 = cos(t70);
t102 = t63 * t82;
t101 = t63 * t83;
t74 = pkin(11) + qJ(6);
t66 = sin(t74);
t100 = t66 * t83;
t68 = cos(t74);
t99 = t68 * t83;
t76 = sin(pkin(11));
t98 = t76 * t83;
t78 = cos(pkin(11));
t97 = t78 * t83;
t96 = t82 * t66;
t95 = t82 * t68;
t94 = t82 * t76;
t93 = t82 * t78;
t81 = -pkin(7) - qJ(2);
t92 = rSges(4,3) - t81;
t69 = cos(t75);
t59 = pkin(3) * t69 + t65;
t73 = -pkin(8) + t81;
t90 = t82 * t59 + t83 * t73;
t89 = rSges(3,3) + qJ(2);
t67 = sin(t75);
t87 = pkin(3) * t67 + t104;
t58 = t83 * t59;
t86 = -t82 * t73 + t58;
t85 = rSges(3,1) * t79 - rSges(3,2) * t77 + pkin(1);
t84 = rSges(4,1) * t69 - rSges(4,2) * t67 + t65;
t64 = pkin(5) * t78 + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t83 - t82 * rSges(2,2)) + g(2) * (t82 * rSges(2,1) + rSges(2,2) * t83) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t77 + rSges(3,2) * t79 + pkin(6)) + (g(1) * t85 - g(2) * t89) * t83 + (g(1) * t89 + g(2) * t85) * t82) - m(4) * (g(3) * (rSges(4,1) * t67 + rSges(4,2) * t69 + t104) + (g(1) * t84 - g(2) * t92) * t83 + (g(1) * t92 + g(2) * t84) * t82) - m(5) * (g(1) * (rSges(5,1) * t101 - t83 * t103 + t58) + g(2) * (-rSges(5,3) * t83 + t90) + g(3) * (rSges(5,1) * t62 + rSges(5,2) * t63 + t87) + (g(1) * (rSges(5,3) - t73) + g(2) * (rSges(5,1) * t63 - t103)) * t82) - m(6) * (g(1) * (pkin(4) * t101 + (t63 * t97 + t94) * rSges(6,1) + (-t63 * t98 + t93) * rSges(6,2) + t86) + g(2) * (pkin(4) * t102 + (t63 * t93 - t98) * rSges(6,1) + (-t63 * t94 - t97) * rSges(6,2) + t90) + g(3) * (-t108 * t63 + t87) + (g(3) * (rSges(6,1) * t78 - rSges(6,2) * t76 + pkin(4)) + t107 * t108) * t62) - m(7) * (g(1) * (t64 * t101 + pkin(5) * t94 + (t63 * t99 + t96) * rSges(7,1) + (-t63 * t100 + t95) * rSges(7,2) + t86) + g(2) * (t64 * t102 - pkin(5) * t98 + (t63 * t95 - t100) * rSges(7,1) + (-t63 * t96 - t99) * rSges(7,2) + t90) + g(3) * (-t109 * t63 + t87) + (g(3) * (rSges(7,1) * t68 - rSges(7,2) * t66 + t64) + t107 * t109) * t62);
U  = t1;

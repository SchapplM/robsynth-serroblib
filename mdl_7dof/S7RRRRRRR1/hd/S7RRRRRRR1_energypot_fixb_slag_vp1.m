% Calculate potential energy for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% rSges [8x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S7RRRRRRR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1),zeros(8,1),zeros(8,3)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_energypot_fixb_slag_vp1: qJ has to be [7x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_energypot_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_energypot_fixb_slag_vp1: m has to be [8x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [8,3]), ...
  'S7RRRRRRR1_energypot_fixb_slag_vp1: rSges has to be [8x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 06:29:07
% EndTime: 2019-03-10 06:29:07
% DurationCPUTime: 0.39s
% Computational Cost: add. (294->109), mult. (689->160), div. (0->0), fcn. (857->14), ass. (0->53)
t112 = rSges(6,3) + pkin(3);
t94 = cos(qJ(2));
t111 = t94 * pkin(2) + pkin(1);
t110 = pkin(4) + rSges(8,3);
t109 = cos(qJ(4));
t87 = sin(qJ(3));
t88 = sin(qJ(2));
t108 = t87 * t88;
t89 = sin(qJ(1));
t107 = t89 * t88;
t106 = t89 * t94;
t95 = cos(qJ(1));
t105 = t95 * t87;
t104 = t95 * t88;
t93 = cos(qJ(3));
t103 = t95 * t93;
t102 = pkin(2) * t107;
t101 = pkin(2) * t104;
t100 = t88 * t109;
t86 = sin(qJ(4));
t74 = t88 * t93 * t86 + t94 * t109;
t99 = t74 * pkin(3) + t111;
t98 = rSges(3,1) * t94 - rSges(3,2) * t88;
t77 = t93 * t106 + t105;
t70 = -t89 * t100 + t77 * t86;
t97 = t70 * pkin(3) - t102;
t79 = t94 * t103 - t89 * t87;
t72 = -t95 * t100 + t79 * t86;
t96 = t72 * pkin(3) - t101;
t92 = cos(qJ(5));
t91 = cos(qJ(6));
t90 = cos(qJ(7));
t85 = sin(qJ(5));
t84 = sin(qJ(6));
t83 = sin(qJ(7));
t78 = -t94 * t105 - t89 * t93;
t76 = -t87 * t106 + t103;
t75 = t93 * t100 - t94 * t86;
t73 = t86 * t104 + t79 * t109;
t71 = t86 * t107 + t77 * t109;
t69 = -t85 * t108 + t75 * t92;
t68 = t92 * t108 + t75 * t85;
t67 = t73 * t92 + t78 * t85;
t66 = t73 * t85 - t78 * t92;
t65 = t71 * t92 + t76 * t85;
t64 = t71 * t85 - t76 * t92;
t63 = t69 * t91 + t74 * t84;
t62 = -t69 * t84 + t74 * t91;
t61 = t67 * t91 + t72 * t84;
t60 = -t67 * t84 + t72 * t91;
t59 = t65 * t91 + t70 * t84;
t58 = -t65 * t84 + t70 * t91;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t95 * rSges(2,1) - t89 * rSges(2,2)) + g(2) * (t89 * rSges(2,1) + t95 * rSges(2,2)) + g(3) * (pkin(1) + rSges(2,3))) - m(3) * (g(1) * (t89 * rSges(3,3) + t98 * t95) + g(2) * (-t95 * rSges(3,3) + t98 * t89) + g(3) * (t88 * rSges(3,1) + t94 * rSges(3,2) + pkin(1))) - m(4) * (g(1) * (t79 * rSges(4,1) + t78 * rSges(4,2)) + g(2) * (t77 * rSges(4,1) + t76 * rSges(4,2)) + g(3) * (t94 * rSges(4,3) + t111) + (g(3) * (rSges(4,1) * t93 - rSges(4,2) * t87) + (g(1) * t95 + g(2) * t89) * (-rSges(4,3) - pkin(2))) * t88) - m(5) * (g(1) * (t73 * rSges(5,1) - t72 * rSges(5,2) + t78 * rSges(5,3) - t101) + g(2) * (t71 * rSges(5,1) - t70 * rSges(5,2) + t76 * rSges(5,3) - t102) + g(3) * (t75 * rSges(5,1) - t74 * rSges(5,2) - rSges(5,3) * t108 + t111)) - m(6) * (g(1) * (t67 * rSges(6,1) - t66 * rSges(6,2) + t112 * t72 - t101) + g(2) * (t65 * rSges(6,1) - t64 * rSges(6,2) + t112 * t70 - t102) + g(3) * (t69 * rSges(6,1) - t68 * rSges(6,2) + t112 * t74 + t111)) - m(7) * (g(1) * (t61 * rSges(7,1) + t60 * rSges(7,2) + t66 * rSges(7,3) + t96) + g(2) * (t59 * rSges(7,1) + t58 * rSges(7,2) + t64 * rSges(7,3) + t97) + g(3) * (t63 * rSges(7,1) + t62 * rSges(7,2) + t68 * rSges(7,3) + t99)) - m(8) * (g(1) * ((t61 * t90 - t66 * t83) * rSges(8,1) + (-t61 * t83 - t66 * t90) * rSges(8,2) + t110 * t60 + t96) + g(2) * ((t59 * t90 - t64 * t83) * rSges(8,1) + (-t59 * t83 - t64 * t90) * rSges(8,2) + t110 * t58 + t97) + g(3) * ((t63 * t90 - t68 * t83) * rSges(8,1) + (-t63 * t83 - t68 * t90) * rSges(8,2) + t110 * t62 + t99));
U  = t1;

% Calculate potential energy for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP4_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:45
% EndTime: 2019-03-09 09:58:45
% DurationCPUTime: 0.51s
% Computational Cost: add. (183->101), mult. (252->122), div. (0->0), fcn. (240->8), ass. (0->40)
t112 = rSges(7,1) + pkin(5);
t111 = rSges(5,3) + pkin(8);
t110 = rSges(7,3) + qJ(6);
t81 = sin(qJ(1));
t84 = cos(qJ(1));
t109 = g(1) * t84 + g(2) * t81;
t80 = sin(qJ(2));
t105 = t80 * pkin(2) + pkin(6);
t104 = t80 * t81;
t103 = t80 * t84;
t77 = qJ(4) + pkin(9);
t72 = cos(t77);
t102 = t81 * t72;
t79 = sin(qJ(4));
t101 = t81 * t79;
t82 = cos(qJ(4));
t100 = t81 * t82;
t83 = cos(qJ(2));
t99 = t81 * t83;
t98 = t82 * t84;
t95 = t84 * pkin(1) + t81 * pkin(7);
t94 = qJ(3) * t80;
t93 = t79 * t103;
t92 = t80 * t101;
t75 = t81 * pkin(1);
t91 = pkin(2) * t99 + t81 * t94 + t75;
t90 = -pkin(4) * t79 - qJ(3);
t89 = t95 + (pkin(2) * t83 + t94) * t84;
t78 = -qJ(5) - pkin(8);
t88 = -t78 * t80 + t105;
t87 = -t84 * pkin(7) + t91;
t70 = pkin(4) * t82 + pkin(3);
t86 = pkin(4) * t93 + t81 * t70 + t89;
t85 = pkin(4) * t92 - t70 * t84 + t87;
t71 = sin(t77);
t61 = t104 * t71 - t72 * t84;
t60 = t102 * t80 + t71 * t84;
t59 = t103 * t71 + t102;
t58 = -t103 * t72 + t71 * t81;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t84 - t81 * rSges(2,2)) + g(2) * (t81 * rSges(2,1) + rSges(2,2) * t84) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t81 * rSges(3,3) + t95) + g(2) * (rSges(3,1) * t99 - rSges(3,2) * t104 + t75) + g(3) * (rSges(3,1) * t80 + rSges(3,2) * t83 + pkin(6)) + (g(1) * (rSges(3,1) * t83 - rSges(3,2) * t80) + g(2) * (-rSges(3,3) - pkin(7))) * t84) - m(4) * (g(1) * (t81 * rSges(4,1) + t89) + g(2) * (-rSges(4,2) * t99 + rSges(4,3) * t104 + t91) + g(3) * (-rSges(4,2) * t80 + (-rSges(4,3) - qJ(3)) * t83 + t105) + (g(1) * (-rSges(4,2) * t83 + rSges(4,3) * t80) + g(2) * (-rSges(4,1) - pkin(7))) * t84) - m(5) * (g(1) * (t81 * pkin(3) + (t93 + t100) * rSges(5,1) + (t80 * t98 - t101) * rSges(5,2) + t89) + g(2) * (-t84 * pkin(3) + (t92 - t98) * rSges(5,1) + (t80 * t100 + t79 * t84) * rSges(5,2) + t87) + g(3) * (t111 * t80 + t105) + (g(3) * (-rSges(5,1) * t79 - rSges(5,2) * t82 - qJ(3)) + t109 * t111) * t83) - m(6) * (g(1) * (t59 * rSges(6,1) - t58 * rSges(6,2) + t86) + g(2) * (t61 * rSges(6,1) + t60 * rSges(6,2) + t85) + g(3) * (rSges(6,3) * t80 + t88) + (g(3) * (-rSges(6,1) * t71 - rSges(6,2) * t72 + t90) + t109 * (rSges(6,3) - t78)) * t83) - m(7) * (g(1) * (t110 * t58 + t112 * t59 + t86) + g(2) * (-t110 * t60 + t112 * t61 + t85) + g(3) * (rSges(7,2) * t80 + t88) + (g(3) * (t110 * t72 - t112 * t71 + t90) + t109 * (rSges(7,2) - t78)) * t83);
U  = t1;

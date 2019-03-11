% Calculate potential energy for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:45:03
% EndTime: 2019-03-09 09:45:03
% DurationCPUTime: 0.44s
% Computational Cost: add. (239->95), mult. (227->115), div. (0->0), fcn. (215->10), ass. (0->41)
t110 = rSges(7,1) + pkin(5);
t109 = rSges(5,3) + pkin(8);
t108 = rSges(7,3) + qJ(6);
t83 = sin(qJ(1));
t86 = cos(qJ(1));
t107 = g(1) * t86 + g(2) * t83;
t104 = rSges(3,3) + pkin(7);
t82 = sin(qJ(2));
t102 = t82 * pkin(2) + pkin(6);
t78 = qJ(2) + pkin(9);
t73 = sin(t78);
t101 = rSges(4,2) * t73;
t75 = cos(t78);
t100 = t83 * t75;
t81 = sin(qJ(4));
t99 = t83 * t81;
t84 = cos(qJ(4));
t98 = t83 * t84;
t97 = t86 * t75;
t96 = t86 * t81;
t95 = t86 * t84;
t85 = cos(qJ(2));
t71 = t85 * pkin(2) + pkin(1);
t80 = -qJ(3) - pkin(7);
t92 = t83 * t71 + t86 * t80;
t65 = t86 * t71;
t91 = -t83 * t80 + t65;
t70 = t84 * pkin(4) + pkin(3);
t79 = -qJ(5) - pkin(8);
t90 = t73 * t70 + t75 * t79 + t102;
t89 = pkin(4) * t99 + t70 * t97 + t91;
t88 = rSges(3,1) * t85 - rSges(3,2) * t82 + pkin(1);
t87 = -pkin(4) * t96 + t70 * t100 + t92;
t77 = qJ(4) + pkin(10);
t74 = cos(t77);
t72 = sin(t77);
t60 = t83 * t72 + t74 * t97;
t59 = t72 * t97 - t83 * t74;
t58 = t74 * t100 - t86 * t72;
t57 = t72 * t100 + t86 * t74;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t86 * rSges(2,1) - t83 * rSges(2,2)) + g(2) * (t83 * rSges(2,1) + t86 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (t82 * rSges(3,1) + t85 * rSges(3,2) + pkin(6)) + (g(1) * t88 - g(2) * t104) * t86 + (g(1) * t104 + g(2) * t88) * t83) - m(4) * (g(1) * (rSges(4,1) * t97 - t86 * t101 + t65) + g(2) * (-t86 * rSges(4,3) + t92) + g(3) * (t73 * rSges(4,1) + t75 * rSges(4,2) + t102) + (g(1) * (rSges(4,3) - t80) + g(2) * (rSges(4,1) * t75 - t101)) * t83) - m(5) * (g(1) * (pkin(3) * t97 + (t75 * t95 + t99) * rSges(5,1) + (-t75 * t96 + t98) * rSges(5,2) + t91) + g(2) * (pkin(3) * t100 + (t75 * t98 - t96) * rSges(5,1) + (-t75 * t99 - t95) * rSges(5,2) + t92) + g(3) * (-t109 * t75 + t102) + (g(3) * (rSges(5,1) * t84 - rSges(5,2) * t81 + pkin(3)) + t107 * t109) * t73) - m(6) * (g(1) * (t60 * rSges(6,1) - t59 * rSges(6,2) + t89) + g(2) * (t58 * rSges(6,1) - t57 * rSges(6,2) + t87) + g(3) * (-t75 * rSges(6,3) + t90) + (g(3) * (rSges(6,1) * t74 - rSges(6,2) * t72) + t107 * (rSges(6,3) - t79)) * t73) - m(7) * (g(1) * (t108 * t59 + t110 * t60 + t89) + g(2) * (t108 * t57 + t110 * t58 + t87) + g(3) * (-t75 * rSges(7,2) + t90) + (g(3) * (t108 * t72 + t110 * t74) + t107 * (rSges(7,2) - t79)) * t73);
U  = t1;

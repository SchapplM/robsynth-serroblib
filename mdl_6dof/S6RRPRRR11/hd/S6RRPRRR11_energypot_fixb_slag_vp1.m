% Calculate potential energy for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR11_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:29:07
% EndTime: 2019-03-09 14:29:07
% DurationCPUTime: 0.55s
% Computational Cost: add. (186->110), mult. (237->135), div. (0->0), fcn. (221->10), ass. (0->37)
t102 = rSges(5,3) + pkin(8);
t78 = -pkin(9) - pkin(8);
t101 = rSges(6,3) - t78;
t100 = rSges(7,3) + pkin(10) - t78;
t74 = sin(qJ(1));
t77 = cos(qJ(1));
t99 = g(1) * t77 + g(2) * t74;
t72 = sin(qJ(4));
t98 = pkin(4) * t72;
t73 = sin(qJ(2));
t94 = t73 * pkin(2) + pkin(6);
t75 = cos(qJ(4));
t61 = t75 * pkin(4) + pkin(3);
t93 = t73 * t74;
t92 = t73 * t77;
t91 = t74 * t75;
t76 = cos(qJ(2));
t90 = t74 * t76;
t89 = t75 * t77;
t71 = qJ(4) + qJ(5);
t64 = qJ(6) + t71;
t60 = cos(t64);
t88 = t77 * t60;
t85 = t77 * pkin(1) + t74 * pkin(7);
t84 = qJ(3) * t73;
t83 = t72 * t93;
t82 = t72 * t92;
t67 = t74 * pkin(1);
t81 = pkin(2) * t90 + t74 * t84 + t67;
t80 = t85 + (pkin(2) * t76 + t84) * t77;
t79 = -t77 * pkin(7) + t81;
t63 = cos(t71);
t62 = sin(t71);
t59 = sin(t64);
t54 = pkin(5) * t62 + t98;
t53 = pkin(5) * t63 + t61;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t77 - rSges(2,2) * t74) + g(2) * (rSges(2,1) * t74 + rSges(2,2) * t77) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t74 + t85) + g(2) * (rSges(3,1) * t90 - rSges(3,2) * t93 + t67) + g(3) * (rSges(3,1) * t73 + rSges(3,2) * t76 + pkin(6)) + (g(1) * (rSges(3,1) * t76 - rSges(3,2) * t73) + g(2) * (-rSges(3,3) - pkin(7))) * t77) - m(4) * (g(1) * (rSges(4,1) * t74 + t80) + g(2) * (-rSges(4,2) * t90 + rSges(4,3) * t93 + t81) + g(3) * (-rSges(4,2) * t73 + (-rSges(4,3) - qJ(3)) * t76 + t94) + (g(1) * (-rSges(4,2) * t76 + rSges(4,3) * t73) + g(2) * (-rSges(4,1) - pkin(7))) * t77) - m(5) * (g(1) * (t74 * pkin(3) + (t82 + t91) * rSges(5,1) + (-t72 * t74 + t73 * t89) * rSges(5,2) + t80) + g(2) * (-t77 * pkin(3) + (t83 - t89) * rSges(5,1) + (t72 * t77 + t73 * t91) * rSges(5,2) + t79) + g(3) * (t102 * t73 + t94) + (g(3) * (-rSges(5,1) * t72 - rSges(5,2) * t75 - qJ(3)) + t99 * t102) * t76) - m(6) * (g(1) * (t74 * t61 + pkin(4) * t82 + (t62 * t92 + t63 * t74) * rSges(6,1) + (-t62 * t74 + t63 * t92) * rSges(6,2) + t80) + g(2) * (-t77 * t61 + pkin(4) * t83 + (t62 * t93 - t63 * t77) * rSges(6,1) + (t62 * t77 + t63 * t93) * rSges(6,2) + t79) + g(3) * (t101 * t73 + t94) + (g(3) * (-rSges(6,1) * t62 - rSges(6,2) * t63 - qJ(3) - t98) + t99 * t101) * t76) - m(7) * (g(1) * (t74 * t53 + t54 * t92 + (t59 * t92 + t60 * t74) * rSges(7,1) + (-t59 * t74 + t73 * t88) * rSges(7,2) + t80) + g(2) * (-t77 * t53 + t54 * t93 + (t59 * t93 - t88) * rSges(7,1) + (t59 * t77 + t60 * t93) * rSges(7,2) + t79) + g(3) * (t100 * t73 + t94) + (g(3) * (-rSges(7,1) * t59 - rSges(7,2) * t60 - qJ(3) - t54) + t99 * t100) * t76);
U  = t1;

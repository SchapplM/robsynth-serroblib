% Calculate potential energy for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR10_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:00
% EndTime: 2019-03-09 09:34:00
% DurationCPUTime: 0.55s
% Computational Cost: add. (186->110), mult. (237->135), div. (0->0), fcn. (221->10), ass. (0->40)
t75 = -pkin(8) - qJ(4);
t106 = rSges(6,3) - t75;
t105 = rSges(7,3) + pkin(9) - t75;
t104 = rSges(5,3) + qJ(4);
t77 = sin(qJ(1));
t79 = cos(qJ(1));
t103 = g(1) * t79 + g(2) * t77;
t73 = sin(pkin(10));
t102 = pkin(4) * t73;
t76 = sin(qJ(2));
t99 = t76 * pkin(2) + pkin(6);
t74 = cos(pkin(10));
t62 = t74 * pkin(4) + pkin(3);
t98 = t76 * t77;
t97 = t76 * t79;
t72 = pkin(10) + qJ(5);
t65 = qJ(6) + t72;
t60 = sin(t65);
t96 = t77 * t60;
t61 = cos(t65);
t95 = t77 * t61;
t63 = sin(t72);
t94 = t77 * t63;
t64 = cos(t72);
t93 = t77 * t64;
t92 = t77 * t73;
t91 = t77 * t74;
t78 = cos(qJ(2));
t90 = t77 * t78;
t87 = t79 * pkin(1) + t77 * pkin(7);
t86 = qJ(3) * t76;
t84 = t73 * t97;
t83 = t76 * t92;
t69 = t77 * pkin(1);
t82 = pkin(2) * t90 + t77 * t86 + t69;
t81 = t87 + (pkin(2) * t78 + t86) * t79;
t80 = -t79 * pkin(7) + t82;
t55 = pkin(5) * t63 + t102;
t54 = pkin(5) * t64 + t62;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t79 - t77 * rSges(2,2)) + g(2) * (t77 * rSges(2,1) + rSges(2,2) * t79) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t77 * rSges(3,3) + t87) + g(2) * (rSges(3,1) * t90 - rSges(3,2) * t98 + t69) + g(3) * (rSges(3,1) * t76 + rSges(3,2) * t78 + pkin(6)) + (g(1) * (rSges(3,1) * t78 - rSges(3,2) * t76) + g(2) * (-rSges(3,3) - pkin(7))) * t79) - m(4) * (g(1) * (t77 * rSges(4,1) + t81) + g(2) * (-rSges(4,2) * t90 + rSges(4,3) * t98 + t82) + g(3) * (-rSges(4,2) * t76 + (-rSges(4,3) - qJ(3)) * t78 + t99) + (g(1) * (-rSges(4,2) * t78 + rSges(4,3) * t76) + g(2) * (-rSges(4,1) - pkin(7))) * t79) - m(5) * (g(1) * (t77 * pkin(3) + (t84 + t91) * rSges(5,1) + (t74 * t97 - t92) * rSges(5,2) + t81) + g(2) * (-t79 * pkin(3) + (-t74 * t79 + t83) * rSges(5,1) + (t73 * t79 + t76 * t91) * rSges(5,2) + t80) + g(3) * (t104 * t76 + t99) + (g(3) * (-rSges(5,1) * t73 - rSges(5,2) * t74 - qJ(3)) + t103 * t104) * t78) - m(6) * (g(1) * (t77 * t62 + pkin(4) * t84 + (t63 * t97 + t93) * rSges(6,1) + (t64 * t97 - t94) * rSges(6,2) + t81) + g(2) * (-t79 * t62 + pkin(4) * t83 + (-t64 * t79 + t76 * t94) * rSges(6,1) + (t63 * t79 + t76 * t93) * rSges(6,2) + t80) + g(3) * (t106 * t76 + t99) + (g(3) * (-rSges(6,1) * t63 - rSges(6,2) * t64 - qJ(3) - t102) + t103 * t106) * t78) - m(7) * (g(1) * (t77 * t54 + t55 * t97 + (t60 * t97 + t95) * rSges(7,1) + (t61 * t97 - t96) * rSges(7,2) + t81) + g(2) * (-t79 * t54 + t55 * t98 + (-t61 * t79 + t76 * t96) * rSges(7,1) + (t60 * t79 + t76 * t95) * rSges(7,2) + t80) + g(3) * (t105 * t76 + t99) + (g(3) * (-rSges(7,1) * t60 - rSges(7,2) * t61 - qJ(3) - t55) + t103 * t105) * t78);
U  = t1;

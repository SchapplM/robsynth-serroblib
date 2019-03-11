% Calculate potential energy for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:14
% EndTime: 2019-03-09 11:47:14
% DurationCPUTime: 0.44s
% Computational Cost: add. (226->99), mult. (212->120), div. (0->0), fcn. (196->10), ass. (0->44)
t108 = rSges(5,3) + pkin(8);
t83 = -pkin(9) - pkin(8);
t107 = rSges(6,3) - t83;
t106 = rSges(7,3) + qJ(6) - t83;
t79 = sin(qJ(1));
t82 = cos(qJ(1));
t105 = g(1) * t82 + g(2) * t79;
t102 = rSges(3,3) + pkin(7);
t78 = sin(qJ(2));
t100 = t78 * pkin(2) + pkin(6);
t80 = cos(qJ(4));
t65 = t80 * pkin(4) + pkin(3);
t74 = qJ(2) + pkin(10);
t67 = sin(t74);
t99 = rSges(4,2) * t67;
t68 = cos(t74);
t98 = t68 * t79;
t97 = t68 * t82;
t75 = qJ(4) + qJ(5);
t69 = sin(t75);
t96 = t69 * t79;
t95 = t69 * t82;
t70 = cos(t75);
t94 = t70 * t79;
t93 = t70 * t82;
t77 = sin(qJ(4));
t92 = t77 * t79;
t91 = t77 * t82;
t90 = t79 * t80;
t89 = t80 * t82;
t81 = cos(qJ(2));
t66 = pkin(2) * t81 + pkin(1);
t76 = -qJ(3) - pkin(7);
t86 = t79 * t66 + t82 * t76;
t63 = t82 * t66;
t85 = -t79 * t76 + t63;
t84 = rSges(3,1) * t81 - rSges(3,2) * t78 + pkin(1);
t61 = pkin(4) * t77 + pkin(5) * t69;
t60 = pkin(5) * t70 + t65;
t59 = t68 * t93 + t96;
t58 = -t68 * t95 + t94;
t57 = t68 * t94 - t95;
t56 = -t68 * t96 - t93;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t82 - rSges(2,2) * t79) + g(2) * (rSges(2,1) * t79 + rSges(2,2) * t82) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t78 + rSges(3,2) * t81 + pkin(6)) + (g(1) * t84 - g(2) * t102) * t82 + (g(1) * t102 + g(2) * t84) * t79) - m(4) * (g(1) * (rSges(4,1) * t97 - t82 * t99 + t63) + g(2) * (-rSges(4,3) * t82 + t86) + g(3) * (rSges(4,1) * t67 + rSges(4,2) * t68 + t100) + (g(1) * (rSges(4,3) - t76) + g(2) * (rSges(4,1) * t68 - t99)) * t79) - m(5) * (g(1) * (pkin(3) * t97 + (t68 * t89 + t92) * rSges(5,1) + (-t68 * t91 + t90) * rSges(5,2) + t85) + g(2) * (pkin(3) * t98 + (t68 * t90 - t91) * rSges(5,1) + (-t68 * t92 - t89) * rSges(5,2) + t86) + g(3) * (-t108 * t68 + t100) + (g(3) * (rSges(5,1) * t80 - rSges(5,2) * t77 + pkin(3)) + t105 * t108) * t67) - m(6) * (g(1) * (t59 * rSges(6,1) + t58 * rSges(6,2) + pkin(4) * t92 + t65 * t97 + t85) + g(2) * (t57 * rSges(6,1) + t56 * rSges(6,2) - pkin(4) * t91 + t65 * t98 + t86) + g(3) * (-t107 * t68 + t100) + (g(3) * (rSges(6,1) * t70 - rSges(6,2) * t69 + t65) + t105 * t107) * t67) - m(7) * (g(1) * (rSges(7,1) * t59 + rSges(7,2) * t58 + t60 * t97 + t61 * t79 + t85) + g(2) * (rSges(7,1) * t57 + rSges(7,2) * t56 + t60 * t98 - t61 * t82 + t86) + g(3) * (-t106 * t68 + t100) + (g(3) * (rSges(7,1) * t70 - rSges(7,2) * t69 + t60) + t105 * t106) * t67);
U  = t1;

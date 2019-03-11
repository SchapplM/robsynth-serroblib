% Calculate potential energy for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP11_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:18
% EndTime: 2019-03-09 12:45:19
% DurationCPUTime: 0.51s
% Computational Cost: add. (176->105), mult. (237->127), div. (0->0), fcn. (221->8), ass. (0->37)
t107 = rSges(5,3) + pkin(8);
t84 = -pkin(9) - pkin(8);
t106 = rSges(6,3) - t84;
t105 = rSges(7,3) + qJ(6) - t84;
t80 = sin(qJ(1));
t83 = cos(qJ(1));
t104 = g(1) * t83 + g(2) * t80;
t78 = sin(qJ(4));
t103 = pkin(4) * t78;
t79 = sin(qJ(2));
t99 = t79 * pkin(2) + pkin(6);
t81 = cos(qJ(4));
t68 = t81 * pkin(4) + pkin(3);
t98 = t79 * t80;
t97 = t79 * t83;
t96 = t80 * t81;
t82 = cos(qJ(2));
t95 = t80 * t82;
t94 = t81 * t83;
t91 = t83 * pkin(1) + t80 * pkin(7);
t90 = qJ(3) * t79;
t89 = t78 * t98;
t88 = t78 * t97;
t73 = t80 * pkin(1);
t87 = pkin(2) * t95 + t80 * t90 + t73;
t86 = t91 + (pkin(2) * t82 + t90) * t83;
t85 = -t83 * pkin(7) + t87;
t77 = qJ(4) + qJ(5);
t70 = cos(t77);
t69 = sin(t77);
t63 = pkin(5) * t69 + t103;
t62 = pkin(5) * t70 + t68;
t61 = t69 * t98 - t70 * t83;
t60 = t69 * t83 + t70 * t98;
t59 = t69 * t97 + t70 * t80;
t58 = -t69 * t80 + t70 * t97;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t83 - rSges(2,2) * t80) + g(2) * (rSges(2,1) * t80 + rSges(2,2) * t83) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t80 * rSges(3,3) + t91) + g(2) * (rSges(3,1) * t95 - rSges(3,2) * t98 + t73) + g(3) * (rSges(3,1) * t79 + rSges(3,2) * t82 + pkin(6)) + (g(1) * (rSges(3,1) * t82 - rSges(3,2) * t79) + g(2) * (-rSges(3,3) - pkin(7))) * t83) - m(4) * (g(1) * (rSges(4,1) * t80 + t86) + g(2) * (-rSges(4,2) * t95 + rSges(4,3) * t98 + t87) + g(3) * (-rSges(4,2) * t79 + (-rSges(4,3) - qJ(3)) * t82 + t99) + (g(1) * (-rSges(4,2) * t82 + rSges(4,3) * t79) + g(2) * (-rSges(4,1) - pkin(7))) * t83) - m(5) * (g(1) * (t80 * pkin(3) + (t88 + t96) * rSges(5,1) + (-t78 * t80 + t79 * t94) * rSges(5,2) + t86) + g(2) * (-t83 * pkin(3) + (t89 - t94) * rSges(5,1) + (t78 * t83 + t79 * t96) * rSges(5,2) + t85) + g(3) * (t107 * t79 + t99) + (g(3) * (-rSges(5,1) * t78 - rSges(5,2) * t81 - qJ(3)) + t104 * t107) * t82) - m(6) * (g(1) * (t59 * rSges(6,1) + t58 * rSges(6,2) + pkin(4) * t88 + t80 * t68 + t86) + g(2) * (t61 * rSges(6,1) + t60 * rSges(6,2) + pkin(4) * t89 - t83 * t68 + t85) + g(3) * (t106 * t79 + t99) + (g(3) * (-rSges(6,1) * t69 - rSges(6,2) * t70 - qJ(3) - t103) + t104 * t106) * t82) - m(7) * (g(1) * (rSges(7,1) * t59 + rSges(7,2) * t58 + t62 * t80 + t63 * t97 + t86) + g(2) * (rSges(7,1) * t61 + rSges(7,2) * t60 - t62 * t83 + t63 * t98 + t85) + g(3) * (t105 * t79 + t99) + (g(3) * (-rSges(7,1) * t69 - rSges(7,2) * t70 - qJ(3) - t63) + t104 * t105) * t82);
U  = t1;

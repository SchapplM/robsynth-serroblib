% Calculate potential energy for
% S6RRPRRP12
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRP12_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP12_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:16
% EndTime: 2019-03-09 12:50:17
% DurationCPUTime: 0.50s
% Computational Cost: add. (183->101), mult. (252->122), div. (0->0), fcn. (240->8), ass. (0->38)
t108 = rSges(7,1) + pkin(5);
t107 = rSges(5,3) + pkin(8);
t106 = rSges(7,3) + qJ(6);
t78 = sin(qJ(1));
t81 = cos(qJ(1));
t105 = g(1) * t81 + g(2) * t78;
t77 = sin(qJ(2));
t101 = t77 * pkin(2) + pkin(6);
t100 = t77 * t78;
t99 = t77 * t81;
t79 = cos(qJ(4));
t98 = t78 * t79;
t80 = cos(qJ(2));
t97 = t78 * t80;
t96 = t79 * t81;
t93 = t81 * pkin(1) + t78 * pkin(7);
t92 = qJ(3) * t77;
t76 = sin(qJ(4));
t91 = t76 * t100;
t90 = t76 * t99;
t73 = t78 * pkin(1);
t89 = pkin(2) * t97 + t78 * t92 + t73;
t88 = -pkin(4) * t76 - qJ(3);
t87 = t93 + (pkin(2) * t80 + t92) * t81;
t82 = -pkin(9) - pkin(8);
t86 = -t77 * t82 + t101;
t85 = -t81 * pkin(7) + t89;
t68 = pkin(4) * t79 + pkin(3);
t84 = pkin(4) * t90 + t78 * t68 + t87;
t83 = pkin(4) * t91 - t81 * t68 + t85;
t75 = qJ(4) + qJ(5);
t70 = cos(t75);
t69 = sin(t75);
t59 = t69 * t100 - t70 * t81;
t58 = t70 * t100 + t69 * t81;
t57 = t69 * t99 + t70 * t78;
t56 = t69 * t78 - t70 * t99;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t81 - rSges(2,2) * t78) + g(2) * (rSges(2,1) * t78 + rSges(2,2) * t81) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t78 + t93) + g(2) * (rSges(3,1) * t97 - rSges(3,2) * t100 + t73) + g(3) * (rSges(3,1) * t77 + rSges(3,2) * t80 + pkin(6)) + (g(1) * (rSges(3,1) * t80 - rSges(3,2) * t77) + g(2) * (-rSges(3,3) - pkin(7))) * t81) - m(4) * (g(1) * (rSges(4,1) * t78 + t87) + g(2) * (-rSges(4,2) * t97 + rSges(4,3) * t100 + t89) + g(3) * (-rSges(4,2) * t77 + (-rSges(4,3) - qJ(3)) * t80 + t101) + (g(1) * (-rSges(4,2) * t80 + rSges(4,3) * t77) + g(2) * (-rSges(4,1) - pkin(7))) * t81) - m(5) * (g(1) * (t78 * pkin(3) + (t90 + t98) * rSges(5,1) + (-t76 * t78 + t77 * t96) * rSges(5,2) + t87) + g(2) * (-t81 * pkin(3) + (t91 - t96) * rSges(5,1) + (t76 * t81 + t77 * t98) * rSges(5,2) + t85) + g(3) * (t107 * t77 + t101) + (g(3) * (-rSges(5,1) * t76 - rSges(5,2) * t79 - qJ(3)) + t105 * t107) * t80) - m(6) * (g(1) * (t57 * rSges(6,1) - t56 * rSges(6,2) + t84) + g(2) * (t59 * rSges(6,1) + t58 * rSges(6,2) + t83) + g(3) * (t77 * rSges(6,3) + t86) + (g(3) * (-rSges(6,1) * t69 - rSges(6,2) * t70 + t88) + t105 * (rSges(6,3) - t82)) * t80) - m(7) * (g(1) * (t106 * t56 + t108 * t57 + t84) + g(2) * (-t106 * t58 + t108 * t59 + t83) + g(3) * (t77 * rSges(7,2) + t86) + (g(3) * (t106 * t70 - t108 * t69 + t88) + t105 * (rSges(7,2) - t82)) * t80);
U  = t1;

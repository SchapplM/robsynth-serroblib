% Calculate potential energy for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:56:56
% EndTime: 2019-03-09 06:56:56
% DurationCPUTime: 0.45s
% Computational Cost: add. (234->93), mult. (172->115), div. (0->0), fcn. (152->12), ass. (0->40)
t100 = rSges(6,3) + pkin(9);
t99 = rSges(7,3) + pkin(10) + pkin(9);
t67 = qJ(1) + pkin(11);
t58 = sin(t67);
t59 = cos(t67);
t98 = g(1) * t59 + g(2) * t58;
t95 = rSges(4,3) + pkin(7);
t69 = qJ(3) + qJ(4);
t61 = sin(t69);
t93 = rSges(5,2) * t61;
t63 = cos(t69);
t92 = t58 * t63;
t70 = sin(qJ(5));
t91 = t58 * t70;
t90 = t59 * t63;
t89 = t59 * t70;
t68 = qJ(5) + qJ(6);
t60 = sin(t68);
t88 = t60 * t63;
t62 = cos(t68);
t87 = t62 * t63;
t86 = t63 * t70;
t73 = cos(qJ(5));
t85 = t63 * t73;
t83 = pkin(6) + qJ(2);
t74 = cos(qJ(3));
t57 = pkin(3) * t74 + pkin(2);
t75 = cos(qJ(1));
t66 = t75 * pkin(1);
t82 = t59 * t57 + t66;
t71 = sin(qJ(3));
t81 = t71 * pkin(3) + t83;
t72 = sin(qJ(1));
t65 = t72 * pkin(1);
t77 = -pkin(8) - pkin(7);
t80 = t58 * t57 + t59 * t77 + t65;
t79 = -t58 * t77 + t82;
t78 = rSges(4,1) * t74 - rSges(4,2) * t71 + pkin(2);
t56 = pkin(5) * t73 + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t75 - rSges(2,2) * t72) + g(2) * (rSges(2,1) * t72 + rSges(2,2) * t75) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t59 - rSges(3,2) * t58 + t66) + g(2) * (rSges(3,1) * t58 + rSges(3,2) * t59 + t65) + g(3) * (rSges(3,3) + t83)) - m(4) * (g(1) * t66 + g(2) * t65 + g(3) * (rSges(4,1) * t71 + rSges(4,2) * t74 + t83) + (g(1) * t78 - g(2) * t95) * t59 + (g(1) * t95 + g(2) * t78) * t58) - m(5) * (g(1) * (rSges(5,1) * t90 - t59 * t93 + t82) + g(2) * (-rSges(5,3) * t59 + t80) + g(3) * (rSges(5,1) * t61 + rSges(5,2) * t63 + t81) + (g(1) * (rSges(5,3) - t77) + g(2) * (rSges(5,1) * t63 - t93)) * t58) - m(6) * (g(1) * (pkin(4) * t90 + (t59 * t85 + t91) * rSges(6,1) + (t58 * t73 - t59 * t86) * rSges(6,2) + t79) + g(2) * (pkin(4) * t92 + (t58 * t85 - t89) * rSges(6,1) + (-t58 * t86 - t59 * t73) * rSges(6,2) + t80) + g(3) * (-t100 * t63 + t81) + (g(3) * (rSges(6,1) * t73 - rSges(6,2) * t70 + pkin(4)) + t98 * t100) * t61) - m(7) * (g(1) * (t56 * t90 + pkin(5) * t91 + (t58 * t60 + t59 * t87) * rSges(7,1) + (t58 * t62 - t59 * t88) * rSges(7,2) + t79) + g(2) * (t56 * t92 - pkin(5) * t89 + (t58 * t87 - t59 * t60) * rSges(7,1) + (-t58 * t88 - t59 * t62) * rSges(7,2) + t80) + g(3) * (-t99 * t63 + t81) + (g(3) * (rSges(7,1) * t62 - rSges(7,2) * t60 + t56) + t98 * t99) * t61);
U  = t1;

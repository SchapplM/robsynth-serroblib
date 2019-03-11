% Calculate potential energy for
% S6RPRRRR3
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR3_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:12
% EndTime: 2019-03-09 07:00:12
% DurationCPUTime: 0.50s
% Computational Cost: add. (245->102), mult. (198->130), div. (0->0), fcn. (182->12), ass. (0->41)
t103 = rSges(5,3) + pkin(8);
t79 = -pkin(9) - pkin(8);
t102 = rSges(6,3) - t79;
t101 = rSges(7,3) + pkin(10) - t79;
t70 = qJ(1) + pkin(11);
t62 = sin(t70);
t63 = cos(t70);
t100 = g(1) * t63 + g(2) * t62;
t76 = cos(qJ(4));
t61 = t76 * pkin(4) + pkin(3);
t74 = sin(qJ(3));
t96 = rSges(4,2) * t74;
t73 = sin(qJ(4));
t95 = t62 * t73;
t77 = cos(qJ(3));
t94 = t62 * t77;
t93 = t63 * t73;
t92 = t63 * t77;
t72 = qJ(4) + qJ(5);
t64 = sin(t72);
t91 = t64 * t77;
t65 = cos(t72);
t90 = t65 * t77;
t89 = t73 * t77;
t88 = t76 * t77;
t54 = pkin(5) * t65 + t61;
t87 = t77 * t54;
t86 = t77 * t61;
t83 = pkin(6) + qJ(2);
t75 = sin(qJ(1));
t67 = t75 * pkin(1);
t82 = t62 * pkin(2) + t67;
t78 = cos(qJ(1));
t69 = t78 * pkin(1);
t81 = t63 * pkin(2) + t62 * pkin(7) + t69;
t80 = -t63 * pkin(7) + t82;
t66 = qJ(6) + t72;
t60 = cos(t66);
t59 = sin(t66);
t55 = t73 * pkin(4) + pkin(5) * t64;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t78 - rSges(2,2) * t75) + g(2) * (rSges(2,1) * t75 + rSges(2,2) * t78) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t63 - rSges(3,2) * t62 + t69) + g(2) * (rSges(3,1) * t62 + rSges(3,2) * t63 + t67) + g(3) * (rSges(3,3) + t83)) - m(4) * (g(1) * (rSges(4,3) * t62 + t81) + g(2) * (rSges(4,1) * t94 - t62 * t96 + t82) + g(3) * (rSges(4,1) * t74 + rSges(4,2) * t77 + t83) + (g(1) * (rSges(4,1) * t77 - t96) + g(2) * (-rSges(4,3) - pkin(7))) * t63) - m(5) * (g(1) * (pkin(3) * t92 + (t63 * t88 + t95) * rSges(5,1) + (t62 * t76 - t63 * t89) * rSges(5,2) + t81) + g(2) * (pkin(3) * t94 + (t62 * t88 - t93) * rSges(5,1) + (-t62 * t89 - t63 * t76) * rSges(5,2) + t80) + g(3) * (-t103 * t77 + t83) + (g(3) * (rSges(5,1) * t76 - rSges(5,2) * t73 + pkin(3)) + t100 * t103) * t74) - m(6) * (g(1) * (t63 * t86 + pkin(4) * t95 + (t62 * t64 + t63 * t90) * rSges(6,1) + (t62 * t65 - t63 * t91) * rSges(6,2) + t81) + g(2) * (t62 * t86 - pkin(4) * t93 + (t62 * t90 - t63 * t64) * rSges(6,1) + (-t62 * t91 - t63 * t65) * rSges(6,2) + t80) + g(3) * (-t102 * t77 + t83) + (g(3) * (rSges(6,1) * t65 - rSges(6,2) * t64 + t61) + t100 * t102) * t74) - m(7) * (g(1) * (t63 * t87 + t62 * t55 + (t59 * t62 + t60 * t92) * rSges(7,1) + (-t59 * t92 + t60 * t62) * rSges(7,2) + t81) + g(2) * (t62 * t87 - t63 * t55 + (-t59 * t63 + t60 * t94) * rSges(7,1) + (-t59 * t94 - t60 * t63) * rSges(7,2) + t80) + g(3) * (-t101 * t77 + t83) + (g(3) * (rSges(7,1) * t60 - rSges(7,2) * t59 + t54) + t100 * t101) * t74);
U  = t1;

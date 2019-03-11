% Calculate potential energy for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
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
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_energypot_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR1_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR1_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:37:54
% EndTime: 2019-03-09 02:37:55
% DurationCPUTime: 0.43s
% Computational Cost: add. (234->93), mult. (172->111), div. (0->0), fcn. (152->12), ass. (0->38)
t95 = rSges(7,3) + pkin(8) + qJ(5);
t66 = qJ(1) + pkin(9);
t57 = sin(t66);
t60 = cos(t66);
t94 = g(1) * t60 + g(2) * t57;
t93 = rSges(6,3) + qJ(5);
t90 = rSges(4,3) + pkin(7);
t65 = qJ(3) + pkin(10);
t56 = sin(t65);
t89 = rSges(5,2) * t56;
t59 = cos(t65);
t88 = t57 * t59;
t67 = sin(pkin(11));
t87 = t57 * t67;
t68 = cos(pkin(11));
t86 = t57 * t68;
t85 = t60 * t59;
t84 = t60 * t67;
t83 = t60 * t68;
t81 = pkin(6) + qJ(2);
t73 = cos(qJ(3));
t54 = t73 * pkin(3) + pkin(2);
t74 = cos(qJ(1));
t63 = t74 * pkin(1);
t80 = t60 * t54 + t63;
t71 = sin(qJ(3));
t78 = t71 * pkin(3) + t81;
t72 = sin(qJ(1));
t62 = t72 * pkin(1);
t69 = -qJ(4) - pkin(7);
t77 = t57 * t54 + t60 * t69 + t62;
t76 = -t57 * t69 + t80;
t75 = rSges(4,1) * t73 - rSges(4,2) * t71 + pkin(2);
t64 = pkin(11) + qJ(6);
t58 = cos(t64);
t55 = sin(t64);
t53 = t68 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t74 * rSges(2,1) - t72 * rSges(2,2)) + g(2) * (t72 * rSges(2,1) + t74 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (t60 * rSges(3,1) - t57 * rSges(3,2) + t63) + g(2) * (t57 * rSges(3,1) + t60 * rSges(3,2) + t62) + g(3) * (rSges(3,3) + t81)) - m(4) * (g(1) * t63 + g(2) * t62 + g(3) * (t71 * rSges(4,1) + t73 * rSges(4,2) + t81) + (g(1) * t75 - g(2) * t90) * t60 + (g(1) * t90 + g(2) * t75) * t57) - m(5) * (g(1) * (rSges(5,1) * t85 - t60 * t89 + t80) + g(2) * (-t60 * rSges(5,3) + t77) + g(3) * (t56 * rSges(5,1) + t59 * rSges(5,2) + t78) + (g(1) * (rSges(5,3) - t69) + g(2) * (rSges(5,1) * t59 - t89)) * t57) - m(6) * (g(1) * (pkin(4) * t85 + (t59 * t83 + t87) * rSges(6,1) + (-t59 * t84 + t86) * rSges(6,2) + t76) + g(2) * (pkin(4) * t88 + (t59 * t86 - t84) * rSges(6,1) + (-t59 * t87 - t83) * rSges(6,2) + t77) + g(3) * (-t93 * t59 + t78) + (g(3) * (rSges(6,1) * t68 - rSges(6,2) * t67 + pkin(4)) + t94 * t93) * t56) - m(7) * (g(1) * (t53 * t85 + pkin(5) * t87 + (t57 * t55 + t58 * t85) * rSges(7,1) + (-t55 * t85 + t57 * t58) * rSges(7,2) + t76) + g(2) * (t53 * t88 - pkin(5) * t84 + (-t60 * t55 + t58 * t88) * rSges(7,1) + (-t55 * t88 - t60 * t58) * rSges(7,2) + t77) + g(3) * (-t95 * t59 + t78) + (g(3) * (rSges(7,1) * t58 - rSges(7,2) * t55 + t53) + t94 * t95) * t56);
U  = t1;

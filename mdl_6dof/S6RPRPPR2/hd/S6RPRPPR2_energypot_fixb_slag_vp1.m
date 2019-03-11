% Calculate potential energy for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRPPR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:10
% EndTime: 2019-03-09 02:41:10
% DurationCPUTime: 0.37s
% Computational Cost: add. (214->88), mult. (159->103), div. (0->0), fcn. (135->10), ass. (0->34)
t89 = rSges(7,3) + pkin(8);
t88 = rSges(4,3) + pkin(7);
t63 = qJ(3) + pkin(10);
t56 = sin(t63);
t64 = qJ(1) + pkin(9);
t59 = cos(t64);
t86 = t56 * t59;
t57 = sin(t64);
t66 = sin(qJ(6));
t85 = t57 * t66;
t69 = cos(qJ(6));
t84 = t57 * t69;
t58 = cos(t63);
t83 = t58 * t59;
t82 = t59 * t66;
t81 = t59 * t69;
t80 = pkin(6) + qJ(2);
t70 = cos(qJ(3));
t55 = pkin(3) * t70 + pkin(2);
t71 = cos(qJ(1));
t62 = t71 * pkin(1);
t79 = t59 * t55 + t62;
t78 = qJ(5) * t56;
t67 = sin(qJ(3));
t77 = t67 * pkin(3) + t80;
t68 = sin(qJ(1));
t61 = t68 * pkin(1);
t65 = -qJ(4) - pkin(7);
t76 = t57 * t55 + t59 * t65 + t61;
t75 = t56 * pkin(4) + t77;
t74 = pkin(4) * t83 + t59 * t78 + t79;
t73 = t76 + (pkin(4) * t58 + t78) * t57;
t72 = rSges(4,1) * t70 - rSges(4,2) * t67 + pkin(2);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t71 - t68 * rSges(2,2)) + g(2) * (t68 * rSges(2,1) + rSges(2,2) * t71) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t59 - rSges(3,2) * t57 + t62) + g(2) * (rSges(3,1) * t57 + rSges(3,2) * t59 + t61) + g(3) * (rSges(3,3) + t80)) - m(4) * (g(1) * t62 + g(2) * t61 + g(3) * (rSges(4,1) * t67 + rSges(4,2) * t70 + t80) + (g(1) * t72 - g(2) * t88) * t59 + (g(1) * t88 + g(2) * t72) * t57) - m(5) * (g(1) * (rSges(5,1) * t83 - rSges(5,2) * t86 + t79) + g(2) * (-rSges(5,3) * t59 + t76) + g(3) * (rSges(5,1) * t56 + rSges(5,2) * t58 + t77) + (g(1) * (rSges(5,3) - t65) + g(2) * (rSges(5,1) * t58 - rSges(5,2) * t56)) * t57) - m(6) * (g(1) * (-rSges(6,2) * t83 + rSges(6,3) * t86 + t74) + g(2) * (-rSges(6,1) * t59 + t73) + g(3) * (-rSges(6,2) * t56 + (-rSges(6,3) - qJ(5)) * t58 + t75) + (g(1) * (rSges(6,1) - t65) + g(2) * (-rSges(6,2) * t58 + rSges(6,3) * t56)) * t57) - m(7) * (g(1) * ((t56 * t82 + t84) * rSges(7,1) + (t56 * t81 - t85) * rSges(7,2) + t74 + (pkin(5) - t65) * t57) + g(2) * (-t59 * pkin(5) + (t56 * t85 - t81) * rSges(7,1) + (t56 * t84 + t82) * rSges(7,2) + t73) + g(3) * (t89 * t56 + t75) + (g(3) * (-rSges(7,1) * t66 - rSges(7,2) * t69 - qJ(5)) + (g(1) * t59 + g(2) * t57) * t89) * t58);
U  = t1;

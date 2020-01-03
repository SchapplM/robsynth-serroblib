% Calculate potential energy for
% S5RRPPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRPPR11_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR11_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:20
% EndTime: 2019-12-31 19:46:20
% DurationCPUTime: 0.42s
% Computational Cost: add. (116->84), mult. (172->107), div. (0->0), fcn. (156->8), ass. (0->30)
t84 = rSges(6,3) + pkin(7) + qJ(4);
t60 = sin(qJ(1));
t62 = cos(qJ(1));
t83 = g(1) * t62 + g(2) * t60;
t82 = rSges(5,3) + qJ(4);
t59 = sin(qJ(2));
t79 = t59 * pkin(2) + pkin(5);
t78 = t59 * t60;
t77 = t59 * t62;
t55 = pkin(8) + qJ(5);
t49 = sin(t55);
t76 = t60 * t49;
t50 = cos(t55);
t75 = t60 * t50;
t56 = sin(pkin(8));
t74 = t60 * t56;
t57 = cos(pkin(8));
t73 = t60 * t57;
t61 = cos(qJ(2));
t72 = t60 * t61;
t70 = t62 * pkin(1) + t60 * pkin(6);
t69 = qJ(3) * t59;
t67 = t56 * t77;
t66 = t59 * t74;
t53 = t60 * pkin(1);
t65 = pkin(2) * t72 + t60 * t69 + t53;
t64 = t70 + (pkin(2) * t61 + t69) * t62;
t63 = -t62 * pkin(6) + t65;
t48 = pkin(4) * t57 + pkin(3);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t62 - t60 * rSges(2,2)) + g(2) * (t60 * rSges(2,1) + rSges(2,2) * t62) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (t60 * rSges(3,3) + t70) + g(2) * (rSges(3,1) * t72 - rSges(3,2) * t78 + t53) + g(3) * (t59 * rSges(3,1) + t61 * rSges(3,2) + pkin(5)) + (g(1) * (rSges(3,1) * t61 - rSges(3,2) * t59) + g(2) * (-rSges(3,3) - pkin(6))) * t62) - m(4) * (g(1) * (t60 * rSges(4,1) + t64) + g(2) * (-rSges(4,2) * t72 + rSges(4,3) * t78 + t65) + g(3) * (-rSges(4,2) * t59 + (-rSges(4,3) - qJ(3)) * t61 + t79) + (g(1) * (-rSges(4,2) * t61 + rSges(4,3) * t59) + g(2) * (-rSges(4,1) - pkin(6))) * t62) - m(5) * (g(1) * (t60 * pkin(3) + (t67 + t73) * rSges(5,1) + (t57 * t77 - t74) * rSges(5,2) + t64) + g(2) * (-t62 * pkin(3) + (-t57 * t62 + t66) * rSges(5,1) + (t56 * t62 + t59 * t73) * rSges(5,2) + t63) + g(3) * (t82 * t59 + t79) + (g(3) * (-rSges(5,1) * t56 - rSges(5,2) * t57 - qJ(3)) + t83 * t82) * t61) - m(6) * (g(1) * (t60 * t48 + pkin(4) * t67 + (t49 * t77 + t75) * rSges(6,1) + (t50 * t77 - t76) * rSges(6,2) + t64) + g(2) * (-t62 * t48 + pkin(4) * t66 + (-t50 * t62 + t59 * t76) * rSges(6,1) + (t49 * t62 + t59 * t75) * rSges(6,2) + t63) + g(3) * (t84 * t59 + t79) + (g(3) * (-rSges(6,1) * t49 - rSges(6,2) * t50 - pkin(4) * t56 - qJ(3)) + t83 * t84) * t61);
U = t1;

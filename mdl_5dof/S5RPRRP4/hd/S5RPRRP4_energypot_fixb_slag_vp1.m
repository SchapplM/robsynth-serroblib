% Calculate potential energy for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:43
% EndTime: 2022-01-23 09:31:43
% DurationCPUTime: 0.41s
% Computational Cost: add. (142->83), mult. (173->101), div. (0->0), fcn. (161->8), ass. (0->35)
t90 = rSges(4,3) + pkin(6);
t73 = pkin(7) + pkin(6);
t89 = rSges(6,3) + qJ(5) + t73;
t70 = sin(qJ(1));
t72 = cos(qJ(1));
t88 = g(1) * t72 + g(2) * t70;
t69 = sin(qJ(3));
t87 = pkin(3) * t69;
t71 = cos(qJ(3));
t86 = pkin(3) * t71;
t68 = cos(pkin(8));
t83 = pkin(2) * t68 + pkin(1);
t67 = sin(pkin(8));
t82 = t67 * pkin(2) + pkin(5);
t81 = rSges(3,2) * t67;
t80 = t68 * t70;
t79 = t68 * t72;
t61 = t70 * qJ(2);
t77 = t72 * pkin(1) + t61;
t76 = rSges(4,1) * t71 - rSges(4,2) * t69;
t75 = t69 * rSges(4,1) + t71 * rSges(4,2);
t74 = t90 * t67 + t76 * t68 + t83;
t66 = qJ(3) + qJ(4);
t63 = t70 * pkin(1);
t59 = cos(t66);
t58 = sin(t66);
t57 = qJ(2) + t87;
t55 = pkin(4) * t58 + t87;
t54 = pkin(4) * t59 + pkin(2) + t86;
t53 = t58 * t70 + t59 * t79;
t52 = -t58 * t79 + t59 * t70;
t51 = -t58 * t72 + t59 * t80;
t50 = -t58 * t80 - t59 * t72;
t49 = t67 * t73 + t68 * t86 + t83;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t72 - rSges(2,2) * t70) + g(2) * (rSges(2,1) * t70 + rSges(2,2) * t72) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t70 + t77) + g(2) * (rSges(3,1) * t80 - t70 * t81 + t63) + g(3) * (rSges(3,1) * t67 + rSges(3,2) * t68 + pkin(5)) + (g(1) * (rSges(3,1) * t68 - t81) + g(2) * (-rSges(3,3) - qJ(2))) * t72) - m(4) * (g(1) * (t75 * t70 + t74 * t72 + t61) + g(2) * ((-qJ(2) - t75) * t72 + t74 * t70) + g(3) * (t76 * t67 - t90 * t68 + t82)) - m(5) * (g(1) * (rSges(5,1) * t53 + rSges(5,2) * t52 + t49 * t72 + t57 * t70) + g(2) * (rSges(5,1) * t51 + rSges(5,2) * t50 + t49 * t70 - t57 * t72) + g(3) * (t82 + (-rSges(5,3) - t73) * t68) + (g(3) * (rSges(5,1) * t59 - rSges(5,2) * t58 + t86) + t88 * rSges(5,3)) * t67) - m(6) * (g(1) * (rSges(6,1) * t53 + rSges(6,2) * t52 + t54 * t79 + t55 * t70 + t77) + g(2) * (rSges(6,1) * t51 + rSges(6,2) * t50 + t54 * t80 + t63 + (-qJ(2) - t55) * t72) + g(3) * (-t89 * t68 + pkin(5)) + (g(3) * (rSges(6,1) * t59 - rSges(6,2) * t58 + t54) + t88 * t89) * t67);
U = t1;

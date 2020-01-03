% Calculate potential energy for
% S5RPRRP10
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRRP10_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP10_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:51
% EndTime: 2019-12-31 18:50:51
% DurationCPUTime: 0.33s
% Computational Cost: add. (141->74), mult. (154->91), div. (0->0), fcn. (138->8), ass. (0->33)
t82 = rSges(5,3) + pkin(7);
t81 = rSges(6,3) + qJ(5) + pkin(7);
t61 = sin(qJ(1));
t63 = cos(qJ(1));
t80 = g(1) * t63 + g(2) * t61;
t56 = sin(pkin(8));
t76 = t56 * pkin(2) + pkin(5);
t55 = pkin(8) + qJ(3);
t52 = sin(t55);
t75 = rSges(4,2) * t52;
t53 = cos(t55);
t74 = t53 * t61;
t73 = t53 * t63;
t60 = sin(qJ(4));
t72 = t60 * t63;
t71 = t61 * t60;
t62 = cos(qJ(4));
t70 = t61 * t62;
t69 = t62 * t63;
t57 = cos(pkin(8));
t49 = pkin(2) * t57 + pkin(1);
t59 = -pkin(6) - qJ(2);
t67 = t61 * t49 + t63 * t59;
t66 = rSges(3,3) + qJ(2);
t48 = t63 * t49;
t65 = -t61 * t59 + t48;
t64 = rSges(3,1) * t57 - rSges(3,2) * t56 + pkin(1);
t51 = pkin(4) * t62 + pkin(3);
t46 = t53 * t69 + t71;
t45 = -t53 * t72 + t70;
t44 = t53 * t70 - t72;
t43 = -t53 * t71 - t69;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t63 - t61 * rSges(2,2)) + g(2) * (t61 * rSges(2,1) + rSges(2,2) * t63) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t56 + rSges(3,2) * t57 + pkin(5)) + (g(1) * t64 - g(2) * t66) * t63 + (g(1) * t66 + g(2) * t64) * t61) - m(4) * (g(1) * (rSges(4,1) * t73 - t63 * t75 + t48) + g(2) * (-rSges(4,3) * t63 + t67) + g(3) * (rSges(4,1) * t52 + rSges(4,2) * t53 + t76) + (g(1) * (rSges(4,3) - t59) + g(2) * (rSges(4,1) * t53 - t75)) * t61) - m(5) * (g(1) * (t46 * rSges(5,1) + t45 * rSges(5,2) + pkin(3) * t73 + t65) + g(2) * (rSges(5,1) * t44 + rSges(5,2) * t43 + pkin(3) * t74 + t67) + g(3) * (-t82 * t53 + t76) + (g(3) * (rSges(5,1) * t62 - rSges(5,2) * t60 + pkin(3)) + t80 * t82) * t52) - m(6) * (g(1) * (t46 * rSges(6,1) + t45 * rSges(6,2) + pkin(4) * t71 + t51 * t73 + t65) + g(2) * (t44 * rSges(6,1) + t43 * rSges(6,2) - pkin(4) * t72 + t51 * t74 + t67) + g(3) * (-t81 * t53 + t76) + (g(3) * (rSges(6,1) * t62 - rSges(6,2) * t60 + t51) + t80 * t81) * t52);
U = t1;

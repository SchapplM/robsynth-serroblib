% Calculate potential energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:43:58
% EndTime: 2019-12-05 17:43:59
% DurationCPUTime: 0.49s
% Computational Cost: add. (152->90), mult. (180->114), div. (0->0), fcn. (168->10), ass. (0->39)
t60 = -pkin(6) - qJ(3);
t88 = rSges(5,3) - t60;
t87 = rSges(6,3) + pkin(7) - t60;
t61 = sin(qJ(1));
t62 = cos(qJ(1));
t86 = -g(2) * t61 + g(3) * t62;
t85 = rSges(4,3) + qJ(3);
t58 = cos(pkin(9));
t46 = t58 * pkin(3) + pkin(2);
t57 = sin(pkin(8));
t82 = rSges(3,2) * t57;
t59 = cos(pkin(8));
t81 = t59 * t61;
t80 = t59 * t62;
t55 = pkin(9) + qJ(4);
t49 = qJ(5) + t55;
t44 = sin(t49);
t79 = t61 * t44;
t45 = cos(t49);
t78 = t61 * t45;
t47 = sin(t55);
t77 = t61 * t47;
t48 = cos(t55);
t76 = t61 * t48;
t56 = sin(pkin(9));
t75 = t61 * t56;
t74 = t61 * t58;
t73 = t62 * t44;
t72 = t62 * t45;
t71 = t62 * t47;
t70 = t62 * t48;
t69 = t62 * t56;
t68 = t62 * t58;
t65 = t62 * pkin(1) + t61 * qJ(2);
t52 = t62 * qJ(2);
t63 = -t61 * pkin(1) + t52;
t43 = t56 * pkin(3) + pkin(4) * t47;
t42 = pkin(4) * t48 + t46;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (-t61 * rSges(2,1) - t62 * rSges(2,2)) + g(3) * (t62 * rSges(2,1) - t61 * rSges(2,2))) - m(3) * (g(1) * (t57 * rSges(3,1) + t59 * rSges(3,2) + pkin(5)) + g(2) * (t62 * rSges(3,3) + t52) + g(3) * (rSges(3,1) * t80 - t62 * t82 + t65) + (g(2) * (-rSges(3,1) * t59 - pkin(1) + t82) + g(3) * rSges(3,3)) * t61) - m(4) * (g(1) * (-t85 * t59 + pkin(5)) + g(2) * (-pkin(2) * t81 + (-t59 * t74 + t69) * rSges(4,1) + (t59 * t75 + t68) * rSges(4,2) + t63) + g(3) * (pkin(2) * t80 + (t59 * t68 + t75) * rSges(4,1) + (-t59 * t69 + t74) * rSges(4,2) + t65) + (g(1) * (rSges(4,1) * t58 - rSges(4,2) * t56 + pkin(2)) + t86 * t85) * t57) - m(5) * (g(1) * (-t88 * t59 + pkin(5)) + g(2) * (-t46 * t81 + pkin(3) * t69 + (-t59 * t76 + t71) * rSges(5,1) + (t59 * t77 + t70) * rSges(5,2) + t63) + g(3) * (t46 * t80 + pkin(3) * t75 + (t59 * t70 + t77) * rSges(5,1) + (-t59 * t71 + t76) * rSges(5,2) + t65) + (g(1) * (rSges(5,1) * t48 - rSges(5,2) * t47 + t46) + t86 * t88) * t57) - m(6) * (g(1) * (-t87 * t59 + pkin(5)) + g(2) * (-t42 * t81 + t62 * t43 + (-t59 * t78 + t73) * rSges(6,1) + (t59 * t79 + t72) * rSges(6,2) + t63) + g(3) * (t42 * t80 + t61 * t43 + (t59 * t72 + t79) * rSges(6,1) + (-t59 * t73 + t78) * rSges(6,2) + t65) + (g(1) * (rSges(6,1) * t45 - rSges(6,2) * t44 + t42) + t86 * t87) * t57);
U = t1;

% Calculate potential energy for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:17
% EndTime: 2019-12-05 15:16:18
% DurationCPUTime: 0.44s
% Computational Cost: add. (150->87), mult. (254->114), div. (0->0), fcn. (266->10), ass. (0->34)
t88 = pkin(7) + pkin(6) + rSges(6,3);
t87 = pkin(6) + rSges(5,3);
t63 = sin(pkin(9));
t64 = sin(pkin(8));
t86 = t63 * t64;
t66 = cos(pkin(8));
t85 = t63 * t66;
t67 = sin(qJ(4));
t84 = t63 * t67;
t69 = cos(qJ(4));
t83 = t63 * t69;
t65 = cos(pkin(9));
t82 = t64 * t65;
t68 = sin(qJ(3));
t81 = t64 * t68;
t70 = cos(qJ(3));
t80 = t64 * t70;
t79 = t66 * t68;
t78 = t66 * t70;
t77 = t66 * pkin(1) + t64 * qJ(2);
t76 = t63 * pkin(2) + qJ(1);
t75 = t66 * t65 * pkin(2) + pkin(5) * t85 + t77;
t60 = t64 * pkin(1);
t74 = pkin(2) * t82 + pkin(5) * t86 - t66 * qJ(2) + t60;
t62 = qJ(4) + qJ(5);
t57 = sin(t62);
t58 = cos(t62);
t73 = rSges(6,1) * t58 - rSges(6,2) * t57 + pkin(4) * t69 + pkin(3);
t72 = t57 * rSges(6,1) + t58 * rSges(6,2) + t67 * pkin(4);
t49 = t65 * t78 + t81;
t48 = t65 * t79 - t80;
t47 = t65 * t80 - t79;
t46 = t65 * t81 + t78;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t66 - rSges(2,2) * t64) + g(2) * (rSges(2,1) * t64 + rSges(2,2) * t66) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t64 + t77) + g(2) * (rSges(3,1) * t82 - rSges(3,2) * t86 + t60) + g(3) * (rSges(3,1) * t63 + rSges(3,2) * t65 + qJ(1)) + (g(1) * (rSges(3,1) * t65 - rSges(3,2) * t63) + g(2) * (-rSges(3,3) - qJ(2))) * t66) - m(4) * (g(1) * (rSges(4,1) * t49 - rSges(4,2) * t48 + rSges(4,3) * t85 + t75) + g(2) * (rSges(4,1) * t47 - rSges(4,2) * t46 + rSges(4,3) * t86 + t74) + g(3) * ((-rSges(4,3) - pkin(5)) * t65 + (rSges(4,1) * t70 - rSges(4,2) * t68) * t63 + t76)) - m(5) * (g(1) * (t49 * pkin(3) + (t49 * t69 + t66 * t84) * rSges(5,1) + (-t49 * t67 + t66 * t83) * rSges(5,2) + t87 * t48 + t75) + g(2) * (t47 * pkin(3) + (t47 * t69 + t64 * t84) * rSges(5,1) + (-t47 * t67 + t64 * t83) * rSges(5,2) + t87 * t46 + t74) + g(3) * ((-t67 * rSges(5,1) - t69 * rSges(5,2) - pkin(5)) * t65 + (t87 * t68 + (t69 * rSges(5,1) - t67 * rSges(5,2) + pkin(3)) * t70) * t63 + t76)) - m(6) * (g(1) * (t88 * t48 + t75) + g(2) * (t88 * t46 + t74) + (g(1) * t66 + g(2) * t64) * t72 * t63 + (g(1) * t49 + g(2) * t47) * t73 + (t76 + (-pkin(5) - t72) * t65 + (t88 * t68 + t73 * t70) * t63) * g(3));
U = t1;

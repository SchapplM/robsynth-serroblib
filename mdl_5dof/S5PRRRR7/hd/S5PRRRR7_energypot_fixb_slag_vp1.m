% Calculate potential energy for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR7_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:42
% EndTime: 2019-12-05 17:11:42
% DurationCPUTime: 0.47s
% Computational Cost: add. (152->90), mult. (180->118), div. (0->0), fcn. (168->10), ass. (0->33)
t82 = rSges(4,3) + pkin(6);
t62 = -pkin(7) - pkin(6);
t81 = rSges(5,3) - t62;
t80 = rSges(6,3) + pkin(8) - t62;
t56 = sin(pkin(9));
t57 = cos(pkin(9));
t79 = g(1) * t57 + g(2) * t56;
t60 = cos(qJ(3));
t46 = t60 * pkin(3) + pkin(2);
t59 = sin(qJ(2));
t75 = rSges(3,2) * t59;
t58 = sin(qJ(3));
t74 = t56 * t58;
t61 = cos(qJ(2));
t73 = t56 * t61;
t72 = t57 * t58;
t71 = t57 * t61;
t70 = t58 * t61;
t69 = t60 * t61;
t55 = qJ(3) + qJ(4);
t48 = cos(t55);
t42 = pkin(4) * t48 + t46;
t68 = t61 * t42;
t67 = t61 * t46;
t64 = t57 * pkin(1) + t56 * pkin(5);
t50 = t56 * pkin(1);
t63 = -t57 * pkin(5) + t50;
t52 = qJ(5) + t55;
t47 = sin(t55);
t45 = cos(t52);
t44 = sin(t52);
t43 = t58 * pkin(3) + pkin(4) * t47;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t57 - rSges(2,2) * t56) + g(2) * (rSges(2,1) * t56 + rSges(2,2) * t57) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,3) * t56 + t64) + g(2) * (rSges(3,1) * t73 - t56 * t75 + t50) + g(3) * (rSges(3,1) * t59 + rSges(3,2) * t61 + qJ(1)) + (g(1) * (rSges(3,1) * t61 - t75) + g(2) * (-rSges(3,3) - pkin(5))) * t57) - m(4) * (g(1) * (pkin(2) * t71 + (t57 * t69 + t74) * rSges(4,1) + (t56 * t60 - t57 * t70) * rSges(4,2) + t64) + g(2) * (pkin(2) * t73 + (t56 * t69 - t72) * rSges(4,1) + (-t56 * t70 - t57 * t60) * rSges(4,2) + t63) + g(3) * (-t82 * t61 + qJ(1)) + (g(3) * (rSges(4,1) * t60 - rSges(4,2) * t58 + pkin(2)) + t79 * t82) * t59) - m(5) * (g(1) * (t57 * t67 + pkin(3) * t74 + (t47 * t56 + t48 * t71) * rSges(5,1) + (-t47 * t71 + t48 * t56) * rSges(5,2) + t64) + g(2) * (t56 * t67 - pkin(3) * t72 + (-t47 * t57 + t48 * t73) * rSges(5,1) + (-t47 * t73 - t48 * t57) * rSges(5,2) + t63) + g(3) * (-t81 * t61 + qJ(1)) + (g(3) * (rSges(5,1) * t48 - rSges(5,2) * t47 + t46) + t79 * t81) * t59) - m(6) * (g(1) * (t57 * t68 + t56 * t43 + (t44 * t56 + t45 * t71) * rSges(6,1) + (-t44 * t71 + t45 * t56) * rSges(6,2) + t64) + g(2) * (t56 * t68 - t57 * t43 + (-t44 * t57 + t45 * t73) * rSges(6,1) + (-t44 * t73 - t45 * t57) * rSges(6,2) + t63) + g(3) * (-t80 * t61 + qJ(1)) + (g(3) * (rSges(6,1) * t45 - rSges(6,2) * t44 + t42) + t79 * t80) * t59);
U = t1;

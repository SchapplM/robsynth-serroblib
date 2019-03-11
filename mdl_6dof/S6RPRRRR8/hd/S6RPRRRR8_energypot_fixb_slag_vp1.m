% Calculate potential energy for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPRRRR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR8_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:40
% EndTime: 2019-03-09 07:19:41
% DurationCPUTime: 0.38s
% Computational Cost: add. (173->83), mult. (180->93), div. (0->0), fcn. (160->10), ass. (0->35)
t55 = qJ(3) + qJ(4);
t49 = cos(t55);
t82 = rSges(7,3) + pkin(10) + pkin(9);
t84 = t82 * t49;
t83 = rSges(6,3) + pkin(9);
t81 = pkin(2) + pkin(6);
t47 = sin(t55);
t80 = g(1) * t47;
t79 = g(2) * t47;
t58 = sin(qJ(1));
t51 = t58 * pkin(1);
t78 = g(2) * t51;
t57 = sin(qJ(3));
t77 = t57 * pkin(3);
t76 = rSges(4,3) + pkin(7);
t63 = -pkin(8) - pkin(7);
t74 = rSges(5,3) - t63;
t61 = cos(qJ(1));
t73 = t61 * pkin(1) + t58 * qJ(2);
t60 = cos(qJ(3));
t72 = t60 * pkin(3) + t81;
t71 = -qJ(2) - t77;
t70 = rSges(4,1) * t57 + rSges(4,2) * t60;
t69 = rSges(5,1) * t47 + rSges(5,2) * t49;
t56 = sin(qJ(5));
t59 = cos(qJ(5));
t68 = rSges(6,1) * t59 - rSges(6,2) * t56 + pkin(4);
t67 = t56 * rSges(6,1) + t59 * rSges(6,2) - t63;
t54 = qJ(5) + qJ(6);
t46 = sin(t54);
t48 = cos(t54);
t66 = rSges(7,1) * t48 - rSges(7,2) * t46 + pkin(5) * t59 + pkin(4);
t65 = g(1) * (t58 * t77 + t73) + t78;
t64 = t46 * rSges(7,1) + t48 * rSges(7,2) + t56 * pkin(5) - t63;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t61 - rSges(2,2) * t58) + g(2) * (rSges(2,1) * t58 + rSges(2,2) * t61) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t61 + rSges(3,3) * t58 + t73) + g(2) * (-rSges(3,2) * t58 + t51 + (-rSges(3,3) - qJ(2)) * t61) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * t73 + t78 + g(3) * (rSges(4,1) * t60 - t57 * rSges(4,2) + t81) + (g(1) * t70 + g(2) * t76) * t58 + (g(1) * t76 + g(2) * (-qJ(2) - t70)) * t61) - m(5) * (g(3) * (rSges(5,1) * t49 - rSges(5,2) * t47 + t72) + (g(1) * t69 + g(2) * t74) * t58 + (g(1) * t74 + g(2) * (-t69 + t71)) * t61 + t65) - m(6) * (g(3) * (t83 * t47 + t72) + (g(2) * t67 + t68 * t80) * t58 + (g(1) * t67 + g(2) * t71 - t68 * t79) * t61 + (g(3) * t68 + (-g(1) * t58 + g(2) * t61) * t83) * t49 + t65) - m(7) * ((-g(1) * t84 + g(2) * t64 + t66 * t80) * t58 + (g(1) * t64 + g(2) * (t71 + t84) - t66 * t79) * t61 + t65 + (t82 * t47 + t66 * t49 + t72) * g(3));
U  = t1;

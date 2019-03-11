% Calculate potential energy for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_energypot_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR2_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:27
% EndTime: 2019-03-09 01:31:27
% DurationCPUTime: 0.32s
% Computational Cost: add. (189->76), mult. (141->82), div. (0->0), fcn. (117->10), ass. (0->33)
t79 = rSges(7,3) + pkin(8);
t53 = pkin(10) + qJ(5);
t46 = sin(t53);
t58 = sin(qJ(6));
t60 = cos(qJ(6));
t64 = rSges(7,1) * t60 - rSges(7,2) * t58 + pkin(5);
t78 = t46 * t64;
t55 = sin(pkin(10));
t77 = pkin(4) * t55;
t57 = -pkin(7) - qJ(4);
t75 = rSges(6,3) - t57;
t74 = pkin(6) + qJ(2);
t54 = qJ(1) + pkin(9);
t47 = sin(t54);
t59 = sin(qJ(1));
t51 = t59 * pkin(1);
t73 = t47 * pkin(2) + t51;
t72 = rSges(5,3) + qJ(4);
t71 = pkin(3) + t74;
t49 = cos(t54);
t61 = cos(qJ(1));
t52 = t61 * pkin(1);
t70 = t49 * pkin(2) + t47 * qJ(3) + t52;
t69 = g(2) * t73;
t68 = -qJ(3) - t77;
t56 = cos(pkin(10));
t67 = t56 * pkin(4) + t71;
t66 = rSges(5,1) * t55 + rSges(5,2) * t56;
t48 = cos(t53);
t65 = rSges(6,1) * t46 + rSges(6,2) * t48;
t63 = rSges(7,1) * t58 + rSges(7,2) * t60 - t57;
t62 = g(1) * (t47 * t77 + t70) + t69;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t61 - t59 * rSges(2,2)) + g(2) * (t59 * rSges(2,1) + rSges(2,2) * t61) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t49 - rSges(3,2) * t47 + t52) + g(2) * (rSges(3,1) * t47 + rSges(3,2) * t49 + t51) + g(3) * (rSges(3,3) + t74)) - m(4) * (g(1) * (-rSges(4,2) * t49 + rSges(4,3) * t47 + t70) + g(2) * (-rSges(4,2) * t47 + (-rSges(4,3) - qJ(3)) * t49 + t73) + g(3) * (rSges(4,1) + t74)) - m(5) * (g(1) * t70 + t69 + g(3) * (rSges(5,1) * t56 - rSges(5,2) * t55 + t71) + (g(1) * t66 + g(2) * t72) * t47 + (g(1) * t72 + g(2) * (-qJ(3) - t66)) * t49) - m(6) * (g(3) * (rSges(6,1) * t48 - rSges(6,2) * t46 + t67) + (g(1) * t65 + g(2) * t75) * t47 + (g(1) * t75 + g(2) * (-t65 + t68)) * t49 + t62) - m(7) * (g(3) * (t79 * t46 + t67) + (g(1) * t78 + g(2) * t63) * t47 + (g(1) * t63 + (t68 - t78) * g(2)) * t49 + (g(3) * t64 + (-g(1) * t47 + g(2) * t49) * t79) * t48 + t62);
U  = t1;

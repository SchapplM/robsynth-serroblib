% Calculate potential energy for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:27
% EndTime: 2019-03-09 02:30:27
% DurationCPUTime: 0.37s
% Computational Cost: add. (128->89), mult. (163->105), div. (0->0), fcn. (143->8), ass. (0->30)
t80 = rSges(6,3) + pkin(8);
t79 = rSges(7,3) + pkin(9) + pkin(8);
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t78 = g(1) * t57 + g(2) * t54;
t77 = pkin(2) + pkin(6);
t52 = sin(qJ(5));
t73 = t54 * t52;
t53 = sin(qJ(4));
t72 = t54 * t53;
t55 = cos(qJ(5));
t71 = t54 * t55;
t70 = t57 * t52;
t69 = t57 * t53;
t68 = t57 * t55;
t48 = t54 * pkin(1);
t66 = t54 * qJ(3) + t48;
t65 = t57 * pkin(1) + t54 * qJ(2);
t64 = pkin(3) + t77;
t63 = t57 * pkin(7) + t66;
t62 = t57 * qJ(3) + t65;
t56 = cos(qJ(4));
t61 = rSges(5,1) * t53 + rSges(5,2) * t56;
t60 = -t54 * pkin(7) + t62;
t59 = -t57 * qJ(2) + t63;
t51 = qJ(5) + qJ(6);
t44 = cos(t51);
t43 = sin(t51);
t42 = t55 * pkin(5) + pkin(4);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t57 * rSges(2,1) - t54 * rSges(2,2)) + g(2) * (t54 * rSges(2,1) + t57 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-t57 * rSges(3,2) + t54 * rSges(3,3) + t65) + g(2) * (-t54 * rSges(3,2) + t48 + (-rSges(3,3) - qJ(2)) * t57) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (t54 * rSges(4,2) + t57 * rSges(4,3) + t62) + g(2) * (t54 * rSges(4,3) + (-rSges(4,2) - qJ(2)) * t57 + t66) + g(3) * (rSges(4,1) + t77)) - m(5) * (g(1) * t62 + g(2) * t63 + g(3) * (t56 * rSges(5,1) - t53 * rSges(5,2) + t64) + (g(1) * t61 + g(2) * (rSges(5,3) - qJ(2))) * t57 + (g(1) * (-rSges(5,3) - pkin(7)) + g(2) * t61) * t54) - m(6) * (g(1) * (pkin(4) * t69 + (t53 * t68 - t73) * rSges(6,1) + (-t52 * t69 - t71) * rSges(6,2) + t60) + g(2) * (pkin(4) * t72 + (t53 * t71 + t70) * rSges(6,1) + (-t52 * t72 + t68) * rSges(6,2) + t59) + g(3) * (t80 * t53 + t64) + (g(3) * (rSges(6,1) * t55 - rSges(6,2) * t52 + pkin(4)) - t78 * t80) * t56) - m(7) * (g(1) * (t42 * t69 - pkin(5) * t73 + (-t54 * t43 + t44 * t69) * rSges(7,1) + (-t43 * t69 - t54 * t44) * rSges(7,2) + t60) + g(2) * (t42 * t72 + pkin(5) * t70 + (t57 * t43 + t44 * t72) * rSges(7,1) + (-t43 * t72 + t57 * t44) * rSges(7,2) + t59) + g(3) * (t79 * t53 + t64) + (g(3) * (rSges(7,1) * t44 - rSges(7,2) * t43 + t42) - t78 * t79) * t56);
U  = t1;

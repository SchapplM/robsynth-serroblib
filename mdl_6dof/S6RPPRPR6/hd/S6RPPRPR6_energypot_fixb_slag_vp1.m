% Calculate potential energy for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRPR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR6_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:32
% EndTime: 2019-03-09 01:50:33
% DurationCPUTime: 0.33s
% Computational Cost: add. (110->84), mult. (150->96), div. (0->0), fcn. (126->6), ass. (0->25)
t73 = rSges(7,3) + pkin(8);
t72 = pkin(2) + pkin(6);
t51 = sin(qJ(4));
t71 = pkin(4) * t51;
t52 = sin(qJ(1));
t54 = cos(qJ(4));
t69 = t52 * t54;
t50 = sin(qJ(6));
t55 = cos(qJ(1));
t68 = t55 * t50;
t53 = cos(qJ(6));
t67 = t55 * t53;
t46 = t52 * pkin(1);
t66 = t52 * qJ(3) + t46;
t65 = t55 * pkin(1) + t52 * qJ(2);
t64 = t54 * qJ(5);
t63 = pkin(3) + t72;
t62 = t55 * pkin(7) + t66;
t61 = t55 * qJ(3) + t65;
t60 = t52 * t71 + t62;
t59 = t55 * t71 + t61;
t58 = t54 * pkin(4) + t51 * qJ(5) + t63;
t57 = rSges(5,1) * t51 + rSges(5,2) * t54;
t56 = -rSges(6,2) * t51 - rSges(6,3) * t54 - t64;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t55 * rSges(2,1) - t52 * rSges(2,2)) + g(2) * (t52 * rSges(2,1) + t55 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-t55 * rSges(3,2) + t52 * rSges(3,3) + t65) + g(2) * (-t52 * rSges(3,2) + t46 + (-rSges(3,3) - qJ(2)) * t55) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (t52 * rSges(4,2) + t55 * rSges(4,3) + t61) + g(2) * (t52 * rSges(4,3) + (-rSges(4,2) - qJ(2)) * t55 + t66) + g(3) * (rSges(4,1) + t72)) - m(5) * (g(1) * t61 + g(2) * t62 + g(3) * (t54 * rSges(5,1) - t51 * rSges(5,2) + t63) + (g(1) * t57 + g(2) * (rSges(5,3) - qJ(2))) * t55 + (g(1) * (-rSges(5,3) - pkin(7)) + g(2) * t57) * t52) - m(6) * (g(1) * t59 + g(2) * t60 + g(3) * (-t54 * rSges(6,2) + t51 * rSges(6,3) + t58) + (g(1) * t56 + g(2) * (rSges(6,1) - qJ(2))) * t55 + (g(1) * (-rSges(6,1) - pkin(7)) + g(2) * t56) * t52) - m(7) * (g(1) * (-t55 * t64 + t59 + (-t68 * rSges(7,1) - t67 * rSges(7,2)) * t54 + (-t53 * rSges(7,1) + t50 * rSges(7,2) - pkin(5) - pkin(7)) * t52) + g(2) * (-t52 * t64 + (-t50 * t69 + t67) * rSges(7,1) + (-t53 * t69 - t68) * rSges(7,2) + t60 + (pkin(5) - qJ(2)) * t55) + g(3) * (t73 * t54 + t58) + (g(3) * (rSges(7,1) * t50 + rSges(7,2) * t53) + (g(1) * t55 + g(2) * t52) * t73) * t51);
U  = t1;

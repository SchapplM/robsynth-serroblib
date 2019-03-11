% Calculate potential energy for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPPRR5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energypot_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_energypot_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR5_energypot_fixb_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:09
% EndTime: 2019-03-09 01:37:10
% DurationCPUTime: 0.33s
% Computational Cost: add. (136->79), mult. (197->94), div. (0->0), fcn. (205->8), ass. (0->27)
t72 = rSges(7,3) + pkin(8);
t54 = cos(qJ(1));
t66 = sin(qJ(1));
t47 = t66 * pkin(1);
t63 = t66 * qJ(3) + t47;
t71 = t63 + (-pkin(3) - qJ(2)) * t54;
t70 = pkin(2) + pkin(6);
t53 = cos(qJ(5));
t69 = t53 * pkin(5);
t68 = -rSges(6,3) - pkin(7);
t50 = sin(qJ(6));
t65 = t50 * t53;
t52 = cos(qJ(6));
t64 = t52 * t53;
t62 = t54 * pkin(1) + t66 * qJ(2);
t61 = sin(pkin(9));
t60 = qJ(4) + t70;
t59 = t54 * qJ(3) + t62;
t58 = t66 * pkin(3) + t59;
t49 = cos(pkin(9));
t41 = t66 * t49 + t54 * t61;
t57 = t41 * pkin(4) + t58;
t51 = sin(qJ(5));
t56 = rSges(6,1) * t53 - rSges(6,2) * t51;
t40 = t49 * t54 - t66 * t61;
t55 = -t40 * pkin(4) + t71;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t54 * rSges(2,1) - t66 * rSges(2,2)) + g(2) * (t66 * rSges(2,1) + t54 * rSges(2,2)) + g(3) * (pkin(6) + rSges(2,3))) - m(3) * (g(1) * (-t54 * rSges(3,2) + t66 * rSges(3,3) + t62) + g(2) * (-t66 * rSges(3,2) + t47 + (-rSges(3,3) - qJ(2)) * t54) + g(3) * (pkin(6) + rSges(3,1))) - m(4) * (g(1) * (t66 * rSges(4,1) + t54 * rSges(4,3) + t59) + g(2) * (t66 * rSges(4,3) + (-rSges(4,1) - qJ(2)) * t54 + t63) + g(3) * (-rSges(4,2) + t70)) - m(5) * (g(1) * (rSges(5,1) * t41 + rSges(5,2) * t40 + t58) + g(2) * (-t40 * rSges(5,1) + t41 * rSges(5,2) + t71) + g(3) * (rSges(5,3) + t60)) - m(6) * (g(1) * t57 + g(2) * t55 + g(3) * (rSges(6,1) * t51 + rSges(6,2) * t53 + t60) + (g(1) * t56 + g(2) * t68) * t41 + (g(1) * t68 - g(2) * t56) * t40) - m(7) * (g(1) * (t41 * t69 - t40 * pkin(7) + (-t40 * t50 + t41 * t64) * rSges(7,1) + (-t40 * t52 - t41 * t65) * rSges(7,2) + t57) + g(2) * (-t40 * t69 - t41 * pkin(7) + (-t40 * t64 - t41 * t50) * rSges(7,1) + (t40 * t65 - t41 * t52) * rSges(7,2) + t55) + g(3) * (-t72 * t53 + t60) + (g(3) * (rSges(7,1) * t52 - rSges(7,2) * t50 + pkin(5)) + (g(1) * t41 - g(2) * t40) * t72) * t51);
U  = t1;

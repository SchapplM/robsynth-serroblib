% Calculate potential energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(5,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energypot_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:04
% EndTime: 2019-07-18 17:17:04
% DurationCPUTime: 0.22s
% Computational Cost: add. (126->72), mult. (144->99), div. (0->0), fcn. (128->10), ass. (0->31)
t53 = cos(qJ(2));
t69 = pkin(1) * t53;
t48 = qJ(2) + qJ(3);
t43 = sin(t48);
t68 = pkin(5) * t43;
t50 = sin(qJ(2));
t67 = t50 * pkin(1) + pkin(4);
t45 = cos(t48);
t51 = sin(qJ(1));
t66 = t51 * t45;
t49 = sin(qJ(4));
t65 = t51 * t49;
t52 = cos(qJ(4));
t64 = t51 * t52;
t54 = cos(qJ(1));
t63 = t54 * t45;
t62 = t54 * t49;
t61 = t54 * t52;
t39 = t51 * t69;
t60 = t51 * t68 + t39;
t40 = t54 * t69;
t59 = t54 * t68 + t40;
t58 = -t45 * pkin(5) + t67;
t57 = g(1) * t54 + g(2) * t51;
t56 = rSges(3,1) * t53 - rSges(3,2) * t50;
t55 = rSges(4,1) * t45 - rSges(4,2) * t43;
t47 = qJ(4) + qJ(5);
t44 = cos(t47);
t42 = sin(t47);
t41 = t52 * pkin(3) + pkin(2);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t54 * rSges(2,1) - t51 * rSges(2,2)) + g(2) * (t51 * rSges(2,1) + t54 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(1) * (t51 * rSges(3,3) + t56 * t54) + g(2) * (-t54 * rSges(3,3) + t56 * t51) + g(3) * (t50 * rSges(3,1) + t53 * rSges(3,2) + pkin(4))) - m(4) * (g(1) * (t51 * rSges(4,3) + t55 * t54 + t40) + g(2) * (-t54 * rSges(4,3) + t55 * t51 + t39) + g(3) * (t43 * rSges(4,1) + t45 * rSges(4,2) + t67)) - m(5) * (g(1) * (pkin(2) * t63 + (t45 * t61 + t65) * rSges(5,1) + (-t45 * t62 + t64) * rSges(5,2) + t59) + g(2) * (pkin(2) * t66 + (t45 * t64 - t62) * rSges(5,1) + (-t45 * t65 - t61) * rSges(5,2) + t60) + g(3) * (-t45 * rSges(5,3) + t58) + (g(3) * (rSges(5,1) * t52 - rSges(5,2) * t49 + pkin(2)) + t57 * rSges(5,3)) * t43) - m(6) * (g(1) * (t41 * t63 + pkin(3) * t65 + (t51 * t42 + t44 * t63) * rSges(6,1) + (-t42 * t63 + t51 * t44) * rSges(6,2) + t59) + g(2) * (t41 * t66 - pkin(3) * t62 + (-t54 * t42 + t44 * t66) * rSges(6,1) + (-t42 * t66 - t54 * t44) * rSges(6,2) + t60) + g(3) * (-t45 * rSges(6,3) + t58) + (g(3) * (rSges(6,1) * t44 - rSges(6,2) * t42 + t41) + t57 * rSges(6,3)) * t43);
U  = t1;

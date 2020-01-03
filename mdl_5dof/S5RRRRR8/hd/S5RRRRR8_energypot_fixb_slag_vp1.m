% Calculate potential energy for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR8_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR8_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:18
% EndTime: 2019-12-31 22:24:19
% DurationCPUTime: 0.38s
% Computational Cost: add. (151->79), mult. (154->99), div. (0->0), fcn. (138->10), ass. (0->31)
t79 = rSges(5,3) + pkin(8);
t78 = rSges(6,3) + pkin(9) + pkin(8);
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t77 = g(1) * t59 + g(2) * t56;
t74 = rSges(3,3) + pkin(6);
t55 = sin(qJ(2));
t72 = t55 * pkin(2) + pkin(5);
t53 = qJ(2) + qJ(3);
t48 = sin(t53);
t71 = rSges(4,2) * t48;
t54 = sin(qJ(4));
t70 = t54 * t56;
t69 = t54 * t59;
t50 = cos(t53);
t68 = t56 * t50;
t57 = cos(qJ(4));
t67 = t56 * t57;
t66 = t59 * t50;
t58 = cos(qJ(2));
t45 = pkin(2) * t58 + pkin(1);
t61 = -pkin(7) - pkin(6);
t64 = t56 * t45 + t59 * t61;
t43 = t59 * t45;
t63 = -t56 * t61 + t43;
t62 = rSges(3,1) * t58 - rSges(3,2) * t55 + pkin(1);
t52 = qJ(4) + qJ(5);
t49 = cos(t52);
t47 = sin(t52);
t44 = pkin(4) * t57 + pkin(3);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t59 - rSges(2,2) * t56) + g(2) * (rSges(2,1) * t56 + rSges(2,2) * t59) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t55 + rSges(3,2) * t58 + pkin(5)) + (g(1) * t62 - g(2) * t74) * t59 + (g(1) * t74 + g(2) * t62) * t56) - m(4) * (g(1) * (rSges(4,1) * t66 - t59 * t71 + t43) + g(2) * (-rSges(4,3) * t59 + t64) + g(3) * (rSges(4,1) * t48 + rSges(4,2) * t50 + t72) + (g(1) * (rSges(4,3) - t61) + g(2) * (rSges(4,1) * t50 - t71)) * t56) - m(5) * (g(1) * (pkin(3) * t66 + (t57 * t66 + t70) * rSges(5,1) + (-t54 * t66 + t67) * rSges(5,2) + t63) + g(2) * (pkin(3) * t68 + (t50 * t67 - t69) * rSges(5,1) + (-t54 * t68 - t57 * t59) * rSges(5,2) + t64) + g(3) * (-t79 * t50 + t72) + (g(3) * (rSges(5,1) * t57 - rSges(5,2) * t54 + pkin(3)) + t77 * t79) * t48) - m(6) * (g(1) * (t44 * t66 + pkin(4) * t70 + (t47 * t56 + t49 * t66) * rSges(6,1) + (-t47 * t66 + t49 * t56) * rSges(6,2) + t63) + g(2) * (t44 * t68 - pkin(4) * t69 + (-t47 * t59 + t49 * t68) * rSges(6,1) + (-t47 * t68 - t49 * t59) * rSges(6,2) + t64) + g(3) * (-t78 * t50 + t72) + (g(3) * (rSges(6,1) * t49 - rSges(6,2) * t47 + t44) + t77 * t78) * t48);
U = t1;

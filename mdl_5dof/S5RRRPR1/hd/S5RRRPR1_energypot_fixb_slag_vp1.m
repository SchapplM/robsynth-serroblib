% Calculate potential energy for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRPR1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:44
% EndTime: 2019-12-05 18:37:44
% DurationCPUTime: 0.18s
% Computational Cost: add. (138->59), mult. (110->66), div. (0->0), fcn. (86->10), ass. (0->28)
t58 = -pkin(7) - pkin(6);
t68 = rSges(3,3) + pkin(6);
t54 = sin(qJ(2));
t67 = t54 * pkin(2) + pkin(5);
t56 = cos(qJ(2));
t45 = t56 * pkin(2) + pkin(1);
t66 = rSges(4,3) - t58;
t52 = -qJ(4) + t58;
t65 = rSges(5,3) - t52;
t64 = rSges(6,3) + pkin(8) - t52;
t53 = qJ(2) + qJ(3);
t47 = sin(t53);
t63 = pkin(3) * t47 + t67;
t48 = cos(t53);
t37 = pkin(3) * t48 + t45;
t46 = pkin(9) + t53;
t62 = rSges(3,1) * t56 - rSges(3,2) * t54 + pkin(1);
t61 = rSges(4,1) * t48 - rSges(4,2) * t47 + t45;
t40 = sin(t46);
t41 = cos(t46);
t60 = rSges(5,1) * t41 - rSges(5,2) * t40 + t37;
t44 = qJ(5) + t46;
t38 = sin(t44);
t39 = cos(t44);
t59 = rSges(6,1) * t39 - rSges(6,2) * t38 + pkin(4) * t41 + t37;
t57 = cos(qJ(1));
t55 = sin(qJ(1));
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t57 - rSges(2,2) * t55) + g(2) * (rSges(2,1) * t55 + rSges(2,2) * t57) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(3) * (rSges(3,1) * t54 + rSges(3,2) * t56 + pkin(5)) + (g(1) * t62 - g(2) * t68) * t57 + (g(1) * t68 + g(2) * t62) * t55) - m(4) * (g(3) * (rSges(4,1) * t47 + rSges(4,2) * t48 + t67) + (g(1) * t61 - g(2) * t66) * t57 + (g(1) * t66 + g(2) * t61) * t55) - m(5) * (g(3) * (rSges(5,1) * t40 + rSges(5,2) * t41 + t63) + (g(1) * t60 - g(2) * t65) * t57 + (g(1) * t65 + g(2) * t60) * t55) - m(6) * (g(3) * (rSges(6,1) * t38 + rSges(6,2) * t39 + pkin(4) * t40 + t63) + (g(1) * t59 - g(2) * t64) * t57 + (g(1) * t64 + g(2) * t59) * t55);
U = t1;

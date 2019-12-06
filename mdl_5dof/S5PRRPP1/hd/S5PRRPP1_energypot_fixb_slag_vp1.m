% Calculate potential energy for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPP1_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP1_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:07
% EndTime: 2019-12-05 16:06:08
% DurationCPUTime: 0.23s
% Computational Cost: add. (137->61), mult. (104->68), div. (0->0), fcn. (80->8), ass. (0->25)
t63 = rSges(6,1) + pkin(4);
t46 = qJ(3) + pkin(8);
t39 = sin(t46);
t41 = cos(t46);
t62 = rSges(5,1) * t41 - rSges(5,2) * t39;
t61 = rSges(6,3) + qJ(5);
t60 = rSges(4,3) + pkin(6);
t57 = pkin(5) + qJ(1);
t51 = cos(qJ(3));
t37 = t51 * pkin(3) + pkin(2);
t45 = pkin(7) + qJ(2);
t40 = cos(t45);
t48 = cos(pkin(7));
t43 = t48 * pkin(1);
t56 = t40 * t37 + t43;
t50 = sin(qJ(3));
t55 = t50 * pkin(3) + t57;
t38 = sin(t45);
t47 = sin(pkin(7));
t42 = t47 * pkin(1);
t49 = -qJ(4) - pkin(6);
t54 = t38 * t37 + t40 * t49 + t42;
t53 = rSges(4,1) * t51 - rSges(4,2) * t50 + pkin(2);
t52 = t61 * t39 + t63 * t41;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t48 * rSges(2,1) - t47 * rSges(2,2)) + g(2) * (t47 * rSges(2,1) + t48 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t40 * rSges(3,1) - t38 * rSges(3,2) + t43) + g(2) * (t38 * rSges(3,1) + t40 * rSges(3,2) + t42) + g(3) * (rSges(3,3) + t57)) - m(4) * (g(1) * t43 + g(2) * t42 + g(3) * (t50 * rSges(4,1) + t51 * rSges(4,2) + t57) + (g(1) * t53 - g(2) * t60) * t40 + (g(1) * t60 + g(2) * t53) * t38) - m(5) * (g(1) * (t62 * t40 + t56) + g(2) * (-t40 * rSges(5,3) + t54) + g(3) * (t39 * rSges(5,1) + t41 * rSges(5,2) + t55) + (g(1) * (rSges(5,3) - t49) + g(2) * t62) * t38) - m(6) * (g(1) * t56 + g(2) * t54 + g(3) * (t63 * t39 - t61 * t41 + t55) + (-g(2) * rSges(6,2) + g(1) * t52) * t40 + (g(1) * (rSges(6,2) - t49) + g(2) * t52) * t38);
U = t1;

% Calculate potential energy for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_energypot_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:52:52
% EndTime: 2019-12-05 18:52:52
% DurationCPUTime: 0.15s
% Computational Cost: add. (100->53), mult. (99->73), div. (0->0), fcn. (79->10), ass. (0->26)
t46 = cos(qJ(3));
t57 = pkin(2) * t46;
t40 = qJ(3) + qJ(4);
t33 = sin(t40);
t56 = t33 * rSges(6,3);
t41 = qJ(1) + qJ(2);
t34 = sin(t41);
t42 = sin(qJ(5));
t55 = t34 * t42;
t45 = cos(qJ(5));
t54 = t34 * t45;
t36 = cos(t41);
t53 = t36 * t42;
t52 = t36 * t45;
t44 = sin(qJ(1));
t38 = t44 * pkin(1);
t51 = t34 * t57 + t38;
t47 = cos(qJ(1));
t39 = t47 * pkin(1);
t50 = t36 * t57 + t39;
t43 = sin(qJ(3));
t49 = rSges(4,1) * t46 - rSges(4,2) * t43;
t35 = cos(t40);
t48 = rSges(5,1) * t35 - rSges(5,2) * t33;
t37 = t43 * pkin(2);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t47 * rSges(2,1) - t44 * rSges(2,2)) + g(2) * (t44 * rSges(2,1) + t47 * rSges(2,2)) + g(3) * rSges(2,3)) - m(3) * (g(1) * (t36 * rSges(3,1) - t34 * rSges(3,2) + t39) + g(2) * (t34 * rSges(3,1) + t36 * rSges(3,2) + t38) + g(3) * rSges(3,3)) - m(4) * (g(1) * (t34 * rSges(4,3) + t36 * t49 + t39) + g(2) * (-t36 * rSges(4,3) + t34 * t49 + t38) + g(3) * (t43 * rSges(4,1) + t46 * rSges(4,2))) - m(5) * (g(1) * (t34 * rSges(5,3) + t36 * t48 + t50) + g(2) * (-t36 * rSges(5,3) + t34 * t48 + t51) + g(3) * (t33 * rSges(5,1) + t35 * rSges(5,2) + t37)) - m(6) * (g(1) * ((t35 * t52 + t55) * rSges(6,1) + (-t35 * t53 + t54) * rSges(6,2) + t36 * t56 + t50) + g(2) * ((t35 * t54 - t53) * rSges(6,1) + (-t35 * t55 - t52) * rSges(6,2) + t34 * t56 + t51) + g(3) * (-t35 * rSges(6,3) + t37 + (rSges(6,1) * t45 - rSges(6,2) * t42) * t33));
U = t1;

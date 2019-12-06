% Calculate potential energy for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:18:59
% EndTime: 2019-12-05 16:18:59
% DurationCPUTime: 0.18s
% Computational Cost: add. (135->59), mult. (97->64), div. (0->0), fcn. (73->10), ass. (0->27)
t63 = rSges(4,3) + pkin(6);
t54 = cos(qJ(3));
t37 = t54 * pkin(3) + pkin(2);
t52 = -qJ(4) - pkin(6);
t62 = rSges(5,3) - t52;
t61 = rSges(6,3) + pkin(7) - t52;
t60 = pkin(5) + qJ(1);
t49 = qJ(3) + pkin(9);
t53 = sin(qJ(3));
t59 = t53 * pkin(3) + t60;
t50 = sin(pkin(8));
t43 = t50 * pkin(1);
t51 = cos(pkin(8));
t44 = t51 * pkin(1);
t58 = g(1) * t44 + g(2) * t43;
t57 = rSges(4,1) * t54 - rSges(4,2) * t53 + pkin(2);
t39 = sin(t49);
t41 = cos(t49);
t56 = rSges(5,1) * t41 - rSges(5,2) * t39 + t37;
t42 = qJ(5) + t49;
t35 = sin(t42);
t36 = cos(t42);
t55 = rSges(6,1) * t36 - rSges(6,2) * t35 + pkin(4) * t41 + t37;
t48 = pkin(8) + qJ(2);
t40 = cos(t48);
t38 = sin(t48);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t51 - rSges(2,2) * t50) + g(2) * (rSges(2,1) * t50 + rSges(2,2) * t51) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (rSges(3,1) * t40 - rSges(3,2) * t38 + t44) + g(2) * (rSges(3,1) * t38 + rSges(3,2) * t40 + t43) + g(3) * (rSges(3,3) + t60)) - m(4) * (g(3) * (t53 * rSges(4,1) + t54 * rSges(4,2) + t60) + (g(1) * t57 - g(2) * t63) * t40 + (g(1) * t63 + g(2) * t57) * t38 + t58) - m(5) * (g(3) * (rSges(5,1) * t39 + rSges(5,2) * t41 + t59) + (g(1) * t56 - g(2) * t62) * t40 + (g(1) * t62 + g(2) * t56) * t38 + t58) - m(6) * (g(3) * (rSges(6,1) * t35 + rSges(6,2) * t36 + pkin(4) * t39 + t59) + (g(1) * t55 - g(2) * t61) * t40 + (g(1) * t61 + g(2) * t55) * t38 + t58);
U = t1;

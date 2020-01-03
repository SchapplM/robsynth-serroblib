% Calculate potential energy for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRP5_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_energypot_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP5_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:48
% EndTime: 2019-12-31 16:28:48
% DurationCPUTime: 0.26s
% Computational Cost: add. (77->60), mult. (125->79), div. (0->0), fcn. (113->6), ass. (0->25)
t66 = rSges(4,3) + pkin(5);
t65 = rSges(5,3) + qJ(4) + pkin(5);
t44 = sin(pkin(6));
t45 = cos(pkin(6));
t64 = g(1) * t45 + g(2) * t44;
t48 = sin(qJ(2));
t60 = rSges(3,2) * t48;
t47 = sin(qJ(3));
t59 = t44 * t47;
t50 = cos(qJ(2));
t58 = t44 * t50;
t57 = t45 * t47;
t56 = t45 * t50;
t55 = t47 * t50;
t49 = cos(qJ(3));
t54 = t49 * t50;
t52 = t45 * pkin(1) + t44 * pkin(4);
t42 = t44 * pkin(1);
t51 = -t45 * pkin(4) + t42;
t40 = pkin(3) * t49 + pkin(2);
t39 = t45 * t54 + t59;
t38 = t44 * t49 - t45 * t55;
t37 = t44 * t54 - t57;
t36 = -t44 * t55 - t45 * t49;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t45 - rSges(2,2) * t44) + g(2) * (rSges(2,1) * t44 + rSges(2,2) * t45) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t44 * rSges(3,3) + t52) + g(2) * (rSges(3,1) * t58 - t44 * t60 + t42) + g(3) * (t48 * rSges(3,1) + rSges(3,2) * t50 + qJ(1)) + (g(1) * (rSges(3,1) * t50 - t60) + g(2) * (-rSges(3,3) - pkin(4))) * t45) - m(4) * (g(1) * (t39 * rSges(4,1) + t38 * rSges(4,2) + pkin(2) * t56 + t52) + g(2) * (t37 * rSges(4,1) + t36 * rSges(4,2) + pkin(2) * t58 + t51) + g(3) * (-t66 * t50 + qJ(1)) + (g(3) * (rSges(4,1) * t49 - rSges(4,2) * t47 + pkin(2)) + t64 * t66) * t48) - m(5) * (g(1) * (t39 * rSges(5,1) + t38 * rSges(5,2) + pkin(3) * t59 + t40 * t56 + t52) + g(2) * (t37 * rSges(5,1) + t36 * rSges(5,2) - pkin(3) * t57 + t40 * t58 + t51) + g(3) * (-t65 * t50 + qJ(1)) + (g(3) * (rSges(5,1) * t49 - rSges(5,2) * t47 + t40) + t64 * t65) * t48);
U = t1;

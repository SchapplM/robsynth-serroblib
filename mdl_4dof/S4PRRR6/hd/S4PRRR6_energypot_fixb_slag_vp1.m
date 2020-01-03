% Calculate potential energy for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PRRR6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR6_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:36
% EndTime: 2019-12-31 16:34:36
% DurationCPUTime: 0.29s
% Computational Cost: add. (87->65), mult. (125->88), div. (0->0), fcn. (113->8), ass. (0->25)
t65 = rSges(4,3) + pkin(5);
t64 = rSges(5,3) + pkin(6) + pkin(5);
t42 = sin(pkin(7));
t43 = cos(pkin(7));
t63 = g(1) * t43 + g(2) * t42;
t45 = sin(qJ(2));
t59 = rSges(3,2) * t45;
t44 = sin(qJ(3));
t58 = t42 * t44;
t47 = cos(qJ(2));
t57 = t42 * t47;
t56 = t43 * t44;
t55 = t43 * t47;
t54 = t44 * t47;
t46 = cos(qJ(3));
t53 = t46 * t47;
t35 = t46 * pkin(3) + pkin(2);
t52 = t47 * t35;
t50 = t43 * pkin(1) + t42 * pkin(4);
t39 = t42 * pkin(1);
t49 = -t43 * pkin(4) + t39;
t41 = qJ(3) + qJ(4);
t37 = cos(t41);
t36 = sin(t41);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t43 * rSges(2,1) - t42 * rSges(2,2)) + g(2) * (t42 * rSges(2,1) + t43 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(1) * (t42 * rSges(3,3) + t50) + g(2) * (rSges(3,1) * t57 - t42 * t59 + t39) + g(3) * (t45 * rSges(3,1) + t47 * rSges(3,2) + qJ(1)) + (g(1) * (rSges(3,1) * t47 - t59) + g(2) * (-rSges(3,3) - pkin(4))) * t43) - m(4) * (g(1) * (pkin(2) * t55 + (t43 * t53 + t58) * rSges(4,1) + (t42 * t46 - t43 * t54) * rSges(4,2) + t50) + g(2) * (pkin(2) * t57 + (t42 * t53 - t56) * rSges(4,1) + (-t42 * t54 - t43 * t46) * rSges(4,2) + t49) + g(3) * (-t65 * t47 + qJ(1)) + (g(3) * (rSges(4,1) * t46 - rSges(4,2) * t44 + pkin(2)) + t63 * t65) * t45) - m(5) * (g(1) * (t43 * t52 + pkin(3) * t58 + (t42 * t36 + t37 * t55) * rSges(5,1) + (-t36 * t55 + t42 * t37) * rSges(5,2) + t50) + g(2) * (t42 * t52 - pkin(3) * t56 + (-t43 * t36 + t37 * t57) * rSges(5,1) + (-t36 * t57 - t43 * t37) * rSges(5,2) + t49) + g(3) * (-t64 * t47 + qJ(1)) + (g(3) * (rSges(5,1) * t37 - rSges(5,2) * t36 + t35) + t63 * t64) * t45);
U = t1;

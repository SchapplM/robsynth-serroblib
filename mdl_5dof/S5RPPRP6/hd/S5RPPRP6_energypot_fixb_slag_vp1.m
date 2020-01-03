% Calculate potential energy for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRP6_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP6_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:54
% EndTime: 2019-12-31 17:54:54
% DurationCPUTime: 0.20s
% Computational Cost: add. (100->61), mult. (110->65), div. (0->0), fcn. (86->6), ass. (0->25)
t58 = rSges(6,1) + pkin(4);
t57 = rSges(6,3) + qJ(5);
t56 = pkin(2) + pkin(5);
t39 = sin(pkin(7));
t55 = pkin(3) * t39;
t42 = sin(qJ(1));
t36 = t42 * pkin(1);
t54 = g(2) * t36;
t41 = -pkin(6) - qJ(3);
t53 = rSges(6,2) - t41;
t52 = rSges(5,3) - t41;
t43 = cos(qJ(1));
t51 = t43 * pkin(1) + t42 * qJ(2);
t50 = rSges(4,3) + qJ(3);
t40 = cos(pkin(7));
t49 = t40 * pkin(3) + t56;
t48 = -qJ(2) - t55;
t47 = rSges(4,1) * t39 + rSges(4,2) * t40;
t38 = pkin(7) + qJ(4);
t32 = sin(t38);
t33 = cos(t38);
t46 = rSges(5,1) * t32 + rSges(5,2) * t33;
t45 = g(1) * (t42 * t55 + t51) + t54;
t44 = t58 * t32 - t57 * t33;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t43 * rSges(2,1) - t42 * rSges(2,2)) + g(2) * (t42 * rSges(2,1) + t43 * rSges(2,2)) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (-t43 * rSges(3,2) + t42 * rSges(3,3) + t51) + g(2) * (-t42 * rSges(3,2) + t36 + (-rSges(3,3) - qJ(2)) * t43) + g(3) * (pkin(5) + rSges(3,1))) - m(4) * (g(1) * t51 + t54 + g(3) * (t40 * rSges(4,1) - t39 * rSges(4,2) + t56) + (g(1) * t47 + g(2) * t50) * t42 + (g(1) * t50 + g(2) * (-qJ(2) - t47)) * t43) - m(5) * (g(3) * (t33 * rSges(5,1) - t32 * rSges(5,2) + t49) + (g(1) * t46 + g(2) * t52) * t42 + (g(1) * t52 + g(2) * (-t46 + t48)) * t43 + t45) - m(6) * (g(3) * (t57 * t32 + t58 * t33 + t49) + (g(1) * t44 + g(2) * t53) * t42 + (g(1) * t53 + g(2) * (-t44 + t48)) * t43 + t45);
U = t1;

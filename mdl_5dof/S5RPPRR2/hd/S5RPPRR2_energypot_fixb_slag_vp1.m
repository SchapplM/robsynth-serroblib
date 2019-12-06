% Calculate potential energy for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPRR2_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_energypot_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR2_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:30
% EndTime: 2019-12-05 17:39:30
% DurationCPUTime: 0.18s
% Computational Cost: add. (102->61), mult. (103->65), div. (0->0), fcn. (79->8), ass. (0->24)
t58 = pkin(2) + pkin(5);
t43 = sin(pkin(8));
t57 = pkin(3) * t43;
t45 = -pkin(6) - qJ(3);
t56 = rSges(5,3) - t45;
t55 = rSges(6,3) + pkin(7) - t45;
t46 = sin(qJ(1));
t47 = cos(qJ(1));
t54 = t47 * pkin(1) + t46 * qJ(2);
t53 = rSges(4,3) + qJ(3);
t42 = pkin(8) + qJ(4);
t44 = cos(pkin(8));
t52 = t44 * pkin(3) + t58;
t51 = rSges(4,1) * t43 + rSges(4,2) * t44;
t36 = qJ(5) + t42;
t32 = sin(t36);
t33 = cos(t36);
t34 = sin(t42);
t50 = rSges(6,1) * t32 + rSges(6,2) * t33 + pkin(4) * t34 + t57;
t39 = t46 * pkin(1);
t49 = g(1) * t54 + g(2) * t39;
t35 = cos(t42);
t48 = rSges(5,1) * t34 + rSges(5,2) * t35 + t57;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (rSges(2,1) * t47 - t46 * rSges(2,2)) + g(2) * (t46 * rSges(2,1) + rSges(2,2) * t47) + g(3) * (pkin(5) + rSges(2,3))) - m(3) * (g(1) * (-rSges(3,2) * t47 + t46 * rSges(3,3) + t54) + g(2) * (-t46 * rSges(3,2) + t39 + (-rSges(3,3) - qJ(2)) * t47) + g(3) * (pkin(5) + rSges(3,1))) - m(4) * (g(3) * (rSges(4,1) * t44 - rSges(4,2) * t43 + t58) + (g(1) * t51 + g(2) * t53) * t46 + (g(1) * t53 + g(2) * (-qJ(2) - t51)) * t47 + t49) - m(5) * (g(3) * (rSges(5,1) * t35 - rSges(5,2) * t34 + t52) + (g(1) * t48 + g(2) * t56) * t46 + (g(1) * t56 + g(2) * (-qJ(2) - t48)) * t47 + t49) - m(6) * (g(3) * (rSges(6,1) * t33 - rSges(6,2) * t32 + pkin(4) * t35 + t52) + (g(1) * t50 + g(2) * t55) * t46 + (g(1) * t55 + g(2) * (-qJ(2) - t50)) * t47 + t49);
U = t1;

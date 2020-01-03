% Calculate potential energy for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPRPR3_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_energypot_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_energypot_fixb_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:52
% EndTime: 2020-01-03 11:35:53
% DurationCPUTime: 0.38s
% Computational Cost: add. (152->66), mult. (105->76), div. (0->0), fcn. (85->10), ass. (0->26)
t58 = rSges(6,3) + pkin(7);
t41 = cos(pkin(9));
t57 = t41 * pkin(4);
t45 = cos(qJ(1));
t56 = t45 * pkin(1);
t42 = sin(qJ(5));
t54 = t41 * t42;
t44 = cos(qJ(5));
t53 = t41 * t44;
t52 = pkin(5) + qJ(2);
t39 = qJ(1) + pkin(8);
t35 = sin(t39);
t43 = sin(qJ(1));
t38 = t43 * pkin(1);
t51 = pkin(2) * t35 + t38;
t50 = -rSges(5,3) - qJ(4);
t49 = pkin(6) + t52;
t37 = qJ(3) + t39;
t33 = sin(t37);
t48 = t33 * pkin(3) + t51;
t36 = cos(t39);
t47 = -pkin(2) * t36 - t56;
t40 = sin(pkin(9));
t46 = rSges(5,1) * t41 - rSges(5,2) * t40;
t34 = cos(t37);
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (pkin(5) + rSges(2,3)) + g(2) * (t43 * rSges(2,1) + t45 * rSges(2,2)) + g(3) * (-t45 * rSges(2,1) + t43 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) + t52) + g(2) * (t35 * rSges(3,1) + t36 * rSges(3,2) + t38) + g(3) * (-t36 * rSges(3,1) + t35 * rSges(3,2) - t56)) - m(4) * (g(1) * (rSges(4,3) + t49) + g(2) * (t33 * rSges(4,1) + t34 * rSges(4,2) + t51) + g(3) * (-t34 * rSges(4,1) + t33 * rSges(4,2) + t47)) - m(5) * (g(1) * (t40 * rSges(5,1) + t41 * rSges(5,2) + t49) + g(2) * t48 + g(3) * t47 + (g(2) * t46 + g(3) * t50) * t33 + (g(2) * t50 + g(3) * (-pkin(3) - t46)) * t34) - m(6) * (g(1) * (-t58 * t41 + t49) + g(2) * (t33 * t57 - t34 * qJ(4) + (t33 * t53 - t34 * t42) * rSges(6,1) + (-t33 * t54 - t34 * t44) * rSges(6,2) + t48) + g(3) * (t47 + (-t42 * rSges(6,1) - t44 * rSges(6,2) - qJ(4)) * t33 + (-t53 * rSges(6,1) + t54 * rSges(6,2) - pkin(3) - t57) * t34) + (g(1) * (rSges(6,1) * t44 - rSges(6,2) * t42 + pkin(4)) + (g(2) * t33 - g(3) * t34) * t58) * t40);
U = t1;

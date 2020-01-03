% Calculate potential energy for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RPRR7_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR7_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:39
% EndTime: 2019-12-31 16:53:40
% DurationCPUTime: 0.24s
% Computational Cost: add. (89->56), mult. (101->72), div. (0->0), fcn. (85->8), ass. (0->24)
t56 = rSges(5,3) + pkin(6);
t38 = sin(pkin(7));
t54 = t38 * pkin(2) + pkin(4);
t37 = pkin(7) + qJ(3);
t34 = sin(t37);
t53 = rSges(4,2) * t34;
t41 = sin(qJ(4));
t42 = sin(qJ(1));
t52 = t42 * t41;
t43 = cos(qJ(4));
t51 = t42 * t43;
t35 = cos(t37);
t44 = cos(qJ(1));
t50 = t44 * t35;
t49 = t44 * t41;
t48 = t44 * t43;
t39 = cos(pkin(7));
t32 = t39 * pkin(2) + pkin(1);
t40 = -pkin(5) - qJ(2);
t47 = t42 * t32 + t44 * t40;
t46 = rSges(3,3) + qJ(2);
t45 = rSges(3,1) * t39 - rSges(3,2) * t38 + pkin(1);
t31 = t44 * t32;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t44 * rSges(2,1) - t42 * rSges(2,2)) + g(2) * (t42 * rSges(2,1) + t44 * rSges(2,2)) + g(3) * (pkin(4) + rSges(2,3))) - m(3) * (g(3) * (t38 * rSges(3,1) + t39 * rSges(3,2) + pkin(4)) + (g(1) * t45 - g(2) * t46) * t44 + (g(1) * t46 + g(2) * t45) * t42) - m(4) * (g(1) * (rSges(4,1) * t50 - t44 * t53 + t31) + g(2) * (-t44 * rSges(4,3) + t47) + g(3) * (t34 * rSges(4,1) + t35 * rSges(4,2) + t54) + (g(1) * (rSges(4,3) - t40) + g(2) * (rSges(4,1) * t35 - t53)) * t42) - m(5) * (g(1) * (pkin(3) * t50 + t31 - t42 * t40 + (t35 * t48 + t52) * rSges(5,1) + (-t35 * t49 + t51) * rSges(5,2)) + g(2) * (t42 * t35 * pkin(3) + (t35 * t51 - t49) * rSges(5,1) + (-t35 * t52 - t48) * rSges(5,2) + t47) + g(3) * (-t56 * t35 + t54) + (g(3) * (rSges(5,1) * t43 - rSges(5,2) * t41 + pkin(3)) + (g(1) * t44 + g(2) * t42) * t56) * t34);
U = t1;

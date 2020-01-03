% Calculate potential energy for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4PPRR4_energypot_fixb_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_energypot_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_energypot_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_energypot_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_energypot_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR4_energypot_fixb_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:23
% EndTime: 2019-12-31 16:18:23
% DurationCPUTime: 0.26s
% Computational Cost: add. (89->56), mult. (101->72), div. (0->0), fcn. (85->8), ass. (0->24)
t56 = rSges(5,3) + pkin(5);
t37 = pkin(7) + qJ(3);
t34 = sin(t37);
t54 = rSges(4,2) * t34;
t39 = sin(pkin(6));
t43 = sin(qJ(4));
t53 = t39 * t43;
t44 = cos(qJ(4));
t52 = t39 * t44;
t35 = cos(t37);
t41 = cos(pkin(6));
t51 = t41 * t35;
t50 = t41 * t43;
t49 = t41 * t44;
t40 = cos(pkin(7));
t33 = t40 * pkin(2) + pkin(1);
t42 = -pkin(4) - qJ(2);
t48 = t39 * t33 + t41 * t42;
t47 = rSges(3,3) + qJ(2);
t38 = sin(pkin(7));
t46 = t38 * pkin(2) + qJ(1);
t45 = rSges(3,1) * t40 - rSges(3,2) * t38 + pkin(1);
t31 = t41 * t33;
t1 = -m(1) * (g(1) * rSges(1,1) + g(2) * rSges(1,2) + g(3) * rSges(1,3)) - m(2) * (g(1) * (t41 * rSges(2,1) - t39 * rSges(2,2)) + g(2) * (t39 * rSges(2,1) + t41 * rSges(2,2)) + g(3) * (qJ(1) + rSges(2,3))) - m(3) * (g(3) * (t38 * rSges(3,1) + t40 * rSges(3,2) + qJ(1)) + (g(1) * t45 - g(2) * t47) * t41 + (g(1) * t47 + g(2) * t45) * t39) - m(4) * (g(1) * (rSges(4,1) * t51 - t41 * t54 + t31) + g(2) * (-t41 * rSges(4,3) + t48) + g(3) * (t34 * rSges(4,1) + t35 * rSges(4,2) + t46) + (g(1) * (rSges(4,3) - t42) + g(2) * (rSges(4,1) * t35 - t54)) * t39) - m(5) * (g(1) * (pkin(3) * t51 + t31 - t39 * t42 + (t35 * t49 + t53) * rSges(5,1) + (-t35 * t50 + t52) * rSges(5,2)) + g(2) * (t39 * t35 * pkin(3) + (t35 * t52 - t50) * rSges(5,1) + (-t35 * t53 - t49) * rSges(5,2) + t48) + g(3) * (-t56 * t35 + t46) + (g(3) * (rSges(5,1) * t44 - rSges(5,2) * t43 + pkin(3)) + (g(1) * t41 + g(2) * t39) * t56) * t34);
U = t1;

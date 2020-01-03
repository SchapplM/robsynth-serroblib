% Calculate Gravitation load on the joints for
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
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:56
% EndTime: 2020-01-03 11:35:58
% DurationCPUTime: 0.22s
% Computational Cost: add. (262->64), mult. (180->87), div. (0->0), fcn. (158->10), ass. (0->37)
t57 = rSges(6,3) + pkin(7);
t56 = -m(5) - m(6);
t34 = qJ(1) + pkin(8);
t31 = qJ(3) + t34;
t27 = sin(t31);
t35 = sin(pkin(9));
t55 = t27 * t35;
t36 = cos(pkin(9));
t54 = t27 * t36;
t28 = cos(t31);
t53 = t28 * t35;
t52 = t28 * t36;
t37 = sin(qJ(5));
t51 = t36 * t37;
t39 = cos(qJ(5));
t50 = t36 * t39;
t49 = t28 * pkin(3) + t27 * qJ(4);
t48 = t27 * rSges(4,1) + t28 * rSges(4,2);
t29 = sin(t34);
t38 = sin(qJ(1));
t32 = t38 * pkin(1);
t47 = pkin(2) * t29 + t32;
t30 = cos(t34);
t40 = cos(qJ(1));
t33 = t40 * pkin(1);
t46 = pkin(2) * t30 + t33;
t45 = t28 * rSges(4,1) - rSges(4,2) * t27;
t7 = -t27 * t39 + t28 * t51;
t8 = t27 * t37 + t28 * t50;
t44 = t8 * rSges(6,1) - t7 * rSges(6,2) + pkin(4) * t52 + t57 * t53 + t49;
t43 = rSges(5,1) * t52 - rSges(5,2) * t53 + t27 * rSges(5,3) + t49;
t23 = t27 * pkin(3);
t5 = -t27 * t51 - t28 * t39;
t6 = t27 * t50 - t28 * t37;
t42 = t6 * rSges(6,1) + t5 * rSges(6,2) + pkin(4) * t54 - qJ(4) * t28 + t57 * t55 + t23;
t41 = -rSges(5,2) * t55 + rSges(5,1) * t54 + t23 + (-rSges(5,3) - qJ(4)) * t28;
t1 = [-m(2) * (g(2) * (rSges(2,1) * t40 - t38 * rSges(2,2)) + g(3) * (t38 * rSges(2,1) + rSges(2,2) * t40)) - m(3) * (g(2) * (rSges(3,1) * t30 - rSges(3,2) * t29 + t33) + g(3) * (rSges(3,1) * t29 + rSges(3,2) * t30 + t32)) - m(4) * (g(2) * (t45 + t46) + g(3) * (t47 + t48)) - m(5) * (g(2) * (t43 + t46) + g(3) * (t41 + t47)) - m(6) * (g(2) * (t44 + t46) + g(3) * (t42 + t47)), (-m(3) - m(4) + t56) * g(1), -m(4) * (g(2) * t45 + g(3) * t48) - m(5) * (g(2) * t43 + g(3) * t41) - m(6) * (g(2) * t44 + g(3) * t42), t56 * (-g(2) * t28 - g(3) * t27), -m(6) * (g(2) * (rSges(6,1) * t5 - rSges(6,2) * t6) + g(3) * (rSges(6,1) * t7 + rSges(6,2) * t8) + g(1) * (-rSges(6,1) * t37 - rSges(6,2) * t39) * t35)];
taug = t1(:);

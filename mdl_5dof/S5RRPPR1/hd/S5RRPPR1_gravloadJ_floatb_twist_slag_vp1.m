% Calculate Gravitation load on the joints for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:43
% EndTime: 2020-01-03 11:55:44
% DurationCPUTime: 0.23s
% Computational Cost: add. (243->57), mult. (150->67), div. (0->0), fcn. (114->10), ass. (0->31)
t57 = rSges(6,3) + pkin(7) + qJ(4);
t32 = pkin(9) + qJ(5);
t25 = sin(t32);
t26 = cos(t32);
t56 = rSges(6,1) * t26 - rSges(6,2) * t25;
t55 = rSges(5,3) + qJ(4);
t35 = cos(pkin(9));
t54 = -rSges(5,2) * sin(pkin(9)) + pkin(3) + rSges(5,1) * t35;
t53 = pkin(4) * t35 + pkin(3) + t56;
t52 = -m(5) - m(6);
t33 = qJ(1) + qJ(2);
t28 = sin(t33);
t29 = cos(t33);
t47 = t28 * rSges(3,1) + t29 * rSges(3,2);
t27 = pkin(8) + t33;
t20 = sin(t27);
t21 = cos(t27);
t23 = pkin(2) * t28;
t46 = t20 * rSges(4,1) + t21 * rSges(4,2) + t23;
t45 = t29 * rSges(3,1) - rSges(3,2) * t28;
t24 = pkin(2) * t29;
t44 = t21 * rSges(4,1) - rSges(4,2) * t20 + t24;
t42 = t55 * t20 + t54 * t21 + t24;
t41 = t53 * t20 - t57 * t21 + t23;
t40 = t57 * t20 + t53 * t21 + t24;
t39 = t54 * t20 - t55 * t21 + t23;
t38 = cos(qJ(1));
t37 = sin(qJ(1));
t31 = t38 * pkin(1);
t30 = t37 * pkin(1);
t1 = [-m(2) * (g(2) * (rSges(2,1) * t38 - t37 * rSges(2,2)) + g(3) * (t37 * rSges(2,1) + rSges(2,2) * t38)) - m(3) * (g(2) * (t31 + t45) + g(3) * (t30 + t47)) - m(4) * (g(2) * (t31 + t44) + g(3) * (t30 + t46)) - m(5) * (g(2) * (t31 + t42) + g(3) * (t30 + t39)) - m(6) * (g(2) * (t31 + t40) + g(3) * (t30 + t41)), -m(3) * (g(2) * t45 + g(3) * t47) - m(4) * (g(2) * t44 + g(3) * t46) - m(5) * (g(2) * t42 + g(3) * t39) - m(6) * (g(2) * t40 + g(3) * t41), (-m(4) + t52) * g(1), t52 * (-g(2) * t21 - g(3) * t20), -m(6) * (g(1) * t56 + (-g(2) * t20 + g(3) * t21) * (rSges(6,1) * t25 + rSges(6,2) * t26))];
taug = t1(:);

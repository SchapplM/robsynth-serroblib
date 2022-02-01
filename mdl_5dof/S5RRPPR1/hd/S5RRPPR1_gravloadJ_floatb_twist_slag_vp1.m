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
% m [6x1]
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:51:19
% EndTime: 2022-01-20 09:51:20
% DurationCPUTime: 0.22s
% Computational Cost: add. (243->57), mult. (150->67), div. (0->0), fcn. (114->10), ass. (0->31)
t53 = pkin(7) + qJ(4) + rSges(6,3);
t26 = pkin(9) + qJ(5);
t20 = sin(t26);
t21 = cos(t26);
t52 = rSges(6,1) * t21 - rSges(6,2) * t20;
t51 = rSges(5,3) + qJ(4);
t29 = cos(pkin(9));
t50 = rSges(5,2) * sin(pkin(9)) - pkin(3) - rSges(5,1) * t29;
t49 = -pkin(4) * t29 - pkin(3) - t52;
t48 = -m(5) - m(6);
t31 = sin(qJ(1));
t47 = pkin(1) * t31;
t27 = qJ(1) + qJ(2);
t23 = sin(t27);
t46 = pkin(2) * t23;
t24 = cos(t27);
t41 = t24 * rSges(3,1) - rSges(3,2) * t23;
t22 = pkin(8) + t27;
t16 = sin(t22);
t17 = cos(t22);
t19 = pkin(2) * t24;
t40 = t17 * rSges(4,1) - rSges(4,2) * t16 + t19;
t39 = -rSges(3,1) * t23 - rSges(3,2) * t24;
t37 = -rSges(4,1) * t16 - rSges(4,2) * t17 - t46;
t36 = t51 * t16 - t50 * t17 + t19;
t35 = t53 * t16 - t49 * t17 + t19;
t34 = t50 * t16 + t51 * t17 - t46;
t33 = t49 * t16 + t53 * t17 - t46;
t32 = cos(qJ(1));
t25 = t32 * pkin(1);
t1 = [-m(2) * (g(1) * (-t31 * rSges(2,1) - rSges(2,2) * t32) + g(2) * (rSges(2,1) * t32 - t31 * rSges(2,2))) - m(3) * (g(1) * (t39 - t47) + g(2) * (t25 + t41)) - m(4) * (g(1) * (t37 - t47) + g(2) * (t25 + t40)) - m(5) * (g(1) * (t34 - t47) + g(2) * (t25 + t36)) - m(6) * (g(1) * (t33 - t47) + g(2) * (t25 + t35)), -m(3) * (g(1) * t39 + g(2) * t41) - m(4) * (g(1) * t37 + g(2) * t40) - m(5) * (g(1) * t34 + g(2) * t36) - m(6) * (g(1) * t33 + g(2) * t35), (-m(4) + t48) * g(3), t48 * (g(1) * t16 - g(2) * t17), -m(6) * (g(3) * t52 + (g(1) * t17 + g(2) * t16) * (-rSges(6,1) * t20 - rSges(6,2) * t21))];
taug = t1(:);

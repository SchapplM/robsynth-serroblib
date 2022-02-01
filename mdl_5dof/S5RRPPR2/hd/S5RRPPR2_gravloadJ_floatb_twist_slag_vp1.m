% Calculate Gravitation load on the joints for
% S5RRPPR2
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:26
% EndTime: 2022-01-20 10:05:27
% DurationCPUTime: 0.22s
% Computational Cost: add. (281->68), mult. (192->89), div. (0->0), fcn. (168->10), ass. (0->39)
t55 = rSges(6,3) + pkin(7);
t54 = -m(5) - m(6);
t29 = qJ(1) + qJ(2);
t26 = sin(t29);
t53 = pkin(2) * t26;
t25 = pkin(8) + t29;
t22 = sin(t25);
t52 = g(1) * t22;
t33 = sin(qJ(1));
t51 = t33 * pkin(1);
t23 = cos(t25);
t30 = sin(pkin(9));
t50 = t23 * t30;
t31 = cos(pkin(9));
t49 = t23 * t31;
t32 = sin(qJ(5));
t48 = t31 * t32;
t34 = cos(qJ(5));
t47 = t31 * t34;
t27 = cos(t29);
t24 = pkin(2) * t27;
t46 = t23 * pkin(3) + t22 * qJ(4) + t24;
t45 = t23 * qJ(4) - t53;
t44 = t27 * rSges(3,1) - t26 * rSges(3,2);
t43 = t23 * rSges(4,1) - t22 * rSges(4,2) + t24;
t5 = t22 * t48 + t23 * t34;
t6 = -t22 * t47 + t23 * t32;
t42 = t6 * rSges(6,1) + t5 * rSges(6,2) + t45;
t41 = -t26 * rSges(3,1) - t27 * rSges(3,2);
t7 = t22 * t34 - t23 * t48;
t8 = t22 * t32 + t23 * t47;
t40 = t8 * rSges(6,1) + t7 * rSges(6,2) + pkin(4) * t49 + t55 * t50 + t46;
t39 = -t22 * rSges(4,1) - t23 * rSges(4,2) - t53;
t38 = rSges(5,1) * t49 - rSges(5,2) * t50 + t22 * rSges(5,3) + t46;
t37 = (-pkin(4) * t31 - t55 * t30 - pkin(3)) * t52;
t36 = t23 * rSges(5,3) + t45 + (-rSges(5,1) * t31 + rSges(5,2) * t30 - pkin(3)) * t22;
t35 = cos(qJ(1));
t28 = t35 * pkin(1);
t1 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - t35 * rSges(2,2)) + g(2) * (t35 * rSges(2,1) - t33 * rSges(2,2))) - m(3) * (g(1) * (t41 - t51) + g(2) * (t28 + t44)) - m(4) * (g(1) * (t39 - t51) + g(2) * (t28 + t43)) - m(5) * (g(1) * (t36 - t51) + g(2) * (t28 + t38)) - m(6) * (g(1) * (t42 - t51) + g(2) * (t28 + t40) + t37), -m(3) * (g(1) * t41 + g(2) * t44) - m(4) * (g(1) * t39 + g(2) * t43) - m(5) * (g(1) * t36 + g(2) * t38) - m(6) * (g(1) * t42 + g(2) * t40 + t37), (-m(4) + t54) * g(3), t54 * (-g(2) * t23 + t52), -m(6) * (g(1) * (t7 * rSges(6,1) - t8 * rSges(6,2)) + g(2) * (-t5 * rSges(6,1) + t6 * rSges(6,2)) + g(3) * (-rSges(6,1) * t32 - rSges(6,2) * t34) * t30)];
taug = t1(:);

% Calculate Gravitation load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:33
% EndTime: 2022-01-20 11:30:34
% DurationCPUTime: 0.23s
% Computational Cost: add. (308->60), mult. (162->72), div. (0->0), fcn. (120->10), ass. (0->37)
t49 = rSges(6,3) + pkin(8);
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t48 = rSges(6,1) * t25 - rSges(6,2) * t23;
t47 = -pkin(4) - t48;
t24 = sin(qJ(1));
t46 = pkin(1) * t24;
t22 = qJ(1) + qJ(2);
t18 = sin(t22);
t45 = pkin(2) * t18;
t20 = qJ(3) + t22;
t16 = sin(t20);
t44 = pkin(3) * t16;
t17 = cos(t20);
t41 = t17 * rSges(4,1) - rSges(4,2) * t16;
t19 = cos(t22);
t40 = t19 * rSges(3,1) - rSges(3,2) * t18;
t14 = pkin(2) * t19;
t39 = t14 + t41;
t15 = pkin(9) + t20;
t10 = sin(t15);
t11 = cos(t15);
t12 = pkin(3) * t17;
t38 = t11 * rSges(5,1) - rSges(5,2) * t10 + t12;
t37 = -rSges(3,1) * t18 - rSges(3,2) * t19;
t36 = -rSges(4,1) * t16 - rSges(4,2) * t17;
t34 = t14 + t38;
t33 = t36 - t45;
t32 = -rSges(5,1) * t10 - rSges(5,2) * t11 - t44;
t31 = t49 * t10 - t47 * t11 + t12;
t30 = t14 + t31;
t29 = t32 - t45;
t28 = t47 * t10 + t49 * t11 - t44;
t27 = t28 - t45;
t26 = cos(qJ(1));
t21 = t26 * pkin(1);
t1 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - rSges(2,2) * t26) + g(2) * (rSges(2,1) * t26 - t24 * rSges(2,2))) - m(3) * (g(1) * (t37 - t46) + g(2) * (t21 + t40)) - m(4) * (g(1) * (t33 - t46) + g(2) * (t21 + t39)) - m(5) * (g(1) * (t29 - t46) + g(2) * (t21 + t34)) - m(6) * (g(1) * (t27 - t46) + g(2) * (t21 + t30)), -m(3) * (g(1) * t37 + g(2) * t40) - m(4) * (g(1) * t33 + g(2) * t39) - m(5) * (g(1) * t29 + g(2) * t34) - m(6) * (g(1) * t27 + g(2) * t30), -m(4) * (g(1) * t36 + g(2) * t41) - m(5) * (g(1) * t32 + g(2) * t38) - m(6) * (g(1) * t28 + g(2) * t31), (-m(5) - m(6)) * g(3), -m(6) * (g(3) * t48 + (g(1) * t11 + g(2) * t10) * (-rSges(6,1) * t23 - rSges(6,2) * t25))];
taug = t1(:);

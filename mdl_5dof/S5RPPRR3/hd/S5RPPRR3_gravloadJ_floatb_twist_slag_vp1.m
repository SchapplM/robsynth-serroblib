% Calculate Gravitation load on the joints for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:08
% EndTime: 2022-01-23 09:14:09
% DurationCPUTime: 0.25s
% Computational Cost: add. (186->52), mult. (139->67), div. (0->0), fcn. (105->10), ass. (0->30)
t16 = pkin(9) + qJ(4);
t10 = cos(t16);
t12 = qJ(5) + t16;
t5 = sin(t12);
t6 = cos(t12);
t30 = t6 * rSges(6,1) - rSges(6,2) * t5;
t40 = pkin(4) * t10 + t30;
t17 = qJ(1) + pkin(8);
t11 = cos(t17);
t9 = sin(t17);
t39 = g(1) * t11 + g(2) * t9;
t21 = sin(qJ(1));
t36 = pkin(1) * t21;
t19 = cos(pkin(9));
t7 = t19 * pkin(3) + pkin(2);
t20 = -pkin(6) - qJ(3);
t34 = rSges(5,3) - t20;
t33 = rSges(6,3) + pkin(7) - t20;
t32 = rSges(4,3) + qJ(3);
t31 = -m(4) - m(5) - m(6);
t29 = -rSges(6,1) * t5 - rSges(6,2) * t6;
t8 = sin(t16);
t28 = rSges(5,1) * t10 - rSges(5,2) * t8;
t22 = cos(qJ(1));
t14 = t22 * pkin(1);
t27 = -g(1) * t36 + g(2) * t14;
t26 = t7 + t40;
t25 = t28 + t7;
t24 = rSges(4,1) * t19 - rSges(4,2) * sin(pkin(9)) + pkin(2);
t1 = [-m(2) * (g(1) * (-t21 * rSges(2,1) - rSges(2,2) * t22) + g(2) * (rSges(2,1) * t22 - t21 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t9 - rSges(3,2) * t11 - t36) + g(2) * (rSges(3,1) * t11 - rSges(3,2) * t9 + t14)) - m(4) * ((-g(1) * t24 + g(2) * t32) * t9 + (g(1) * t32 + g(2) * t24) * t11 + t27) - m(5) * ((-g(1) * t25 + g(2) * t34) * t9 + (g(1) * t34 + g(2) * t25) * t11 + t27) - m(6) * ((-g(1) * t26 + g(2) * t33) * t9 + (g(1) * t33 + g(2) * t26) * t11 + t27), (-m(3) + t31) * g(3), t31 * (g(1) * t9 - g(2) * t11), (-m(5) * t28 - m(6) * t40) * g(3) + t39 * (-m(5) * (-rSges(5,1) * t8 - rSges(5,2) * t10) - m(6) * (-pkin(4) * t8 + t29)), -m(6) * (g(3) * t30 + t39 * t29)];
taug = t1(:);

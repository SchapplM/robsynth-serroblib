% Calculate Gravitation load on the joints for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:33:58
% EndTime: 2022-01-20 10:33:59
% DurationCPUTime: 0.22s
% Computational Cost: add. (285->57), mult. (152->69), div. (0->0), fcn. (112->10), ass. (0->35)
t48 = rSges(6,3) + pkin(8);
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t47 = rSges(6,1) * t25 - rSges(6,2) * t23;
t46 = -pkin(4) - t47;
t24 = sin(qJ(1));
t45 = pkin(1) * t24;
t22 = qJ(1) + qJ(2);
t19 = sin(t22);
t44 = pkin(2) * t19;
t18 = pkin(9) + t22;
t15 = cos(t18);
t20 = cos(t22);
t16 = pkin(2) * t20;
t41 = pkin(3) * t15 + t16;
t17 = qJ(4) + t18;
t11 = sin(t17);
t12 = cos(t17);
t40 = t12 * rSges(5,1) - rSges(5,2) * t11;
t39 = t20 * rSges(3,1) - rSges(3,2) * t19;
t14 = sin(t18);
t38 = t15 * rSges(4,1) - rSges(4,2) * t14 + t16;
t37 = -pkin(3) * t14 - t44;
t36 = -rSges(3,1) * t19 - rSges(3,2) * t20;
t35 = -rSges(5,1) * t11 - rSges(5,2) * t12;
t33 = t40 + t41;
t32 = t48 * t11 - t46 * t12;
t31 = -rSges(4,1) * t14 - rSges(4,2) * t15 - t44;
t30 = t46 * t11 + t48 * t12;
t29 = t32 + t41;
t28 = t35 + t37;
t27 = t30 + t37;
t26 = cos(qJ(1));
t21 = t26 * pkin(1);
t1 = [-m(2) * (g(1) * (-t24 * rSges(2,1) - rSges(2,2) * t26) + g(2) * (rSges(2,1) * t26 - t24 * rSges(2,2))) - m(3) * (g(1) * (t36 - t45) + g(2) * (t21 + t39)) - m(4) * (g(1) * (t31 - t45) + g(2) * (t21 + t38)) - m(5) * (g(1) * (t28 - t45) + g(2) * (t21 + t33)) - m(6) * (g(1) * (t27 - t45) + g(2) * (t21 + t29)), -m(3) * (g(1) * t36 + g(2) * t39) - m(4) * (g(1) * t31 + g(2) * t38) - m(5) * (g(1) * t28 + g(2) * t33) - m(6) * (g(1) * t27 + g(2) * t29), (-m(4) - m(5) - m(6)) * g(3), -m(5) * (g(1) * t35 + g(2) * t40) - m(6) * (g(1) * t30 + g(2) * t32), -m(6) * (g(3) * t47 + (g(1) * t12 + g(2) * t11) * (-rSges(6,1) * t23 - rSges(6,2) * t25))];
taug = t1(:);

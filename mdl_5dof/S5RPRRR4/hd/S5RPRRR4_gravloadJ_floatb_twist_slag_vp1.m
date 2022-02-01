% Calculate Gravitation load on the joints for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:28
% EndTime: 2022-01-23 09:34:29
% DurationCPUTime: 0.18s
% Computational Cost: add. (266->54), mult. (140->66), div. (0->0), fcn. (102->10), ass. (0->33)
t45 = rSges(6,3) + pkin(8);
t22 = sin(qJ(5));
t24 = cos(qJ(5));
t44 = rSges(6,1) * t24 - rSges(6,2) * t22;
t43 = -pkin(4) - t44;
t23 = sin(qJ(1));
t42 = pkin(1) * t23;
t21 = qJ(1) + pkin(9);
t19 = qJ(3) + t21;
t14 = sin(t19);
t41 = pkin(3) * t14;
t18 = cos(t21);
t25 = cos(qJ(1));
t20 = t25 * pkin(1);
t38 = pkin(2) * t18 + t20;
t15 = cos(t19);
t37 = t15 * rSges(4,1) - rSges(4,2) * t14;
t16 = qJ(4) + t19;
t11 = sin(t16);
t12 = cos(t16);
t36 = t12 * rSges(5,1) - rSges(5,2) * t11;
t10 = pkin(3) * t15;
t35 = t10 + t36;
t17 = sin(t21);
t34 = -pkin(2) * t17 - t42;
t33 = -rSges(4,1) * t14 - rSges(4,2) * t15;
t32 = -rSges(5,1) * t11 - rSges(5,2) * t12;
t30 = t45 * t11 - t43 * t12;
t29 = t32 - t41;
t28 = t10 + t30;
t27 = t43 * t11 + t45 * t12;
t26 = t27 - t41;
t1 = [-m(2) * (g(1) * (-t23 * rSges(2,1) - rSges(2,2) * t25) + g(2) * (rSges(2,1) * t25 - t23 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t17 - rSges(3,2) * t18 - t42) + g(2) * (rSges(3,1) * t18 - rSges(3,2) * t17 + t20)) - m(4) * (g(1) * (t33 + t34) + g(2) * (t37 + t38)) - m(5) * (g(1) * (t29 + t34) + g(2) * (t35 + t38)) - m(6) * (g(1) * (t26 + t34) + g(2) * (t28 + t38)), (-m(3) - m(4) - m(5) - m(6)) * g(3), -m(4) * (g(1) * t33 + g(2) * t37) - m(5) * (g(1) * t29 + g(2) * t35) - m(6) * (g(1) * t26 + g(2) * t28), -m(5) * (g(1) * t32 + g(2) * t36) - m(6) * (g(1) * t27 + g(2) * t30), -m(6) * (g(3) * t44 + (g(1) * t12 + g(2) * t11) * (-rSges(6,1) * t22 - rSges(6,2) * t24))];
taug = t1(:);

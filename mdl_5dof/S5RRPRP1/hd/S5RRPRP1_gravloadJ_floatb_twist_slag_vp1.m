% Calculate Gravitation load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:33
% EndTime: 2022-01-20 10:19:35
% DurationCPUTime: 0.28s
% Computational Cost: add. (244->59), mult. (166->71), div. (0->0), fcn. (127->8), ass. (0->30)
t52 = rSges(6,1) + pkin(4);
t26 = sin(qJ(4));
t28 = cos(qJ(4));
t51 = -t26 * rSges(6,2) + t52 * t28;
t50 = pkin(7) + rSges(5,3);
t49 = qJ(5) + pkin(7) + rSges(6,3);
t47 = t28 * rSges(5,1) - t26 * rSges(5,2);
t46 = -pkin(3) - t47;
t45 = -pkin(3) - t51;
t24 = qJ(1) + qJ(2);
t20 = sin(t24);
t44 = pkin(2) * t20;
t27 = sin(qJ(1));
t43 = t27 * pkin(1);
t21 = cos(t24);
t38 = t21 * rSges(3,1) - rSges(3,2) * t20;
t19 = pkin(8) + t24;
t15 = sin(t19);
t16 = cos(t19);
t17 = pkin(2) * t21;
t37 = t16 * rSges(4,1) - rSges(4,2) * t15 + t17;
t36 = -rSges(3,1) * t20 - rSges(3,2) * t21;
t35 = -rSges(4,1) * t15 - rSges(4,2) * t16 - t44;
t34 = t50 * t15 - t46 * t16 + t17;
t33 = t49 * t15 - t45 * t16 + t17;
t32 = t46 * t15 + t50 * t16 - t44;
t31 = t45 * t15 + t49 * t16 - t44;
t29 = cos(qJ(1));
t23 = t29 * pkin(1);
t1 = [-m(2) * (g(1) * (-t27 * rSges(2,1) - t29 * rSges(2,2)) + g(2) * (t29 * rSges(2,1) - t27 * rSges(2,2))) - m(3) * (g(1) * (t36 - t43) + g(2) * (t23 + t38)) - m(4) * (g(1) * (t35 - t43) + g(2) * (t23 + t37)) - m(5) * (g(1) * (t32 - t43) + g(2) * (t23 + t34)) - m(6) * (g(1) * (t31 - t43) + g(2) * (t23 + t33)), -m(3) * (g(1) * t36 + g(2) * t38) - m(4) * (g(1) * t35 + g(2) * t37) - m(5) * (g(1) * t32 + g(2) * t34) - m(6) * (g(1) * t31 + g(2) * t33), (-m(4) - m(5) - m(6)) * g(3), (-m(5) * t47 - m(6) * t51) * g(3) + (g(1) * t16 + g(2) * t15) * (-m(5) * (-rSges(5,1) * t26 - rSges(5,2) * t28) - m(6) * (-rSges(6,2) * t28 - t52 * t26)), -m(6) * (g(1) * t15 - g(2) * t16)];
taug = t1(:);

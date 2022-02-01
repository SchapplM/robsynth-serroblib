% Calculate Gravitation load on the joints for
% S5RRPRR4
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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:03
% EndTime: 2022-01-20 10:48:05
% DurationCPUTime: 0.27s
% Computational Cost: add. (272->61), mult. (177->73), div. (0->0), fcn. (135->10), ass. (0->35)
t31 = cos(qJ(4));
t27 = qJ(4) + qJ(5);
t21 = sin(t27);
t23 = cos(t27);
t43 = t23 * rSges(6,1) - rSges(6,2) * t21;
t58 = t31 * pkin(4) + t43;
t57 = pkin(7) + rSges(5,3);
t56 = pkin(8) + pkin(7) + rSges(6,3);
t29 = sin(qJ(4));
t55 = t31 * rSges(5,1) - t29 * rSges(5,2);
t54 = -pkin(3) - t55;
t53 = -pkin(3) - t58;
t28 = qJ(1) + qJ(2);
t20 = pkin(9) + t28;
t16 = sin(t20);
t17 = cos(t20);
t52 = g(1) * t17 + g(2) * t16;
t22 = sin(t28);
t51 = pkin(2) * t22;
t30 = sin(qJ(1));
t48 = t30 * pkin(1);
t24 = cos(t28);
t44 = t24 * rSges(3,1) - t22 * rSges(3,2);
t18 = pkin(2) * t24;
t42 = t17 * rSges(4,1) - t16 * rSges(4,2) + t18;
t41 = -t22 * rSges(3,1) - t24 * rSges(3,2);
t40 = -rSges(6,1) * t21 - rSges(6,2) * t23;
t39 = -t16 * rSges(4,1) - t17 * rSges(4,2) - t51;
t38 = t57 * t16 - t54 * t17 + t18;
t37 = t56 * t16 - t53 * t17 + t18;
t36 = t54 * t16 + t57 * t17 - t51;
t35 = t53 * t16 + t56 * t17 - t51;
t32 = cos(qJ(1));
t26 = t32 * pkin(1);
t1 = [-m(2) * (g(1) * (-t30 * rSges(2,1) - t32 * rSges(2,2)) + g(2) * (t32 * rSges(2,1) - t30 * rSges(2,2))) - m(3) * (g(1) * (t41 - t48) + g(2) * (t26 + t44)) - m(4) * (g(1) * (t39 - t48) + g(2) * (t26 + t42)) - m(5) * (g(1) * (t36 - t48) + g(2) * (t26 + t38)) - m(6) * (g(1) * (t35 - t48) + g(2) * (t26 + t37)), -m(3) * (g(1) * t41 + g(2) * t44) - m(4) * (g(1) * t39 + g(2) * t42) - m(5) * (g(1) * t36 + g(2) * t38) - m(6) * (g(1) * t35 + g(2) * t37), (-m(4) - m(5) - m(6)) * g(3), (-m(5) * t55 - m(6) * t58) * g(3) + t52 * (-m(5) * (-rSges(5,1) * t29 - rSges(5,2) * t31) - m(6) * (-pkin(4) * t29 + t40)), -m(6) * (g(3) * t43 + t52 * t40)];
taug = t1(:);

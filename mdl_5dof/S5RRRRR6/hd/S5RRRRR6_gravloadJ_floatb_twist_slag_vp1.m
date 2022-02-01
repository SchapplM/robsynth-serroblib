% Calculate Gravitation load on the joints for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:40
% EndTime: 2022-01-20 12:07:41
% DurationCPUTime: 0.36s
% Computational Cost: add. (328->68), mult. (245->84), div. (0->0), fcn. (194->10), ass. (0->42)
t36 = qJ(3) + qJ(4);
t30 = cos(t36);
t32 = qJ(5) + t36;
t25 = sin(t32);
t26 = cos(t32);
t55 = t26 * rSges(6,1) - t25 * rSges(6,2);
t54 = pkin(4) * t30 + t55;
t75 = rSges(4,3) + pkin(7);
t28 = sin(t36);
t56 = t30 * rSges(5,1) - t28 * rSges(5,2);
t42 = -pkin(8) - pkin(7);
t74 = pkin(9) - t42 + rSges(6,3);
t73 = -t42 + rSges(5,3);
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t72 = t40 * rSges(4,1) - t38 * rSges(4,2);
t71 = -pkin(2) - t72;
t33 = t40 * pkin(3);
t27 = t33 + pkin(2);
t70 = -t27 - t54;
t69 = -t27 - t56;
t51 = -rSges(6,1) * t25 - rSges(6,2) * t26;
t68 = -pkin(4) * t28 + t51;
t37 = qJ(1) + qJ(2);
t29 = sin(t37);
t31 = cos(t37);
t67 = g(1) * t31 + g(2) * t29;
t63 = t38 * pkin(3);
t39 = sin(qJ(1));
t62 = t39 * pkin(1);
t57 = t31 * rSges(3,1) - t29 * rSges(3,2);
t53 = -t29 * rSges(3,1) - t31 * rSges(3,2);
t52 = -rSges(5,1) * t28 - rSges(5,2) * t30;
t50 = t75 * t29 - t71 * t31;
t49 = t70 * t29 + t74 * t31;
t48 = t71 * t29 + t75 * t31;
t47 = t73 * t29 - t69 * t31;
t46 = t74 * t29 - t70 * t31;
t45 = t69 * t29 + t73 * t31;
t41 = cos(qJ(1));
t34 = t41 * pkin(1);
t1 = [-m(2) * (g(1) * (-t39 * rSges(2,1) - t41 * rSges(2,2)) + g(2) * (t41 * rSges(2,1) - t39 * rSges(2,2))) - m(3) * (g(1) * (t53 - t62) + g(2) * (t34 + t57)) - m(4) * (g(1) * (t48 - t62) + g(2) * (t34 + t50)) - m(5) * (g(1) * (t45 - t62) + g(2) * (t34 + t47)) - m(6) * (g(1) * (t49 - t62) + g(2) * (t34 + t46)), -m(3) * (g(1) * t53 + g(2) * t57) - m(4) * (g(1) * t48 + g(2) * t50) - m(5) * (g(1) * t45 + g(2) * t47) - m(6) * (g(1) * t49 + g(2) * t46), (-m(4) * t72 - m(5) * (t33 + t56) - m(6) * (t33 + t54)) * g(3) + t67 * (-m(4) * (-rSges(4,1) * t38 - rSges(4,2) * t40) - m(5) * (t52 - t63) - m(6) * (-t63 + t68)), (-m(5) * t56 - m(6) * t54) * g(3) + t67 * (-m(5) * t52 - m(6) * t68), -m(6) * (g(3) * t55 + t67 * t51)];
taug = t1(:);

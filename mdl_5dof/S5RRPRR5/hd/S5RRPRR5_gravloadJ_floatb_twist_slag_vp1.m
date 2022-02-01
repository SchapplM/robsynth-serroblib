% Calculate Gravitation load on the joints for
% S5RRPRR5
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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:12
% EndTime: 2022-01-20 11:02:13
% DurationCPUTime: 0.28s
% Computational Cost: add. (276->61), mult. (199->76), div. (0->0), fcn. (157->10), ass. (0->36)
t35 = pkin(9) + qJ(4);
t28 = cos(t35);
t29 = qJ(5) + t35;
t23 = sin(t29);
t24 = cos(t29);
t51 = t24 * rSges(6,1) - t23 * rSges(6,2);
t69 = pkin(4) * t28 + t51;
t39 = -pkin(7) - qJ(3);
t68 = pkin(8) - t39 + rSges(6,3);
t67 = -t39 + rSges(5,3);
t27 = sin(t35);
t66 = t28 * rSges(5,1) - t27 * rSges(5,2);
t65 = qJ(3) + rSges(4,3);
t38 = cos(pkin(9));
t25 = t38 * pkin(3) + pkin(2);
t64 = -t25 - t69;
t63 = rSges(4,2) * sin(pkin(9)) - pkin(2) - rSges(4,1) * t38;
t62 = -t25 - t66;
t36 = qJ(1) + qJ(2);
t30 = sin(t36);
t31 = cos(t36);
t61 = g(1) * t31 + g(2) * t30;
t40 = sin(qJ(1));
t58 = t40 * pkin(1);
t52 = t31 * rSges(3,1) - t30 * rSges(3,2);
t50 = -t30 * rSges(3,1) - t31 * rSges(3,2);
t49 = -rSges(6,1) * t23 - rSges(6,2) * t24;
t48 = t65 * t30 - t63 * t31;
t47 = t64 * t30 + t68 * t31;
t46 = t63 * t30 + t65 * t31;
t45 = t67 * t30 - t62 * t31;
t44 = t68 * t30 - t64 * t31;
t43 = t62 * t30 + t67 * t31;
t41 = cos(qJ(1));
t33 = t41 * pkin(1);
t1 = [-m(2) * (g(1) * (-t40 * rSges(2,1) - t41 * rSges(2,2)) + g(2) * (t41 * rSges(2,1) - t40 * rSges(2,2))) - m(3) * (g(1) * (t50 - t58) + g(2) * (t33 + t52)) - m(4) * (g(1) * (t46 - t58) + g(2) * (t33 + t48)) - m(5) * (g(1) * (t43 - t58) + g(2) * (t33 + t45)) - m(6) * (g(1) * (t47 - t58) + g(2) * (t33 + t44)), -m(3) * (g(1) * t50 + g(2) * t52) - m(4) * (g(1) * t46 + g(2) * t48) - m(5) * (g(1) * t43 + g(2) * t45) - m(6) * (g(1) * t47 + g(2) * t44), (-m(4) - m(5) - m(6)) * (g(1) * t30 - g(2) * t31), (-m(5) * t66 - m(6) * t69) * g(3) + t61 * (-m(5) * (-rSges(5,1) * t27 - rSges(5,2) * t28) - m(6) * (-pkin(4) * t27 + t49)), -m(6) * (g(3) * t51 + t61 * t49)];
taug = t1(:);

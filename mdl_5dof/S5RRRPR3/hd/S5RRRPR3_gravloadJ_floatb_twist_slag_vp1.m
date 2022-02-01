% Calculate Gravitation load on the joints for
% S5RRRPR3
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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:28
% EndTime: 2022-01-20 11:42:29
% DurationCPUTime: 0.30s
% Computational Cost: add. (290->66), mult. (218->81), div. (0->0), fcn. (173->10), ass. (0->39)
t36 = qJ(3) + pkin(9);
t28 = sin(t36);
t29 = cos(t36);
t41 = cos(qJ(3));
t33 = t41 * pkin(3);
t74 = t29 * rSges(5,1) - t28 * rSges(5,2) + t33;
t30 = qJ(5) + t36;
t22 = sin(t30);
t23 = cos(t30);
t52 = t23 * rSges(6,1) - t22 * rSges(6,2);
t73 = pkin(4) * t29 + t33 + t52;
t72 = rSges(4,3) + pkin(7);
t38 = -qJ(4) - pkin(7);
t71 = pkin(8) - t38 + rSges(6,3);
t70 = -t38 + rSges(5,3);
t39 = sin(qJ(3));
t69 = t41 * rSges(4,1) - t39 * rSges(4,2);
t67 = -pkin(2) - t69;
t66 = -pkin(2) - t73;
t65 = -pkin(2) - t74;
t37 = qJ(1) + qJ(2);
t31 = sin(t37);
t32 = cos(t37);
t64 = g(1) * t32 + g(2) * t31;
t61 = t39 * pkin(3);
t40 = sin(qJ(1));
t60 = t40 * pkin(1);
t53 = t32 * rSges(3,1) - t31 * rSges(3,2);
t51 = -t31 * rSges(3,1) - t32 * rSges(3,2);
t50 = -rSges(6,1) * t22 - rSges(6,2) * t23;
t49 = t72 * t31 - t67 * t32;
t48 = t66 * t31 + t71 * t32;
t47 = t67 * t31 + t72 * t32;
t46 = t70 * t31 - t65 * t32;
t45 = t71 * t31 - t66 * t32;
t44 = t65 * t31 + t70 * t32;
t42 = cos(qJ(1));
t34 = t42 * pkin(1);
t1 = [-m(2) * (g(1) * (-t40 * rSges(2,1) - t42 * rSges(2,2)) + g(2) * (t42 * rSges(2,1) - t40 * rSges(2,2))) - m(3) * (g(1) * (t51 - t60) + g(2) * (t34 + t53)) - m(4) * (g(1) * (t47 - t60) + g(2) * (t34 + t49)) - m(5) * (g(1) * (t44 - t60) + g(2) * (t34 + t46)) - m(6) * (g(1) * (t48 - t60) + g(2) * (t34 + t45)), -m(3) * (g(1) * t51 + g(2) * t53) - m(4) * (g(1) * t47 + g(2) * t49) - m(5) * (g(1) * t44 + g(2) * t46) - m(6) * (g(1) * t48 + g(2) * t45), (-m(4) * t69 - m(5) * t74 - m(6) * t73) * g(3) + t64 * (-m(4) * (-rSges(4,1) * t39 - rSges(4,2) * t41) - m(5) * (-rSges(5,1) * t28 - rSges(5,2) * t29 - t61) - m(6) * (-pkin(4) * t28 + t50 - t61)), (-m(5) - m(6)) * (g(1) * t31 - g(2) * t32), -m(6) * (g(3) * t52 + t64 * t50)];
taug = t1(:);

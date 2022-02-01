% Calculate Gravitation load on the joints for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:36
% EndTime: 2022-01-23 09:20:37
% DurationCPUTime: 0.23s
% Computational Cost: add. (262->65), mult. (180->86), div. (0->0), fcn. (158->10), ass. (0->37)
t52 = rSges(6,3) + pkin(7);
t51 = -m(5) - m(6);
t32 = sin(qJ(1));
t50 = pkin(1) * t32;
t28 = qJ(1) + pkin(8);
t26 = qJ(3) + t28;
t22 = sin(t26);
t49 = g(1) * t22;
t23 = cos(t26);
t29 = sin(pkin(9));
t48 = t23 * t29;
t30 = cos(pkin(9));
t47 = t23 * t30;
t31 = sin(qJ(5));
t46 = t30 * t31;
t33 = cos(qJ(5));
t45 = t30 * t33;
t44 = t23 * pkin(3) + t22 * qJ(4);
t25 = cos(t28);
t34 = cos(qJ(1));
t27 = t34 * pkin(1);
t43 = pkin(2) * t25 + t27;
t16 = t23 * qJ(4);
t5 = t22 * t46 + t23 * t33;
t6 = -t22 * t45 + t23 * t31;
t42 = t6 * rSges(6,1) + t5 * rSges(6,2) + t16;
t41 = t23 * rSges(4,1) - rSges(4,2) * t22;
t24 = sin(t28);
t40 = -pkin(2) * t24 - t50;
t39 = -rSges(4,1) * t22 - rSges(4,2) * t23;
t7 = t22 * t33 - t23 * t46;
t8 = t22 * t31 + t23 * t45;
t38 = t8 * rSges(6,1) + t7 * rSges(6,2) + pkin(4) * t47 + t52 * t48 + t44;
t37 = rSges(5,1) * t47 - rSges(5,2) * t48 + t22 * rSges(5,3) + t44;
t36 = t23 * rSges(5,3) + t16 + (-rSges(5,1) * t30 + rSges(5,2) * t29 - pkin(3)) * t22;
t35 = (-pkin(4) * t30 - t52 * t29 - pkin(3)) * t49;
t1 = [-m(2) * (g(1) * (-t32 * rSges(2,1) - rSges(2,2) * t34) + g(2) * (rSges(2,1) * t34 - t32 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t24 - rSges(3,2) * t25 - t50) + g(2) * (rSges(3,1) * t25 - rSges(3,2) * t24 + t27)) - m(4) * (g(1) * (t39 + t40) + g(2) * (t41 + t43)) - m(5) * (g(1) * (t36 + t40) + g(2) * (t37 + t43)) - m(6) * (g(1) * (t40 + t42) + g(2) * (t38 + t43) + t35), (-m(3) - m(4) + t51) * g(3), -m(4) * (g(1) * t39 + g(2) * t41) - m(5) * (g(1) * t36 + g(2) * t37) - m(6) * (g(1) * t42 + g(2) * t38 + t35), t51 * (-g(2) * t23 + t49), -m(6) * (g(1) * (rSges(6,1) * t7 - rSges(6,2) * t8) + g(2) * (-rSges(6,1) * t5 + rSges(6,2) * t6) + g(3) * (-rSges(6,1) * t31 - rSges(6,2) * t33) * t29)];
taug = t1(:);

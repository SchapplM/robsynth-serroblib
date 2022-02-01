% Calculate Gravitation load on the joints for
% S5RPRRR5
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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:34
% EndTime: 2022-01-20 09:48:36
% DurationCPUTime: 0.26s
% Computational Cost: add. (253->56), mult. (165->70), div. (0->0), fcn. (125->10), ass. (0->33)
t30 = cos(qJ(4));
t27 = qJ(4) + qJ(5);
t22 = sin(t27);
t23 = cos(t27);
t41 = t23 * rSges(6,1) - rSges(6,2) * t22;
t56 = t30 * pkin(4) + t41;
t55 = pkin(7) + rSges(5,3);
t54 = pkin(8) + pkin(7) + rSges(6,3);
t28 = sin(qJ(4));
t53 = rSges(5,1) * t30 - rSges(5,2) * t28;
t52 = -pkin(3) - t53;
t51 = -pkin(3) - t56;
t26 = qJ(1) + pkin(9);
t21 = qJ(3) + t26;
t16 = sin(t21);
t17 = cos(t21);
t50 = g(1) * t17 + g(2) * t16;
t29 = sin(qJ(1));
t47 = t29 * pkin(1);
t20 = cos(t26);
t31 = cos(qJ(1));
t25 = t31 * pkin(1);
t43 = pkin(2) * t20 + t25;
t42 = t17 * rSges(4,1) - rSges(4,2) * t16;
t19 = sin(t26);
t40 = -pkin(2) * t19 - t47;
t39 = -rSges(4,1) * t16 - rSges(4,2) * t17;
t38 = -rSges(6,1) * t22 - rSges(6,2) * t23;
t37 = t55 * t16 - t52 * t17;
t36 = t52 * t16 + t55 * t17;
t35 = t54 * t16 - t51 * t17;
t34 = t51 * t16 + t54 * t17;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t29 - rSges(2,2) * t31) + g(2) * (rSges(2,1) * t31 - rSges(2,2) * t29)) - m(3) * (g(1) * (-rSges(3,1) * t19 - rSges(3,2) * t20 - t47) + g(2) * (rSges(3,1) * t20 - rSges(3,2) * t19 + t25)) - m(4) * (g(1) * (t39 + t40) + g(2) * (t42 + t43)) - m(5) * (g(1) * (t36 + t40) + g(2) * (t37 + t43)) - m(6) * (g(1) * (t34 + t40) + g(2) * (t35 + t43)), (-m(3) - m(4) - m(5) - m(6)) * g(3), -m(4) * (g(1) * t39 + g(2) * t42) - m(5) * (g(1) * t36 + g(2) * t37) - m(6) * (g(1) * t34 + g(2) * t35), (-m(5) * t53 - m(6) * t56) * g(3) + t50 * (-m(5) * (-rSges(5,1) * t28 - rSges(5,2) * t30) - m(6) * (-pkin(4) * t28 + t38)), -m(6) * (g(3) * t41 + t50 * t38)];
taug = t1(:);

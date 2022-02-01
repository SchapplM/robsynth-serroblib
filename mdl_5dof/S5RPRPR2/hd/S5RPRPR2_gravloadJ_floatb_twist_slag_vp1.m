% Calculate Gravitation load on the joints for
% S5RPRPR2
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:48
% EndTime: 2022-01-23 09:18:49
% DurationCPUTime: 0.20s
% Computational Cost: add. (224->52), mult. (138->64), div. (0->0), fcn. (104->10), ass. (0->29)
t51 = pkin(7) + qJ(4) + rSges(6,3);
t25 = pkin(9) + qJ(5);
t19 = sin(t25);
t21 = cos(t25);
t50 = rSges(6,1) * t21 - rSges(6,2) * t19;
t49 = rSges(5,3) + qJ(4);
t28 = cos(pkin(9));
t48 = rSges(5,2) * sin(pkin(9)) - pkin(3) - rSges(5,1) * t28;
t47 = -pkin(4) * t28 - pkin(3) - t50;
t46 = -m(5) - m(6);
t30 = sin(qJ(1));
t45 = pkin(1) * t30;
t26 = qJ(1) + pkin(8);
t22 = cos(t26);
t31 = cos(qJ(1));
t24 = t31 * pkin(1);
t40 = pkin(2) * t22 + t24;
t23 = qJ(3) + t26;
t16 = sin(t23);
t17 = cos(t23);
t39 = t17 * rSges(4,1) - rSges(4,2) * t16;
t20 = sin(t26);
t38 = -pkin(2) * t20 - t45;
t37 = -rSges(4,1) * t16 - rSges(4,2) * t17;
t35 = t49 * t16 - t48 * t17;
t34 = t48 * t16 + t49 * t17;
t33 = t51 * t16 - t47 * t17;
t32 = t47 * t16 + t51 * t17;
t1 = [-m(2) * (g(1) * (-t30 * rSges(2,1) - rSges(2,2) * t31) + g(2) * (rSges(2,1) * t31 - t30 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t20 - rSges(3,2) * t22 - t45) + g(2) * (rSges(3,1) * t22 - rSges(3,2) * t20 + t24)) - m(4) * (g(1) * (t37 + t38) + g(2) * (t39 + t40)) - m(5) * (g(1) * (t34 + t38) + g(2) * (t35 + t40)) - m(6) * (g(1) * (t32 + t38) + g(2) * (t33 + t40)), (-m(3) - m(4) + t46) * g(3), -m(4) * (g(1) * t37 + g(2) * t39) - m(5) * (g(1) * t34 + g(2) * t35) - m(6) * (g(1) * t32 + g(2) * t33), t46 * (g(1) * t16 - g(2) * t17), -m(6) * (g(3) * t50 + (g(1) * t17 + g(2) * t16) * (-rSges(6,1) * t19 - rSges(6,2) * t21))];
taug = t1(:);

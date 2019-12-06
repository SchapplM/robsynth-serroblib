% Calculate Gravitation load on the joints for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:04
% EndTime: 2019-12-05 17:06:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (255->47), mult. (126->57), div. (0->0), fcn. (90->8), ass. (0->29)
t40 = rSges(6,3) + pkin(8);
t21 = sin(qJ(5));
t22 = cos(qJ(5));
t39 = rSges(6,1) * t22 - t21 * rSges(6,2);
t38 = -pkin(4) - t39;
t20 = pkin(9) + qJ(2);
t17 = sin(t20);
t37 = pkin(2) * t17;
t19 = qJ(3) + t20;
t14 = sin(t19);
t36 = pkin(3) * t14;
t15 = cos(t19);
t33 = t15 * rSges(4,1) - rSges(4,2) * t14;
t16 = qJ(4) + t19;
t11 = sin(t16);
t12 = cos(t16);
t32 = t12 * rSges(5,1) - rSges(5,2) * t11;
t10 = pkin(3) * t15;
t31 = t10 + t32;
t30 = -rSges(4,1) * t14 - rSges(4,2) * t15;
t29 = -rSges(5,1) * t11 - rSges(5,2) * t12;
t27 = t40 * t11 - t38 * t12;
t26 = t29 - t36;
t25 = t10 + t27;
t24 = t38 * t11 + t40 * t12;
t23 = t24 - t36;
t18 = cos(t20);
t13 = pkin(2) * t18;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t17 - rSges(3,2) * t18) + g(2) * (rSges(3,1) * t18 - rSges(3,2) * t17)) - m(4) * (g(1) * (t30 - t37) + g(2) * (t13 + t33)) - m(5) * (g(1) * (t26 - t37) + g(2) * (t13 + t31)) - m(6) * (g(1) * (t23 - t37) + g(2) * (t13 + t25)), -m(4) * (g(1) * t30 + g(2) * t33) - m(5) * (g(1) * t26 + g(2) * t31) - m(6) * (g(1) * t23 + g(2) * t25), -m(5) * (g(1) * t29 + g(2) * t32) - m(6) * (g(1) * t24 + g(2) * t27), -m(6) * (g(3) * t39 + (g(1) * t12 + g(2) * t11) * (-rSges(6,1) * t21 - rSges(6,2) * t22))];
taug = t1(:);

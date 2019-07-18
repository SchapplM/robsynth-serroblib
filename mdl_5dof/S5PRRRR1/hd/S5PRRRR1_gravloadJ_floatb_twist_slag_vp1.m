% Calculate Gravitation load on the joints for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:20
% EndTime: 2019-07-18 13:28:21
% DurationCPUTime: 0.28s
% Computational Cost: add. (129->60), mult. (195->93), div. (0->0), fcn. (162->8), ass. (0->34)
t14 = qJ(3) + qJ(4);
t12 = sin(t14);
t13 = cos(t14);
t24 = t13 * rSges(5,1) - t12 * rSges(5,2);
t19 = cos(qJ(3));
t40 = pkin(2) * t19;
t46 = -t24 - t40;
t15 = sin(qJ(5));
t36 = rSges(6,2) * t15;
t45 = rSges(6,3) * t13 + t12 * t36;
t17 = sin(qJ(2));
t20 = cos(qJ(2));
t44 = g(1) * t20 + g(3) * t17;
t43 = t45 * t17;
t42 = t45 * t20;
t16 = sin(qJ(3));
t41 = pkin(2) * t16;
t18 = cos(qJ(5));
t37 = rSges(6,1) * t18;
t34 = t12 * rSges(6,3);
t32 = t15 * t20;
t31 = t17 * t15;
t30 = t17 * t18;
t29 = t18 * t20;
t28 = t12 * t37;
t25 = rSges(4,1) * t19 - rSges(4,2) * t16;
t23 = -rSges(5,1) * t12 - rSges(5,2) * t13;
t22 = -t34 + (t36 - t37) * t13;
t10 = t20 * t40;
t4 = t13 * t29 + t31;
t3 = -t13 * t32 + t30;
t2 = -t13 * t30 + t32;
t1 = t13 * t31 + t29;
t5 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-t17 * rSges(3,1) - rSges(3,2) * t20) + g(3) * (rSges(3,1) * t20 - t17 * rSges(3,2))) - m(4) * (g(1) * (rSges(4,3) * t20 - t25 * t17) + g(3) * (t17 * rSges(4,3) + t25 * t20)) - m(5) * (g(3) * t10 + (g(1) * rSges(5,3) + g(3) * t24) * t20 + (g(3) * rSges(5,3) + g(1) * t46) * t17) - m(6) * (g(1) * (rSges(6,1) * t2 + rSges(6,2) * t1 + (-t34 - t40) * t17) + g(3) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t20 * t34 + t10)), -m(6) * (g(1) * t42 + g(3) * t43) + (m(4) * t25 - m(5) * t46 - m(6) * (t22 - t40)) * g(2) + t44 * (-m(4) * (-rSges(4,1) * t16 - rSges(4,2) * t19) - m(5) * (t23 - t41) - m(6) * (-t28 - t41)), -m(5) * (-g(2) * t24 + t44 * t23) - m(6) * (g(1) * (-t20 * t28 + t42) + g(2) * t22 + g(3) * (-t17 * t28 + t43)), -m(6) * (g(1) * (rSges(6,1) * t3 - rSges(6,2) * t4) + g(3) * (-rSges(6,1) * t1 + rSges(6,2) * t2) + g(2) * (rSges(6,1) * t15 + rSges(6,2) * t18) * t12)];
taug  = t5(:);

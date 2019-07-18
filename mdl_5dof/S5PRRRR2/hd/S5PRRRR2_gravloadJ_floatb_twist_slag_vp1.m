% Calculate Gravitation load on the joints for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:11
% EndTime: 2019-07-18 13:30:12
% DurationCPUTime: 0.17s
% Computational Cost: add. (165->45), mult. (120->57), div. (0->0), fcn. (84->8), ass. (0->27)
t37 = rSges(6,3) + pkin(6);
t17 = sin(qJ(5));
t19 = cos(qJ(5));
t36 = -rSges(6,1) * t19 + rSges(6,2) * t17;
t16 = qJ(2) + qJ(3);
t12 = sin(t16);
t35 = pkin(3) * t12;
t18 = sin(qJ(2));
t34 = t18 * pkin(2);
t13 = cos(t16);
t31 = t13 * rSges(4,1) - rSges(4,2) * t12;
t14 = qJ(4) + t16;
t10 = sin(t14);
t11 = cos(t14);
t30 = t11 * rSges(5,1) - t10 * rSges(5,2);
t9 = pkin(3) * t13;
t29 = t30 + t9;
t28 = -rSges(4,1) * t12 - rSges(4,2) * t13;
t27 = -t10 * rSges(5,1) - t11 * rSges(5,2);
t25 = t36 * t10 + t37 * t11;
t24 = t37 * t10 - t36 * t11;
t23 = t24 + t9;
t22 = t27 - t35;
t21 = t25 - t35;
t20 = cos(qJ(2));
t15 = t20 * pkin(2);
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-t18 * rSges(3,1) - rSges(3,2) * t20) + g(2) * (rSges(3,1) * t20 - t18 * rSges(3,2))) - m(4) * (g(1) * (t28 - t34) + g(2) * (t15 + t31)) - m(5) * (g(1) * (t22 - t34) + g(2) * (t15 + t29)) - m(6) * (g(1) * (t21 - t34) + g(2) * (t15 + t23)), -m(4) * (g(1) * t28 + g(2) * t31) - m(5) * (g(1) * t22 + g(2) * t29) - m(6) * (g(1) * t21 + g(2) * t23), -m(5) * (g(1) * t27 + g(2) * t30) - m(6) * (g(1) * t25 + g(2) * t24), -m(6) * (-g(3) * t36 + (g(1) * t11 + g(2) * t10) * (-rSges(6,1) * t17 - rSges(6,2) * t19))];
taug  = t1(:);

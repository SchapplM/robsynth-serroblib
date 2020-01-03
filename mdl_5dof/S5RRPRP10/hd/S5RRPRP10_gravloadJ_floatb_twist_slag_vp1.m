% Calculate Gravitation load on the joints for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:25
% EndTime: 2019-12-31 20:09:27
% DurationCPUTime: 0.51s
% Computational Cost: add. (171->106), mult. (349->148), div. (0->0), fcn. (322->6), ass. (0->45)
t53 = rSges(6,1) + pkin(4);
t47 = rSges(5,3) + pkin(7);
t39 = rSges(6,3) + qJ(5) + pkin(7);
t22 = cos(qJ(1));
t19 = sin(qJ(1));
t51 = g(2) * t19;
t28 = g(1) * t22 + t51;
t21 = cos(qJ(2));
t50 = g(3) * t21;
t20 = cos(qJ(4));
t49 = t20 * pkin(4);
t13 = t21 * pkin(2);
t18 = sin(qJ(2));
t46 = t18 * t22;
t17 = sin(qJ(4));
t45 = t19 * t17;
t44 = t19 * t20;
t43 = t21 * rSges(4,2);
t42 = t21 * t22;
t41 = t22 * t17;
t40 = t22 * t20;
t11 = t18 * qJ(3);
t38 = t11 + t13;
t37 = t22 * pkin(1) + t19 * pkin(6);
t36 = qJ(3) * t21;
t35 = -pkin(2) - t47;
t34 = t18 * t17 * pkin(4);
t33 = -pkin(2) - t39;
t32 = pkin(2) * t42 + t22 * t11 + t37;
t31 = -pkin(1) - t11;
t30 = g(1) * t35;
t29 = g(1) * t33;
t27 = t21 * rSges(3,1) - t18 * rSges(3,2);
t25 = rSges(5,1) * t17 + rSges(5,2) * t20;
t2 = t18 * t40 - t45;
t4 = t18 * t44 + t41;
t24 = rSges(6,2) * t20 + t53 * t17;
t6 = t19 * t36;
t8 = t22 * t36;
t23 = g(1) * t8 + g(2) * t6 + g(3) * t38;
t14 = t22 * pkin(6);
t10 = pkin(3) + t49;
t5 = -t18 * t45 + t40;
t3 = t18 * t41 + t44;
t1 = [-m(2) * (g(1) * (-t19 * rSges(2,1) - t22 * rSges(2,2)) + g(2) * (t22 * rSges(2,1) - t19 * rSges(2,2))) - m(3) * (g(1) * (t22 * rSges(3,3) + t14) + g(2) * (rSges(3,1) * t42 - rSges(3,2) * t46 + t37) + (g(1) * (-pkin(1) - t27) + g(2) * rSges(3,3)) * t19) - m(4) * (g(1) * (t22 * rSges(4,1) + t14) + g(2) * (-rSges(4,2) * t42 + rSges(4,3) * t46 + t32) + (g(1) * (-t18 * rSges(4,3) - t13 + t31 + t43) + g(2) * rSges(4,1)) * t19) - m(5) * (g(1) * (t5 * rSges(5,1) - t4 * rSges(5,2) + t22 * pkin(3) + t14) + g(2) * (t3 * rSges(5,1) + t2 * rSges(5,2) + t47 * t42 + t32) + (g(2) * pkin(3) + g(1) * t31 + t21 * t30) * t19) - m(6) * (g(1) * (t5 * rSges(6,1) - t4 * rSges(6,2) + t14) + g(2) * (t3 * rSges(6,1) + t2 * rSges(6,2) + t32) + (g(1) * t10 + g(2) * (t39 * t21 + t34)) * t22 + (g(1) * (t31 - t34) + g(2) * t10 + t21 * t29) * t19), -m(3) * (g(3) * t27 + t28 * (-rSges(3,1) * t18 - rSges(3,2) * t21)) - m(4) * (g(1) * (rSges(4,3) * t42 + t8) + g(2) * (t19 * t21 * rSges(4,3) + t6) + g(3) * (t38 - t43) + (g(3) * rSges(4,3) + t28 * (rSges(4,2) - pkin(2))) * t18) - m(5) * ((g(3) * t47 + t28 * t25) * t21 + (g(3) * t25 + t22 * t30 + t35 * t51) * t18 + t23) - m(6) * ((g(3) * t39 + t28 * t24) * t21 + (g(3) * t24 + t22 * t29 + t33 * t51) * t18 + t23), (-m(4) - m(5) - m(6)) * (t28 * t18 - t50), -m(5) * (g(1) * (t2 * rSges(5,1) - t3 * rSges(5,2)) + g(2) * (t4 * rSges(5,1) + t5 * rSges(5,2))) - m(6) * (g(1) * (-t3 * rSges(6,2) + t53 * t2) + g(2) * (t5 * rSges(6,2) + t53 * t4)) + (-m(5) * (-rSges(5,1) * t20 + rSges(5,2) * t17) - m(6) * (-rSges(6,1) * t20 + rSges(6,2) * t17 - t49)) * t50, -m(6) * (g(3) * t18 + t28 * t21)];
taug = t1(:);

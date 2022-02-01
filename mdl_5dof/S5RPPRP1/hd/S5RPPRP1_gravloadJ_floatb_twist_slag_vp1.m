% Calculate Gravitation load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:13
% EndTime: 2022-01-23 09:12:14
% DurationCPUTime: 0.32s
% Computational Cost: add. (182->70), mult. (194->98), div. (0->0), fcn. (179->8), ass. (0->28)
t37 = rSges(6,1) + pkin(4);
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t36 = rSges(4,1) * t15 - rSges(4,2) * t14;
t35 = pkin(3) * t15 + (rSges(5,3) + pkin(6)) * t14;
t17 = sin(qJ(4));
t33 = pkin(4) * t17;
t13 = qJ(1) + pkin(7);
t10 = sin(t13);
t32 = g(1) * t10;
t11 = cos(t13);
t31 = g(2) * t11;
t18 = sin(qJ(1));
t30 = t18 * pkin(1);
t26 = t15 * t17;
t19 = cos(qJ(4));
t25 = t15 * t19;
t24 = -m(4) - m(5) - m(6);
t20 = cos(qJ(1));
t12 = t20 * pkin(1);
t23 = t11 * pkin(2) + t10 * qJ(3) + t12;
t22 = t11 * qJ(3) - t30;
t3 = t10 * t19 - t11 * t26;
t1 = t10 * t26 + t11 * t19;
t21 = t15 * (t19 * pkin(4) + pkin(3)) + (rSges(6,3) + qJ(5) + pkin(6)) * t14;
t4 = t10 * t17 + t11 * t25;
t2 = -t10 * t25 + t11 * t17;
t5 = [-m(2) * (g(1) * (-t18 * rSges(2,1) - t20 * rSges(2,2)) + g(2) * (t20 * rSges(2,1) - t18 * rSges(2,2))) - m(3) * (g(1) * (-t10 * rSges(3,1) - t11 * rSges(3,2) - t30) + g(2) * (t11 * rSges(3,1) - t10 * rSges(3,2) + t12)) - m(4) * (g(1) * (t11 * rSges(4,3) + t22) + g(2) * (t36 * t11 + t23) + (g(1) * (-pkin(2) - t36) + g(2) * rSges(4,3)) * t10) - m(5) * (g(1) * (t2 * rSges(5,1) + t1 * rSges(5,2) + t22) + g(2) * (t4 * rSges(5,1) + t3 * rSges(5,2) + t23) + t35 * t31 + (-pkin(2) - t35) * t32) - m(6) * (g(1) * (t2 * rSges(6,1) + t1 * rSges(6,2) + t22) + g(2) * (t4 * rSges(6,1) + t3 * rSges(6,2) + t23) + (g(1) * t33 + g(2) * t21) * t11 + (g(1) * (-pkin(2) - t21) + g(2) * t33) * t10), (-m(3) + t24) * g(3), t24 * (-t31 + t32), -m(5) * (g(1) * (t3 * rSges(5,1) - rSges(5,2) * t4) + g(2) * (-rSges(5,1) * t1 + rSges(5,2) * t2)) - m(6) * (g(1) * (-t4 * rSges(6,2) + t37 * t3) + g(2) * (t2 * rSges(6,2) - t37 * t1)) + (-m(5) * (-rSges(5,1) * t17 - rSges(5,2) * t19) - m(6) * (-rSges(6,1) * t17 - rSges(6,2) * t19 - t33)) * g(3) * t14, -m(6) * (-g(3) * t15 + (g(1) * t11 + g(2) * t10) * t14)];
taug = t5(:);

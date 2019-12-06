% Calculate Gravitation load on the joints for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:34
% EndTime: 2019-12-05 15:10:35
% DurationCPUTime: 0.30s
% Computational Cost: add. (135->58), mult. (322->89), div. (0->0), fcn. (355->8), ass. (0->36)
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t30 = rSges(6,3) + qJ(5);
t40 = rSges(6,1) + pkin(4);
t27 = t30 * t23 + t40 * t25;
t19 = sin(pkin(8));
t26 = cos(qJ(3));
t35 = t19 * t26;
t41 = g(3) * pkin(6) * t35;
t39 = rSges(6,2) + pkin(6);
t38 = rSges(5,3) + pkin(6);
t37 = t19 * t23;
t36 = t19 * t25;
t20 = sin(pkin(7));
t24 = sin(qJ(3));
t34 = t20 * t24;
t33 = t20 * t26;
t22 = cos(pkin(7));
t32 = t22 * t24;
t31 = t22 * t26;
t29 = -m(3) - m(4) - m(5) - m(6);
t28 = rSges(5,1) * t25 - rSges(5,2) * t23;
t21 = cos(pkin(8));
t12 = -t21 * t23 + t25 * t35;
t11 = t21 * t25 + t23 * t35;
t10 = t21 * t31 + t34;
t9 = -t21 * t32 + t33;
t8 = t21 * t33 - t32;
t7 = -t21 * t34 - t31;
t6 = t9 * pkin(3);
t5 = t7 * pkin(3);
t4 = t10 * t25 + t22 * t37;
t3 = t10 * t23 - t22 * t36;
t2 = t20 * t37 + t8 * t25;
t1 = -t20 * t36 + t8 * t23;
t13 = [(-m(2) + t29) * g(3), t29 * (g(1) * t20 - g(2) * t22), -m(4) * (g(1) * (t9 * rSges(4,1) - t10 * rSges(4,2)) + g(2) * (t7 * rSges(4,1) - t8 * rSges(4,2))) - m(5) * (g(1) * (t38 * t10 + t28 * t9 + t6) + g(2) * (t28 * t7 + t38 * t8 + t5) + t41) - m(6) * (g(1) * (t39 * t10 + t27 * t9 + t6) + g(2) * (t27 * t7 + t39 * t8 + t5) + t41) + ((m(4) * rSges(4,2) - m(5) * rSges(5,3) - m(6) * rSges(6,2)) * t26 + (m(4) * rSges(4,1) - m(5) * (-pkin(3) - t28) - m(6) * (-pkin(3) - t27)) * t24) * g(3) * t19, -m(5) * (g(1) * (-t3 * rSges(5,1) - t4 * rSges(5,2)) + g(2) * (-t1 * rSges(5,1) - t2 * rSges(5,2)) + g(3) * (-t11 * rSges(5,1) - t12 * rSges(5,2))) - m(6) * (g(1) * (-t40 * t3 + t30 * t4) + g(2) * (-t40 * t1 + t30 * t2) + g(3) * (-t40 * t11 + t30 * t12)), -m(6) * (g(1) * t3 + g(2) * t1 + g(3) * t11)];
taug = t13(:);

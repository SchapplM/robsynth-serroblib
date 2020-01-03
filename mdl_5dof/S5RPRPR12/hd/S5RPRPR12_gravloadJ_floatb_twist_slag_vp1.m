% Calculate Gravitation load on the joints for
% S5RPRPR12
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:26
% EndTime: 2019-12-31 18:29:27
% DurationCPUTime: 0.37s
% Computational Cost: add. (236->73), mult. (258->101), div. (0->0), fcn. (233->10), ass. (0->38)
t15 = pkin(8) + qJ(3);
t11 = sin(t15);
t13 = cos(t15);
t16 = sin(pkin(9));
t18 = cos(pkin(9));
t29 = rSges(5,1) * t18 - rSges(5,2) * t16 + pkin(3);
t33 = rSges(5,3) + qJ(4);
t48 = t11 * t33 + t29 * t13;
t45 = rSges(6,3) + pkin(7) + qJ(4);
t22 = sin(qJ(1));
t23 = cos(qJ(1));
t44 = g(1) * t23 + g(2) * t22;
t14 = pkin(9) + qJ(5);
t10 = sin(t14);
t12 = cos(t14);
t8 = t18 * pkin(4) + pkin(3);
t43 = m(5) * t29 + m(6) * (rSges(6,1) * t12 - rSges(6,2) * t10 + t8) + m(4) * rSges(4,1);
t19 = cos(pkin(8));
t9 = t19 * pkin(2) + pkin(1);
t6 = t23 * t9;
t42 = g(2) * t6;
t41 = -m(5) - m(6);
t37 = t22 * t13;
t36 = t23 * t13;
t21 = -pkin(6) - qJ(2);
t35 = rSges(4,3) - t21;
t34 = rSges(3,3) + qJ(2);
t32 = pkin(4) * t16 - t21;
t31 = t13 * rSges(4,1) - t11 * rSges(4,2);
t30 = rSges(3,1) * t19 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t27 = t16 * rSges(5,1) + t18 * rSges(5,2) - t21;
t26 = t11 * t45 + t13 * t8;
t25 = m(4) * rSges(4,2) - m(5) * t33 - m(6) * t45;
t5 = t22 * t10 + t12 * t36;
t4 = -t10 * t36 + t22 * t12;
t3 = t23 * t10 - t12 * t37;
t2 = t10 * t37 + t23 * t12;
t1 = [-m(2) * (g(1) * (-t22 * rSges(2,1) - t23 * rSges(2,2)) + g(2) * (t23 * rSges(2,1) - t22 * rSges(2,2))) - m(3) * ((g(1) * t34 + g(2) * t30) * t23 + (-g(1) * t30 + g(2) * t34) * t22) - m(4) * (t42 + (g(1) * t35 + g(2) * t31) * t23 + (g(1) * (-t31 - t9) + g(2) * t35) * t22) - m(5) * (t42 + (g(1) * t27 + t48 * g(2)) * t23 + (g(2) * t27 + (-t9 - t48) * g(1)) * t22) - m(6) * (g(1) * (t3 * rSges(6,1) + t2 * rSges(6,2)) + g(2) * (t5 * rSges(6,1) + t4 * rSges(6,2) + t6) + (g(1) * t32 + g(2) * t26) * t23 + (g(1) * (-t26 - t9) + g(2) * t32) * t22), (-m(3) - m(4) + t41) * (g(1) * t22 - g(2) * t23), (t25 * t11 - t13 * t43) * g(3) + t44 * (t11 * t43 + t25 * t13), t41 * (-g(3) * t13 + t11 * t44), -m(6) * (g(1) * (t4 * rSges(6,1) - t5 * rSges(6,2)) + g(2) * (-t2 * rSges(6,1) + t3 * rSges(6,2)) + g(3) * (-rSges(6,1) * t10 - rSges(6,2) * t12) * t11)];
taug = t1(:);

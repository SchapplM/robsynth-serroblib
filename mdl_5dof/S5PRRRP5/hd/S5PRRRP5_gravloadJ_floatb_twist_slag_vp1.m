% Calculate Gravitation load on the joints for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:46
% EndTime: 2019-12-05 16:47:47
% DurationCPUTime: 0.37s
% Computational Cost: add. (208->78), mult. (300->113), div. (0->0), fcn. (278->8), ass. (0->40)
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t57 = g(1) * t24 + g(2) * t23;
t29 = -pkin(7) - pkin(6);
t22 = qJ(3) + qJ(4);
t18 = sin(t22);
t19 = cos(t22);
t28 = cos(qJ(2));
t46 = t23 * t28;
t10 = t24 * t18 - t19 * t46;
t9 = -t18 * t46 - t24 * t19;
t56 = t9 * rSges(6,1) + t10 * rSges(6,2);
t55 = t9 * rSges(5,1) + t10 * rSges(5,2);
t45 = t24 * t28;
t11 = -t18 * t45 + t23 * t19;
t12 = -t23 * t18 - t19 * t45;
t54 = t11 * rSges(6,1) + t12 * rSges(6,2);
t53 = t11 * rSges(5,1) + t12 * rSges(5,2);
t52 = pkin(4) * t18;
t26 = sin(qJ(2));
t49 = g(3) * t26;
t25 = sin(qJ(3));
t48 = t25 * pkin(3);
t47 = rSges(4,3) + pkin(6);
t44 = t25 * t28;
t27 = cos(qJ(3));
t43 = t27 * t28;
t42 = rSges(5,3) - t29;
t41 = rSges(6,3) + qJ(5) - t29;
t20 = t27 * pkin(3);
t15 = pkin(4) * t19 + t20;
t39 = -rSges(5,1) * t18 - rSges(5,2) * t19;
t38 = -rSges(6,1) * t18 - rSges(6,2) * t19;
t37 = t27 * rSges(4,1) - t25 * rSges(4,2) + pkin(2);
t36 = t23 * t27 - t24 * t44;
t35 = -t23 * t44 - t24 * t27;
t34 = rSges(5,1) * t19 - rSges(5,2) * t18 + pkin(2) + t20;
t33 = rSges(6,1) * t19 - rSges(6,2) * t18 + pkin(2) + t15;
t14 = -t48 - t52;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(3) * (t28 * rSges(3,1) - t26 * rSges(3,2)) + t57 * (-rSges(3,1) * t26 - rSges(3,2) * t28)) - m(4) * (g(3) * (t47 * t26 + t37 * t28) + t57 * (-t37 * t26 + t47 * t28)) - m(5) * (g(3) * (t42 * t26 + t34 * t28) + t57 * (-t34 * t26 + t42 * t28)) - m(6) * (g(3) * (t41 * t26 + t33 * t28) + t57 * (-t33 * t26 + t41 * t28)), -m(4) * (g(1) * (t36 * rSges(4,1) + (-t23 * t25 - t24 * t43) * rSges(4,2)) + g(2) * (t35 * rSges(4,1) + (-t23 * t43 + t24 * t25) * rSges(4,2))) - m(5) * (g(1) * (pkin(3) * t36 + t53) + g(2) * (pkin(3) * t35 + t55)) - m(6) * (g(1) * (t14 * t45 + t23 * t15 + t54) + g(2) * (t14 * t46 - t24 * t15 + t56)) + (-m(4) * (-rSges(4,1) * t25 - rSges(4,2) * t27) - m(5) * (t39 - t48) - m(6) * (t14 + t38)) * t49, -m(5) * (g(1) * t53 + g(2) * t55) - m(6) * (g(1) * (pkin(4) * t11 + t54) + g(2) * (pkin(4) * t9 + t56)) + (-m(5) * t39 - m(6) * (t38 - t52)) * t49, -m(6) * (-g(3) * t28 + t57 * t26)];
taug = t1(:);

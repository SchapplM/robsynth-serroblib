% Calculate Gravitation load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:21:59
% EndTime: 2019-12-05 18:22:01
% DurationCPUTime: 0.33s
% Computational Cost: add. (244->62), mult. (166->73), div. (0->0), fcn. (127->8), ass. (0->32)
t51 = rSges(6,1) + pkin(4);
t20 = sin(qJ(4));
t22 = cos(qJ(4));
t50 = -t20 * rSges(6,2) + t51 * t22;
t49 = rSges(5,3) + pkin(7);
t48 = rSges(6,3) + qJ(5) + pkin(7);
t47 = -pkin(3) - t50;
t37 = t22 * rSges(5,1);
t46 = -pkin(3) - t37;
t18 = qJ(1) + qJ(2);
t15 = sin(t18);
t43 = pkin(2) * t15;
t16 = cos(t18);
t42 = pkin(2) * t16;
t21 = sin(qJ(1));
t41 = t21 * pkin(1);
t23 = cos(qJ(1));
t40 = t23 * pkin(1);
t39 = t20 * rSges(5,2);
t35 = -t16 * rSges(3,1) + t15 * rSges(3,2);
t14 = pkin(8) + t18;
t11 = sin(t14);
t12 = cos(t14);
t33 = t11 * t39 + t49 * t12 - t43;
t32 = -t15 * rSges(3,1) - t16 * rSges(3,2);
t30 = -t12 * rSges(4,1) + t11 * rSges(4,2) - t42;
t28 = -t11 * rSges(4,1) - t12 * rSges(4,2) - t43;
t27 = -t42 + (t39 + t46) * t12;
t26 = (-g(2) * t49 + g(3) * t46) * t11;
t25 = t47 * t11 + t48 * t12 - t43;
t24 = -t48 * t11 + t47 * t12 - t42;
t1 = [-m(2) * (g(2) * (-t23 * rSges(2,1) + t21 * rSges(2,2)) + g(3) * (-t21 * rSges(2,1) - t23 * rSges(2,2))) - m(3) * (g(2) * (t35 - t40) + g(3) * (t32 - t41)) - m(4) * (g(2) * (t30 - t40) + g(3) * (t28 - t41)) - m(5) * (g(2) * (t27 - t40) + g(3) * (t33 - t41) + t26) - m(6) * (g(2) * (t24 - t40) + g(3) * (t25 - t41)), -m(3) * (g(2) * t35 + g(3) * t32) - m(4) * (g(2) * t30 + g(3) * t28) - m(5) * (g(2) * t27 + g(3) * t33 + t26) - m(6) * (g(2) * t24 + g(3) * t25), (-m(4) - m(5) - m(6)) * g(1), (-m(5) * (t37 - t39) - m(6) * t50) * g(1) + (-g(2) * t11 + g(3) * t12) * (m(5) * (rSges(5,1) * t20 + rSges(5,2) * t22) + m(6) * (rSges(6,2) * t22 + t51 * t20)), -m(6) * (g(2) * t12 + g(3) * t11)];
taug = t1(:);

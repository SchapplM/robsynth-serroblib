% Calculate Gravitation load on the joints for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:14
% EndTime: 2019-12-05 17:59:16
% DurationCPUTime: 0.39s
% Computational Cost: add. (148->74), mult. (188->95), div. (0->0), fcn. (148->6), ass. (0->35)
t44 = rSges(6,1) + pkin(4);
t14 = qJ(3) + qJ(4);
t8 = sin(t14);
t9 = cos(t14);
t22 = -t9 * rSges(6,2) - t44 * t8;
t15 = sin(qJ(3));
t39 = pkin(3) * t15;
t45 = t22 - t39;
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t43 = g(1) * t16 - g(2) * t18;
t42 = pkin(4) * t9;
t19 = -pkin(7) - pkin(6);
t41 = rSges(5,2) * t8;
t40 = rSges(6,2) * t8;
t17 = cos(qJ(3));
t38 = pkin(3) * t17;
t35 = t16 * t9;
t34 = t18 * t8;
t33 = t18 * t9;
t30 = rSges(4,3) + pkin(6);
t29 = rSges(5,3) - t19;
t28 = rSges(6,3) + qJ(5) - t19;
t27 = t18 * pkin(1) + t16 * qJ(2);
t26 = -rSges(5,1) * t8 - rSges(5,2) * t9;
t24 = rSges(4,1) * t15 + rSges(4,2) * t17;
t11 = t18 * qJ(2);
t21 = g(1) * t11 + g(2) * t27;
t20 = -t26 + t39;
t7 = rSges(5,2) * t34;
t6 = rSges(6,2) * t34;
t5 = rSges(5,1) * t35;
t4 = rSges(6,1) * t35;
t2 = t38 + t42;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t16 - rSges(2,2) * t18) + g(2) * (rSges(2,1) * t18 - rSges(2,2) * t16)) - m(3) * (g(1) * (rSges(3,3) * t18 + t11 + (rSges(3,2) - pkin(1)) * t16) + g(2) * (-rSges(3,2) * t18 + rSges(3,3) * t16 + t27)) - m(4) * ((g(1) * t24 + g(2) * t30) * t18 + (g(1) * (-pkin(1) - t30) + g(2) * t24) * t16 + t21) - m(5) * ((g(1) * t20 + g(2) * t29) * t18 + (g(1) * (-pkin(1) - t29) + g(2) * t20) * t16 + t21) - m(6) * ((-g(1) * t45 + g(2) * t28) * t18 + (g(1) * (-pkin(1) - t28) - g(2) * t45) * t16 + t21), (-m(3) - m(4) - m(5) - m(6)) * t43, -m(4) * (-g(3) * t24 + t43 * (rSges(4,1) * t17 - rSges(4,2) * t15)) - m(5) * (g(1) * (t5 + (t38 - t41) * t16) + g(2) * (t7 + (-rSges(5,1) * t9 - t38) * t18) - g(3) * t20) - m(6) * (g(1) * (t4 + (t2 - t40) * t16) + g(2) * (t6 + (-rSges(6,1) * t9 - t2) * t18) + g(3) * t45), -m(5) * (g(1) * (-t16 * t41 + t5) + g(2) * (-rSges(5,1) * t33 + t7) + g(3) * t26) - m(6) * (g(1) * (t4 + (-t40 + t42) * t16) + g(2) * (-t33 * t44 + t6) + g(3) * t22), -m(6) * (g(1) * t18 + g(2) * t16)];
taug = t1(:);

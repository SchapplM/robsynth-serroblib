% Calculate Gravitation load on the joints for
% S5PRRRP2
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:37
% EndTime: 2019-12-05 16:41:38
% DurationCPUTime: 0.24s
% Computational Cost: add. (238->54), mult. (160->66), div. (0->0), fcn. (125->6), ass. (0->28)
t38 = rSges(6,1) + pkin(4);
t22 = sin(qJ(4));
t31 = rSges(6,3) + qJ(5);
t45 = t31 * t22;
t23 = cos(qJ(4));
t34 = t22 * rSges(5,2);
t44 = rSges(5,1) * t23 - t34;
t21 = pkin(8) + qJ(2);
t20 = qJ(3) + t21;
t16 = sin(t20);
t17 = cos(t20);
t43 = g(1) * t17 + g(2) * t16;
t42 = t38 * t23 + t45;
t18 = sin(t21);
t41 = pkin(2) * t18;
t35 = t17 * t23;
t13 = t17 * pkin(7);
t33 = t17 * rSges(6,2) + t13;
t32 = t17 * pkin(3) + t16 * pkin(7);
t30 = t17 * rSges(4,1) - rSges(4,2) * t16;
t29 = -rSges(4,1) * t16 - rSges(4,2) * t17;
t28 = t16 * rSges(6,2) + t17 * t45 + t38 * t35 + t32;
t27 = rSges(5,1) * t35 + t16 * rSges(5,3) - t17 * t34 + t32;
t26 = t17 * rSges(5,3) + t13 + (-pkin(3) - t44) * t16;
t25 = g(1) * (-pkin(3) - t42) * t16;
t19 = cos(t21);
t15 = pkin(2) * t19;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t18 - rSges(3,2) * t19) + g(2) * (rSges(3,1) * t19 - rSges(3,2) * t18)) - m(4) * (g(1) * (t29 - t41) + g(2) * (t15 + t30)) - m(5) * (g(1) * (t26 - t41) + g(2) * (t15 + t27)) - m(6) * (g(1) * (t33 - t41) + g(2) * (t15 + t28) + t25), -m(4) * (g(1) * t29 + g(2) * t30) - m(5) * (g(1) * t26 + g(2) * t27) - m(6) * (g(1) * t33 + g(2) * t28 + t25), (-m(5) * t44 - m(6) * t42) * g(3) + t43 * (-m(5) * (-rSges(5,1) * t22 - rSges(5,2) * t23) - m(6) * (-t38 * t22 + t31 * t23)), -m(6) * (-g(3) * t23 + t43 * t22)];
taug = t1(:);

% Calculate Gravitation load on the joints for
% S5PRRRP1
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:47
% EndTime: 2019-12-05 16:39:49
% DurationCPUTime: 0.20s
% Computational Cost: add. (214->53), mult. (140->63), div. (0->0), fcn. (105->6), ass. (0->26)
t42 = pkin(7) + rSges(5,3);
t41 = qJ(5) + pkin(7) + rSges(6,3);
t25 = cos(qJ(4));
t24 = sin(qJ(4));
t34 = t24 * rSges(5,2);
t40 = rSges(5,1) * t25 - t34;
t33 = t24 * rSges(6,2);
t39 = rSges(6,1) * t25 - t33;
t22 = pkin(8) + qJ(2);
t18 = sin(t22);
t38 = pkin(2) * t18;
t20 = qJ(3) + t22;
t16 = cos(t20);
t35 = t16 * t25;
t15 = sin(t20);
t32 = t16 * rSges(4,1) - rSges(4,2) * t15;
t31 = -rSges(4,1) * t15 - rSges(4,2) * t16;
t30 = rSges(5,1) * t35 + (pkin(3) - t34) * t16 + t42 * t15;
t29 = t42 * t16 + (-pkin(3) - t40) * t15;
t21 = t25 * pkin(4);
t17 = t21 + pkin(3);
t28 = rSges(6,1) * t35 + (t17 - t33) * t16 + t41 * t15;
t27 = t41 * t16 + (-t17 - t39) * t15;
t19 = cos(t22);
t14 = pkin(2) * t19;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t18 - rSges(3,2) * t19) + g(2) * (rSges(3,1) * t19 - rSges(3,2) * t18)) - m(4) * (g(1) * (t31 - t38) + g(2) * (t14 + t32)) - m(5) * (g(1) * (t29 - t38) + g(2) * (t14 + t30)) - m(6) * (g(1) * (t27 - t38) + g(2) * (t14 + t28)), -m(4) * (g(1) * t31 + g(2) * t32) - m(5) * (g(1) * t29 + g(2) * t30) - m(6) * (g(1) * t27 + g(2) * t28), (-m(5) * t40 - m(6) * (t21 + t39)) * g(3) + (g(1) * t16 + g(2) * t15) * (-m(5) * (-rSges(5,1) * t24 - rSges(5,2) * t25) - m(6) * (-rSges(6,2) * t25 + (-rSges(6,1) - pkin(4)) * t24)), -m(6) * (g(1) * t15 - g(2) * t16)];
taug = t1(:);

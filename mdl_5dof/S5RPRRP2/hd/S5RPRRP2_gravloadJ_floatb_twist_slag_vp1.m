% Calculate Gravitation load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:27:52
% EndTime: 2022-01-23 09:27:53
% DurationCPUTime: 0.24s
% Computational Cost: add. (225->62), mult. (154->78), div. (0->0), fcn. (117->8), ass. (0->30)
t45 = pkin(7) + rSges(5,3);
t44 = qJ(5) + pkin(7) + rSges(6,3);
t26 = sin(qJ(1));
t43 = pkin(1) * t26;
t27 = cos(qJ(4));
t42 = rSges(5,1) * t27;
t41 = rSges(6,1) * t27;
t23 = qJ(1) + pkin(8);
t20 = qJ(3) + t23;
t15 = sin(t20);
t25 = sin(qJ(4));
t40 = t15 * t25;
t16 = cos(t20);
t39 = t16 * t25;
t38 = t16 * t27;
t19 = cos(t23);
t28 = cos(qJ(1));
t22 = t28 * pkin(1);
t37 = pkin(2) * t19 + t22;
t36 = t16 * rSges(4,1) - rSges(4,2) * t15;
t18 = sin(t23);
t35 = -pkin(2) * t18 - t43;
t34 = -rSges(4,1) * t15 - rSges(4,2) * t16;
t33 = rSges(5,1) * t38 - rSges(5,2) * t39 + t16 * pkin(3) + t45 * t15;
t32 = rSges(5,2) * t40 + (-pkin(3) - t42) * t15 + t45 * t16;
t21 = t27 * pkin(4);
t17 = t21 + pkin(3);
t31 = rSges(6,1) * t38 - rSges(6,2) * t39 + t44 * t15 + t16 * t17;
t30 = rSges(6,2) * t40 + (-t17 - t41) * t15 + t44 * t16;
t1 = [-m(2) * (g(1) * (-t26 * rSges(2,1) - rSges(2,2) * t28) + g(2) * (rSges(2,1) * t28 - t26 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t18 - rSges(3,2) * t19 - t43) + g(2) * (rSges(3,1) * t19 - rSges(3,2) * t18 + t22)) - m(4) * (g(1) * (t34 + t35) + g(2) * (t36 + t37)) - m(5) * (g(1) * (t32 + t35) + g(2) * (t33 + t37)) - m(6) * (g(1) * (t30 + t35) + g(2) * (t31 + t37)), (-m(3) - m(4) - m(5) - m(6)) * g(3), -m(4) * (g(1) * t34 + g(2) * t36) - m(5) * (g(1) * t32 + g(2) * t33) - m(6) * (g(1) * t30 + g(2) * t31), (-m(5) * (-rSges(5,2) * t25 + t42) - m(6) * (-rSges(6,2) * t25 + t21 + t41)) * g(3) + (g(1) * t16 + g(2) * t15) * (-m(5) * (-rSges(5,1) * t25 - rSges(5,2) * t27) - m(6) * (-rSges(6,2) * t27 + (-rSges(6,1) - pkin(4)) * t25)), -m(6) * (g(1) * t15 - g(2) * t16)];
taug = t1(:);

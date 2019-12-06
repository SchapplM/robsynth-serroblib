% Calculate Gravitation load on the joints for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:22
% EndTime: 2019-12-05 15:16:23
% DurationCPUTime: 0.27s
% Computational Cost: add. (161->59), mult. (305->102), div. (0->0), fcn. (328->10), ass. (0->34)
t43 = rSges(5,3) + pkin(6);
t17 = sin(pkin(9));
t18 = sin(pkin(8));
t42 = t17 * t18;
t20 = cos(pkin(8));
t41 = t17 * t20;
t21 = sin(qJ(4));
t40 = t17 * t21;
t23 = cos(qJ(4));
t39 = t17 * t23;
t24 = cos(qJ(3));
t38 = t17 * t24;
t22 = sin(qJ(3));
t37 = t18 * t22;
t36 = t18 * t24;
t35 = t20 * t22;
t34 = t20 * t24;
t33 = rSges(6,3) + pkin(7) + pkin(6);
t32 = -m(3) - m(4) - m(5) - m(6);
t19 = cos(pkin(9));
t8 = t19 * t36 - t35;
t31 = t18 * t39 - t8 * t21;
t30 = t23 * rSges(5,1) - t21 * rSges(5,2) + pkin(3);
t10 = t19 * t34 + t37;
t29 = -t10 * t21 + t20 * t39;
t28 = -t19 * t23 - t21 * t38;
t16 = qJ(4) + qJ(5);
t14 = sin(t16);
t15 = cos(t16);
t27 = rSges(6,1) * t15 - rSges(6,2) * t14 + t23 * pkin(4) + pkin(3);
t26 = m(6) * (g(1) * ((-t10 * t14 + t15 * t41) * rSges(6,1) + (-t10 * t15 - t14 * t41) * rSges(6,2)) + g(2) * ((-t8 * t14 + t15 * t42) * rSges(6,1) + (-t14 * t42 - t8 * t15) * rSges(6,2)) + g(3) * ((-t14 * t38 - t19 * t15) * rSges(6,1) + (t19 * t14 - t15 * t38) * rSges(6,2)));
t9 = -t19 * t35 + t36;
t7 = -t19 * t37 - t34;
t1 = [(-m(2) + t32) * g(3), t32 * (g(1) * t18 - g(2) * t20), -m(4) * (g(1) * (t9 * rSges(4,1) - t10 * rSges(4,2)) + g(2) * (t7 * rSges(4,1) - t8 * rSges(4,2))) - m(5) * (g(1) * (t43 * t10 + t30 * t9) + g(2) * (t30 * t7 + t43 * t8)) - m(6) * (g(1) * (t33 * t10 + t27 * t9) + g(2) * (t27 * t7 + t33 * t8)) + ((m(4) * rSges(4,2) - m(5) * t43 - m(6) * t33) * t24 + (m(4) * rSges(4,1) + m(5) * t30 + m(6) * t27) * t22) * g(3) * t17, -m(5) * (g(1) * (t29 * rSges(5,1) + (-t10 * t23 - t20 * t40) * rSges(5,2)) + g(2) * (t31 * rSges(5,1) + (-t18 * t40 - t8 * t23) * rSges(5,2)) + g(3) * (t28 * rSges(5,1) + (t19 * t21 - t23 * t38) * rSges(5,2))) - t26 - m(6) * (g(1) * t29 + g(2) * t31 + g(3) * t28) * pkin(4), -t26];
taug = t1(:);

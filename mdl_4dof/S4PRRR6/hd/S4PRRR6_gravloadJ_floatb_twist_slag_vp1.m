% Calculate Gravitation load on the joints for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:40
% EndTime: 2019-12-31 16:34:41
% DurationCPUTime: 0.22s
% Computational Cost: add. (104->45), mult. (175->69), div. (0->0), fcn. (159->8), ass. (0->22)
t11 = sin(qJ(3));
t13 = cos(qJ(3));
t8 = qJ(3) + qJ(4);
t6 = sin(t8);
t7 = cos(t8);
t31 = m(4) * (t13 * rSges(4,1) - t11 * rSges(4,2) + pkin(2)) + m(5) * (rSges(5,1) * t7 - rSges(5,2) * t6 + pkin(3) * t13 + pkin(2)) + m(3) * rSges(3,1);
t10 = cos(pkin(7));
t14 = cos(qJ(2));
t9 = sin(pkin(7));
t26 = t14 * t9;
t30 = (-t10 * t7 - t6 * t26) * rSges(5,1) + (t10 * t6 - t7 * t26) * rSges(5,2);
t25 = t10 * t14;
t29 = (-t6 * t25 + t7 * t9) * rSges(5,1) + (-t7 * t25 - t6 * t9) * rSges(5,2);
t12 = sin(qJ(2));
t27 = g(3) * t12;
t24 = t11 * t14;
t23 = t13 * t14;
t22 = -rSges(5,1) * t6 - rSges(5,2) * t7;
t20 = -t10 * t13 - t9 * t24;
t19 = -t10 * t24 + t13 * t9;
t17 = m(3) * rSges(3,2) - m(4) * (rSges(4,3) + pkin(5)) - m(5) * (rSges(5,3) + pkin(6) + pkin(5));
t1 = [(-m(2) - m(3) - m(4) - m(5)) * g(3), (t17 * t12 - t31 * t14) * g(3) + (g(1) * t10 + g(2) * t9) * (t31 * t12 + t17 * t14), -m(4) * (g(1) * (t19 * rSges(4,1) + (-t10 * t23 - t11 * t9) * rSges(4,2)) + g(2) * (t20 * rSges(4,1) + (t10 * t11 - t9 * t23) * rSges(4,2))) - m(5) * (g(1) * (t19 * pkin(3) + t29) + g(2) * (t20 * pkin(3) + t30)) + (-m(4) * (-rSges(4,1) * t11 - rSges(4,2) * t13) - m(5) * (-pkin(3) * t11 + t22)) * t27, -m(5) * (g(1) * t29 + g(2) * t30 + t22 * t27)];
taug = t1(:);

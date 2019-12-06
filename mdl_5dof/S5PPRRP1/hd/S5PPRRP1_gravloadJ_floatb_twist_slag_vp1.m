% Calculate Gravitation load on the joints for
% S5PPRRP1
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:39
% EndTime: 2019-12-05 15:06:40
% DurationCPUTime: 0.21s
% Computational Cost: add. (133->37), mult. (170->55), div. (0->0), fcn. (148->6), ass. (0->21)
t29 = rSges(6,1) + pkin(4);
t10 = sin(pkin(7));
t11 = cos(pkin(7));
t28 = g(1) * t11 + g(2) * t10;
t13 = sin(qJ(4));
t14 = cos(qJ(4));
t27 = m(5) * (rSges(5,1) * t14 - rSges(5,2) * t13 + pkin(3)) + m(6) * (-rSges(6,2) * t13 + t29 * t14 + pkin(3)) + m(4) * rSges(4,1);
t23 = t10 * t13;
t22 = t10 * t14;
t21 = t11 * t13;
t20 = t11 * t14;
t19 = -m(3) - m(4) - m(5) - m(6);
t9 = pkin(8) + qJ(3);
t8 = cos(t9);
t3 = -t8 * t21 + t22;
t1 = -t8 * t23 - t20;
t16 = m(4) * rSges(4,2) - m(5) * (rSges(5,3) + pkin(6)) - m(6) * (rSges(6,3) + qJ(5) + pkin(6));
t7 = sin(t9);
t4 = -t8 * t20 - t23;
t2 = -t8 * t22 + t21;
t5 = [(-m(2) + t19) * g(3), t19 * (g(1) * t10 - g(2) * t11), (t16 * t7 - t27 * t8) * g(3) + t28 * (t16 * t8 + t27 * t7), -m(5) * (g(1) * (rSges(5,1) * t3 + rSges(5,2) * t4) + g(2) * (rSges(5,1) * t1 + rSges(5,2) * t2)) - m(6) * (g(1) * (t4 * rSges(6,2) + t29 * t3) + g(2) * (t2 * rSges(6,2) + t29 * t1)) + (-m(5) * (-rSges(5,1) * t13 - rSges(5,2) * t14) - m(6) * (-rSges(6,2) * t14 - t29 * t13)) * g(3) * t7, -m(6) * (-g(3) * t8 + t28 * t7)];
taug = t5(:);

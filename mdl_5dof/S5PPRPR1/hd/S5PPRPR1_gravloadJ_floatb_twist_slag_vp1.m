% Calculate Gravitation load on the joints for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PPRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:53
% EndTime: 2019-12-05 15:00:53
% DurationCPUTime: 0.17s
% Computational Cost: add. (131->33), mult. (140->47), div. (0->0), fcn. (117->8), ass. (0->17)
t11 = sin(pkin(7));
t13 = cos(pkin(7));
t27 = g(1) * t13 + g(2) * t11;
t12 = cos(pkin(9));
t8 = pkin(9) + qJ(5);
t4 = sin(t8);
t6 = cos(t8);
t26 = m(5) * (rSges(5,1) * t12 - rSges(5,2) * sin(pkin(9)) + pkin(3)) + m(6) * (rSges(6,1) * t6 - rSges(6,2) * t4 + pkin(4) * t12 + pkin(3)) + m(4) * rSges(4,1);
t25 = -m(5) - m(6);
t9 = pkin(8) + qJ(3);
t7 = cos(t9);
t21 = t11 * t7;
t20 = t13 * t7;
t19 = -m(3) - m(4) + t25;
t16 = m(4) * rSges(4,2) - m(5) * (rSges(5,3) + qJ(4)) - m(6) * (rSges(6,3) + pkin(6) + qJ(4));
t5 = sin(t9);
t1 = [(-m(2) + t19) * g(3), t19 * (g(1) * t11 - g(2) * t13), (t16 * t5 - t26 * t7) * g(3) + t27 * (t16 * t7 + t26 * t5), t25 * (-g(3) * t7 + t27 * t5), -m(6) * (g(1) * ((t11 * t6 - t4 * t20) * rSges(6,1) + (-t11 * t4 - t6 * t20) * rSges(6,2)) + g(2) * ((-t13 * t6 - t4 * t21) * rSges(6,1) + (t13 * t4 - t6 * t21) * rSges(6,2)) + g(3) * (-rSges(6,1) * t4 - rSges(6,2) * t6) * t5)];
taug = t1(:);

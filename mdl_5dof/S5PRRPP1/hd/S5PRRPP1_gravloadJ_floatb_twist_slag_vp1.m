% Calculate Gravitation load on the joints for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:11
% EndTime: 2019-12-05 16:06:13
% DurationCPUTime: 0.22s
% Computational Cost: add. (181->52), mult. (149->60), div. (0->0), fcn. (117->6), ass. (0->26)
t9 = pkin(7) + qJ(2);
t6 = cos(t9);
t30 = g(2) * t6;
t29 = rSges(6,1) + pkin(4);
t4 = sin(t9);
t28 = g(1) * t6 + g(2) * t4;
t27 = rSges(6,3) + qJ(5);
t10 = qJ(3) + pkin(8);
t5 = sin(t10);
t7 = cos(t10);
t15 = t27 * t5 + t29 * t7;
t13 = cos(qJ(3));
t8 = t13 * pkin(3);
t3 = t8 + pkin(2);
t25 = t3 * t30;
t23 = -m(5) - m(6);
t12 = sin(qJ(3));
t22 = pkin(3) * t12;
t21 = rSges(4,3) + pkin(6);
t11 = -qJ(4) - pkin(6);
t20 = rSges(6,2) - t11;
t19 = rSges(5,3) - t11;
t18 = rSges(5,1) * t7 - rSges(5,2) * t5;
t17 = rSges(4,1) * t13 - t12 * rSges(4,2);
t16 = pkin(2) + t17;
t1 = [(-m(2) - m(3) - m(4) + t23) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t4 - rSges(3,2) * t6) + g(2) * (rSges(3,1) * t6 - rSges(3,2) * t4)) - m(4) * ((g(1) * t21 + g(2) * t16) * t6 + (-g(1) * t16 + g(2) * t21) * t4) - m(5) * (t25 + (g(1) * t19 + g(2) * t18) * t6 + (g(1) * (-t18 - t3) + g(2) * t19) * t4) - m(6) * (t25 + (g(1) * t20 + g(2) * t15) * t6 + (g(1) * (-t15 - t3) + g(2) * t20) * t4), (-m(4) * t17 - m(5) * (t18 + t8) - m(6) * (t8 + t15)) * g(3) + t28 * (-m(4) * (-rSges(4,1) * t12 - rSges(4,2) * t13) - m(5) * (-rSges(5,1) * t5 - rSges(5,2) * t7 - t22) - m(6) * (t27 * t7 - t29 * t5 - t22)), t23 * (g(1) * t4 - t30), -m(6) * (-g(3) * t7 + t28 * t5)];
taug = t1(:);
